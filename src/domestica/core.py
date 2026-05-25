import logging
from pathlib import Path
from typing import List, Optional, Dict, Any

from dnachisel import NoSolutionError

from domestica import io_utils
from domestica import AA_analysis
from domestica import codon_opt
from domestica import vector_assembly
from domestica import idt
from domestica.models import PipelineConfig, ResultRow

def ingest(config: PipelineConfig) -> List[Dict[str, str]]:
    """Reads input sequences from the specified path."""
    logging.info(f"Reading input from {config.input_path}...")
    return io_utils.read_input(config.input_path, config.name_col, config.seq_col)

def analyze(record: Dict[str, str], config: PipelineConfig) -> ResultRow:
    """Performs biochemical analysis on a single protein record."""
    record_id = record["id"]
    protein_seq = record["sequence"]

    row_data = {
        "Name": record_id,
        "AA_seq": protein_seq,
        "N_tag": config.n_tag,
        "C_tag": config.c_tag,
        "insert_name": config.vector_path.name if config.vector_path else "no_vector.gb",
    }

    if config.params:
        analyzed_seq = config.n_tag + protein_seq + config.c_tag
        calculated_metrics = AA_analysis.calculate_selected_params(
            sequence=analyzed_seq,
            requested_params=config.params,
            ph=config.ph
        )
        row_data.update(calculated_metrics)
        if config.n_tag or config.c_tag:
            row_data["AA_Seq_Final"] = analyzed_seq

    return ResultRow(**row_data)

from Bio import BiopythonParserWarning
import warnings

def optimize(row: ResultRow, config: PipelineConfig, idt_user_info: Optional[Dict] = None) -> ResultRow:
    """Performs codon optimization on a single protein record."""
    if not config.optimize:
        return row

    vector_path = config.vector_path
    if not vector_path:
        vector_path = Path(__file__).parent / "specifications" / "no_vector.gb"
        logging.info(f"No vector provided. Defaulting to: {vector_path}")
        row.insert_name = vector_path.name

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=BiopythonParserWarning)
            vector_record = vector_assembly.load_vector_record(vector_path)
        naive_record = vector_assembly.make_naive_vector_record_by_seq(vector_record, row.aa_seq)
        candidate_solutions = []

        optimizer = codon_opt.OPTIMIZERS_REGISTRY["dnachisel"]
        evaluator = idt.EVALUATORS_REGISTRY["idt"]

        for n in range(config.nstruct):
            logging.info(f"Generating optimized structure {n + 1}/{config.nstruct} for {row.name}")
            try:
                optimized_vector_solution = optimizer.optimize(naive_record)
                optimized_vector_record = optimized_vector_solution.to_record()

                fragment_to_synthesize = None
                for feature in optimized_vector_record.features:
                    if feature.type == "domesticator" and feature.qualifiers.get('label') == ["synthesize"]:
                        fragment_to_synthesize = feature.extract(optimized_vector_record.seq)
                        break

                if fragment_to_synthesize is None:
                    fragment_to_synthesize = optimized_vector_record.seq
                    logging.warning(f"No 'synthesize' feature found in vector. Outputting full vector sequence for {row.name}.")

                sequence_str = str(fragment_to_synthesize)

                if not config.skip_idt:
                    try:
                        eval_score, issues = evaluator.evaluate(
                            sequence_str,
                            user_info=idt_user_info,
                            kind=config.idt_type,
                            verbose=False # Or take from config if added
                        )
                        # Original code used print for IDT details
                        logging.debug(f"SOLUTION {n + 1}: Total Score: {eval_score}")
                    except Exception as e:
                        logging.error(f"Failed to query IDT complexity: {e}")
                        eval_score = float('inf') # Ensure it's not chosen if it fails

                    if eval_score < config.idt_threshold:
                        candidate_solutions.append({"seq": sequence_str, "score": eval_score})
                        break
                else:
                    eval_score = optimized_vector_solution.objectives_evaluations().scores_sum()

                candidate_solutions.append({"seq": sequence_str, "score": eval_score})
            except NoSolutionError:
                logging.warning(f"No valid codon optimization solution found on attempt {n + 1}/{config.nstruct}.")
                continue
            except Exception as e:
                logging.error(f"Optimization attempt {n + 1}/{config.nstruct} failed: {e}")
                continue

        if candidate_solutions:
            best_solution = min(candidate_solutions, key=lambda x: x["score"])
            row.dna_seq = best_solution["seq"]

            if not config.skip_idt:
                row.idt_score = best_solution["score"]
                logging.info(f"Best structure for {row.name} selected with IDT Score: {best_solution['score']}")
            else:
                row.dnachisel_score = best_solution["score"]
                logging.info(f"Best structure for {row.name} selected with Dnachisel Score: {best_solution['score']}")

    except Exception as e:
        logging.error(f"Codon optimization failed for {row.name}: {e}")

    return row

def export(results: List[ResultRow], config: PipelineConfig):
    """Writes the results to the specified output path."""
    logging.info(f"Writing {len(results)} records to {config.output_path}...")

    # Convert ResultRows back to dictionaries for pandas, using field aliases
    data = [row.model_dump(by_alias=True) for row in results]
    import pandas as pd
    df = pd.DataFrame(data)

    if config.out_cols:
        df = df.reindex(columns=config.out_cols)

    df.to_excel(config.output_path, index=False)
    logging.info("Pipeline completed successfully!")

from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn

def run_pipeline_with_config(config: PipelineConfig):
    """Executes the full pipeline using a PipelineConfig object."""
    idt_user_info = None
    if not config.skip_idt:
        user_info_file = idt.use_dir(config.idt_credentials_dir)
        idt_user_info = idt.get_user_info(user_info_file)

    records = ingest(config)
    results = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("[cyan]Processing sequences...", total=len(records))
        for record in records:
            progress.update(task, description=f"Processing: {record['id']}")
            row = analyze(record, config)
            row = optimize(row, config, idt_user_info)
            results.append(row)
            progress.advance(task)

    export(results, config)

def run_pipeline(
        input_path: Path,
        output_path: Path,
        params: list[str],
        optimize: bool,
        vector_path: Optional[Path],
        nstruct: int,
        skip_idt: bool,
        ph: float,
        name_col: str,
        seq_col: str,
        idt_type: str,
        idt_credentials_dir: str,
        idt_threshold: float,
        n_tag: str = "",
        c_tag: str = "",
        out_cols: Optional[list[str]] = None
):
    """Legacy entry point for run_pipeline."""
    config = PipelineConfig(
        input_path=input_path,
        output_path=output_path,
        params=params,
        optimize=optimize,
        vector_path=vector_path,
        nstruct=nstruct,
        skip_idt=skip_idt,
        ph=ph,
        name_col=name_col,
        seq_col=seq_col,
        idt_type=idt_type,
        idt_credentials_dir=Path(idt_credentials_dir),
        idt_threshold=idt_threshold,
        n_tag=n_tag,
        c_tag=c_tag,
        out_cols=out_cols
    )
    run_pipeline_with_config(config)
