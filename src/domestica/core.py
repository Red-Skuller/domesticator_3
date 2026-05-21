import logging
from pathlib import Path
import pandas as pd
from typing import Optional

from dnachisel import NoSolutionError

from domestica import io_utils
from domestica import AA_analysis
from domestica import codon_opt
from domestica import vector_assembly
from domestica import idt

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
    if not skip_idt:
        user_info_file = idt.use_dir(idt_credentials_dir)
        idt_user_info = idt.get_user_info(user_info_file)

    logging.info(f"Reading input from {input_path}...")
    records = io_utils.read_input(input_path, name_col, seq_col)

    results = []

    for record in records:
        record_id = record["id"]
        protein_seq = record["sequence"]
        logging.info(f"Processing: {record_id}")
        # Pre-populate row data with the new variables to support dynamic column selection
        row_data = {
            "Name": record_id,
            "AA_seq": protein_seq,
            "N_tag": n_tag,
            "C_tag": c_tag,
            "insert_name": Path(vector_path).name,
        }

        # 1. Protein Parameters
        if params:
            analyzed_seq = n_tag + protein_seq + c_tag
            calculated_metrics = AA_analysis.calculate_selected_params(
                sequence=analyzed_seq,
                requested_params=params,
                ph=ph
            )
            row_data.update(calculated_metrics)
            if n_tag or c_tag:
                row_data["AA_Seq_Final"] = analyzed_seq

        # 2. Codon Optimization
        if optimize:
            if not vector_path:
                #Empty fallback vector file without flanking sequences and default constraints
                vector_path = Path(__file__).parent / "specifications" / "no_vector.gb"
                logging.info(f"No vector provided. Defaulting to: {vector_path}")
            try:
                vector_record = vector_assembly.load_vector_record(vector_path)
                naive_record = vector_assembly.make_naive_vector_record_by_seq(vector_record, protein_seq)
                candidate_solutions = []
                for n in range(nstruct):
                    logging.info(f"Generating optimized structure {n + 1}/{nstruct} for {record_id}")
                    try:
                        optimized_vector_solution = codon_opt.optimize_naive_record(naive_record)
                        optimized_vector_record = optimized_vector_solution.to_record()
                        fragment_to_synthesize = None
                        for feature in optimized_vector_record.features:
                            # Extract exact sub-fragment to order if tagged by domesticator flag
                            if feature.type == "domesticator" and feature.qualifiers.get('label') == ["synthesize"]:
                                fragment_to_synthesize = feature.extract(optimized_vector_record.seq)
                                break

                        # Fallback in case the vector lacks a domesticator tag
                        if fragment_to_synthesize is None:
                            fragment_to_synthesize = optimized_vector_record.seq
                            logging.warning(
                                f"No 'synthesize' feature found in vector. Outputting full vector sequence for {record_id}.")
                        sequence_str = str(fragment_to_synthesize)

                        # IDT Complexity Check
                        if not skip_idt:
                            try:
                                eval_score, issues = idt.get_complexity_score(sequence_str, idt_user_info, kind=idt_type)

                                print(f"SOLUTION {n + 1}:")
                                for issue in issues:
                                    print(f"{issue['Score']} {issue['Name']}")
                                print(f"Total Score: {eval_score}")
                            except Exception as e:
                                logging.error(f"Failed to query IDT complexity: {e}")
                            if eval_score < idt_threshold:
                                candidate_solutions.append({
                                    "seq": sequence_str,
                                    "score": eval_score,
                                    "record": optimized_vector_record  # <-- ADDED
                                })
                                break
                        else:
                            # Use dnachisel's internal objective evaluation score
                            eval_score = optimized_vector_solution.objectives_evaluations().scores_sum()

                        candidate_solutions.append({
                            "seq": sequence_str,
                            "score": eval_score,
                            "record": optimized_vector_record  # <-- ADDED
                        })
                    except NoSolutionError:
                        # Only skip THIS iteration, not the whole protein
                        logging.warning(f"No valid codon optimization solution found on attempt {n + 1}/{nstruct}.")
                        continue
                    except Exception as e:
                        logging.error(f"Optimization attempt {n + 1}/{nstruct} failed: {e}")
                        continue
                if candidate_solutions:
                    best_solution = min(candidate_solutions, key=lambda x: x["score"])
                    row_data["DNA_seq"] = best_solution["seq"]
                    row_data["DNA_record"] = best_solution["record"]  # <-- ADDED

                    # Log and store the corresponding score type
                    if not skip_idt:
                        row_data["IDT_Score"] = best_solution["score"]
                        logging.info(
                            f"Best structure for {record_id} selected with IDT Score: {best_solution['score']}")
                    else:
                        row_data["Dnachisel_Score"] = best_solution["score"]
                        logging.info(
                            f"Best structure for {record_id} selected with Dnachisel Score: {best_solution['score']}")
            except Exception as e:
                logging.error(f"Codon optimization failed for {record_id}: {e}")

        results.append(row_data)

    # 3. Write Output
    logging.info(f"Writing {len(results)} records to {output_path}...")
    io_utils.write_output(results, output_path, out_cols)
    logging.info("Pipeline completed successfully!")