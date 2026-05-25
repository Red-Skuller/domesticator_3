import typer
import yaml
import logging
from pathlib import Path
from typing import List, Optional
from rich.logging import RichHandler
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn

from domestica.core import run_pipeline_with_config
from domestica.models import PipelineConfig

app = typer.Typer(help="Protein Designer Pipeline: Excel/FASTA -> Params -> DNA Opt -> Excel")
console = Console()

def setup_logging(verbose: bool):
    """Sets up terminal output using Rich."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, console=console)]
    )

@app.command()
def run(
    input_path: Optional[Path] = typer.Option(None, "--input", "-i", help="Input fasta (.fasta) or Excel file (.xlsx)"),
    output_path: Optional[Path] = typer.Option(None, "--output", "-o", help="Output Excel file (.xlsx)"),
    params: Optional[List[str]] = typer.Option(None, "--params", help="Select which protein parameters to calculate."),
    optimize: Optional[bool] = typer.Option(None, "--optimize", help="Perform codon optimization"),
    vector: Optional[Path] = typer.Option(None, "--vector", "-v", help="Genbank vector file for insertion"),
    nstruct: Optional[int] = typer.Option(None, "--nstruct", "-n", help="Number of optimized DNA structures"),
    ph: Optional[float] = typer.Option(None, "--ph", help="pH for net charge calculation"),
    name_col: Optional[str] = typer.Option(None, "--name-col", help="Column header for protein names"),
    seq_col: Optional[str] = typer.Option(None, "--seq-col", help="Column header for protein sequences"),
    idt_credentials_dir: Optional[str] = typer.Option(None, "--idt-credentials-dir", help="Path to IDT API credentials"),
    skip_idt: Optional[bool] = typer.Option(None, "--skip-idt", help="Skip IDT complexity checking"),
    idt_type: Optional[str] = typer.Option(None, "--idt-type", help="Type of sequence to query IDT for"),
    idt_threshold: Optional[float] = typer.Option(None, "--idt-threshold", help="IDT score threshold"),
    n_tag: Optional[str] = typer.Option(None, "--n-tag", help="N-terminal tag for analysis"),
    c_tag: Optional[str] = typer.Option(None, "--c-tag", help="C-terminal tag for analysis"),
    out_cols: Optional[List[str]] = typer.Option(None, "--out-cols", help="Order and selection of output columns"),
    config_file: Optional[Path] = typer.Option(None, "--config", "-c", help="Path to a YAML configuration file"),
    verbose: bool = typer.Option(False, "--verbose", help="Enable verbose logging")
):
    setup_logging(verbose)

    config_data = {}
    if config_file:
        logging.info(f"Loading configuration from {config_file}")
        with open(config_file, "r") as f:
            config_data = yaml.safe_load(f)

    # Merge CLI options into config_data, overriding if not None
    cli_options = {
        "input_path": input_path,
        "output_path": output_path,
        "params": params,
        "optimize": optimize,
        "vector_path": vector,
        "nstruct": nstruct,
        "ph": ph,
        "name_col": name_col,
        "seq_col": seq_col,
        "idt_credentials_dir": Path(idt_credentials_dir) if idt_credentials_dir else None,
        "skip_idt": skip_idt,
        "idt_type": idt_type,
        "idt_threshold": idt_threshold,
        "n_tag": n_tag,
        "c_tag": c_tag,
        "out_cols": out_cols,
    }
    for key, value in cli_options.items():
        if value is not None:
            config_data[key] = value

    try:
        config = PipelineConfig(**config_data)
    except Exception as e:
        logging.error(f"Configuration error: {e}")
        raise typer.Exit(code=1)

    if not config.params and not config.optimize:
        logging.warning("Neither --params nor --optimize was selected. The output will just mirror the input.")

    logging.info("Starting pipeline...")
    run_pipeline_with_config(config)

def main():
    app()

if __name__ == "__main__":
    main()
