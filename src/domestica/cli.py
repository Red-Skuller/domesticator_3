import os
import sys
import logging
import concurrent.futures
from pathlib import Path
from typing import Optional

import typer

from domestica.schema import PipelineConfig, SequenceRecord, ResultRow
from domestica.io import parse_input_file, write_output_file
from domestica.optimizer import optimize_sequence
from domestica.vendors.base import get_evaluator

app = typer.Typer(name="domestica", add_completion=False)
logger = logging.getLogger(__name__)

_worker_evaluator = None


def _init_worker(vendor_target: Optional[str], product: str, verbose: bool) -> None:
    global _worker_evaluator
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        stream=sys.stdout
    )
    logger.debug("Worker process initialized. Vendor: %s, Product: %s, Verbose: %s", vendor_target, product, verbose)
    if vendor_target:
        _worker_evaluator = get_evaluator(vendor_target, product)
    else:
        _worker_evaluator = None


def _worker_task(record: SequenceRecord, template_path: Path, min_length: int) -> ResultRow:
    global _worker_evaluator
    logger.info("Starting processing for record ID: %s", record.record_id)
    try:
        evaluator_func = _worker_evaluator.evaluate if _worker_evaluator else None

        naive, opt_seq, accepted, score, opt_rec = optimize_sequence(
            protein_sequence=record.protein_sequence,
            template_path=template_path,
            min_length=min_length,
            evaluator=evaluator_func
        )
        logger.info("Successfully optimized record ID: %s. Accepted: %s, Score: %s", record.record_id, accepted, score)
        return ResultRow(
            record_id=record.record_id, protein_sequence=record.protein_sequence,
            naive_dna_sequence=naive, optimized_sequence=opt_seq, vendor_score=score,
            accepted=accepted, status="SUCCESS", optimized_record=opt_rec
        )
    except Exception as e:
        logger.exception("Execution failed for record ID: %s due to an unhandled exception.", record.record_id)
        return ResultRow(
            record_id=record.record_id, protein_sequence=record.protein_sequence,
            status="FAILED", error_message=str(e)
        )


@app.command()
def optimize(
        input_path: Optional[Path] = typer.Argument(None, help="Path to inputs.", exists=True),
        output_path: Path = typer.Argument(..., help="Path for outputs."),
        template_path: Path = typer.Option(..., "--template", "-t", exists=True),
        vendor: Optional[str] = typer.Option(None, "--vendor", "-v"),
        product: str = typer.Option("eblocks", "--product", "-p"),
        max_workers: int = typer.Option(max(1, (os.cpu_count() or 2) - 1), "--workers", "-w"),
        min_length: int = typer.Option(300, "--min-length", "-m", help="Minimum sequence length in base pairs."),
        verbose: bool = typer.Option(False, "--verbose")
) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        stream=sys.stdout
    )

    logger.info("Initializing optimization pipeline workflow.")
    config = PipelineConfig(vendor_target=vendor, product=product, max_workers=max_workers, min_length=min_length)
    logger.debug("Pipeline Configuration: %s", config.model_dump())

    if input_path is not None:
        logger.debug("Input path provided: %s", input_path)
        records = parse_input_file(input_path)
    else:
        logger.info("No input path provided. Transitioning optimization execution target directly to base template.")
        records = [SequenceRecord(record_id="optimized_template", protein_sequence=None)]

    results = []

    logger.info("Spawning ProcessPoolExecutor with %d workers.", config.max_workers)
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=config.max_workers,
            initializer=_init_worker,
            initargs=(config.vendor_target, config.product, verbose)
    ) as executor:
        future_to_record = {
            executor.submit(_worker_task, r, template_path, config.min_length): r for r in records
        }
        for future in concurrent.futures.as_completed(future_to_record):
            rec = future_to_record[future]
            try:
                res = future.result()
                results.append(res)
                logger.debug("Retrieved future result execution status for record %s: %s", rec.record_id, res.status)
            except Exception:
                logger.exception("Critical error retrieving process future result for record %s.", rec.record_id)

    write_output_file(results, output_path)
    logger.info("Optimization process pipeline complete.")

    logger.info("Optimization Summary:")
    header = f"{'Record ID':<25} | {'Status':<10} | {'Vendor Accepted':<15} | {'Vendor Score':<15} | {'Length':<10}"
    logger.info(header)
    logger.info("-" * len(header))

    for res in results:
        seq_len = str(len(res.optimized_sequence)) if res.optimized_sequence else "N/A"
        score = f"{res.vendor_score:.2f}" if res.vendor_score is not None else "N/A"

        row_str = f"{str(res.record_id):<25} | {str(res.status):<10} | {str(res.accepted):<15} | {score:<15} | {seq_len:<10}"
        logger.info(row_str)

if __name__ == "__main__":
    app()