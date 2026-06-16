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


def _init_worker(vendor_target: Optional[str], product: str) -> None:
    global _worker_evaluator
    if vendor_target:
        _worker_evaluator = get_evaluator(vendor_target, product)
    else:
        _worker_evaluator = None


def _worker_task(record: SequenceRecord, template_path: Path) -> ResultRow:
    global _worker_evaluator
    try:
        evaluator_func = _worker_evaluator.evaluate if _worker_evaluator else None

        naive, opt_seq, accepted, score, opt_rec = optimize_sequence(
            protein_sequence=record.protein_sequence,
            template_path=template_path,
            evaluator=evaluator_func
        )
        return ResultRow(
            record_id=record.record_id, protein_sequence=record.protein_sequence,
            naive_dna_sequence=naive, optimized_sequence=opt_seq, vendor_score=score,
            accepted=accepted, status="SUCCESS", optimized_record=opt_rec
        )
    except Exception as e:
        return ResultRow(
            record_id=record.record_id, protein_sequence=record.protein_sequence,
            status="FAILED", error_message=str(e)
        )


@app.command()
def optimize(
        input_path: Optional[Path] = typer.Argument(None,
                                                    help="Path to inputs. If omitted, the template itself is optimized.",
                                                    exists=True),
        output_path: Path = typer.Argument(..., help="Path for outputs."),
        template_path: Path = typer.Option(..., "--template", "-t", exists=True),
        vendor: Optional[str] = typer.Option(None, "--vendor", "-v"),
        product: str = typer.Option("eblocks", "--product", "-p"),
        max_workers: int = typer.Option(max(1, (os.cpu_count() or 2) - 1), "--workers", "-w"),
        verbose: bool = typer.Option(False, "--verbose")
) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        stream=sys.stdout
    )

    config = PipelineConfig(vendor_target=vendor, product=product, max_workers=max_workers)

    if input_path is not None:
        records = parse_input_file(input_path)
    else:
        # Optimization execution target directly transitions to the base template file
        records = [SequenceRecord(record_id="optimized_template", protein_sequence=None)]

    results = []

    with concurrent.futures.ProcessPoolExecutor(
            max_workers=config.max_workers,
            initializer=_init_worker,
            initargs=(config.vendor_target, config.product)
    ) as executor:
        future_to_record = {
            executor.submit(_worker_task, r, template_path): r for r in records
        }
        for future in concurrent.futures.as_completed(future_to_record):
            results.append(future.result())
    print(results)
    write_output_file(results, output_path)


if __name__ == "__main__":
    app()