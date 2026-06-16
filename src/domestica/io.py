import logging
import re
from pathlib import Path
from typing import List
import pandas as pd
from Bio import SeqIO
from domestica.schema import SequenceRecord, ResultRow

logger = logging.getLogger(__name__)


def parse_input_file(file_path: Path) -> List[SequenceRecord]:
    ext = file_path.suffix.lower()
    records = []
    logger.info("Parsing input sequence data from file: %s", file_path)

    if ext in (".fasta", ".fa"):
        logger.debug("Identified format as FASTA.")
        for record in SeqIO.parse(file_path, "fasta"):
            seq_str = str(record.seq).upper().strip()
            records.append(SequenceRecord(
                record_id=record.id,
                protein_sequence=seq_str if seq_str else None,
                metadata={"description": record.description}
            ))
    elif ext in (".xlsx", ".xls"):
        logger.debug("Identified format as Excel spreadsheet.")
        df = pd.read_excel(file_path)
        name_col = df.columns[0]
        seq_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
        logger.debug("Parsing columns: ID Column='%s', Sequence Column='%s'", name_col, seq_col)
        for index, row in df.iterrows():
            val = row.get(seq_col)
            prot_seq = str(val).upper().strip() if pd.notna(val) else None
            records.append(SequenceRecord(
                record_id=str(row.get(name_col, f"seq_{index}")),
                protein_sequence=prot_seq if prot_seq else None,
                metadata=row.to_dict()
            ))
    else:
        logger.error("Failed parsing input file: Unsupported format standard '%s'", ext)
        raise ValueError(f"Unsupported input format: {ext}")

    logger.info("Parsed %d total structural sequence records from file.", len(records))
    return records


def write_output_file(results: List[ResultRow], destination: Path) -> None:
    if not results:
        logger.warning("Execution generated empty result set. Omitting output operations.")
        return
    ext = destination.suffix.lower()
    logger.info("Writing output data to destination path: %s", destination)

    if ext in (".gb", ".genbank"):
        logger.debug("Formatting output as GenBank sequence data records.")
        records_to_write = []
        for row in results:
            if row.status == "SUCCESS" and row.optimized_record is not None:
                rec = row.optimized_record
                rec.id = row.record_id
                rec.name = re.sub(r'[^a-zA-Z0-9_]', '_', row.record_id)[:16]
                records_to_write.append(rec)

        if records_to_write:
            with open(destination, "w") as f:
                SeqIO.write(records_to_write, f, "genbank")
            logger.info("Successfully recorded %d sequence structures to GenBank file.", len(records_to_write))
        else:
            logger.warning("No successful optimized records available to serialize to GenBank format.")

    elif ext in (".xlsx", ".xls"):
        logger.debug("Formatting output matrix as Excel file dataframe tables.")
        df = pd.DataFrame([row.model_dump() for row in results])
        df.to_excel(destination, index=False, engine="openpyxl")
        logger.info("Successfully wrote dataset matrix to Excel format.")
    elif ext == ".csv":
        logger.debug("Formatting output matrix as CSV tabular flat-file.")
        df = pd.DataFrame([row.model_dump() for row in results])
        df.to_csv(destination, index=False)
        logger.info("Successfully wrote dataset matrix to CSV format.")
    else:
        logger.error("Failed processing serialization: Unsupported destination file format standard '%s'", ext)
        raise ValueError(f"Unsupported output format: {ext}")