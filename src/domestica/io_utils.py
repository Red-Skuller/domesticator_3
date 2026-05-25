from pathlib import Path
import pandas as pd
from Bio import SeqIO
import re

_AA_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]")


def clean_seq(seq: str) -> str:
    """Return an uppercase protein sequence containing only the 20 standard AAs.

    Args:
        seq: The input sequence string.

    Returns:
        The cleaned, uppercase sequence.
    """
    if seq is None: return ""
    s = str(seq).upper()
    s = re.sub(r"\s+", "", s)
    s = s.replace("-", "").replace("*", "")
    return _AA_RE.sub("", s)


from typing import Callable, List, Dict

def read_fasta(path: Path, name_col: str, seq_col: str) -> List[Dict[str, str]]:
    """Reads a FASTA file and returns a list of sequence records.

    Args:
        path: Path to the FASTA file.
        name_col: Unused for FASTA but kept for signature consistency.
        seq_col: Unused for FASTA but kept for signature consistency.

    Returns:
        List of dictionaries with 'id' and 'sequence'.
    """
    records = []
    for rec in SeqIO.parse(path, "fasta"):
        records.append({"id": rec.id, "sequence": clean_seq(str(rec.seq))})
    return records


def read_xlsx(path: Path, name_col: str, seq_col: str) -> List[Dict[str, str]]:
    """Reads an Excel file and returns a list of sequence records.

    Args:
        path: Path to the Excel file.
        name_col: Column header for names.
        seq_col: Column header for sequences.

    Returns:
        List of dictionaries with 'id' and 'sequence'.
    """
    df = pd.read_excel(path)
    # Assuming the standard columns are Name and Sequence or fallback to indexes
    actual_name_col = name_col if name_col in df.columns else df.columns[0]
    actual_seq_col = seq_col if seq_col in df.columns else df.columns[1]

    records = []
    for _, row in df.iterrows():
        records.append({
            "id": str(row[actual_name_col]),
            "sequence": clean_seq(str(row[actual_seq_col]))
        })
    return records

READERS_REGISTRY: Dict[str, Callable[[Path, str, str], List[Dict[str, str]]]] = {
    ".fasta": read_fasta,
    ".fa": read_fasta,
    ".xlsx": read_xlsx,
}


def register_reader(ext: str, func: Callable[[Path, str, str], List[Dict[str, str]]]):
    """Register a new file reader for a given extension.

    Args:
        ext: File extension including the dot (e.g., '.fasta').
        func: Reader function.
    """
    READERS_REGISTRY[ext.lower()] = func


def read_input(path: Path, name_col: str = "Name", seq_col: str = "Sequence") -> List[Dict[str, str]]:
    """Reads input using the appropriate reader from the registry.

    Args:
        path: Path to the input file.
        name_col: Column header for names (for Excel).
        seq_col: Column header for sequences (for Excel).

    Returns:
        List of dictionaries with 'id' and 'sequence'.
    """
    ext = path.suffix.lower()
    reader = READERS_REGISTRY.get(ext)
    if not reader:
        raise ValueError(f"Unsupported file extension: {ext}")
    return reader(path, name_col, seq_col)