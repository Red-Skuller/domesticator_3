from pathlib import Path
import pandas as pd
from Bio import SeqIO
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

_AA_RE = re.compile(r"[^ACDEFGHIKLMNPQRSTVWY]")


def clean_seq(seq: str) -> str:
    """Return an uppercase protein sequence containing only the 20 standard AAs."""
    if seq is None: return ""
    s = str(seq).upper()
    s = re.sub(r"\s+", "", s)
    s = s.replace("-", "").replace("*", "")
    return _AA_RE.sub("", s)


def read_input(path: Path, name_col="Name", seq_col="Sequence") -> list[dict]:
    """Reads FASTA or XLSX and returns a list of dicts: [{'id': ..., 'sequence': ...}]"""
    ext = path.suffix.lower()
    records = []

    if ext in [".fasta", ".fa"]:
        for rec in SeqIO.parse(path, "fasta"):
            records.append({"id": rec.id, "sequence": clean_seq(str(rec.seq))})

    elif ext in [".xlsx"]:
        df = pd.read_excel(path)
        # Assuming the standard columns are Name and Sequence or fallback to indexes
        name_col = name_col if name_col in df.columns else df.columns[0]
        seq_col = seq_col if seq_col in df.columns else df.columns[1]

        for _, row in df.iterrows():
            records.append({
                "id": str(row[name_col]),
                "sequence": clean_seq(str(row[seq_col]))
            })
    else:
        raise ValueError(f"Unsupported file extension: {ext}")

    return records


def write_output(results: list[dict], output_path: Path, out_cols: list[str] = None):
    """Write pipeline results to .xlsx, .fasta, or .gb based on file extension."""
    ext = output_path.suffix.lower()

    if ext == ".xlsx":
        # Remove the 'DNA_record' object before passing to pandas to avoid serialization errors
        clean_results = [{k: v for k, v in row.items() if k != 'DNA_record'} for row in results]
        df = pd.DataFrame(clean_results)
        if out_cols:
            df = df.reindex(columns=out_cols)
        df.to_excel(output_path, index=False)

    elif ext in [".fasta", ".fa", ".gb", ".gbk"]:
        records = []
        for row in results:
            safe_id = str(row.get("Name", "unknown"))[:16]

            # If generating a GenBank file AND we preserved the fully annotated vector record
            if ext in [".gb", ".gbk"] and "DNA_record" in row:
                rec = row["DNA_record"]
                rec.id = safe_id
                rec.name = safe_id
                rec.annotations["molecule_type"] = "DNA"
                records.append(rec)

            # Fallback for FASTA (or if optimization was skipped/failed)
            else:
                seq_str = row.get("DNA_seq") or row.get("AA_Seq_Final") or row.get("AA_seq")
                if not seq_str:
                    continue

                rec = SeqRecord(
                    Seq(seq_str),
                    id=safe_id,
                    name=safe_id,
                    description="Domestica Pipeline Output"
                )
                if ext in [".gb", ".gbk"]:
                    rec.annotations["molecule_type"] = "DNA" if "DNA_seq" in row else "protein"

                records.append(rec)

        fmt = "genbank" if ext in [".gb", ".gbk"] else "fasta"
        SeqIO.write(records, output_path, fmt)
    else:
        raise ValueError(f"Unsupported output extension: {ext}")