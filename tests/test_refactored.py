import pytest
from pathlib import Path
from domestica.models import PipelineConfig, ResultRow
from domestica.core import ingest, analyze, export
import pandas as pd

def test_pipeline_config_validation():
    config = PipelineConfig(
        input_path=Path("input.fasta"),
        output_path=Path("output.xlsx"),
        params=["mw", "pi"]
    )
    assert config.input_path == Path("input.fasta")
    assert "mw" in config.params

def test_result_row_alias():
    row = ResultRow(
        Name="Test",
        AA_seq="MACT",
        MW=123.4
    )
    assert row.name == "Test"
    assert row.aa_seq == "MACT"

    data = row.model_dump(by_alias=True)
    assert data["Name"] == "Test"
    assert data["AA_seq"] == "MACT"

def test_analyze_logic():
    config = PipelineConfig(
        input_path=Path("dummy.fasta"),
        output_path=Path("dummy.xlsx"),
        params=["length"]
    )
    record = {"id": "Prot1", "sequence": "MACT"}
    row = analyze(record, config)
    assert row.name == "Prot1"
    assert row.length == 4

def test_ingest_fasta(tmp_path):
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(">Seq1\nMACT\n>Seq2\nMW")
    config = PipelineConfig(
        input_path=fasta_path,
        output_path=Path("out.xlsx")
    )
    records = ingest(config)
    assert len(records) == 2
    assert records[0]["id"] == "Seq1"
    assert records[0]["sequence"] == "MACT"

def test_export_xlsx(tmp_path):
    output_path = tmp_path / "output.xlsx"
    config = PipelineConfig(
        input_path=Path("in.fasta"),
        output_path=output_path,
        out_cols=["Name", "AA_seq", "length"]
    )
    results = [
        ResultRow(Name="P1", AA_seq="MACT", length=4),
        ResultRow(Name="P2", AA_seq="MW", length=2)
    ]
    export(results, config)

    assert output_path.exists()
    df = pd.read_excel(output_path)
    assert len(df) == 2
    assert list(df.columns) == ["Name", "AA_seq", "length"]
