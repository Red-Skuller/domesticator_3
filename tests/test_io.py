from pathlib import Path
from domestica.io import parse_input_file, write_output_file
from domestica.schema import ResultRow
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def test_parse_input_file_fasta(input_fasta: Path):
    """Test FASTA file parsing."""
    records = parse_input_file(input_fasta)
    assert len(records) == 2
    assert records[0].record_id == "test1"
    assert records[0].protein_sequence == "RYANKIVLER"
    assert records[1].record_id == "test2"
    assert records[1].protein_sequence == "HQTDQGSALESMAN"


def test_parse_input_file_xlsx(input_xlsx: Path):
    """Test Excel file parsing."""
    records = parse_input_file(input_xlsx)
    assert len(records) == 2
    assert records[0].record_id == "test1"
    assert records[0].protein_sequence == "RYANKIVLER"
    assert records[1].record_id == "test2"
    assert records[1].protein_sequence == "HQTDQGSALESMAN"


def test_write_output_file_csv(tmp_path: Path):
    """Test writing outputs to a CSV file."""
    out_path = tmp_path / "output.csv"
    results = [ResultRow(record_id="test1", status="SUCCESS", optimized_sequence="ATGC")]
    write_output_file(results, out_path)
    assert out_path.exists()
    content = out_path.read_text()
    assert "test1" in content
    assert "ATGC" in content


def test_write_output_file_genbank(tmp_path: Path):
    """Test writing outputs to GenBank files."""
    out_path = tmp_path / "output.gb"
    rec = SeqRecord(Seq("ATGC"), id="seq_1",annotations={"molecule_type": "DNA"})
    results = [ResultRow(record_id="seq_1", status="SUCCESS", optimized_record=rec)]
    write_output_file(results, out_path)

    expected_file = tmp_path / "output_seq_1.gb"
    assert expected_file.exists()
    content = expected_file.read_text()
    assert "LOCUS" in content