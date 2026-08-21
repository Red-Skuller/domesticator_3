import pytest
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from domestica.optimizer import _get_codon_table, _insert_into_template, optimize_sequence


@pytest.mark.parametrize(
    "template_fixture",
    ["template_gb_longer_insert", "template_gb_shorter_insert"]
)
def test_get_codon_table_parameterized(template_fixture: str, request: pytest.FixtureRequest) -> None:
    """Verify codon table extraction for both longer and shorter insert templates."""
    template_path = request.getfixturevalue(template_fixture)
    records = list(SeqIO.parse(template_path, "genbank"))
    record = records[0]
    table_name = _get_codon_table(record)
    assert isinstance(table_name, str)
    assert len(table_name) > 0


@pytest.mark.parametrize(
    "template_fixture",
    ["template_gb_longer_insert", "template_gb_shorter_insert"]
)
def test_insert_into_template_and_annotations(template_fixture: str, request: pytest.FixtureRequest) -> None:
    """Verify sequence insertion, feature mapping, and annotation transfer for both templates."""
    template_path = request.getfixturevalue(template_fixture)
    records = list(SeqIO.parse(template_path, "genbank"))
    record = records[0]
    insert_dna = "ATGGCTCGTAAATAA" * 20

    merged_record = _insert_into_template(record, insert_dna)

    assert isinstance(merged_record, SeqRecord)
    assert insert_dna in str(merged_record.seq)

    # Verify that all top-level annotations from the template are transferred and preserved
    for key, value in record.annotations.items():
        assert key in merged_record.annotations
        assert merged_record.annotations[key] == value

    # Verify features are transferred and re-attached correctly
    assert len(merged_record.features) == len(record.features)


@pytest.mark.parametrize(
    "template_fixture",
    ["template_gb_longer_insert", "template_gb_shorter_insert"]
)
def test_optimize_sequence_parameterized(template_fixture: str, request: pytest.FixtureRequest) -> None:
    """Verify full sequence optimization execution with both template types."""
    template_path = request.getfixturevalue(template_fixture)
    protein_seq = "MAR"
    naive_dna, optimized_dna, is_accepted, score, merged_record = optimize_sequence(
        protein_sequence=protein_seq,
        template_path=template_path,
        min_length=150
    )

    assert naive_dna is not None
    assert optimized_dna is not None
    assert isinstance(is_accepted, bool)
    assert isinstance(merged_record, SeqRecord)
    assert len(optimized_dna) >= 150


def test_optimize_sequence_with_evaluator(template_gb_longer_insert: Path) -> None:
    """Verify custom evaluator function integration correctly sets acceptance and score flags."""
    def mock_evaluator(sequence: str):
        return True, 4.5

    naive_dna, optimized_dna, is_accepted, score, merged_record = optimize_sequence(
        protein_sequence="MAR",
        template_path=template_gb_longer_insert,
        min_length=150,
        evaluator=mock_evaluator
    )

    assert is_accepted is True
    assert score == 4.5
    assert optimized_dna is not None