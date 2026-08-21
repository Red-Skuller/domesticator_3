import pytest
from pydantic import ValidationError
from domestica.schema import SequenceRecord, ResultRow, PipelineConfig


def test_sequence_record_defaults() -> None:
    record = SequenceRecord(record_id="seq_001")
    assert record.record_id == "seq_001"
    assert record.protein_sequence is None
    assert record.metadata == {}


def test_sequence_record_full() -> None:
    record = SequenceRecord(
        record_id="seq_002",
        protein_sequence="MKTAYIAK",
        metadata={"source": "uniprot"}
    )
    assert record.record_id == "seq_002"
    assert record.protein_sequence == "MKTAYIAK"
    assert record.metadata["source"] == "uniprot"


def test_result_row_defaults() -> None:
    row = ResultRow(record_id="seq_001", status="PENDING")
    assert row.record_id == "seq_001"
    assert row.protein_sequence is None
    assert row.naive_dna_sequence is None
    assert row.optimized_sequence is None
    assert row.vendor_score is None
    assert row.accepted is False
    assert row.status == "PENDING"
    assert row.error_message is None
    assert row.optimized_record is None


def test_result_row_optimized_record_exclusion() -> None:
    dummy_seq_record = object()
    row = ResultRow(
        record_id="seq_001",
        status="SUCCESS",
        optimized_record=dummy_seq_record
    )

    dumped_data = row.model_dump()

    # Ensure optimized_record is excluded from standard serializations
    assert "optimized_record" not in dumped_data
    assert dumped_data["status"] == "SUCCESS"
    assert row.optimized_record is dummy_seq_record


def test_pipeline_config_defaults() -> None:
    config = PipelineConfig()
    assert config.vendor_target is None
    assert config.product == "eblocks"
    assert config.max_workers == 2
    assert config.min_length == 300


def test_pipeline_config_field_validation() -> None:
    # Valid boundary values
    config = PipelineConfig(max_workers=1, min_length=0)
    assert config.max_workers == 1
    assert config.min_length == 0

    config_max = PipelineConfig(max_workers=64)
    assert config_max.max_workers == 64

    # Invalid max_workers (below minimum)
    with pytest.raises(ValidationError):
        PipelineConfig(max_workers=0)

    # Invalid max_workers (above maximum)
    with pytest.raises(ValidationError):
        PipelineConfig(max_workers=65)

    # Invalid min_length (negative value)
    with pytest.raises(ValidationError):
        PipelineConfig(min_length=-1)


def test_pipeline_config_environment_variables(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("DOMESTICA_VENDOR_TARGET", "idt")
    monkeypatch.setenv("DOMESTICA_PRODUCT", "gblocks")
    monkeypatch.setenv("DOMESTICA_MAX_WORKERS", "16")
    monkeypatch.setenv("DOMESTICA_MIN_LENGTH", "450")

    config = PipelineConfig()
    assert config.vendor_target == "idt"
    assert config.product == "gblocks"
    assert config.max_workers == 16
    assert config.min_length == 450