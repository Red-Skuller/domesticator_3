import pytest
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

@pytest.fixture(autouse=True)
def mock_env_vars(monkeypatch: pytest.MonkeyPatch) -> None:
    """Set default environment variables for Pydantic settings."""
    monkeypatch.setenv("DOMESTICA_IDT_CLIENT_ID", "test_idt_id")
    monkeypatch.setenv("DOMESTICA_IDT_CLIENT_SECRET", "test_idt_secret")
    monkeypatch.setenv("DOMESTICA_IDT_USERNAME", "test_idt_user")
    monkeypatch.setenv("DOMESTICA_IDT_PASSWORD", "test_idt_pass")
    monkeypatch.setenv("DOMESTICA_THERMOFISHER_CLIENT_ID", "test_tf_id")
    monkeypatch.setenv("DOMESTICA_THERMOFISHER_CLIENT_SECRET", "test_tf_secret")

@pytest.fixture
def template_gb_longer_insert() -> Path:
    file_path = Path(__file__).parent / "data" / "gg_insert_longer.gb"
    return file_path

@pytest.fixture
def template_gb_shorter_insert() -> Path:
    file_path = Path(__file__).parent / "data" / "gg_insert_shorter.gb"
    return file_path

@pytest.fixture
def input_fasta() -> Path:
    file_path = Path(__file__).parent / "data" / "test.protein.fasta"
    return file_path

@pytest.fixture
def input_xlsx() -> Path:
    file_path = Path(__file__).parent / "data" / "input.xlsx"
    return file_path

@pytest.fixture
def input_pdb() -> Path:
    file_path = Path(__file__).parent / "data" / "actor_9x_het3_003_005_006_charge_styr.pdb"
    return file_path