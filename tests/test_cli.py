# File: /home/lukah/Projects/domestica/tests/test_cli.py

import pytest
import logging
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock, call
from typer.testing import CliRunner
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from domestica.cli import app, _worker_task, _init_worker, _worker_evaluator
from domestica.schema import SequenceRecord, ResultRow, PipelineConfig
from domestica.optimizer import optimize_sequence

# Create a test runner
runner = CliRunner()


# Test fixtures and helpers
@pytest.fixture
def sample_sequence_records():
    """Create sample SequenceRecord objects."""
    return [
        SequenceRecord(
            record_id="test_seq_1",
            protein_sequence="MKTIIALSYIFCLVFADYKDDDDK",
            metadata={"description": "Test sequence 1"}
        ),
        SequenceRecord(
            record_id="test_seq_2",
            protein_sequence="MVSKGEELFTGVVPILVELDGDVNGHK",
            metadata={"description": "Test sequence 2"}
        ),
    ]


# Test basic CLI command execution
class TestCLIBasicExecution:

    def test_optimize_with_fasta_input(self, input_fasta, tmp_path):
        """Test optimization with FASTA input."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1",
                "--min-length", "300"
            ]
        )

        assert result.exit_code == 0, f"CLI failed with error: {result.output}"
        assert output_path.exists()

        # Check CSV output
        df = pd.read_csv(output_path)
        assert len(df) > 0
        assert 'record_id' in df.columns
        assert 'status' in df.columns
        assert all(row['status'] == 'SUCCESS' for _, row in df.iterrows())

    def test_optimize_with_xlsx_input(self, input_xlsx, tmp_path):
        """Test optimization with Excel input."""
        output_path = tmp_path / "output.xlsx"

        result = runner.invoke(
            app,
            [
                "--input", str(input_xlsx),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code == 0, f"CLI failed with error: {result.output}"
        assert output_path.exists()

        # Check Excel output
        df = pd.read_excel(output_path)
        assert len(df) > 0
        assert 'record_id' in df.columns
        assert 'status' in df.columns

    def test_optimize_without_input(self, tmp_path):
        """Test optimization without input (template-only mode)."""
        output_path = tmp_path / "output.csv"
        template_path = Path(__file__).parent / "data" / "gg_insert_longer.gb"

        result = runner.invoke(
            app,
            [
                "--output", str(output_path),
                "--template", str(template_path),
                "--workers", "1"
            ]
        )

        assert result.exit_code == 0, f"CLI failed with error: {result.output}"
        assert output_path.exists()

        # Verify output has the template record
        df = pd.read_csv(output_path)
        assert len(df) == 1
        assert df.iloc[0]['record_id'] == 'optimized_template'

    def test_optimize_with_verbose_flag(self, input_fasta, tmp_path):
        """Test optimization with verbose logging."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1",
                "--verbose"
            ]
        )

        assert result.exit_code == 0, f"CLI failed with error: {result.output}"
        assert output_path.exists()

    def test_optimize_with_vendor(self, input_fasta, tmp_path):
        """Test optimization with vendor evaluator."""
        output_path = tmp_path / "output.csv"

        # Since we can't actually call vendor APIs in tests,
        # this should fail gracefully but still complete
        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--vendor", "idt",
                "--workers", "1"
            ]
        )

        # The process should still complete even if vendor is unavailable
        assert result.exit_code == 0
        assert output_path.exists()

    def test_optimize_with_short_flags(self, input_fasta, tmp_path):
        """Test optimization with short flag names."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "-i", str(input_fasta),
                "-o", str(output_path),
                "-t", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "-w", "1"
            ]
        )

        assert result.exit_code == 0, f"CLI failed with error: {result.output}"
        assert output_path.exists()


# Test different output formats
class TestOutputFormats:

    def test_genbank_output(self, input_fasta, tmp_path):
        """Test GenBank output format."""
        output_path = tmp_path / "output.gb"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code == 0
        # GenBank output creates multiple files with suffixes
        genbank_files = list(tmp_path.glob("output_*.gb"))
        assert len(genbank_files) > 0

        # Verify each GenBank file is valid
        for gb_file in genbank_files:
            records = list(SeqIO.parse(gb_file, "genbank"))
            assert len(records) == 1

    def test_csv_output(self, input_fasta, tmp_path):
        """Test CSV output format."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code == 0
        df = pd.read_csv(output_path)
        assert 'record_id' in df.columns
        assert 'protein_sequence' in df.columns
        assert 'optimized_sequence' in df.columns
        assert 'vendor_score' in df.columns
        assert 'accepted' in df.columns
        assert 'status' in df.columns

    def test_excel_output(self, input_fasta, tmp_path):
        """Test Excel output format."""
        output_path = tmp_path / "output.xlsx"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code == 0
        df = pd.read_excel(output_path)
        assert 'record_id' in df.columns
        assert 'status' in df.columns


# Test error handling and edge cases
class TestErrorHandling:

    def test_invalid_input_format(self, input_pdb, tmp_path):
        """Test with invalid input format (PDB)."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_pdb),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code != 0

    def test_nonexistent_template(self, input_fasta, tmp_path):
        """Test with nonexistent template file."""
        output_path = tmp_path / "output.csv"
        nonexistent_template = tmp_path / "nonexistent.gb"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(nonexistent_template),
                "--workers", "1"
            ]
        )

        assert result.exit_code != 0

    def test_invalid_output_format(self, input_fasta, tmp_path):
        """Test with unsupported output format."""
        output_path = tmp_path / "output.xyz"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code != 0

    def test_invalid_workers_count(self, input_fasta, tmp_path):
        """Test with invalid workers count."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "100"  # Invalid: exceeds 64
            ]
        )

        assert result.exit_code != 0

    def test_invalid_min_length(self, input_fasta, tmp_path):
        """Test with invalid min_length (negative)."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--min-length", "-100",  # Invalid
                "--workers", "1"
            ]
        )

        assert result.exit_code != 0

    def test_missing_template_option(self, input_fasta, tmp_path):
        """Test without required template option."""
        output_path = tmp_path / "output.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--workers", "1"
            ]
        )

        assert result.exit_code != 0

    def test_missing_output_option(self, input_fasta, tmp_path):
        """Test without required output option."""
        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code != 0

    def test_nonexistent_input_file(self, tmp_path):
        """Test with nonexistent input file."""
        output_path = tmp_path / "output.csv"
        nonexistent_input = tmp_path / "nonexistent.fasta"

        result = runner.invoke(
            app,
            [
                "--input", str(nonexistent_input),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code != 0


# Test worker functions
class TestWorkerFunctions:

    def test_worker_task_success(self, template_gb_longer_insert, sample_sequence_records):
        """Test worker task with successful optimization."""
        record = sample_sequence_records[0]

        # Mock the optimize_sequence function
        mock_record = MagicMock(spec=SeqRecord)
        with patch('domestica.cli.optimize_sequence') as mock_opt:
            mock_opt.return_value = (
                "ATGCTAGCTA",  # naive
                "ATGCTAGCTAGCTA",  # optimized
                True,  # accepted
                9.5,  # score
                mock_record  # record
            )

            result = _worker_task(record, template_gb_longer_insert, 300)

            assert isinstance(result, ResultRow)
            assert result.record_id == record.record_id
            assert result.status == "SUCCESS"
            assert result.accepted == True
            assert result.vendor_score == 9.5
            assert result.optimized_sequence is not None
            assert result.optimized_record == mock_record

    def test_worker_task_failure(self, template_gb_longer_insert, sample_sequence_records):
        """Test worker task with failed optimization."""
        record = sample_sequence_records[1]

        # Mock the optimize_sequence function to raise an exception
        with patch('domestica.cli.optimize_sequence', side_effect=Exception("Test error")):
            result = _worker_task(record, template_gb_longer_insert, 300)

            assert isinstance(result, ResultRow)
            assert result.record_id == record.record_id
            assert result.status == "FAILED"
            assert result.error_message == "Test error"

    def test_worker_task_with_none_sequence(self, template_gb_longer_insert):
        """Test worker task with None protein sequence."""
        record = SequenceRecord(record_id="test_none", protein_sequence=None)

        with patch('domestica.cli.optimize_sequence') as mock_opt:
            mock_opt.return_value = (
                None,  # naive
                "ATGCTAGCTAGCTA",  # optimized
                True,  # accepted
                None,  # score
                MagicMock(spec=SeqRecord)  # record
            )

            result = _worker_task(record, template_gb_longer_insert, 300)

            assert isinstance(result, ResultRow)
            assert result.record_id == "test_none"
            assert result.protein_sequence is None
            assert result.status == "SUCCESS"

    def test_worker_init_with_vendor(self):
        """Test worker initialization with vendor."""
        import domestica.cli as cli_module

        with patch('domestica.cli.get_evaluator') as mock_get_eval:
            mock_evaluator = MagicMock()
            mock_get_eval.return_value = mock_evaluator

            cli_module._init_worker("idt", "eblocks", False)

            assert cli_module._worker_evaluator == mock_evaluator
            mock_get_eval.assert_called_once_with("idt", "eblocks")

    def test_worker_init_without_vendor(self):
        """Test worker initialization without vendor."""
        import domestica.cli as cli_module

        cli_module._init_worker(None, "eblocks", False)
        assert cli_module._worker_evaluator is None


# Test with different template files
class TestTemplateVariations:

    def test_optimize_with_longer_insert_template(self, input_fasta, tmp_path, template_gb_longer_insert):
        """Test optimization with longer insert template."""
        output_path = tmp_path / "output_longer.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(template_gb_longer_insert),
                "--workers", "1"
            ]
        )

        assert result.exit_code == 0
        assert output_path.exists()

        # Verify results
        df = pd.read_csv(output_path)
        assert len(df) > 0
        assert all(row['status'] == 'SUCCESS' for _, row in df.iterrows())

    def test_optimize_with_shorter_insert_template(self, input_fasta, tmp_path, template_gb_shorter_insert):
        """Test optimization with shorter insert template."""
        output_path = tmp_path / "output_shorter.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(template_gb_shorter_insert),
                "--workers", "1"
            ]
        )

        assert result.exit_code == 0
        assert output_path.exists()

        # Verify results (may have failures due to template constraints)
        df = pd.read_csv(output_path)
        assert len(df) > 0


# Test concurrent execution
class TestConcurrentExecution:

    def test_multiple_workers(self, input_fasta, tmp_path):
        """Test optimization with multiple workers."""
        output_path = tmp_path / "output_multi.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "2"
            ]
        )

        assert result.exit_code == 0, f"CLI failed with error: {result.output}"
        assert output_path.exists()

        # Verify output contains results for all input sequences
        df = pd.read_csv(output_path)
        assert len(df) > 0

    def test_single_worker(self, input_fasta, tmp_path):
        """Test optimization with single worker."""
        output_path = tmp_path / "output_single.csv"

        result = runner.invoke(
            app,
            [
                "--input", str(input_fasta),
                "--output", str(output_path),
                "--template", str(Path(__file__).parent / "data" / "gg_insert_longer.gb"),
                "--workers", "1"
            ]
        )

        assert result.exit_code == 0
        assert output_path.exists()


# Test configuration validation
class TestConfiguration:

    def test_pipeline_config_defaults(self):
        """Test PipelineConfig with default values."""
        config = PipelineConfig()
        assert config.max_workers == 2
        assert config.min_length == 300
        assert config.product == "eblocks"
        assert config.vendor_target is None

    def test_pipeline_config_custom(self):
        """Test PipelineConfig with custom values."""
        config = PipelineConfig(
            vendor_target="idt",
            product="gblocks",
            max_workers=4,
            min_length=500
        )
        assert config.vendor_target == "idt"
        assert config.product == "gblocks"
        assert config.max_workers == 4
        assert config.min_length == 500

    def test_pipeline_config_validation(self):
        """Test PipelineConfig validation constraints."""
        with pytest.raises(Exception):
            PipelineConfig(max_workers=100)  # Exceeds maximum of 64

        with pytest.raises(Exception):
            PipelineConfig(max_workers=0)  # Below minimum of 1

        with pytest.raises(Exception):
            PipelineConfig(min_length=-1)  # Below minimum of 0

    def test_pipeline_config_env_vars(self):
        """Test PipelineConfig with environment variables."""
        with patch.dict('os.environ', {
            'DOMESTICA_MAX_WORKERS': '8',
            'DOMESTICA_MIN_LENGTH': '500',
            'DOMESTICA_PRODUCT': 'gblocks',
            'DOMESTICA_VENDOR_TARGET': 'idt'
        }, clear=False):
            config = PipelineConfig()
            assert config.max_workers == 8
            assert config.min_length == 500
            assert config.product == 'gblocks'
            assert config.vendor_target == 'idt'