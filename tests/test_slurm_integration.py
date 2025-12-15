"""Integration tests for SLURM modules to improve coverage."""

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from click.testing import CliRunner

from mtr_analysis.cli_slurm import (
    _generate_mutation_scripts,
    _generate_rna_map_scripts,
    _handle_job_submission,
    _print_pipeline_instructions,
    _submit_pipeline,
    slurm,
)
from mtr_analysis.slurm.job import SlurmConfig
from mtr_analysis.slurm.runner import SlurmRunner


@pytest.fixture
def cli_runner() -> CliRunner:
    """Create a Click test runner."""
    return CliRunner()


@pytest.fixture
def temp_dir() -> Path:
    """Create temporary directory."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


class TestGenerateRnaMapScripts:
    """Tests for _generate_rna_map_scripts function."""

    def test_generates_scripts(self, temp_dir: Path) -> None:
        """Test generating RNA-MaP scripts."""
        df = pd.DataFrame(
            {
                "construct": ["sample_t0", "sample_t15min"],
                "barcode_seq": ["ACGT", "TTTT"],
            }
        )
        config = SlurmConfig()
        scripts = _generate_rna_map_scripts(df, "data", str(temp_dir), config)
        assert len(scripts) == 2


class TestGenerateMutationScripts:
    """Tests for _generate_mutation_scripts function."""

    def test_generates_scripts(self, temp_dir: Path) -> None:
        """Test generating mutation analysis scripts."""
        dirs = [str(temp_dir / "dir1"), str(temp_dir / "dir2")]
        config = SlurmConfig()
        scripts = _generate_mutation_scripts(dirs, str(temp_dir), config)
        assert len(scripts) == 2


class TestHandleJobSubmission:
    """Tests for _handle_job_submission function."""

    def test_no_submit(self, temp_dir: Path, capsys) -> None:
        """Test without submitting jobs."""
        script_path = temp_dir / "test.sh"
        script_path.write_text("#!/bin/bash\necho test")

        _handle_job_submission(
            [script_path], str(temp_dir), dry_run=False, submit=False
        )
        captured = capsys.readouterr()
        assert "Use --submit" in captured.out

    def test_dry_run(self, temp_dir: Path, capsys) -> None:
        """Test with dry run."""
        script_path = temp_dir / "test.sh"
        script_path.write_text("#!/bin/bash\necho test")

        _handle_job_submission([script_path], str(temp_dir), dry_run=True, submit=False)
        captured = capsys.readouterr()
        assert "DRY RUN" in captured.out

    def test_submit(self, temp_dir: Path) -> None:
        """Test with actual submission (mocked)."""
        script_path = temp_dir / "test.sh"
        script_path.write_text("#!/bin/bash\necho test")

        with patch.object(SlurmRunner, "submit_job", return_value="12345"):
            _handle_job_submission(
                [script_path], str(temp_dir), dry_run=False, submit=True
            )


class TestPrintPipelineInstructions:
    """Tests for _print_pipeline_instructions function."""

    def test_prints_stages(self, capsys) -> None:
        """Test printing pipeline instructions."""
        scripts = {
            "rna_map": ["script1.sh", "script2.sh"],
            "mutations": ["mut1.sh"],
            "aggregation": ["agg.sh"],
            "fitting": ["fit.sh"],
        }
        _print_pipeline_instructions("slurm_scripts", scripts)
        captured = capsys.readouterr()
        assert "Pipeline stages" in captured.out
        assert "RNA-MaP jobs: 2" in captured.out


class TestSubmitPipeline:
    """Tests for _submit_pipeline function."""

    def test_submit_pipeline_dry_run(self, temp_dir: Path) -> None:
        """Test submitting pipeline in dry run mode."""
        # Create mock script files
        for name in ["rna1.sh", "rna2.sh", "mut1.sh", "agg.sh", "fit.sh"]:
            (temp_dir / name).write_text("#!/bin/bash\necho test")

        scripts = {
            "rna_map": [temp_dir / "rna1.sh", temp_dir / "rna2.sh"],
            "mutations": [temp_dir / "mut1.sh"],
            "aggregation": [temp_dir / "agg.sh"],
            "fitting": [temp_dir / "fit.sh"],
        }
        _submit_pipeline(scripts, str(temp_dir), dry_run=True)


class TestSlurmRunnerExtended:
    """Extended tests for SlurmRunner."""

    def test_print_summary_empty(self, temp_dir: Path, capsys) -> None:
        """Test print_summary with no jobs."""
        runner = SlurmRunner(script_dir=temp_dir, dry_run=True)
        runner.print_summary()
        captured = capsys.readouterr()
        assert "No jobs submitted" in captured.out

    def test_print_summary_with_jobs(self, temp_dir: Path, capsys) -> None:
        """Test print_summary with jobs."""
        runner = SlurmRunner(script_dir=temp_dir, dry_run=True)
        script = temp_dir / "test.sh"
        script.write_text("#!/bin/bash")
        runner.submit_job(script)
        runner.print_summary()
        captured = capsys.readouterr()
        assert "Submitted 1 jobs" in captured.out


class TestSetupMutationsWithData:
    """Tests for setup-mutations with data directories."""

    def test_with_data_dirs(self, cli_runner: CliRunner) -> None:
        """Test setup-mutations with existing data directories."""
        with cli_runner.isolated_filesystem():
            Path("data").mkdir()
            Path("data/sample1").mkdir()
            Path("data/sample2").mkdir()

            result = cli_runner.invoke(
                slurm, ["setup-mutations", "--script-dir", "scripts"]
            )
            assert result.exit_code == 0
            assert "Generated" in result.output


class TestSetupAggregationWithData:
    """Tests for setup-aggregation command."""

    def test_generates_script(self, cli_runner: CliRunner) -> None:
        """Test that aggregation script is generated."""
        with cli_runner.isolated_filesystem():
            Path("data").mkdir()
            Path("data/sample1").mkdir()

            result = cli_runner.invoke(
                slurm,
                ["setup-aggregation", "--script-dir", "scripts"],
            )
            assert result.exit_code == 0
            assert "Generated" in result.output
            assert Path("scripts/aggregate_mutations.sh").exists()
