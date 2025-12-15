"""Tests for CLI SLURM module."""

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from mtr_analysis.cli import cli
from mtr_analysis.cli_slurm import (
    _build_slurm_config,
    _load_constructs,
    slurm,
)


@pytest.fixture
def cli_runner() -> CliRunner:
    """Create a Click test runner."""
    return CliRunner()


class TestSlurmCliGroup:
    """Tests for slurm CLI command group."""

    def test_slurm_help(self, cli_runner: CliRunner) -> None:
        """Test slurm help command."""
        result = cli_runner.invoke(slurm, ["--help"])
        assert result.exit_code == 0
        assert "SLURM job management" in result.output

    def test_slurm_subcommands(self, cli_runner: CliRunner) -> None:
        """Test that slurm subcommands are available."""
        result = cli_runner.invoke(slurm, ["--help"])
        assert "setup-rna-map" in result.output
        assert "setup-mutations" in result.output
        assert "setup-aggregation" in result.output
        assert "setup-fitting" in result.output
        assert "setup-full-pipeline" in result.output


class TestBuildSlurmConfig:
    """Tests for _build_slurm_config function."""

    def test_basic_config(self) -> None:
        """Test building basic config."""
        config = _build_slurm_config(
            time="02:00:00",
            memory="8G",
        )

        assert config.time == "02:00:00"
        assert config.memory == "8G"

    def test_config_defaults(self) -> None:
        """Test building config with defaults."""
        config = _build_slurm_config(
            time="01:00:00",
            memory="4G",
        )

        assert config.time == "01:00:00"
        assert config.memory == "4G"
        assert config.cpus == 1


class TestLoadConstructs:
    """Tests for _load_constructs function."""

    def test_loads_and_filters(self, temp_dir: Path) -> None:
        """Test loading and filtering constructs."""
        csv_path = temp_dir / "constructs.csv"
        df = pd.DataFrame(
            {
                "code": ["C01HP", "OTHER", "C01HP"],
                "construct": ["mtr1_t0", "other_t0", "mtr1_t30min"],
                "barcode_seq": ["ACGT", "GGGG", "TTTT"],
            }
        )
        df.to_csv(csv_path, index=False)

        result = _load_constructs(str(csv_path))

        assert len(result) == 2
        assert "time" in result.columns


class TestSetupRnaMapCommand:
    """Tests for setup-rna-map command."""

    def test_help(self, cli_runner: CliRunner) -> None:
        """Test setup-rna-map help."""
        result = cli_runner.invoke(slurm, ["setup-rna-map", "--help"])
        assert result.exit_code == 0
        assert "DATA_CSV" in result.output
        assert "--time" in result.output
        assert "--dry-run" in result.output

    def test_generates_scripts(self, cli_runner: CliRunner, temp_dir: Path) -> None:
        """Test that scripts are generated."""
        # Create test data CSV
        csv_path = temp_dir / "data.csv"
        df = pd.DataFrame(
            {
                "code": ["C01HP"],
                "construct": ["mtr1_t0"],
                "barcode_seq": ["ACGT"],
            }
        )
        df.to_csv(csv_path, index=False)

        script_dir = temp_dir / "scripts"

        result = cli_runner.invoke(
            slurm,
            [
                "setup-rna-map",
                str(csv_path),
                "--script-dir",
                str(script_dir),
            ],
        )

        assert result.exit_code == 0
        assert "Generated" in result.output


class TestSetupMutationsCommand:
    """Tests for setup-mutations command."""

    def test_help(self, cli_runner: CliRunner) -> None:
        """Test setup-mutations help."""
        result = cli_runner.invoke(slurm, ["setup-mutations", "--help"])
        assert result.exit_code == 0
        assert "--time" in result.output
        assert "--submit" in result.output

    def test_no_data_dirs(self, cli_runner: CliRunner) -> None:
        """Test with no data directories."""
        with cli_runner.isolated_filesystem():
            result = cli_runner.invoke(slurm, ["setup-mutations"])
            assert "No data directories found" in result.output


class TestSetupAggregationCommand:
    """Tests for setup-aggregation command."""

    def test_help(self, cli_runner: CliRunner) -> None:
        """Test setup-aggregation help."""
        result = cli_runner.invoke(slurm, ["setup-aggregation", "--help"])
        assert result.exit_code == 0
        assert "--output" in result.output


class TestSetupFittingCommand:
    """Tests for setup-fitting command."""

    def test_help(self, cli_runner: CliRunner) -> None:
        """Test setup-fitting help."""
        result = cli_runner.invoke(slurm, ["setup-fitting", "--help"])
        assert result.exit_code == 0
        assert "--min-info-count" in result.output
        assert "--plot" in result.output

    def test_generates_script(self, cli_runner: CliRunner, temp_dir: Path) -> None:
        """Test that fitting script is generated."""
        script_dir = temp_dir / "scripts"

        result = cli_runner.invoke(
            slurm,
            [
                "setup-fitting",
                "--script-dir",
                str(script_dir),
                "--min-info-count",
                "500",
            ],
        )

        assert result.exit_code == 0
        assert "Generated" in result.output


class TestSetupFullPipelineCommand:
    """Tests for setup-full-pipeline command."""

    def test_help(self, cli_runner: CliRunner) -> None:
        """Test setup-full-pipeline help."""
        result = cli_runner.invoke(slurm, ["setup-full-pipeline", "--help"])
        assert result.exit_code == 0
        assert "DATA_CSV" in result.output
        assert "--output-dir" in result.output
        assert "--submit" in result.output

    def test_generates_pipeline(self, cli_runner: CliRunner, temp_dir: Path) -> None:
        """Test that full pipeline scripts are generated."""
        # Create test data CSV
        csv_path = temp_dir / "data.csv"
        df = pd.DataFrame(
            {
                "code": ["C01HP", "C01HP"],
                "construct": ["mtr1_t0", "mtr1_t15min"],
                "barcode_seq": ["ACGT", "TTTT"],
            }
        )
        df.to_csv(csv_path, index=False)

        script_dir = temp_dir / "scripts"

        result = cli_runner.invoke(
            slurm,
            [
                "setup-full-pipeline",
                str(csv_path),
                "--script-dir",
                str(script_dir),
            ],
        )

        assert result.exit_code == 0
        assert "Setting up full analysis pipeline" in result.output
        assert "Pipeline stages" in result.output


class TestSlurmCliIntegration:
    """Integration tests for SLURM CLI through main cli."""

    def test_slurm_accessible_from_main_cli(self, cli_runner: CliRunner) -> None:
        """Test that slurm commands are accessible from main CLI."""
        result = cli_runner.invoke(cli, ["slurm", "--help"])
        assert result.exit_code == 0
        assert "SLURM job management" in result.output

    def test_slurm_subcommand_from_main_cli(self, cli_runner: CliRunner) -> None:
        """Test accessing slurm subcommand from main CLI."""
        result = cli_runner.invoke(cli, ["slurm", "setup-fitting", "--help"])
        assert result.exit_code == 0
        assert "--min-info-count" in result.output
