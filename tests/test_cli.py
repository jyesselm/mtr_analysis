"""Tests for CLI module."""

from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from mtr_analysis.cli import (
    DEFAULT_SEQUENCE,
    _create_mutation_dataframe,
    _load_and_filter_data,
    _load_and_filter_fractions,
    cli,
)
from mtr_analysis.mutations import MutationResult


@pytest.fixture
def cli_runner() -> CliRunner:
    """Create a Click test runner."""
    return CliRunner()


class TestCli:
    """Tests for main CLI group."""

    def test_cli_help(self, cli_runner: CliRunner) -> None:
        """Test CLI help command."""
        result = cli_runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "MTR ribozyme analysis toolkit" in result.output


class TestLoadAndFilterData:
    """Tests for _load_and_filter_data function."""

    def test_filters_by_code(self, temp_dir: Path) -> None:
        """Test that data is filtered by code."""
        csv_path = temp_dir / "test_data.csv"
        df = pd.DataFrame(
            {
                "code": ["C01HP", "OTHER", "C01HP"],
                "construct": ["mtr1_t0", "other_t0", "mtr1_t15min"],
                "barcode_seq": ["ACGT", "GGGG", "TTTT"],
            }
        )
        df.to_csv(csv_path, index=False)

        result = _load_and_filter_data(str(csv_path))

        assert len(result) == 2
        assert all(result["code"] == "C01HP")

    def test_adds_time_column(self, temp_dir: Path) -> None:
        """Test that time column is added."""
        csv_path = temp_dir / "test_data.csv"
        df = pd.DataFrame(
            {
                "code": ["C01HP"],
                "construct": ["mtr1_t30min"],
                "barcode_seq": ["ACGT"],
            }
        )
        df.to_csv(csv_path, index=False)

        result = _load_and_filter_data(str(csv_path))

        assert "time" in result.columns
        assert result["time"].iloc[0] == 30


class TestCreateMutationDataframe:
    """Tests for _create_mutation_dataframe function."""

    def test_creates_correct_columns(self) -> None:
        """Test that DataFrame has correct columns."""
        result = MutationResult(
            mutation_counts={"WT": 50, "A42G": 10},
            info_counts={"WT": 100, "A42G": 50},
            position_histogram=[0] * 100,
            total_reads=150,
        )

        df = _create_mutation_dataframe(result)

        assert list(df.columns) == ["mut", "mut_count", "info_count", "mut_fraction"]

    def test_computes_fractions(self) -> None:
        """Test that mutation fractions are computed correctly."""
        result = MutationResult(
            mutation_counts={"WT": 50},
            info_counts={"WT": 100},
            position_histogram=[0] * 100,
            total_reads=100,
        )

        df = _create_mutation_dataframe(result)

        assert df[df["mut"] == "WT"]["mut_fraction"].iloc[0] == 0.5

    def test_sorted_by_fraction(self) -> None:
        """Test that results are sorted by mutation fraction."""
        result = MutationResult(
            mutation_counts={"WT": 10, "A42G": 50},
            info_counts={"WT": 100, "A42G": 100},
            position_histogram=[0] * 100,
            total_reads=200,
        )

        df = _create_mutation_dataframe(result)

        # Higher fraction should be first
        assert df.iloc[0]["mut"] == "A42G"


class TestLoadAndFilterFractions:
    """Tests for _load_and_filter_fractions function."""

    def test_filters_by_info_count(self, temp_dir: Path) -> None:
        """Test filtering by minimum info count."""
        csv_path = temp_dir / "fractions.csv"
        df = pd.DataFrame(
            {
                "mut": ["WT", "A42G", "C55T"],
                "info_count": [1000, 500, 2000],
                "mut_fraction": [0.1, 0.2, 0.3],
                "time": [0, 0, 0],
            }
        )
        df.to_csv(csv_path, index=False)

        result = _load_and_filter_fractions(str(csv_path), min_count=1000)

        assert len(result) == 2
        assert "A42G" not in result["mut"].values

    def test_sorts_by_time(self, temp_dir: Path) -> None:
        """Test that results are sorted by time."""
        csv_path = temp_dir / "fractions.csv"
        df = pd.DataFrame(
            {
                "mut": ["WT", "WT", "WT"],
                "info_count": [1000, 1000, 1000],
                "mut_fraction": [0.1, 0.2, 0.3],
                "time": [60, 15, 30],
            }
        )
        df.to_csv(csv_path, index=False)

        result = _load_and_filter_fractions(str(csv_path), min_count=100)

        assert list(result["time"]) == [15, 30, 60]


class TestDefaultSequence:
    """Tests for DEFAULT_SEQUENCE constant."""

    def test_sequence_length(self) -> None:
        """Test that default sequence has expected length."""
        # The sequence should be long enough for analysis
        assert len(DEFAULT_SEQUENCE) > 100

    def test_sequence_valid_characters(self) -> None:
        """Test that sequence contains only valid nucleotides."""
        valid_chars = {"A", "C", "G", "T"}
        assert all(char in valid_chars for char in DEFAULT_SEQUENCE)


class TestCliCommands:
    """Tests for CLI commands (integration-style)."""

    def test_run_group_exists(self, cli_runner: CliRunner) -> None:
        """Test that run command group is available."""
        result = cli_runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        assert "Pipeline step commands" in result.output

    def test_run_single_rna_map_help(self, cli_runner: CliRunner) -> None:
        """Test run single-rna-map help."""
        result = cli_runner.invoke(cli, ["run", "single-rna-map", "--help"])
        assert result.exit_code == 0
        assert "BARCODE_SEQ" in result.output

    def test_get_mutation_fractions_no_data(self, cli_runner: CliRunner) -> None:
        """Test run get-mutations with no data directories."""
        with cli_runner.isolated_filesystem():
            result = cli_runner.invoke(cli, ["run", "get-mutations"])
            assert "No data directories found" in result.output

    def test_process_single_dir_help(self, cli_runner: CliRunner) -> None:
        """Test run process-dir help."""
        result = cli_runner.invoke(cli, ["run", "process-dir", "--help"])
        assert result.exit_code == 0
        assert "DIR_PATH" in result.output

    def test_fit_mut_fractions_help(self, cli_runner: CliRunner) -> None:
        """Test run fit help."""
        result = cli_runner.invoke(cli, ["run", "fit", "--help"])
        assert result.exit_code == 0
        assert "--min-info-count" in result.output
        assert "--plot" in result.output

    def test_slurm_group_exists(self, cli_runner: CliRunner) -> None:
        """Test that slurm command group is available."""
        result = cli_runner.invoke(cli, ["slurm", "--help"])
        assert result.exit_code == 0
        assert "SLURM job management" in result.output
