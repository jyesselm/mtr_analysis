"""Integration tests for CLI module to improve coverage."""

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from mtr_analysis.cli import (
    _fit_all_mutations,
    _fit_single_mutation,
    _generate_mutation_plot,
    _process_directory,
    _save_combined_results,
    _save_kinetics_results,
    cli,
)
from mtr_analysis.fitting import FitResult


@pytest.fixture
def cli_runner() -> CliRunner:
    """Create a Click test runner."""
    return CliRunner()


class TestFitAllMutations:
    """Tests for _fit_all_mutations function."""

    def test_skips_small_groups(self) -> None:
        """Test that groups with < 5 points are skipped."""
        df = pd.DataFrame(
            {
                "mut": ["A1G"] * 4,  # Only 4 points
                "time": [0, 15, 30, 60],
                "mut_fraction": [0.0, 0.1, 0.2, 0.3],
            }
        )
        results = _fit_all_mutations(df, generate_plots=False)
        assert len(results) == 0

    def test_fits_valid_group(self) -> None:
        """Test that valid groups are fitted."""
        times = [0, 15, 30, 60, 120, 240]
        values = [0.0, 0.15, 0.25, 0.4, 0.55, 0.65]
        df = pd.DataFrame(
            {
                "mut": ["A1G"] * 6,
                "time": times,
                "mut_fraction": values,
            }
        )
        results = _fit_all_mutations(df, generate_plots=False)
        assert len(results) == 1
        assert results[0][0] == "A1G"


class TestFitSingleMutation:
    """Tests for _fit_single_mutation function."""

    def test_fit_success(self) -> None:
        """Test successful fitting."""
        times = [0, 15, 30, 60, 120, 240]
        values = [0.0, 0.15, 0.25, 0.4, 0.55, 0.65]
        df = pd.DataFrame({"time": times, "mut_fraction": values})
        result = _fit_single_mutation("A1G", df, generate_plots=False)
        assert result is not None
        assert result[0] == "A1G"

    def test_fit_with_plots(self) -> None:
        """Test fitting with plot generation."""
        times = [0, 15, 30, 60, 120, 240]
        values = [0.0, 0.15, 0.25, 0.4, 0.55, 0.65]
        df = pd.DataFrame({"time": times, "mut_fraction": values})
        with patch("mtr_analysis.cli._generate_mutation_plot") as mock_plot:
            result = _fit_single_mutation("A1G", df, generate_plots=True)
            assert result is not None
            mock_plot.assert_called_once()


class TestGenerateMutationPlot:
    """Tests for _generate_mutation_plot function."""

    def test_creates_plot(self, temp_dir: Path) -> None:
        """Test that plot is generated and saved."""
        times = np.array([0, 15, 30, 60, 120])
        values = np.array([0.0, 0.2, 0.4, 0.5, 0.6])
        fit_result = FitResult(y_max=0.7, k=0.02, k_std=0.005)

        plot_dir = temp_dir / "plots"
        plot_dir.mkdir()

        with (
            patch("mtr_analysis.cli.save_plot"),
            patch("mtr_analysis.cli.create_kinetics_plot") as mock_create,
        ):
            mock_create.return_value = "mock_fig"
            _generate_mutation_plot("A1G", times, values, fit_result)
            mock_create.assert_called_once()


class TestSaveCombinedResults:
    """Tests for _save_combined_results function."""

    def test_saves_csv(self, cli_runner: CliRunner) -> None:
        """Test that combined results are saved."""
        df1 = pd.DataFrame({"mut": ["A1G"], "count": [10]})
        df2 = pd.DataFrame({"mut": ["C2T"], "count": [20]})

        with cli_runner.isolated_filesystem():
            _save_combined_results([df1, df2])
            result_df = pd.read_csv("all_mut_fractions.csv")
            assert len(result_df) == 2


class TestSaveKineticsResults:
    """Tests for _save_kinetics_results function."""

    def test_saves_csv(self, cli_runner: CliRunner) -> None:
        """Test that kinetics results are saved."""
        results = [
            ["A1G", 0.8, 0.02, 0.005],
            ["C2T", 0.6, 0.01, 0.003],
        ]

        with cli_runner.isolated_filesystem():
            _save_kinetics_results(results, "kinetics.csv")
            df = pd.read_csv("kinetics.csv")
            assert len(df) == 2
            assert list(df.columns) == ["mut", "y_max", "k", "k_std"]


class TestAggregateMutations:
    """Tests for aggregate-mutations command."""

    def test_no_files(self, cli_runner: CliRunner) -> None:
        """Test with no mutation fraction files."""
        with cli_runner.isolated_filesystem():
            Path("data").mkdir()
            Path("data/sample1").mkdir()
            result = cli_runner.invoke(cli, ["aggregate-mutations"])
            assert "No mutation fraction files found" in result.output

    def test_aggregates_files(self, cli_runner: CliRunner) -> None:
        """Test aggregating multiple files."""
        with cli_runner.isolated_filesystem():
            Path("data").mkdir()
            Path("data/sample_t15min").mkdir()
            Path("data/sample_t30min").mkdir()

            df1 = pd.DataFrame({"mut": ["A1G"], "count": [10]})
            df2 = pd.DataFrame({"mut": ["C2T"], "count": [20]})
            df1.to_csv("data/sample_t15min/mut_fractions.csv", index=False)
            df2.to_csv("data/sample_t30min/mut_fractions.csv", index=False)

            result = cli_runner.invoke(cli, ["aggregate-mutations"])
            assert "Aggregated 2 files" in result.output


class TestFitMutFractions:
    """Tests for fit-mut-fractions command."""

    def test_file_not_found(self, cli_runner: CliRunner) -> None:
        """Test with missing input file."""
        with cli_runner.isolated_filesystem():
            result = cli_runner.invoke(
                cli, ["fit-mut-fractions", "--input-file", "nonexistent.csv"]
            )
            assert result.exit_code != 0

    def test_fit_with_data(self, cli_runner: CliRunner) -> None:
        """Test fitting with valid data."""
        with cli_runner.isolated_filesystem():
            # Create test data
            times = [0, 15, 30, 60, 120, 240]
            values = [0.0, 0.1, 0.2, 0.35, 0.5, 0.6]
            df = pd.DataFrame(
                {
                    "mut": ["A1G"] * 6,
                    "time": times,
                    "mut_fraction": values,
                    "info_count": [2000] * 6,
                }
            )
            df.to_csv("all_mut_fractions.csv", index=False)

            result = cli_runner.invoke(cli, ["fit-mut-fractions"])
            assert result.exit_code == 0 or "Saved kinetics results" in result.output


class TestProcessDirectory:
    """Tests for _process_directory function."""

    def test_process_directory(self, temp_dir: Path) -> None:
        """Test processing a directory with bitvector data."""
        # Create mock directory structure
        data_dir = temp_dir / "data" / "sample_t30min"
        output_dir = data_dir / "output" / "BitVector_Files"
        output_dir.mkdir(parents=True)

        # Create mock bitvector file
        sequence = "A" * 120
        bitvector_content = "header1\nheader2\nheader3\nread1\t" + "." * 120 + "\t0\n"
        (output_dir / "mtr1_mut_lib_wt_bitvectors.txt").write_text(bitvector_content)

        df = _process_directory(str(data_dir), sequence)
        assert "mut" in df.columns
        assert "time" in df.columns


@pytest.fixture
def temp_dir() -> Path:
    """Create temporary directory."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)
