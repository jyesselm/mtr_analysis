"""Tests for plotting module."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mtr_analysis.fitting import FitResult
from mtr_analysis.plotting import (
    _configure_axes,
    _plot_data_points,
    _plot_fit_curve,
    create_kinetics_plot,
    save_plot,
)


class TestCreateKineticsPlot:
    """Tests for create_kinetics_plot function."""

    def test_returns_figure(self) -> None:
        """Test that function returns a matplotlib Figure."""
        times = np.array([0, 15, 30, 60, 120])
        values = np.array([0.0, 0.2, 0.4, 0.6, 0.7])
        fit_result = FitResult(y_max=0.8, k=0.02, k_std=0.005)

        fig = create_kinetics_plot(times, values, fit_result, "TestMutation")

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_figure_has_axes(self) -> None:
        """Test that figure contains axes."""
        times = np.array([0, 15, 30, 60])
        values = np.array([0.0, 0.2, 0.4, 0.5])
        fit_result = FitResult(y_max=0.6, k=0.03, k_std=0.01)

        fig = create_kinetics_plot(times, values, fit_result, "Test")

        assert len(fig.axes) == 1
        plt.close(fig)


class TestSavePlot:
    """Tests for save_plot function."""

    def test_saves_file(self, temp_dir: Path) -> None:
        """Test that plot is saved to file."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        output_path = temp_dir / "test_plot.png"

        save_plot(fig, output_path, dpi=100)

        assert output_path.exists()
        assert output_path.stat().st_size > 0

    def test_closes_figure(self, temp_dir: Path) -> None:
        """Test that figure is closed after saving."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        output_path = temp_dir / "test_plot.png"

        save_plot(fig, output_path)

        # Figure should be closed (not in list of figures)
        assert fig not in plt.get_fignums()


class TestPlotDataPoints:
    """Tests for _plot_data_points helper function."""

    def test_scatter_plot_created(self) -> None:
        """Test that scatter plot is created."""
        fig, ax = plt.subplots()
        times = np.array([0, 10, 20])
        values = np.array([0.1, 0.2, 0.3])

        _plot_data_points(ax, times, values)

        # Check that children were added to axes
        assert len(ax.collections) > 0
        plt.close(fig)


class TestPlotFitCurve:
    """Tests for _plot_fit_curve helper function."""

    def test_line_plot_created(self) -> None:
        """Test that line plot is created."""
        fig, ax = plt.subplots()
        times = np.array([0, 10, 20, 30])
        fit_result = FitResult(y_max=0.5, k=0.05, k_std=0.01)

        _plot_fit_curve(ax, times, fit_result)

        # Check that line was added
        assert len(ax.lines) > 0
        plt.close(fig)


class TestConfigureAxes:
    """Tests for _configure_axes helper function."""

    def test_labels_set(self) -> None:
        """Test that axis labels are set correctly."""
        fig, ax = plt.subplots()
        fit_result = FitResult(y_max=0.8, k=0.02, k_std=0.005)

        _configure_axes(ax, "TestMut", fit_result)

        assert ax.get_xlabel() == "Time (min)"
        assert ax.get_ylabel() == "Mutational Fraction"
        plt.close(fig)

    def test_title_contains_mutation_name(self) -> None:
        """Test that title contains mutation name."""
        fig, ax = plt.subplots()
        fit_result = FitResult(y_max=0.8, k=0.02, k_std=0.005)

        _configure_axes(ax, "A42G", fit_result)

        assert "A42G" in ax.get_title()
        plt.close(fig)

    def test_title_contains_parameters(self) -> None:
        """Test that title contains fit parameters."""
        fig, ax = plt.subplots()
        fit_result = FitResult(y_max=0.8, k=0.02, k_std=0.005)

        _configure_axes(ax, "Test", fit_result)

        title = ax.get_title()
        assert "Ymax" in title
        assert "k:" in title
        plt.close(fig)
