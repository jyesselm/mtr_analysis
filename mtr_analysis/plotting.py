"""
Plotting utilities for mutation kinetics visualization.

This module provides functions for creating publication-quality plots
of mutation fraction time series and curve fits.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mtr_analysis.fitting import FitResult, monoexponential


def create_kinetics_plot(
    times: NDArray[np.floating],
    values: NDArray[np.floating],
    fit_result: FitResult,
    mutation_name: str,
) -> Figure:
    """
    Create a kinetics plot with data and fitted curve.

    Args:
        times: Time values (x-axis).
        values: Mutation fraction values (y-axis).
        fit_result: Fitted parameters.
        mutation_name: Name for the plot title.

    Returns:
        Matplotlib Figure object.
    """
    fig, ax = plt.subplots()
    _plot_data_points(ax, times, values)
    _plot_fit_curve(ax, times, fit_result)
    _configure_axes(ax, mutation_name, fit_result)
    return fig


def _plot_data_points(ax: plt.Axes, times: NDArray, values: NDArray) -> None:
    """
    Plot experimental data points.

    Args:
        ax: Matplotlib axes.
        times: Time values.
        values: Observed values.
    """
    ax.scatter(times, values)


def _plot_fit_curve(ax: plt.Axes, times: NDArray, fit_result: FitResult) -> None:
    """
    Plot fitted monoexponential curve.

    Args:
        ax: Matplotlib axes.
        times: Original time values (used for range).
        fit_result: Fitted parameters.
    """
    t_fit = np.linspace(min(times), max(times), 1000)
    y_fit = monoexponential(t_fit, fit_result.y_max, fit_result.k)
    ax.plot(t_fit, y_fit, color="black")


def _configure_axes(ax: plt.Axes, mutation_name: str, fit_result: FitResult) -> None:
    """
    Configure axes labels, scale, and title.

    Args:
        ax: Matplotlib axes.
        mutation_name: Name for the title.
        fit_result: Fitted parameters for the title.
    """
    ax.set_xscale("symlog")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Mutational Fraction")
    title = (
        f"{mutation_name} "
        f"Ymax: {fit_result.y_max:.2f} "
        f"k: {fit_result.k:.5f} ± {fit_result.k_std:.5f}"
    )
    ax.set_title(title)


def save_plot(fig: Figure, output_path: Path, dpi: int = 200) -> None:
    """
    Save figure to file and close it.

    Args:
        fig: Matplotlib Figure to save.
        output_path: Destination file path.
        dpi: Resolution in dots per inch.
    """
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
