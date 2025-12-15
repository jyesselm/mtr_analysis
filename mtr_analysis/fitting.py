"""
Curve fitting for kinetic data.

This module provides monoexponential fitting with bootstrap error estimation
for analyzing mutation fraction time series.
"""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import OptimizeWarning, curve_fit


@dataclass
class FitResult:
    """Results from monoexponential curve fitting."""

    y_max: float
    k: float
    k_std: float


def monoexponential(t: NDArray[np.floating], y_max: float, k: float) -> NDArray:
    """
    Monoexponential growth model.

    Args:
        t: Time values.
        y_max: Maximum asymptotic value.
        k: Rate constant.

    Returns:
        Model predictions at each time point.
    """
    return y_max * (1 - np.exp(-k * t))


def fit_monoexponential(
    times: NDArray[np.floating],
    values: NDArray[np.floating],
    n_bootstrap: int = 1000,
    seed: int = 42,
) -> FitResult:
    """
    Fit monoexponential model with bootstrap error estimation.

    Performs curve fitting to estimate parameters, then uses bootstrap
    resampling to estimate uncertainty in the rate constant.

    Args:
        times: Time values (independent variable).
        values: Observed values (dependent variable).
        n_bootstrap: Number of bootstrap iterations.
        seed: Random seed for reproducibility.

    Returns:
        FitResult containing fitted parameters and error estimate.
    """
    y_max, k = _fit_curve(times, values)
    k_std = _bootstrap_k_error(times, values, n_bootstrap, seed)
    return FitResult(y_max=y_max, k=k, k_std=k_std)


def _fit_curve(
    times: NDArray[np.floating], values: NDArray[np.floating]
) -> tuple[float, float]:
    """
    Perform single curve fit.

    Args:
        times: Time values.
        values: Observed values.

    Returns:
        Tuple of (y_max, k) parameters.
    """
    initial_guesses = [1.0, 0.1]
    popt, _ = curve_fit(
        monoexponential, times, values, p0=initial_guesses, maxfev=10000
    )
    return float(popt[0]), float(popt[1])


def _bootstrap_k_error(
    times: NDArray[np.floating],
    values: NDArray[np.floating],
    n_bootstrap: int,
    seed: int,
) -> float:
    """
    Estimate rate constant error using bootstrap resampling.

    Args:
        times: Time values.
        values: Observed values.
        n_bootstrap: Number of bootstrap samples.
        seed: Random seed.

    Returns:
        Standard deviation of rate constant estimates.
    """
    rng = np.random.default_rng(seed=seed)
    k_values = []
    initial_guesses = [1.0, 0.1]

    for _ in range(n_bootstrap):
        k_estimate = _single_bootstrap_iteration(times, values, rng, initial_guesses)
        if k_estimate is not None:
            k_values.append(k_estimate)

    if not k_values:
        return 0.0
    return float(np.std(k_values))


def _single_bootstrap_iteration(
    times: NDArray[np.floating],
    values: NDArray[np.floating],
    rng: np.random.Generator,
    initial_guesses: list,
) -> float | None:
    """
    Perform single bootstrap iteration.

    Args:
        times: Original time values.
        values: Original observed values.
        rng: Random number generator.
        initial_guesses: Initial parameter guesses.

    Returns:
        Fitted rate constant, or None if fit failed.
    """
    indices = rng.integers(0, len(times), len(times))
    times_sample = times[indices]
    values_sample = values[indices]
    try:
        popt, _ = curve_fit(
            monoexponential,
            times_sample,
            values_sample,
            p0=initial_guesses,
            maxfev=10000,
        )
        return float(popt[1])
    except (RuntimeError, OptimizeWarning):
        return None
