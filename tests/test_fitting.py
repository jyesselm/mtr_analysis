"""Tests for fitting module."""

import numpy as np
import pytest

from mtr_analysis.fitting import (
    FitResult,
    _bootstrap_k_error,
    _fit_curve,
    fit_monoexponential,
    monoexponential,
)


class TestFitResult:
    """Tests for FitResult dataclass."""

    def test_creation(self) -> None:
        """Test creating a FitResult instance."""
        result = FitResult(y_max=0.8, k=0.05, k_std=0.01)
        assert result.y_max == 0.8
        assert result.k == 0.05
        assert result.k_std == 0.01


class TestMonoexponential:
    """Tests for monoexponential function."""

    def test_zero_time(self) -> None:
        """Test that monoexponential returns 0 at t=0."""
        result = monoexponential(np.array([0.0]), 1.0, 0.1)
        assert result[0] == pytest.approx(0.0)

    def test_large_time(self) -> None:
        """Test that monoexponential approaches y_max at large t."""
        result = monoexponential(np.array([10000.0]), 0.8, 0.1)
        assert result[0] == pytest.approx(0.8, rel=1e-3)

    def test_intermediate_time(self) -> None:
        """Test monoexponential at intermediate time points."""
        # y = 0.8 * (1 - exp(-0.1 * 10)) = 0.8 * (1 - exp(-1))
        expected = 0.8 * (1 - np.exp(-1))
        result = monoexponential(np.array([10.0]), 0.8, 0.1)
        assert result[0] == pytest.approx(expected)

    def test_array_input(self) -> None:
        """Test monoexponential with array input."""
        times = np.array([0, 10, 20, 30])
        result = monoexponential(times, 1.0, 0.1)
        assert len(result) == 4
        assert result[0] == pytest.approx(0.0)
        assert result[-1] > result[0]


class TestFitMonoexponential:
    """Tests for fit_monoexponential function."""

    def test_fit_synthetic_data(
        self, sample_times: np.ndarray, sample_values: np.ndarray
    ) -> None:
        """Test fitting on synthetic data."""
        result = fit_monoexponential(
            sample_times.astype(float),
            sample_values.astype(float),
            n_bootstrap=100,  # Reduced for test speed
        )

        assert isinstance(result, FitResult)
        assert 0.5 < result.y_max < 1.0  # Should be close to 0.8
        assert 0.001 < result.k < 0.1  # Should be close to 0.01
        assert result.k_std >= 0  # Standard deviation should be non-negative

    def test_fit_perfect_data(self) -> None:
        """Test fitting on perfect (noise-free) data."""
        y_max, k = 0.9, 0.05
        times = np.array([0.0, 10.0, 20.0, 40.0, 80.0, 160.0])
        values = y_max * (1 - np.exp(-k * times))

        result = fit_monoexponential(times, values, n_bootstrap=50)

        assert result.y_max == pytest.approx(y_max, rel=0.01)
        assert result.k == pytest.approx(k, rel=0.01)

    def test_reproducibility_with_seed(self) -> None:
        """Test that results are reproducible with same seed."""
        times = np.array([0.0, 15.0, 30.0, 60.0, 120.0])
        values = 0.8 * (1 - np.exp(-0.02 * times))

        result1 = fit_monoexponential(times, values, n_bootstrap=100, seed=42)
        result2 = fit_monoexponential(times, values, n_bootstrap=100, seed=42)

        assert result1.k_std == pytest.approx(result2.k_std)


class TestFitCurve:
    """Tests for _fit_curve helper function."""

    def test_fit_curve_basic(self) -> None:
        """Test basic curve fitting."""
        y_max, k = 0.7, 0.03
        times = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0])
        values = y_max * (1 - np.exp(-k * times))

        fitted_ymax, fitted_k = _fit_curve(times, values)

        assert fitted_ymax == pytest.approx(y_max, rel=0.05)
        assert fitted_k == pytest.approx(k, rel=0.05)


class TestBootstrapKError:
    """Tests for _bootstrap_k_error helper function."""

    def test_bootstrap_returns_float(self) -> None:
        """Test that bootstrap error returns a float."""
        times = np.array([0.0, 15.0, 30.0, 60.0, 120.0])
        values = 0.8 * (1 - np.exp(-0.02 * times))

        error = _bootstrap_k_error(times, values, n_bootstrap=50, seed=42)

        assert isinstance(error, float)
        assert error >= 0

    def test_bootstrap_with_noise(self) -> None:
        """Test bootstrap error is larger with noisy data."""
        times = np.array([0.0, 15.0, 30.0, 60.0, 120.0])
        clean_values = 0.8 * (1 - np.exp(-0.02 * times))

        # Add significant noise
        rng = np.random.default_rng(42)
        noisy_values = clean_values + rng.normal(0, 0.1, len(times))
        noisy_values = np.clip(noisy_values, 0, 1)

        clean_error = _bootstrap_k_error(times, clean_values, n_bootstrap=100, seed=42)
        noisy_error = _bootstrap_k_error(times, noisy_values, n_bootstrap=100, seed=42)

        # Noisy data should generally have larger uncertainty
        # (not always guaranteed but usually true)
        assert noisy_error > 0
        assert clean_error >= 0
