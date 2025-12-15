"""Pytest fixtures and configuration for mtr_analysis tests."""

import tempfile
from collections.abc import Generator
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_sequence() -> str:
    """Return a sample RNA sequence for testing."""
    return "GGAAGATCGAGTAGATCAAAGGAGGCTGACCGACCCCCCGAGCTTCGGCTCGGGGACAACTA"


@pytest.fixture
def sample_bitvector_content() -> str:
    """Return sample bitvector file content for testing."""
    # Header lines
    header = "header1\nheader2\nheader3\n"
    # Data lines: bitvector string with mutations, then tab, then mutation count
    # Position 0 has 'A' mutation in first read
    # Position 5 has 'G' mutation in second read
    seq_len = 120
    bitvector1 = "." * seq_len
    bitvector2 = "." * 5 + "G" + "." * (seq_len - 6)
    bitvector3 = "." * 86 + "T" + "." * (seq_len - 87)  # mutation at target pos

    data = f"read1\t{bitvector1}\t0\nread2\t{bitvector2}\t1\nread3\t{bitvector3}\t1\n"
    return header + data


@pytest.fixture
def sample_times() -> np.ndarray:
    """Return sample time points for fitting tests."""
    return np.array([0, 15, 30, 60, 120, 240])


@pytest.fixture
def sample_values() -> np.ndarray:
    """Return sample mutation fraction values for fitting tests."""
    # Simulate monoexponential growth: y = 0.8 * (1 - exp(-0.01 * t))
    times = np.array([0, 15, 30, 60, 120, 240])
    y_max, k = 0.8, 0.01
    values = y_max * (1 - np.exp(-k * times))
    # Add small noise
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 0.01, len(times))
    return np.clip(values + noise, 0, 1)
