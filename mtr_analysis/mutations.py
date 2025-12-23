"""
Mutation processing and analysis.

This module provides utilities for parsing mutation data from bitvector
files and computing mutation statistics.
"""

from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path


@dataclass
class MutationResult:
    """Results from processing mutation data."""

    mutation_counts: dict[str, int]
    info_counts: dict[str, int]
    position_histogram: list[int]
    total_reads: int


def process_mutations(sequence: str, path: Path | str) -> MutationResult:
    """
    Process mutation data from a bitvector file.

    Reads a bitvector file, filters reads by mutation count, and computes
    mutation frequencies at a target position.

    Args:
        sequence: Reference RNA sequence.
        path: Path to the bitvector file.

    Returns:
        MutationResult containing counts and histogram data.
    """
    lines = _read_bitvector_lines(path)
    return _analyze_mutations(sequence, lines)


def _read_bitvector_lines(path: Path | str) -> Iterator[str]:
    """
    Read and preprocess bitvector file lines.

    Yields lines one at a time to avoid loading entire file into memory.

    Args:
        path: Path to the bitvector file.

    Yields:
        Data lines (header lines skipped).
    """
    with open(path) as f:
        # Skip first 3 header lines
        for _ in range(3):
            next(f, None)
        # Yield remaining lines one at a time
        for line in f:
            yield line.strip()


def _analyze_mutations(sequence: str, lines: Iterator[str]) -> MutationResult:
    """
    Analyze mutation data from bitvector lines.

    Args:
        sequence: Reference RNA sequence.
        lines: Iterator of preprocessed bitvector data lines.

    Returns:
        MutationResult with computed statistics.
    """
    mutation_counts: dict[str, int] = defaultdict(int)
    info_counts: dict[str, int] = defaultdict(int)
    histogram = [0] * len(sequence)
    total_reads = 0
    target_position = 86

    for line in lines:
        result = _process_single_read(line, sequence, target_position)
        if result is None:
            continue
        total_reads += 1
        mutation_name, is_mutated, positions = result
        _update_histogram(histogram, positions)
        mutation_counts[mutation_name] += int(is_mutated)
        info_counts[mutation_name] += 1

    return MutationResult(
        mutation_counts=dict(mutation_counts),
        info_counts=dict(info_counts),
        position_histogram=histogram,
        total_reads=total_reads,
    )


def _process_single_read(
    line: str, sequence: str, target_position: int
) -> tuple[str, bool, list[int]] | None:
    """
    Process a single read from the bitvector file.

    Args:
        line: Tab-separated line from bitvector file.
        sequence: Reference sequence.
        target_position: Position to check for mutation.

    Returns:
        Tuple of (mutation_name, is_target_mutated, mutation_positions),
        or None if read should be filtered.
    """
    fields = line.split("\t")
    mutation_count = int(fields[-1])
    if mutation_count > 3:
        return None
    bitvector = fields[1]
    if bitvector[target_position] == ".":
        return None
    positions = _find_mutation_positions(bitvector)
    mutation_name = _build_mutation_name(sequence, positions, target_position)
    is_mutated = target_position in [p for p, _ in positions]
    return mutation_name, is_mutated, [p for p, _ in positions]


def _find_mutation_positions(bitvector: str) -> list[tuple[int, str]]:
    """
    Find positions of mutations in a bitvector string.

    Args:
        bitvector: String where mutations are marked with nucleotide letters.

    Returns:
        List of (position, nucleotide) tuples for each mutation.
    """
    nucleotides = {"T", "C", "A", "G"}
    return [(i, char) for i, char in enumerate(bitvector) if char in nucleotides]


def _build_mutation_name(
    sequence: str, positions: list[tuple[int, str]], target_position: int
) -> str:
    """
    Build a mutation name string from mutation positions.

    Args:
        sequence: Reference sequence.
        positions: List of (position, nucleotide) tuples.
        target_position: Position to exclude from the name.

    Returns:
        Mutation name like 'A42G_C55T' or 'WT' for wild type.
    """
    parts = []
    for pos, char in positions:
        if pos == target_position:
            continue
        original = sequence[pos]
        parts.append(f"{original}{pos + 1}{char}")
    if not parts:
        return "WT"
    return "_".join(parts)


def _update_histogram(histogram: list[int], positions: list[int]) -> None:
    """
    Update position histogram with mutation positions.

    Args:
        histogram: List to update in place.
        positions: Positions to increment.
    """
    for pos in positions:
        histogram[pos] += 1


def compute_mutation_fractions(
    mutation_counts: dict[str, int], info_counts: dict[str, int]
) -> dict[str, float]:
    """
    Compute mutation fractions from counts.

    Args:
        mutation_counts: Number of mutations per variant.
        info_counts: Total reads per variant.

    Returns:
        Dictionary mapping variant names to mutation fractions.
    """
    fractions = {}
    for name in mutation_counts:
        if info_counts[name] > 0:
            fractions[name] = mutation_counts[name] / info_counts[name]
        else:
            fractions[name] = 0.0
    return fractions
