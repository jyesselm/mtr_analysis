"""Tests for mutations module."""

from pathlib import Path

from mtr_analysis.mutations import (
    MutationResult,
    _build_mutation_name,
    _find_mutation_positions,
    _read_bitvector_lines,
    _update_histogram,
    compute_mutation_fractions,
    process_mutations,
)


class TestMutationResult:
    """Tests for MutationResult dataclass."""

    def test_creation(self) -> None:
        """Test creating a MutationResult instance."""
        result = MutationResult(
            mutation_counts={"WT": 10, "A42G": 5},
            info_counts={"WT": 100, "A42G": 50},
            position_histogram=[0] * 100,
            total_reads=150,
        )
        assert result.mutation_counts["WT"] == 10
        assert result.info_counts["A42G"] == 50
        assert result.total_reads == 150


class TestProcessMutations:
    """Tests for process_mutations function."""

    def test_process_mutations(
        self, temp_dir: Path, sample_bitvector_content: str
    ) -> None:
        """Test processing mutations from a bitvector file."""
        # Write test file
        bitvector_file = temp_dir / "test_bitvectors.txt"
        bitvector_file.write_text(sample_bitvector_content)

        # Use a sequence long enough for the target position
        sequence = "A" * 120
        result = process_mutations(sequence, bitvector_file)

        assert isinstance(result, MutationResult)
        assert result.total_reads >= 0


class TestReadBitvectorLines:
    """Tests for _read_bitvector_lines helper function."""

    def test_read_and_skip_header(self, temp_dir: Path) -> None:
        """Test that header lines are skipped."""
        content = "header1\nheader2\nheader3\ndata1\ndata2\n"
        test_file = temp_dir / "test.txt"
        test_file.write_text(content)

        lines = list(_read_bitvector_lines(test_file))

        assert len(lines) == 2
        assert lines[0] == "data1"
        assert lines[1] == "data2"


class TestFindMutationPositions:
    """Tests for _find_mutation_positions helper function."""

    def test_no_mutations(self) -> None:
        """Test bitvector with no mutations."""
        bitvector = "." * 50
        positions = _find_mutation_positions(bitvector)
        assert positions == []

    def test_single_mutation(self) -> None:
        """Test bitvector with single mutation."""
        bitvector = "." * 10 + "A" + "." * 39
        positions = _find_mutation_positions(bitvector)
        assert positions == [(10, "A")]

    def test_multiple_mutations(self) -> None:
        """Test bitvector with multiple mutations."""
        bitvector = "A" + "." * 8 + "G" + "." * 40
        positions = _find_mutation_positions(bitvector)
        assert positions == [(0, "A"), (9, "G")]

    def test_all_nucleotides(self) -> None:
        """Test that all nucleotide types are detected."""
        bitvector = "ATCG" + "." * 46
        positions = _find_mutation_positions(bitvector)
        assert len(positions) == 4
        assert {char for _, char in positions} == {"A", "T", "C", "G"}


class TestBuildMutationName:
    """Tests for _build_mutation_name helper function."""

    def test_wild_type(self) -> None:
        """Test that empty positions returns WT."""
        name = _build_mutation_name("ACGT", [], target_position=2)
        assert name == "WT"

    def test_single_mutation(self) -> None:
        """Test building name for single mutation."""
        sequence = "ACGTACGT"
        positions = [(0, "T")]  # A->T at position 0
        name = _build_mutation_name(sequence, positions, target_position=5)
        assert name == "A1T"

    def test_multiple_mutations(self) -> None:
        """Test building name for multiple mutations."""
        sequence = "ACGTACGT"
        positions = [(0, "G"), (2, "A")]  # A->G at 0, G->A at 2
        name = _build_mutation_name(sequence, positions, target_position=5)
        assert name == "A1G_G3A"

    def test_target_position_excluded(self) -> None:
        """Test that target position is excluded from name."""
        sequence = "ACGTACGT"
        positions = [(0, "T"), (5, "G")]
        name = _build_mutation_name(sequence, positions, target_position=5)
        assert name == "A1T"  # Position 5 excluded


class TestUpdateHistogram:
    """Tests for _update_histogram helper function."""

    def test_update_single_position(self) -> None:
        """Test updating histogram at single position."""
        histogram = [0, 0, 0, 0, 0]
        _update_histogram(histogram, [2])
        assert histogram == [0, 0, 1, 0, 0]

    def test_update_multiple_positions(self) -> None:
        """Test updating histogram at multiple positions."""
        histogram = [0, 0, 0, 0, 0]
        _update_histogram(histogram, [0, 2, 4])
        assert histogram == [1, 0, 1, 0, 1]

    def test_update_empty_positions(self) -> None:
        """Test that empty positions list doesn't change histogram."""
        histogram = [1, 2, 3]
        _update_histogram(histogram, [])
        assert histogram == [1, 2, 3]


class TestComputeMutationFractions:
    """Tests for compute_mutation_fractions function."""

    def test_basic_fractions(self) -> None:
        """Test computing basic mutation fractions."""
        mutation_counts = {"WT": 50, "A42G": 20}
        info_counts = {"WT": 100, "A42G": 100}

        fractions = compute_mutation_fractions(mutation_counts, info_counts)

        assert fractions["WT"] == 0.5
        assert fractions["A42G"] == 0.2

    def test_zero_info_count(self) -> None:
        """Test handling of zero info count."""
        mutation_counts = {"WT": 10}
        info_counts = {"WT": 0}

        fractions = compute_mutation_fractions(mutation_counts, info_counts)

        assert fractions["WT"] == 0.0

    def test_zero_mutation_count(self) -> None:
        """Test handling of zero mutation count."""
        mutation_counts = {"WT": 0}
        info_counts = {"WT": 100}

        fractions = compute_mutation_fractions(mutation_counts, info_counts)

        assert fractions["WT"] == 0.0
