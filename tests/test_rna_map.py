"""Tests for rna_map module."""

from pathlib import Path

from mtr_analysis.rna_map import (
    RnaMapConfig,
    build_rna_map_command,
)


class TestRnaMapConfig:
    """Tests for RnaMapConfig dataclass."""

    def test_creation(self) -> None:
        """Test creating an RnaMapConfig instance."""
        config = RnaMapConfig(
            params_file=Path("/path/to/params.yml"),
            dotbracket_file=Path("/path/to/dotbracket.csv"),
            fasta_file=Path("/path/to/sequence.fasta"),
            fastq1=Path("/path/to/read1.fastq.gz"),
            fastq2=Path("/path/to/read2.fastq.gz"),
        )
        assert config.params_file == Path("/path/to/params.yml")
        assert config.fasta_file == Path("/path/to/sequence.fasta")


class TestBuildRnaMapCommand:
    """Tests for build_rna_map_command function."""

    def test_command_structure(self) -> None:
        """Test that command has correct structure."""
        config = RnaMapConfig(
            params_file=Path("/params.yml"),
            dotbracket_file=Path("/dotb.csv"),
            fasta_file=Path("/seq.fasta"),
            fastq1=Path("/r1.fastq.gz"),
            fastq2=Path("/r2.fastq.gz"),
        )

        cmd = build_rna_map_command(config)

        assert cmd[0] == "rna-map"
        assert "-fa" in cmd
        assert "-fq1" in cmd
        assert "-fq2" in cmd
        assert "--dot-bracket" in cmd
        assert "--param-file" in cmd

    def test_command_paths(self) -> None:
        """Test that paths are correctly included in command."""
        config = RnaMapConfig(
            params_file=Path("/path/params.yml"),
            dotbracket_file=Path("/path/dotb.csv"),
            fasta_file=Path("/path/seq.fasta"),
            fastq1=Path("/path/r1.fastq.gz"),
            fastq2=Path("/path/r2.fastq.gz"),
        )

        cmd = build_rna_map_command(config)

        assert "/path/seq.fasta" in cmd
        assert "/path/r1.fastq.gz" in cmd
        assert "/path/r2.fastq.gz" in cmd
        assert "/path/dotb.csv" in cmd
        assert "/path/params.yml" in cmd

    def test_command_order(self) -> None:
        """Test that command arguments are in correct order."""
        config = RnaMapConfig(
            params_file=Path("/params.yml"),
            dotbracket_file=Path("/dotb.csv"),
            fasta_file=Path("/seq.fasta"),
            fastq1=Path("/r1.fastq.gz"),
            fastq2=Path("/r2.fastq.gz"),
        )

        cmd = build_rna_map_command(config)

        # Check that each flag is followed by its value
        fa_idx = cmd.index("-fa")
        assert cmd[fa_idx + 1] == "/seq.fasta"

        fq1_idx = cmd.index("-fq1")
        assert cmd[fq1_idx + 1] == "/r1.fastq.gz"

        fq2_idx = cmd.index("-fq2")
        assert cmd[fq2_idx + 1] == "/r2.fastq.gz"
