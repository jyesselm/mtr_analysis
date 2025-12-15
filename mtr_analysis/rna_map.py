"""
RNA-MaP execution utilities.

This module handles running the rna-map tool for processing
sequencing data from RNA structure probing experiments.
"""

import subprocess
from contextlib import ExitStack
from dataclasses import dataclass
from importlib.resources import as_file, files
from pathlib import Path


@dataclass
class RnaMapConfig:
    """Configuration for RNA-MaP execution."""

    params_file: Path
    dotbracket_file: Path
    fasta_file: Path
    fastq1: Path
    fastq2: Path


def get_default_resource_paths() -> tuple[Path, Path, Path]:
    """
    Get paths to default resource files.

    Returns:
        Tuple of (params_path, dotbracket_path, fasta_path).
    """
    data_files = files("mtr_analysis.data")
    with ExitStack() as stack:
        params = stack.enter_context(as_file(data_files.joinpath("params.yml")))
        dotb = stack.enter_context(as_file(data_files.joinpath("C01HP.csv")))
        fa = stack.enter_context(as_file(data_files.joinpath("C01HP.fasta")))
        return Path(params), Path(dotb), Path(fa)


def build_rna_map_command(config: RnaMapConfig) -> list[str]:
    """
    Build command line arguments for rna-map.

    Args:
        config: RNA-MaP configuration.

    Returns:
        List of command line arguments.
    """
    return [
        "rna-map",
        "-fa",
        str(config.fasta_file),
        "-fq1",
        str(config.fastq1),
        "-fq2",
        str(config.fastq2),
        "--dot-bracket",
        str(config.dotbracket_file),
        "--param-file",
        str(config.params_file),
    ]


def run_rna_map(
    config: RnaMapConfig, check: bool = True
) -> subprocess.CompletedProcess:
    """
    Execute rna-map with the given configuration.

    Args:
        config: RNA-MaP configuration.
        check: Whether to raise on non-zero exit code.

    Returns:
        Completed process result.
    """
    cmd = build_rna_map_command(config)
    return subprocess.run(cmd, check=check)


def run_rna_map_for_barcode(barcode_seq: str) -> subprocess.CompletedProcess:
    """
    Run RNA-MaP for a specific barcode using default resources.

    Args:
        barcode_seq: Barcode sequence identifier.

    Returns:
        Completed process result.
    """
    with ExitStack() as stack:
        data_files = files("mtr_analysis.data")
        params = stack.enter_context(as_file(data_files.joinpath("params.yml")))
        dotb = stack.enter_context(as_file(data_files.joinpath("C01HP.csv")))
        fa = stack.enter_context(as_file(data_files.joinpath("C01HP.fasta")))

        config = RnaMapConfig(
            params_file=Path(params),
            dotbracket_file=Path(dotb),
            fasta_file=Path(fa),
            fastq1=Path(f"demultiplexed/{barcode_seq}/test_R2.fastq.gz"),
            fastq2=Path(f"demultiplexed/{barcode_seq}/test_R1.fastq.gz"),
        )
        return run_rna_map(config)
