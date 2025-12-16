"""
SLURM job templates for MTR analysis tasks.

This module provides pre-configured job templates for common
analysis tasks like RNA-MaP processing and mutation analysis.
"""

from pathlib import Path

from mtr_analysis.slurm.job import SlurmConfig, SlurmJob


def create_rna_map_job(
    construct: str,
    barcode_seq: str,
    output_dir: str,
    config: SlurmConfig | None = None,
) -> SlurmJob:
    """
    Create a SLURM job for RNA-MaP analysis.

    Args:
        construct: Construct name/identifier.
        barcode_seq: Barcode sequence for demultiplexing.
        output_dir: Base output directory.
        config: SLURM configuration (uses defaults if None).

    Returns:
        Configured SlurmJob for RNA-MaP analysis.
    """
    if config is None:
        config = SlurmConfig(time="02:00:00", memory="8G")
    commands = _build_rna_map_commands(construct, barcode_seq, output_dir)
    return SlurmJob(name=f"rna_map_{construct}", commands=commands, config=config)


def _build_rna_map_commands(
    construct: str, barcode_seq: str, output_dir: str
) -> list[str]:
    """Build shell commands for RNA-MaP job."""
    return [
        f"# Process construct: {construct}",
        "",
        "# Clean up any previous run",
        "rm -rf log output input",
        "",
        "# Create output directory",
        f"mkdir -p {output_dir}/{construct}",
        "",
        "# Run RNA-MaP",
        f"mtr-analysis run single-rna-map {barcode_seq}",
        "",
        "# Move results",
        f"mv output {output_dir}/{construct}/",
        "",
        'echo "Completed processing {construct}"',
    ]


def create_mutation_job(
    data_dir: str,
    sequence: str,  # noqa: ARG001 - kept for API consistency
    config: SlurmConfig | None = None,
) -> SlurmJob:
    """
    Create a SLURM job for mutation fraction analysis.

    Args:
        data_dir: Directory containing RNA-MaP output.
        sequence: Reference RNA sequence.
        config: SLURM configuration (uses defaults if None).

    Returns:
        Configured SlurmJob for mutation analysis.
    """
    if config is None:
        config = SlurmConfig(time="00:30:00", memory="4G")
    dir_name = Path(data_dir).name
    commands = _build_mutation_commands(data_dir)
    return SlurmJob(name=f"mutations_{dir_name}", commands=commands, config=config)


def _build_mutation_commands(data_dir: str) -> list[str]:
    """Build shell commands for mutation analysis job."""
    return [
        f"# Analyze mutations in: {data_dir}",
        "",
        f"mtr-analysis run process-dir {data_dir}",
        "",
        f'echo "Completed mutation analysis for {data_dir}"',
    ]


def create_aggregation_job(
    data_dirs: list[str],  # noqa: ARG001 - kept for API consistency
    output_file: str,
    data_dir: str = "data",
    config: SlurmConfig | None = None,
) -> SlurmJob:
    """
    Create a SLURM job for aggregating mutation results.

    Args:
        data_dirs: List of directories with mutation results.
        output_file: Path to output aggregated CSV.
        data_dir: Base data directory.
        config: SLURM configuration (uses defaults if None).

    Returns:
        Configured SlurmJob for aggregation.
    """
    if config is None:
        config = SlurmConfig(time="00:15:00", memory="2G")
    commands = [
        "# Aggregate mutation fractions",
        "",
        f"mtr-analysis run aggregate --output {output_file} --data-dir {data_dir}",
        "",
        f'echo "Aggregation complete: {output_file}"',
    ]
    return SlurmJob(name="aggregate_mutations", commands=commands, config=config)


def create_fitting_job(
    input_file: str,  # noqa: ARG001 - kept for API consistency
    output_file: str,
    min_info_count: int = 1000,
    generate_plots: bool = False,
    config: SlurmConfig | None = None,
) -> SlurmJob:
    """
    Create a SLURM job for curve fitting.

    Args:
        input_file: Path to aggregated mutation fractions CSV.
        output_file: Path for kinetics output CSV.
        min_info_count: Minimum read count filter.
        generate_plots: Whether to generate plots.
        config: SLURM configuration (uses defaults if None).

    Returns:
        Configured SlurmJob for curve fitting.
    """
    if config is None:
        config = SlurmConfig(time="00:30:00", memory="4G")
    cmd = f"mtr-analysis run fit --min-info-count {min_info_count}"
    if generate_plots:
        cmd += " --plot"
    commands = [
        "# Fit mutation kinetics",
        "",
        cmd,
        "",
        f'echo "Fitting complete: {output_file}"',
    ]
    return SlurmJob(name="fit_kinetics", commands=commands, config=config)
