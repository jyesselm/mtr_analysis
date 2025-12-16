"""
SLURM job templates for MTR analysis tasks.

This module provides pre-configured job templates for common
analysis tasks like RNA-MaP processing and mutation analysis.
"""

from pathlib import Path

from mtr_analysis.slurm.job import SlurmConfig, SlurmJob


def _get_absolute_path(path: str | None) -> str | None:
    """Convert a path to absolute path if provided."""
    if path is None:
        return None
    return str(Path(path).resolve())


def create_rna_map_job(
    construct: str,
    barcode_seq: str,
    output_dir: str,
    config: SlurmConfig | None = None,
    config_file: str | None = None,
) -> SlurmJob:
    """
    Create a SLURM job for RNA-MaP analysis.

    Args:
        construct: Construct name/identifier.
        barcode_seq: Barcode sequence for demultiplexing.
        output_dir: Base output directory.
        config: SLURM configuration (uses defaults if None).
        config_file: Path to config file to pass to mtr-analysis.

    Returns:
        Configured SlurmJob for RNA-MaP analysis.
    """
    if config is None:
        config = SlurmConfig(time="02:00:00", memory="8G")
    # Convert paths to absolute
    abs_output_dir = _get_absolute_path(output_dir)
    abs_config_file = _get_absolute_path(config_file)
    commands = _build_rna_map_commands(construct, barcode_seq, abs_output_dir, abs_config_file)
    return SlurmJob(name=f"rna_map_{construct}", commands=commands, config=config)


def _build_rna_map_commands(
    construct: str, barcode_seq: str, output_dir: str, config_file: str | None = None
) -> list[str]:
    """Build shell commands for RNA-MaP job."""
    # Build the mtr-analysis command with optional config
    cmd = f"mtr-analysis run single-rna-map {barcode_seq}"
    if config_file:
        cmd += f" --config {config_file}"

    # Create a unique run directory for this job to avoid file collisions
    run_dir = f"{output_dir}/run/{construct}"

    return [
        f"# Process construct: {construct}",
        "",
        "# Create unique run directory to avoid file collisions",
        f"mkdir -p {run_dir}",
        f"cd {run_dir}",
        "",
        "# Clean up any previous run in this directory",
        "rm -rf log output input",
        "",
        "# Create output directory",
        f"mkdir -p {output_dir}/{construct}",
        "",
        "# Run RNA-MaP",
        cmd,
        "",
        "# Move results",
        f"mv output {output_dir}/{construct}/",
        "",
        f'echo "Completed processing {construct}"',
    ]


def create_mutation_job(
    data_dir: str,
    sequence: str,  # noqa: ARG001 - kept for API consistency
    config: SlurmConfig | None = None,
    config_file: str | None = None,
) -> SlurmJob:
    """
    Create a SLURM job for mutation fraction analysis.

    Args:
        data_dir: Directory containing RNA-MaP output.
        sequence: Reference RNA sequence.
        config: SLURM configuration (uses defaults if None).
        config_file: Path to config file to pass to mtr-analysis.

    Returns:
        Configured SlurmJob for mutation analysis.
    """
    if config is None:
        config = SlurmConfig(time="00:30:00", memory="4G")
    # Convert paths to absolute
    abs_data_dir = _get_absolute_path(data_dir)
    abs_config_file = _get_absolute_path(config_file)
    dir_name = Path(data_dir).name
    commands = _build_mutation_commands(abs_data_dir, abs_config_file)
    return SlurmJob(name=f"mutations_{dir_name}", commands=commands, config=config)


def _build_mutation_commands(data_dir: str, config_file: str | None = None) -> list[str]:
    """Build shell commands for mutation analysis job."""
    cmd = f"mtr-analysis run process-dir {data_dir}"
    if config_file:
        cmd += f" --config {config_file}"

    return [
        f"# Analyze mutations in: {data_dir}",
        "",
        cmd,
        "",
        f'echo "Completed mutation analysis for {data_dir}"',
    ]


def create_aggregation_job(
    data_dirs: list[str],  # noqa: ARG001 - kept for API consistency
    output_file: str,
    data_dir: str = "data",
    config: SlurmConfig | None = None,
    config_file: str | None = None,
) -> SlurmJob:
    """
    Create a SLURM job for aggregating mutation results.

    Args:
        data_dirs: List of directories with mutation results.
        output_file: Path to output aggregated CSV.
        data_dir: Base data directory.
        config: SLURM configuration (uses defaults if None).
        config_file: Path to config file to pass to mtr-analysis.

    Returns:
        Configured SlurmJob for aggregation.
    """
    if config is None:
        config = SlurmConfig(time="00:15:00", memory="2G")
    # Convert paths to absolute
    abs_data_dir = _get_absolute_path(data_dir)
    abs_output_file = _get_absolute_path(output_file)
    abs_config_file = _get_absolute_path(config_file)

    cmd = f"mtr-analysis run aggregate --output {abs_output_file} --data-dir {abs_data_dir}"
    if abs_config_file:
        cmd += f" --config {abs_config_file}"

    commands = [
        "# Aggregate mutation fractions",
        "",
        cmd,
        "",
        f'echo "Aggregation complete: {abs_output_file}"',
    ]
    return SlurmJob(name="aggregate_mutations", commands=commands, config=config)


def create_fitting_job(
    input_file: str,  # noqa: ARG001 - kept for API consistency
    output_file: str,
    min_info_count: int = 1000,
    generate_plots: bool = False,
    config: SlurmConfig | None = None,
    config_file: str | None = None,
) -> SlurmJob:
    """
    Create a SLURM job for curve fitting.

    Args:
        input_file: Path to aggregated mutation fractions CSV.
        output_file: Path for kinetics output CSV.
        min_info_count: Minimum read count filter.
        generate_plots: Whether to generate plots.
        config: SLURM configuration (uses defaults if None).
        config_file: Path to config file to pass to mtr-analysis.

    Returns:
        Configured SlurmJob for curve fitting.
    """
    if config is None:
        config = SlurmConfig(time="00:30:00", memory="4G")
    # Convert paths to absolute
    abs_config_file = _get_absolute_path(config_file)

    cmd = f"mtr-analysis run fit --min-info-count {min_info_count}"
    if generate_plots:
        cmd += " --plot"
    if abs_config_file:
        cmd += f" --config {abs_config_file}"
    commands = [
        "# Fit mutation kinetics",
        "",
        cmd,
        "",
        f'echo "Fitting complete: {output_file}"',
    ]
    return SlurmJob(name="fit_kinetics", commands=commands, config=config)
