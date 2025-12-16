"""
SLURM-specific CLI commands for MTR analysis.

This module provides commands for generating and submitting SLURM jobs
to parallelize analysis across HPC clusters.
"""

import glob
from pathlib import Path

import click
import pandas as pd

from mtr_analysis.config import Config, ConfigError, load_config
from mtr_analysis.slurm.job import SlurmConfig
from mtr_analysis.slurm.runner import SlurmRunner
from mtr_analysis.slurm.templates import (
    create_aggregation_job,
    create_fitting_job,
    create_mutation_job,
    create_rna_map_job,
)
from mtr_analysis.time_parser import get_minutes_from_dir_name


def _load_config_file(config_path: str | None) -> Config | None:
    """Load config file if provided."""
    if config_path is None:
        return None
    try:
        return load_config(config_path)
    except (ConfigError, FileNotFoundError) as e:
        raise click.ClickException(f"Config error: {e}")


def _write_submit_file(
    script_paths: list, script_dir: str, stage_name: str, stage_num: int
) -> Path:
    """Write a README_SUBMIT file with sbatch commands for all scripts."""
    submit_file = Path(script_dir) / f"README_SUBMIT_{stage_num}_{stage_name.upper()}.sh"
    lines = [
        "#!/bin/bash",
        f"# Stage {stage_num}: {stage_name}",
        f"# Generated submission script for {len(script_paths)} job(s)",
        "# Usage: source " + str(submit_file),
        "",
    ]
    for path in script_paths:
        lines.append(f"sbatch {path}")
    lines.append("")
    lines.append(f'echo "Submitted {len(script_paths)} {stage_name} job(s)"')

    submit_file.write_text("\n".join(lines))
    print(f"Generated submit file: {submit_file}")
    return submit_file


@click.group()
def slurm() -> None:
    """SLURM job management commands."""
    pass


@slurm.command()
@click.argument("data_csv", type=click.Path(exists=True))
@click.option("--config", "-c", "config_file", type=click.Path(exists=True), help="Config file path")
@click.option("--output-dir", default="data", help="Output directory")
@click.option("--script-dir", default="slurm_scripts", help="Script output directory")
@click.option("--time", default=None, help="Time limit per job (default: 02:00:00)")
@click.option("--memory", default=None, help="Memory per job (default: 8G)")
@click.option("--extra-commands", default=None, help="Extra commands to run (e.g., module load)")
@click.option("--dry-run", is_flag=True, help="Generate scripts without submitting")
@click.option("--submit", is_flag=True, help="Submit jobs after generating")
def setup_rna_map(
    data_csv: str,
    config_file: str | None,
    output_dir: str,
    script_dir: str,
    time: str | None,
    memory: str | None,
    extra_commands: str | None,
    dry_run: bool,
    submit: bool,
) -> None:
    """
    Generate SLURM jobs for RNA-MaP analysis.

    Creates one job per construct, allowing parallel processing.
    Use --config to load settings from a config file.
    """
    cfg = _load_config_file(config_file)
    slurm_config = _build_slurm_config(time, memory, extra_commands, cfg, stage="rna_map")
    df = _load_constructs(data_csv)
    script_paths = _generate_rna_map_scripts(df, output_dir, script_dir, slurm_config, config_file)
    _write_submit_file(script_paths, script_dir, "rna_map", 1)
    _handle_job_submission(script_paths, script_dir, dry_run, submit)


def _build_slurm_config(
    time: str | None,
    memory: str | None,
    extra_commands: str | None = None,
    cfg: Config | None = None,
    stage: str | None = None,
) -> SlurmConfig:
    """Build SLURM configuration from CLI options and/or config file.

    CLI options override config file values. Stage-specific config values
    are used when available (e.g., rna_map_time for the rna_map stage).
    """
    # Determine defaults based on stage and config
    if cfg is not None:
        slurm_cfg = cfg.slurm
        # Get stage-specific defaults from config
        if stage == "rna_map":
            default_time = slurm_cfg.rna_map_time
            default_memory = slurm_cfg.rna_map_memory
        elif stage == "mutations":
            default_time = slurm_cfg.mutations_time
            default_memory = slurm_cfg.mutations_memory
        elif stage == "aggregation":
            default_time = slurm_cfg.aggregation_time
            default_memory = slurm_cfg.aggregation_memory
        elif stage == "fitting":
            default_time = slurm_cfg.fitting_time
            default_memory = slurm_cfg.fitting_memory
        else:
            default_time = slurm_cfg.time
            default_memory = slurm_cfg.memory
        default_extra = slurm_cfg.extra_commands
    else:
        # Fallback defaults when no config
        default_time = "01:00:00"
        default_memory = "4G"
        default_extra = None

    # CLI options override config/defaults
    final_time = time if time is not None else default_time
    final_memory = memory if memory is not None else default_memory
    final_extra = extra_commands if extra_commands is not None else default_extra

    return SlurmConfig(time=final_time, memory=final_memory, extra_commands=final_extra)


def _load_constructs(data_csv: str) -> pd.DataFrame:
    """Load and filter construct data from CSV."""
    df = pd.read_csv(data_csv)
    df = df.query("code == 'C01HP'")
    df["time"] = df["construct"].apply(get_minutes_from_dir_name)
    return df


def _generate_rna_map_scripts(
    df: pd.DataFrame, output_dir: str, script_dir: str, config: SlurmConfig,
    config_file: str | None = None
) -> list:
    """Generate SLURM scripts for each construct."""
    script_path = Path(script_dir)
    script_path.mkdir(parents=True, exist_ok=True)
    config.output_dir = script_path / "logs"
    config.output_dir.mkdir(parents=True, exist_ok=True)

    scripts = []
    for _, row in df.iterrows():
        job = create_rna_map_job(
            construct=row["construct"],
            barcode_seq=row["barcode_seq"],
            output_dir=output_dir,
            config=config,
            config_file=config_file,
        )
        path = job.write_script(script_path)
        scripts.append(path)
        print(f"Generated: {path}")
    return scripts


def _handle_job_submission(
    script_paths: list, script_dir: str, dry_run: bool, submit: bool
) -> None:
    """Handle job submission based on CLI flags."""
    print(f"\nGenerated {len(script_paths)} job scripts in {script_dir}/")
    if submit or dry_run:
        runner = SlurmRunner(script_dir=Path(script_dir), dry_run=dry_run)
        runner.submit_parallel_jobs(script_paths)
        runner.print_summary()
    elif not submit:
        print("Use --submit to submit jobs or --dry-run to preview.")


@slurm.command()
@click.option("--config", "-c", "config_file", type=click.Path(exists=True), help="Config file path")
@click.option("--data-dir", default="data", help="Data directory containing constructs")
@click.option("--script-dir", default="slurm_scripts", help="Script output directory")
@click.option("--time", default=None, help="Time limit per job (default: 00:30:00)")
@click.option("--memory", default=None, help="Memory per job (default: 4G)")
@click.option("--extra-commands", default=None, help="Extra commands to run (e.g., module load)")
@click.option("--dry-run", is_flag=True, help="Generate scripts without submitting")
@click.option("--submit", is_flag=True, help="Submit jobs after generating")
def setup_mutations(
    config_file: str | None,
    data_dir: str,
    script_dir: str,
    time: str | None,
    memory: str | None,
    extra_commands: str | None,
    dry_run: bool,
    submit: bool,
) -> None:
    """
    Generate SLURM jobs for mutation analysis.

    Creates one job per data directory.
    Use --config to load settings from a config file.
    """
    cfg = _load_config_file(config_file)
    slurm_config = _build_slurm_config(time, memory, extra_commands, cfg, stage="mutations")
    dirs = glob.glob(f"{data_dir}/*")
    if not dirs:
        print("No data directories found.")
        return
    script_paths = _generate_mutation_scripts(dirs, script_dir, slurm_config, config_file)
    _write_submit_file(script_paths, script_dir, "mutations", 2)
    _handle_job_submission(script_paths, script_dir, dry_run, submit)


def _generate_mutation_scripts(
    dirs: list, script_dir: str, config: SlurmConfig, config_file: str | None = None
) -> list:
    """Generate SLURM scripts for mutation analysis."""
    script_path = Path(script_dir)
    script_path.mkdir(parents=True, exist_ok=True)
    config.output_dir = script_path / "logs"
    config.output_dir.mkdir(parents=True, exist_ok=True)

    scripts = []
    for dir_path in dirs:
        job = create_mutation_job(data_dir=dir_path, sequence="", config=config, config_file=config_file)
        path = job.write_script(script_path)
        scripts.append(path)
        print(f"Generated: {path}")
    return scripts


@slurm.command()
@click.option("--config", "-c", "config_file", type=click.Path(exists=True), help="Config file path")
@click.option("--data-dir", default="data", help="Data directory containing constructs")
@click.option("--script-dir", default="slurm_scripts", help="Script output directory")
@click.option("--time", default=None, help="Time limit (default: 00:15:00)")
@click.option("--memory", default=None, help="Memory (default: 2G)")
@click.option("--output", default="all_mut_fractions.csv", help="Output file")
@click.option("--extra-commands", default=None, help="Extra commands to run (e.g., module load)")
@click.option("--dry-run", is_flag=True, help="Preview without submitting")
@click.option("--submit", is_flag=True, help="Submit job")
def setup_aggregation(
    config_file: str | None,
    data_dir: str,
    script_dir: str,
    time: str | None,
    memory: str | None,
    output: str,
    extra_commands: str | None,
    dry_run: bool,
    submit: bool,
) -> None:
    """Generate SLURM job for aggregating mutation results.

    Use --config to load settings from a config file.
    """
    cfg = _load_config_file(config_file)
    slurm_config = _build_slurm_config(time, memory, extra_commands, cfg, stage="aggregation")
    dirs = glob.glob(f"{data_dir}/*")
    job = create_aggregation_job(data_dirs=dirs, output_file=output, data_dir=data_dir, config=slurm_config, config_file=config_file)
    script_path = Path(script_dir)
    script_path.mkdir(parents=True, exist_ok=True)
    slurm_config.output_dir = script_path / "logs"
    slurm_config.output_dir.mkdir(parents=True, exist_ok=True)
    path = job.write_script(script_path)
    print(f"Generated: {path}")
    _write_submit_file([path], script_dir, "aggregation", 3)
    _handle_job_submission([path], script_dir, dry_run, submit)


@slurm.command()
@click.option("--config", "-c", "config_file", type=click.Path(exists=True), help="Config file path")
@click.option("--script-dir", default="slurm_scripts", help="Script output directory")
@click.option("--time", default=None, help="Time limit (default: 00:30:00)")
@click.option("--memory", default=None, help="Memory (default: 4G)")
@click.option("--min-info-count", default=1000, help="Minimum read count filter")
@click.option("--plot", is_flag=True, help="Generate plots")
@click.option("--extra-commands", default=None, help="Extra commands to run (e.g., module load)")
@click.option("--dry-run", is_flag=True, help="Preview without submitting")
@click.option("--submit", is_flag=True, help="Submit job")
def setup_fitting(
    config_file: str | None,
    script_dir: str,
    time: str | None,
    memory: str | None,
    min_info_count: int,
    plot: bool,
    extra_commands: str | None,
    dry_run: bool,
    submit: bool,
) -> None:
    """Generate SLURM job for curve fitting.

    Use --config to load settings from a config file.
    """
    cfg = _load_config_file(config_file)
    slurm_config = _build_slurm_config(time, memory, extra_commands, cfg, stage="fitting")
    job = create_fitting_job(
        input_file="all_mut_fractions.csv",
        output_file="mut_kinetics.csv",
        min_info_count=min_info_count,
        generate_plots=plot,
        config=slurm_config,
        config_file=config_file,
    )
    script_path = Path(script_dir)
    script_path.mkdir(parents=True, exist_ok=True)
    slurm_config.output_dir = script_path / "logs"
    slurm_config.output_dir.mkdir(parents=True, exist_ok=True)
    path = job.write_script(script_path)
    print(f"Generated: {path}")
    _write_submit_file([path], script_dir, "fitting", 4)
    _handle_job_submission([path], script_dir, dry_run, submit)


@slurm.command()
@click.argument("data_csv", type=click.Path(exists=True))
@click.option("--config", "-c", "config_file", type=click.Path(exists=True), help="Config file path")
@click.option("--output-dir", default="data", help="Output directory")
@click.option("--script-dir", default="slurm_scripts", help="Script output directory")
@click.option("--min-info-count", default=1000, help="Minimum read count filter")
@click.option("--plot", is_flag=True, help="Generate plots")
@click.option("--extra-commands", default=None, help="Extra commands to run (e.g., module load)")
@click.option("--dry-run", is_flag=True, help="Preview without submitting")
@click.option("--submit", is_flag=True, help="Submit all jobs")
def setup_full_pipeline(
    data_csv: str,
    config_file: str | None,
    output_dir: str,
    script_dir: str,
    min_info_count: int,
    plot: bool,
    extra_commands: str | None,
    dry_run: bool,
    submit: bool,
) -> None:
    """
    Generate complete SLURM pipeline for full analysis.

    Creates jobs for: RNA-MaP -> Mutations -> Aggregation -> Fitting
    Jobs are set up with appropriate dependencies.
    Use --config to load settings from a config file.
    """
    cfg = _load_config_file(config_file)
    script_path = Path(script_dir)
    script_path.mkdir(parents=True, exist_ok=True)
    log_dir = script_path / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    print("Setting up full analysis pipeline...")
    print(f"Scripts will be written to: {script_dir}/")

    # Generate all scripts
    df = _load_constructs(data_csv)
    all_scripts = _generate_full_pipeline_scripts(
        df, output_dir, script_path, log_dir, min_info_count, plot, extra_commands, cfg, config_file
    )

    # Generate README_SUBMIT files for each stage
    _write_submit_file(all_scripts["rna_map"], script_dir, "rna_map", 1)
    _write_submit_file(all_scripts["mutations"], script_dir, "mutations", 2)
    _write_submit_file(all_scripts["aggregation"], script_dir, "aggregation", 3)
    _write_submit_file(all_scripts["fitting"], script_dir, "fitting", 4)

    print(f"\nGenerated {len(all_scripts)} total job scripts")
    _print_pipeline_instructions(script_dir, all_scripts)

    if submit or dry_run:
        _submit_pipeline(all_scripts, script_dir, dry_run)


def _generate_full_pipeline_scripts(
    df: pd.DataFrame,
    output_dir: str,
    script_path: Path,
    log_dir: Path,
    min_info_count: int,
    plot: bool,
    extra_commands: str | None = None,
    cfg: Config | None = None,
    config_file: str | None = None,
) -> dict[str, list[Path]]:
    """Generate all scripts for full pipeline."""
    scripts: dict[str, list[Path]] = {
        "rna_map": [],
        "mutations": [],
        "aggregation": [],
        "fitting": [],
    }

    # Build stage-specific configs using config file values if available
    rna_slurm = _build_slurm_config(None, None, extra_commands, cfg, stage="rna_map")
    rna_slurm.output_dir = log_dir

    mut_slurm = _build_slurm_config(None, None, extra_commands, cfg, stage="mutations")
    mut_slurm.output_dir = log_dir

    agg_slurm = _build_slurm_config(None, None, extra_commands, cfg, stage="aggregation")
    agg_slurm.output_dir = log_dir

    fit_slurm = _build_slurm_config(None, None, extra_commands, cfg, stage="fitting")
    fit_slurm.output_dir = log_dir

    # RNA-MaP jobs
    for _, row in df.iterrows():
        job = create_rna_map_job(
            row["construct"], row["barcode_seq"], output_dir, rna_slurm, config_file
        )
        scripts["rna_map"].append(job.write_script(script_path))

    # Mutation jobs
    for _, row in df.iterrows():
        dir_path = f"{output_dir}/{row['construct']}"
        job = create_mutation_job(dir_path, "", mut_slurm, config_file)
        scripts["mutations"].append(job.write_script(script_path))

    # Aggregation job
    dirs = [f"{output_dir}/{row['construct']}" for _, row in df.iterrows()]
    job = create_aggregation_job(dirs, "all_mut_fractions.csv", output_dir, agg_slurm, config_file)
    scripts["aggregation"].append(job.write_script(script_path))

    # Fitting job
    job = create_fitting_job(
        "all_mut_fractions.csv", "mut_kinetics.csv", min_info_count, plot, fit_slurm, config_file
    )
    scripts["fitting"].append(job.write_script(script_path))

    return scripts


def _print_pipeline_instructions(script_dir: str, scripts: dict) -> None:
    """Print instructions for manual submission."""
    print("\nPipeline stages:")
    print(f"  1. RNA-MaP jobs: {len(scripts['rna_map'])} scripts")
    print(f"  2. Mutation jobs: {len(scripts['mutations'])} scripts")
    print(f"  3. Aggregation: {len(scripts['aggregation'])} script")
    print(f"  4. Fitting: {len(scripts['fitting'])} script")
    print("\nTo submit each stage (after previous stage completes):")
    print(f"  source {script_dir}/README_SUBMIT_1_RNA_MAP.sh")
    print(f"  source {script_dir}/README_SUBMIT_2_MUTATIONS.sh")
    print(f"  source {script_dir}/README_SUBMIT_3_AGGREGATION.sh")
    print(f"  source {script_dir}/README_SUBMIT_4_FITTING.sh")


def _submit_pipeline(scripts: dict, script_dir: str, dry_run: bool) -> None:
    """Submit pipeline jobs with dependencies."""
    runner = SlurmRunner(script_dir=Path(script_dir), dry_run=dry_run)

    # Submit RNA-MaP jobs in parallel
    rna_map_ids = runner.submit_parallel_jobs(scripts["rna_map"])

    # Submit mutation jobs with dependency on all RNA-MaP jobs
    dependency = ",".join(rna_map_ids) if rna_map_ids else None
    for script in scripts["mutations"]:
        runner.submit_job(script, dependency=dependency)

    # Submit aggregation after mutations
    mut_ids = list(runner.submitted_jobs.keys())[-len(scripts["mutations"]) :]
    agg_dep = ",".join(mut_ids) if mut_ids else None
    agg_id = runner.submit_job(scripts["aggregation"][0], dependency=agg_dep)

    # Submit fitting after aggregation
    runner.submit_job(scripts["fitting"][0], dependency=agg_id)

    runner.print_summary()
