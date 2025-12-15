"""
Command-line interface for MTR analysis.

This module provides CLI commands for running RNA-MaP analysis,
computing mutation fractions, and fitting kinetic data.
"""

import glob
import os
from pathlib import Path

import click
import pandas as pd

from mtr_analysis.fitting import fit_monoexponential
from mtr_analysis.mutations import compute_mutation_fractions, process_mutations
from mtr_analysis.plotting import create_kinetics_plot, save_plot
from mtr_analysis.rna_map import run_rna_map_for_barcode
from mtr_analysis.time_parser import get_minutes_from_dir_name

# Default sequence used in analysis
DEFAULT_SEQUENCE = (
    "GGAAGATCGAGTAGATCAAAGGAGGCTGACCGACCCCCCGAGCTTCGGCTCGGGGACAACTA"
    "GACATACAGTATCTTCGGATACTGAGCCTCCACAAAGAAACAACAACAACAAC"
)


@click.group()
def cli() -> None:
    """MTR ribozyme analysis toolkit."""
    pass


@cli.command()
@click.argument("data_csv", type=click.Path(exists=True))
@click.option("--data-dir", default="data", help="Output directory for data")
def run_rna_map(data_csv: str, data_dir: str) -> None:
    """
    Run RNA-MaP on all constructs in data CSV.

    DATA_CSV: Path to CSV file with construct and barcode information.
    """
    df = _load_and_filter_data(data_csv)
    _print_construct_info(df)
    _process_all_constructs(df, data_dir)


def _load_and_filter_data(data_csv: str) -> pd.DataFrame:
    """Load CSV and filter for C01HP constructs."""
    df = pd.read_csv(data_csv)
    df = df.query("code == 'C01HP'")
    df["time"] = df["construct"].apply(get_minutes_from_dir_name)
    return df


def _print_construct_info(df: pd.DataFrame) -> None:
    """Print construct information to console."""
    print("Constructs:")
    print(df[["construct", "barcode_seq", "time"]].to_string(index=False))


def _process_all_constructs(df: pd.DataFrame, data_dir: str) -> None:
    """Process all constructs and save results."""
    os.makedirs(data_dir, exist_ok=True)
    data = []
    for _, row in df.iterrows():
        result = _process_single_construct(row, data_dir)
        data.append(result)
    df_summary = pd.DataFrame(data)
    print(df_summary)


def _process_single_construct(row: pd.Series, data_dir: str) -> dict:
    """Process a single construct and return summary."""
    os.system("rm -rf log output input")
    construct_dir = f"{data_dir}/{row['construct']}"
    os.makedirs(construct_dir, exist_ok=True)
    run_rna_map_for_barcode(row["barcode_seq"])
    os.system(f"mv output {construct_dir}")
    return _extract_summary(row["construct"], row["time"], construct_dir)


def _extract_summary(construct: str, time: int, construct_dir: str) -> dict:
    """Extract summary statistics from RNA-MaP output."""
    summary_path = f"{construct_dir}/output/BitVector_Files/summary.csv"
    df_summary = pd.read_csv(summary_path)
    row_data = df_summary.iloc[0]
    return {
        "construct": construct,
        "time": time,
        "reads": row_data["reads"],
        "aligned": row_data["aligned"],
    }


@cli.command()
@click.argument("barcode_seq")
def run_single_rna_map(barcode_seq: str) -> None:
    """
    Run RNA-MaP for a single barcode.

    BARCODE_SEQ: Barcode sequence identifier.
    """
    run_rna_map_for_barcode(barcode_seq)
    print(f"Completed RNA-MaP for barcode: {barcode_seq}")


@cli.command()
@click.option("--sequence", default=DEFAULT_SEQUENCE, help="Reference sequence")
@click.option("--data-dir", default="data", help="Data directory containing constructs")
def get_mutation_fractions(sequence: str, data_dir: str) -> None:
    """Compute mutation fractions for all data directories."""
    dirs = glob.glob(f"{data_dir}/*")
    if not dirs:
        print("No data directories found.")
        return
    dfs = [_process_directory(d, sequence) for d in dirs]
    _save_combined_results(dfs)


def _process_directory(dir_path: str, sequence: str) -> pd.DataFrame:
    """Process a single directory and save results."""
    bitvector_path = f"{dir_path}/output/BitVector_Files/mtr1_mut_lib_wt_bitvectors.txt"
    result = process_mutations(sequence, bitvector_path)
    df = _create_mutation_dataframe(result)
    df.to_csv(f"{dir_path}/mut_fractions.csv", index=False)
    df["time"] = get_minutes_from_dir_name(dir_path)
    return df


def _create_mutation_dataframe(result) -> pd.DataFrame:
    """Create DataFrame from mutation processing result."""
    fractions = compute_mutation_fractions(result.mutation_counts, result.info_counts)
    data = [
        [name, result.mutation_counts[name], result.info_counts[name], frac]
        for name, frac in fractions.items()
    ]
    df = pd.DataFrame(data, columns=["mut", "mut_count", "info_count", "mut_fraction"])
    return df.sort_values(by="mut_fraction", ascending=False)


def _save_combined_results(dfs: list) -> None:
    """Combine and save all mutation fraction results."""
    df = pd.concat(dfs)
    df.to_csv("all_mut_fractions.csv", index=False)
    print("Saved combined results to all_mut_fractions.csv")


@cli.command()
@click.argument("dir_path")
@click.option("--sequence", default=DEFAULT_SEQUENCE, help="Reference sequence")
def process_single_dir(dir_path: str, sequence: str) -> None:
    """
    Process mutations for a single directory.

    DIR_PATH: Path to directory containing RNA-MaP output.
    """
    _process_directory(dir_path, sequence)
    print(f"Processed: {dir_path}")


@cli.command()
@click.option("--output", default="all_mut_fractions.csv", help="Output file path")
@click.option("--data-dir", default="data", help="Data directory containing constructs")
def aggregate_mutations(output: str, data_dir: str) -> None:
    """Aggregate mutation fractions from all data directories."""
    dirs = glob.glob(f"{data_dir}/*")
    dfs = []
    for dir_path in dirs:
        csv_path = f"{dir_path}/mut_fractions.csv"
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            df["time"] = get_minutes_from_dir_name(dir_path)
            dfs.append(df)
    if dfs:
        combined = pd.concat(dfs)
        combined.to_csv(output, index=False)
        print(f"Aggregated {len(dfs)} files to {output}")
    else:
        print("No mutation fraction files found.")


@cli.command()
@click.option("--min-info-count", type=int, default=1000, help="Minimum read count")
@click.option("--plot", is_flag=True, help="Generate plots")
@click.option("--input-file", default="all_mut_fractions.csv", help="Input CSV")
@click.option("--output-file", default="mut_kinetics.csv", help="Output CSV")
def fit_mut_fractions(
    min_info_count: int, plot: bool, input_file: str, output_file: str
) -> None:
    """Fit monoexponential curves to mutation fraction data."""
    if plot:
        os.makedirs("plots", exist_ok=True)
    df = _load_and_filter_fractions(input_file, min_info_count)
    results = _fit_all_mutations(df, plot)
    _save_kinetics_results(results, output_file)


def _load_and_filter_fractions(path: str, min_count: int) -> pd.DataFrame:
    """Load mutation fractions and apply filters."""
    df = pd.read_csv(path)
    df = df.query(f"info_count >= {min_count}")
    return df.sort_values(by="time")


def _fit_all_mutations(df: pd.DataFrame, generate_plots: bool) -> list:
    """Fit curves for all mutations with sufficient data."""
    results = []
    for mutation, group in df.groupby("mut"):
        group = group.sort_values(by="time")
        if len(group) < 5:
            continue
        fit_result = _fit_single_mutation(str(mutation), group, generate_plots)
        if fit_result:
            results.append(fit_result)
    return results


def _fit_single_mutation(
    mutation: str, group: pd.DataFrame, generate_plots: bool
) -> list | None:
    """Fit curve for a single mutation and optionally plot."""
    import numpy as np

    times = np.asarray(group["time"].values, dtype=np.float64)
    values = np.asarray(group["mut_fraction"].values, dtype=np.float64)
    try:
        result = fit_monoexponential(times, values)
    except RuntimeError:
        return None
    if generate_plots:
        _generate_mutation_plot(mutation, times, values, result)
    return [mutation, result.y_max, result.k, result.k_std]


def _generate_mutation_plot(mutation: str, times, values, result) -> None:
    """Generate and save plot for a mutation."""
    fig = create_kinetics_plot(times, values, result, mutation)
    save_plot(fig, Path(f"plots/{mutation}.png"))


def _save_kinetics_results(results: list, output_file: str) -> None:
    """Save fitting results to CSV."""
    df = pd.DataFrame(results, columns=["mut", "y_max", "k", "k_std"])
    df = df.sort_values(by="k", ascending=False)
    df.to_csv(output_file, index=False)
    print(f"Saved kinetics results to {output_file}")


# Import and register SLURM commands (must be at bottom to avoid circular imports)
from mtr_analysis.cli_slurm import slurm  # noqa: E402

cli.add_command(slurm)


if __name__ == "__main__":
    cli()
