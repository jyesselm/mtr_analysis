import pandas as pd
import click
import os
import glob
import json
from collections import defaultdict
import matplotlib.pyplot as plt
import re
import numpy as np
from scipy.optimize import curve_fit


def get_mins_fm_dir_name(dir_name):
    """
    Extracts the number of minutes from a directory name.
    Handles both 'min' and 'hr' (converts hours to minutes).
    Examples:
        'data/mtr1_mut_lib_t15min' -> 15
        'data/mtr1_mut_lib_t4hr' -> 240
        'data/mtr1_mut_lib_t0' -> 0
    """
    import re

    if dir_name.count("/"):
        base = os.path.basename(dir_name)
    else:
        base = dir_name
    match = re.search(r"_t(\d+)(min|hr)?", base)
    if match:
        num = int(match.group(1))
        unit = match.group(2)
        if unit == "hr":
            return num * 60
        else:
            return num
    # Handle case like 'data/mtr1_mut_lib_t0'
    match = re.search(r"_t(\d+)$", base)
    if match:
        return int(match.group(1))
    raise ValueError(f"Could not parse minutes from directory name: {dir_name}")


def process_mutations(sequence, path):
    with open(path, "r") as f:
        lines = []
        for i, line in enumerate(f):
            lines.append(line.strip())
    lines = lines[3:]
    muts = defaultdict(int)
    infos = defaultdict(int)
    histo = [0] * len(sequence)
    count = 0
    max_index = 86
    for l in lines:
        spl = l.split("\t")
        n_mut = int(spl[-1])
        if n_mut > 3:
            continue
        # Find all positions of T, C, A, G in spl[1]
        if spl[1][max_index] == ".":
            continue
        count += 1
        positions = []
        for i, char in enumerate(spl[1]):
            if char in ["T", "C", "A", "G"]:
                positions.append((i, char))
                histo[i] += 1
        is_mutated = False
        name = ""
        for pos, char in positions:
            if pos == max_index:
                is_mutated = True
            else:
                name += sequence[pos] + str(pos + 1) + char + "_"
        if name == "":
            name = "WT"
        else:
            name = name[:-1]
        muts[name] += int(is_mutated)
        infos[name] += 1

    return muts, infos


# Define the monoexponential model
def monoexponential(t, Ymax, k):
    return Ymax * (1 - np.exp(-k * t))


def fit_monoexponential(X, Y):
    # Curve fitting for best estimate
    initial_guesses = [1.0, 0.1]
    popt, pcov = curve_fit(monoexponential, X, Y, p0=initial_guesses, maxfev=10000)
    Ymax, k = popt

    # Bootstrapping
    n_bootstrap = 1000
    k_values = []

    rng = np.random.default_rng(seed=42)  # Reproducibility

    for _ in range(n_bootstrap):
        indices = rng.integers(0, len(X), len(X))  # resample with replacement
        X_sample = X[indices]
        Y_sample = Y[indices]

        try:
            popt_sample, _ = curve_fit(
                monoexponential, X_sample, Y_sample, p0=initial_guesses, maxfev=10000
            )
            k_values.append(popt_sample[1])  # store only k
        except RuntimeError:
            continue  # skip if fit fails

    # Calculate bootstrap error for k
    k_values = np.array(k_values)
    k_std = np.std(k_values)
    return Ymax, k, k_std


@click.group()
def cli():
    pass


@cli.command()
@click.argument("data_csv", type=click.Path(exists=True))
def run_rna_map(data_csv):
    """
    Run RNA-MaP on the data.
    """
    df = pd.read_csv(data_csv)
    df = df.query("code == 'C01HP'")
    df["time"] = df["construct"].apply(get_mins_fm_dir_name)
    data = []
    print("Constructs:")
    print(df[["construct", "barcode_seq", "time"]].to_string(index=False))
    os.makedirs("data", exist_ok=True)
    for i, row in df.iterrows():
        os.system("rm -rf log output input")
        os.makedirs(f"data/{row['construct']}", exist_ok=True)
        cmd = f"rna-map -fa $SEQPATH/fastas/C01HP.fasta -fq1 demultiplexed/{row['barcode_seq']}/test_R2.fastq.gz -fq2 demultiplexed/{row['barcode_seq']}/test_R1.fastq.gz --dot-bracket $SEQPATH/rna/C01HP.csv --param-file params.yml"
        os.system(cmd)
        os.system(f"mv output data/{row['construct']}")
        df_summary = pd.read_csv(
            f"data/{row['construct']}/output/BitVector_Files/summary.csv"
        )
        sum_row = df_summary.iloc[0]
        data.append(
            {
                "construct": row["construct"],
                "time": row["time"],
                "reads": sum_row["reads"],
                "aligned": sum_row["aligned"],
            }
        )
    df_sum = pd.DataFrame(data)
    print(df_sum)


@cli.command()
def get_mutation_fractions():
    dirs = glob.glob("data/*")
    sequence = "GGAAGATCGAGTAGATCAAAGGAGGCTGACCGACCCCCCGAGCTTCGGCTCGGGGACAACTAGACATACAGTATCTTCGGATACTGAGCCTCCACAAAGAAACAACAACAACAAC"
    dfs = []
    for dir in dirs:
        path = f"{dir}/output/BitVector_Files/mtr1_mut_lib_wt_bitvectors.txt"
        muts, infos = process_mutations(sequence, path)
        keys = muts.keys()
        data = []
        for k in keys:
            data.append([k, muts[k], infos[k], muts[k] / infos[k]])
        df = pd.DataFrame(
            data, columns=["mut", "mut_count", "info_count", "mut_fraction"]
        )
        df = df.sort_values(by="mut_fraction", ascending=False)
        df.to_csv(f"{dir}/mut_fractions.csv", index=False)
        df["time"] = get_mins_fm_dir_name(dir)
        dfs.append(df)
    df = pd.concat(dfs)
    df.to_csv("all_mut_fractions.csv", index=False)


@cli.command()
@cli.option("--min-info-count", type=int, default=1000)
@cli.option("--plot", is_flag=True)
def fit_mut_fractions(min_info_count, plot):
    os.makedirs("plots", exist_ok=True)
    df = pd.read_csv("all_mut_fractions.csv")
    df.query(f"info_count >= {min_info_count}", inplace=True)
    df.sort_values(by="time", inplace=True)
    data = []
    for mut, g in df.groupby("mut"):
        g.sort_values(by="time", inplace=True)
        if len(g) < 5:
            continue
        Ymax, k, k_std = fit_monoexponential(g["time"].values, g["mut_fraction"].values)
        if plot:
            fig, ax = plt.subplots()
            ax.scatter(g["time"].values, g["mut_fraction"].values)
            ax.set_xscale("symlog")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Mutational Fraction")
            ax.set_title(f"{mut} Ymax: {Ymax:.2f} k: {k:.5f} ± {k_std:.5f}")
            t_fit = np.linspace(min(g["time"].values), max(g["time"].values), 1000)
            Y_fit = monoexponential(t_fit, Ymax, k)
            ax.plot(t_fit, Y_fit, color="black")
            fig.savefig(f"plots/{mut}.png", dpi=200)
            plt.close(fig)
        data.append([mut, Ymax, k, k_std])
    df = pd.DataFrame(data, columns=["mut", "y_max", "k", "k_std"])
    df.sort_values(by="k", ascending=False, inplace=True)
    df.to_csv("mut_kinetics.csv", index=False)


if __name__ == "__main__":
    cli()
