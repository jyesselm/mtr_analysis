# mtr_analysis

[![Test Example](https://github.com/jyesselm/mtr_analysis/actions/workflows/test-example.yml/badge.svg)](https://github.com/jyesselm/mtr_analysis/actions/workflows/test-example.yml)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)

Analysis scripts for processing MTR ribozyme mutation data.

## Install

```shell
git clone https://github.com/jyesselm/mtr_analysis
cd mtr_analysis
conda env create -f environment.yml
conda activate mtr-analysis
pip install -e ".[dev]"  # includes dev dependencies for testing
```

## Usage

There are two ways to run the analysis pipeline:

1. **Local execution** - Run all steps sequentially on a single machine
2. **SLURM execution** - Parallelize processing across an HPC cluster

---

## Configuration

The pipeline uses a YAML configuration file to manage all settings. This is the recommended way to configure your analysis.

### Quick Start

```shell
# Generate an example config file
mtr-analysis config init --output config.yml

# Edit config.yml to set your paths (must be absolute!)
# Then validate your configuration
mtr-analysis config validate config.yml
```

### Configuration File Structure

```yaml
# All paths must be absolute (start with / or ~)
paths:
  demultiplex_dir: /path/to/demultiplexed    # Required: FASTQ input directory
  data_dir: /path/to/output/data              # Required: RNA-MaP output directory
  output_dir: /path/to/output                 # Optional: aggregated results
  plots_dir: /path/to/plots                   # Optional: kinetics plots

fastq:
  read1_pattern: test_R1.fastq.gz             # FASTQ filename pattern
  read2_pattern: test_R2.fastq.gz

sequence:
  reference_sequence: GGAAGATCGAG...          # Reference RNA sequence
  construct_filter: C01HP                      # Filter constructs by code
  target_position: 86                          # 0-indexed mutation position

mutation:
  mutation_count_filter: 3                     # Max mutations per read
  min_info_count: 1000                         # Min reads for fitting
  min_data_points: 5                           # Min time points for fitting

fitting:
  n_bootstrap: 1000                            # Bootstrap iterations
  random_seed: 42                              # For reproducibility

output:
  mutation_fractions_file: all_mut_fractions.csv
  kinetics_file: mut_kinetics.csv
  generate_plots: false

slurm:
  enabled: false                               # Enable for HPC clusters
  time: "01:00:00"                             # Default time limit
  memory: 4G                                   # Default memory
  email: user@example.com                      # Job notifications
  email_type: FAIL                             # NONE, BEGIN, END, FAIL, ALL
```

### Config Commands

```shell
# Generate example config with documentation
mtr-analysis config init --output config.yml

# Validate configuration (checks paths, options, etc.)
mtr-analysis config validate config.yml

# Display all parsed configuration values
mtr-analysis config show config.yml
```

### Validation

The config system validates:

| Check | Description |
|-------|-------------|
| Absolute paths | All paths must be absolute (not relative) |
| Path existence | `demultiplex_dir` must exist |
| Valid nucleotides | Sequence must contain only A, C, G, T, U |
| SLURM time format | Must be `HH:MM:SS` |
| SLURM memory format | Must be number + unit (e.g., `4G`, `500M`) |
| SLURM email type | Must be `NONE`, `BEGIN`, `END`, `FAIL`, or `ALL` |
| Numeric bounds | Positive values where required |

Example error messages:
```
Error: Path for 'paths.demultiplex_dir' must be absolute. Got relative path: 'data'
Error: Invalid email_type 'INVALID'. Must be one of: ALL, BEGIN, END, FAIL, NONE
Error: Invalid time '1hour'. Must be in format HH:MM:SS
```

---

## Local Execution

### Step 1: Run RNA-MaP to get bitvector files

This is the longest step and will take > 1 hour.

```shell
# Download the demultiplexed data and data.csv from $NRDSTOR/run_name
mtr-analysis run-rna-map data.csv
```

Expected output:
```
            construct  barcode_seq  time
  mtr1_mut_lib_t15min CTGCGTGCAAAC    15
  mtr1_mut_lib_t60min GCAAATGTGCTA    60
 mtr1_mut_lib_t180min AAGGACCACTGG   180
 mtr1_mut_lib_t420min CGGGCACGGCGG   420
mtr1_mut_lib_t1440min CTAGCAATGTGA  1440

# Summary at the end:
               construct  time    reads  aligned
0    mtr1_mut_lib_t15min    15  8372027    98.38
1    mtr1_mut_lib_t60min    60  7073223    98.88
2   mtr1_mut_lib_t180min   180  5069385    98.63
3   mtr1_mut_lib_t420min   420  9167912    98.58
4  mtr1_mut_lib_t1440min  1440      654    97.25
```

### Step 2: Get mutation fractions for each construct

```shell
mtr-analysis get-mutation-fractions
```

This creates:
- `data/<construct>/mut_fractions.csv` for each construct
- `all_mut_fractions.csv` combining all constructs

### Step 3: Fit mutation fractions to monoexponential model

```shell
mtr-analysis fit-mut-fractions --plot
```

This creates:
- `mut_kinetics.csv` with y_max, k, and k_std for each mutation
- `plots/<mutation>.png` with fitted curves (if `--plot` is used)

---

## SLURM Execution (HPC Clusters)

For large datasets, use SLURM to parallelize processing across cluster nodes.

### Quick Start: Full Pipeline

Generate and submit the complete pipeline with one command:

```shell
mtr-analysis slurm setup-full-pipeline data.csv --submit
```

This creates jobs with proper dependencies:
1. RNA-MaP jobs (parallel, one per construct)
2. Mutation analysis jobs (parallel, after RNA-MaP completes)
3. Aggregation job (after all mutation jobs complete)
4. Fitting job (after aggregation completes)

### Step-by-Step SLURM Setup

For more control, set up each stage separately:

#### Step 1: Generate RNA-MaP jobs

```shell
# Generate job scripts (one per construct)
mtr-analysis slurm setup-rna-map data.csv --script-dir slurm_scripts

# Preview what would be submitted
mtr-analysis slurm setup-rna-map data.csv --dry-run

# Generate and submit immediately
mtr-analysis slurm setup-rna-map data.csv --submit
```

#### Step 2: Generate mutation analysis jobs

```shell
# After RNA-MaP jobs complete, generate mutation jobs
mtr-analysis slurm setup-mutations --script-dir slurm_scripts --submit
```

#### Step 3: Generate aggregation job

```shell
mtr-analysis slurm setup-aggregation --script-dir slurm_scripts --submit
```

#### Step 4: Generate fitting job

```shell
mtr-analysis slurm setup-fitting --script-dir slurm_scripts --plot --submit
```

### SLURM Options

All SLURM commands support these options:

| Option | Default | Description |
|--------|---------|-------------|
| `--script-dir` | `slurm_scripts` | Directory for generated scripts |
| `--time` | varies | Time limit per job |
| `--memory` | varies | Memory per job |
| `--dry-run` | False | Preview without submitting |
| `--submit` | False | Submit jobs after generating |

### Manual Submission with Dependencies

If you prefer manual control, generate scripts and submit with dependencies:

```shell
# Generate all scripts
mtr-analysis slurm setup-full-pipeline data.csv

# Submit RNA-MaP jobs and capture job IDs
RNA_JOBS=$(for f in slurm_scripts/rna_map_*.sh; do
    sbatch $f | awk '{print $4}'
done | paste -sd,)

# Submit mutation jobs after RNA-MaP completes
for f in slurm_scripts/mutations_*.sh; do
    sbatch --dependency=afterok:$RNA_JOBS $f
done

# Continue with aggregation and fitting...
```

### Example: Full Pipeline on HPC

```shell
# 1. Transfer data to cluster
scp -r demultiplexed/ data.csv user@cluster:/scratch/project/

# 2. SSH to cluster and set up environment
ssh user@cluster
cd /scratch/project
module load conda
conda activate mtr-analysis

# 3. Generate and submit full pipeline
mtr-analysis slurm setup-full-pipeline data.csv --submit

# 4. Monitor jobs
squeue -u $USER

# 5. After completion, retrieve results
# Results will be in: mut_kinetics.csv, plots/
```

---

## CLI Reference

### Main Commands

| Command | Description |
|---------|-------------|
| `run-rna-map` | Run RNA-MaP on all constructs |
| `run-single-rna-map` | Run RNA-MaP for a single barcode |
| `get-mutation-fractions` | Compute mutation fractions for all directories |
| `process-single-dir` | Process mutations for a single directory |
| `aggregate-mutations` | Aggregate mutation results from all directories |
| `fit-mut-fractions` | Fit monoexponential curves to mutation data |

### Config Commands

| Command | Description |
|---------|-------------|
| `config init` | Generate an example configuration file |
| `config validate` | Validate a configuration file |
| `config show` | Display all parsed configuration values |

### SLURM Commands

| Command | Description |
|---------|-------------|
| `slurm setup-rna-map` | Generate SLURM jobs for RNA-MaP (parallel) |
| `slurm setup-mutations` | Generate SLURM jobs for mutation analysis |
| `slurm setup-aggregation` | Generate SLURM job for aggregation |
| `slurm setup-fitting` | Generate SLURM job for curve fitting |
| `slurm setup-full-pipeline` | Generate complete pipeline with dependencies |

---

## Development

```shell
# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run linting
ruff check mtr_analysis/ tests/

# Run type checking
mypy mtr_analysis/

# Format code
ruff format mtr_analysis/ tests/
```
