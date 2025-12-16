"""
Configuration management for MTR analysis.

This module provides YAML-based configuration loading with validation
for paths, options, and parameter constraints.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


class ConfigError(Exception):
    """Raised when configuration validation fails."""

    pass


@dataclass
class PathsConfig:
    """Configuration for input/output paths."""

    demultiplex_dir: Path
    data_dir: Path
    output_dir: Path | None = None
    plots_dir: Path | None = None

    def __post_init__(self) -> None:
        """Convert strings to Paths and validate."""
        self.demultiplex_dir = Path(self.demultiplex_dir)
        self.data_dir = Path(self.data_dir)
        if self.output_dir is not None:
            self.output_dir = Path(self.output_dir)
        if self.plots_dir is not None:
            self.plots_dir = Path(self.plots_dir)


@dataclass
class FastqConfig:
    """Configuration for FASTQ file naming."""

    read1_pattern: str = "test_R1.fastq.gz"
    read2_pattern: str = "test_R2.fastq.gz"

    def __post_init__(self) -> None:
        """Validate FASTQ patterns."""
        if not self.read1_pattern.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
            raise ConfigError(
                f"Invalid read1_pattern '{self.read1_pattern}': "
                "must end with .fastq, .fastq.gz, .fq, or .fq.gz"
            )
        if not self.read2_pattern.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
            raise ConfigError(
                f"Invalid read2_pattern '{self.read2_pattern}': "
                "must end with .fastq, .fastq.gz, .fq, or .fq.gz"
            )


@dataclass
class SequenceConfig:
    """Configuration for sequence analysis."""

    reference_sequence: str
    construct_filter: str = "C01HP"
    target_position: int = 86

    def __post_init__(self) -> None:
        """Validate sequence configuration."""
        # Remove whitespace from sequence (YAML multiline strings may add spaces)
        self.reference_sequence = "".join(self.reference_sequence.split())

        valid_nucleotides = set("ACGTU")
        invalid_chars = set(self.reference_sequence.upper()) - valid_nucleotides
        if invalid_chars:
            raise ConfigError(
                f"Invalid nucleotides in reference_sequence: {invalid_chars}"
            )
        if self.target_position < 0:
            raise ConfigError(
                f"target_position must be >= 0, got {self.target_position}"
            )
        if self.target_position >= len(self.reference_sequence):
            raise ConfigError(
                f"target_position ({self.target_position}) exceeds sequence length "
                f"({len(self.reference_sequence)})"
            )


@dataclass
class MutationConfig:
    """Configuration for mutation analysis."""

    mutation_count_filter: int = 3
    min_info_count: int = 1000
    min_data_points: int = 5

    def __post_init__(self) -> None:
        """Validate mutation configuration."""
        if self.mutation_count_filter < 0:
            raise ConfigError(
                f"mutation_count_filter must be >= 0, got {self.mutation_count_filter}"
            )
        if self.min_info_count < 1:
            raise ConfigError(
                f"min_info_count must be >= 1, got {self.min_info_count}"
            )
        if self.min_data_points < 2:
            raise ConfigError(
                f"min_data_points must be >= 2, got {self.min_data_points}"
            )


@dataclass
class FittingConfig:
    """Configuration for curve fitting."""

    n_bootstrap: int = 1000
    random_seed: int = 42
    max_iterations: int = 10000
    initial_y_max: float = 1.0
    initial_k: float = 0.1

    def __post_init__(self) -> None:
        """Validate fitting configuration."""
        if self.n_bootstrap < 1:
            raise ConfigError(f"n_bootstrap must be >= 1, got {self.n_bootstrap}")
        if self.max_iterations < 100:
            raise ConfigError(
                f"max_iterations must be >= 100, got {self.max_iterations}"
            )
        if self.initial_y_max <= 0:
            raise ConfigError(
                f"initial_y_max must be > 0, got {self.initial_y_max}"
            )
        if self.initial_k <= 0:
            raise ConfigError(f"initial_k must be > 0, got {self.initial_k}")


@dataclass
class OutputConfig:
    """Configuration for output files."""

    mutation_fractions_file: str = "all_mut_fractions.csv"
    kinetics_file: str = "mut_kinetics.csv"
    generate_plots: bool = False

    def __post_init__(self) -> None:
        """Validate output configuration."""
        if not self.mutation_fractions_file.endswith(".csv"):
            raise ConfigError(
                f"mutation_fractions_file must end with .csv, "
                f"got '{self.mutation_fractions_file}'"
            )
        if not self.kinetics_file.endswith(".csv"):
            raise ConfigError(
                f"kinetics_file must end with .csv, got '{self.kinetics_file}'"
            )


VALID_SLURM_EMAIL_TYPES = frozenset({"NONE", "BEGIN", "END", "FAIL", "ALL"})
VALID_SLURM_TIME_PATTERN = r"^\d{1,3}:\d{2}:\d{2}$"
VALID_SLURM_MEMORY_PATTERN = r"^\d+[KMGT]?$"


@dataclass
class SlurmConfig:
    """Configuration for SLURM job submission."""

    enabled: bool = False
    time: str = "01:00:00"
    memory: str = "4G"
    cpus: int = 1
    job_name_prefix: str = "mtr"
    email: str | None = None
    email_type: str = "FAIL"
    script_dir: Path | None = None
    log_dir: Path | None = None

    # Extra commands to run at the beginning of each job (e.g., module load, conda activate)
    extra_commands: str | None = None

    # Stage-specific time overrides
    rna_map_time: str = "02:00:00"
    rna_map_memory: str = "8G"
    mutations_time: str = "00:30:00"
    mutations_memory: str = "4G"
    aggregation_time: str = "00:15:00"
    aggregation_memory: str = "2G"
    fitting_time: str = "00:30:00"
    fitting_memory: str = "4G"

    def __post_init__(self) -> None:
        """Convert paths and validate configuration."""
        import re

        if self.script_dir is not None:
            self.script_dir = Path(self.script_dir)
        if self.log_dir is not None:
            self.log_dir = Path(self.log_dir)

        # Validate email type
        if self.email_type.upper() not in VALID_SLURM_EMAIL_TYPES:
            raise ConfigError(
                f"Invalid email_type '{self.email_type}'. "
                f"Must be one of: {', '.join(sorted(VALID_SLURM_EMAIL_TYPES))}"
            )

        # Validate time formats
        time_fields = [
            ("time", self.time),
            ("rna_map_time", self.rna_map_time),
            ("mutations_time", self.mutations_time),
            ("aggregation_time", self.aggregation_time),
            ("fitting_time", self.fitting_time),
        ]
        for name, value in time_fields:
            if not re.match(VALID_SLURM_TIME_PATTERN, value):
                raise ConfigError(
                    f"Invalid {name} '{value}'. Must be in format HH:MM:SS or H:MM:SS"
                )

        # Validate memory formats
        memory_fields = [
            ("memory", self.memory),
            ("rna_map_memory", self.rna_map_memory),
            ("mutations_memory", self.mutations_memory),
            ("aggregation_memory", self.aggregation_memory),
            ("fitting_memory", self.fitting_memory),
        ]
        for name, value in memory_fields:
            if not re.match(VALID_SLURM_MEMORY_PATTERN, value, re.IGNORECASE):
                raise ConfigError(
                    f"Invalid {name} '{value}'. Must be a number optionally "
                    "followed by K, M, G, or T (e.g., '4G', '500M')"
                )

        # Validate cpus
        if self.cpus < 1:
            raise ConfigError(f"cpus must be >= 1, got {self.cpus}")


# Default reference sequence
DEFAULT_SEQUENCE = (
    "GGAAGATCGAGTAGATCAAAGGAGGCTGACCGACCCCCCGAGCTTCGGCTCGGGGACAACTA"
    "GACATACAGTATCTTCGGATACTGAGCCTCCACAAAGAAACAACAACAACAAC"
)


@dataclass
class Config:
    """Main configuration container for MTR analysis."""

    paths: PathsConfig
    fastq: FastqConfig = field(default_factory=FastqConfig)
    sequence: SequenceConfig = field(
        default_factory=lambda: SequenceConfig(reference_sequence=DEFAULT_SEQUENCE)
    )
    mutation: MutationConfig = field(default_factory=MutationConfig)
    fitting: FittingConfig = field(default_factory=FittingConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    slurm: SlurmConfig = field(default_factory=SlurmConfig)


def _validate_absolute_path(path: Path, field_name: str) -> None:
    """
    Validate that a path is absolute.

    Args:
        path: Path to validate.
        field_name: Name of the field for error messages.

    Raises:
        ConfigError: If path is not absolute.
    """
    if not path.is_absolute():
        raise ConfigError(
            f"Path for '{field_name}' must be absolute. "
            f"Got relative path: '{path}'. "
            f"Please use an absolute path like '/home/user/data' or expand with ~"
        )


def _expand_path(path_str: str) -> Path:
    """
    Expand user home directory and return absolute path.

    Args:
        path_str: Path string, possibly with ~ for home directory.

    Returns:
        Expanded absolute Path.
    """
    return Path(path_str).expanduser().resolve()


def validate_config(config: Config) -> None:
    """
    Validate the entire configuration.

    Performs validation that requires cross-field checks or
    validation that couldn't be done in __post_init__.

    Args:
        config: Configuration to validate.

    Raises:
        ConfigError: If validation fails.
    """
    # Validate required absolute paths
    _validate_absolute_path(config.paths.demultiplex_dir, "paths.demultiplex_dir")
    _validate_absolute_path(config.paths.data_dir, "paths.data_dir")

    if config.paths.output_dir is not None:
        _validate_absolute_path(config.paths.output_dir, "paths.output_dir")
    if config.paths.plots_dir is not None:
        _validate_absolute_path(config.paths.plots_dir, "paths.plots_dir")

    # Validate SLURM paths if SLURM is enabled
    if config.slurm.enabled:
        if config.slurm.script_dir is not None:
            _validate_absolute_path(config.slurm.script_dir, "slurm.script_dir")
        if config.slurm.log_dir is not None:
            _validate_absolute_path(config.slurm.log_dir, "slurm.log_dir")

    # Warn if demultiplex_dir doesn't exist
    if not config.paths.demultiplex_dir.exists():
        raise ConfigError(
            f"demultiplex_dir does not exist: {config.paths.demultiplex_dir}"
        )


def load_config(path: str | Path) -> Config:
    """
    Load configuration from a YAML file.

    Args:
        path: Path to the YAML configuration file.

    Returns:
        Validated Config object.

    Raises:
        ConfigError: If the file cannot be read or validation fails.
        FileNotFoundError: If the config file doesn't exist.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")

    try:
        with open(path) as f:
            data = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ConfigError(f"Invalid YAML in configuration file: {e}")

    if data is None:
        raise ConfigError("Configuration file is empty")

    return _parse_config(data)


def _parse_config(data: dict[str, Any]) -> Config:
    """
    Parse configuration dictionary into Config object.

    Args:
        data: Dictionary from YAML parsing.

    Returns:
        Validated Config object.

    Raises:
        ConfigError: If required fields are missing or validation fails.
    """
    # Validate required sections
    if "paths" not in data:
        raise ConfigError("Missing required section: 'paths'")

    paths_data = data["paths"]
    if "demultiplex_dir" not in paths_data:
        raise ConfigError("Missing required field: 'paths.demultiplex_dir'")
    if "data_dir" not in paths_data:
        raise ConfigError("Missing required field: 'paths.data_dir'")

    # Build paths config with expanded paths
    paths_config = PathsConfig(
        demultiplex_dir=_expand_path(paths_data["demultiplex_dir"]),
        data_dir=_expand_path(paths_data["data_dir"]),
        output_dir=_expand_path(paths_data["output_dir"])
        if paths_data.get("output_dir")
        else None,
        plots_dir=_expand_path(paths_data["plots_dir"])
        if paths_data.get("plots_dir")
        else None,
    )

    # Build fastq config
    fastq_data = data.get("fastq", {})
    fastq_config = FastqConfig(
        read1_pattern=fastq_data.get("read1_pattern", "test_R1.fastq.gz"),
        read2_pattern=fastq_data.get("read2_pattern", "test_R2.fastq.gz"),
    )

    # Build sequence config
    seq_data = data.get("sequence", {})
    sequence_config = SequenceConfig(
        reference_sequence=seq_data.get("reference_sequence", DEFAULT_SEQUENCE),
        construct_filter=seq_data.get("construct_filter", "C01HP"),
        target_position=seq_data.get("target_position", 86),
    )

    # Build mutation config
    mut_data = data.get("mutation", {})
    mutation_config = MutationConfig(
        mutation_count_filter=mut_data.get("mutation_count_filter", 3),
        min_info_count=mut_data.get("min_info_count", 1000),
        min_data_points=mut_data.get("min_data_points", 5),
    )

    # Build fitting config
    fit_data = data.get("fitting", {})
    fitting_config = FittingConfig(
        n_bootstrap=fit_data.get("n_bootstrap", 1000),
        random_seed=fit_data.get("random_seed", 42),
        max_iterations=fit_data.get("max_iterations", 10000),
        initial_y_max=fit_data.get("initial_y_max", 1.0),
        initial_k=fit_data.get("initial_k", 0.1),
    )

    # Build output config
    out_data = data.get("output", {})
    output_config = OutputConfig(
        mutation_fractions_file=out_data.get(
            "mutation_fractions_file", "all_mut_fractions.csv"
        ),
        kinetics_file=out_data.get("kinetics_file", "mut_kinetics.csv"),
        generate_plots=out_data.get("generate_plots", False),
    )

    # Build SLURM config
    slurm_data = data.get("slurm", {})
    slurm_config = SlurmConfig(
        enabled=slurm_data.get("enabled", False),
        time=slurm_data.get("time", "01:00:00"),
        memory=slurm_data.get("memory", "4G"),
        cpus=slurm_data.get("cpus", 1),
        job_name_prefix=slurm_data.get("job_name_prefix", "mtr"),
        email=slurm_data.get("email"),
        email_type=slurm_data.get("email_type", "FAIL"),
        script_dir=_expand_path(slurm_data["script_dir"])
        if slurm_data.get("script_dir")
        else None,
        log_dir=_expand_path(slurm_data["log_dir"])
        if slurm_data.get("log_dir")
        else None,
        extra_commands=slurm_data.get("extra_commands"),
        rna_map_time=slurm_data.get("rna_map_time", "02:00:00"),
        rna_map_memory=slurm_data.get("rna_map_memory", "8G"),
        mutations_time=slurm_data.get("mutations_time", "00:30:00"),
        mutations_memory=slurm_data.get("mutations_memory", "4G"),
        aggregation_time=slurm_data.get("aggregation_time", "00:15:00"),
        aggregation_memory=slurm_data.get("aggregation_memory", "2G"),
        fitting_time=slurm_data.get("fitting_time", "00:30:00"),
        fitting_memory=slurm_data.get("fitting_memory", "4G"),
    )

    # Create and validate config
    config = Config(
        paths=paths_config,
        fastq=fastq_config,
        sequence=sequence_config,
        mutation=mutation_config,
        fitting=fitting_config,
        output=output_config,
        slurm=slurm_config,
    )

    validate_config(config)
    return config


def create_example_config(output_path: str | Path) -> None:
    """
    Create an example configuration file with documentation.

    Args:
        output_path: Where to write the example config.
    """
    example = '''# MTR Analysis Configuration File
# All paths must be absolute (start with / or ~)

# =============================================================================
# PATHS - Required section
# =============================================================================
paths:
  # Directory containing demultiplexed FASTQ files
  # Structure: demultiplex_dir/{barcode_seq}/read1.fastq.gz
  demultiplex_dir: /path/to/demultiplexed

  # Output directory for RNA-MaP results and analysis
  data_dir: /path/to/output/data

  # Optional: Custom output directory for aggregated results
  # If not specified, uses current working directory
  # output_dir: /path/to/output

  # Optional: Directory for kinetics plots
  # If not specified, uses 'plots' subdirectory of output_dir
  # plots_dir: /path/to/plots

# =============================================================================
# FASTQ - File naming patterns
# =============================================================================
fastq:
  # Pattern for Read 1 FASTQ files within barcode directories
  read1_pattern: test_R1.fastq.gz

  # Pattern for Read 2 FASTQ files within barcode directories
  read2_pattern: test_R2.fastq.gz

# =============================================================================
# SEQUENCE - Reference sequence and analysis settings
# =============================================================================
sequence:
  # Reference RNA sequence for mutation analysis
  # Default is the MTR ribozyme sequence
  reference_sequence: >-
    GGAAGATCGAGTAGATCAAAGGAGGCTGACCGACCCCCCGAGCTTCGGCTCGGGGACAACTA
    GACATACAGTATCTTCGGATACTGAGCCTCCACAAAGAAACAACAACAACAAC

  # Code to filter constructs in the data CSV (e.g., "C01HP")
  construct_filter: C01HP

  # 0-indexed position to analyze for mutations
  target_position: 86

# =============================================================================
# MUTATION - Mutation filtering parameters
# =============================================================================
mutation:
  # Maximum mutations allowed per read (reads with more are filtered out)
  mutation_count_filter: 3

  # Minimum read count required for a mutation to be included in fitting
  min_info_count: 1000

  # Minimum number of time points required for curve fitting
  min_data_points: 5

# =============================================================================
# FITTING - Curve fitting parameters
# =============================================================================
fitting:
  # Number of bootstrap iterations for error estimation
  n_bootstrap: 1000

  # Random seed for reproducibility
  random_seed: 42

  # Maximum iterations for curve_fit optimization
  max_iterations: 10000

  # Initial guess for y_max parameter (maximum asymptotic value)
  initial_y_max: 1.0

  # Initial guess for k parameter (rate constant)
  initial_k: 0.1

# =============================================================================
# OUTPUT - Output file configuration
# =============================================================================
output:
  # Filename for aggregated mutation fractions
  mutation_fractions_file: all_mut_fractions.csv

  # Filename for kinetics results
  kinetics_file: mut_kinetics.csv

  # Whether to generate kinetics plots
  generate_plots: false

# =============================================================================
# SLURM - HPC cluster configuration (optional)
# =============================================================================
slurm:
  # Enable SLURM job submission
  enabled: false

  # Default time limit (HH:MM:SS)
  time: "01:00:00"

  # Default memory allocation
  memory: 4G

  # CPUs per task
  cpus: 1

  # Prefix for job names
  job_name_prefix: mtr

  # Email for job notifications (optional)
  # email: user@example.com

  # When to send email notifications
  # Valid values: NONE, BEGIN, END, FAIL, ALL
  email_type: FAIL

  # Directory for SLURM scripts (optional, must be absolute if specified)
  # script_dir: /path/to/slurm_scripts

  # Directory for SLURM logs (optional, must be absolute if specified)
  # log_dir: /path/to/slurm_logs

  # Extra commands to run at the start of each SLURM job
  # Use this for environment setup like module loads and conda activation
  # Commands are separated by newlines in the YAML string
  # extra_commands: |
  #   module load anaconda3
  #   conda activate mtr_env

  # Stage-specific resource allocations
  rna_map_time: "02:00:00"
  rna_map_memory: 8G
  mutations_time: "00:30:00"
  mutations_memory: 4G
  aggregation_time: "00:15:00"
  aggregation_memory: 2G
  fitting_time: "00:30:00"
  fitting_memory: 4G
'''
    Path(output_path).write_text(example)
