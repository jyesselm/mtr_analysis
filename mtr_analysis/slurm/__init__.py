"""
SLURM job management utilities.

This package provides tools for creating and submitting SLURM jobs
to parallelize analysis tasks across HPC clusters.
"""

from mtr_analysis.slurm.job import SlurmConfig, SlurmJob
from mtr_analysis.slurm.runner import SlurmRunner
from mtr_analysis.slurm.templates import create_mutation_job, create_rna_map_job

__all__ = [
    "SlurmJob",
    "SlurmConfig",
    "SlurmRunner",
    "create_rna_map_job",
    "create_mutation_job",
]
