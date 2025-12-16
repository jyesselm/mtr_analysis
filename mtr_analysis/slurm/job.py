"""
SLURM job definition and configuration.

This module defines the data structures and generation utilities
for SLURM batch job scripts.
"""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class SlurmConfig:
    """Configuration for SLURM job submission."""

    time: str = "01:00:00"
    memory: str = "4G"
    cpus: int = 1
    output_dir: Path = field(default_factory=lambda: Path("slurm_logs"))
    job_name_prefix: str = "mtr"
    email: str | None = None
    email_type: str = "FAIL"
    extra_commands: str | None = None


@dataclass
class SlurmJob:
    """Represents a single SLURM job."""

    name: str
    commands: list[str]
    config: SlurmConfig
    dependencies: list[str] = field(default_factory=list)

    def generate_script(self) -> str:
        """
        Generate the SLURM batch script content.

        Returns:
            Complete SLURM script as a string.
        """
        lines = self._generate_header()
        lines.extend(self._generate_environment())
        lines.extend(self.commands)
        return "\n".join(lines)

    def _generate_header(self) -> list[str]:
        """Generate SLURM header directives."""
        header = [
            "#!/bin/bash",
            f"#SBATCH --job-name={self.config.job_name_prefix}_{self.name}",
            f"#SBATCH --time={self.config.time}",
            f"#SBATCH --mem={self.config.memory}",
            f"#SBATCH --cpus-per-task={self.config.cpus}",
            f"#SBATCH --output={self.config.output_dir}/{self.name}_%j.out",
            f"#SBATCH --error={self.config.output_dir}/{self.name}_%j.err",
        ]
        if self.config.email:
            header.append(f"#SBATCH --mail-user={self.config.email}")
            header.append(f"#SBATCH --mail-type={self.config.email_type}")
        return header

    def _generate_environment(self) -> list[str]:
        """Generate environment setup commands."""
        lines = [
            "",
            "# Exit on error",
            "set -e",
            "",
            "# Print job info",
            'echo "Job started at $(date)"',
            'echo "Running on node: $(hostname)"',
            "",
        ]
        # Add extra commands if provided (e.g., module load, conda activate)
        if self.config.extra_commands:
            lines.append("# Environment setup")
            for cmd in self.config.extra_commands.strip().split("\n"):
                cmd = cmd.strip()
                if cmd:
                    lines.append(cmd)
            lines.append("")
        return lines

    def write_script(self, output_dir: Path) -> Path:
        """
        Write the job script to a file.

        Args:
            output_dir: Directory to write the script.

        Returns:
            Path to the written script file.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        script_path = output_dir / f"{self.name}.sh"
        script_path.write_text(self.generate_script())
        return script_path


def create_job(
    name: str,
    commands: list[str],
    config: SlurmConfig | None = None,
) -> SlurmJob:
    """
    Create a SLURM job with the given configuration.

    Args:
        name: Job name identifier.
        commands: Shell commands to execute.
        config: SLURM configuration (uses defaults if None).

    Returns:
        Configured SlurmJob instance.
    """
    if config is None:
        config = SlurmConfig()
    return SlurmJob(name=name, commands=commands, config=config)
