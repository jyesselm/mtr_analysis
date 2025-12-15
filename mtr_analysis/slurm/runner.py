"""
SLURM job submission and execution management.

This module handles the submission, monitoring, and coordination
of SLURM jobs for the MTR analysis pipeline.
"""

import subprocess
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class JobStatus:
    """Status information for a submitted job."""

    job_id: str
    name: str
    state: str = "PENDING"


@dataclass
class SlurmRunner:
    """Manages SLURM job submission and tracking."""

    script_dir: Path
    dry_run: bool = False
    submitted_jobs: dict[str, JobStatus] = field(default_factory=dict)

    def submit_job(self, script_path: Path, dependency: str | None = None) -> str:
        """
        Submit a job script to SLURM.

        Args:
            script_path: Path to the SLURM script.
            dependency: Optional job ID to depend on.

        Returns:
            Job ID of the submitted job.
        """
        cmd = self._build_submit_command(script_path, dependency)
        if self.dry_run:
            return self._handle_dry_run(script_path, cmd)
        return self._execute_submit(script_path, cmd)

    def _build_submit_command(
        self, script_path: Path, dependency: str | None
    ) -> list[str]:
        """Build sbatch command with optional dependency."""
        cmd = ["sbatch"]
        if dependency:
            cmd.extend(["--dependency", f"afterok:{dependency}"])
        cmd.append(str(script_path))
        return cmd

    def _handle_dry_run(self, script_path: Path, cmd: list[str]) -> str:
        """Handle dry run mode by printing command."""
        print(f"[DRY RUN] Would submit: {' '.join(cmd)}")
        fake_id = f"dry_run_{script_path.stem}"
        self._record_job(fake_id, script_path.stem)
        return fake_id

    def _execute_submit(self, script_path: Path, cmd: list[str]) -> str:
        """Execute sbatch and parse job ID."""
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        job_id = self._parse_job_id(result.stdout)
        self._record_job(job_id, script_path.stem)
        return job_id

    def _parse_job_id(self, output: str) -> str:
        """Parse job ID from sbatch output."""
        # Output format: "Submitted batch job 123456"
        return output.strip().split()[-1]

    def _record_job(self, job_id: str, name: str) -> None:
        """Record submitted job for tracking."""
        self.submitted_jobs[job_id] = JobStatus(job_id=job_id, name=name)

    def submit_job_chain(self, script_paths: list[Path]) -> list[str]:
        """
        Submit jobs as a dependency chain.

        Each job waits for the previous one to complete.

        Args:
            script_paths: Ordered list of script paths.

        Returns:
            List of submitted job IDs.
        """
        job_ids = []
        prev_id = None
        for script_path in script_paths:
            job_id = self.submit_job(script_path, dependency=prev_id)
            job_ids.append(job_id)
            prev_id = job_id
        return job_ids

    def submit_parallel_jobs(self, script_paths: list[Path]) -> list[str]:
        """
        Submit jobs to run in parallel.

        Args:
            script_paths: List of script paths to submit.

        Returns:
            List of submitted job IDs.
        """
        return [self.submit_job(path) for path in script_paths]

    def get_job_status(self, job_id: str) -> JobStatus | None:
        """
        Get the status of a submitted job.

        Args:
            job_id: Job ID to query.

        Returns:
            JobStatus if found, None otherwise.
        """
        return self.submitted_jobs.get(job_id)

    def print_summary(self) -> None:
        """Print summary of all submitted jobs."""
        if not self.submitted_jobs:
            print("No jobs submitted.")
            return
        print(f"\nSubmitted {len(self.submitted_jobs)} jobs:")
        for job_id, status in self.submitted_jobs.items():
            print(f"  {job_id}: {status.name}")
