"""Tests for SLURM modules."""

from pathlib import Path

from mtr_analysis.slurm.job import SlurmConfig, SlurmJob, create_job
from mtr_analysis.slurm.runner import JobStatus, SlurmRunner
from mtr_analysis.slurm.templates import (
    create_aggregation_job,
    create_fitting_job,
    create_mutation_job,
    create_rna_map_job,
)


class TestSlurmConfig:
    """Tests for SlurmConfig dataclass."""

    def test_default_values(self) -> None:
        """Test default configuration values."""
        config = SlurmConfig()
        assert config.time == "01:00:00"
        assert config.memory == "4G"
        assert config.cpus == 1
        assert config.extra_commands is None

    def test_custom_values(self) -> None:
        """Test custom configuration values."""
        config = SlurmConfig(
            time="04:00:00",
            memory="16G",
            cpus=4,
        )
        assert config.time == "04:00:00"
        assert config.memory == "16G"
        assert config.cpus == 4

    def test_extra_commands(self) -> None:
        """Test extra_commands configuration."""
        config = SlurmConfig(
            extra_commands="module load anaconda3\nconda activate mtr_env",
        )
        assert config.extra_commands == "module load anaconda3\nconda activate mtr_env"


class TestSlurmJob:
    """Tests for SlurmJob class."""

    def test_generate_script(self) -> None:
        """Test generating a SLURM script."""
        config = SlurmConfig()
        job = SlurmJob(
            name="test_job",
            commands=["echo 'Hello World'"],
            config=config,
        )

        script = job.generate_script()

        assert "#!/bin/bash" in script
        assert "#SBATCH --job-name=mtr_test_job" in script
        assert "echo 'Hello World'" in script

    def test_generate_script_with_email(self) -> None:
        """Test that email settings are included when specified."""
        config = SlurmConfig(email="user@example.com", email_type="ALL")
        job = SlurmJob(name="test", commands=["pwd"], config=config)

        script = job.generate_script()

        assert "#SBATCH --mail-user=user@example.com" in script
        assert "#SBATCH --mail-type=ALL" in script

    def test_generate_script_with_extra_commands(self) -> None:
        """Test that extra_commands are included in script."""
        config = SlurmConfig(
            extra_commands="module load anaconda3\nconda activate mtr_env"
        )
        job = SlurmJob(name="test", commands=["echo 'test'"], config=config)

        script = job.generate_script()

        assert "# Environment setup" in script
        assert "module load anaconda3" in script
        assert "conda activate mtr_env" in script

    def test_generate_script_without_extra_commands(self) -> None:
        """Test that script works without extra_commands."""
        config = SlurmConfig()
        job = SlurmJob(name="test", commands=["echo 'test'"], config=config)

        script = job.generate_script()

        assert "# Environment setup" not in script
        assert "echo 'test'" in script

    def test_write_script(self, temp_dir: Path) -> None:
        """Test writing script to file."""
        config = SlurmConfig()
        job = SlurmJob(name="write_test", commands=["ls -la"], config=config)

        script_path = job.write_script(temp_dir)

        assert script_path.exists()
        assert script_path.name == "write_test.sh"
        content = script_path.read_text()
        assert "ls -la" in content


class TestCreateJob:
    """Tests for create_job factory function."""

    def test_create_with_defaults(self) -> None:
        """Test creating job with default config."""
        job = create_job(name="test", commands=["echo test"])

        assert job.name == "test"
        assert job.commands == ["echo test"]
        assert job.config.time == "01:00:00"

    def test_create_with_custom_config(self) -> None:
        """Test creating job with custom config."""
        config = SlurmConfig(time="02:00:00")
        job = create_job(name="test", commands=["pwd"], config=config)

        assert job.config.time == "02:00:00"


class TestJobStatus:
    """Tests for JobStatus dataclass."""

    def test_creation(self) -> None:
        """Test creating a JobStatus instance."""
        status = JobStatus(job_id="12345", name="test_job")
        assert status.job_id == "12345"
        assert status.name == "test_job"
        assert status.state == "PENDING"


class TestSlurmRunner:
    """Tests for SlurmRunner class."""

    def test_dry_run_submit(self, temp_dir: Path) -> None:
        """Test job submission in dry run mode."""
        runner = SlurmRunner(script_dir=temp_dir, dry_run=True)

        # Create a test script
        script_path = temp_dir / "test.sh"
        script_path.write_text("#!/bin/bash\necho test")

        job_id = runner.submit_job(script_path)

        assert "dry_run" in job_id
        assert len(runner.submitted_jobs) == 1

    def test_submit_parallel_jobs_dry_run(self, temp_dir: Path) -> None:
        """Test parallel job submission in dry run mode."""
        runner = SlurmRunner(script_dir=temp_dir, dry_run=True)

        # Create test scripts
        scripts = []
        for i in range(3):
            script_path = temp_dir / f"job_{i}.sh"
            script_path.write_text(f"#!/bin/bash\necho job {i}")
            scripts.append(script_path)

        job_ids = runner.submit_parallel_jobs(scripts)

        assert len(job_ids) == 3
        assert len(runner.submitted_jobs) == 3

    def test_submit_job_chain_dry_run(self, temp_dir: Path) -> None:
        """Test sequential job submission in dry run mode."""
        runner = SlurmRunner(script_dir=temp_dir, dry_run=True)

        # Create test scripts
        scripts = []
        for i in range(2):
            script_path = temp_dir / f"chain_{i}.sh"
            script_path.write_text(f"#!/bin/bash\necho step {i}")
            scripts.append(script_path)

        job_ids = runner.submit_job_chain(scripts)

        assert len(job_ids) == 2

    def test_get_job_status(self, temp_dir: Path) -> None:
        """Test retrieving job status."""
        runner = SlurmRunner(script_dir=temp_dir, dry_run=True)

        script_path = temp_dir / "status_test.sh"
        script_path.write_text("#!/bin/bash\necho test")

        job_id = runner.submit_job(script_path)
        status = runner.get_job_status(job_id)

        assert status is not None
        assert status.job_id == job_id

    def test_get_nonexistent_job_status(self, temp_dir: Path) -> None:
        """Test retrieving status for nonexistent job."""
        runner = SlurmRunner(script_dir=temp_dir, dry_run=True)
        status = runner.get_job_status("nonexistent")
        assert status is None


class TestSlurmTemplates:
    """Tests for SLURM job templates."""

    def test_create_rna_map_job(self) -> None:
        """Test creating RNA-MaP job."""
        job = create_rna_map_job(
            construct="test_construct",
            barcode_seq="ACGT",
            output_dir="data",
        )

        assert "rna_map_test_construct" in job.name
        script = job.generate_script()
        assert "test_construct" in script
        # Verify new command structure
        assert "mtr-analysis run single-rna-map" in script

    def test_create_rna_map_job_with_extra_commands(self) -> None:
        """Test creating RNA-MaP job with extra_commands."""
        config = SlurmConfig(
            time="02:00:00",
            memory="8G",
            extra_commands="module load anaconda3",
        )
        job = create_rna_map_job(
            construct="test_construct",
            barcode_seq="ACGT",
            output_dir="data",
            config=config,
        )

        script = job.generate_script()
        assert "module load anaconda3" in script

    def test_create_mutation_job(self) -> None:
        """Test creating mutation analysis job."""
        job = create_mutation_job(
            data_dir="data/test_dir",
            sequence="ACGT",
        )

        assert "mutations_" in job.name
        script = job.generate_script()
        assert "data/test_dir" in script
        # Verify new command structure
        assert "mtr-analysis run process-dir" in script

    def test_create_aggregation_job(self) -> None:
        """Test creating aggregation job."""
        job = create_aggregation_job(
            data_dirs=["data/dir1", "data/dir2"],
            output_file="results.csv",
        )

        assert "aggregate" in job.name
        script = job.generate_script()
        assert "results.csv" in script
        # Verify new command structure
        assert "mtr-analysis run aggregate" in script

    def test_create_fitting_job(self) -> None:
        """Test creating fitting job."""
        job = create_fitting_job(
            input_file="mutations.csv",
            output_file="kinetics.csv",
            min_info_count=500,
            generate_plots=True,
        )

        assert "fit_kinetics" in job.name
        script = job.generate_script()
        assert "--min-info-count 500" in script
        assert "--plot" in script
        # Verify new command structure
        assert "mtr-analysis run fit" in script

    def test_fitting_job_without_plots(self) -> None:
        """Test creating fitting job without plots."""
        job = create_fitting_job(
            input_file="mutations.csv",
            output_file="kinetics.csv",
            generate_plots=False,
        )

        script = job.generate_script()
        assert "--plot" not in script
