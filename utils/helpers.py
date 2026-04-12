#!/usr/bin/env python3

from pathlib import Path
from typing import List
import shutil
import subprocess


# general helper functions


def vprint(verbose: bool, *args, **kwargs) -> None:
    """Print only when verbose mode is enabled.

    Args:
        verbose: Whether printing is enabled.
        *args: Positional arguments passed to print().
        **kwargs: Keyword arguments passed to print().
    """
    if verbose:
        print(*args, **kwargs)


def require_file(path: Path, label: str) -> None:
    """Check that a required file exists.

    Args:
        path: Path to the required file.
        label: Human-readable label used in the error message.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def require_executable(name: str, label: str) -> None:
    """Check that a required executable is available on PATH.

    Args:
        name: Name of the executable to look up.
        label: Human-readable label used in the error message.

    Raises:
        FileNotFoundError: If the executable is not found on PATH.
    """
    if shutil.which(name) is None:
        raise FileNotFoundError(f"{label} not found on PATH: {name}")


def run_command(cmd: List[str], verbose: bool = True) -> None:
    """Run a shell command and raise an error if it fails.

    Args:
        cmd: Command and arguments to execute.
        verbose: Whether to print the command before running it.

    Raises:
        subprocess.CalledProcessError: If the command exits with a non-zero
            status.
    """
    vprint(verbose, "Running:")
    vprint(verbose, "  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def ensure_dir(path: Path) -> None:
    """Create a directory if it does not already exist.

    Args:
        path: Directory path to create.
    """
    path.mkdir(parents=True, exist_ok=True)


def remove_empty_tmp_dirs(results_dir: Path, verbose: bool = True) -> None:
    """Remove empty directories named ``tmp`` under a results directory.

    Args:
        results_dir: Top-level results directory to scan.
        verbose: Whether to print routine progress messages.

    Raises:
        FileNotFoundError: If ``results_dir`` does not exist.
        NotADirectoryError: If ``results_dir`` is not a directory.
    """
    require_file(results_dir, "Results directory")
    if not results_dir.is_dir():
        raise NotADirectoryError(f"Results directory is not a directory: {results_dir}")

    removed_any = False

    for path in sorted(results_dir.rglob("tmp"), key=lambda p: len(p.parts), reverse=True):
        if not path.is_dir():
            continue

        try:
            next(path.iterdir())
        except StopIteration:
            path.rmdir()
            removed_any = True
            print(f"Removed empty tmp directory: {path}")
        else:
            vprint(verbose, f"Keeping non-empty tmp directory: {path}")

    if not removed_any:
        vprint(verbose, f"No empty tmp directories found under: {results_dir}")


# HALPER related helper functions


def should_use_existing_halper_output(
    ortholog_bed: Path,
    verbose: bool = True,
) -> bool:
    """Ask whether to reuse an existing HALPER output file.

    If the HALPER ortholog BED file already exists, the user is prompted to
    decide whether to reuse it. Pressing Enter defaults to yes.

    Args:
        ortholog_bed: Path to the expected HALPER ortholog output file.
        verbose: Whether to print informational messages after the prompt.

    Returns:
        True if the existing output should be reused.
        False if a new HALPER job should be submitted instead.
    """
    if not ortholog_bed.exists():
        return False

    while True:
        response = input(
            f"Existing HALPER output found:\n"
            f"  {ortholog_bed}\n"
            f"Use this existing output? [Y/n]: "
        ).strip().lower()

        if response in ("", "y", "yes"):
            vprint(verbose, f"Using existing HALPER output: {ortholog_bed}")
            return True

        if response in ("n", "no"):
            vprint(verbose, f"Will regenerate HALPER output: {ortholog_bed}")
            return False

        print("Please enter y, yes, n, no, or press Enter for yes.")


def submit_or_run_job(
    use_sbatch: bool,
    script_path: Path,
    verbose: bool = True,
) -> None:
    """Submit a batch script with sbatch or run it directly with bash.

    Args:
        use_sbatch: Whether to submit the script through SLURM.
        script_path: Path to the script to submit or run.
        verbose: Whether to print the command before running it.

    Raises:
        subprocess.CalledProcessError: If the submission or execution command
            fails.
    """
    if use_sbatch:
        run_command(["sbatch", str(script_path)], verbose=verbose)
    else:
        run_command(["bash", str(script_path)], verbose=verbose)


def submit_halper_job_or_exit(
    use_sbatch: bool,
    script_path: Path,
    verbose: bool = True,
) -> None:
    """Submit or run a HALPER job, then stop the program.

    Args:
        use_sbatch: Whether to submit through SLURM.
        script_path: Path to the job script.
        verbose: Whether to print the submission command.

    Raises:
        SystemExit: Always raised after the job is submitted or started.
    """
    submit_or_run_job(use_sbatch, script_path, verbose=verbose)
    raise SystemExit(
        "HALPER job submitted/launched. Re-run the pipeline after HALPER finishes."
    )


# bedtools related helper functions


def run_bedtools_to_file(
    cmd: List[str],
    output_path: Path,
    verbose: bool = True,
) -> None:
    """Run a bedtools command and write stdout to an output file safely.

    Output is first written to a temporary file and then moved into place only
    if the command succeeds.

    Args:
        cmd: Command and arguments to execute.
        output_path: Destination file for command output.
        verbose: Whether to print the command before running it.

    Raises:
        subprocess.CalledProcessError: If the BEDTools command fails.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    try:
        with open(tmp_path, "w") as fout:
            vprint(verbose, "Running:")
            vprint(verbose, "  " + " ".join(cmd))
            subprocess.run(cmd, stdout=fout, check=True)
        tmp_path.replace(output_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()