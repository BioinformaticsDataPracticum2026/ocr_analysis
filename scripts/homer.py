#!/usr/bin/env python3

from pathlib import Path
from typing import Dict, List, Optional

import utils.helpers as helpers


"""
Run HOMER motif enrichment analyses for cross-species OCR sets.

This wrapper prepares and submits one HOMER job per:
- dataset
- analysis mode

Main biological questions:
1. How many motifs occur by chance?
   - compare default HOMER background vs matched background BEDs
   - default_background vs matched_background

2. Do we expect motifs to be randomly distributed in open chromatin regions?
   - compare full-region runs (-size given) vs centered fixed-width runs (-size 200)
   - matched_background vs centered_200
"""


def homer_results_exist(output_dir: Path) -> bool:
    """Check whether a HOMER output directory already contains results.

    This is used to decide whether a HOMER run can be skipped. The check is
    intentionally simple and looks for the standard HOMER output subdirectories
    created for known or de novo motif results.

    Args:
        output_dir: Path to the expected HOMER output directory.

    Returns:
        True if the directory exists and contains either ``knownResults`` or
        ``homerResults``. False otherwise.
    """
    if not output_dir.exists():
        return False

    known_dir = output_dir / "knownResults"
    homer_dir = output_dir / "homerResults"

    return known_dir.exists() or homer_dir.exists()


def build_homer_jobs(config: dict) -> List[Dict[str, object]]:
    """Build the list of HOMER jobs to prepare and submit.

    This function defines the input BED files, matched background BED files,
    genome FASTA files, preparsed output directories, and analysis modes for
    all HOMER runs. Cross-species BED files are paired with the FASTA of the
    target-species coordinate system.

    Args:
        config: Pipeline configuration dictionary loaded from YAML.

    Returns:
        A list of dictionaries, where each dictionary describes one HOMER job.
        Each job includes the input BED, background BED, genome FASTA,
        preparsed directory, output directory, script path, log path, and
        analysis settings.
    """
    results_dir = Path(config["project"]["output_dir"])

    cross_species_dir = results_dir / "bedtools" / "cross_species_ep"
    promoter_enhancer_dir = results_dir / "bedtools" / "promoter_enhancer"
    homer_dir = results_dir / "homer"
    preparsed_root = homer_dir / "preparsed"

    species_1_genome_fasta = Path(config["annotations"]["species_1_genome_fasta"])
    species_2_genome_fasta = Path(config["annotations"]["species_2_genome_fasta"])

    datasets = [
        {
            # Human OCRs mapped into mouse, retained if open in mouse,
            # then classified as enhancer-like.
            # These BED coordinates are in mouse space, so use species_2_genome_fasta.
            # The matched background is the broader human pancreas enhancer set.
            "label": "human_open_in_mouse_enhancers",
            "input_bed": cross_species_dir / "human_open_in_mouse_enhancers.bed",
            "background_bed": promoter_enhancer_dir / "human_pancreas_enhancers.bed",
            "genome_fasta": species_2_genome_fasta,
        },
        {
            # Human OCRs mapped into mouse, retained if open in mouse,
            # then classified as promoter-like.
            # These BED coordinates are in mouse space, so use species_2_genome_fasta.
            # The matched background is the broader human pancreas promoter set.
            "label": "human_open_in_mouse_promoters",
            "input_bed": cross_species_dir / "human_open_in_mouse_promoters.bed",
            "background_bed": promoter_enhancer_dir / "human_pancreas_promoters.bed",
            "genome_fasta": species_2_genome_fasta,
        },
        {
            # Mouse OCRs mapped into human, retained if open in human,
            # then classified as enhancer-like.
            # These BED coordinates are in human space, so use species_1_genome_fasta.
            # The matched background is the broader mouse pancreas enhancer set.
            "label": "mouse_open_in_human_enhancers",
            "input_bed": cross_species_dir / "mouse_open_in_human_enhancers.bed",
            "background_bed": promoter_enhancer_dir / "mouse_pancreas_enhancers.bed",
            "genome_fasta": species_1_genome_fasta,
        },
        {
            # Mouse OCRs mapped into human, retained if open in human,
            # then classified as promoter-like.
            # These BED coordinates are in human space, so use species_1_genome_fasta.
            # The matched background is the broader mouse pancreas promoter set.
            "label": "mouse_open_in_human_promoters",
            "input_bed": cross_species_dir / "mouse_open_in_human_promoters.bed",
            "background_bed": promoter_enhancer_dir / "mouse_pancreas_promoters.bed",
            "genome_fasta": species_1_genome_fasta,
        },
    ]

    analysis_modes = [
        {
            # 1. How many motifs occur in random regions of the genome
            #    (HOMER findMotifsGenome.pl default settings)?
            "mode": "default_background",
            "size_arg": "given",
            "use_background": False,
        },
        {
            # 2. How many motifs occur in control sequences
            #    (HOMER findMotifsGenome.pl -bg [background regions])?
            "mode": "matched_background",
            "size_arg": "given",
            "use_background": True,
        },
        {
            # 3. Do motifs appear more enriched in centered fixed-width windows
            #    than across the full OCR?
            "mode": "centered_200",
            "size_arg": "200",
            "use_background": True,
        },
    ]

    jobs: List[Dict[str, object]] = []

    for dataset in datasets:
        for mode in analysis_modes:
            output_dir = homer_dir / dataset["label"] / mode["mode"]
            script_dir = output_dir / "tmp"
            log_dir = output_dir / "logs"

            jobs.append(
                {
                    "label": dataset["label"],
                    "mode": mode["mode"],
                    "input_bed": dataset["input_bed"],
                    "background_bed": dataset["background_bed"],
                    "genome_fasta": dataset["genome_fasta"],
                    "preparsed_dir": preparsed_root / dataset["genome_fasta"].stem,
                    "size_arg": mode["size_arg"],
                    "use_background": mode["use_background"],
                    "output_dir": output_dir,
                    "script_path": script_dir / f"{dataset['label']}__{mode['mode']}.sbatch.sh",
                    "log_path": log_dir / f"{dataset['label']}__{mode['mode']}.log",
                }
            )

    return jobs


def write_homer_sbatch_script(
    script_path: Path,
    log_path: Path,
    partition: str,
    time_limit: str,
    nodes: int,
    ntasks: int,
    mem: str,
    allocation_id: str,
    input_bed: Path,
    genome_fasta: Path,
    preparsed_dir: Path,
    output_dir: Path,
    size_arg: str,
    background_bed: Optional[Path],
    job_name: str,
) -> None:
    """Write an sbatch script for a single HOMER job.

    The generated script loads HOMER on the cluster, validates that
    ``findMotifsGenome.pl`` is available, and runs motif enrichment using a BED
    file, a genome FASTA file, a writable preparsed directory, and an optional
    matched background BED file.

    Args:
        script_path: Path where the sbatch script should be written.
        log_path: Path to the SLURM log file for this job.
        partition: SLURM partition name.
        time_limit: Walltime limit for the job, formatted for sbatch.
        nodes: Number of nodes to request.
        ntasks: Number of CPUs/tasks to request and pass to HOMER with ``-p``.
        mem: Memory request string for sbatch, such as ``"64000MB"``.
        allocation_id: SLURM allocation or account ID.
        input_bed: Input BED file for motif analysis.
        genome_fasta: Genome FASTA file used by HOMER to extract sequences.
        preparsed_dir: Writable directory where HOMER can store preparsed genome
            files for custom FASTA-based runs.
        output_dir: HOMER output directory for this run.
        size_arg: HOMER ``-size`` argument, such as ``"given"`` or ``"200"``.
        background_bed: Optional matched background BED file. If provided, it is
            passed to HOMER with ``-bg``.
        job_name: Job name used in the sbatch header and log messages.

    Returns:
        None. The function writes the sbatch script to disk and makes it
        executable.
    """
    helpers.ensure_dir(script_path.parent)
    helpers.ensure_dir(log_path.parent)
    helpers.ensure_dir(output_dir)
    helpers.ensure_dir(preparsed_dir)

    sbatch_lines = [
        "#!/bin/bash",
        f"#SBATCH -p {partition}",
        f"#SBATCH -N {nodes}",
        f"#SBATCH --ntasks={ntasks}",
        f"#SBATCH -t {time_limit}",
        f"#SBATCH -o {log_path}",
        f"#SBATCH -J {job_name}",
    ]

    if mem:
        sbatch_lines.append(f"#SBATCH --mem={mem}")

    if allocation_id:
        sbatch_lines.append(f"#SBATCH -A {allocation_id}")

    sbatch_header = "\n".join(sbatch_lines)

    bg_arg = ""
    if background_bed is not None:
        bg_arg = f' \\\n  -bg "{background_bed}"'

    script_text = f"""{sbatch_header}

set -euo pipefail

echo "Starting HOMER job: {job_name}"
echo "Host: $(hostname)"
echo "Date: $(date)"

module load homer

FIND_MOTIFS="$(which findMotifsGenome.pl)"
if [ -z "$FIND_MOTIFS" ]; then
    echo "ERROR: findMotifsGenome.pl not found on PATH after module load homer"
    exit 1
fi

echo "FIND_MOTIFS: $FIND_MOTIFS"
echo "GENOME_FASTA: {genome_fasta}"
echo "PREPARSED_DIR: {preparsed_dir}"
echo "PATH: $PATH"
echo "ntasks: {ntasks}"

which findMotifsGenome.pl
mkdir -p "{preparsed_dir}"

"$FIND_MOTIFS" \\
  "{input_bed}" \\
  "{genome_fasta}" \\
  "{output_dir}" \\
  -size {size_arg} \\
  -p {ntasks} \\
  -preparsedDir "{preparsed_dir}"{bg_arg}

echo "Finished HOMER job: {job_name}"
echo "Date: $(date)"
"""

    script_path.write_text(script_text)
    script_path.chmod(0o755)


def run_homer(config: dict) -> Dict[str, Path]:
    """Prepare and submit or run HOMER motif enrichment jobs.

    This function validates required files and executables, builds the set of
    HOMER jobs, writes one sbatch script per job, and either submits or runs
    each job depending on the cluster configuration. If any job is submitted or
    launched, the function exits so the pipeline can be rerun after HOMER
    finishes.

    Args:
        config: Pipeline configuration dictionary loaded from YAML.

    Returns:
        A dictionary mapping each HOMER run name to its output directory.

    Raises:
        FileNotFoundError: If a required BED file, FASTA file, or executable is
            missing.
        SystemExit: If one or more HOMER jobs are submitted or launched in this
            run.
    """
    verbose = bool(config.get("project", {}).get("verbose", False))

    helpers.require_executable("findMotifsGenome.pl", "HOMER findMotifsGenome.pl")

    cluster = config.get("cluster", {})
    use_sbatch = bool(cluster.get("use_sbatch", False))
    partition = cluster.get("partition", "RM-shared")
    time_limit = cluster.get("time", "08:00:00")
    nodes = int(cluster.get("nodes", 1))
    ntasks = int(cluster.get("ntasks", 1))
    mem = str(cluster.get("mem", ""))
    allocation_id = str(cluster.get("allocation_id", ""))

    jobs = build_homer_jobs(config)

    outputs: Dict[str, Path] = {}
    submitted_any_job = False

    print("=" * 80)
    print("Preparing HOMER motif enrichment jobs")
    print("=" * 80)

    helpers.vprint(verbose, f"use_sbatch: {use_sbatch}")
    helpers.vprint(verbose, f"partition: {partition}")
    helpers.vprint(verbose, f"time: {time_limit}")
    helpers.vprint(verbose, f"nodes: {nodes}")
    helpers.vprint(verbose, f"ntasks: {ntasks}")
    helpers.vprint(verbose, f"mem: {mem}")
    helpers.vprint(verbose, f"allocation_id: {allocation_id}")

    for job in jobs:
        label = str(job["label"])
        mode = str(job["mode"])
        size_arg = str(job["size_arg"])
        use_background = bool(job["use_background"])

        input_bed = Path(job["input_bed"])
        background_bed = Path(job["background_bed"])
        genome_fasta = Path(job["genome_fasta"])
        preparsed_dir = Path(job["preparsed_dir"])
        output_dir = Path(job["output_dir"])
        script_path = Path(job["script_path"])
        log_path = Path(job["log_path"])

        helpers.require_file(input_bed, f"{label} input BED")
        helpers.require_file(genome_fasta, f"{label} genome FASTA")
        if use_background:
            helpers.require_file(background_bed, f"{label} background BED")

        run_name = f"{label}__{mode}"
        outputs[run_name] = output_dir

        helpers.vprint(verbose, f"Preparing HOMER job: {run_name}")
        helpers.vprint(verbose, f"  input: {input_bed}")
        helpers.vprint(verbose, f"  genome FASTA: {genome_fasta}")
        helpers.vprint(verbose, f"  preparsed dir: {preparsed_dir}")
        helpers.vprint(
            verbose,
            f"  background: {background_bed if use_background else 'HOMER default'}",
        )
        helpers.vprint(verbose, f"  size: {size_arg}")
        helpers.vprint(verbose, f"  output: {output_dir}")

        if homer_results_exist(output_dir):
            print(f"Skipping HOMER job; output already exists: {run_name}")
            continue

        write_homer_sbatch_script(
            script_path=script_path,
            log_path=log_path,
            partition=partition,
            time_limit=time_limit,
            nodes=nodes,
            ntasks=ntasks,
            mem=mem,
            allocation_id=allocation_id,
            input_bed=input_bed,
            genome_fasta=genome_fasta,
            preparsed_dir=preparsed_dir,
            output_dir=output_dir,
            size_arg=size_arg,
            background_bed=background_bed if use_background else None,
            job_name=f"homer_{run_name}",
        )

        print(f"Submitting/running HOMER job: {run_name}")
        helpers.submit_or_run_job(use_sbatch, script_path, verbose=verbose)
        submitted_any_job = True

    if submitted_any_job:
        raise SystemExit(
            "Submitted HOMER job(s). Re-run the pipeline after HOMER finishes."
        )

    print("HOMER outputs already exist; continuing pipeline.")
    return outputs
