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
   -> compare default HOMER background vs matched background BEDs

2. Do we expect motifs to be randomly distributed in open chromatin regions?
   -> compare full-region runs (-size given) vs centered fixed-width runs (-size 200)

Expected upstream inputs:
- results/bedtools/cross_species_ep/human_open_in_mouse_enhancers.bed
- results/bedtools/cross_species_ep/human_open_in_mouse_promoters.bed
- results/bedtools/cross_species_ep/mouse_open_in_human_enhancers.bed
- results/bedtools/cross_species_ep/mouse_open_in_human_promoters.bed

Matched background sets:
- results/bedtools/promoter_enhancer/human_pancreas_enhancers.bed
- results/bedtools/promoter_enhancer/human_pancreas_promoters.bed
- results/bedtools/promoter_enhancer/mouse_pancreas_enhancers.bed
- results/bedtools/promoter_enhancer/mouse_pancreas_promoters.bed
"""


def homer_results_exist(output_dir: Path) -> bool:
    """Return True if a HOMER output directory looks complete enough to reuse."""
    if not output_dir.exists():
        return False

    known_dir = output_dir / "knownResults"
    homer_dir = output_dir / "homerResults"

    return known_dir.exists() or homer_dir.exists()


def build_homer_jobs(config: dict) -> List[Dict[str, object]]:
    """Build all HOMER jobs to run."""
    results_dir = Path(config["project"]["output_dir"])

    cross_species_dir = results_dir / "bedtools" / "cross_species_ep"
    promoter_enhancer_dir = results_dir / "bedtools" / "promoter_enhancer"
    homer_dir = results_dir / "homer"

    datasets = [
        {
            "label": "human_open_in_mouse_enhancers",
            "input_bed": cross_species_dir / "human_open_in_mouse_enhancers.bed",
            "background_bed": promoter_enhancer_dir / "human_pancreas_enhancers.bed",
            "genome": "hg38",
        },
        {
            "label": "human_open_in_mouse_promoters",
            "input_bed": cross_species_dir / "human_open_in_mouse_promoters.bed",
            "background_bed": promoter_enhancer_dir / "human_pancreas_promoters.bed",
            "genome": "hg38",
        },
        {
            "label": "mouse_open_in_human_enhancers",
            "input_bed": cross_species_dir / "mouse_open_in_human_enhancers.bed",
            "background_bed": promoter_enhancer_dir / "mouse_pancreas_enhancers.bed",
            "genome": "mm10",
        },
        {
            "label": "mouse_open_in_human_promoters",
            "input_bed": cross_species_dir / "mouse_open_in_human_promoters.bed",
            "background_bed": promoter_enhancer_dir / "mouse_pancreas_promoters.bed",
            "genome": "mm10",
        },
    ]

    analysis_modes = [
        {
            "mode": "default_background",
            "size_arg": "given",
            "use_background": False,
        },
        {
            "mode": "matched_background",
            "size_arg": "given",
            "use_background": True,
        },
        {
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
                    "genome": dataset["genome"],
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
    genome: str,
    output_dir: Path,
    size_arg: str,
    background_bed: Optional[Path],
    job_name: str,
) -> None:
    """Write an sbatch script for one HOMER job."""
    helpers.ensure_dir(script_path.parent)
    helpers.ensure_dir(log_path.parent)
    helpers.ensure_dir(output_dir)

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
        bg_arg = f' -bg "{background_bed}"'

    script_text = f"""{sbatch_header}

set -euo pipefail

echo "Starting HOMER job: {job_name}"
echo "Host: $(hostname)"
echo "Date: $(date)"

module load homer

findMotifsGenome.pl \\
  "{input_bed}" \\
  {genome} \\
  "{output_dir}" \\
  -size {size_arg}{bg_arg}

echo "Finished HOMER job: {job_name}"
echo "Date: $(date)"
"""

    script_path.write_text(script_text)
    script_path.chmod(0o755)


def run_homer(config: dict) -> Dict[str, Path]:
    """Prepare and submit/run HOMER jobs.

    If one or more jobs are submitted/launched, this function exits so the user
    can re-run the pipeline after the jobs finish.
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
        genome = str(job["genome"])
        size_arg = str(job["size_arg"])
        use_background = bool(job["use_background"])

        input_bed = Path(job["input_bed"])
        background_bed = Path(job["background_bed"])
        output_dir = Path(job["output_dir"])
        script_path = Path(job["script_path"])
        log_path = Path(job["log_path"])

        helpers.require_file(input_bed, f"{label} input BED")
        helpers.require_file(background_bed, f"{label} background BED")

        run_name = f"{label}__{mode}"
        outputs[run_name] = output_dir

        helpers.vprint(verbose, f"Preparing HOMER job: {run_name}")
        helpers.vprint(verbose, f"  input: {input_bed}")
        helpers.vprint(verbose, f"  genome: {genome}")
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
            genome=genome,
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