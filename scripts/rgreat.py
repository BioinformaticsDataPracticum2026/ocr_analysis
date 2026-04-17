#!/usr/bin/env python3

from pathlib import Path
from typing import Dict

import utils.helpers as helpers

"""
run rGREAT on the final open/closed peaks BED files.

expected inputs:
- human_peaks_closed_in_mouse.bed   -> genome hg38
- mouse_peaks_closed_in_mouse.bed   -> genome mm10
- human_peaks_open_in_mouse.bed   -> shared
- mouse_peaks_open_in_human.bed   -> shared
"""


def run_rgreat(config: dict) -> Dict[str, Path]:
    """Run rGREAT for the four cross-species promoter/enhancer BED sets.

    Args:
        config: Pipeline configuration dictionary.

    Returns:
        Dictionary mapping analysis labels to their output directories.

    Raises:
        FileNotFoundError: If required directories, BED files, or the R script
            are missing.
        NotADirectoryError: If an expected directory path exists but is not a
            directory.
    """
    verbose = bool(config.get("project", {}).get("verbose", False))

    helpers.require_executable("Rscript", "Rscript")

    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]
    species_1_genome = config["rGreat"]["species_1_genome"]
    species_2_genome = config["rGreat"]["species_2_genome"]

    output_dir = Path(config["project"]["output_dir"])
    bedtools_dir = output_dir / "bedtools"
    open_closed_dir = bedtools_dir / "open_closed"
    rgreat_dir = output_dir / "rgreat"

    r_script = Path("scripts") / "rgreat_online.R"
    helpers.require_file(r_script, "rGREAT R script")

    if not bedtools_dir.exists():
        raise FileNotFoundError(
            f"BEDTools output directory not found: {bedtools_dir}\n"
            "Run the BEDTools steps before running rGREAT."
        )
    if not bedtools_dir.is_dir():
        raise NotADirectoryError(f"BEDTools path is not a directory: {bedtools_dir}")

    if not open_closed_dir.exists():
        raise FileNotFoundError(
            f"Cross-species promoter/enhancer directory not found: {open_closed_dir}\n"
            "Run cross-species promoter/enhancer classification before rGREAT."
        )
    if not open_closed_dir.is_dir():
        raise NotADirectoryError(
            f"Cross-species promoter/enhancer path is not a directory: {open_closed_dir}"
        )

    bed_files = sorted(open_closed_dir.glob("*.bed"))
    if not bed_files:
        raise FileNotFoundError(
            f"No BED files found in: {open_closed_dir}\n"
            "Run cross-species promoter/enhancer classification before rGREAT."
        )

    helpers.ensure_dir(rgreat_dir)

    jobs = [
        {
            "label": f"{species_1_name}_peaks_open_in_{species_2_name}",
            "bed_file": open_closed_dir / f"{species_1_name}_peaks_open_in_{species_2_name}.bed",
            "genome": f"{species_2_genome}", 
        },
        {
            "label": f"{species_1_name}_peaks_closed_in_{species_2_name}",
            "bed_file": open_closed_dir / f"{species_1_name}_peaks_closed_in_{species_2_name}.bed",
            "genome": f"{species_2_genome}",
        },
        {
            "label": f"{species_2_name}_peaks_open_in_{species_1_name}",
            "bed_file": open_closed_dir / f"{species_2_name}_peaks_open_in_{species_1_name}.bed",
            "genome": f"{species_1_genome}",
        },
        {
            "label": f"{species_2_name}_peaks_closed_in_{species_1_name}",
            "bed_file": open_closed_dir / f"{species_2_name}_peaks_closed_in_{species_1_name}.bed",
            "genome": f"{species_1_genome}",
        },
        {
            "label": f"{species_1_name}_pancreas_peaks",
            "bed_file": bedtools_dir / f"{species_1_name}_pancreas_peaks.bed",
            "genome": f"{species_1_genome}",
        },
        {
            "label": f"{species_2_name}_pancreas_peaks",
            "bed_file": bedtools_dir / f"{species_2_name}_pancreas_peaks.bed",
            "genome": f"{species_2_genome}",
        },
    ]

    helpers.vprint(verbose, f"BEDTools directory: {bedtools_dir}")
    helpers.vprint(verbose, f"Open-Closed EP directory: {open_closed_dir}")
    helpers.vprint(verbose, f"Found {len(bed_files)} BED file(s) in {open_closed_dir}")
    helpers.vprint(verbose, f"rGREAT output directory: {rgreat_dir}")
    helpers.vprint(verbose, f"Using rGREAT script: {r_script}")
    helpers.vprint(verbose, f"Configured rGREAT jobs: {len(jobs)}")

    outputs: Dict[str, Path] = {}

    print("Running rGREAT analyses")
    for job in jobs:
        label = job["label"]
        bed_file = job["bed_file"]
        genome = job["genome"]
        run_output_dir = rgreat_dir / label

        helpers.require_file(bed_file, f"rGREAT input BED for {label}")
        helpers.ensure_dir(run_output_dir)

        metadata_file = run_output_dir / f"{label}.metadata.txt"

        helpers.vprint(verbose, f"Preparing rGREAT job: {label}")
        helpers.vprint(verbose, f"  BED file: {bed_file}")
        helpers.vprint(verbose, f"  Genome: {genome}")
        helpers.vprint(verbose, f"  Output dir: {run_output_dir}")
        helpers.vprint(verbose, f"  Metadata file: {metadata_file}")

        if metadata_file.exists():
            print(f"Skipping existing rGREAT run: {label}")
            outputs[label] = run_output_dir
            continue

        print(f"Running rGREAT for: {label}")
        helpers.run_command(
            [
                "Rscript",
                str(r_script),
                str(bed_file),
                genome,
                str(run_output_dir),
            ],
            verbose=verbose,
        )

        outputs[label] = run_output_dir
        print(f"Wrote: {run_output_dir}")

    print("rGREAT analysis complete")
    return outputs