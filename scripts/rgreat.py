#!/usr/bin/env python3

from pathlib import Path
from typing import Dict

import utils.helpers as helpers

"""
run rGREAT on the final cross-species promoter/enhancer BED files.

expected inputs:
- human_open_in_mouse_promoters.bed   -> genome mm10
- human_open_in_mouse_enhancers.bed   -> genome mm10
- mouse_open_in_human_promoters.bed   -> genome hg38
- mouse_open_in_human_enhancers.bed   -> genome hg38
"""


def run_rgreat(config: dict) -> Dict[str, Path]:
    """Run rGREAT for the four cross-species promoter/enhancer BED sets.

    Args:
        config: Pipeline configuration dictionary.

    Returns:
        Dictionary mapping analysis labels to their output directories.

    Raises:
        FileNotFoundError: If required BED files or the R script are missing.
    """
    verbose = bool(config.get("project", {}).get("verbose", False))

    helpers.require_executable("Rscript", "Rscript")

    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]

    output_dir = Path(config["project"]["output_dir"])
    cross_ep_dir = output_dir / "bedtools" / "cross_species_ep"
    rgreat_dir = output_dir / "rgreat"

    helpers.ensure_dir(rgreat_dir)

    r_script = Path("scripts") / "rgreat_online.R"
    helpers.require_file(r_script, "rGREAT R script")

    jobs = [
        {
            "label": f"{species_1_name}_open_in_{species_2_name}_promoters",
            "bed_file": cross_ep_dir / f"{species_1_name}_open_in_{species_2_name}_promoters.bed",
            "genome": "mm10",
        },
        {
            "label": f"{species_1_name}_open_in_{species_2_name}_enhancers",
            "bed_file": cross_ep_dir / f"{species_1_name}_open_in_{species_2_name}_enhancers.bed",
            "genome": "mm10",
        },
        {
            "label": f"{species_2_name}_open_in_{species_1_name}_promoters",
            "bed_file": cross_ep_dir / f"{species_2_name}_open_in_{species_1_name}_promoters.bed",
            "genome": "hg38",
        },
        {
            "label": f"{species_2_name}_open_in_{species_1_name}_enhancers",
            "bed_file": cross_ep_dir / f"{species_2_name}_open_in_{species_1_name}_enhancers.bed",
            "genome": "hg38",
        },
    ]

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
