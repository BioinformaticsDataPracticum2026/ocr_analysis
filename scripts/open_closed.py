#!/usr/bin/env python3

from typing import Dict
from pathlib import Path

import utils.helpers as helpers

"""
Here, we classify OCRs as open or closed in the target species.
The idea is:
- if a HALPER-mapped peak overlaps a real peak in the target species, it is "open".
- if it does not overlap any real peak in the target species, it is "closed".

The overlap is determined using BEDTools intersect:
- human_to_mouse_halper.bed vs mouse_pancreas_peaks.bed
- mouse_to_human_halper.bed vs human_pancreas_peaks.bed

human_peaks_open_in_mouse.bed
- human peaks whose mapped ortholog overlaps a mouse open chromatin peak
human_peaks_closed_in_mouse.bed
- human peaks whose mapped ortholog does not overlap a mouse open chromatin peak
mouse_peaks_open_in_human.bed
- mouse peaks whose mapped ortholog overlaps a human open chromatin peak
mouse_peaks_closed_in_human.bed
- mouse peaks whose mapped ortholog does not overlap a human open chromatin peak
"""


def classify_open_closed(
    mapped_bed: Path,
    target_peak_bed: Path,
    open_output: Path,
    closed_output: Path,
    verbose: bool = True,
) -> None:
    """Classify mapped peaks as open or closed in the target species.

    A mapped peak is labeled as open if it overlaps a real peak in the target
    species, and closed if it does not overlap any target peak.

    Args:
        mapped_bed: BED file of lifted-over or mapped peaks.
        target_peak_bed: BED file of real peaks in the target species.
        open_output: Output BED file for mapped peaks classified as open.
        closed_output: Output BED file for mapped peaks classified as closed.
        verbose: Whether to print detailed command information.

    Raises:
        FileNotFoundError: If either input BED file does not exist.
    """
    helpers.require_file(mapped_bed, "Mapped BED file")
    helpers.require_file(target_peak_bed, "Target peak BED file")

    helpers.run_bedtools_to_file(
        [
            "bedtools",
            "intersect",
            "-a",
            str(mapped_bed),
            "-b",
            str(target_peak_bed),
            "-u",
        ],
        open_output,
        verbose=verbose,
    )

    helpers.run_bedtools_to_file(
        [
            "bedtools",
            "intersect",
            "-a",
            str(mapped_bed),
            "-b",
            str(target_peak_bed),
            "-v",
        ],
        closed_output,
        verbose=verbose,
    )


def run_open_closed(config: dict) -> Dict[str, Path]:
    """Run open/closed peak classification for configured species comparisons.

    This function checks for required inputs, builds expected BED file paths,
    runs open/closed classification for forward liftover, and optionally runs
    the reverse direction if bidirectional liftover is enabled.

    Args:
        config: Pipeline configuration dictionary containing project,
            comparison, and parameter settings.

    Returns:
        A dictionary mapping result labels to output BED file paths.

    Raises:
        FileNotFoundError: If BEDTools is not available on PATH.
        KeyError: If required configuration fields are missing.
    """
    verbose = bool(config.get("project", {}).get("verbose", False))

    helpers.require_executable("bedtools", "BEDTools")

    results_dir = Path(config["project"]["output_dir"]) / "bedtools" / "open_closed"
    helpers.ensure_dir(results_dir)

    organ = config["comparison"]["organ"]
    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]

    bedtools_dir = Path(config["project"]["output_dir"]) / "bedtools"

    species_1_peak_bed = bedtools_dir / f"{species_1_name}_{organ}_peaks.bed"
    species_2_peak_bed = bedtools_dir / f"{species_2_name}_{organ}_peaks.bed"

    species_1_to_species_2_halper_bed = Path(
        config["halper_outputs"]["species_1_to_species_2"]
    )
    species_2_to_species_1_halper_bed = Path(
        config["halper_outputs"]["species_2_to_species_1"]
    )

    outputs: Dict[str, Path] = {}

    # species 1 -> species 2
    if species_1_to_species_2_halper_bed.exists():
        s1_open_in_s2 = results_dir / f"{species_1_name}_peaks_open_in_{species_2_name}.bed"
        s1_closed_in_s2 = results_dir / f"{species_1_name}_peaks_closed_in_{species_2_name}.bed"

        print(f"Classifying {species_1_name} peaks in {species_2_name}")
        classify_open_closed(
            mapped_bed=species_1_to_species_2_halper_bed,
            target_peak_bed=species_2_peak_bed,
            open_output=s1_open_in_s2,
            closed_output=s1_closed_in_s2,
            verbose=verbose,
        )

        outputs["species_1_open_in_species_2"] = s1_open_in_s2
        outputs["species_1_closed_in_species_2"] = s1_closed_in_s2
        print(f"Wrote: {s1_open_in_s2}")
        print(f"Wrote: {s1_closed_in_s2}")
    else:
        print(f"Skipping missing mapped BED: {species_1_to_species_2_halper_bed}")

    # species 2 -> species 1
    bidirectional = bool(config.get("parameters", {}).get("run_bidirectional_liftover", False))
    if bidirectional:
        if species_2_to_species_1_halper_bed.exists():
            s2_open_in_s1 = results_dir / f"{species_2_name}_peaks_open_in_{species_1_name}.bed"
            s2_closed_in_s1 = results_dir / f"{species_2_name}_peaks_closed_in_{species_1_name}.bed"

            print(f"Classifying {species_2_name} peaks in {species_1_name}")
            classify_open_closed(
                mapped_bed=species_2_to_species_1_halper_bed,
                target_peak_bed=species_1_peak_bed,
                open_output=s2_open_in_s1,
                closed_output=s2_closed_in_s1,
                verbose=verbose,
            )

            outputs["species_2_open_in_species_1"] = s2_open_in_s1
            outputs["species_2_closed_in_species_1"] = s2_closed_in_s1
            print(f"Wrote: {s2_open_in_s1}")
            print(f"Wrote: {s2_closed_in_s1}")
        else:
            print(f"Skipping missing mapped BED: {species_2_to_species_1_halper_bed}")
    else:
        print("Bidirectional liftover disabled; skipping reverse open/closed analysis")

    print("Open/closed classification complete")
    return outputs