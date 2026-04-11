#!/usr/bin/env python3

from typing import Dict, List
from pathlib import Path
import gzip
import sys
import subprocess

import utils.helpers as helpers

"""
Here, we classify OCRs as promoter-like or enhancer-like within each species.
The idea is:
- if a peak is close to a TSS, it is "promoter-like"
- if a peak is farther from a TSS, it is "enhancer-like"

The distance is determined using BEDTools closest:
- human_pancreas_peaks.bed vs human TSS BED
- mouse_pancreas_peaks.bed vs mouse TSS BED

By default:
- distance <= promoter_max_distance -> promoter-like
- distance > promoter_max_distance  -> enhancer-like
"""


def split_by_tss_distance(
    annotated_bed: Path,
    promoter_output: Path,
    enhancer_output: Path,
    promoter_max_distance: int,
) -> None:
    """Split a TSS-annotated BED file into promoter-like and enhancer-like peaks.

    This function expects the last column of each input line to contain the
    distance to the closest TSS. Peaks with distance less than or equal to
    ``promoter_max_distance`` are written to the promoter output; all others
    are written to the enhancer output.

    Args:
        annotated_bed: BED file annotated with TSS distances in the last column.
        promoter_output: Output BED file for promoter-like peaks.
        enhancer_output: Output BED file for enhancer-like peaks.
        promoter_max_distance: Maximum distance from a TSS for a peak to be
            classified as promoter-like.

    Raises:
        FileNotFoundError: If the annotated BED file does not exist.
        ValueError: If the last column of a non-empty line is not an integer.
    """
    helpers.require_file(annotated_bed, "TSS-annotated BED file")

    promoter_output.parent.mkdir(parents=True, exist_ok=True)
    enhancer_output.parent.mkdir(parents=True, exist_ok=True)

    promoter_tmp = promoter_output.with_suffix(promoter_output.suffix + ".tmp")
    enhancer_tmp = enhancer_output.with_suffix(enhancer_output.suffix + ".tmp")

    try:
        with open(annotated_bed, "r") as fin, \
             open(promoter_tmp, "w") as promoter_fout, \
             open(enhancer_tmp, "w") as enhancer_fout:

            for line in fin:
                if not line.strip():
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 1:
                    continue

                try:
                    distance = int(fields[-1])
                except ValueError:
                    raise ValueError(
                        f"Last column is not an integer distance in file: {annotated_bed}\n"
                        f"Problem line: {line}"
                    )

                if distance <= promoter_max_distance:
                    promoter_fout.write(line)
                else:
                    enhancer_fout.write(line)

        promoter_tmp.replace(promoter_output)
        enhancer_tmp.replace(enhancer_output)

    finally:
        if promoter_tmp.exists():
            promoter_tmp.unlink()
        if enhancer_tmp.exists():
            enhancer_tmp.unlink()


def sort_bed(input_bed: Path, output_bed: Path) -> None:
    """Sort a BED file using BEDTools.

    Args:
        input_bed: Input BED file to sort.
        output_bed: Output path for the sorted BED file.

    Raises:
        FileNotFoundError: If the input BED file does not exist.
        subprocess.CalledProcessError: If the BEDTools sort command fails.
    """
    helpers.require_file(input_bed, "BED file to sort")
    helpers.run_bedtools_to_file(
        [
            "bedtools",
            "sort",
            "-i",
            str(input_bed),
        ],
        output_bed,
    )


def classify_promoter_enhancer(
    peak_bed: Path,
    tss_bed: Path,
    annotated_output: Path,
    promoter_output: Path,
    enhancer_output: Path,
    promoter_max_distance: int,
) -> None:
    """Classify peaks as promoter-like or enhancer-like using closest TSS distance.

    This function sorts the peak BED and TSS BED files, annotates each peak with
    the distance to the nearest TSS using ``bedtools closest``, and then splits
    the annotated peaks into promoter-like and enhancer-like sets.

    Args:
        peak_bed: BED file of peaks to classify.
        tss_bed: BED file of transcription start sites.
        annotated_output: Output path for the TSS-annotated BED file.
        promoter_output: Output BED file for promoter-like peaks.
        enhancer_output: Output BED file for enhancer-like peaks.
        promoter_max_distance: Maximum distance from a TSS for a peak to be
            classified as promoter-like.

    Raises:
        FileNotFoundError: If the peak BED or TSS BED file does not exist.
        subprocess.CalledProcessError: If a BEDTools command fails.
        ValueError: If the annotated output contains a malformed distance column.
    """
    helpers.require_file(peak_bed, "Peak BED file")
    helpers.require_file(tss_bed, "TSS BED file")

    sorted_peak_bed = annotated_output.parent / f"{peak_bed.stem}.sorted.bed"
    sorted_tss_bed = annotated_output.parent / f"{tss_bed.stem}.sorted.bed"

    print(f"Sorting peak BED: {peak_bed}")
    sort_bed(peak_bed, sorted_peak_bed)

    print(f"Sorting TSS BED: {tss_bed}")
    sort_bed(tss_bed, sorted_tss_bed)

    helpers.run_bedtools_to_file(
        [
            "bedtools",
            "closest",
            "-a",
            str(sorted_peak_bed),
            "-b",
            str(sorted_tss_bed),
            "-d",
            "-t",
            "first",
        ],
        annotated_output,
    )

    split_by_tss_distance(
        annotated_bed=annotated_output,
        promoter_output=promoter_output,
        enhancer_output=enhancer_output,
        promoter_max_distance=promoter_max_distance,
    )


def run_promoter_enhancer(config: dict) -> Dict[str, Path]:
    """Run promoter/enhancer classification for both species in the config.

    This function reads peak and TSS file paths from the configuration, creates
    output directories, runs promoter/enhancer classification for each species,
    and returns the generated file paths.

    Args:
        config: Configuration dictionary containing project, comparison, and
            annotation settings.

    Returns:
        Dictionary mapping output labels to generated BED file paths.

    Raises:
        FileNotFoundError: If BEDTools or required input files are missing.
        KeyError: If required configuration keys are missing.
        subprocess.CalledProcessError: If a BEDTools command fails.
        ValueError: If an annotated BED file contains an invalid distance value.
    """
    helpers.require_executable("bedtools", "BEDTools")

    results_dir = Path(config["project"]["output_dir"]) / "bedtools" / "promoter_enhancer"
    tmp_dir = results_dir / "tmp"
    helpers.ensure_dir(results_dir)
    helpers.ensure_dir(tmp_dir)

    organ = config["comparison"]["organ"]
    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]

    bedtools_dir = Path(config["project"]["output_dir"]) / "bedtools"

    species_1_peak_bed = bedtools_dir / f"{species_1_name}_{organ}_peaks.bed"
    species_2_peak_bed = bedtools_dir / f"{species_2_name}_{organ}_peaks.bed"

    species_1_tss_bed = Path(config["annotations"]["species_1_tss_file"])
    species_2_tss_bed = Path(config["annotations"]["species_2_tss_file"])

    promoter_max_distance = int(
        config.get("annotations", {}).get("promoter_max_distance", 5000)
    )

    outputs: Dict[str, Path] = {}

    # species 1
    s1_annotated = tmp_dir / f"{species_1_name}_{organ}_tss_annotated.bed"
    s1_promoters = results_dir / f"{species_1_name}_{organ}_promoters.bed"
    s1_enhancers = results_dir / f"{species_1_name}_{organ}_enhancers.bed"

    print(f"Classifying {species_1_name} {organ} peaks as promoter-like or enhancer-like")
    classify_promoter_enhancer(
        peak_bed=species_1_peak_bed,
        tss_bed=species_1_tss_bed,
        annotated_output=s1_annotated,
        promoter_output=s1_promoters,
        enhancer_output=s1_enhancers,
        promoter_max_distance=promoter_max_distance,
    )

    outputs["species_1_tss_annotated"] = s1_annotated
    outputs["species_1_promoters"] = s1_promoters
    outputs["species_1_enhancers"] = s1_enhancers
    print(f"Wrote: {s1_promoters}")
    print(f"Wrote: {s1_enhancers}")

    # species 2
    s2_annotated = tmp_dir / f"{species_2_name}_{organ}_tss_annotated.bed"
    s2_promoters = results_dir / f"{species_2_name}_{organ}_promoters.bed"
    s2_enhancers = results_dir / f"{species_2_name}_{organ}_enhancers.bed"

    print(f"Classifying {species_2_name} {organ} peaks as promoter-like or enhancer-like")
    classify_promoter_enhancer(
        peak_bed=species_2_peak_bed,
        tss_bed=species_2_tss_bed,
        annotated_output=s2_annotated,
        promoter_output=s2_promoters,
        enhancer_output=s2_enhancers,
        promoter_max_distance=promoter_max_distance,
    )

    outputs["species_2_tss_annotated"] = s2_annotated
    outputs["species_2_promoters"] = s2_promoters
    outputs["species_2_enhancers"] = s2_enhancers
    print(f"Wrote: {s2_promoters}")
    print(f"Wrote: {s2_enhancers}")

    print("Promoter/enhancer classification complete")
    return outputs
