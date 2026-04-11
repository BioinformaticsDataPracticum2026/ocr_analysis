#!/usr/bin/env python3

from typing import Dict, List
from pathlib import Path
import shutil
import subprocess


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


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def require_executable(name: str, label: str) -> None:
    if shutil.which(name) is None:
        raise FileNotFoundError(f"{label} not found on PATH: {name}")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def run_bedtools_to_file(cmd: List[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    try:
        with open(tmp_path, "w") as fout:
            print("Running:")
            print("  " + " ".join(cmd))
            subprocess.run(cmd, stdout=fout, check=True)
        tmp_path.replace(output_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def split_by_tss_distance(
    annotated_bed: Path,
    promoter_output: Path,
    enhancer_output: Path,
    promoter_max_distance: int,
) -> None:
    require_file(annotated_bed, "TSS-annotated BED file")

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
    require_file(input_bed, "BED file to sort")
    run_bedtools_to_file(
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
    require_file(peak_bed, "Peak BED file")
    require_file(tss_bed, "TSS BED file")

    sorted_peak_bed = annotated_output.parent / f"{peak_bed.stem}.sorted.bed"
    sorted_tss_bed = annotated_output.parent / f"{tss_bed.stem}.sorted.bed"

    print(f"Sorting peak BED: {peak_bed}")
    sort_bed(peak_bed, sorted_peak_bed)

    print(f"Sorting TSS BED: {tss_bed}")
    sort_bed(tss_bed, sorted_tss_bed)

    run_bedtools_to_file(
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
    require_executable("bedtools", "BEDTools")

    results_dir = Path(config["project"]["output_dir"]) / "bedtools" / "promoter_enhancer"
    ensure_dir(results_dir)

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
    s1_annotated = results_dir / f"{species_1_name}_{organ}_tss_annotated.bed"
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
    print(f"Wrote: {s1_annotated}")
    print(f"Wrote: {s1_promoters}")
    print(f"Wrote: {s1_enhancers}")

    # species 2
    s2_annotated = results_dir / f"{species_2_name}_{organ}_tss_annotated.bed"
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
    print(f"Wrote: {s2_annotated}")
    print(f"Wrote: {s2_promoters}")
    print(f"Wrote: {s2_enhancers}")

    print("Promoter/enhancer classification complete")
    return outputs