#!/usr/bin/env python3

from typing import List
from pathlib import Path
import shutil
import subprocess


"""
Here, we classify OCRs as open or closed in the target species.
The idea is: 
- if a HALPER-mapped peak overlaps a real peak in the target species, it is "open".
- If it does not overlap any real peak in the target species, it is "closed".

The overlap is determined using BEDTools intersect:
- human_to_mouse_halper.bed vs mouse_pancreas_peaks.bed
- mouse_to_human_halper.bed vs human_pancreas_peaks.bed
"""

def require_executable(name: str, label: str) -> None:
    if shutil.which(name) is None:
        raise FileNotFoundError(f"{label} not found on PATH: {name}")

def run_command(cmd: List[str]) -> None:
    print("Running:")
    print("  " + " ".join(str(x) for x in cmd))
    subprocess.run(cmd, check=True)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def classify_open_closed(
    mapped_bed: Path,
    target_peak_bed: Path,
    open_output: Path,
    closed_output: Path,
) -> None:
    ensure_dir(open_output.parent)
    ensure_dir(closed_output.parent)

    if not mapped_bed.exists():
        raise FileNotFoundError(f"Mapped BED file not found: {mapped_bed}")
    if not target_peak_bed.exists():
        raise FileNotFoundError(f"Target peak BED file not found: {target_peak_bed}")

    # Open in target species: mapped peak overlaps a real peak in target
    with open(open_output, "w") as fout:
        subprocess.run(
            [
                "bedtools",
                "intersect",
                "-a",
                str(mapped_bed),
                "-b",
                str(target_peak_bed),
                "-u",
            ],
            stdout=fout,
            check=True,
        )

    # Closed in target species: mapped peak does not overlap any real peak in target
    with open(closed_output, "w") as fout:
        subprocess.run(
            [
                "bedtools",
                "intersect",
                "-a",
                str(mapped_bed),
                "-b",
                str(target_peak_bed),
                "-v",
            ],
            stdout=fout,
            check=True,
        )


def run_open_closed(config: dict) -> dict:
    require_executable("bedtools", "BEDTools")
    
    results_dir = Path(config["project"]["output_dir"]) / "bedtools" / "open_closed"
    ensure_dir(results_dir)

    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]

    bedtools_dir = Path(config["project"]["output_dir"]) / "bedtools"

    species_1_peak_bed = bedtools_dir / f"{species_1_name}_{config['comparison']['organ']}_peaks.bed"
    species_2_peak_bed = bedtools_dir / f"{species_2_name}_{config['comparison']['organ']}_peaks.bed"

    species_1_to_species_2_halper_bed = bedtools_dir / f"{species_1_name}_to_{species_2_name}_halper.bed"
    species_2_to_species_1_halper_bed = bedtools_dir / f"{species_2_name}_to_{species_1_name}_halper.bed"

    outputs = {}

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
