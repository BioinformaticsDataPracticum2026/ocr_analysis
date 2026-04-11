#!/usr/bin/env python3

from pathlib import Path
import gzip
import shutil
import yaml

'''
This script cleans the files from HALPER into proper BED format.
It only reformats existing peak/HALPER files into simple BED3 files. 

Here's how the files are processed:
- Human_to_Mouse.orthologs.bed from HALPER > human_to_mouse_halper.bed
- Mouse_to_Human.orthologs.bed from HALPER > mouse_to_human_halper.bed
- idr.optimal_peak.narrowPeak.gz -> human/mouse_pancreas_peaks.bed
'''

def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def write_bed3(input_path: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open_maybe_gzip(input_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"File has fewer than 3 columns: {input_path}")

            fout.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\n")


def run_bedtools_preprocess(config) -> dict:
    bedtools_dir = Path(config["project"]["output_dir"]) / "bedtools"
    tmp_dir = bedtools_dir / "tmp"

    bedtools_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]
    organ = config["comparison"]["organ"]

    species_1_peak_out = bedtools_dir / f"{species_1_name}_{organ}_peaks.bed"
    species_2_peak_out = bedtools_dir / f"{species_2_name}_{organ}_peaks.bed"

    print(f"Preprocessing peak file: {config['peaks']['species_1_peak_file']}")
    write_bed3(Path(config["peaks"]["species_1_peak_file"]), species_1_peak_out)
    print(f"Wrote: {species_1_peak_out}")

    print(f"Preprocessing peak file: {config['peaks']['species_2_peak_file']}")
    write_bed3(Path(config["peaks"]["species_2_peak_file"]), species_2_peak_out)
    print(f"Wrote: {species_2_peak_out}")

    outputs = {
        "output_dir": bedtools_dir,
        "tmp_dir": tmp_dir,
        "species_1_peak_bed3": species_1_peak_out,
        "species_2_peak_bed3": species_2_peak_out,
    }

    bidirectional = bool(
        config.get("parameters", {}).get("run_bidirectional_liftover", False)
    )
    halper_outputs = config.get("halper_outputs", {})

    if halper_outputs.get("species_1_to_species_2"):
        s1_to_s2_in = Path(halper_outputs["species_1_to_species_2"])
        s1_to_s2_out = bedtools_dir / f"{species_1_name}_to_{species_2_name}_halper.bed"

        if s1_to_s2_in.exists():
            print(f"Preprocessing HALPER file: {s1_to_s2_in}")
            write_bed3(s1_to_s2_in, s1_to_s2_out)
            print(f"Wrote: {s1_to_s2_out}")
            outputs["species_1_to_species_2_halper_bed3"] = s1_to_s2_out
        else:
            print(f"Skipping missing HALPER file: {s1_to_s2_in}")

    if bidirectional and halper_outputs.get("species_2_to_species_1"):
        s2_to_s1_in = Path(halper_outputs["species_2_to_species_1"])
        s2_to_s1_out = bedtools_dir / f"{species_2_name}_to_{species_1_name}_halper.bed"

        if s2_to_s1_in.exists():
            print(f"Preprocessing HALPER file: {s2_to_s1_in}")
            write_bed3(s2_to_s1_in, s2_to_s1_out)
            print(f"Wrote: {s2_to_s1_out}")
            outputs["species_2_to_species_1_halper_bed3"] = s2_to_s1_out
        else:
            print(f"Skipping missing HALPER file: {s2_to_s1_in}")
    elif not bidirectional:
        print("Bidirectional liftover disabled; skipping species_2_to_species_1 HALPER preprocessing")

    print("BEDTools preprocessing complete")
    print(f"Output directory: {bedtools_dir}")

    return outputs