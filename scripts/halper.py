#!/usr/bin/env python3

from pathlib import Path
from typing import List
import gzip
import shutil
import subprocess


def require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")


def require_executable(name: str, label: str) -> None:
    if shutil.which(name) is None:
        raise FileNotFoundError(f"{label} not found on PATH: {name}")


def run_command(cmd: List[str]) -> None:
    print("Running:")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def decompress_if_needed(input_path: Path, output_path: Path) -> Path:
    if input_path.suffix != ".gz":
        return input_path

    with gzip.open(input_path, "rt") as fin, open(output_path, "w") as fout:
        for line in fin:
            fout.write(line)

    return output_path


def make_bed4_from_narrowpeak(input_path: Path, output_path: Path, prefix: str) -> None:
    seen = set()

    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for i, line in enumerate(fin, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Invalid narrowPeak line in {input_path}: {line}")

            chrom = fields[0]
            start = fields[1]
            end = fields[2]

            if len(fields) >= 4 and fields[3] not in ("", "."):
                name = fields[3]
            else:
                name = f"{prefix}_peak_{i:06d}"

            if name in seen:
                name = f"{name}_{i:06d}"
            seen.add(name)

            fout.write(f"{chrom}\t{start}\t{end}\t{name}\n")


def make_summit_bed_from_narrowpeak(input_path: Path, output_path: Path, prefix: str) -> None:
    seen = set()

    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for i, line in enumerate(fin, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise ValueError(f"Expected 10 columns in narrowPeak file: {input_path}")

            chrom = fields[0]
            start = int(fields[1])
            summit_offset = int(fields[9])

            if len(fields) >= 4 and fields[3] not in ("", "."):
                name = fields[3]
            else:
                name = f"{prefix}_peak_{i:06d}"

            if name in seen:
                name = f"{name}_{i:06d}"
            seen.add(name)

            summit = start + summit_offset
            fout.write(f"{chrom}\t{summit}\t{summit + 1}\t{name}\n")


def run_halper_one_direction(
    peak_file: Path,
    hal_file: Path,
    source_species: str,
    target_species: str,
    output_dir: Path,
    orthologfind_py: Path,
    min_len: int,
    max_len: int,
    protect_dist: int,
) -> None:
    ensure_dir(output_dir)
    tmp_dir = output_dir / "tmp"
    ensure_dir(tmp_dir)

    prefix = f"{source_species}_to_{target_species}"

    peak_plain = tmp_dir / f"{prefix}.input.narrowPeak"
    query_bed = tmp_dir / f"{prefix}.query.bed"
    summit_bed = tmp_dir / f"{prefix}.summits.bed"
    mapped_bed = tmp_dir / f"{prefix}.mapped.bed"
    mapped_summits_bed = tmp_dir / f"{prefix}.mapped.summits.bed"
    ortholog_bed = output_dir / f"{prefix}.orthologs.bed"

    require_file(peak_file, "Peak file")
    require_file(hal_file, "HAL file")
    require_file(orthologfind_py, "orthologFind.py")
    require_executable("halLiftover", "halLiftover")

    peak_plain_path = decompress_if_needed(peak_file, peak_plain)

    make_bed4_from_narrowpeak(peak_plain_path, query_bed, prefix)
    make_summit_bed_from_narrowpeak(peak_plain_path, summit_bed, prefix)

    run_command([
        "halLiftover",
        "--bedType", "4",
        str(hal_file),
        source_species,
        str(query_bed),
        target_species,
        str(mapped_bed),
    ])

    run_command([
        "halLiftover",
        str(hal_file),
        source_species,
        str(summit_bed),
        target_species,
        str(mapped_summits_bed),
    ])

    run_command([
        "python",
        str(orthologfind_py),
        "-max_len", str(max_len),
        "-min_len", str(min_len),
        "-protect_dist", str(protect_dist),
        "-qFile", str(query_bed),
        "-tFile", str(mapped_bed),
        "-sFile", str(mapped_summits_bed),
        "-oFile", str(ortholog_bed),
        "-mult_keepone",
    ])

    if not ortholog_bed.exists():
        raise RuntimeError(f"HALPER finished but output file was not created: {ortholog_bed}")


def run_halper(config: dict) -> None:
    peak_1 = Path(config["peaks"]["species_1_peak_file"])
    peak_2 = Path(config["peaks"]["species_2_peak_file"])
    hal_file = Path(config["alignments"]["hal_file"])

    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]
    species_1_hal = config["comparison"]["species_1_hal_name"]
    species_2_hal = config["comparison"]["species_2_hal_name"]

    orthologfind_py = Path(config["tools"]["orthologfind_py"])
    output_root = Path(config["project"]["output_dir"]) / "halper"

    min_len = int(config["parameters"].get("min_peak_length", 50))
    max_len = int(config["parameters"].get("max_peak_length", 1000))
    protect_dist = int(config["parameters"].get("protect_dist", 5))
    bidirectional = bool(config["parameters"].get("run_bidirectional_liftover", True))

    print("Running HALPER...")
    print(f"  species 1 name: {species_1_name}")
    print(f"  species 2 name: {species_2_name}")
    print(f"  species 1 HAL name: {species_1_hal}")
    print(f"  species 2 HAL name: {species_2_hal}")
    print(f"  species 1 peak file: {peak_1}")
    print(f"  species 2 peak file: {peak_2}")
    print(f"  HAL file: {hal_file}")
    print(f"  orthologFind.py: {orthologfind_py}")
    print(f"  output root: {output_root}")

    run_halper_one_direction(
        peak_file=peak_1,
        hal_file=hal_file,
        source_species=species_1_hal,
        target_species=species_2_hal,
        output_dir=output_root / f"{species_1_name}_to_{species_2_name}",
        orthologfind_py=orthologfind_py,
        min_len=min_len,
        max_len=max_len,
        protect_dist=protect_dist,
    )

    if bidirectional:
        run_halper_one_direction(
            peak_file=peak_2,
            hal_file=hal_file,
            source_species=species_2_hal,
            target_species=species_1_hal,
            output_dir=output_root / f"{species_2_name}_to_{species_1_name}",
            orthologfind_py=orthologfind_py,
            min_len=min_len,
            max_len=max_len,
            protect_dist=protect_dist,
        )