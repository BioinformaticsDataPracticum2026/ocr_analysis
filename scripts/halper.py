#!/usr/bin/env python3

from typing import Dict
from pathlib import Path
import gzip
import sys

import utils.helpers as helpers


def decompress_if_needed(input_path: Path, output_path: Path) -> Path:
    """Decompress a gzipped file if needed.

    If the input file does not end with ``.gz``, this function returns the
    original input path unchanged. Otherwise, it writes the decompressed
    contents to ``output_path`` and returns that path.

    Args:
        input_path: Input file path, possibly gzipped.
        output_path: Output path for the decompressed file.

    Returns:
        Path to the plain-text file that should be used downstream.
    """
    if input_path.suffix != ".gz":
        return input_path

    with gzip.open(input_path, "rt") as fin, open(output_path, "w") as fout:
        for line in fin:
            fout.write(line)

    return output_path


def make_bed4_from_narrowpeak(input_path: Path, output_path: Path, prefix: str) -> None:
    """Convert a narrowPeak file to a BED4 file.

    The output contains chromosome, start, end, and a unique peak name. If the
    input already contains a usable name in column 4, that name is kept;
    otherwise a name is generated from ``prefix``.

    Args:
        input_path: Path to the input narrowPeak file.
        output_path: Path to the output BED4 file.
        prefix: Prefix used when generating peak names.

    Raises:
        ValueError: If an input line has fewer than 3 columns.
    """
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
    """Create a 1-bp summit BED file from a narrowPeak file.

    The summit position is computed as ``start + summit_offset``, where the
    summit offset is taken from column 10 of the narrowPeak file. The output
    contains chromosome, summit start, summit end, and a unique peak name.

    Args:
        input_path: Path to the input narrowPeak file.
        output_path: Path to the output summit BED file.
        prefix: Prefix used when generating peak names.

    Raises:
        ValueError: If an input line has fewer than 10 columns.
    """
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


def write_sbatch_script(
    script_path: Path,
    log_path: Path,
    partition: str,
    time_limit: str,
    nodes: int,
    ntasks: int,
    mem: str,
    allocation_id: str,
    hal_file: Path,
    source_species: str,
    target_species: str,
    query_bed: Path,
    summit_bed: Path,
    mapped_bed: Path,
    mapped_summits_bed: Path,
    orthologfind_py: Path,
    ortholog_bed: Path,
    min_len: int,
    max_len: int,
    protect_dist: int,
) -> None:
    """Write an sbatch script for one HALPER liftover job.

    The generated script runs ``halLiftover`` on both peak intervals and summit
    intervals, then runs ``orthologFind.py`` to post-process the mapped peaks.

    Args:
        script_path: Path where the sbatch script will be written.
        log_path: Path for the SLURM output log.
        partition: SLURM partition name.
        time_limit: SLURM walltime limit.
        nodes: Number of SLURM nodes to request.
        ntasks: Number of SLURM tasks to request.
        mem: Memory request string, such as ``"16G"``. If empty, no ``--mem``
            line is added.
        allocation_id: SLURM allocation/account ID. If empty, no ``-A`` line is
            added.
        hal_file: Path to the HAL alignment file.
        source_species: Source species name used in the HAL file.
        target_species: Target species name used in the HAL file.
        query_bed: BED4 file of query peaks.
        summit_bed: BED file of 1-bp query summits.
        mapped_bed: Output path for lifted query peaks.
        mapped_summits_bed: Output path for lifted query summits.
        orthologfind_py: Path to the ``orthologFind.py`` script.
        ortholog_bed: Output path for HALPER ortholog peaks.
        min_len: Minimum allowed peak length for ortholog filtering.
        max_len: Maximum allowed peak length for ortholog filtering.
        protect_dist: Protection distance parameter passed to
            ``orthologFind.py``.
    """
    helpers.ensure_dir(script_path.parent)
    helpers.ensure_dir(log_path.parent)

    sbatch_lines = [
        "#!/bin/bash",
        f"#SBATCH -p {partition}",
        f"#SBATCH -N {nodes}",
        f"#SBATCH --ntasks={ntasks}",
        f"#SBATCH -t {time_limit}",
        f"#SBATCH -o {log_path}",
        f"#SBATCH -J halper_{source_species}_to_{target_species}",
    ]

    if mem:
        sbatch_lines.append(f"#SBATCH --mem={mem}")

    if allocation_id:
        sbatch_lines.append(f"#SBATCH -A {allocation_id}")

    sbatch_header = "\n".join(sbatch_lines)

    script_text = f"""{sbatch_header}

set -euo pipefail

echo "Starting HALPER job: {source_species} -> {target_species}"
echo "Host: $(hostname)"
echo "Date: $(date)"

halLiftover --bedType 4 "{hal_file}" "{source_species}" "{query_bed}" "{target_species}" "{mapped_bed}"
halLiftover "{hal_file}" "{source_species}" "{summit_bed}" "{target_species}" "{mapped_summits_bed}"

"{sys.executable}" "{orthologfind_py}" \\
  -max_len "{max_len}" \\
  -min_len "{min_len}" \\
  -protect_dist "{protect_dist}" \\
  -qFile "{query_bed}" \\
  -tFile "{mapped_bed}" \\
  -sFile "{mapped_summits_bed}" \\
  -oFile "{ortholog_bed}" \\
  -mult_keepone

echo "Finished HALPER job: {source_species} -> {target_species}"
echo "Date: $(date)"
"""
    script_path.write_text(script_text)
    script_path.chmod(0o755)


def prepare_halper_one_direction(
    peak_file: Path,
    hal_file: Path,
    source_species: str,
    target_species: str,
    output_dir: Path,
) -> Dict[str, Path]:
    """Prepare all intermediate files for one HALPER mapping direction.

    This function creates working directories, validates required inputs,
    decompresses the peak file if needed, and generates both the BED4 peak file
    and the summit BED file used by the HALPER workflow.

    Args:
        peak_file: Input narrowPeak file for the source species.
        hal_file: HAL alignment file.
        source_species: Source species name used in the HAL file.
        target_species: Target species name used in the HAL file.
        output_dir: Output directory for this mapping direction.

    Returns:
        Dictionary containing paths to generated intermediate files, final
        ortholog output, batch script, and log file.

    Raises:
        FileNotFoundError: If the peak file or HAL file does not exist.
        ValueError: If the input narrowPeak file is malformed.
    """
    helpers.ensure_dir(output_dir)
    tmp_dir = output_dir / "tmp"
    logs_dir = output_dir / "logs"
    helpers.ensure_dir(tmp_dir)
    helpers.ensure_dir(logs_dir)

    prefix = f"{source_species}_to_{target_species}"

    peak_plain = tmp_dir / f"{prefix}.input.narrowPeak"
    query_bed = tmp_dir / f"{prefix}.query.bed"
    summit_bed = tmp_dir / f"{prefix}.summits.bed"
    mapped_bed = tmp_dir / f"{prefix}.mapped.bed"
    mapped_summits_bed = tmp_dir / f"{prefix}.mapped.summits.bed"
    ortholog_bed = output_dir / f"{prefix}.orthologs.bed"
    batch_script = tmp_dir / f"{prefix}.sbatch.sh"
    batch_log = logs_dir / f"{prefix}.log"

    helpers.require_file(peak_file, "Peak file")
    helpers.require_file(hal_file, "HAL file")

    peak_plain_path = decompress_if_needed(peak_file, peak_plain)
    make_bed4_from_narrowpeak(peak_plain_path, query_bed, prefix)
    make_summit_bed_from_narrowpeak(peak_plain_path, summit_bed, prefix)

    return {
        "query_bed": query_bed,
        "summit_bed": summit_bed,
        "mapped_bed": mapped_bed,
        "mapped_summits_bed": mapped_summits_bed,
        "ortholog_bed": ortholog_bed,
        "batch_script": batch_script,
        "batch_log": batch_log,
    }


def run_halper(config: dict) -> None:
    """Run the HALPER workflow for one or two liftover directions.

    This function reads paths and parameters from the configuration, validates
    required tools, prepares intermediate BED files, writes batch scripts, and
    either submits or runs the jobs for forward and optional reverse liftover.

    If one or more HALPER jobs are submitted or launched in this run, the
    program exits afterward so downstream steps do not continue before HALPER
    finishes.

    Args:
        config: Configuration dictionary containing peak files, alignment paths,
            tool paths, comparison metadata, workflow parameters, and optional
            cluster settings.

    Raises:
        FileNotFoundError: If a required file or executable is missing.
        KeyError: If required configuration keys are missing.
        ValueError: If a narrowPeak file is malformed.
        SystemExit: If one or more HALPER jobs are submitted or launched.
    """
    verbose = bool(config.get("project", {}).get("verbose", False))

    peak_1 = Path(config["peaks"]["species_1_peak_file"])
    peak_2 = Path(config["peaks"]["species_2_peak_file"])
    hal_file = Path(config["alignments"]["hal_file"])
    orthologfind_py = Path(config["tools"]["orthologfind_py"])

    helpers.require_file(orthologfind_py, "orthologFind.py")
    helpers.require_executable("halLiftover", "halLiftover")

    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]
    species_1_hal = config["comparison"]["species_1_hal_name"]
    species_2_hal = config["comparison"]["species_2_hal_name"]

    output_root = Path(config["project"]["output_dir"]) / "halper"
    min_len = int(config["parameters"].get("min_peak_length", 50))
    max_len = int(config["parameters"].get("max_peak_length", 1000))
    protect_dist = int(config["parameters"].get("protect_dist", 5))
    bidirectional = bool(config["parameters"].get("run_bidirectional_liftover", True))

    cluster = config.get("cluster", {})
    use_sbatch = bool(cluster.get("use_sbatch", False))
    partition = cluster.get("partition", "RM-shared")
    time_limit = cluster.get("time", "08:00:00")
    nodes = int(cluster.get("nodes", 1))
    ntasks = int(cluster.get("ntasks", 1))
    mem = cluster.get("mem")
    allocation_id = cluster.get("allocation_id")

    helpers.vprint(verbose, "Running HALPER...")
    helpers.vprint(verbose, f"  species 1 name: {species_1_name}")
    helpers.vprint(verbose, f"  species 2 name: {species_2_name}")
    helpers.vprint(verbose, f"  species 1 HAL name: {species_1_hal}")
    helpers.vprint(verbose, f"  species 2 HAL name: {species_2_hal}")
    helpers.vprint(verbose, f"  species 1 peak file: {peak_1}")
    helpers.vprint(verbose, f"  species 2 peak file: {peak_2}")
    helpers.vprint(verbose, f"  HAL file: {hal_file}")
    helpers.vprint(verbose, f"  orthologFind.py: {orthologfind_py}")
    helpers.vprint(verbose, f"  output root: {output_root}")
    helpers.vprint(verbose, f"  use sbatch: {use_sbatch}")
    helpers.vprint(verbose, f"  partition: {partition}")
    helpers.vprint(verbose, f"  time: {time_limit}")
    helpers.vprint(verbose, f"  nodes: {nodes}")
    helpers.vprint(verbose, f"  ntasks: {ntasks}")
    helpers.vprint(verbose, f"  mem: {mem}")
    helpers.vprint(verbose, f"  allocation_id: {allocation_id}")

    submitted_any_job = False

    prepared_1 = prepare_halper_one_direction(
        peak_file=peak_1,
        hal_file=hal_file,
        source_species=species_1_hal,
        target_species=species_2_hal,
        output_dir=output_root / f"{species_1_name}_to_{species_2_name}",
    )

    if helpers.should_use_existing_halper_output(prepared_1["ortholog_bed"], verbose=verbose):
        print(f"Skipping HALPER job for {species_1_name} -> {species_2_name}")
    else:
        write_sbatch_script(
            script_path=prepared_1["batch_script"],
            log_path=prepared_1["batch_log"],
            partition=partition,
            time_limit=time_limit,
            nodes=nodes,
            ntasks=ntasks,
            mem=mem,
            allocation_id=allocation_id,
            hal_file=hal_file,
            source_species=species_1_hal,
            target_species=species_2_hal,
            query_bed=prepared_1["query_bed"],
            summit_bed=prepared_1["summit_bed"],
            mapped_bed=prepared_1["mapped_bed"],
            mapped_summits_bed=prepared_1["mapped_summits_bed"],
            orthologfind_py=orthologfind_py,
            ortholog_bed=prepared_1["ortholog_bed"],
            min_len=min_len,
            max_len=max_len,
            protect_dist=protect_dist,
        )
        helpers.submit_or_run_job(use_sbatch, prepared_1["batch_script"], verbose=verbose)
        submitted_any_job = True

    if bidirectional:
        prepared_2 = prepare_halper_one_direction(
            peak_file=peak_2,
            hal_file=hal_file,
            source_species=species_2_hal,
            target_species=species_1_hal,
            output_dir=output_root / f"{species_2_name}_to_{species_1_name}",
        )

        if helpers.should_use_existing_halper_output(prepared_2["ortholog_bed"], verbose=verbose):
            print(f"Skipping HALPER job for {species_2_name} -> {species_1_name}")
        else:
            write_sbatch_script(
                script_path=prepared_2["batch_script"],
                log_path=prepared_2["batch_log"],
                partition=partition,
                time_limit=time_limit,
                nodes=nodes,
                ntasks=ntasks,
                mem=mem,
                allocation_id=allocation_id,
                hal_file=hal_file,
                source_species=species_2_hal,
                target_species=species_1_hal,
                query_bed=prepared_2["query_bed"],
                summit_bed=prepared_2["summit_bed"],
                mapped_bed=prepared_2["mapped_bed"],
                mapped_summits_bed=prepared_2["mapped_summits_bed"],
                orthologfind_py=orthologfind_py,
                ortholog_bed=prepared_2["ortholog_bed"],
                min_len=min_len,
                max_len=max_len,
                protect_dist=protect_dist,
            )
            helpers.submit_or_run_job(use_sbatch, prepared_2["batch_script"], verbose=verbose)
            submitted_any_job = True

    if submitted_any_job:
        raise SystemExit(
            "Submitted HALPER job(s). Re-run the pipeline after HALPER finishes."
        )

    print("HALPER outputs already exist; continuing pipeline.")