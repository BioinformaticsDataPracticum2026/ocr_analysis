#!/usr/bin/env python3

from pathlib import Path
from typing import Dict, List, Optional

import utils.helpers as helpers

"""
Get a small summary table for the bedtools outputs.
"""


def count_bed_lines(path: Path) -> int:
    """Count non-empty lines in a BED file.

    Args:
        path: BED file path.

    Returns:
        Number of non-empty lines in the file.

    Raises:
        FileNotFoundError: If the BED file does not exist.
    """
    helpers.require_file(path, "BED file")

    count = 0
    with open(path, "r") as fin:
        for line in fin:
            if line.strip():
                count += 1
    return count


def get_row(
    rows: List[dict],
    category: str,
    source_species: str,
    target_species: str,
) -> Optional[dict]:
    """Get the first summary row matching a category and species pair.

    Args:
        rows: Summary row dictionaries.
        category: Summary category to match.
        source_species: Source species name.
        target_species: Target species name.

    Returns:
        Matching row dictionary if found, otherwise None.
    """
    for row in rows:
        if (
            row["category"] == category
            and row["source_species"] == source_species
            and row["target_species"] == target_species
        ):
            return row
    return None


def build_summary_rows(config: dict) -> List[dict]:
    """Build summary rows for key BED outputs in the pipeline.

    Args:
        config: Pipeline configuration dictionary.

    Returns:
        List of row dictionaries for the summary table.
    """
    verbose = bool(config.get("project", {}).get("verbose", False))

    output_dir = Path(config["project"]["output_dir"])
    organ = config["comparison"]["organ"]
    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]

    bedtools_dir = output_dir / "bedtools"
    pe_dir = bedtools_dir / "promoter_enhancer"
    oc_dir = bedtools_dir / "open_closed"
    cross_ep_dir = bedtools_dir / "cross_species_ep"

    rows: List[dict] = []

    def add_row(
        category: str,
        path: Path,
        source_species: str,
        target_species: str,
        coordinate_species: str,
        parent_category: Optional[str] = None,
        parent_count: Optional[int] = None,
    ) -> Optional[int]:
        """Add one summary row if the file exists.

        Args:
            category: Summary category label.
            path: BED file path to summarize.
            source_species: Species the peaks originally came from.
            target_species: Species being compared against.
            coordinate_species: Species whose genome coordinates the BED uses.
            parent_category: Parent category label for percent calculations.
            parent_count: Parent count for percent calculations.

        Returns:
            Count of rows in the BED file if it exists, otherwise None.
        """
        if not path.exists():
            print(f"Skipping missing summary input: {path}")
            return None

        count = count_bed_lines(path)
        percent_of_parent = ""

        if parent_count is not None:
            percent_of_parent = f"{(count / parent_count * 100):.2f}" if parent_count > 0 else "0.00"

        row = {
            "category": category,
            "source_species": source_species,
            "target_species": target_species,
            "coordinate_species": coordinate_species,
            "count": count,
            "parent_category": parent_category or "",
            "percent_of_parent": percent_of_parent,
            "file_path": str(path),
        }
        rows.append(row)
        helpers.vprint(verbose, f"Added summary row: {category} ({count})")
        return count

    # Native peak BEDs from preprocessing
    species_1_peak_bed = bedtools_dir / f"{species_1_name}_{organ}_peaks.bed"
    species_2_peak_bed = bedtools_dir / f"{species_2_name}_{organ}_peaks.bed"

    s1_native_total = add_row(
        category="native_peaks",
        path=species_1_peak_bed,
        source_species=species_1_name,
        target_species=species_1_name,
        coordinate_species=species_1_name,
    )
    s2_native_total = add_row(
        category="native_peaks",
        path=species_2_peak_bed,
        source_species=species_2_name,
        target_species=species_2_name,
        coordinate_species=species_2_name,
    )

    # Native promoter/enhancer
    s1_promoters = pe_dir / f"{species_1_name}_{organ}_promoters.bed"
    s1_enhancers = pe_dir / f"{species_1_name}_{organ}_enhancers.bed"
    s2_promoters = pe_dir / f"{species_2_name}_{organ}_promoters.bed"
    s2_enhancers = pe_dir / f"{species_2_name}_{organ}_enhancers.bed"

    add_row(
        category="native_promoters",
        path=s1_promoters,
        source_species=species_1_name,
        target_species=species_1_name,
        coordinate_species=species_1_name,
        parent_category="native_peaks",
        parent_count=s1_native_total,
    )
    add_row(
        category="native_enhancers",
        path=s1_enhancers,
        source_species=species_1_name,
        target_species=species_1_name,
        coordinate_species=species_1_name,
        parent_category="native_peaks",
        parent_count=s1_native_total,
    )
    add_row(
        category="native_promoters",
        path=s2_promoters,
        source_species=species_2_name,
        target_species=species_2_name,
        coordinate_species=species_2_name,
        parent_category="native_peaks",
        parent_count=s2_native_total,
    )
    add_row(
        category="native_enhancers",
        path=s2_enhancers,
        source_species=species_2_name,
        target_species=species_2_name,
        coordinate_species=species_2_name,
        parent_category="native_peaks",
        parent_count=s2_native_total,
    )

    # Open/closed
    s1_open_in_s2 = oc_dir / f"{species_1_name}_peaks_open_in_{species_2_name}.bed"
    s1_closed_in_s2 = oc_dir / f"{species_1_name}_peaks_closed_in_{species_2_name}.bed"
    s2_open_in_s1 = oc_dir / f"{species_2_name}_peaks_open_in_{species_1_name}.bed"
    s2_closed_in_s1 = oc_dir / f"{species_2_name}_peaks_closed_in_{species_1_name}.bed"

    s1_open_total = add_row(
        category="open_in_target",
        path=s1_open_in_s2,
        source_species=species_1_name,
        target_species=species_2_name,
        coordinate_species=species_2_name,
    )
    add_row(
        category="closed_in_target",
        path=s1_closed_in_s2,
        source_species=species_1_name,
        target_species=species_2_name,
        coordinate_species=species_2_name,
    )

    s2_open_total = add_row(
        category="open_in_target",
        path=s2_open_in_s1,
        source_species=species_2_name,
        target_species=species_1_name,
        coordinate_species=species_1_name,
    )
    add_row(
        category="closed_in_target",
        path=s2_closed_in_s1,
        source_species=species_2_name,
        target_species=species_1_name,
        coordinate_species=species_1_name,
    )

    # Cross-species promoter/enhancer on open peaks
    s1_open_promoters = cross_ep_dir / f"{species_1_name}_open_in_{species_2_name}_promoters.bed"
    s1_open_enhancers = cross_ep_dir / f"{species_1_name}_open_in_{species_2_name}_enhancers.bed"
    s2_open_promoters = cross_ep_dir / f"{species_2_name}_open_in_{species_1_name}_promoters.bed"
    s2_open_enhancers = cross_ep_dir / f"{species_2_name}_open_in_{species_1_name}_enhancers.bed"

    add_row(
        category="open_in_target_promoters",
        path=s1_open_promoters,
        source_species=species_1_name,
        target_species=species_2_name,
        coordinate_species=species_2_name,
        parent_category="open_in_target",
        parent_count=s1_open_total,
    )
    add_row(
        category="open_in_target_enhancers",
        path=s1_open_enhancers,
        source_species=species_1_name,
        target_species=species_2_name,
        coordinate_species=species_2_name,
        parent_category="open_in_target",
        parent_count=s1_open_total,
    )
    add_row(
        category="open_in_target_promoters",
        path=s2_open_promoters,
        source_species=species_2_name,
        target_species=species_1_name,
        coordinate_species=species_1_name,
        parent_category="open_in_target",
        parent_count=s2_open_total,
    )
    add_row(
        category="open_in_target_enhancers",
        path=s2_open_enhancers,
        source_species=species_2_name,
        target_species=species_1_name,
        coordinate_species=species_1_name,
        parent_category="open_in_target",
        parent_count=s2_open_total,
    )

    return rows


def write_summary_csv(rows: List[dict], output_csv: Path) -> None:
    """Write summary rows to CSV.

    Args:
        rows: Summary row dictionaries.
        output_csv: Output CSV path.
    """
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    headers = [
        "category",
        "source_species",
        "target_species",
        "coordinate_species",
        "count",
        "parent_category",
        "percent_of_parent",
        "file_path",
    ]

    with open(output_csv, "w") as fout:
        fout.write(",".join(headers) + "\n")
        for row in rows:
            fout.write(",".join(str(row[h]) for h in headers) + "\n")


def print_key_summary(rows: List[dict], config: dict) -> None:
    """Print a compact human-readable summary of key BEDTools results.

    Args:
        rows: Summary row dictionaries.
        config: Pipeline configuration dictionary.
    """
    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]
    organ = config["comparison"]["organ"]

    print("\nKey BEDTools summary")
    print("-" * 80)
    print(f"Organ: {organ}")
    print(f"Species pair: {species_1_name} vs {species_2_name}")
    print()

    # Native peak composition
    for species_name in [species_1_name, species_2_name]:
        native_peaks = get_row(rows, "native_peaks", species_name, species_name)
        native_promoters = get_row(rows, "native_promoters", species_name, species_name)
        native_enhancers = get_row(rows, "native_enhancers", species_name, species_name)

        if native_peaks:
            print(f"{species_name} native peaks: {native_peaks['count']}")
        if native_promoters:
            pct = native_promoters["percent_of_parent"]
            print(f"  promoters: {native_promoters['count']} ({pct}% of native peaks)")
        if native_enhancers:
            pct = native_enhancers["percent_of_parent"]
            print(f"  enhancers: {native_enhancers['count']} ({pct}% of native peaks)")
        if native_peaks or native_promoters or native_enhancers:
            print()

    # Cross-species open/closed and open promoter/enhancer composition
    comparisons = [
        (species_1_name, species_2_name),
        (species_2_name, species_1_name),
    ]

    for source_species, target_species in comparisons:
        open_row = get_row(rows, "open_in_target", source_species, target_species)
        closed_row = get_row(rows, "closed_in_target", source_species, target_species)
        open_promoters = get_row(rows, "open_in_target_promoters", source_species, target_species)
        open_enhancers = get_row(rows, "open_in_target_enhancers", source_species, target_species)

        if not any([open_row, closed_row, open_promoters, open_enhancers]):
            continue

        print(f"{source_species} peaks in {target_species}:")
        if open_row:
            print(f"  open in target: {open_row['count']}")
        if closed_row:
            print(f"  closed in target: {closed_row['count']}")
        if open_promoters:
            pct = open_promoters["percent_of_parent"]
            print(f"  open peaks that are promoter-like: {open_promoters['count']} ({pct}% of open peaks)")
        if open_enhancers:
            pct = open_enhancers["percent_of_parent"]
            print(f"  open peaks that are enhancer-like: {open_enhancers['count']} ({pct}% of open peaks)")
        print()

    print("-" * 80)


def run_bedtools_summary(config: dict) -> Path:
    """Create the final BEDTools summary CSV.

    Args:
        config: Pipeline configuration dictionary.

    Returns:
        Path to the generated summary CSV.
    """
    verbose = bool(config.get("project", {}).get("verbose", False))

    summary_dir = Path(config["project"]["output_dir"]) / "bedtools" / "summary"
    helpers.ensure_dir(summary_dir)

    output_csv = summary_dir / "bedtools_summary.csv"

    helpers.vprint(verbose, "Building BEDTools summary rows")
    rows = build_summary_rows(config)

    write_summary_csv(rows, output_csv)
    print_key_summary(rows, config)

    print(f"Wrote: {output_csv}")
    print("BEDTools summary complete")
    return output_csv