#!/usr/bin/env python3

from typing import Dict
from pathlib import Path

import utils.helpers as helpers
from scripts.promoter_enhancer import classify_promoter_enhancer


def run_cross_species_ep(config: dict, open_closed_outputs: Dict[str, Path]) -> Dict[str, Path]:
    verbose = bool(config.get("project", {}).get("verbose", False))

    helpers.require_executable("bedtools", "BEDTools")

    results_dir = Path(config["project"]["output_dir"]) / "bedtools" / "cross_species_ep"
    tmp_dir = results_dir / "tmp"
    helpers.ensure_dir(results_dir)
    helpers.ensure_dir(tmp_dir)

    organ = config["comparison"]["organ"]
    species_1_name = config["comparison"]["species_1_name"]
    species_2_name = config["comparison"]["species_2_name"]

    species_1_tss_bed = Path(config["annotations"]["species_1_tss_file"])
    species_2_tss_bed = Path(config["annotations"]["species_2_tss_file"])

    promoter_max_distance = int(
        config.get("annotations", {}).get("promoter_max_distance", 5000)
    )

    outputs: Dict[str, Path] = {}

    # species 1 mapped into species 2, then classified with species 2 TSS
    s1_open_in_s2 = open_closed_outputs.get("species_1_open_in_species_2")
    if s1_open_in_s2 and s1_open_in_s2.exists():
        annotated = tmp_dir / f"{species_1_name}_open_in_{species_2_name}_tss_annotated.bed"
        promoters = results_dir / f"{species_1_name}_open_in_{species_2_name}_promoters.bed"
        enhancers = results_dir / f"{species_1_name}_open_in_{species_2_name}_enhancers.bed"

        print(f"Classifying {species_1_name} peaks open in {species_2_name} as promoter-like or enhancer-like")
        classify_promoter_enhancer(
            peak_bed=s1_open_in_s2,
            tss_bed=species_2_tss_bed,
            annotated_output=annotated,
            promoter_output=promoters,
            enhancer_output=enhancers,
            promoter_max_distance=promoter_max_distance,
            verbose=verbose,
        )

        outputs["species_1_open_in_species_2_tss_annotated"] = annotated
        outputs["species_1_open_in_species_2_promoters"] = promoters
        outputs["species_1_open_in_species_2_enhancers"] = enhancers
        print(f"Wrote: {promoters}")
        print(f"Wrote: {enhancers}")

    # species 2 mapped into species 1, then classified with species 1 TSS
    s2_open_in_s1 = open_closed_outputs.get("species_2_open_in_species_1")
    if s2_open_in_s1 and s2_open_in_s1.exists():
        annotated = tmp_dir / f"{species_2_name}_open_in_{species_1_name}_tss_annotated.bed"
        promoters = results_dir / f"{species_2_name}_open_in_{species_1_name}_promoters.bed"
        enhancers = results_dir / f"{species_2_name}_open_in_{species_1_name}_enhancers.bed"

        print(f"Classifying {species_2_name} peaks open in {species_1_name} as promoter-like or enhancer-like")
        classify_promoter_enhancer(
            peak_bed=s2_open_in_s1,
            tss_bed=species_1_tss_bed,
            annotated_output=annotated,
            promoter_output=promoters,
            enhancer_output=enhancers,
            promoter_max_distance=promoter_max_distance,
            verbose=verbose,
        )

        outputs["species_2_open_in_species_1_tss_annotated"] = annotated
        outputs["species_2_open_in_species_1_promoters"] = promoters
        outputs["species_2_open_in_species_1_enhancers"] = enhancers
        print(f"Wrote: {promoters}")
        print(f"Wrote: {enhancers}")

    print("Cross-species promoter/enhancer classification complete")
    return outputs