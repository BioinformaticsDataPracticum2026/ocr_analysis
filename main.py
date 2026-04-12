#!/usr/bin/env python3

from pathlib import Path
import argparse
import sys

from utils.check_dependencies import check_dependencies


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run cross-species OCR pipeline"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("./config.yaml"),
        help="Path to config YAML file",
    )
    parser.add_argument(
        "--skip-halper",
        action="store_true",
        help="Skip HALPER step",
    )
    args = parser.parse_args()

    # check Python modules + basic executables first
    check_dependencies(skip_halper=args.skip_halper)

    # import only after dependency check passes
    from utils.config import load_config
    from utils.check_dependencies import check_config_dependencies
    from scripts.halper import run_halper
    from scripts.bedtools_preprocess import run_bedtools_preprocess
    from scripts.open_closed import run_open_closed
    from scripts.promoter_enhancer import run_promoter_enhancer
    from scripts.cross_species_ep import run_cross_species_ep

    config = load_config(args.config)

    # now check files/paths that depend on config
    check_config_dependencies(config, skip_halper=args.skip_halper)

    print("=" * 80)
    print("Starting pipeline")
    print(f"Python executable: {sys.executable}")
    print(f"Config file: {args.config}")
    print("=" * 80)

    if not args.skip_halper:
        print("\n[1/5] Running HALPER")
        run_halper(config)
    else:
        print("\n[1/5] Skipping HALPER")

    print("\n[2/5] Running BEDTools preprocessing")
    run_bedtools_preprocess(config)

    print("\n[3/5] Running open/closed classification")
    open_closed_outputs = run_open_closed(config)

    print("\n[4/5] Running promoter/enhancer classification")
    run_promoter_enhancer(config)

    print("\n[5/5] Running cross-species promoter/enhancer classification")
    run_cross_species_ep(config, open_closed_outputs)

    

    print("\n" + "=" * 80)
    print("Pipeline complete")


if __name__ == "__main__":
    main()