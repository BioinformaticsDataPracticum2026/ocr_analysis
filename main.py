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
    parser.add_argument(
        "--rgreat-only",
        action="store_true",
        help="Skip directly to rGREAT using existing BED outputs",
    )
    args = parser.parse_args()

    # if running rGREAT only, HALPER is implicitly skipped.
    skip_halper = args.skip_halper or args.rgreat_only

    # check Python modules + basic executables first
    check_dependencies(skip_halper=skip_halper)

    # import only after dependency check passes
    from utils.config import load_config
    from utils.check_dependencies import check_config_dependencies
    from scripts.halper import run_halper
    from scripts.bedtools_preprocess import run_bedtools_preprocess
    from scripts.open_closed import run_open_closed
    from scripts.promoter_enhancer import run_promoter_enhancer
    from scripts.cross_species_ep import run_cross_species_ep
    from scripts.bedtools_summary import run_bedtools_summary
    from scripts.rgreat import run_rgreat

    config = load_config(args.config)

    # now check files/paths that depend on config
    check_config_dependencies(config, skip_halper=skip_halper)

    print("=" * 80)
    print("Starting pipeline")
    print(f"Python executable: {sys.executable}")
    print(f"Config file: {args.config}")
    print("=" * 80)

    if args.rgreat_only:
        print("\n[1/1] Running rGREAT only")
        run_rgreat(config)
        print("\n" + "=" * 80)
        print("Pipeline complete")
        return

    if not args.skip_halper:
        print("\n[1/7] Running HALPER")
        run_halper(config)
    else:
        print("\n[1/7] Skipping HALPER")

    print("\n[2/7] Running BEDTools preprocessing")
    run_bedtools_preprocess(config)

    print("\n[3/7] Running open/closed classification")
    open_closed_outputs = run_open_closed(config)

    print("\n[4/7] Running promoter/enhancer classification")
    run_promoter_enhancer(config)

    print("\n[5/7] Running cross-species promoter/enhancer classification")
    run_cross_species_ep(config, open_closed_outputs)

    print("\n[6/7] Writing BEDTools summary")
    run_bedtools_summary(config)

    print("\n[7/7] Running rGREAT")
    run_rgreat(config)

    print("\n" + "=" * 80)
    print("Pipeline complete")


if __name__ == "__main__":
    main()