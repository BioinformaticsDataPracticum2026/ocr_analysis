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
        "--skip-bedtools-preprocess",
        action="store_true",
        help="Skip BEDTools preprocessing step",
    )
    args = parser.parse_args()

    # check Python modules + basic executables first
    check_dependencies()

    # import only after dependency check passes
    from utils.config import load_config
    from scripts.halper import run_halper
    from scripts.bedtools_preprocess import run_bedtools_preprocess
    from utils.check_dependencies import check_config_dependencies

    config = load_config(args.config)

    # now check files/paths that depend on config
    check_config_dependencies(config)

    print("=" * 80)
    print("Starting pipeline")
    print(f"Python executable: {sys.executable}")
    print(f"Config file: {args.config}")
    print("=" * 80)

    if not args.skip_halper:
        print("\n[1/2] Running HALPER")
        run_halper(config)
    else:
        print("\n[1/2] Skipping HALPER")

    if not args.skip_bedtools_preprocess:
        print("\n[2/2] Running BEDTools preprocessing")
        run_bedtools_preprocess(config)
    else:
        print("\n[2/2] Skipping BEDTools preprocessing")

    print("\n" + "=" * 80)
    print("Pipeline complete")


if __name__ == "__main__":
    main()