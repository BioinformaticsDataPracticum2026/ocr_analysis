#!/usr/bin/env python3

from pathlib import Path
import argparse
import sys

from utils.config import load_config
from utils.check_dependency import check_dependencies
from scripts.halper import run_halper


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

    config = load_config(args.config)
    check_dependencies(config)

    print("=" * 80)
    print("Starting pipeline")
    print(f"Python executable: {sys.executable}")
    print(f"Config file: {args.config}")
    print("=" * 80)

    if not args.skip_halper:
        run_halper(config)
    else:
        print("Skipping HALPER")

    print("=" * 80)
    print("Pipeline complete")


if __name__ == "__main__":
    main()