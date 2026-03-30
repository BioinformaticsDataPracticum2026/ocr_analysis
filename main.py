#!/usr/bin/env python3

from pathlib import Path
import argparse
import sys

def check_dependencies() -> None:
    try:
        import yaml
    except ModuleNotFoundError:
        project_root = Path(__file__).resolve().parent
        venv_python = project_root / ".venv" / "bin" / "python"

        message = f"""
Missing dependency: PyYAML

This project needs the package:
  pyyaml   (imported in Python as: yaml)

Recommended setup:
  python3 -m venv .venv
  source .venv/bin/activate
  python -m pip install pyyaml

Then run:
  python main.py --config config/config.yaml

If you already created the virtual environment, make sure you activated it first.

Expected venv python:
  {venv_python}
"""
        raise SystemExit(message.strip())

def load_config(config_path: Path) -> dict:
    import yaml

    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def run_halper(config: dict) -> None:
    print("Running HALPER...")
    print(f"  species 1: {config['comparison']['species_1']}")
    print(f"  species 2: {config['comparison']['species_2']}")
    print(f"  organ: {config['comparison']['organ']}")
    print(f"  hal file: {config['alignments']['hal_file']}")
    print(f"  species 1 peaks: {config['peaks']['species_1_peak_file']}")
    print(f"  species 2 peaks: {config['peaks']['species_2_peak_file']}")
    # HALPER code goes here
    # call halper with halper.py

def run_bedtools(config: dict) -> None:
    print("Running bedtools annotation...")
    print(f"  species 1 TSS: {config['genome_info']['species_1_tss_file']}")
    print(f"  species 2 TSS: {config['genome_info']['species_2_tss_file']}")
    # bedtools code goes here
    # call bedtools with bedtools.py

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
        "--skip-bedtools",
        action="store_true",
        help="Skip bedtools step",
    )
    args = parser.parse_args()

    check_dependencies()
    config = load_config(args.config)

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

    if not args.skip_bedtools:
        run_bedtools(config)
    else:
        print("Skipping bedtools")

    print("=" * 80)
    print("Pipeline complete")


if __name__ == "__main__":
    main()
