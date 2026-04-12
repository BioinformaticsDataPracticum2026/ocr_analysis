from pathlib import Path
import shutil
from typing import Optional, Union


def environment_help() -> str:
    """Return a short message for restoring the expected cluster environment.

    Returns:
        A human-readable setup message suggesting module loading and virtual
        environment activation.
    """
    return (
        "You may be in a fresh cluster login session.\n"
        "Recommended setup:\n"
        "  module load python\n"
        "  source .venv/bin/activate\n"
    )


def check_python_package(import_name: str, install_name: Optional[str] = None) -> None:
    """Exit if a required Python package cannot be imported.

    Args:
        import_name: Module name used in Python import statements.
        install_name: Package name the user should install with pip. If not
            provided, ``import_name`` is used.

    Raises:
        SystemExit: If the package is missing.
    """
    try:
        __import__(import_name)
    except ModuleNotFoundError:
        if install_name is None:
            install_name = import_name
        raise SystemExit(
            f"Missing Python package: {install_name}\n\n"
            f"{environment_help()}\n"
            f"Or install just this package:\n"
            f"  python -m pip install {install_name}"
        )


def check_executable(executable_name: str) -> None:
    """Exit if a required executable is not available on PATH.

    Args:
        executable_name: Name of the executable to look up.

    Raises:
        SystemExit: If the executable is not found.
    """
    if shutil.which(executable_name) is None:
        raise SystemExit(
            f"Missing executable on PATH: {executable_name}\n\n"
            f"{environment_help()}"
        )


def check_file_exists(path_str: Union[str, Path], label: str) -> None:
    """Exit if a required file does not exist.

    Args:
        path_str: Path to the required file.
        label: Human-readable label describing the file.

    Raises:
        SystemExit: If the file does not exist.
    """
    path = Path(path_str)
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def check_dependencies(skip_halper: bool = False) -> None:
    """Validate required Python packages and executables.

    This function checks core Python dependencies first. It always checks for
    BEDTools because downstream analysis depends on it. It only checks for
    HALPER-related executables when HALPER is not being skipped.

    Args:
        skip_halper: Whether the HALPER step will be skipped.

    Raises:
        SystemExit: If any required package or executable is missing.
    """
    check_python_package("yaml", "pyyaml")
    check_python_package("numpy", "numpy")
    check_python_package("matplotlib", "matplotlib")

    check_executable("bedtools")

    if not skip_halper:
        check_executable("halLiftover")


def check_config_dependencies(config: dict, skip_halper: bool = False) -> None:
    """Validate required config-driven input files.

    When HALPER is enabled, this checks the HAL file, orthologFind script, and
    source peak files needed for mapping. When HALPER is skipped, this instead
    checks that the expected HALPER output BED files already exist. In all
    cases, it checks the native peak files and TSS annotation files needed for
    downstream BEDTools steps.

    Args:
        config: Parsed pipeline configuration dictionary.
        skip_halper: Whether the HALPER step will be skipped.

    Raises:
        SystemExit: If any required configured file is missing.
        KeyError: If expected config keys are missing.
    """
    check_file_exists(config["peaks"]["species_1_peak_file"], "species 1 peak file")
    check_file_exists(config["peaks"]["species_2_peak_file"], "species 2 peak file")
    check_file_exists(
        config["annotations"]["species_1_tss_file"],
        "species 1 TSS file",
    )
    check_file_exists(
        config["annotations"]["species_2_tss_file"],
        "species 2 TSS file",
    )

    if skip_halper:
        check_file_exists(
            config["halper_outputs"]["species_1_to_species_2"],
            "species 1 to species 2 HALPER output BED",
        )

        bidirectional = bool(
            config.get("parameters", {}).get("run_bidirectional_liftover", False)
        )
        if bidirectional:
            check_file_exists(
                config["halper_outputs"]["species_2_to_species_1"],
                "species 2 to species 1 HALPER output BED",
            )
    else:
        check_file_exists(config["alignments"]["hal_file"], "HAL file")
        check_file_exists(config["tools"]["orthologfind_py"], "orthologFind.py")