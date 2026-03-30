from pathlib import Path
import shutil


def check_python_package(import_name: str, install_name: str = None) -> None:
    try:
        __import__(import_name)
    except ModuleNotFoundError:
        if install_name is None:
            install_name = import_name
        raise SystemExit(
            f"Missing Python package: {install_name}\n"
            f"Try:\n"
            f"  python -m pip install {install_name}"
        )


def check_executable(executable_name: str) -> None:
    if shutil.which(executable_name) is None:
        raise SystemExit(
            f"Missing executable on PATH: {executable_name}"
        )


def check_file_exists(path_str: str, label: str) -> None:
    path = Path(path_str)
    if not path.exists():
        raise SystemExit(f"Missing {label}: {path}")


def check_dependencies(config: dict) -> None:
    # Python packages
    check_python_package("yaml", "pyyaml")
    check_python_package("numpy", "numpy")
    check_python_package("matplotlib", "matplotlib")

    # Executables
    check_executable("halLiftover")

    # Important files
    check_file_exists(config["alignments"]["hal_file"], "HAL file")
    check_file_exists(config["tools"]["orthologfind_py"], "orthologFind.py")
    check_file_exists(config["peaks"]["species_1_peak_file"], "species 1 peak file")
    check_file_exists(config["peaks"]["species_2_peak_file"], "species 2 peak file")