#!/usr/bin/env python3
"""Run mean_field for matching flow results found in NickelCUT's data tree."""

import argparse
import math
import re
import shlex
import subprocess
import tempfile
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_DATA_ROOT = (PROJECT_ROOT / "../../data/nickel_cut").resolve()
EXTRACTED_CHANNELS = "lowest_ROD_extraced_channels.bin"


def read_config(path):
    values = {}
    lines = path.read_text().splitlines(keepends=True)
    for line_number, line in enumerate(lines, start=1):
        fields = shlex.split(line, comments=True)
        if not fields:
            continue
        if fields[0] in values:
            raise ValueError(f"{path}:{line_number}: duplicate parameter {fields[0]}")
        if len(fields) < 2:
            raise ValueError(f"{path}:{line_number}: missing value for {fields[0]}")
        values[fields[0]] = fields[1]

    for key in ("U_0", "T", "tprime", "E_F", "output_dir"):
        if key not in values:
            raise ValueError(f"{path}: required parameter {key} is missing")
    for key in ("U_0", "T", "tprime"):
        try:
            values[key] = float(values[key])
        except ValueError as error:
            raise ValueError(f"{path}: {key} must be numeric") from error
    return values, lines


def lattice_size():
    header = (PROJECT_ROOT / "sources/L.hpp").read_text()
    match = re.search(r"static constexpr int L\s*=\s*(\d+)\s*;", header)
    if not match:
        raise ValueError("Could not determine L from sources/L.hpp")
    return int(match.group(1))


def same_number(path_value, target):
    try:
        return math.isclose(float(path_value), target, rel_tol=1e-12, abs_tol=1e-12)
    except ValueError:
        return False


def matching_flow_results(binary_root, parameters, size):
    results = []
    for state_file in binary_root.rglob(EXTRACTED_CHANNELS):
        components = {}
        for part in state_file.parent.parts:
            if "=" in part:
                key, value = part.split("=", 1)
                components[key] = value

        if components.get("L") != str(size):
            continue
        if not all(
            same_number(components.get(key, ""), parameters[key])
            for key in ("T", "U_0", "tprime")
        ):
            continue
        if state_file.stat().st_size == 0:
            continue
        try:
            fermi_energy = float(components["E_F"])
        except (KeyError, ValueError):
            continue
        results.append((fermi_energy, state_file))

    return sorted(results, key=lambda result: result[0])


def config_for_fermi_energy(lines, fermi_energy):
    updated = []
    replacements = 0
    for line in lines:
        if shlex.split(line, comments=True)[:1] == ["E_F"]:
            line, count = re.subn(
                r"^(\s*E_F\s+)\S+",
                lambda match: f"{match.group(1)}{fermi_energy}",
                line,
                count=1,
            )
            replacements += count
        updated.append(line)
    if replacements != 1:
        raise ValueError("Expected exactly one E_F entry in the config")
    return "".join(updated)


def main():
    parser = argparse.ArgumentParser(
        description="Run mean_field for every matching saved NickelCUT flow result."
    )
    parser.add_argument(
        "config",
        nargs="?",
        type=Path,
        default=PROJECT_ROOT / "params/mean_field_basic.config",
        help="input parameter file (default: params/mean_field_basic.config)",
    )
    parser.add_argument(
        "--mean-field",
        type=Path,
        default=PROJECT_ROOT / "build/default/apps/mean_field",
        help="mean_field executable (default: build/default/apps/mean_field)",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help="compiled NickelCUT data root (default: ../../data/nickel_cut)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="list matching flow results without launching mean_field",
    )
    args = parser.parse_args()

    config_path = args.config.resolve()
    data_root = args.data_root.resolve()
    parameters, lines = read_config(config_path)
    binary_root = data_root / parameters["output_dir"] / "binaries"
    matches = matching_flow_results(binary_root, parameters, lattice_size())

    if not matches:
        print(f"No matching flow results found under {binary_root}")
        return 0

    if not args.dry_run and not args.mean_field.is_file():
        parser.error(f"mean_field executable not found: {args.mean_field}")

    failures = 0
    with tempfile.TemporaryDirectory(prefix="nickelcut_mean_field_") as temp_dir:
        for fermi_energy, state_file in matches:
            print(f"E_F={fermi_energy:g}: {state_file}")
            if args.dry_run:
                continue

            run_config = Path(temp_dir) / f"E_F_{fermi_energy:g}.config"
            run_config.write_text(config_for_fermi_energy(lines, fermi_energy))
            completed = subprocess.run(
                [str(args.mean_field.resolve()), str(run_config)],
                cwd=PROJECT_ROOT,
                check=False,
            )
            if completed.returncode != 0:
                print(f"mean_field failed for E_F={fermi_energy:g} (exit {completed.returncode})")
                failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())