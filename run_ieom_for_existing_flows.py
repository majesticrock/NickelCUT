#!/usr/bin/env python3
"""Run susceptibility_ieom for matching flow results in NickelCUT's data tree."""

import argparse
import math
import re
import shlex
import subprocess
import tempfile
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_DATA_ROOT = (PROJECT_ROOT / "../../data/nickel_cut").resolve()
EXTRACTED_CHANNELS = "final_flow_state.bin"

# Parameters that can be matched against the directory hierarchy and varied.
NUMERIC_PARAMETERS = ("U_0", "T", "tprime", "E_F")
REQUIRED_PARAMETERS = (*NUMERIC_PARAMETERS, "output_dir")


def read_config(path):
    values = {}
    lines = path.read_text().splitlines(keepends=True)

    for line_number, line in enumerate(lines, start=1):
        fields = shlex.split(line, comments=True)
        if not fields:
            continue

        if fields[0] in values:
            raise ValueError(
                f"{path}:{line_number}: duplicate parameter {fields[0]}"
            )
        if len(fields) < 2:
            raise ValueError(
                f"{path}:{line_number}: missing value for {fields[0]}"
            )

        values[fields[0]] = fields[1]

    for key in REQUIRED_PARAMETERS:
        if key not in values:
            raise ValueError(f"{path}: required parameter {key} is missing")

    for key in NUMERIC_PARAMETERS:
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
        return math.isclose(
            float(path_value),
            target,
            rel_tol=1e-12,
            abs_tol=1e-12,
        )
    except (TypeError, ValueError):
        return False


def directory_parameters(state_file):
    """Extract key=value parameters from the file's parent directory hierarchy."""
    components = {}

    for part in state_file.parent.parts:
        if "=" in part:
            key, value = part.split("=", 1)
            components[key] = value

    return components


def matching_flow_results(binary_root, parameters, size, varied_parameters):
    """Find results matching every fixed parameter and return their parameter values."""
    results = []

    if not binary_root.is_dir():
        return results

    for state_file in binary_root.rglob(EXTRACTED_CHANNELS):
        if state_file.stat().st_size == 0:
            continue

        components = directory_parameters(state_file)

        if components.get("L") != str(size):
            continue

        # Every numeric parameter that is not varied must match the input config.
        fixed_parameters = (
            set(NUMERIC_PARAMETERS) - set(varied_parameters)
        )

        if not all(
            same_number(components.get(key, ""), parameters[key])
            for key in fixed_parameters
        ):
            continue

        # Each varied parameter must be present and numeric in the directory path.
        try:
            varied_values = {
                key: float(components[key])
                for key in varied_parameters
            }
        except (KeyError, ValueError):
            continue

        results.append((varied_values, state_file))

    # Sort lexicographically by the requested parameter order.
    return sorted(
        results,
        key=lambda result: tuple(
            result[0][key] for key in varied_parameters
        ),
    )


def config_for_parameters(lines, varied_values):
    """Replace the requested parameter values while preserving other config lines."""
    updated = []
    replacements = {key: 0 for key in varied_values}

    for line in lines:
        fields = shlex.split(line, comments=True)
        if not fields or fields[0] not in varied_values:
            updated.append(line)
            continue

        key = fields[0]
        value = varied_values[key]

        # Replace only the value, preserving whitespace, trailing comments,
        # and the rest of the original line.
        pattern = rf"^(\s*{re.escape(key)}\s+)\S+"
        new_line, count = re.subn(
            pattern,
            lambda match: f"{match.group(1)}{value:.17g}",
            line,
            count=1,
        )

        if count != 1:
            raise ValueError(f"Could not update parameter {key}")

        replacements[key] += count
        updated.append(new_line)

    invalid = [key for key, count in replacements.items() if count != 1]
    if invalid:
        raise ValueError(
            "Expected exactly one config entry for: " + ", ".join(invalid)
        )

    return "".join(updated)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Run susceptibility_ieom for saved NickelCUT flow results, "
            "varying one or more parameters."
        )
    )
    parser.add_argument(
        "config",
        nargs="?",
        type=Path,
        default=PROJECT_ROOT / "params/mean_field_basic.config",
        help="input parameter file (default: params/mean_field_basic.config)",
    )
    parser.add_argument(
        "--susceptibility_ieom",
        type=Path,
        default=PROJECT_ROOT / "build/default/apps/susceptibility_ieom",
        help="susceptibility_ieom executable",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help="compiled NickelCUT data root (default: ../../data/nickel_cut)",
    )
    parser.add_argument(
        "--vary",
        nargs="+",
        choices=NUMERIC_PARAMETERS,
        default=["E_F"],
        metavar="PARAMETER",
        help=(
            "parameter(s) to vary. Choices: U_0, T, tprime, E_F "
            "(default: E_F). Example: --vary E_F T"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="list matching flow results without launching susceptibility_ieom",
    )
    args = parser.parse_args()

    # Reject duplicate parameter names to avoid ambiguous output/config handling.
    if len(set(args.vary)) != len(args.vary):
        parser.error("--vary parameters must not be repeated")

    config_path = args.config.resolve()
    data_root = args.data_root.resolve()

    try:
        parameters, lines = read_config(config_path)
        size = lattice_size()
        binary_root = data_root / parameters["output_dir"] / "binaries"

        matches = matching_flow_results(
            binary_root,
            parameters,
            size,
            args.vary,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))

    if not matches:
        print(f"No matching flow results found under {binary_root}")
        return 0

    executable = args.susceptibility_ieom.resolve()
    if not args.dry_run and not executable.is_file():
        parser.error(f"susceptibility_ieom executable not found: {executable}")

    failures = 0

    with tempfile.TemporaryDirectory(
        prefix="nickelcut_susceptibility_ieom_"
    ) as temp_dir:
        for varied_values, state_file in matches:
            description = ", ".join(
                f"{key}={varied_values[key]:g}" for key in args.vary
            )
            print(f"{description}: {state_file}")

            if args.dry_run:
                continue

            # Include all varied values in the filename to prevent collisions.
            filename_parts = [
                f"{key}_{varied_values[key]:.17g}" for key in args.vary
            ]
            run_config = Path(temp_dir) / (
                "_".join(filename_parts) + ".config"
            )

            run_config.write_text(
                config_for_parameters(lines, varied_values)
            )

            completed = subprocess.run(
                [str(executable), str(run_config)],
                cwd=PROJECT_ROOT,
                check=False,
            )

            if completed.returncode != 0:
                print(
                    f"susceptibility_ieom failed for {description} "
                    f"(exit {completed.returncode})"
                )
                failures += 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())