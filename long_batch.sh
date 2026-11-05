#!/bin/bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_CONFIG="${PROJECT_ROOT}/params/long_queue.config"
INPUT_FILE="${PROJECT_ROOT}/params/for_auto.txt"
SLURM_TEMPLATE="${PROJECT_ROOT}/slurm/flow_long.slurm"

if (( $# > 2 )); then
  echo "Usage: $0 [L [parameter]]"
  exit 1
fi

if (( $# >= 1 )); then
  L="$1"
else
  read -r -p "Lattice size L: " L
fi
if [[ ! "${L}" =~ ^[0-9]+$ ]] || ((10#${L} < 2 || 10#${L} % 2 != 0)); then
  echo "L must be a positive even integer."
  exit 1
fi
L=$((10#${L}))

BUILD_DIR="${BUILD_DIR:-${PROJECT_ROOT}/build/cascadelake/L${L}}"
FLOW_EXECUTABLE="${BUILD_DIR}/apps/flow"

if [[ ! -x "${FLOW_EXECUTABLE}" ]]; then
  echo "WARNING: No flow executable found for L=${L}: ${FLOW_EXECUTABLE}"
  echo "Build it first with: make cascadelake L=${L}"
  echo "No jobs were generated or submitted."
  exit 1
fi

if [[ ! -f "${BUILD_DIR}/CMakeCache.txt" ]]; then
  echo "WARNING: Cannot verify the lattice size of ${FLOW_EXECUTABLE}."
  echo "No jobs were generated or submitted."
  exit 1
fi

BUILT_L="$(sed -n 's/^NICKELCUT_L:STRING=//p' "${BUILD_DIR}/CMakeCache.txt")"
if [[ "${BUILT_L}" != "${L}" ]]; then
  echo "WARNING: ${FLOW_EXECUTABLE} was configured for L=${BUILT_L:-unknown}, not L=${L}."
  echo "No jobs were generated or submitted."
  exit 1
fi

if [[ ! -f "${BASE_CONFIG}" ]]; then
  echo "Missing NickelCUT base config: ${BASE_CONFIG}"
  exit 1
fi

if [[ ! -f "${SLURM_TEMPLATE}" ]]; then
  echo "Missing Slurm template: ${SLURM_TEMPLATE}"
  exit 1
fi

if [[ ! -f "${INPUT_FILE}" ]]; then
  echo "Missing parameter list: ${INPUT_FILE}"
  echo "Expected one value per line, for example:"
  echo "  -0.4"
  echo "  -0.3"
  exit 1
fi

readarray -t NEW_VALUES < "${INPUT_FILE}"

# Parameter names used by NickelCUT configs.
declare -A TOKENS=(
  [1]="U_0"
  [U]="U_0"
  [u]="U_0"
  [2]="T"
  [T]="T"
  [t]="T"
  [3]="filling"
  [F]="filling"
  [f]="filling"
  [4]="tprime"
  [p]="tprime"
  [P]="tprime"
)

if (( $# == 2 )); then
  choice="$2"
else
  echo "Select the parameter to vary:"
  echo "1 / U: local interaction"
  echo "2 / T: temperature"
  echo "3 / F: filling"
  echo "4 / p: next-nearest hopping"
  read -r -p "Enter your choice: " choice
fi
TOKEN="${TOKENS[$choice]:-}"

if [[ -z "${TOKEN}" ]]; then
  echo "Invalid choice. Exiting."
  exit 1
fi

CURRENT_TIME="$(date +"%Y%m%d_%H%M%S")"
OUTPUT_DIR="${PROJECT_ROOT}/auto_generated_${CURRENT_TIME}"
mkdir -p "${OUTPUT_DIR}"

for NEW_VALUE in "${NEW_VALUES[@]}"; do
  NEW_VALUE="${NEW_VALUE%$'\r'}"
  [[ -z "${NEW_VALUE}" ]] && continue

  NEW_NAME="$(printf '%s' "${NEW_VALUE}" | tr -cs '[:alnum:]._-' '_')"
  NEW_NAME="${NEW_NAME// /_}"
  NEW_NAME="${NEW_NAME//-/_}"

  if ! awk -v token="${TOKEN}" -v new_value="${NEW_VALUE}" '
      $1 == token { print token, new_value; found = 1; next }
      { print }
      END { exit found ? 0 : 1 }
    ' "${BASE_CONFIG}" > "${OUTPUT_DIR}/${NEW_NAME}.config"; then
    echo "Parameter ${TOKEN} was not found in ${BASE_CONFIG}."
    exit 1
  fi

  slurm_path="${OUTPUT_DIR}/${NEW_NAME}.slurm"
  sed \
      -e "s|#SBATCH --job-name=flow|#SBATCH --job-name=flow_${NEW_NAME}_${CURRENT_TIME}|" \
      -e "s|#SBATCH --output=/home/althueser/phd/cpp/NickelCUT/output_long.txt|#SBATCH --output=${PROJECT_ROOT}/output_${CURRENT_TIME}_${NEW_NAME}.txt|" \
      -e "s|^#SBATCH --partition=.*|#SBATCH --partition=long|" \
      -e "s|^[[:space:]]*\\./build/cascadelake/apps/flow .*|${FLOW_EXECUTABLE} ${OUTPUT_DIR}/${NEW_NAME}.config|" \
      "${SLURM_TEMPLATE}" > "${slurm_path}"

  echo "Submitting ${slurm_path}"
  sbatch "${slurm_path}"
done
