#!/bin/bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_CONFIG="${PROJECT_ROOT}/params/med_queue.config"
INPUT_FILE="${PROJECT_ROOT}/params/for_auto.txt"
SLURM_TEMPLATE="${PROJECT_ROOT}/slurm/flow_med.slurm"
BUILD_DIR="${BUILD_DIR:-${PROJECT_ROOT}/build/cascadelake}"

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
  [3]="E_F"
  [E]="E_F"
  [e]="E_F"
  [4]="tprime"
  [p]="tprime"
  [P]="tprime"
)

echo "Select the parameter to vary:"
echo "1 / U: local interaction"
echo "2 / T: hopping"
echo "3 / E: Fermi energy"
echo "4 / p: next-nearest hopping"
read -p "Enter your choice: " choice
TOKEN="${TOKENS[$choice]}"

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
      -e "s|#SBATCH --output=/home/althueser/phd/cpp/NickelCUT/output_med.txt|#SBATCH --output=${PROJECT_ROOT}/output_${CURRENT_TIME}_${NEW_NAME}.txt|" \
      -e "s|^#SBATCH --partition=.*|#SBATCH --partition=med|" \
      -e "s|./build/cascadelake/flow params/med_queue.config|${BUILD_DIR}/flow ${OUTPUT_DIR}/${NEW_NAME}.config|" \
      "${SLURM_TEMPLATE}" > "${slurm_path}"

  echo "Submitting ${slurm_path}"
  #sbatch "${slurm_path}"
done
