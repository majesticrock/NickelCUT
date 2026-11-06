#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------
# 0.  Resolve paths that are independent of the current working dir
# -----------------------------------------------------------------
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_CONFIG="${PROJECT_ROOT}/params/long_queue.config"
INPUT_FILE="${PROJECT_ROOT}/params/for_auto.txt"
SLURM_TEMPLATE="${PROJECT_ROOT}/slurm/flow_template.slurm"

# -----------------------------------------------------------------
# 1.  Argument handling (optional: Lattice size  L  and parameter)
# -----------------------------------------------------------------
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

# -----------------------------------------------------------------
# 2.  Build‑directory handling – we keep the old “cascadelake” default
# -----------------------------------------------------------------
# The actual binary that will be used is resolved *inside* the SLURM script,
# therefore we only need a *representative* path here to test that the
# executable exists for *at least one* architecture.
BUILD_DIR_DEFAULT="${PROJECT_ROOT}/build/cascadelake/L${L}"
BUILD_DIR="${BUILD_DIR:-${BUILD_DIR_DEFAULT}}"
FLOW_EXECUTABLE="${BUILD_DIR}/apps/flow"

if [[ ! -x "${FLOW_EXECUTABLE}" ]]; then
  echo "WARNING: No flow executable found for L=${L}: ${FLOW_EXECUTABLE}"
  echo "Build it first with e.g.:"
  echo "   make cascadelake L=${L}   # for CascadeLake"
  echo "   make icelake    L=${L}   # for IceLake"
  echo "No jobs were generated or submitted."
  exit 1
fi

# -----------------------------------------------------------------
# 3.  Verify that the binary has been configured for the requested L
# -----------------------------------------------------------------
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

# -----------------------------------------------------------------
# 4.  Sanity checks for the ancillary files we need
# -----------------------------------------------------------------
[[ -f "${BASE_CONFIG}" ]]   || { echo "Missing NickelCUT base config: ${BASE_CONFIG}";   exit 1; }
[[ -f "${SLURM_TEMPLATE}" ]] || { echo "Missing Slurm template: ${SLURM_TEMPLATE}"; exit 1; }
[[ -f "${INPUT_FILE}" ]]    || {
  echo "Missing parameter list: ${INPUT_FILE}"
  echo "Expected one value per line, e.g.:"
  echo "  -0.4"
  echo "  -0.3"
  exit 1
}

readarray -t NEW_VALUES < "${INPUT_FILE}"

# -----------------------------------------------------------------
# 5.  Token table – maps the user‑choice to the name used in the config
# -----------------------------------------------------------------
declare -A TOKENS=(
  [1]="U_0"  [U]="U_0"  [u]="U_0"
  [2]="T"    [T]="T"    [t]="T"
  [3]="filling" [F]="filling" [f]="filling"
  [4]="tprime"  [p]="tprime"   [P]="tprime"
)

# -----------------------------------------------------------------
# 6.  Choose which parameter we want to sweep
# -----------------------------------------------------------------
if (( $# == 2 )); then
  choice="$2"
else
  echo "Select the parameter to vary:"
  echo " 1 / U : local interaction"
  echo " 2 / T : temperature"
  echo " 3 / F : filling"
  echo " 4 / p : next‑nearest hopping"
  read -r -p "Enter the number or the letter: " choice
fi
TOKEN="${TOKENS[$choice]:-}"

if [[ -z "${TOKEN}" ]]; then
  echo "Invalid choice. Exiting."
  exit 1
fi

# -----------------------------------------------------------------
# 7.  Output folder for the *generated* configs & slurm files
# -----------------------------------------------------------------
CURRENT_TIME="$(date +"%Y%m%d_%H%M%S")"
OUTPUT_DIR="${PROJECT_ROOT}/auto_generated_${CURRENT_TIME}"
mkdir -p "${OUTPUT_DIR}"

# -----------------------------------------------------------------
# 8.  Loop over every value in the list and create a job script
# -----------------------------------------------------------------
for NEW_VALUE in "${NEW_VALUES[@]}"; do
  # Strip possible Windows line‑end
  NEW_VALUE="${NEW_VALUE%$'\r'}"
  [[ -z "${NEW_VALUE}" ]] && continue

  # Create a filesystem‑friendly name from the value
  NEW_NAME="$(printf '%s' "${NEW_VALUE}" | tr -cs '[:alnum:]._-' '_')"
  NEW_NAME="${NEW_NAME// /_}"
  NEW_NAME="${NEW_NAME//-/_}"

  # -----------------------------------------------------------------
  # 8a. Produce a *modified* NickelCUT config where the chosen token
  #     gets the new value
  # -----------------------------------------------------------------
  if ! awk -v token="${TOKEN}" -v new_value="${NEW_VALUE}" '
        $1 == token { print token, new_value; found = 1; next }
        { print }
        END { exit found ? 0 : 1 }
      ' "${BASE_CONFIG}" > "${OUTPUT_DIR}/${NEW_NAME}.config"; then
    echo "Parameter ${TOKEN} was not found in ${BASE_CONFIG}."
    exit 1
  fi

  # -----------------------------------------------------------------
  # 8b. Build a *full* SLURM script from the generic template
  # -----------------------------------------------------------------
  slurm_path="${OUTPUT_DIR}/${NEW_NAME}.slurm"

  sed \
      -e "s|#SBATCH --job-name=flow|#SBATCH --job-name=flow_${NEW_NAME}_${CURRENT_TIME}|" \
      -e "s|#SBATCH --output=/home/althueser/phd/cpp/NickelCUT/output_long.txt|#SBATCH --output=${PROJECT_ROOT}/output_${CURRENT_TIME}_${NEW_NAME}.txt|" \
      -e "s|^#SBATCH --partition=.*|#SBATCH --partition=long|" \
      -e "s|#PLACEHOLDER_EXEC|${FLOW_EXECUTABLE} ${OUTPUT_DIR}/${NEW_NAME}.config|" \
      "${SLURM_TEMPLATE}" > "${slurm_path}"

  echo "Submitting ${slurm_path}"
  sbatch "${slurm_path}"
done