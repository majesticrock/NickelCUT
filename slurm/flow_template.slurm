#!/bin/bash
#SBATCH --job-name=flow               # will be overwritten by the generator
#SBATCH --output=/home/althueser/phd/cpp/NickelCUT/output_long.txt   # overwritten
#SBATCH --time=48:00:00               # (long) – overwritten for med queue
#SBATCH --partition=long              # overwritten for med queue
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0                       # all RAM on the node
#SBATCH --hint=nomultithread
#SBATCH --constraint="CascadeLake|IceLake"   # accept both micro‑architectures
#SBATCH --mail-user=joshua.althueser@tu-dortmund.de
#SBATCH --mail-type=END

# -------------------------------------------------------
# 1.  Move to the project root (same for long & med)
# -------------------------------------------------------
cd /home/althueser/phd/cpp/NickelCUT/ || { echo "Cannot cd to project dir"; exit 1; }

# Lattice size is filled in by the batch-script generator.
L=__LATTICE_SIZE__

# -------------------------------------------------------
# 2.  Helper: pick the correct binary at *runtime*
# -------------------------------------------------------
CPU_INFO="$(awk -F: '
    /^(vendor_id|cpu family|model)[[:space:]]*:/ {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
        if ($1 ~ /^vendor_id/) vendor = $2
        else if ($1 ~ /^cpu family/) family = $2
        else if ($1 ~ /^model/) model = $2
    }
    END { print vendor, family, model }
' /proc/cpuinfo)"
read -r CPU_VENDOR CPU_FAMILY CPU_MODEL <<< "$CPU_INFO"

case "${CPU_VENDOR}:${CPU_FAMILY}:${CPU_MODEL}" in
    GenuineIntel:6:85)  ARCH_DIR="cascadelake" ;;
    GenuineIntel:6:106) ARCH_DIR="icelake" ;;
    *)
        echo "ERROR: Unsupported CPU (vendor=$CPU_VENDOR family=$CPU_FAMILY model=$CPU_MODEL)"
        exit 1
        ;;
esac

# L is encoded in the build path; the generator writes the correct
# directory name into the placeholder below.
BUILD_DIR="./build/${ARCH_DIR}/L${L}"
EXECUTABLE="${BUILD_DIR}/apps/flow"

if [[ ! -x "${EXECUTABLE}" ]]; then
    echo "ERROR: Executable ${EXECUTABLE} not found or not executable."
    exit 1
fi

# -------------------------------------------------------
# 3.  OpenMP – use *all* physical cores allocated on the node
# -------------------------------------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
echo "Running on $(hostname) with $OMP_NUM_THREADS threads (CPU model $CPU_MODEL, arch $ARCH_DIR, L=$L)"

# -------------------------------------------------------
# 4.  Run the program – the actual config file is supplied
#    at generation time.
# -------------------------------------------------------
"$EXECUTABLE" "__CONFIG_FILE__"