#!/bin/bash
#SBATCH --job-name=flow
#SBATCH --output=/home/althueser/phd/cpp/NickelCUT/output_long.txt
#SBATCH --time=48:00:00
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=40G
#SBATCH --hint=nomultithread
#SBATCH --constraint="CascadeLake|IceLake"
#SBATCH --mail-user=joshua.althueser@tu-dortmund.de
#SBATCH --mail-type=END

# -------------------------------------------------
# 0. Go to the project directory
# -------------------------------------------------
cd /home/althueser/phd/cpp/NickelCUT/ || { echo "Cannot cd!"; exit 1; }

date
echo "--- START ---"

# -------------------------------------------------
# 1. Determine which CPU family we got
# -------------------------------------------------
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
    *) echo "ERROR: Unsupported CPU (vendor=$CPU_VENDOR family=$CPU_FAMILY model=$CPU_MODEL)"; exit 1 ;;
esac

# -------------------------------------------------
# 2. Pick the correct binary (adjust L‑folder if needed)
# -------------------------------------------------
LAYER="L1"                       # <-- adapt if you have several builds
EXEC="./build/${ARCH_DIR}/${LAYER}/apps/flow"

if [[ ! -x "$EXEC" ]]; then
    echo "ERROR: Executable $EXEC not found."
    exit 1
fi

# -------------------------------------------------
# 3. OpenMP threads = all allocated physical cores
# -------------------------------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
echo "Running on $(hostname) – $OMP_NUM_THREADS cores (CPU model $CPU_MODEL, arch $ARCH_DIR)"
echo "Executable = $EXEC"

# -------------------------------------------------
# 4. Launch the program
# -------------------------------------------------
"$EXEC" params/long_queue.config

echo "--- END ---"
date
echo "$(whoami) is leaving from $(hostname) ..."
echo