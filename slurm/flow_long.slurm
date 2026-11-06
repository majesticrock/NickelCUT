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
echo "SLURM constraint reported: $SLURM_CONSTRAINT"
case "$SLURM_CONSTRAINT" in
    CascadeLake) ARCH_DIR="cascadelake" ;;
    IceLake)    ARCH_DIR="icelake"    ;;
    *) echo "ERROR: Unexpected CPU family"; exit 1 ;;
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
echo "Running on $(hostname) – $OMP_NUM_THREADS cores"
echo "Executable = $EXEC"

# -------------------------------------------------
# 4. Launch the program
# -------------------------------------------------
"$EXEC" params/long_queue.config

echo "--- END ---"
date
echo "$(whoami) is leaving from $(hostname) ..."
echo