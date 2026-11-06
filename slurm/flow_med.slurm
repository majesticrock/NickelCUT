#!/bin/bash
#SBATCH --job-name=flow
#SBATCH --output=/home/althueser/phd/cpp/NickelCUT/output_med.txt
#SBATCH --time=8:00:00
#SBATCH --partition=med
#SBATCH --nodes=1                 # exactly one node
#SBATCH --exclusive               # give us the whole node
#SBATCH --mem=40G
#SBATCH --hint=nomultithread      # deactivate hyper‑threading (optional but recommended)
#SBATCH --constraint="CascadeLake|IceLake"   # accept either family
#SBATCH --mail-user=joshua.althueser@tu-dortmund.de
#SBATCH --mail-type=END

# --------------------------------------------------------------------
# 1.  Go to the project directory
# --------------------------------------------------------------------
cd /home/althueser/phd/cpp/NickelCUT/ || { echo "Cannot cd to project dir"; exit 1; }

date
echo "--- START ---"

# --------------------------------------------------------------------
# 2.  Find out which CPU family the allocated node belongs to
# --------------------------------------------------------------------
echo "SLURM reports constraint: $SLURM_CONSTRAINT"

case "$SLURM_CONSTRAINT" in
    CascadeLake)
        ARCH_DIR="cascadelake"
        ;;
    IceLake)
        ARCH_DIR="icelake"
        ;;
    *)
        echo "ERROR: Unexpected CPU family '$SLURM_CONSTRAINT'."
        exit 1
        ;;
esac

# --------------------------------------------------------------------
# 3.  Choose the right executable (the L‑folder is identical for both builds)
# --------------------------------------------------------------------
# If you have several L‑folders (e.g. L1, L2, …) change the value of LAYER
LAYER="L1"                                    # <-- adapt if necessary
EXEC="./build/${ARCH_DIR}/${LAYER}/apps/flow"

if [[ ! -x "$EXEC" ]]; then
    echo "ERROR: Executable $EXEC not found or not executable."
    exit 1
fi

# --------------------------------------------------------------------
# 4.  Set OpenMP thread count to the *exact* number of cores allocated
# --------------------------------------------------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
echo "Running $EXEC with OMP_NUM_THREADS=$OMP_NUM_THREADS on $(hostname)"

# --------------------------------------------------------------------
# 5.  Launch the program
# --------------------------------------------------------------------
"$EXEC" params/med_queue.config

echo "--- END ---"
date
echo "$(whoami) is leaving from $(hostname) ..."
echo