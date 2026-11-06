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

# -------------------------------------------------------
# 2.  Helper: pick the correct binary at *runtime*
# -------------------------------------------------------
case "$SLURM_CONSTRAINT" in
    CascadeLake) ARCH_DIR="cascadelake" ;;
    IceLake)    ARCH_DIR="icelake"    ;;
    *)
        echo "ERROR: Unexpected SLURM constraint '$SLURM_CONSTRAINT'"
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
echo "Running on $(hostname) with $OMP_NUM_THREADS threads (arch $ARCH_DIR, L=$L)"

# -------------------------------------------------------
# 4.  Run the program – the actual config file is supplied
#    as the *first* argument by the generator.
# -------------------------------------------------------
"$EXECUTABLE" "$1"