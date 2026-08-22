#!/usr/bin/env bash
#SBATCH --job-name=bufdeeph-ufd
#SBATCH --partition=normal
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --array=1-9

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$SLURM_SUBMIT_DIR}"
INPUT_DIR="${INPUT_DIR:?Set INPUT_DIR to the DeepH data directory}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$PROJECT_DIR/output/unfolding}"
NUM_PARTS="${NUM_PARTS:-9}"
M_TRANS="${M_TRANS:-3.89}"
N_GVEC="${N_GVEC:-9}"
NUM_BAND="${NUM_BAND:-5200}"
MAX_ITER="${MAX_ITER:-1500}"
NUM_POINTS="${NUM_POINTS:-25}"
FERMI_SHIFT="${FERMI_SHIFT:--2.25}"
MID_Z="${MID_Z:-18.0}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-1}"

PART_INDEX="${SLURM_ARRAY_TASK_ID}"
PART_DIR="$OUTPUT_ROOT/part${PART_INDEX}"
mkdir -p "$PART_DIR"

cd "$PROJECT_DIR"
date > "$PART_DIR/time.dat"

srun --cpu-bind=cores julia --threads="$JULIA_NUM_THREADS" src/band_unfolding.jl \
  --input_dir "$INPUT_DIR" \
  --output_dir "$PART_DIR" \
  --output-file data.txt \
  --m-trans "$M_TRANS" \
  --n-gvec "$N_GVEC" \
  --num-band "$NUM_BAND" \
  --max-iter "$MAX_ITER" \
  --num-points "$NUM_POINTS" \
  --fermi-shift "$FERMI_SHIFT" \
  --mid-z "$MID_Z" \
  --part-index "$PART_INDEX" \
  --num-parts "$NUM_PARTS"

date >> "$PART_DIR/time.dat"
