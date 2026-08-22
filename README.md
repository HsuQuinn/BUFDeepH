# Band Unfolding for Moiré Systems
Author: Hsu Quinn
Version: 3.1.0

BUFDeepH implements band unfolding for honeycomb moire superlattices. The core solver is written in Julia and is designed around sparse DeepH Hamiltonians, overlap matrices, and k-path batch execution on Slurm clusters.

## Features

- Band unfolding along the K->Gamma->M path.
- Sparse Hamiltonian and overlap handling for DeepH output.
- Parameterized moire transformation matrices, band count, Arpack iteration budget, Fermi shift, and layer z cutoff.
- Reuses the preprocessed sparse matrices across all selected k points.
- Supports direct point ranges and Slurm job-array partitioning.

## Requirements

Install the following Julia packages:
```julia
using Pkg
Pkg.add(["DelimitedFiles", "LinearAlgebra", "JSON", "HDF5", "ArgParse", "SparseArrays", "Pardiso", "Arpack", "LinearMaps", "JLD", "PyPlot", "Statistics"])
```

## Usage

Prepare an input directory with DeepH output files, then run:

```bash
julia src/band_unfolding.jl \
  --input_dir ./test/9.43/deephdata \
  --output_dir ./output \
  --m-trans 9.43 \
  --num-band 1500 \
  --max-iter 1000 \
  --num-points 25 \
  --fermi-shift 0.0
```

The default configuration follows the latest 3.89-degree experiment:

```bash
julia src/band_unfolding.jl \
  -i /path/to/deephdata \
  -o ./output/part1 \
  --m-trans 3.89 \
  --n-gvec 9 \
  --num-band 5200 \
  --max-iter 1500 \
  --num-points 25 \
  --fermi-shift -2.25
```

To run only a direct 1-based k-point range:

```bash
julia src/band_unfolding.jl -i /path/to/deephdata -o ./output \
  --point-start 45 --point-end 49
```

To split the k path into partitions:

```bash
julia src/band_unfolding.jl -i /path/to/deephdata -o ./output/part9 \
  --part-index 9 --num-parts 9
```

The script writes `data.txt` by default. Each record is a three-line triplet:

```text
idx_p:<local point index>
energy
weight
```

## Slurm Array

Use `scripts/slurm_band_unfolding_array.sh` for k-path partitioning on a cluster:

```bash
export INPUT_DIR=/path/to/deephdata
export OUTPUT_ROOT=/path/to/output/unfolding_389
export NUM_PARTS=9
export M_TRANS=3.89
sbatch --array=1-9 scripts/slurm_band_unfolding_array.sh
```

For a 9.43-degree case:

```bash
export INPUT_DIR=/path/to/deephdata
export OUTPUT_ROOT=/path/to/output/unfolding_943
export NUM_PARTS=5
export M_TRANS=9.43
export N_GVEC=4
export NUM_BAND=1500
export MAX_ITER=1000
export FERMI_SHIFT=0.0
sbatch --array=1-5 scripts/slurm_band_unfolding_array.sh
```

## Input Files

- `rlat.dat`: Reciprocal lattice vectors.
- `orbital_types.dat`: Orbital types for each site.
- `site_positions.dat`: Atomic positions.
- `hamiltonians_pred.h5`: Predicted Hamiltonians.
- `overlaps.h5`: Overlap matrices.

## Output

- `data.txt`: Unfolded energy-weight triplets for the selected k points.
- `sparse_matrix.jld`: Cached sparse matrices in the input directory after the first run.

## Notes

- For Slurm batches, split k points with `--part-index` and `--num-parts`.
- `--m-trans` accepts `3.89`, `9.43`, or a custom 3x3 matrix such as `"8 -9 0; 9 17 0; 0 0 1"`.
- Tune `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and `JULIA_NUM_THREADS` for the cluster and linear algebra backend.
- The code currently supports honeycomb lattices but can be adapted for other lattice types.

## Changelog

- Version: 3.1.0
   1. Formalized the 2025-11-19 experimental unfolding script as `src/band_unfolding.jl`.
   2. Added CLI parameters for moire matrix, band count, Arpack iteration budget, point ranges, Fermi shift, and batch partitions.
   3. Added Slurm job-array support through `scripts/slurm_band_unfolding_array.sh`.
   4. Moved sparse-matrix preprocessing outside the selected k-point loop.
- Version: latest
   1. Rewrote `src/band_unfolding.jl` using Method1 (the new implementation uses Method1 as the main unfolding engine).
   2. Added `test/9.43/` as an example case demonstrating the new unfolding algorithm.

