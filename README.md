# Band Unfolding for Moiré Systems
Author: Hsu Quinn
Date: April 2025

This project implements band unfolding for honeycomb lattices, specifically designed for Moiré superlattices. The code is written in Julia and supports sparse matrix computations and band structure visualization.

## Features

- Band unfolding for honeycomb lattices.
- Efficient sparse matrix handling.
- Visualization of energy bands along the K->Γ->M path.

## Requirements

Install the following Julia packages:
```julia
using Pkg
Pkg.add(["DelimitedFiles", "LinearAlgebra", "JSON", "HDF5", "ArgParse", "SparseArrays", "Pardiso", "Arpack", "LinearMaps", "JLD", "PyPlot", "Statistics"])
```

## Usage

1. **Prepare Input Data**: Place required files (`rlat.dat`, `orbital_types.dat`, `site_positions.dat`, etc.) in the input directory.
2. **Run the Script**:
   ```bash
   julia src/band_unfolding.jl --input_dir ./test/deephdata --output_dir ./output
   ```
3. **View Results**: The output directory will contain `BandUnFolDing.png`, showing the unfolded band structure.

## Input Files

- `rlat.dat`: Reciprocal lattice vectors.
- `orbital_types.dat`: Orbital types for each site.
- `site_positions.dat`: Atomic positions.
- `hamiltonians_pred.h5`: Predicted Hamiltonians.
- `overlaps.h5`: Overlap matrices.

## Output

- `BandUnFolDing.png`: Band structure visualization.

## Notes

- Ensure `OMP_NUM_THREADS` and `JULIA_NUM_THREADS` are set for optimal performance.
- The code currently supports honeycomb lattices but can be adapted for other lattice types.