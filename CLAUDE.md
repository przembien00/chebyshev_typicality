# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This codebase implements the **Chebyshev Expansion Algorithm for Typicality (CET)** to compute dynamical spin correlations in Heisenberg spin systems at thermal equilibrium. It uses quantum typicality to avoid full diagonalization of large Hilbert spaces.

The project is a hybrid C++/Python research codebase:
- **C++ core** (`Algorithm/`): high-performance MPI-parallelized simulation
- **Python scripts**: data analysis and visualization of HDF5 output

## Build

```bash
# From repo root
cmake -S . -B build
cmake --build build
```

Requires: MPI, HDF5, Boost, Blaze (linear algebra), lambda-lanczos.

Produces: `executable_DOUBLE.out` and `executable_rc_DOUBLE.out`.

## Running Simulations

```bash
mpirun -n <cores> ./executable_DOUBLE.out \
  --spinmodel=ISO \
  --srcfile=Couplings/<coupling_file>.h5 \
  --beta=<inverse_temperature> \
  --num_TimePoints=<N> \
  --rescale=<factor>
```

HPC batch submission uses `array_job.sh` or `numvec_array_job.sh`.

## C++ Architecture

```
Algorithm/
├── main.cpp                  # MPI parallel main loop, typicality sampling
├── Parameter_Space/          # CLI argument parsing and parameter validation
├── Hamiltonians/             # Quantum Hamiltonian operators
├── Functions/                # Core physics: Chebyshev evolution, spin ops, correlations
├── Types/                    # State, correlation, and tensor data structures
└── Storage_Concept/          # HDF5 I/O for saving correlations and metadata
```

The main loop in `main.cpp`:
1. Initializes random quantum states (typicality samples)
2. Applies Chebyshev expansion for thermal state preparation (`e^{-βH/2}`)
3. Evolves spin operators via Chebyshev time evolution
4. Computes spin-spin correlations (xx, xy, zz, etc.) including statistical errors
5. Writes results to HDF5

## Data Format

Output lives in `Data/` as HDF5 files. Key stored quantities:
- `num_TimePoints`, `beta`, system size `N`, rescale factors
- Real and imaginary spin correlations per component
- Statistical errors from typicality sampling

Coupling matrices (lattice definitions) are HDF5 files in `Couplings/`. Generate them with `make_couplings*.py` scripts.

## Python Visualization

Each `Plot_*.py` script reads one or more HDF5 files from `Data/` and produces matplotlib figures. Scripts use a shared `.mplconfig/` for consistent styling.

Key scripts:
- `Plot_sdmft.py` — compare CET vs spin-DMFT results
- `Plot_h_field.py` / `Plot_h.py` — magnetic field parameter studies
- `Plot_CspinDMFT_sites.py` — site-resolved spin-DMFT correlations

## Supported Physics

- Spin models: isotropic Heisenberg (ISO) with customizable coupling matrices
- Lattice types: square (nearest-neighbor), chain, frustrated, disordered (random couplings)
- Observables: imaginary-time spin correlations S(τ) for all components
- Also supports full exact diagonalization mode for small systems (validation)
