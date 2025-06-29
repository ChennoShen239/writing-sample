# Pessimism Default Model (Fortran + MATLAB)

This project implements a pessimism-based sovereign default model using Fortran for computation and MATLAB for analysis.

## Structure

- `src/` - Fortran source files (.f90)
- `matlab/` - MATLAB processing and visualization scripts
- `results/` - Model simulation results and binary files
- `Makefile` - Build configuration
- `*.bat`, `*.sh` - Compilation scripts

## Building the Model

Use the provided Makefile or compilation scripts:
```bash
make
# or
./compile.sh
```

## Running Simulations

After compilation, run the model executable and process results with MATLAB scripts in the `matlab/` directory.

## Results

Binary result files (.bin) and processed data (.mat) are stored in the `results/` directory.