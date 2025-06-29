# Heuristic Default Model (Julia)

This project implements a heuristic default model using Julia for economic analysis.

## Structure

- `src/` - Julia source files
  - `full_info/` - Full information model implementation
  - `baseline_info/` - Baseline model with lambda information
- `results/` - Model results and visualizations
- `scripts/` - Utility scripts for plotting and analysis

## Running the Model

The main model files are:
- `src/full_info/full_info.jl` - Full information model
- `src/baseline_info/lambda_info.jl` - Lambda information model
- `src/baseline_info/lambda_info_calibration.jl` - Model calibration

## Results

Generated figures and results are stored in the `results/` directory, organized by model type (full_info, lambda_info, etc.).