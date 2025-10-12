# Project Plan: Separate Two Research Projects

## Problem Analysis
The repository currently contains two distinct research projects that are intermingled:

1. **Heuristic Default Model (Julia Project)**:
   - Location: Root directory, `/full_info/`, `/baseline_info/` folders
   - Main files: `full_info.jl`, `lambda_info.jl`, `lambda_info_calibration.jl`
   - Results: `Result/` folder with visualization plots (Figs/)
   - Language: Julia with plotting capabilities

2. **PRO Default Model (Fortran + MATLAB Project)**:
   - Location: `/long_term/defaultModel/` folder
   - Source files: Fortran `.f90` files in `src/` folder
   - Processing: MATLAB `.m` files for results analysis
   - Compiled binaries and simulation results

## Proposed Solution
Create two separate project directories to clearly delineate the research:

### Directory Structure After Separation:

```
/
├── README.md (updated with project descriptions)
├── heuristic-default-model/     # Julia project
│   ├── src/
│   │   ├── full_info.jl
│   │   ├── lambda_info.jl
│   │   └── lambda_info_calibration.jl
│   ├── results/
│   │   └── [current Result/ contents]
│   ├── scripts/
│   │   └── [current scripts/ contents if Julia-related]
│   └── README.md
├── pro-default-model/     # Fortran + MATLAB project
│   ├── src/
│   │   └── [current long_term/defaultModel/src/ contents]
│   ├── matlab/
│   │   └── [MATLAB processing files]
│   ├── results/
│   │   └── [current long_term/defaultModel/results/ contents]
│   ├── Makefile
│   ├── compile scripts
│   └── README.md
└── paper/                       # Shared paper materials
    └── [current paper/ contents]
```

## Todo Items
- [ ] Create `heuristic-default-model/` directory structure
- [ ] Move Julia files to heuristic project directory
- [ ] Create `pro-default-model/` directory structure  
- [ ] Move Fortran/MATLAB files to PRO project directory
- [ ] Update file paths in source files to reflect new structure
- [ ] Create individual README files for each project
- [ ] Update main README.md to describe both projects
- [ ] Test that both projects still compile/run correctly
- [ ] Clean up any remaining orphaned files

## Key Considerations
- Preserve all existing functionality - no code changes
- Maintain relative paths within each project
- Keep paper materials accessible to both projects
- Follow principle of simplicity - clear separation without complexity
- Ensure both projects remain independently runnable

## Review

### Changes Made
Successfully separated the two research projects into distinct directories:

1. **Created heuristic-default-model/ directory**: Contains Julia-based heuristic default model
   - Moved `full_info/`, `baseline_info/` to `src/`
   - Moved `Result/` to `results/`
   - Moved `scripts/` directory
   - Created project-specific README.md

2. **Created pro-default-model/ directory**: Contains Fortran + MATLAB PRO model
   - Copied Fortran source files from `long_term/defaultModel/src/` to `src/`
   - Copied MATLAB files to `matlab/` directory
   - Copied simulation results to `results/`
   - Copied build files (Makefile, compilation scripts)
   - Created project-specific README.md

3. **Updated main README.md**: Now clearly describes both projects with separate sections and links to individual project documentation

4. **Cleaned up repository**: Removed original directories (`long_term/`, `results/`) and orphaned files (gumbel plots) after successful separation

### Final Structure
```
/
├── README.md (updated)
├── heuristic-default-model/     # Julia project
│   ├── src/
│   ├── results/
│   ├── scripts/
│   └── README.md
├── pro-default-model/     # Fortran + MATLAB project
│   ├── src/
│   ├── matlab/
│   ├── results/
│   ├── Makefile
│   └── README.md
└── paper/                       # Shared paper materials
```

Both projects are now clearly separated while maintaining their functionality and preserving all original code and results.
Successfully merged two extension subsections into one cohesive subsection with three paragraphs:

1. **Removed duplicate subsection headers**: Eliminated the first "Theoretical Extension" subsection header and integrated content under the unified "Theoretical Extensions" header.

2. **Created coherent three-paragraph structure**:
   - **Paragraph 1: Ramsey Planning Under PRO** - Welfare analysis showing persistent loss even under optimal fiscal policy
   - **Paragraph 2: Endogenous Belief Formation** - Bayesian learning with negativity bias leading to persistent PRO beliefs
   - **Paragraph 3: Optimal Policy Communication** - Strategic transparency choice to mitigate PRO-constrained beliefs

3. **Added smooth transitions**: Included connecting sentences between paragraphs to create logical flow from welfare analysis → belief formation → communication strategy.

4. **Preserved mathematical rigor**: All propositions, equations, and appendix references remain completely intact.

### Final Structure
The merged subsection now provides a comprehensive view of three theoretical extensions while maintaining academic rigor and logical coherence. The structure follows the principle of simplicity with minimal changes to achieve the integration goal.

### Compilation Fix
Fixed LaTeX compilation error on line 2383 by properly formatting the transition operator equation:
- Changed `T_\mathcal{H}` to `T_{\mathcal{H}}` to ensure proper subscript grouping
- Document now compiles successfully generating a 65-page PDF
- All mathematical notation and references remain intact
