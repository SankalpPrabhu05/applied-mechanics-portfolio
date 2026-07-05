# FEM Structures (TME245) — Chalmers University of Technology

Academic year 2024/2025, Study Period 3
Course: TME245 - Finite Element Method: Structures

**Authors (Group 27):**
- Sankalp Manojkumar Prabhu ([sankalpm@chalmers.se](mailto:sankalpm@chalmers.se))
- Paulo Liborio ([librio@chalmers.se](mailto:librio@chalmers.se))

---

## Overview

This repository contains two computer assignments (CAs) completed for TME245, covering
nonlinear large-deformation analysis and linear/buckling analysis of plate structures,
each combining commercial FE software (Abaqus) with custom Matlab implementations of the
underlying element formulations.

| Assignment | Topic | Tools |
|---|---|---|
| [CA1 — Large Deformation Analysis](./CA1-large-deformation) | Nonlinear hyperelastic (Yeoh) analysis of a rubber profile under uniaxial tension/compression | Abaqus, Matlab (6-node TRIA element, Newton-Raphson) |
| [CA2 — Kirchhoff Plate Analysis](./CA2-kirchhoff-plate) | Linear elastic plate bending/membrane analysis and linearized pre-buckling of a pool wall section | Matlab (Kirchhoff plate element, generalized eigenvalue buckling) |

Each assignment folder contains its own README with a more detailed summary, the full
report PDF, source code, and (where available) result figures.

## Dependencies and Code Availability

- **CALFEM**: Several scripts call functions from [CALFEM](https://github.com/CALFEM/calfem-matlab),
  an open-source Matlab toolbox for FE analysis developed at Lund University. These calls
  are included as-is in the scripts; CALFEM itself is not bundled in this repository and
  should be installed separately if you want to attempt to run the code.
- **Course-provided routines**: Some functions used in these assignments (e.g. mesh
  generation, load application helpers) were provided by the course instructor and are
  **not included** in this repository, pending permission from the instructor to publish
  them. As a result, none of the scripts here will run standalone, even with CALFEM
  installed — see each CA's README for the specific list of excluded functions.
- Only functions and scripts that appear explicitly in the report appendices are included
  here, since these represent work authored or reproduced by the group for the assignment.

## Repository Structure

```
FEM_Structures_TME245/
├── README.md                       (this file)
├── CA1-large-deformation/
│   ├── README.md
│   ├── report/
│   ├── src/
│   │   ├── scripts/                (runnable analysis scripts)
│   │   └── functions/               (element routines, material models)
│   ├── data/
│   ├── results/
│   └── problem_statements/
└── CA2-kirchhoff-plate/
    ├── README.md
    ├── report/
    ├── src/
    │   ├── scripts/
    │   └── functions/
    ├── data/
    ├── results/
    └── problem_statements/
```

## License

This repository is released under the MIT License, unless otherwise noted (e.g. for
third-party CALFEM code, which retains its own license).
