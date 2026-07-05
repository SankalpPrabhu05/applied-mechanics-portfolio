# Applied Mechanics Portfolio

FEM, CAE, and applied mechanics coursework from my MSc in Applied Mechanics at Chalmers University of Technology.

Each course is organized as its own top-level folder, with computer assignments (CAs) broken out further where a course has more than one. Every assignment folder includes the full report, source code (where sharable), and — where applicable — a note on any course-provided/third-party dependencies that are excluded from this repo.

## Courses

| Course | Topic | Tools |
|---|---|---|
| [`Solid_Mechanics_TME235/`](./Solid_Mechanics_TME235) | Beam theory, axisymmetric elasticity, plate bending, hyperelasticity — solved analytically, numerically, and validated against Abaqus | Python (SymPy, NumPy, Pandas, Matplotlib), Abaqus, CALFEM |
| [`FEM_Structures_TME245/`](./FEM_Structures_TME245) | Nonlinear large-deformation (Yeoh hyperelastic) analysis and Kirchhoff plate bending/buckling | Abaqus, MATLAB (custom TRIA6 element, Newton-Raphson, generalized eigenvalue buckling) |
| [`FEM_Solid_TME250/`](./FEM_Solid_TME250) | Incompressible elasticity (CST vs. Taylor-Hood mixed elements), Lagrange-multiplier contact mechanics (quasi-static + implicit/explicit dynamics), and adaptive mesh refinement (a-posteriori error estimation) | MATLAB, CALFEM |

## Repository Structure

Each course folder generally follows:

```
<Course_Name>/
├── README.md
├── CA<n>-<topic>/              (or a flat structure, if the course has one assignment)
│   ├── README.md
│   ├── report/                 full write-up (PDF)
│   ├── src/
│   │   ├── scripts/            runnable analysis scripts
│   │   └── functions/          element routines, helper functions
│   ├── data/                   (where applicable)
│   ├── results/                key output plots/figures
│   └── problem_statements/     original assignment brief
```

Structure is adapted per course/assignment based on how similar the individual CAs are to each other — see each course's own README for its exact layout and any deviations.

## Dependencies and Code Availability

- **CALFEM**: several scripts call functions from [CALFEM](https://github.com/CALFEM/calfem-matlab), an open-source MATLAB toolbox for FE analysis developed at Lund University. These calls are included as-is; CALFEM itself is not bundled here and should be installed separately to run the code.
- **Course-provided routines**: some functions (e.g. mesh generators, load-application helpers) were provided by course instructors and are **not included** in this repository, pending instructor permission to publish. As a result, most scripts here will not run standalone, even with CALFEM installed — each assignment's README lists the specific excluded functions.
- Only functions/scripts that appear explicitly in the report appendices are included, since these represent work authored (or explicitly reproduced with attribution) for the assignment.
- Some assignments were completed with a co-author (noted per course/CA); solo work is noted where applicable.

## License

This repository is released under the MIT License, unless otherwise noted (e.g. for third-party CALFEM code, which retains its own license).

## Author

Sankalp Manojkumar Prabhu — MSc Applied Mechanics, Chalmers University of Technology
[sankalpm@chalmers.se](mailto:sankalpm@chalmers.se)
