# CA1 — Incompressible Elasticity (TME250)

Solo assignment. Report: `report/TME250_CA1.pdf`.

## Summary

- **Task 1**: analytical traction–displacement (T–U) relations for two limiting boundary
  cases (ϵ11 = 0 and σ11 = 0) under plane strain, as functions of Poisson's ratio ν.
- **Task 2**: CST (Constant Strain Triangle) elements under uniform traction / prescribed
  displacement loading, with ν swept from 0.3 toward 0.5. Shows the CST element's
  volumetric locking and checkerboard pressure pattern near incompressibility, benchmarked
  against the Task 1 analytical bounds, with a mesh convergence study.
- **Task 3**: Taylor-Hood mixed element (quadratic displacement / linear pressure) for the
  same two load cases, demonstrating stable, smooth pressure fields all the way to
  ν = 0.5 (fully incompressible), in contrast to Task 2.

## File Structure

```
CA1-large-deformation/
├── README.md
├── report/
│   └── TME250_CA1.pdf
├── src/
│   ├── scripts/
│   │   ├── Task3_Taylor_hood_element.m               # Taylor-Hood mixed element: mesh convergence + pressure field for nu -> 0.5
│   │   └── Task2_CST_element.m                       # CST element: mesh convergence + checkerboard pressure field for nu -> 0.5
│   └── functions/
│       ├── Be_Taylor_Hood_implementation.m           # Symbolic derivation of Taylor-Hood element formulation
│       ├── Be_Taylor_hood.m                          # Taylor-Hood element formulation
│       ├── Loadvector_interpolation.m                # Unit edge load vector for linear (CST) boundary segments
│       ├── Loadvector_interpolation_Taylor_hood.m    # Unit edge load vector for quadratic (TRIA6) boundary segments
│       └── Taylor_hood.m                             # TRIA6/P1 mixed element stiffness routine (nearly/fully incompressible elasticity)
├── results/                                          # Exported figures (PNG)
└── problem_statements/                        
```
## Dependencies / what's excluded

- `hooke.m`, `plante.m`, `plants.m` are **CALFEM** functions (open-source) — referenced
  directly, no redistribution restriction.
- `mesh.m` (custom mesh generator) is **course-provided** and is excluded from this repo
  pending instructor permission. Code here will not run standalone without it.
- Only functions explicitly shown in the report appendices are included in `src/`; any
  other course-provided routines are likewise excluded.
