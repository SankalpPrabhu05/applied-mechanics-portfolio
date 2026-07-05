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

## Contents

- `src/scripts/` — main driver scripts for Task 2 (CST) and Task 3 (Taylor-Hood)
- `src/functions/` — `Loadvector_interpolation.m` (linear edge load interpolation),
  `Taylor_hood.m` (Taylor-Hood element stiffness routine), and its symbolic-derivation
  supporting script
- `problem_statements/` — assignment brief (add if available)

## Dependencies / what's excluded

- `hooke.m`, `plante.m`, `plants.m` are **CALFEM** functions (open-source) — referenced
  directly, no redistribution restriction.
- `mesh.m` (custom mesh generator) is **course-provided** and is excluded from this repo
  pending instructor permission. Code here will not run standalone without it.
- Only functions explicitly shown in the report appendices are included in `src/`; any
  other course-provided routines are likewise excluded.
