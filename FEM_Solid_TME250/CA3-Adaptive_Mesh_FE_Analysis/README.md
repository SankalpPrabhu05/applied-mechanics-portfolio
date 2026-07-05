# CA2 — Contact Mechanics (TME250)

Solo assignment. Report: `report/TME250_CA2.pdf`.

## Summary

- **Task 1 — Analytical (manufactured) solution:** A displacement field satisfying
  homogeneous Dirichlet boundary conditions on a 40×20 mm rectangular domain is chosen
  symbolically. Strains, stresses, and the corresponding body force field are derived,
  and the exact strain energy norm (‖u‖ₐ = 19.117) is computed for use as a reference
  in later tasks.
- **Task 2 — Convergence study:** Relative error in energy norm vs. number of degrees
  of freedom is tracked over 11 uniform refinement levels for both linear (3-node CST)
  and quadratic (6-node) triangular elements, confirming the expected convergence rates
  (q ≈ 1 for linear, q ≈ 2 for quadratic).
- **Task 3 — Error estimator / effectivity index:** For each mesh, the linear solution
  is compared against a quadratic solution on the same topology (linear solution
  projected via `lin2quad`) to estimate the local/global error. The effectivity index
  η = (estimated error)/(exact error) is tracked with refinement and shown to approach 1.
- **Task 4 — Adaptive mesh refinement:** An iterative, threshold-based adaptive
  refinement algorithm (refine elements whose estimated error exceeds 50% of the
  maximum element error) is implemented and compared against uniform refinement,
  showing the adaptive strategy reaches a given error tolerance with fewer degrees
  of freedom.
## Contents

- `src/scripts/` — main driver scripts, one per scenario (quasi-static, implicit dynamic
  without/with contact, explicit dynamic)
- `src/functions/` — `gap_lin.m` (linearized gap function / contact normal computation)
- `problem_statements/` — assignment brief (add if available)

## Dependencies / what's excluded

- `hooke.m`, `plante.m`, `plant6e.m` are **CALFEM** functions (open-source) — referenced
  directly, no redistribution restriction.
- `grid2CALFEM.m`, `refinegrid.m`, `lin2quad.m`, `eldraw2.m` (custom mesh generator) is **course-provided** and is excluded from this repo
  pending instructor permission. Code here will not run standalone without it.
- Only functions explicitly shown in the report appendices are included in `src/`; any
  other course-provided routines are likewise excluded.
