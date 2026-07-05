# CA2 — Contact Mechanics (TME250)

Solo assignment. Report: `report/TME250_CA2.pdf`.

## Summary

- **Task 1**: quasi-static contact between a rigid cylinder and an elastic block, using
  Lagrange multipliers with active-set iteration to resolve contact nodes at each load
  increment. Reports force–indentation response and von Mises stress/displacement fields.
- **Task 2**: dynamic version of the same contact problem — an elastic block with initial
  velocity impacting a stationary cylinder — solved three ways: (i) implicit dynamics
  without contact (sanity check), (ii) implicit dynamics with contact, (iii) explicit
  dynamics with contact, including a deliberately unstable large-timestep case and a
  CFL-compliant stable case. Compares energy conservation, contact force history, and
  stability behavior across schemes.

## Contents

- `src/scripts/` — main driver scripts, one per scenario (quasi-static, implicit dynamic
  without/with contact, explicit dynamic)
- `src/functions/` — `gap_lin.m` (linearized gap function / contact normal computation)
- `problem_statements/` — assignment brief (add if available)

## Dependencies / what's excluded

- `hooke.m`, `plante.m`, `plants.m` are **CALFEM** functions (open-source) — referenced
  directly, no redistribution restriction.
- `mesh.m` (custom mesh generator) is **course-provided** and is excluded from this repo
  pending instructor permission. Code here will not run standalone without it.
- Only functions explicitly shown in the report appendices are included in `src/`; any
  other course-provided routines are likewise excluded.
