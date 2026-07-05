# CA2 — Contact Mechanics (TME250)

Report: `report/TME250_CA2.pdf`.

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

## File Structure

```
CA2-Contact_Mechanics/
├── README.md
├── report/
│   └── TME250_CA2.pdf
├── src/
│   ├── Scenario_1_Quassi_Static_Contact.m       # Quasi-static cylinder indentation via Lagrange-multiplier active-set contact
│   ├── Scenario_2_Implicit_without_Contact.m    # Implicit dynamics sanity check: free oscillation, no contact, energy conservation
│   ├── Scenario_2_Implicit_with_Contact.m       # Implicit dynamics + active-set contact: impact, rebound, multi-cycle contact-separation
│   ├── Scenario_2_Explicit_Dynamics.m           # Explicit dynamics with contact: unstable large-dt vs. CFL-stable small-dt cases
│   └── gap_lin.m                                # Linearized gap function and contact normal for cylinder-node contact
├── results/                                     # Exported figures (PNG)
└── problem_statements/                        
```

## Dependencies / what's excluded

- `hooke.m`, `plante.m`, `plants.m` are **CALFEM** functions (open-source) — referenced
  directly, no redistribution restriction.
- `mesh.m` (custom mesh generator) is **course-provided** and is excluded from this repo
  pending instructor permission. Code here will not run standalone without it.
- Only functions explicitly shown in the report appendices are included in `src/`; any
  other course-provided routines are likewise excluded.
