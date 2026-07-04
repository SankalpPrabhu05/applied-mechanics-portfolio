# Solid Mechanics — Computer Assignments (TME235, Chalmers)

**Course:** Mechanics of Solids, TME235 | Chalmers University of Technology
**Tools:** Python (SymPy, NumPy, Pandas, Matplotlib), Abaqus, CALFEM
**Collaborator:** Tapasi Himanth Karthik

## Summary
Four computer assignments covering beam theory, axisymmetric elasticity, plate bending, and hyperelasticity — each solved analytically (symbolic derivation in Python), numerically (custom FEM implementation), and validated against Abaqus.

## Assignments

**CA1 — Cantilever Beam Under Point Load**
Deflection and stress in a clamped beam carrying a point load at its free end, solved via Euler-Bernoulli and Timoshenko beam theory, and compared against Abaqus and a custom CALFEM-based FE model for two beam lengths (2 m and 0.25 m). Includes a mesh convergence study and analysis of stress concentration at the clamped boundary, with design recommendations (fillets) to reduce it.

**CA2 — Axisymmetric Disc Under Pressure**
Analytical solution (via symbolic integration of the governing ODE) and a from-scratch FEM implementation (element stiffness assembly in SymPy/NumPy) for radial displacement and stress in a pressurized disc, extended to a variable-thickness case.

**CA3 — Axisymmetric Plate Bending**
Kirchhoff plate bending theory solved analytically (MATLAB/Python) and via Abaqus, for a plate free at the inner radius and simply supported at the outer radius. Investigates the effect of halving plate thickness and mesh refinement on maximum deflection, von Mises stress, and stress concentration.

**CA4 — Beam on Hyperelastic (Rubber) Support**
Steel beam supported by a compressible neo-Hookean rubber block, modeled in Abaqus under both tension (P = +10 kN) and compression (P = -10 kN) loading. Validates the small-strain (beam) / large-strain (rubber) modeling assumption via principal strain contours, and compares deflection against an Euler-Bernoulli beam model without rubber support.

## Key Results
- CA1: FEM (Abaqus/CALFEM) and analytical beam models agreed closely for the long beam (L = 2 m); Timoshenko theory matched FEM more closely than Euler-Bernoulli for the short, thick beam (L = 0.25 m) due to shear effects
- CA2: radial displacement and stress from the custom FEM code matched the analytical solution closely (0.044 mm displacement, 135.6 MPa stress at the inner boundary) for constant thickness
- CA3: von Mises stress increased from 5.2 GPa to 5.7 GPa when mesh was refined from 1 mm to 0.1 mm elements, illustrating stress singularity/ mesh-dependence at a geometric discontinuity
- CA4: the rubber support reduced beam-tip deflection by ~73% compared to an unsupported cantilever (64 mm vs. 235 mm under a 10 kN load)

## Repository Structure
'''
solid-mechanics-tme235/
├── report/               Full write-ups with methodology, plots, and discussion
├── src/                  Python scripts for each assignment (CA1-CA4)
├── data/                 Abaqus-exported CSV/TXT data used for comparison plots
├── results/              Key output plots and Abaqus contour screenshots
└── problem_statements/   Original assignment briefs
'''
## Notes
Co-authored with Tapasi Himanth Karthik as part of coursework requirements for TME235, Chalmers University of Technology (Autumn 2024).
