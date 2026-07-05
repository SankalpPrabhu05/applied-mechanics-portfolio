# CA1 — Large Deformation Analysis

Part of [TME245 - Finite Element Method: Structures](../README.md), Chalmers University
of Technology, SP3 2024/2025.

**Authors:** Sankalp Manojkumar Prabhu, Paulo Liborio

## Summary

This assignment investigates the large-deformation behavior of a rubber square profile
(with three rectangular slots) under uniaxial tension and compression, using the nearly
incompressible **Yeoh hyperelastic material model**.

**Task 1 — Abaqus:** The rubber profile was modeled in Abaqus using CPE6M (6-node
modified quadratic plane strain triangle) elements, with NLGEOM enabled for large
deformation. Two load cases were analyzed:
- **LC1 (Tension):** +20 mm prescribed displacement
- **LC2 (Compression):** −15 mm prescribed displacement

Von Mises stress, longitudinal strain, and reaction force were evaluated, along with a
mesh convergence study for LC1. Results showed higher peak Von Mises stress in tension
(20.37 MPa) than compression (15.03 MPa), while reaction force was higher in compression
— consistent with polymer chain stretching/stiffening under tension versus lateral
expansion (incompressibility) under compression.

**Task 2 — Matlab:**
- **2A:** Linear elastic verification of the same geometry/BCs using a Newton-Raphson
  solver (sanity check before introducing hyperelasticity).
- **2B:** Comparison of Yeoh vs. Neo-Hookean stress-stretch response, showing the two
  models agree closely for small stretches (0.9 < λ < 1.1) but diverge at larger
  deformations, since Neo-Hooke is a first-order approximation and cannot capture
  strain-stiffening.
- **2C:** Verification of a custom 6-node triangular (TRIA6) large-deformation element
  routine against a reference deformation gradient/stress test case.
- **2D:** Full nonlinear rubber profile analysis using the TRIA6 element with the Yeoh
  model, compared against the Abaqus results (Task 1) and evaluated for mesh dependency.

## Key Results

| Load Case | Mesh | Reaction Force [N] | Max Von Mises [MPa] |
|---|---|---|---|
| Compressive | Coarse | 19189.03 | 11.42 |
| Compressive | Fine | 19175.06 | 12.10 |
| Tensile | Coarse | -16167.14 | 17.62 |
| Tensile | Fine | -16165.55 | 17.57 |

Full discussion, including mesh convergence and comparison against Abaqus, is in the
report (`report/TME245_CA1_Group27.pdf`).

## Folder Structure

```
CA1-large-deformation/
├── README.md
├── report/
│   └── TME245_CA1.pdf
├── src/
│   ├── scripts/
│   │   ├── CA1_T2A.m                           # Newton-Raphson linear elastic verification
│   │   ├── CA1_T2B.m                           # Stress-stretch model comparison
│   │   ├── CA1_T2C.m                           # TRIA6 element routine verification
│   │   └── CA1_T2D_Coarse_Compression.m        # Full nonlinear rubber profile analysis
│   └── functions/
│       ├── NH_func.m                           # Neo-Hookean stress/tangent (symbolic-derived)
│       ├── Yeoh_func_S33.m                     # Yeoh stress/tangent (symbolic-derived)
│       ├── Be0_TRIA6_func.m                    # Symbolic derivation of TRIA6 B0-matrix
│       └── TRIA6_LDef.m                        # TRIA6 large-deformation element routine
├── results/                                    # Exported figures (PNG)
└── problem_statements/                        
```

## Dependencies / What Won't Run Standalone

- **CALFEM** functions used in these scripts (e.g. `eldraw2`, `eldisp2`, `plante`,
  `plantf`, `plants`, `hooke`, `extract_ed`) are part of the open-source
  [CALFEM](https://github.com/CALFEM/calfem-matlab) Matlab toolbox and are not included
  here — install CALFEM separately if attempting to run this code.
- **Course-provided files** (e.g. mesh/topology generation such as
  `topology_coarse_3node.mat`, `topology_medium_6node.mat`, and any associated mesh-import
  routines) were supplied by the course instructor and are excluded. Without these, the scripts will not run standalone.
- Only the functions and scripts explicitly listed in the report's appendices (A and B)
  are included in `src/`.
