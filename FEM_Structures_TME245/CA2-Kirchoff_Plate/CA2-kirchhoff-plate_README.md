# CA2 — Kirchhoff Plate Analysis

Part of [TME245 - Finite Element Method: Structures](../README.md), Chalmers University
of Technology, SP3 2024/2025.

**Authors:** Sankalp Manojkumar Prabhu, Paulo Liborio

## Summary

This assignment analyzes a pool wall section (2 m × 2 m, 10 mm thick, clamped on three
edges with the top edge free) using a custom **Kirchhoff plate element** implemented in
Matlab, subjected to two combined loads:
- A hydrostatic pressure from water inside the pool, q(y) = ρgȳ (increasing linearly
  with depth)
- An in-plane traction on part of the top edge, representing a person of mass 150 kg
  standing on a platform

**Task 1 — Linear Elastic Plate Analysis:**
- Derivation of the weak form (decoupled membrane and bending problems) and boundary
  condition treatment.
- Implementation and verification of the Kirchhoff plate element routine
  (`Full_KQuad_Func`), computing membrane stiffness, bending stiffness, and both external
  force contributions.
- Full FE solution for in-plane and out-of-plane displacements (max out-of-plane
  displacement ≈ 13.1 mm at the lower-middle section).
- In-plane stress computation at Gauss points, Von Mises equivalent stress, and Factor
  of Safety (FoS) evaluated at the top, middle, and bottom surfaces through the
  thickness. Minimum FoS ≈ 3.47 at the outer surfaces (bending-dominated), confirming
  the structure meets the FS ≥ 2 yield strength requirement.

**Task 2 — Linearized Buckling Analysis:**
- Computation of the membrane stress state from the in-plane solution, used to build the
  geometric stiffness matrix (`Kirch_G`) via a 3×3 Gauss integration scheme (required
  because the Kirchhoff B-matrix is inherently second-order).
- Solution of the generalized eigenvalue problem `Kww_FF z = -λ Gr_FF z` to find the
  critical buckling load multiplier. Since the geometric stiffness was built from the
  actual applied loads (not a unit reference load), the resulting eigenvalue directly
  represents "how many times the current load could increase before buckling" — found to
  be λ₁ ≈ 132.86, i.e. an effective buckling factor of safety of ≈ 133.
- Visualization of the first buckling eigenmode (mode shape only — the eigenvector is
  unscaled/arbitrary in magnitude, showing the pattern of instability rather than a
  physical displacement).

## Key Results

| z position | Max Von Mises [MPa] | Min FoS |
|---|---|---|
| −h/2 (bottom) | 129.67 | 3.47 |
| 0 (mid-plane) | 0.13 | 3396 |
| h/2 (top) | 129.74 | 3.47 |

**Buckling load multiplier:** λ₁ ≈ 132.86 (safe against instability under linear
elastic/geometric assumptions, no imperfections).

Full discussion, derivations, and figures are in the report
(`report/TME245_CA2_Group27.pdf`).

## Folder Structure

```
CA2-kirchhoff-plate/
├── README.md
├── report/
│   └── TME245_CA2_Group27.pdf
├── src/
│   ├── scripts/
│   │   ├── task1c_solve_displacements.m    # Full FE solve: in-plane + out-of-plane
│   │   ├── task1d_stress_function.m        # In-plane stress + Von Mises/FoS (K_inplane_S)
│   │   ├── task2a_geometric_stiffness.m    # Geometric stiffness matrix (Kirch_G)
│   │   └── task2b_buckling_problem.m       # Full linearized pre-buckling analysis
│   └── functions/
│       ├── kirchhoff_element_routine.m     # Full_KQuad_Func — membrane + bending element
│       ├── kirchhoff_bmatrix_derivation.m  # Symbolic derivation: Nk_func, Be_kirch_func,
│       │                                    #   Bast_kirchoff_func
│       └── quad_bmatrix_derivation.m       # Symbolic derivation: Be_Quad_func (membrane B0)
├── data/                                   # Mesh/topology files, if available
├── results/                                # Exported figures (PNG)
└── problem_statements/                     # Not available for this assignment
```

## Dependencies / What Won't Run Standalone

- **CALFEM** functions (e.g. `hooke`, `extract_ed`) used in these scripts are part of the
  open-source [CALFEM](https://github.com/CALFEM/calfem-matlab) Matlab toolbox and are not
  included here — install CALFEM separately if attempting to run this code.
- **Course-provided files**, including the mesh generation routine (`rectMesh`) and the
  pressure-load helper (`pressure_load`), were supplied by the course instructor and are
  excluded pending permission to publish. Without these, the scripts will not run
  standalone.
- Only the functions and scripts explicitly listed in the report's appendices (A and B)
  are included in `src/`.
