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
│   └── TME245_CA2.pdf
├── src/
│   ├── CA2_Task_1_Kirchoff_Plate.m      # Task 1B-1D: element routine verification,
│   │                                    #   displacement solve, in-plane stress + Von
│   │                                    #   Mises/FoS computation
│   ├── CA2_Task_2_Buckling_Problem.m    # Task 2A-2B: geometric stiffness assembly +
│   │                                    #   generalized eigenvalue buckling solve
│   ├── Full_KQuad_Func.m                # Kirchhoff plate element routine (membrane +
│   │                                    #   bending stiffness, external forces)
│   ├── K_inplane_S.m                    # In-plane stress computation (Gauss point
│   │                                    #   averaged, Voigt form)
│   ├── Kirch_G.m                        # Geometric stiffness matrix (3x3 Gauss scheme)
│   ├── Be_Quad_Symb2Func.m              # Symbolic derivation -> membrane B0-matrix
│   │                                    #   function (Be_Quad_func)
│   └── Bstar_KQuad_Symb2Func.m          # Symbolic derivation -> Kirchhoff bending
│                                        #   B*-matrix function (Bast_kirchoff_func,
│                                        #   Nk_func, Be_kirch_func)
├── results/                             # Exported figures (PNG)
└── problem_statements/               
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
