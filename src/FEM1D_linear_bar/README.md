# MAE 5036: 1D Bar FEM

A teaching example for MAE 5036, Advanced Computational Solid Mechanics II. The program solves the axial displacement and stress of a linearly elastic bar using two-node linear finite elements.

## Problem

The bar has length $L=10$, Young's modulus $E=100$, cross-sectional area $A=1$, distributed load $q(x)=x^2$, and end load $P=5$.

$$
\frac{d}{dx}\left(EA\frac{du}{dx}\right) + x^2 = 0,
\qquad u(0)=0,
\qquad EAu_{,x}(L)=P.
$$

## Files

| File | Purpose |
| --- | --- |
| `main.py` | Runs the FEM analysis and saves results. |
| `fem1d_lib.py` | FEM functions for meshing, integration, assembly, boundary conditions, postprocessing, and plotting. |
| `../FEM1D_linear_bar.ipynb` | Step-by-step notebook version with a convergence study. |

## Requirements

```sh
python -m pip install numpy matplotlib
```

## Run

From this folder, run:

```sh
python main.py
```

## Results

For a mesh with `number_element_in_one_direction = 2`, the program creates:

```text
Results_1d_bar_2_elements/
├── nodal_displacements.txt
├── displacement.png
└── stress.png
```

The output folder name changes automatically with the number of elements. `nodal_displacements.txt` contains the coordinate and FEM displacement at every node.

## Learning Objectives

- Build a uniform 1D finite-element mesh.
- Assemble the global stiffness matrix and force vector.
- Apply displacement and traction boundary conditions.
- Solve for nodal displacements and postprocess axial stress.
- Compare FEM results with the analytical solution.

Copyright (c) 2026 Jiarui Wang.
