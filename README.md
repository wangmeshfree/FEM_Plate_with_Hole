# MAE 5036: 1D Bar Finite Element Analysis

Teaching material for MAE 5036, Advanced Computational Solid Mechanics II. This example solves the axial response of a linearly elastic bar using two-node linear finite elements and Gauss integration.

## Problem

The bar has length $L=10$, Young's modulus $E=100$, cross-sectional area $A=1$, a distributed axial load $q(x)=x^2$, and an end force $P=5$.

$$
\frac{d}{dx}\left(EA\frac{du}{dx}\right) + x^2 = 0,
\qquad u(0)=0,
\qquad EAu_{,x}(L)=P.
$$

The code computes the nodal displacements, axial stress, and the support reaction. It also compares the FEM result with the analytical solution:

$$
u(x) = \frac{\left(P + L^3/3\right)x - x^4/12}{EA}.
$$

## Files

| Path | Purpose |
| --- | --- |
| [src/FEM1D_linear_bar/main.py](src/FEM1D_linear_bar/main.py) | Runs the 1D bar analysis and saves results. |
| [src/FEM1D_linear_bar/fem1d_lib.py](src/FEM1D_linear_bar/fem1d_lib.py) | Reusable FEM functions: mesh generation, element integration, assembly, boundary conditions, and plotting. |
| [src/FEM1D_linear_bar.ipynb](src/FEM1D_linear_bar.ipynb) | Step-by-step notebook version, including a mesh-convergence study. |

## Requirements

- Python 3
- NumPy
- Matplotlib

Install the Python packages with:

```sh
python -m pip install numpy matplotlib
```

## Run the Example

From the repository root, run:

```sh
.venv/bin/python src/FEM1D_linear_bar/main.py
```

To work through the calculations interactively, open [src/FEM1D_linear_bar.ipynb](src/FEM1D_linear_bar.ipynb) in VS Code or Jupyter and run its cells in order.

## Results

For `number_element_in_one_direction = 2`, the script creates:

```text
src/FEM1D_linear_bar/Results_1d_bar_2_elements/
├── nodal_displacements.txt
├── displacement.png
└── stress.png
```

The output folder name updates automatically when the number of elements changes. The displacement file contains the nodal coordinate and FEM displacement in two columns.

## Learning Objectives

- Construct a uniform 1D finite-element mesh and element connectivity.
- Evaluate linear shape functions and their derivatives.
- Assemble the global stiffness matrix and equivalent nodal load vector.
- Apply essential and traction boundary conditions.
- Solve for nodal displacements and postprocess axial stress.
- Compare FEM results with the analytical solution.

## Author

Jiarui Wang

Copyright (c) 2026.
