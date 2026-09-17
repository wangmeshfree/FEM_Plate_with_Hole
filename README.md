# MAE 5036: Advanced Computational Solid Mechanics II

This repository contains finite element method examples for MAE 5036. The runnable Python programs are in `src`; the Jupyter notebooks are step-by-step materials for classroom instruction.

## Repository Structure

| Path | Purpose |
| --- | --- | --- |
| [src/FEM1D_linear_bar](src/FEM1D_linear_bar) | Runnable 1D linear-bar FEM program. See its [local guide](src/FEM1D_linear_bar/README.md). |
| [src/FEM1D_linear_bar.ipynb](src/FEM1D_linear_bar.ipynb) | Classroom notebook for the 1D bar example. |
| [src/FEM2D_linear_plate_with_hole.ipynb](src/FEM2D_linear_plate_with_hole.ipynb) | Classroom notebook for the plate-with-hole example. |

## Requirements

- Python 3
- NumPy
- Matplotlib

```sh
python -m pip install numpy matplotlib
```

Students are encouraged to install Python and learn to use it independently. You may use a different programming language for your work, provided you can implement and explain your solution.

The program creates a mesh-specific results folder under `src/FEM1D_linear_bar` containing nodal displacements plus displacement and stress plots.

## Classroom Notebooks

Open the notebooks in VS Code or Jupyter for the guided derivations, intermediate FEM calculations, plots, and convergence study used during class instruction.

## Author

Jiarui Wang

Copyright (c) 2026.
