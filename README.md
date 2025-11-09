# Finite Element Analysis - Plate with Circular Hole under Tension

A comprehensive Python implementation of finite element method (FEM) analysis for studying stress distribution in a plate with a circular hole under tensile loading. This project includes mesh convergence studies and comparison with analytical solutions.


## Overview

This project implements a 2D finite element solver for the classical problem of a plate with a circular hole subjected to uniaxial tension. The solution demonstrates:

- **Stress concentration** around the circular hole
- **Convergence characteristics** of Q4 (4-node quadrilateral) elements
- **Energy norm error** calculation by comparing with analytical solutions
- **Mesh refinement studies** to validate numerical accuracy

## Features

- ✅ **Q4 Element**: Four-node isoparametric quadrilateral elements
- ✅ **Plane Strain Formulation**: Material stiffness matrix for plane strain problems
- ✅ **Gauss Integration**: 2×2 Gauss quadrature for element stiffness matrices
- ✅ **Energy Norm Error**: Rigorous error estimation compared to exact solution
- ✅ **Convergence Analysis**: Mesh refinement studies with convergence rate calculation
- ✅ **Visualization**: Mesh plots and stress contour visualizations

## 📐 Theory

### Problem Setup

A square plate with a circular hole at the center is subjected to uniform tensile stress. Due to symmetry, only one quarter of the plate is modeled:

- **Elastic Modulus (E)**: 1000
- **Poisson's Ratio (ν)**: 0.30
- **Hole Radius (R)**: 1.0
- **Plate Size (L)**: 5.0
- **Applied Stress (Tx)**: 1.0

### Analytical Solution

The exact stress field in polar coordinates is given by:

$$\sigma_{rr} = \frac{T_x}{2}\left(1 - \frac{R^2}{r^2}\right) + \frac{T_x}{2}\left(1 - 4\frac{R^2}{r^2} + 3\frac{R^4}{r^4}\right)\cos(2\theta)$$

$$\sigma_{\theta\theta} = \frac{T_x}{2}\left(1 + \frac{R^2}{r^2}\right) - \frac{T_x}{2}\left(1 + 3\frac{R^4}{r^4}\right)\cos(2\theta)$$

$$\sigma_{r\theta} = -\frac{T_x}{2}\left(1 + 2\frac{R^2}{r^2} - 3\frac{R^4}{r^4}\right)\sin(2\theta)$$

At the hole boundary (r = R, θ = π/2), the maximum stress concentration factor is **3.0**.


### Customization

You can modify key parameters in the code:

```python
# Material properties
E = 1.0e3      # Elastic modulus
nu = 0.30      # Poisson's ratio

# Geometry
R = 1.0        # Hole radius
L = 5.0        # Plate size

# Mesh refinement levels
refinements = [10, 15, 20]  # Number of elements in one direction
```

## 📊 Results

### Mesh Convergence

The implementation demonstrates **first-order convergence on the energy error norm** (O(h)) for Q4 elements:

| Mesh Size (h) | Energy Norm Error | Convergence Rate |
|---------------|-------------------|------------------|
| 0.4000        | ~0.XXX           | -                |
| 0.2667        | ~0.XXX           | ~1.0             |
| 0.2000        | ~0.XXX           | ~1.0             |

### Stress Distribution

The analysis produces three stress component contours:
- ** $$\sigma_{11}$$**: Normal stress in x-direction
- ** $$\sigma_{22}$$**: Normal stress in y-direction  
- ** $$\sigma_{12}$$**: Shear stress

Maximum stress occurs at the hole boundary, confirming the theoretical stress concentration factor of 3.0.



## 📚 Mathematical Background

### Finite Element Formulation

The weak form of the equilibrium equation:

$$\int_\Omega \delta\varepsilon^T \sigma \ d\Omega = \int_{\Omega} \delta u^T b \ d\Omega + \int_{\Gamma_t} \delta u^T t \ d\Gamma$$

Where:
- **ε**: Strain tensor
- **σ**: Stress tensor
- **t**: Traction vector
- **Ω**: Domain
- **Γₜ**: Traction boundary

### Element Stiffness Matrix

$$K_e = \int_{\Omega_e} B^T D B \, d\Omega$$

Where:
- **B**: Strain-displacement matrix
- **D**: Material stiffness matrix (plane strain)
- **Ωₑ**: Element domain

### Energy Norm Error

$$\|e\|_E = \sqrt{\int_\Omega (\varepsilon_{FEM} - \varepsilon_{exact})^T D (\varepsilon_{FEM} - \varepsilon_{exact}) \, d\Omega}$$


## 📖 References

1. **Timoshenko, S. P., & Goodier, J. N.** (1970). "Theory of Elasticity"
2. **Hughes, T. J. R.** (2000). "The Finite Element Method: Linear Static and Dynamic Finite Element Analysis"
3. **Zienkiewicz, O. C., & Taylor, R. L.** (2000). "The Finite Element Method"


## 👤 Author

Jiarui Wang & Huanyang Hou @ Southern University of Science and Technology- [@wangmeshfree](https://github.com/wangmeshfree)


**Note**: This is an educational implementation. For production finite element analysis, consider using established libraries like FEniCS, deal.II, or commercial software.
