# ******************************************************************************
# *                                                                            *
# *        MAE 5036 - Advanced Computational Solid Mechanics II                *
# *                                                                            *
# *                  1D Bar Finite Element Analysis                            *
# *                             JIARUI WANG                                    *
# *                         Copyright (c) 2026                                 *
# *                                                                            *
# ******************************************************************************
import numpy as np
import matplotlib.pyplot as plt

COLORS = [
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#6B7280",
    "#332288",
    "#882255",
    "#44AA99",
]


def get_body_force(x):
    """Return the axial load per unit length at position x."""
    return x**2.0


def exact_stress(x, length, end_load, area):
    """Return the analytical axial stress at position x."""
    axial_force = end_load + length**3 / 3 - x**3 / 3
    return axial_force / area


def exact_displacement(x, youngs_modulus, area, length, end_load):
    """Return the analytical axial displacement at position x."""
    numerator = (end_load + length**3 / 3) * x - x**4 / 12
    return numerator / (youngs_modulus * area)


def get_mesh_bar(nelem_1D, L, P):
    nnode = nelem_1D + 1
    coords = np.linspace(0.0, L, nnode)
    tolerance = 1.0e-10

    connect = np.zeros((nelem_1D, 2), dtype=int)
    for i in range(nelem_1D):
        connect[i, 0] = i
        connect[i, 1] = i + 1

    fixed_nodes = np.array([0], dtype=int)
    fixed_values = np.array([0.0], dtype=float)

    nodal_forces = np.zeros(nnode, dtype=float)
    nodal_forces[-1] = P

    return nnode, coords, connect, fixed_nodes, fixed_values, nodal_forces


def shapefunction_1d_two_nodes(xi):
    """Return 2-node bar shape functions and their derivatives with respect to xi."""
    N = np.array([
        0.5 * (1 - xi),
        0.5 * (1 + xi)
    ])
    dNdxi = 0.5 * np.array([
        -1,
        1
    ])
    return N, dNdxi
        

def get_element_k_and_f(
    coord,
    E,
    A,
    gauss_points,
    gauss_weights,
    get_body_force,
):
    """Return the 1D linear bar element stiffness matrix and body-force vector."""
    ke = np.zeros((2 , 2))
    fe = np.zeros(2)

    for xi, weight_xi in zip(gauss_points, gauss_weights):

        N, dNdxi = shapefunction_1d_two_nodes(xi)

        # Isoparametric mapping: x(xi) = N(xi) @ coord.
        x_gauss = N @ coord
        J = dNdxi @ coord  # dx/dxi
        dNdx = dNdxi / J
        
        ke += E * A * np.outer(dNdx, dNdx) * J * weight_xi
        fe += N * get_body_force(x_gauss) * J * weight_xi


    return ke, fe


def get_gauss_integration_points(number_gauss_points_in_one_direction):
    """Return Gauss-Legendre points and weights on the interval [-1, 1]."""
    if not isinstance(number_gauss_points_in_one_direction, int):
        raise TypeError("The number of Gauss points must be an integer.")
    if number_gauss_points_in_one_direction < 1:
        raise ValueError("The number of Gauss points must be at least 1.")

    return np.polynomial.legendre.leggauss(number_gauss_points_in_one_direction)

def get_global_stiffness_and_force(
    total_nodes,
    coords,
    total_element,
    elements,
    E,
    A,
    gauss_points,
    gauss_weights,
    get_body_force,
):
    """Assemble the global stiffness matrix K and body-force vector F."""

    K = np.zeros((total_nodes,total_nodes))
    F = np.zeros(total_nodes)

    for e in range(total_element):
        # Node IDs in the e-th element, e.g. [0, 1] 
        element = elements[e]

        # node coordinates for the e-th element
        element_coordinates = coords[element]

        # calculate the elment k and f for the e-th element
        ke, fe = get_element_k_and_f(
            element_coordinates,
            E,
            A,
            gauss_points,
            gauss_weights,
            get_body_force,
        )

        for i in range(2):
            # global node id
            row = element[i]
            F[row] += fe[i]
            for j in range(2):
                # global node id for the b-th local node
                column = element[j]
                # local node id for the b-th local node
                K[row, column] += ke[i, j]

    return K, F


def apply_traction_boundary_conditions(F, nodal_forces):
    """Apply nodal forces to the global force vector F."""
    F += nodal_forces
    return F


def apply_essential_boundary_conditions(K, F, fixed_nodes, fixed_values):
    """Apply prescribed nodal displacements to K and F."""
    for node, value in zip(fixed_nodes, fixed_values):
        F -= K[:, node] * value
        K[node, :] = 0.0
        K[:, node] = 0.0
        K[node, node] = 1.0
        F[node] = value

    return K, F


def compute_element_displacement_strain_and_stress(
    coords,
    elements,
    U,
    E,
    number_sample_points=6,
):
    """Return elementwise samples without averaging shared boundaries."""
    if number_sample_points < 2:
        raise ValueError("Each element needs at least two sample points.")

    xi_samples = np.linspace(-1.0, 1.0, number_sample_points)
    total_samples = len(elements) * number_sample_points
    sample_coords = np.zeros(total_samples)
    sample_strain = np.zeros(total_samples)
    sample_stress = np.zeros(total_samples)
    disp_sample = np.zeros(total_samples)

    sample_index = 0
    for element in elements:
        element_coords = coords[element]
        element_displacements = U[element]

        for xi in xi_samples:
            N, dNdxi = shapefunction_1d_two_nodes(xi)

            # Use the same isoparametric mapping as in Part 6.
            x_sample = N @ element_coords
            disp_sample[sample_index] = N @ element_displacements
            J = dNdxi @ element_coords  # dx/dxi
            if J <= 0.0:
                raise ValueError(
                    "Element Jacobian must be positive."
                )
            dNdx = dNdxi / J

            sample_coords[sample_index] = x_sample
            sample_strain[sample_index] = (
                dNdx @ element_displacements
            )
            sample_stress[sample_index] = (
                E * sample_strain[sample_index]
            )
            sample_index += 1

    return sample_coords, disp_sample, sample_strain, sample_stress



def calculate_displacement_and_energy_norm_error(
    coordinates,
    elements,
    displacements,
    youngs_modulus,
    area,
    length,
    end_load,
):
    """Return the L2 displacement error and the energy-norm error."""
    gauss_points, gauss_weights = get_gauss_integration_points(4)
    displacement_error_squared = 0.0
    energy_error_squared = 0.0
    for element in elements:
        element_coordinates = coordinates[element]
        element_displacements = displacements[element]
        for xi, weight in zip(gauss_points, gauss_weights):
            N, dNdxi = shapefunction_1d_two_nodes(xi)
            x = N @ element_coordinates
            jacobian = dNdxi @ element_coordinates
            dNdx = dNdxi / jacobian
            displacement_fem = N @ element_displacements
            displacement_exact = exact_displacement(
                x,
                youngs_modulus,
                area,
                length,
                end_load,
            )
            displacement_error = displacement_fem - displacement_exact

            strain_fem = dNdx @ element_displacements
            strain_exact = exact_stress(x, length, end_load, area) / youngs_modulus
            strain_error = strain_fem - strain_exact

            displacement_error_squared += displacement_error**2 * jacobian * weight
            energy_error_squared += (
                youngs_modulus * area * strain_error**2 * jacobian * weight
            )
    return np.sqrt(displacement_error_squared), np.sqrt(energy_error_squared)


def plot_element_variables(
    sample_x,
    sample_y,
    exact_x,
    exact_y,
    fem_node_coords,
    label_input,
    title,
    output_path=None,
):
    """Plot elementwise FEM data and optionally save the figure."""

    colors = [
        '#0072B2',  # blue
        '#D55E00',  # vermilion
        '#009E73',  # green
        '#CC79A7',  # magenta
        '#E69F00',  # gold
        '#56B4E9',  # sky blue
        '#6B7280',  # gray
        '#332288',  # indigo
        '#882255',  # wine
        '#44AA99',  # teal
    ]

    plt.rcParams['font.family'] = 'serif' # Can be 'sans-serif', 'monospace', 'cursive', 'fantasy'
    plt.rcParams['font.style'] = 'normal' # Can be 'normal', 'italic', or 'oblique'
    plt.rcParams['font.weight'] = 'bold'  # Can be 'normal', 'bold', 'heavy', 'light', etc.
    plt.rcParams['font.size'] = 18        # Font size in points


    plt.figure(figsize=(10, 7.5))

    plt.plot(
        sample_x,
        sample_y,
        '-',
        color=colors[1],
        linewidth=2.5,
        label=label_input,
    )
    plt.plot(
        exact_x,
        exact_y,
        '--',
        linewidth=1.5,
        color=colors[0],
        label='Analytical',
    )
    plt.plot(
        fem_node_coords,
        np.zeros_like(fem_node_coords),
        'o',
        markersize=10,
        color=colors[6],
        label='Nodes',
    )
    plt.xlabel('x')
    plt.ylabel(label_input)
    plt.title(title)

    plt.legend()

    plt.tight_layout()
    if output_path is not None:
        plt.savefig(output_path, dpi=300)
    plt.show()
