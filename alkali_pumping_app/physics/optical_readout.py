"""Weak-probe Stokes readout operators for RF-driven ground-state coherences.

The probe is a detector only: it does not enter the optical-pumping generator,
add a light shift, or broaden the magnetic resonance.  The transmitted Stokes
signal is linearized both in the vapor polarizability and in the RF-induced
density-matrix perturbation.  This keeps the scalar, orientation, and alignment
readout contributions additive at the complex-susceptibility level.
"""

from math import pi

import numpy as np

from .atomic_polarizability import (
    calculate_atomic_polarizability_sweep,
    mathur_quadrupole_scale,
)
from .polarization import polarization_ellipse_vector, transverse_basis_for_k
from .rf_response import _spin_operator


_LAB_AXES = ("x", "y", "z")
_AXIS_INDEX = {axis: index for index, axis in enumerate(_LAB_AXES)}


def _levi_civita(i, j, k):
    if len({i, j, k}) < 3:
        return 0.0
    return 1.0 if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1.0


def _manifold_projector(ground_states, target_F):
    size = len(ground_states)
    projector = np.zeros((size, size), dtype=complex)
    for index, state in enumerate(ground_states):
        if np.isclose(float(state["F"]), float(target_F)):
            projector[index, index] = 1.0
    return projector


def _electronic_spin_projection_factor(atom, target_F):
    """Return the coefficient in P_F J P_F = coefficient * P_F F P_F."""
    I = float(atom["I"])
    J = float(atom["ground"]["J"])
    F = float(target_F)
    denominator = 2.0 * F * (F + 1.0)
    if denominator <= 0.0:
        return 0.0
    return (
        F * (F + 1.0) + J * (J + 1.0) - I * (I + 1.0)
    ) / denominator


def _ij_eigenvalue(atom, target_F):
    I = float(atom["I"])
    J = float(atom["ground"]["J"])
    F = float(target_F)
    return 0.5 * (
        F * (F + 1.0) - I * (I + 1.0) - J * (J + 1.0)
    )


def _single_detuning_polarizability(
    atom,
    line,
    detuning_MHz,
    n2_pressure_torr,
    temperature_C,
    n2_width_MHz_per_torr,
    n2_shift_MHz_per_torr,
):
    # The public sweep routine requires at least two samples. Repeating the
    # requested detuning avoids introducing an artificial finite difference.
    sweep = calculate_atomic_polarizability_sweep(
        atom=atom,
        line=line,
        detunings_MHz=np.array([detuning_MHz, detuning_MHz], dtype=float),
        n2_pressure_torr=n2_pressure_torr,
        temperature_C=temperature_C,
        n2_width_MHz_per_torr=n2_width_MHz_per_torr,
        n2_shift_MHz_per_torr=n2_shift_MHz_per_torr,
    )
    return {
        "alpha_eq": complex(sweep["alpha_eq"][0]),
        "alpha_hfs": complex(sweep["alpha_hfs"][0]),
        "alpha_gt": {F: complex(values[0]) for F, values in sweep["alpha_gt"].items()},
        "alpha_br": {F: complex(values[0]) for F, values in sweep["alpha_br"].items()},
        "lorentz_fwhm_MHz": float(sweep["lorentz_fwhm_MHz"]),
        "doppler_fwhm_MHz": float(sweep["doppler_fwhm_MHz"]),
    }


def polarizability_tensor_operators(
    atom,
    ground_states,
    q_axis,
    response,
):
    """Return scalar, orientation, and alignment polarizability operators.

    Each returned array has shape ``(3, 3, N, N)``. The first two indices are
    laboratory optical-polarization coordinates; the final two are ground-state
    density-matrix coordinates.  Diagonal hyperfine response functions are used,
    consistently with the Atomic polarizability page and the absence of ground-
    hyperfine optical coherences in the current model.

    The alignment operator uses the same raw Cartesian convention as the
    existing atomic response, Q_ij=(F_i F_j+F_j F_i)/2, with its scalar trace
    removed.  Mathur's normalized rank-2 polarizability is converted to that
    convention, including its opposite sign.
    """
    size = len(ground_states)
    tensors = {
        name: np.zeros((3, 3, size, size), dtype=complex)
        for name in ("scalar", "orientation", "alignment")
    }
    manifolds = sorted({float(state["F"]) for state in ground_states})

    for F in manifolds:
        projector = _manifold_projector(ground_states, F)
        spin = {
            axis: _spin_operator(ground_states, q_axis, axis, F)
            for axis in _LAB_AXES
        }

        scalar_coefficient = (
            response["alpha_eq"] + response["alpha_hfs"] * _ij_eigenvalue(atom, F)
        )
        for axis_index in range(3):
            tensors["scalar"][axis_index, axis_index] += (
                scalar_coefficient * projector
            )

        electronic_factor = _electronic_spin_projection_factor(atom, F)
        alpha_gt = response["alpha_gt"][F]
        for i in range(3):
            for j in range(3):
                for k, axis in enumerate(_LAB_AXES):
                    epsilon = _levi_civita(i, k, j)
                    if epsilon:
                        tensors["orientation"][i, j] += (
                            1j * alpha_gt * epsilon * electronic_factor * spin[axis]
                        )

        alpha_br = response["alpha_br"][F] * mathur_quadrupole_scale(F)
        for i, first_axis in enumerate(_LAB_AXES):
            for j, second_axis in enumerate(_LAB_AXES):
                quadrupole = 0.5 * (
                    spin[first_axis] @ spin[second_axis]
                    + spin[second_axis] @ spin[first_axis]
                )
                if i == j:
                    quadrupole -= F * (F + 1.0) * projector / 3.0
                tensors["alignment"][i, j] += alpha_br * quadrupole

    return tensors


def _transverse_tensor(tensor, k_axis):
    first, second = transverse_basis_for_k(k_axis)
    basis = np.asarray([first.real, second.real], dtype=float)
    return np.einsum("ai,ijmn,bj->abmn", basis, tensor, basis)


def _stokes_matrices():
    return {
        "S0": np.array([[1.0, 0.0], [0.0, 1.0]], dtype=complex),
        "S1": np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
        "S2": np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
        # This convention gives S3=+1 for sigma+ as defined in polarization.py.
        "S3": np.array([[0.0, -1j], [1j, 0.0]], dtype=complex),
    }


def _linear_stokes_operators(transverse_tensor, input_jones, coupling):
    """Return Hermitian operators for the first-order changes in S0..S3."""
    operators = {}
    for name, stokes_matrix in _stokes_matrices().items():
        left = np.conjugate(input_jones) @ stokes_matrix
        kernel = 1j * float(coupling) * np.einsum(
            "a,b,abmn->mn", left, input_jones, transverse_tensor
        )
        operator = kernel + np.conjugate(kernel.T)
        operators[name] = 0.5 * (operator + np.conjugate(operator.T))
    return operators


def _derived_signal_operators(stokes_operators, input_stokes):
    s0, s1, s2, s3 = (float(value) for value in input_stokes)
    if s0 <= 0.0:
        raise ValueError("The probe input Stokes intensity must be positive.")
    normalized = np.array([s1 / s0, s2 / s0, s3 / s0], dtype=float)
    result = {"transmission": stokes_operators["S0"] / s0}
    for index, name in enumerate(("s1", "s2", "s3"), start=1):
        result[name] = (
            stokes_operators[f"S{index}"]
            - normalized[index - 1] * stokes_operators["S0"]
        ) / s0

    linear_norm_sq = normalized[0] ** 2 + normalized[1] ** 2
    if linear_norm_sq > 1e-14:
        result["rotation"] = 0.5 * (
            normalized[0] * result["s2"] - normalized[1] * result["s1"]
        ) / linear_norm_sq
        rotation_available = True
    else:
        result["rotation"] = np.zeros_like(result["s1"])
        rotation_available = False

    ellipticity_denominator = max(0.0, 1.0 - normalized[2] ** 2) ** 0.5
    if ellipticity_denominator > 1e-7:
        result["ellipticity"] = 0.5 * result["s3"] / ellipticity_denominator
        ellipticity_available = True
    else:
        result["ellipticity"] = np.zeros_like(result["s3"])
        ellipticity_available = False
    return result, {
        "rotation_available": rotation_available,
        "ellipticity_available": ellipticity_available,
    }


def weak_probe_readout_operators(
    *,
    atom,
    ground_states,
    q_axis,
    line,
    detuning_MHz,
    k_axis,
    azimuth_deg,
    ellipticity_deg,
    path_length_cm,
    density_cm3,
    n2_pressure_torr,
    temperature_C,
    n2_width_MHz_per_torr,
    n2_shift_MHz_per_torr,
    include_scalar=False,
    include_orientation=True,
    include_alignment=True,
):
    """Build additive weak-probe readout operators for six optical signals."""
    path_length_cm = float(path_length_cm)
    density_cm3 = float(density_cm3)
    if path_length_cm < 0.0:
        raise ValueError("Probe path length must be nonnegative.")
    if density_cm3 < 0.0:
        raise ValueError("Alkali density must be nonnegative.")

    response = _single_detuning_polarizability(
        atom,
        line,
        detuning_MHz,
        n2_pressure_torr,
        temperature_C,
        n2_width_MHz_per_torr,
        n2_shift_MHz_per_torr,
    )
    tensor_parts = polarizability_tensor_operators(
        atom, ground_states, q_axis, response
    )

    input_lab = polarization_ellipse_vector(k_axis, azimuth_deg, ellipticity_deg)
    first, second = transverse_basis_for_k(k_axis)
    input_jones = np.array(
        [np.vdot(first, input_lab), np.vdot(second, input_lab)], dtype=complex
    )
    stokes_matrices = _stokes_matrices()
    input_stokes = np.array(
        [float(np.real(np.vdot(input_jones, matrix @ input_jones))) for matrix in stokes_matrices.values()]
    )

    wavelength_cm = float(atom[f"lambda_{line}_nm"]) * 1e-7
    # Mathur polarizabilities are volumes in Gaussian-cgs units: n-1=2*pi*N*alpha.
    coupling = 4.0 * pi**2 * density_cm3 * path_length_cm / wavelength_cm
    component_operators = {}
    component_info = {}
    for component, tensor in tensor_parts.items():
        transverse = _transverse_tensor(tensor, k_axis)
        stokes_ops = _linear_stokes_operators(transverse, input_jones, coupling)
        signal_ops, availability = _derived_signal_operators(stokes_ops, input_stokes)
        component_operators[component] = signal_ops
        component_info[component] = availability

    enabled = {
        "scalar": bool(include_scalar),
        "orientation": bool(include_orientation),
        "alignment": bool(include_alignment),
    }
    total = {
        signal: sum(
            (
                component_operators[component][signal]
                for component in ("scalar", "orientation", "alignment")
                if enabled[component]
            ),
            np.zeros_like(next(iter(component_operators["scalar"].values()))),
        )
        for signal in component_operators["scalar"]
    }
    operators = {"total": total, **component_operators}
    return operators, {
        "line": line,
        "detuning_MHz": float(detuning_MHz),
        "k_axis": k_axis,
        "azimuth_deg": float(azimuth_deg),
        "ellipticity_deg": float(ellipticity_deg),
        "path_length_cm": path_length_cm,
        "density_cm3": density_cm3,
        "coupling": float(coupling),
        "input_stokes": input_stokes,
        "enabled_contributions": enabled,
        "lorentz_fwhm_MHz": response["lorentz_fwhm_MHz"],
        "doppler_fwhm_MHz": response["doppler_fwhm_MHz"],
        "signal_availability": component_info["orientation"],
        "model": "first-order weak-probe Mathur polarizability readout",
    }
