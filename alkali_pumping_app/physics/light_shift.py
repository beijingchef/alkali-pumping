"""Detuning sweeps of diagonal and full-manifold AC-Stark shifts."""

from math import pi

import numpy as np

from .angular_momentum import (
    build_excited_states,
    build_ground_states,
    dipole_amplitude,
)
from .optical_pumping import optical_rate_scale_from_intensity
from .polarization import (
    spherical_components_from_lab,
    spherical_weights_from_lab,
    tensor_geometry_E20,
    transverse_basis_for_k,
)
from .spectroscopy import (
    complex_voigt_response_relative,
    doppler_fwhm_MHz,
    transition_shift_MHz,
)


def stokes_from_lab_vector(E_lab, k_axis):
    """Return normalized beam-frame Stokes parameters for a Jones vector."""
    vector = np.asarray(E_lab, dtype=complex)
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise ValueError("Polarization vector has zero norm.")
    vector = vector / norm
    e1, e2 = transverse_basis_for_k(k_axis)
    a = np.vdot(e1, vector)
    b = np.vdot(e2, vector)
    return np.array(
        [
            abs(a) ** 2 - abs(b) ** 2,
            2.0 * np.real(np.conj(a) * b),
            2.0 * np.imag(np.conj(a) * b),
        ],
        dtype=float,
    )


def _dipole_couplings(atom, line, ground_states, excited_states, E_lab, q_axis):
    """Return <e|d.epsilon|g> without the common electronic reduced element."""
    spherical = spherical_components_from_lab(E_lab, q_axis)
    I = atom["I"]
    Jg = atom["ground"]["J"]
    Je = atom[line]["Jp"]
    coupling = np.zeros((len(excited_states), len(ground_states)), dtype=complex)
    for ei, excited in enumerate(excited_states):
        for gi, ground in enumerate(ground_states):
            q_float = float(excited["m"]) - float(ground["m"])
            q = int(round(q_float))
            if abs(q_float - q) > 1e-9 or q not in spherical:
                continue
            coupling[ei, gi] = spherical[q] * dipole_amplitude(
                I,
                Jg,
                Je,
                ground["F"],
                ground["m"],
                excited["F"],
                excited["m"],
                q,
            )
    return coupling


def _component_decomposition(ground_states, diagonal_hz):
    """Decompose every detuning into scalar/vector/tensor state contributions."""
    diagonal_hz = np.asarray(diagonal_hz, dtype=float)
    point_count, state_count = diagonal_hz.shape
    scalar = np.zeros_like(diagonal_hz)
    vector = np.zeros_like(diagonal_hz)
    tensor = np.zeros_like(diagonal_hz)
    residual = np.zeros_like(diagonal_hz)
    coefficient_rows = []

    for F in sorted({float(state["F"]) for state in ground_states}):
        indices = np.array(
            [
                index
                for index, state in enumerate(ground_states)
                if np.isclose(float(state["F"]), F)
            ],
            dtype=int,
        )
        m = np.array([float(ground_states[index]["m"]) for index in indices])
        basis = np.column_stack(
            [
                np.ones_like(m),
                m,
                3.0 * m**2 - F * (F + 1.0),
            ]
        )
        coefficients = np.linalg.pinv(basis) @ diagonal_hz[:, indices].T
        fitted = (basis @ coefficients).T
        scalar[:, indices] = coefficients[0, :, None]
        vector[:, indices] = coefficients[1, :, None] * m[None, :]
        tensor[:, indices] = coefficients[2, :, None] * basis[None, :, 2]
        residual[:, indices] = diagonal_hz[:, indices] - fitted
        coefficient_rows.append(
            {
                "F": F,
                "scalar": coefficients[0],
                "vector_per_m": coefficients[1],
                # This is the actual m=0 tensor frequency shift used in the
                # Mathur convention discussed alongside this feature.
                "tensor_m0": -F * (F + 1.0) * coefficients[2],
                "residual_max": np.max(np.abs(diagonal_hz[:, indices] - fitted), axis=1),
            }
        )

    return {
        "scalar": scalar,
        "vector": vector,
        "tensor": tensor,
        "residual": residual,
        "coefficients": coefficient_rows,
        "point_count": point_count,
        "state_count": state_count,
    }


def manifold_light_shift_eigenvalues(
    ground_states,
    hamiltonian_hz_per_uW_cm2,
    bare_zeeman_hz,
    intensity_uW_cm2=1.0,
):
    """Diagonalize each F block at the requested physical light intensity."""
    matrices = np.asarray(hamiltonian_hz_per_uW_cm2, dtype=complex)
    bare = np.asarray(bare_zeeman_hz, dtype=float)
    intensity = float(intensity_uW_cm2)
    if intensity < 0.0:
        raise ValueError("Light intensity must be nonnegative.")
    if matrices.ndim != 3 or matrices.shape[1:] != (len(ground_states), len(ground_states)):
        raise ValueError("Hamiltonian sweep shape does not match the ground states.")
    if bare.shape != (len(ground_states),):
        raise ValueError("bare_zeeman_hz must contain one value per ground state.")

    rows = []
    for F in sorted({float(state["F"]) for state in ground_states}):
        indices = np.array(
            [
                index
                for index, state in enumerate(ground_states)
                if np.isclose(float(state["F"]), F)
            ],
            dtype=int,
        )
        bare_sorted = np.sort(bare[indices])
        shifts = np.empty((len(matrices), len(indices)), dtype=float)
        for point_index, matrix in enumerate(matrices):
            block = (
                intensity * matrix[np.ix_(indices, indices)]
                + np.diag(bare[indices])
            )
            shifts[point_index] = np.linalg.eigvalsh(block) - bare_sorted
        rows.append({"F": F, "shifts": shifts})
    return rows


def calculate_light_shift_sweep(
    atom,
    line,
    detunings_MHz,
    E_lab,
    k_axis,
    q_axis,
    n2_pressure_torr,
    temperature_C,
    n2_width_MHz_per_torr,
    n2_shift_MHz_per_torr,
    bare_zeeman_hz=None,
):
    """Calculate per-intensity light shifts over a laser-detuning array.

    All returned shifts are normalized to 1 microW/cm^2.  Diagonal state
    contributions are reported in the selected |F,m> basis.  Separately, each
    F-manifold block is diagonalized so arbitrary pure polarization remains
    meaningful when Raman terms mix m states.  Couplings between different
    ground hyperfine manifolds are omitted under the usual secular hyperfine
    approximation.
    """
    detunings = np.asarray(detunings_MHz, dtype=float)
    if detunings.ndim != 1 or len(detunings) < 2:
        raise ValueError("detunings_MHz must be a one-dimensional sweep.")
    if line not in ("D1", "D2"):
        raise ValueError("line must be D1 or D2.")

    vector = np.asarray(E_lab, dtype=complex)
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise ValueError("Polarization vector has zero norm.")
    vector = vector / norm

    ground_states = build_ground_states(atom)
    excited_states = build_excited_states(atom, line)
    coupling = _dipole_couplings(
        atom, line, ground_states, excited_states, vector, q_axis
    )

    pressure_width = float(n2_width_MHz_per_torr) * float(n2_pressure_torr)
    pressure_shift = float(n2_shift_MHz_per_torr) * float(n2_pressure_torr)
    lorentz_fwhm = float(atom[line]["gamma_nat_MHz"]) + pressure_width
    doppler_fwhm = doppler_fwhm_MHz(atom, line, float(temperature_C))
    transition_centers = np.array(
        [
            [transition_shift_MHz(ground, excited) + pressure_shift for excited in excited_states]
            for ground in ground_states
        ],
        dtype=float,
    )
    response = complex_voigt_response_relative(
        detunings[:, None, None] - transition_centers[None, :, :],
        lorentz_fwhm,
        doppler_fwhm,
    )
    absorption = np.maximum(0.0, np.real(response))
    dispersion = np.imag(response)
    rate_scale = optical_rate_scale_from_intensity(
        atom=atom,
        line=line,
        intensity_uW_cm2=1.0,
        n2_pressure_torr=n2_pressure_torr,
        temperature_C=temperature_C,
        n2_width_MHz_per_torr=n2_width_MHz_per_torr,
    )

    state_count = len(ground_states)
    hamiltonian_hz = np.zeros(
        (len(detunings), state_count, state_count), dtype=complex
    )
    for gi, ground in enumerate(ground_states):
        for gj, other in enumerate(ground_states):
            if not np.isclose(float(ground["F"]), float(other["F"])):
                continue
            amplitude_product = np.conj(coupling[:, gj]) * coupling[:, gi]
            symmetrized_dispersion = 0.25 * (
                dispersion[:, gi, :] + dispersion[:, gj, :]
            )
            hamiltonian_hz[:, gj, gi] = (
                rate_scale
                * np.sum(symmetrized_dispersion * amplitude_product[None, :], axis=1)
                / (2.0 * pi)
            )
    hamiltonian_hz = 0.5 * (
        hamiltonian_hz + np.swapaxes(hamiltonian_hz.conj(), 1, 2)
    )

    diagonal_hz = np.real(np.diagonal(hamiltonian_hz, axis1=1, axis2=2))
    components = _component_decomposition(ground_states, diagonal_hz)
    scattering_rate = rate_scale * np.sum(
        absorption * (np.abs(coupling.T)[None, :, :] ** 2), axis=2
    )

    if bare_zeeman_hz is None:
        bare = np.zeros(state_count, dtype=float)
    else:
        bare = np.asarray(bare_zeeman_hz, dtype=float)
        if bare.shape != (state_count,):
            raise ValueError("bare_zeeman_hz must contain one value per ground state.")

    eigenvalue_rows = manifold_light_shift_eigenvalues(
        ground_states, hamiltonian_hz, bare, intensity_uW_cm2=1.0
    )

    spherical_weights = spherical_weights_from_lab(vector, q_axis)
    populated_components = sum(value > 1e-10 for value in spherical_weights.values())
    return {
        "detunings_MHz": detunings,
        "ground_states": ground_states,
        "excited_states": excited_states,
        "hamiltonian_hz_per_uW_cm2": hamiltonian_hz,
        "diagonal_hz_per_uW_cm2": diagonal_hz,
        "components_hz_per_uW_cm2": components,
        "eigenvalues_hz_per_uW_cm2": eigenvalue_rows,
        "bare_zeeman_hz": bare,
        "scattering_s_inv_per_uW_cm2": scattering_rate,
        "spherical_weights": spherical_weights,
        "stokes": stokes_from_lab_vector(vector, k_axis),
        "E20": tensor_geometry_E20(vector, q_axis),
        "diagonal_in_selected_basis": populated_components <= 1,
        "lorentz_fwhm_MHz": lorentz_fwhm,
        "doppler_fwhm_MHz": doppler_fwhm,
        "pressure_shift_MHz": pressure_shift,
    }
