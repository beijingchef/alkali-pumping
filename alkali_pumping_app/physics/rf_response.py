"""Weak-RF adjacent-coherence response calculations."""

from math import sqrt

import numpy as np

from .spin_exchange import coupled_basis_amplitudes

# ============================================================
# 7c. Weak linearly polarized RF response
# ============================================================

def lab_axis_in_local_frame(q_axis, lab_axis):
    """Map a laboratory axis to the local frame whose z axis is q_axis.

    This uses the same right-handed local frames as local_components():
      q=z: (x_local,y_local,z_local) = (x,y,z)
      q=x: (x_local,y_local,z_local) = (y,z,x)
      q=y: (x_local,y_local,z_local) = (z,x,y)
    """
    mappings = {
        "z": {"x": "x", "y": "y", "z": "z"},
        "x": {"x": "z", "y": "x", "z": "y"},
        "y": {"x": "y", "y": "z", "z": "x"},
    }
    return mappings[q_axis][lab_axis]


def quadrupole_axes_for_rf_axis(rf_axis):
    """Return the cyclic off-diagonal quadrupole axes for one RF axis."""
    axes = {
        "x": ("y", "z"),
        "y": ("z", "x"),
        "z": ("x", "y"),
    }
    try:
        return axes[str(rf_axis).lower()]
    except KeyError as exc:
        raise ValueError(f"Unsupported RF axis: {rf_axis}") from exc


def rf_observable_display_label(observable, rf_axis):
    """Return the laboratory-frame symbol shown for an RF observable."""
    if observable == "Q_ij":
        first, second = quadrupole_axes_for_rf_axis(rf_axis)
        return f"Q_{first}{second}"
    return str(observable)


def adjacent_spin_matrix_elements(local_axis, ladder_coefficient):
    """Return O_ab and O_ba for a=|F,m>, b=|F,m-1>."""
    C = float(ladder_coefficient)
    if local_axis == "x":
        return 0.5 * C, 0.5 * C
    if local_axis == "y":
        return -0.5j * C, 0.5j * C
    if local_axis == "z":
        return 0.0j, 0.0j
    raise ValueError(local_axis)


def weak_rf_observable_susceptibility(
    frequencies_hz,
    ground_states,
    populations,
    adjacent_transition_hz,
    gamma_op,
    gamma_er,
    gamma_se,
    q_axis,
    rf_axis,
    observable,
    target_F=None,
):
    """Return the zero-drive RF susceptibility of one hyperfine manifold.

    Only adjacent coherences within ``target_F`` contribute. If ``target_F``
    is omitted, the upper ground hyperfine manifold (largest F) is used. The
    returned observable is therefore an upper-manifold spin or alignment
    component

        <O_+> = Tr[rho_{F+} O_+],

    not the coherent sum of the two ground hyperfine spin observables.
    For ``observable='Q_ij'``, the off-diagonal alignment operator is

        Q_ij = (F_i F_j + F_j F_i) / 2,

    with (i,j) selected cyclically from the RF axis: x -> yz, y -> zx,
    and z -> xy.

    The applied field is linearly polarized along the selected laboratory
    axis i = x, y, or z:

        H_rf / hbar = Omega_rf F_i cos(omega t).

    The density matrix is expanded directly to first order,

        rho = rho_0 + Omega_rf rho^(1) + O(Omega_rf^2),

    and this function calculates rho^(1) without assigning any finite numerical
    value to Omega_rf. Let chi_plus be the coefficient of exp(-i omega t).
    The real time-domain susceptibility is written as

        d<F_{+,i}(t)>/dOmega_rf
        = X(omega) cos(omega t) + Y(omega) sin(omega t),

    with

        X = 2 Re[chi_plus],
        Y = 2 Im[chi_plus],
        A = sqrt(X^2 + Y^2) = 2 |chi_plus|.

    The function returns A, X, and Y within the linear weak-drive
    approximation.

    Every adjacent coherence is treated independently with its local rate

        Gamma_m = Gamma_OP,m + Gamma_ER,m + Gamma_SE,m.

    Both Fourier components of the real cosine field are retained. This is
    important when a negative signed transition frequency contributes to a
    positive-frequency resonance through the counter-rotating component.
    """
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    p = np.asarray(populations, dtype=float)
    transition_hz = np.asarray(adjacent_transition_hz, dtype=float)
    gamma_op = np.asarray(gamma_op, dtype=float)
    gamma_er = np.asarray(gamma_er, dtype=float)
    gamma_se = np.asarray(gamma_se, dtype=float)

    susceptibility_amplitude = np.zeros_like(frequencies_hz, dtype=float)
    susceptibility_in_phase = np.zeros_like(frequencies_hz, dtype=float)
    susceptibility_quadrature = np.zeros_like(frequencies_hz, dtype=float)
    susceptibility_phasor = np.zeros_like(frequencies_hz, dtype=complex)

    if len(ground_states) == 0:
        return (
            susceptibility_amplitude,
            susceptibility_in_phase,
            susceptibility_quadrature,
            {
                "used_transitions": 0,
                "nonpositive_linewidths": 0,
                "target_F": np.nan,
            },
        )

    if target_F is None:
        target_F = max(float(state["F"]) for state in ground_states)
    target_F = float(target_F)

    drive_local_axis = lab_axis_in_local_frame(q_axis, rf_axis)
    drive_operator = _spin_operator(
        ground_states, q_axis, rf_axis, target_F
    )
    observable_operator, observable_info = _rf_observable_operator(
        ground_states,
        q_axis,
        rf_axis,
        observable,
        target_F,
    )

    state_index = {
        (float(state["F"]), float(state["m"])): idx
        for idx, state in enumerate(ground_states)
    }

    omega = 2.0 * np.pi * frequencies_hz
    used_transitions = 0
    nonpositive_linewidths = 0

    for a_idx, state in enumerate(ground_states):
        F = float(state["F"])
        if not np.isclose(F, target_F):
            continue
        m = float(state["m"])
        b_idx = state_index.get((F, m - 1.0))
        if b_idx is None or not np.isfinite(transition_hz[a_idx]):
            continue

        C2 = F * (F + 1.0) - m * (m - 1.0)
        if C2 <= 0:
            continue
        C = sqrt(C2)

        drive_ab = drive_operator[a_idx, b_idx]
        drive_ba = drive_operator[b_idx, a_idx]
        obs_ab = observable_operator[a_idx, b_idx]
        obs_ba = observable_operator[b_idx, a_idx]
        if abs(drive_ab) == 0.0 or (abs(obs_ab) == 0.0 and abs(obs_ba) == 0.0):
            continue

        gamma = gamma_op[a_idx] + gamma_er[a_idx] + gamma_se[a_idx]
        if not np.isfinite(gamma):
            continue
        if gamma <= 0.0:
            nonpositive_linewidths += 1
            gamma = 1e-12

        D_m = p[a_idx] - p[b_idx]
        omega_m = 2.0 * np.pi * transition_hz[a_idx]

        # Coefficients of exp(-i omega t) in d(rho_ab)/dOmega_rf and
        # d(rho_ba)/dOmega_rf at Omega_rf=0. The factor 1/2 is the Fourier
        # amplitude of cos(omega t); no finite drive amplitude is introduced.
        drho_ab_plus_domega_rf = (
            1j * 0.5 * drive_ab * D_m
            / (gamma + 1j * (omega_m - omega))
        )
        drho_ba_plus_domega_rf = (
            -1j * 0.5 * drive_ba * D_m
            / (gamma + 1j * (-omega_m - omega))
        )

        susceptibility_phasor += (
            obs_ba * drho_ab_plus_domega_rf
            + obs_ab * drho_ba_plus_domega_rf
        )
        used_transitions += 1

    # With the RF reference chosen as cos(omega t),
    #
    #   d<O(t)>/dOmega_rf = X cos(omega t) + Y sin(omega t),
    #
    # where chi_plus is the coefficient of exp(-i omega t). Therefore
    # X = 2 Re(chi_plus), Y = 2 Im(chi_plus), and A = 2 |chi_plus|.
    susceptibility_in_phase = 2.0 * np.real(susceptibility_phasor)
    susceptibility_quadrature = 2.0 * np.imag(susceptibility_phasor)
    susceptibility_amplitude = np.hypot(
        susceptibility_in_phase,
        susceptibility_quadrature,
    )
    return (
        susceptibility_amplitude,
        susceptibility_in_phase,
        susceptibility_quadrature,
        {
            "used_transitions": used_transitions,
            "nonpositive_linewidths": nonpositive_linewidths,
            "drive_local_axis": drive_local_axis,
            "target_F": target_F,
            **observable_info,
        },
    )


def weak_rf_matrix_susceptibility(
    frequencies_hz,
    ground_states,
    populations,
    adjacent_transition_hz,
    gamma_op,
    gamma_er,
    gamma_se,
    q_axis,
    rf_axis,
    readout_operator,
    target_F=None,
):
    """Return weak-RF susceptibility for an arbitrary Hermitian readout.

    This is the matrix-readout counterpart of
    :func:`weak_rf_observable_susceptibility`.  The RF source remains restricted
    to adjacent coherences of ``target_F``; the readout matrix may represent an
    atomic observable or a weak optical measurement of those coherences.
    """
    frequencies_hz = np.asarray(frequencies_hz, dtype=float)
    populations = np.asarray(populations, dtype=float)
    transition_hz = np.asarray(adjacent_transition_hz, dtype=float)
    gamma_op = np.asarray(gamma_op, dtype=float)
    gamma_er = np.asarray(gamma_er, dtype=float)
    gamma_se = np.asarray(gamma_se, dtype=float)
    readout = np.asarray(readout_operator, dtype=complex)
    state_count = len(ground_states)
    if readout.shape != (state_count, state_count):
        raise ValueError("readout_operator must match the ground-state basis size.")

    phasor = np.zeros_like(frequencies_hz, dtype=complex)
    if state_count == 0:
        zeros = np.zeros_like(frequencies_hz, dtype=float)
        return zeros, zeros.copy(), zeros.copy(), {"used_transitions": 0}
    if target_F is None:
        target_F = max(float(state["F"]) for state in ground_states)
    target_F = float(target_F)
    drive = _spin_operator(ground_states, q_axis, rf_axis, target_F)
    state_index = {
        (float(state["F"]), float(state["m"])): index
        for index, state in enumerate(ground_states)
    }
    omega = 2.0 * np.pi * frequencies_hz
    used_transitions = 0
    nonpositive_linewidths = 0
    for a, state in enumerate(ground_states):
        F = float(state["F"])
        if not np.isclose(F, target_F):
            continue
        m = float(state["m"])
        b = state_index.get((F, m - 1.0))
        if b is None or not np.isfinite(transition_hz[a]):
            continue
        drive_ab = drive[a, b]
        drive_ba = drive[b, a]
        readout_ab = readout[a, b]
        readout_ba = readout[b, a]
        if abs(drive_ab) == 0.0 or (
            abs(readout_ab) == 0.0 and abs(readout_ba) == 0.0
        ):
            continue
        gamma = gamma_op[a] + gamma_er[a] + gamma_se[a]
        if not np.isfinite(gamma):
            continue
        if gamma <= 0.0:
            nonpositive_linewidths += 1
            gamma = 1e-12
        population_difference = populations[a] - populations[b]
        omega_m = 2.0 * np.pi * transition_hz[a]
        drho_ab = (
            0.5j * drive_ab * population_difference
            / (gamma + 1j * (omega_m - omega))
        )
        drho_ba = (
            -0.5j * drive_ba * population_difference
            / (gamma + 1j * (-omega_m - omega))
        )
        phasor += readout_ba * drho_ab + readout_ab * drho_ba
        used_transitions += 1

    in_phase = 2.0 * np.real(phasor)
    quadrature = 2.0 * np.imag(phasor)
    return (
        np.hypot(in_phase, quadrature),
        in_phase,
        quadrature,
        {
            "used_transitions": used_transitions,
            "nonpositive_linewidths": nonpositive_linewidths,
            "target_F": target_F,
        },
    )


def _density_reductions(amplitudes, density_matrix):
    """Return nuclear and electron reduced matrices in the uncoupled basis."""
    nuclear = np.einsum(
        "ais,ab,bjs->ij", amplitudes, density_matrix, amplitudes
    )
    electron = np.einsum(
        "ais,ab,bit->st", amplitudes, density_matrix, amplitudes
    )
    return nuclear, electron


def _recombine_density(amplitudes, nuclear, electron):
    return np.einsum(
        "ais,ij,st,bjt->ab", amplitudes, nuclear, electron, amplitudes
    )


def _ordered_coherences(state_count):
    return [(a, b) for a in range(state_count) for b in range(state_count) if a != b]


def _spin_operator(ground_states, q_axis, lab_axis, target_F):
    """Return one laboratory spin component within the selected F manifold."""
    size = len(ground_states)
    operator = np.zeros((size, size), dtype=complex)
    local_axis = lab_axis_in_local_frame(q_axis, lab_axis)
    state_index = {
        (float(state["F"]), float(state["m"])): index
        for index, state in enumerate(ground_states)
    }
    if local_axis == "z":
        for index, state in enumerate(ground_states):
            if np.isclose(float(state["F"]), float(target_F)):
                operator[index, index] = float(state["m"])
        return operator
    for a, state in enumerate(ground_states):
        F = float(state["F"])
        if not np.isclose(F, float(target_F)):
            continue
        m = float(state["m"])
        b = state_index.get((F, m - 1.0))
        if b is None:
            continue
        coefficient = sqrt(F * (F + 1.0) - m * (m - 1.0))
        ab, ba = adjacent_spin_matrix_elements(local_axis, coefficient)
        operator[a, b] = ab
        operator[b, a] = ba
    return operator


def quadrupole_operator(ground_states, q_axis, first_lab_axis, second_lab_axis, target_F):
    """Return Q_ij=(F_i F_j+F_j F_i)/2 in one F manifold."""
    first = _spin_operator(
        ground_states, q_axis, first_lab_axis, target_F
    )
    second = _spin_operator(
        ground_states, q_axis, second_lab_axis, target_F
    )
    return 0.5 * (first @ second + second @ first)


def _rf_observable_operator(ground_states, q_axis, rf_axis, observable, target_F):
    if observable == "Q_ij":
        first_axis, second_axis = quadrupole_axes_for_rf_axis(rf_axis)
        operator = quadrupole_operator(
            ground_states,
            q_axis,
            first_axis,
            second_axis,
            target_F,
        )
        return operator, {
            "observable_lab_axes": (first_axis, second_axis),
            "observable_local_axes": (
                lab_axis_in_local_frame(q_axis, first_axis),
                lab_axis_in_local_frame(q_axis, second_axis),
            ),
            "observable_label": rf_observable_display_label(observable, rf_axis),
        }

    observable_lab_axis = observable[-1].lower()
    operator = _spin_operator(
        ground_states, q_axis, observable_lab_axis, target_F
    )
    return operator, {
        "observable_lab_axis": observable_lab_axis,
        "observable_local_axis": lab_axis_in_local_frame(
            q_axis, observable_lab_axis
        ),
        "observable_label": rf_observable_display_label(observable, rf_axis),
    }


def _coherence_exchange_derivatives(target, partner):
    """Build target/self and partner-to-target density-map derivatives."""
    amplitudes_target = coupled_basis_amplitudes(
        target["atom"], target["ground_states"]
    )
    amplitudes_partner = coupled_basis_amplitudes(
        partner["atom"], partner["ground_states"]
    )
    rho_target = np.diag(np.asarray(target["population"], dtype=complex))
    rho_partner = np.diag(np.asarray(partner["population"], dtype=complex))
    nuclear_target, electron_target = _density_reductions(
        amplitudes_target, rho_target
    )
    _nuclear_partner, electron_partner = _density_reductions(
        amplitudes_partner, rho_partner
    )
    target_pairs = _ordered_coherences(len(target["ground_states"]))
    partner_pairs = _ordered_coherences(len(partner["ground_states"]))
    target_index = {pair: index for index, pair in enumerate(target_pairs)}

    er = np.zeros((len(target_pairs), len(target_pairs)), dtype=complex)
    self_exchange = np.zeros_like(er)
    cross_target = np.zeros_like(er)
    cross_partner = np.zeros(
        (len(target_pairs), len(partner_pairs)), dtype=complex
    )
    unpolarized_electron = 0.5 * np.eye(2, dtype=complex)

    def project(output, column, matrix):
        for pair, row in target_index.items():
            matrix[row, column] = output[pair]

    for column, (c, d) in enumerate(target_pairs):
        perturbation = np.zeros_like(rho_target)
        perturbation[c, d] = 1.0
        delta_nuclear, delta_electron = _density_reductions(
            amplitudes_target, perturbation
        )
        project(
            _recombine_density(
                amplitudes_target, delta_nuclear, unpolarized_electron
            ),
            column,
            er,
        )
        project(
            _recombine_density(amplitudes_target, delta_nuclear, electron_target)
            + _recombine_density(amplitudes_target, nuclear_target, delta_electron),
            column,
            self_exchange,
        )
        project(
            _recombine_density(amplitudes_target, delta_nuclear, electron_partner),
            column,
            cross_target,
        )

    for column, (c, d) in enumerate(partner_pairs):
        perturbation = np.zeros_like(rho_partner)
        perturbation[c, d] = 1.0
        _delta_nuclear, delta_electron = _density_reductions(
            amplitudes_partner, perturbation
        )
        project(
            _recombine_density(amplitudes_target, nuclear_target, delta_electron),
            column,
            cross_partner,
        )
    return target_pairs, er, self_exchange, cross_target, cross_partner


def _local_coherence_generator(result, pairs, er, self_exchange, cross_target):
    size = len(pairs)
    identity = np.eye(size, dtype=complex)
    state_frequency = (
        np.asarray(result["nu_LS"], dtype=float)
        + np.asarray(result["df_pop"]["nu_B"], dtype=float)
    )
    g_op = np.asarray(result["df_pop"]["G_OP"], dtype=float)
    diagonal = np.empty(size, dtype=complex)
    for index, (a, b) in enumerate(pairs):
        angular_frequency = 2.0 * np.pi * (state_frequency[a] - state_frequency[b])
        optical_decay = 0.5 * (g_op[a] + g_op[b])
        diagonal[index] = -optical_decay - 1j * angular_frequency
    return (
        np.diag(diagonal)
        + float(result["R_ER"]) * (er - identity)
        + float(result["R_SE_self"]) * (self_exchange - identity)
        + float(result["R_SE_cross"]) * (cross_target - identity)
    )


def _drive_and_observable_vectors(result, pairs, total_size, offset):
    target_F = result["rf_upper_F"]
    drive = _spin_operator(
        result["ground_states"], result["q_axis"], result["rf_axis"], target_F
    )
    observable, _observable_info = _rf_observable_operator(
        result["ground_states"],
        result["q_axis"],
        result["rf_axis"],
        result["rf_observable"],
        target_F,
    )
    population = np.asarray(result["population"], dtype=float)
    source = np.zeros(total_size, dtype=complex)
    readout = np.zeros(total_size, dtype=complex)
    for local_index, (a, b) in enumerate(pairs):
        index = offset + local_index
        source[index] = 0.5j * drive[a, b] * (population[a] - population[b])
        readout[index] = observable[b, a]
    return source, readout, int(np.count_nonzero(np.abs(source) > 1e-15))


def _drive_vector(result, pairs, total_size, offset):
    """Return the cosine-RF source vector for one species block."""
    drive = _spin_operator(
        result["ground_states"], result["q_axis"], result["rf_axis"],
        result["rf_upper_F"],
    )
    population = np.asarray(result["population"], dtype=float)
    source = np.zeros(total_size, dtype=complex)
    for local_index, (a, b) in enumerate(pairs):
        source[offset + local_index] = (
            0.5j * drive[a, b] * (population[a] - population[b])
        )
    used = int(np.count_nonzero(np.abs(source) > 1e-15))
    return source, used


def _matrix_readout_vector(operator, pairs, total_size, offset):
    """Embed Tr(rho O) as a readout row in one coherence block."""
    operator = np.asarray(operator, dtype=complex)
    readout = np.zeros(total_size, dtype=complex)
    for local_index, (a, b) in enumerate(pairs):
        readout[offset + local_index] = operator[b, a]
    return readout


def _resolvent_modes(generator):
    eigenvalues, eigenvectors = np.linalg.eig(generator)
    inverse_eigenvectors = np.linalg.inv(eigenvectors)
    return eigenvalues, eigenvectors, inverse_eigenvectors


def _resolvent_susceptibility_from_modes(
    modes, source, readout, frequencies_hz
):
    eigenvalues, eigenvectors, inverse_eigenvectors = modes
    source_modes = inverse_eigenvectors @ source
    readout_modes = readout @ eigenvectors
    omega = 2.0 * np.pi * np.asarray(frequencies_hz, dtype=float)
    denominator = -1j * omega[:, None] - eigenvalues[None, :]
    denominator[np.abs(denominator) < 1e-12] = 1e-12
    phasor = np.sum(
        readout_modes[None, :] * source_modes[None, :] / denominator,
        axis=1,
    )
    in_phase = 2.0 * np.real(phasor)
    quadrature = 2.0 * np.imag(phasor)
    return np.hypot(in_phase, quadrature), in_phase, quadrature


def _resolvent_susceptibility(generator, source, readout, frequencies_hz):
    return _resolvent_susceptibility_from_modes(
        _resolvent_modes(generator), source, readout, frequencies_hz
    )


def coupled_weak_rf_observable_susceptibilities(result_A, result_B):
    """Solve A-only and B-only weak drives with full SE coherence feedback.

    For the A result, RF-A drives only Alkali A and the readout operator belongs
    only to A; RF-B is zero.  The B result is defined conversely.  The shared
    block generator retains self- and cross-species spin-exchange transfer, so
    the undriven alkali can still feed coherence back through collisions.
    """
    if not result_A["light_shift_available"] or not result_B["light_shift_available"]:
        return None
    pairs_A, er_A, self_A, cross_A, cross_AB = _coherence_exchange_derivatives(
        result_A, result_B
    )
    pairs_B, er_B, self_B, cross_B, cross_BA = _coherence_exchange_derivatives(
        result_B, result_A
    )
    local_A = _local_coherence_generator(result_A, pairs_A, er_A, self_A, cross_A)
    local_B = _local_coherence_generator(result_B, pairs_B, er_B, self_B, cross_B)
    generator = np.block([
        [local_A, float(result_A["R_SE_cross"]) * cross_AB],
        [float(result_B["R_SE_cross"]) * cross_BA, local_B],
    ])
    total_size = generator.shape[0]
    source_A, readout_A, used_A = _drive_and_observable_vectors(
        result_A, pairs_A, total_size, 0
    )
    source_B, readout_B, used_B = _drive_and_observable_vectors(
        result_B, pairs_B, total_size, len(pairs_A)
    )
    response_A = _resolvent_susceptibility(
        generator, source_A, readout_A, result_A["rf_frequencies_hz"]
    )
    response_B = _resolvent_susceptibility(
        generator, source_B, readout_B, result_B["rf_frequencies_hz"]
    )
    return {
        "A": (*response_A, {"used_transitions": used_A, "coupled": True}),
        "B": (*response_B, {"used_transitions": used_B, "coupled": True}),
    }


def coupled_weak_rf_matrix_readouts(
    result_A, result_B, readouts_A, readouts_B
):
    """Evaluate arbitrary species-local readouts with coupled SE dynamics.

    ``readouts_A`` and ``readouts_B`` map stable names to Hermitian ground-state
    operators. RF-A drives only A and is detected only by Probe-A; RF-B is the
    corresponding independent diagonal channel. The common coupled resolvent is
    diagonalized once and reused for every Stokes/rank readout.
    """
    if not result_A["light_shift_available"] or not result_B["light_shift_available"]:
        return None
    pairs_A, er_A, self_A, cross_A, cross_AB = _coherence_exchange_derivatives(
        result_A, result_B
    )
    pairs_B, er_B, self_B, cross_B, cross_BA = _coherence_exchange_derivatives(
        result_B, result_A
    )
    local_A = _local_coherence_generator(result_A, pairs_A, er_A, self_A, cross_A)
    local_B = _local_coherence_generator(result_B, pairs_B, er_B, self_B, cross_B)
    generator = np.block([
        [local_A, float(result_A["R_SE_cross"]) * cross_AB],
        [float(result_B["R_SE_cross"]) * cross_BA, local_B],
    ])
    modes = _resolvent_modes(generator)
    total_size = generator.shape[0]
    source_A, used_A = _drive_vector(result_A, pairs_A, total_size, 0)
    source_B, used_B = _drive_vector(
        result_B, pairs_B, total_size, len(pairs_A)
    )

    def evaluate(readouts, result, pairs, offset, source, used):
        responses = {}
        for name, operator in readouts.items():
            readout = _matrix_readout_vector(operator, pairs, total_size, offset)
            response = _resolvent_susceptibility_from_modes(
                modes, source, readout, result["rf_frequencies_hz"]
            )
            responses[name] = (
                *response,
                {"used_transitions": used, "coupled": True},
            )
        return responses

    return {
        "A": evaluate(readouts_A, result_A, pairs_A, 0, source_A, used_A),
        "B": evaluate(
            readouts_B, result_B, pairs_B, len(pairs_A), source_B, used_B
        ),
    }


def largest_abs_Dm_relaxation_reference(df_pop, target_F=None):
    """Return the relaxation rate at the largest |D_m| in one F manifold.

    By default, the upper ground hyperfine manifold (largest F in ``df_pop``)
    is used. Here D_m = P_m - P_(m-1), as displayed in the Zeeman table, and
    the reference transition is selected by the largest absolute magnitude
    |D_m| within that same manifold.
    Rows without an adjacent lower-m partner are excluded. The total local
    relaxation rate is

        Gamma_m = Gamma_OP,m + Gamma_ER,m + Gamma_SE,m.
    """
    required = {
        "F", "m", "population_difference",
        "Gamma_OP", "Gamma_ER", "Gamma_SE",
    }
    if not required.issubset(df_pop.columns):
        return {
            "available": False,
            "reason": "required columns are missing",
        }

    if target_F is None:
        finite_F_values = df_pop.loc[
            np.isfinite(df_pop["F"].to_numpy(dtype=float)), "F"
        ].to_numpy(dtype=float)
        if finite_F_values.size == 0:
            return {
                "available": False,
                "reason": "no finite hyperfine manifold was found",
            }
        target_F = float(np.max(finite_F_values))
    target_F = float(target_F)

    candidate = df_pop.loc[:, list(required)].copy()
    candidate = candidate.loc[np.isclose(candidate["F"].to_numpy(dtype=float), target_F)]
    finite_D = np.isfinite(candidate["population_difference"].to_numpy(dtype=float))
    candidate = candidate.loc[finite_D]
    if candidate.empty:
        return {
            "available": False,
            "reason": f"no adjacent transition in F={target_F:g} has a finite D_m",
        }

    # Select the adjacent transition with the largest population-difference magnitude.
    row_index = candidate["population_difference"].abs().idxmax()
    row = df_pop.loc[row_index]
    gamma_total = float(row["Gamma_OP"] + row["Gamma_ER"] + row["Gamma_SE"])

    if not np.isfinite(gamma_total) or gamma_total <= 0.0:
        return {
            "available": False,
            "reason": "the selected total relaxation rate is not positive and finite",
            "F": float(row["F"]),
            "m": float(row["m"]),
            "D_m": float(row["population_difference"]),
            "Gamma_m": gamma_total,
        }

    return {
        "available": True,
        "F": float(row["F"]),
        "m": float(row["m"]),
        "D_m": float(row["population_difference"]),
        "Gamma_m": gamma_total,
        "Gamma_OP": float(row["Gamma_OP"]),
        "Gamma_ER": float(row["Gamma_ER"]),
        "Gamma_SE": float(row["Gamma_SE"]),
    }
