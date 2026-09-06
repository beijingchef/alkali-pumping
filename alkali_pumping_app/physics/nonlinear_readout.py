"""Self-consistent physical-pump Stokes readout in the weak-RF limit.

The ordinary probe is a detector only.  A physical pump is different: its
RF-induced Stokes change rotates the same vector/tensor light-shift Hamiltonian
that creates the atomic response.  This module closes that feedback loop while
retaining the application's first-order RF density-matrix approximation.

The longitudinal coordinate is normalized to the selected path length.  The
weak full-cell response supplies the inhomogeneous propagation term. Fractional
intensity and normalized Stokes changes drive both the dispersive light-shift
commutator and the dissipative optical-pumping Liouvillian. Rank-1 and rank-2
optical readouts are propagated either separately as counterfactual diagnostics
or together as the physical coherent total. The latter is not the sum of the
two nonlinear counterfactual solutions.
"""

import numpy as np
from scipy.linalg import expm

from .optical_pumping import (
    build_optical_L,
    optical_pumping_stokes_liouvillian_sources,
    optical_rate_scale_from_intensity,
)
from .polarization import transverse_basis_for_k
from .rf_response import (
    _spin_operator,
    build_single_species_liouvillian_response_context,
    weak_liouvillian_source_matrix_readouts,
)


_SIGNALS = ("transmission", "s1", "s2", "s3", "rotation", "ellipticity")
_FEEDBACK_COLUMNS = {"transmission": 0, "s1": 1, "s2": 2, "s3": 3}


def _axis_name(vector):
    values = np.abs(np.asarray(vector, dtype=complex))
    index = int(np.argmax(values))
    if not np.isclose(values[index], 1.0) or np.count_nonzero(values > 1e-12) != 1:
        raise ValueError("The pump transverse basis must follow laboratory axes.")
    return ("x", "y", "z")[index]


def _rank_coefficients_from_circular_reference(result, beam, common):
    """Return per-manifold vector and linear-tensor shift coefficients.

    A sigma+ calculation with quantization along k supplies the intrinsic
    vector coefficient directly.  Its rank-2 coefficient is -1/2 of the
    coefficient for a linearly polarized field, because the symmetric
    polarization dyadic of circular light is half the transverse projector.
    Returned coefficients are angular frequencies in rad/s.
    """
    atom = result["atom"]
    states = result["ground_states"]
    line = beam["line"]
    rate_scale = optical_rate_scale_from_intensity(
        atom=atom,
        line=line,
        intensity_uW_cm2=float(beam["intensity"]),
        n2_pressure_torr=common["n2_pressure_torr"],
        temperature_C=common["temperature_C"],
        n2_width_MHz_per_torr=result["n2_coeffs"][line]["width"],
    )
    _generator, info = build_optical_L(
        atom=atom,
        line=line,
        ground_states=states,
        detuning_MHz=float(beam["detuning"]),
        pump_rate_s=rate_scale,
        selected_transition=None,
        k_axis=beam["k_axis"],
        pol="sigma+",
        q_axis=beam["k_axis"],
        n2_pressure_torr=common["n2_pressure_torr"],
        temperature_C=common["temperature_C"],
        n2_width_MHz_per_torr=result["n2_coeffs"][line]["width"],
        n2_shift_MHz_per_torr=result["n2_coeffs"][line]["shift"],
        normalize_to_selected_total=False,
    )
    shift_hz = np.sum(info["light_shift_ge_angular"], axis=1) / (2.0 * np.pi)
    coefficients = {}
    for F in sorted({float(state["F"]) for state in states}):
        indices = np.array(
            [
                index
                for index, state in enumerate(states)
                if np.isclose(float(state["F"]), F)
            ],
            dtype=int,
        )
        m = np.array([float(states[index]["m"]) for index in indices])
        tensor_basis = 3.0 * m**2 - F * (F + 1.0)
        design = np.column_stack([np.ones_like(m), m, tensor_basis])
        scalar, vector, circular_tensor = np.linalg.lstsq(
            design, shift_hz[indices], rcond=None
        )[0]
        coefficients[F] = {
            "scalar_rad_s": float(2.0 * np.pi * scalar),
            "vector_rad_s": float(2.0 * np.pi * vector),
            "linear_tensor_rad_s": float(-4.0 * np.pi * circular_tensor),
        }
    return coefficients


def pump_stokes_light_shift_drives(result, beam, common):
    """Return d(H/hbar) for fractional intensity and normalized Stokes."""
    states = result["ground_states"]
    size = len(states)
    e1, e2 = transverse_basis_for_k(beam["k_axis"])
    axis_1 = _axis_name(e1)
    axis_2 = _axis_name(e2)
    coefficients = _rank_coefficients_from_circular_reference(
        result, beam, common
    )
    drives = {
        name: np.zeros((size, size), dtype=complex)
        for name in ("s1", "s2", "s3")
    }
    for F, values in coefficients.items():
        first = _spin_operator(states, result["q_axis"], axis_1, F)
        second = _spin_operator(states, result["q_axis"], axis_2, F)
        propagation = _spin_operator(
            states, result["q_axis"], beam["k_axis"], F
        )
        tensor = values["linear_tensor_rad_s"]
        drives["s1"] += 1.5 * tensor * (first @ first - second @ second)
        drives["s2"] += 1.5 * tensor * (first @ second + second @ first)
        drives["s3"] += values["vector_rad_s"] * propagation
    drives = {
        name: 0.5 * (operator + np.conjugate(operator.T))
        for name, operator in drives.items()
    }
    scalar = np.zeros((size, size), dtype=complex)
    for F, values in coefficients.items():
        for index, state in enumerate(states):
            if np.isclose(float(state["F"]), F):
                scalar[index, index] = values["scalar_rad_s"]
    input_stokes = np.asarray(result["probe_info"]["input_stokes"], dtype=float)
    normalized = input_stokes[1:] / input_stokes[0]
    drives["transmission"] = (
        scalar
        + normalized[0] * drives["s1"]
        + normalized[1] * drives["s2"]
        + normalized[2] * drives["s3"]
    )
    return {
        name: drives[name]
        for name in ("transmission", "s1", "s2", "s3")
    }, coefficients


def _complex_response(response):
    return (
        np.asarray(response["in_phase"], dtype=float)
        + 1j * np.asarray(response["quadrature"], dtype=float)
    )


def propagate_stokes_feedback(base_by_signal, feedback_by_signal_and_stokes):
    """Integrate a uniform full-cell response with distributed pump feedback.

    The independent coordinate is t=z/L.  If ``b`` is the ordinary weak
    full-cell response and ``K`` is the full-cell response to a spatially
    uniform normalized-Stokes perturbation, propagation obeys

        dy/dt = b + K @ (y_T, y_s1, y_s2, y_s3),  y(0)=0.

    The matrix exponential integrates this constant-coefficient equation
    without a slice-count or fixed-point convergence parameter.
    """
    sample = next(iter(base_by_signal.values()))
    frequency_count = np.asarray(sample).size
    base = np.column_stack([base_by_signal[name] for name in _SIGNALS])
    feedback = np.empty(
        (frequency_count, len(_SIGNALS), len(_FEEDBACK_COLUMNS)),
        dtype=complex,
    )
    for signal_index, signal in enumerate(_SIGNALS):
        for coordinate, column in _FEEDBACK_COLUMNS.items():
            values = feedback_by_signal_and_stokes[signal].get(coordinate)
            feedback[:, signal_index, column] = (
                np.zeros(frequency_count, dtype=complex)
                if values is None
                else values
            )

    output = np.empty_like(base)
    spectral_radius = np.empty(frequency_count, dtype=float)
    feedback_indices = [
        _SIGNALS.index(name) for name in _FEEDBACK_COLUMNS
    ]
    for frequency_index in range(frequency_count):
        propagation = np.zeros(
            (len(_SIGNALS) + 1, len(_SIGNALS) + 1), dtype=complex
        )
        for column, signal_index in enumerate(feedback_indices):
            propagation[:-1, signal_index] = feedback[frequency_index, :, column]
        propagation[:-1, -1] = base[frequency_index]
        initial = np.zeros(len(_SIGNALS) + 1, dtype=complex)
        initial[-1] = 1.0
        output[frequency_index] = (expm(propagation) @ initial)[:-1]
        loop = feedback[frequency_index, feedback_indices, :]
        spectral_radius[frequency_index] = float(
            np.max(np.abs(np.linalg.eigvals(loop)))
        )
    return {
        signal: output[:, index]
        for index, signal in enumerate(_SIGNALS)
    }, spectral_radius


def _unavailable_response(response):
    unavailable = {}
    for component, signals in response.items():
        unavailable[component] = {}
        for signal, values in signals.items():
            nan = np.full_like(np.asarray(values["amplitude"], dtype=float), np.nan)
            unavailable[component][signal] = {
                "amplitude": nan,
                "in_phase": nan.copy(),
                "quadrature": nan.copy(),
                "info": {**values["info"], "nonlinear_available": False},
            }
    return unavailable


def apply_physical_pump_readout(result, common):
    """Replace one weak linked-probe response by physical-pump propagation."""
    probe = result.get("probe", {})
    if probe.get("mode") != "nonlinear":
        return result

    pump_name = probe.get("pump_name")
    diagnostic = next(
        (
            (beam, info)
            for beam, info in result.get("diagnostics", ())
            if beam.get("name") == pump_name
        ),
        None,
    )
    result["probe_weak_response"] = result["probe_response"]
    if diagnostic is None or float(probe.get("pump_intensity_uW_cm2", 0.0)) <= 0.0:
        result["probe_response"] = _unavailable_response(result["probe_response"])
        result["probe_info"] = {
            **result["probe_info"],
            "mode": "nonlinear physical pump",
            "pump_name": pump_name,
            "nonlinear_available": False,
            "nonlinear_reason": "The selected physical pump has zero intensity.",
        }
        return result

    beam, _info = diagnostic
    drives, rank_coefficients = pump_stokes_light_shift_drives(
        result, beam, common
    )
    input_stokes = np.asarray(result["probe_info"]["input_stokes"], dtype=float)
    normalized_stokes = input_stokes[1:] / input_stokes[0]
    rho_0 = np.diag(np.asarray(result["population"], dtype=complex))
    dissipative_sources = optical_pumping_stokes_liouvillian_sources(
        atom=result["atom"],
        line=beam["line"],
        ground_states=result["ground_states"],
        detuning_MHz=float(beam["detuning"]),
        intensity_uW_cm2=float(beam["intensity"]),
        k_axis=beam["k_axis"],
        q_axis=result["q_axis"],
        input_stokes=normalized_stokes,
        density_matrix=rho_0,
        n2_pressure_torr=common["n2_pressure_torr"],
        temperature_C=common["temperature_C"],
        n2_width_MHz_per_torr=result["n2_coeffs"][beam["line"]]["width"],
        n2_shift_MHz_per_torr=result["n2_coeffs"][beam["line"]]["shift"],
    )
    liouvillian_sources = {}
    source_diagnostics = {}
    for coordinate, drive in drives.items():
        dispersive = -1j * (drive @ rho_0 - rho_0 @ drive)
        dissipative = dissipative_sources[coordinate]
        source = dispersive + dissipative
        source = 0.5 * (source + np.conjugate(source.T))
        source -= np.trace(source) * np.eye(len(source), dtype=complex) / len(source)
        liouvillian_sources[coordinate] = source
        source_diagnostics[coordinate] = {
            "dispersive_frobenius_norm_s_inv": float(np.linalg.norm(dispersive)),
            "dissipative_frobenius_norm_s_inv": float(np.linalg.norm(dissipative)),
            "total_frobenius_norm_s_inv": float(np.linalg.norm(source)),
        }

    flat_readouts = {
        f"{component}:{signal}": operator
        for component, operators in result["probe_readout_operators"].items()
        for signal, operator in operators.items()
    }
    response_context = build_single_species_liouvillian_response_context(result)
    feedback_responses = {}
    used_transitions = 0
    population_source_coordinates = []
    for coordinate, source in liouvillian_sources.items():
        responses = weak_liouvillian_source_matrix_readouts(
            result,
            source,
            flat_readouts,
            context=response_context,
        )
        feedback_responses[coordinate] = responses
        for response in responses.values():
            used_transitions = max(
                used_transitions, int(response[3]["used_transitions"])
            )
            if response[3]["population_source"]:
                population_source_coordinates.append(coordinate)

    nonlinear_response = {}
    feedback_diagnostics = {}
    for component, operators in result["probe_readout_operators"].items():
        base_by_signal = {
            signal: _complex_response(result["probe_response"][component][signal])
            for signal in _SIGNALS
        }
        feedback_by_signal_and_stokes = {
            signal: {} for signal in _SIGNALS
        }
        for coordinate in _FEEDBACK_COLUMNS:
            responses = feedback_responses[coordinate]
            for signal in _SIGNALS:
                response = responses[f"{component}:{signal}"]
                feedback_by_signal_and_stokes[signal][coordinate] = (
                    np.asarray(response[1], dtype=float)
                    + 1j * np.asarray(response[2], dtype=float)
                )
        propagated, spectral_radius = propagate_stokes_feedback(
            base_by_signal, feedback_by_signal_and_stokes
        )
        nonlinear_response[component] = {}
        for signal, complex_values in propagated.items():
            base_info = result["probe_response"][component][signal]["info"]
            in_phase = np.real(complex_values)
            quadrature = np.imag(complex_values)
            nonlinear_response[component][signal] = {
                "amplitude": np.hypot(in_phase, quadrature),
                "in_phase": in_phase,
                "quadrature": quadrature,
                "info": {
                    **base_info,
                    "nonlinear": True,
                    "physical_pump": pump_name,
                    "used_feedback_transitions": used_transitions,
                },
            }
        feedback_diagnostics[component] = {
            "spectral_radius": spectral_radius,
            "max_spectral_radius": float(np.max(spectral_radius)),
        }

    result["probe_response"] = nonlinear_response
    result["probe_info"] = {
        **result["probe_info"],
        "mode": "nonlinear physical pump",
        "pump_name": pump_name,
        "pump_intensity_uW_cm2": float(probe["pump_intensity_uW_cm2"]),
        "nonlinear_available": True,
        "feedback_coordinates": tuple(_FEEDBACK_COLUMNS),
        "feedback_contributions": (
            "dispersive light shift",
            "dissipative optical pumping",
        ),
        "rank_shift_coefficients": rank_coefficients,
        "liouvillian_source_diagnostics": source_diagnostics,
        "population_source_coordinates": tuple(
            dict.fromkeys(population_source_coordinates)
        ),
        "feedback_diagnostics": feedback_diagnostics,
        "model": (
            "distributed weak-RF dispersive+dissipative pump Liouvillian readout"
        ),
    }
    return result
