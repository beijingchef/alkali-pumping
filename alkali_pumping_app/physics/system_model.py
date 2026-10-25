"""High-level construction of single- and dual-alkali model results."""

import numpy as np
import pandas as pd
from copy import deepcopy

from .angular_momentum import build_ground_states
from .constants import (
    ATOMS,
    alkali_vapor_density_cm3,
    cross_spin_exchange_rate_info,
    ground_zeeman_shifts_hz,
    self_spin_exchange_rate_info,
    upper_larmor_frequency_from_field_nT,
)
from .multi_species import (
    coupled_population_jacobian,
    cross_spin_exchange_adjacent_coherence_self_relaxation_rates,
    steady_state_two_species,
)
from .nonlinear_readout import apply_physical_pump_readout
from .observables import (
    add_adjacent_optical_relaxation_columns,
    add_nu_m_column,
    add_population_difference_column,
    add_total_relaxation_columns,
    optical_Lambda_fractional_rates,
    population_by_F,
    steady_state_from_L,
    total_G_OP_by_ground_state,
)
from .optical_pumping import (
    build_optical_L,
    decompose_nu_LS_components,
    optical_rate_scale_from_intensity,
    total_nu_LS_from_diagnostics,
)
from .optical_readout import manifold_probe_readout_operators, weak_probe_readout_operators
from .rf_response import (
    atomic_rf_readout_operators,
    build_rf_response_contexts,
    hyperfine_manifolds,
    largest_abs_Dm_relaxation_reference,
    magnetic_rf_operator,
    weak_rf_matrix_readouts,
)
from .spin_exchange import (
    build_ER_matrix,
    build_spin_exchange_matrix,
    er_adjacent_coherence_self_relaxation_rates,
    er_population_fractional_relaxation_rates,
    spin_exchange_adjacent_coherence_self_relaxation_rates,
    spin_exchange_population_fractional_relaxation_rates,
    spin_exchange_population_jacobian,
    steady_state_with_spin_exchange,
)
from .spectroscopy import line_center_frequency_MHz
from .uniform_relaxation import uniform_relaxation_diagnostics, uniform_relaxation_generator


def _flatten_probe_operators(operators):
    return {
        f"{component}:{signal}": operator
        for component, signals in operators.items()
        for signal, operator in signals.items()
    }


def _nested_probe_responses(flat_responses):
    nested = {}
    for name, response in flat_responses.items():
        component, signal = name.split(":", 1)
        amplitude, in_phase, quadrature, info = response
        nested.setdefault(component, {})[signal] = {
            "amplitude": amplitude,
            "in_phase": in_phase,
            "quadrature": quadrature,
            "info": info,
        }
    return nested


def resolve_alkali_densities(atom_A_name, atom_B_name, temperature_C, mode, ratio):
    """Resolve A/B number densities for the selected mixture convention."""
    saturated_A = alkali_vapor_density_cm3(atom_A_name, temperature_C)
    if atom_B_name in (None, "None") or atom_B_name == atom_A_name:
        return saturated_A, 0.0
    if mode == "Relative concentration":
        liquid_ratio = max(0.0, float(ratio))
        liquid_fraction_A = 1.0 / (1.0 + liquid_ratio)
        liquid_fraction_B = 1.0 - liquid_fraction_A
        saturated_B = alkali_vapor_density_cm3(atom_B_name, temperature_C)
        density_A = liquid_fraction_A * saturated_A
        density_B = liquid_fraction_B * saturated_B
    else:
        density_A = saturated_A
        density_B = alkali_vapor_density_cm3(atom_B_name, temperature_C)
    return float(density_A), float(density_B)


def _nearest_line_and_detuning(atom, absolute_frequency_MHz):
    lines = ("D1", "D2")
    line = min(
        lines,
        key=lambda candidate: abs(
            float(absolute_frequency_MHz)
            - line_center_frequency_MHz(atom, candidate)
        ),
    )
    return line, float(absolute_frequency_MHz) - line_center_frequency_MHz(atom, line)


def build_species_linear_model(config, all_beams, common):
    """Build optical and ER generators for one species before spin exchange."""
    atom_name = config["atom_name"]
    atom = ATOMS[atom_name]
    ground_states = build_ground_states(atom)
    size = len(ground_states)
    optical = np.zeros((size, size), dtype=float)
    diagnostics = []
    species_beams = []
    light_shift_diagnostics = []

    for source_beam in all_beams:
        beam = dict(source_beam)
        line, detuning = _nearest_line_and_detuning(
            atom, beam["absolute_frequency_MHz"]
        )
        beam["line"] = line
        beam["detuning"] = detuning
        beam["q_axis"] = config["q_axis"]
        if beam["target_label"] != config["label"]:
            beam["selected_transition"] = None
            beam["transition_label"] = "off-target"
        species_beams.append(beam)

        if float(beam["intensity"]) <= 0.0:
            continue
        rate_scale = optical_rate_scale_from_intensity(
            atom=atom,
            line=line,
            intensity_uW_cm2=beam["intensity"],
            n2_pressure_torr=common["n2_pressure_torr"],
            temperature_C=common["temperature_C"],
            n2_width_MHz_per_torr=config["n2_coeffs"][line]["width"],
        )
        generator, info = build_optical_L(
            atom=atom,
            line=line,
            ground_states=ground_states,
            detuning_MHz=detuning,
            pump_rate_s=rate_scale,
            selected_transition=beam.get("selected_transition"),
            k_axis=beam["k_axis"],
            pol=beam["pol"],
            q_axis=config["q_axis"],
            n2_pressure_torr=common["n2_pressure_torr"],
            temperature_C=common["temperature_C"],
            n2_width_MHz_per_torr=config["n2_coeffs"][line]["width"],
            n2_shift_MHz_per_torr=config["n2_coeffs"][line]["shift"],
            normalize_to_selected_total=False,
        )
        optical += generator
        diagnostics.append((beam, info))
        # A pump remains in the optical-pumping generator for any physically
        # possible off-resonant absorption, but its AC-Stark shift is reported
        # only for the alkali species that the pump targets.  In particular,
        # the other species' beam geometry must not make this species' local
        # light-shift Hamiltonian appear non-diagonal and blank its table.
        if beam["target_label"] == config["label"]:
            light_shift_diagnostics.append((beam, info))

    er_map = build_ER_matrix(atom, ground_states)
    # Omitted rates retain the previous low-level API behavior. New-session
    # defaults are supplied explicitly by each page's versioned settings.
    R_u = float(config.get("R_u", 0.0))
    uniform = uniform_relaxation_generator(size, R_u)
    relaxation = float(config["R_ER"]) * (er_map - np.eye(size)) + uniform
    linear = optical + relaxation
    self_rate_info = self_spin_exchange_rate_info(
        atom_name,
        common["temperature_C"],
        config["density_cm3"],
    )
    return {
        **config,
        "atom": atom,
        "ground_states": ground_states,
        "size": size,
        "beams": species_beams,
        "diagnostics": diagnostics,
        "light_shift_diagnostics": light_shift_diagnostics,
        "L_optical": optical,
        "M_ER": er_map,
        "R_u": R_u,
        "L_uniform": uniform,
        "L_relaxation": relaxation,
        "L_linear": linear,
        "self_rate_info": self_rate_info,
        "R_SE_self": self_rate_info["rate_s"],
    }


def _finalize_species(
    model,
    population,
    self_map,
    cross_map,
    partner_electron,
    R_SE_cross,
    jacobian_local,
    common,
):
    atom = model["atom"]
    states = model["ground_states"]
    R_self = model["R_SE_self"]
    R_ER = float(model["R_ER"])

    nu_LS, light_shift_available = total_nu_LS_from_diagnostics(
        states,
        [beam for beam, _info in model["light_shift_diagnostics"]],
        model["light_shift_diagnostics"],
    )
    components = decompose_nu_LS_components(states, nu_LS)
    static_field_aligned = model["q_axis"] == common["static_field_axis"]
    field_parallel_nT = common["static_field_nT"] if static_field_aligned else 0.0
    bias_hz = upper_larmor_frequency_from_field_nT(
        model["atom_name"], field_parallel_nT
    )
    nu_B, bias_info = ground_zeeman_shifts_hz(
        model["atom_name"], atom, states, bias_hz
    )
    G_OP = total_G_OP_by_ground_state(states, model["diagnostics"])
    Lambda = optical_Lambda_fractional_rates(model["L_optical"], population, G_OP)
    G_ER = er_population_fractional_relaxation_rates(model["M_ER"], population, R_ER)
    Gamma_ER = er_adjacent_coherence_self_relaxation_rates(atom, states, R_ER)
    G_u, Gamma_u = uniform_relaxation_diagnostics(states, population, model["R_u"])
    G_SE_self = spin_exchange_population_fractional_relaxation_rates(
        self_map, population, R_self
    )
    G_SE_cross = spin_exchange_population_fractional_relaxation_rates(
        cross_map, population, R_SE_cross
    )
    Gamma_SE_self = spin_exchange_adjacent_coherence_self_relaxation_rates(
        atom,
        states,
        population,
        model["electron_marginal"],
        R_self,
    )
    Gamma_SE_cross = cross_spin_exchange_adjacent_coherence_self_relaxation_rates(
        atom,
        states,
        partner_electron,
        R_SE_cross,
    )
    G_SE = G_SE_self + G_SE_cross
    Gamma_SE = Gamma_SE_self + Gamma_SE_cross

    df = pd.DataFrame({
        "F": [state["F"] for state in states],
        "m": [state["m"] for state in states],
        "population": population,
        "nu_VS": components["vector"],
        "nu_TS": components["tensor"],
        "nu_LS": nu_LS,
        "nu_B": nu_B,
        "Lambda": Lambda,
        "G_OP": G_OP,
        "G_ER": G_ER,
        "Gamma_ER": Gamma_ER,
        "G_u": G_u,
        "Gamma_u": Gamma_u,
        "G_SE_self": G_SE_self,
        "Gamma_SE_self": Gamma_SE_self,
        "G_SE_cross": G_SE_cross,
        "Gamma_SE_cross": Gamma_SE_cross,
        "G_SE": G_SE,
        "Gamma_SE": Gamma_SE,
    })
    df = add_population_difference_column(df)
    df = add_nu_m_column(df)
    df = add_adjacent_optical_relaxation_columns(df)
    df = add_total_relaxation_columns(df)

    probe = model.get("probe", {
        "line": "D1",
        "detuning_MHz": 1000.0,
        "k_axis": "x",
        "azimuth_deg": 45.0,
        "ellipticity_deg": 0.0,
        "path_length_cm": 1.0,
        "include_scalar": False,
        "include_orientation": True,
        "include_alignment": True,
    })
    probe_operators, probe_info = weak_probe_readout_operators(
        atom=atom,
        ground_states=states,
        q_axis=model["q_axis"],
        line=probe["line"],
        detuning_MHz=probe["detuning_MHz"],
        k_axis=probe["k_axis"],
        azimuth_deg=probe["azimuth_deg"],
        ellipticity_deg=probe["ellipticity_deg"],
        path_length_cm=probe["path_length_cm"],
        density_cm3=model["density_cm3"],
        n2_pressure_torr=common["n2_pressure_torr"],
        temperature_C=common["temperature_C"],
        n2_width_MHz_per_torr=model["n2_coeffs"][probe["line"]]["width"],
        n2_shift_MHz_per_torr=model["n2_coeffs"][probe["line"]]["shift"],
        include_scalar=probe["include_scalar"],
        include_orientation=probe["include_orientation"],
        include_alignment=probe["include_alignment"],
    )

    upper_F = max(float(state["F"]) for state in states)
    frequencies = model["rf_frequencies_hz"]

    df_F = population_by_F(df)
    df["hyperfine_population"] = np.nan
    for _, row in df_F.iterrows():
        mask = (
            np.isclose(df["F"].to_numpy(dtype=float), float(row["F"]))
            & np.isclose(df["m"].to_numpy(dtype=float), 0.0)
        )
        df.loc[mask, "hyperfine_population"] = float(row["population"])

    model.update({
        "population": population,
        "self_map": self_map,
        "cross_map": cross_map,
        "R_SE_cross": float(R_SE_cross),
        "J_population": jacobian_local,
        "electron_marginal": model["electron_marginal"],
        "partner_electron_marginal": np.asarray(partner_electron, dtype=float),
        "nu_LS": nu_LS,
        "light_shift_available": light_shift_available,
        "bias_larmor_hz": bias_hz,
        "static_field_aligned": static_field_aligned,
        "field_parallel_nT": field_parallel_nT,
        "bias_info": bias_info,
        "df_pop": df,
        "df_F": df_F,
        "rf_upper_F": upper_F,
        "rf_frequencies_hz": frequencies,
        "probe_readout_operators": probe_operators,
        "probe_info": probe_info,
        "rf_relaxation_reference": largest_abs_Dm_relaxation_reference(
            df, target_F="all"
        ),
    })
    return model


def compute_alkali_system(species_A_config, species_B_config, all_beams, common):
    """Compute one active alkali or a coupled A/B system."""
    model_A = build_species_linear_model(species_A_config, all_beams, common)
    model_B = (
        build_species_linear_model(species_B_config, all_beams, common)
        if species_B_config is not None
        else None
    )

    if model_B is None:
        R_self = model_A["R_SE_self"]
        if R_self > 0.0:
            p_A, solve = steady_state_with_spin_exchange(
                model_A["L_linear"],
                model_A["atom"],
                model_A["ground_states"],
                R_self,
            )
            self_map = solve["M_SE"]
            electron_A = solve["electron_marginal"]
        else:
            p_A = steady_state_from_L(model_A["L_linear"])
            self_map, electron_A = build_spin_exchange_matrix(
                model_A["atom"], model_A["ground_states"], p_A
            )
            solve = {
                "iterations": 0,
                "converged": True,
                "residual": float(np.max(np.abs(model_A["L_linear"] @ p_A))),
                "mirror_symmetry_enforced": False,
            }
        model_A["electron_marginal"] = electron_A
        cross_map = np.eye(model_A["size"])
        J_A = model_A["L_linear"] + spin_exchange_population_jacobian(
            model_A["atom"], model_A["ground_states"], p_A, R_self
        )
        result_A = _finalize_species(
            model_A,
            p_A,
            self_map,
            cross_map,
            np.array([0.5, 0.5]),
            0.0,
            J_A,
            common,
        )
        _evaluate_system_responses(result_A, None, common, J_A)
        return {"A": result_A, "B": None, "J_coupled": J_A, "solve": solve}

    cross_A_info = cross_spin_exchange_rate_info(
        model_A["atom_name"],
        model_B["atom_name"],
        common["temperature_C"],
        model_B["density_cm3"],
    )
    cross_B_info = cross_spin_exchange_rate_info(
        model_B["atom_name"],
        model_A["atom_name"],
        common["temperature_C"],
        model_A["density_cm3"],
    )
    R_A_cross = cross_A_info["rate_s"]
    R_B_cross = cross_B_info["rate_s"]
    p_A, p_B, solve = steady_state_two_species(
        model_A["L_linear"],
        model_A["atom"],
        model_A["ground_states"],
        model_A["R_SE_self"],
        R_A_cross,
        model_B["L_linear"],
        model_B["atom"],
        model_B["ground_states"],
        model_B["R_SE_self"],
        R_B_cross,
    )
    model_A["electron_marginal"] = solve["A"]["electron_marginal"]
    model_B["electron_marginal"] = solve["B"]["electron_marginal"]
    full_J, blocks = coupled_population_jacobian(
        model_A["L_linear"],
        model_A["atom"],
        model_A["ground_states"],
        p_A,
        model_A["R_SE_self"],
        R_A_cross,
        model_B["L_linear"],
        model_B["atom"],
        model_B["ground_states"],
        p_B,
        model_B["R_SE_self"],
        R_B_cross,
    )
    result_A = _finalize_species(
        model_A,
        p_A,
        solve["A"]["self_map"],
        solve["A"]["cross_map"],
        solve["B"]["electron_marginal"],
        R_A_cross,
        blocks["AA"],
        common,
    )
    result_B = _finalize_species(
        model_B,
        p_B,
        solve["B"]["self_map"],
        solve["B"]["cross_map"],
        solve["A"]["electron_marginal"],
        R_B_cross,
        blocks["BB"],
        common,
    )
    result_A["cross_rate_info"] = cross_A_info
    result_B["cross_rate_info"] = cross_B_info
    _evaluate_system_responses(result_A, result_B, common, full_J)
    return {
        "A": result_A,
        "B": result_B,
        "J_coupled": full_J,
        "J_blocks": blocks,
        "solve": solve,
    }


def _attach_local_response(result, response_config, common):
    """Attach RF and probe calculations to an already-solved species result."""
    result.update(response_config)
    states = result["ground_states"]
    probe = response_config["probe"]
    operators, probe_info = weak_probe_readout_operators(
        atom=result["atom"], ground_states=states, q_axis=result["q_axis"],
        line=probe["line"], detuning_MHz=probe["detuning_MHz"],
        k_axis=probe["k_axis"], azimuth_deg=probe["azimuth_deg"],
        ellipticity_deg=probe["ellipticity_deg"], path_length_cm=probe["path_length_cm"],
        density_cm3=result["density_cm3"], n2_pressure_torr=common["n2_pressure_torr"],
        temperature_C=common["temperature_C"],
        n2_width_MHz_per_torr=result["n2_coeffs"][probe["line"]]["width"],
        n2_shift_MHz_per_torr=result["n2_coeffs"][probe["line"]]["shift"],
        include_scalar=probe["include_scalar"],
        include_orientation=probe["include_orientation"],
        include_alignment=probe["include_alignment"],
    )
    result.update({
        "probe_readout_operators": operators, "probe_info": probe_info,
    })


def compute_alkali_responses(equilibrium, response_A, response_B, common):
    """Reuse a population equilibrium while recalculating RF/probe responses."""
    system = deepcopy(equilibrium)
    result_A, result_B = system["A"], system["B"]
    _attach_local_response(result_A, response_A, common)
    if result_B is not None:
        _attach_local_response(result_B, response_B, common)
    _evaluate_system_responses(result_A, result_B, common, system["J_coupled"])
    return system


def _response_values(response):
    amplitude, in_phase, quadrature, info = response
    return {"amplitude": amplitude, "in_phase": in_phase,
            "quadrature": quadrature, "info": info}


def _evaluate_system_responses(result_A, result_B, common, population_generator):
    """Use identical dynamics for every atomic moment and optical readout."""
    results = {"A": result_A}
    if result_B is not None:
        results["B"] = result_B
    available = all(
        result["light_shift_available"]
        and (result["static_field_aligned"] or abs(common["static_field_nT"]) == 0.0)
        for result in results.values()
    )
    contexts = (build_rf_response_contexts(result_A, result_B, population_generator)
                if available else {})
    for label, result in results.items():
        result.pop("probe_weak_response", None)
        result.pop("probe_weak_manifold_response", None)
        result["rf_manifolds"] = hyperfine_manifolds(result)
        result["rf_available"] = available
        _, drive_info = magnetic_rf_operator(result)
        result["rf_drive_info"] = drive_info
        operators = result["probe_readout_operators"]
        manifold_operators = manifold_probe_readout_operators(result["ground_states"], operators)
        result["probe_manifold_operators"] = manifold_operators
        readouts = {f"atomic:{name}": op for name, op in atomic_rf_readout_operators(result).items()}
        readouts.update({f"probe:{name}": op for name, op in _flatten_probe_operators(operators).items()})
        for manifold, parts in manifold_operators.items():
            readouts.update({f"{manifold}:{name}": op for name, op in _flatten_probe_operators(parts).items()})
        if available:
            responses = weak_rf_matrix_readouts(result, readouts, context=contexts[label])
        else:
            nan = np.full_like(result["rf_frequencies_hz"], np.nan, dtype=float)
            responses = {name: (nan.copy(), nan.copy(), nan.copy(), {
                "used_transitions": 0, "available": False, **drive_info,
            }) for name in readouts}
        result["rf_response"] = {
            name: _response_values(responses[f"atomic:{name}"])
            for name in result["rf_manifolds"]
        }
        # These scalar aliases refer explicitly to Fa for callers of the old
        # result structure. The application displays the manifold dictionary.
        for key in ("amplitude", "in_phase", "quadrature", "info"):
            result[f"rf_{key}"] = result["rf_response"]["Fa"][key]
        result["probe_response"] = _nested_probe_responses({
            name.removeprefix("probe:"): response
            for name, response in responses.items() if name.startswith("probe:")
        })
        result["probe_manifold_response"] = {
            manifold: _nested_probe_responses({
                name.removeprefix(f"{manifold}:"): response
                for name, response in responses.items() if name.startswith(f"{manifold}:")
            }) for manifold in result["rf_manifolds"]
        }
        result["probe_info"].update({"rf_available": available, **drive_info})
        if not available and result.get("probe", {}).get("mode") == "physical":
            result["probe_info"].update({
                "mode": "uniform-state physical pump", "physical_pump_available": False,
                "pump_name": result["probe"].get("pump_name"),
                "physical_pump_reason": "The coupled atomic RF response is unavailable.",
            })
        if available:
            apply_physical_pump_readout(result, common, context=contexts[label])
