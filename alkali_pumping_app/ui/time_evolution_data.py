"""Plot and export data for population dynamics; units match visible curves."""

import numpy as np
import pandas as pd


def zeeman_plot_frame(time_s, species, manifold):
    """Individual Zeeman states only; manifold totals remain in complete CSVs."""
    states = species["ground_states"]
    values = sorted({state["F"] for state in states})
    F = values[-1] if manifold == "Fa" else values[0]
    return pd.DataFrame({
        "Time (s)": np.tile(time_s, sum(state["F"] == F for state in states)),
        "Population": np.concatenate([species["population"][:, i]
                                       for i, state in enumerate(states) if state["F"] == F]),
        "Series": np.concatenate([np.repeat(f"m = {state['m']:g}", len(time_s))
                                   for state in states if state["F"] == F]),
    })


def series_frame(time_s, series, y_label):
    return pd.DataFrame({"Time (s)": np.tile(time_s, len(series)),
                         y_label: np.concatenate(list(series.values())),
                         "Series": np.repeat(list(series), len(time_s))})


def plot_export_dataframe(frame, y_label, *, species, label, title, value_unit, t1_s, t2_s):
    """Export exactly the displayed curves, with one numeric column per curve.

    The frame already contains the plotted units (including mrad conversion).
    Keep its series order and include event metadata for standalone graph files.
    """
    series = list(frame["Series"].drop_duplicates())
    values = frame.pivot(index="Time (s)", columns="Series", values=y_label)
    values = values.reindex(columns=series).sort_index().rename_axis(columns=None)
    values = values.rename_axis("time_s").reset_index()
    values["pumps_on"] = values["time_s"] < t1_s
    values["pump_off_time_s"] = t1_s
    values["total_time_s"] = t2_s
    values["alkali"] = label
    values["atom"] = species["atom_name"]
    values["quantization_axis"] = species["q_axis"]
    values["zeeman_quantization_axis"] = species.get("zeeman_axis", species["q_axis"])
    values["moment_quantization_axis"] = species.get("moment_axis", species["q_axis"])
    moment_axis = species.get("moment_axis", species["q_axis"])
    values["alignment_u_axis"], values["alignment_v_axis"] = {
        "x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y"),
    }[moment_axis]
    values["observable"] = title
    values["value_units"] = value_unit
    return values


def evolution_export_dataframe(evolution, label, probe_signals, config, probe, common, beams):
    species = evolution[label]
    moment_axis = species.get("moment_axis", species["q_axis"])
    times = evolution["time_s"]
    columns = {
        "time_s": times, "pumps_on": times < evolution["t1_s"],
        "pump_off_time_s": evolution["t1_s"], "total_time_s": evolution["t2_s"],
        "atom": species["atom_name"], "quantization_axis": species["q_axis"],
        "zeeman_quantization_axis": species.get("zeeman_axis", species["q_axis"]),
        "moment_quantization_axis": moment_axis,
        "alignment_u_axis": {"x": "y", "y": "z", "z": "x"}[moment_axis],
        "alignment_v_axis": {"x": "z", "y": "x", "z": "y"}[moment_axis],
        "density_cm-3": species["density_cm3"], "R_ER_s-1": config["R_ER"],
        "R_u_s-1": species.get("R_u", 0.0),
        "R_SE_self_s-1": species["R_SE_self"], "R_SE_cross_s-1": species["R_SE_cross"],
        "initial_state": evolution["initial_state"], "RF_enabled": False,
        "density_model": evolution["model"], "probe_model": probe_signals["model"],
        "sampling_limited": evolution["sampling_limited"],
        "temperature_C": common["temperature_C"], "N2_pressure_Torr": common["n2_pressure_torr"],
        "static_field_theta_deg": common["static_field_theta_deg"],
        "static_field_phi_deg": common["static_field_phi_deg"], "static_field_nT": common["static_field_nT"],
        "probe_line": probe["line"], "probe_detuning_MHz": probe["detuning_MHz"],
        "probe_direction": probe["k_axis"], "probe_path_length_cm": probe["path_length_cm"],
        "probe_input_azimuth_deg": probe["azimuth_deg"], "probe_input_ellipticity_deg": probe["ellipticity_deg"],
    }
    for beam in beams:
        name = beam["name"]
        columns.update({f"{name}_on_intensity_uW_cm-2": beam["intensity"],
                        f"{name}_off_intensity_uW_cm-2": 0.0,
                        f"{name}_absolute_frequency_MHz": beam["absolute_frequency_MHz"],
                        f"{name}_direction": beam["k_axis"], f"{name}_polarization": beam["pol"]})
    for index, state in enumerate(species["ground_states"]):
        columns[f"population_F{state['F']:g}_m{state['m']:g}"] = species["population"][:, index]
        for other_index, other in enumerate(species["ground_states"]):
            if other_index <= index or other["F"] != state["F"]:
                continue
            coherence = species["density_matrix"][:, index, other_index]
            name = f"rho_F{state['F']:g}_m{state['m']:g}_m{other['m']:g}"
            columns[name + "_real"] = coherence.real
            columns[name + "_imag"] = coherence.imag
    for name, values in species["moments"].items():
        suffix = "_hbar2" if name.startswith("alignment") else "_hbar" if name.startswith("orientation") else ""
        columns[name + suffix] = values
    for signal in ("rotation", "ellipticity"):
        columns[signal + "_rad"] = probe_signals[signal + "_rad"]
        columns[signal + "_mrad"] = probe_signals[signal + "_rad"] * 1000
    columns["transmission_Iout_over_Iin"] = probe_signals["transmission"]
    columns["absorption_1_minus_Iout_over_Iin"] = probe_signals["absorption"]
    return pd.DataFrame(columns)


def section_export_dataframe(plots, *, species, label, section, t1_s, t2_s):
    """One time row and one qualified column per displayed section curve.

    Plot titles disambiguate manifold/observable curves, and each column names
    its plotted units. Shared event, species and axis metadata appears once.
    Undefined optical samples remain empty instead of being dropped.
    """
    if not plots:
        return pd.DataFrame(columns=["time_s", "section"])
    curves = []
    for plot in plots:
        table = plot_export_dataframe(
            plot["frame"], plot["y_label"], species=species, label=label,
            title=plot["title"], value_unit=plot["value_unit"], t1_s=t1_s, t2_s=t2_s,
        ).set_index("time_s")
        names = list(plot["frame"]["Series"].drop_duplicates())
        curves.append(table[names].rename(columns={
            name: f"{plot['title']} | {name} [{plot['value_unit']}]" for name in names
        }))
    metadata = table.drop(columns=[*names, "observable", "value_units"])
    metadata["section"] = section
    return pd.concat([*curves, metadata], axis=1).reset_index()


def probe_scan_export_dataframe(spectrum, *, species, label, probe, config, common, beams, steady,
                                selected_signals=None):
    """Stationary scan in plotted units; frequency is relative to bare D center."""
    columns = {
        "detuning_MHz": spectrum["detuning_MHz"],
        "alkali": label, "atom": species["atom_name"], "reference_line": probe["line"],
        "detuning_reference": "zero-pressure fine-structure D-line center",
        "pressure_shift_MHz": spectrum["pressure_shift_MHz"],
        "pumps_on": True, "stationary_residual_s-1": steady["stationary_residual_s-1"],
        "steady_state_convergence_time_s": steady["convergence_time_s"],
        "density_cm-3": species["density_cm3"], "temperature_C": common["temperature_C"],
        "N2_pressure_Torr": common["n2_pressure_torr"],
        "static_field_theta_deg": common["static_field_theta_deg"],
        "static_field_phi_deg": common["static_field_phi_deg"], "static_field_nT": common["static_field_nT"],
        "probe_direction": probe["k_axis"], "probe_path_length_cm": probe["path_length_cm"],
        "probe_input_azimuth_deg": probe["azimuth_deg"], "probe_input_ellipticity_deg": probe["ellipticity_deg"],
        "R_ER_s-1": config["R_ER"], "R_u_s-1": config.get("R_u", 0.),
        "probe_width_MHz_per_Torr": config["n2_coeffs"][probe["line"]]["width"],
        "probe_shift_MHz_per_Torr": config["n2_coeffs"][probe["line"]]["shift"],
        "model": "pump-on coherent steady state; non-perturbing Jones frequency scan",
    }
    signal_columns = {
        'Optical rotation': ('rotation_mrad', 1000 * spectrum['rotation_rad']),
        'Ellipticity': ('ellipticity_mrad', 1000 * spectrum['ellipticity_rad']),
        'Absorption': ('absorption_1_minus_Iout_over_Iin', spectrum['absorption']),
    }
    selected = signal_columns if selected_signals is None else selected_signals
    # Keep selected plot order in the CSV too, after the frequency column.
    columns = {'detuning_MHz': columns.pop('detuning_MHz'),
               **{signal_columns[name][0]: signal_columns[name][1] for name in selected}, **columns}
    population = species.get('population')
    if population is None:
        population = np.diag(species['density_matrix'][0]).real
    population = np.asarray(population).reshape(-1)
    columns['population_quantization_axis'] = species.get('zeeman_axis', species['q_axis'])
    for state, value in zip(species['ground_states'], population):
        columns[f"population_F{state['F']:g}_m{state['m']:g}"] = float(value)
    for partner_label in ("A", "B"):
        partner = steady[partner_label]
        if partner is not None:
            columns[f"alkali_{partner_label}_atom"] = partner["atom_name"]
            columns[f"alkali_{partner_label}_density_cm-3"] = partner["density_cm3"]
            for rate in ("R_ER", "R_u", "R_SE_self", "R_SE_cross"):
                columns[f"alkali_{partner_label}_{rate}_s-1"] = partner[rate]
    for beam in beams:
        name = beam["name"]
        columns.update({f"{name}_intensity_uW_cm-2": beam["intensity"],
                        f"{name}_frequency_MHz": beam["absolute_frequency_MHz"],
                        f"{name}_direction": beam["k_axis"], f"{name}_polarization": beam["pol"]})
    return pd.DataFrame(columns)


def steady_density_matrix_export_dataframe(species, *, label):
    """Export every complex matrix element in the stored basis, per total atom."""
    rho = np.asarray(species["density_matrix"])[0]
    states = species["ground_states"]
    rows, columns = np.indices(rho.shape)
    rows, columns = rows.ravel(), columns.ravel()
    return pd.DataFrame({
        "row_index": rows, "column_index": columns,
        "row_F": [states[i]["F"] for i in rows],
        "row_m": [states[i]["m"] for i in rows],
        "column_F": [states[i]["F"] for i in columns],
        "column_m": [states[i]["m"] for i in columns],
        "rho_real": rho.real.ravel(), "rho_imag": rho.imag.ravel(),
        "alkali": label, "atom": species["atom_name"],
        "quantization_axis": species["q_axis"],
        "normalization": "per total atom", "pumps_on": True,
    })
