"""Condition-file serialization and Streamlit session-state helpers."""

import json
from datetime import datetime

import streamlit as st

from ..version import CONDITION_SCHEMA_VERSION


RF_FIELDS = (
    "axis",
    "observable",
    "frequency_lower_hz",
    "frequency_upper_hz",
    "show_amplitude",
    "show_in_phase",
    "add_pi_in_phase",
    "show_quadrature",
    "add_pi_quadrature",
    "relaxation_normalized",
    "density_factor",
)


def _rf_condition_keys(label):
    return tuple(f"rf_{field}_{label}" for field in RF_FIELDS)


RF_CONDITION_KEYS = (*_rf_condition_keys("A"), *_rf_condition_keys("B"))

PROBE_FIELDS = (
    "source",
    "line",
    "transition",
    "det_rel",
    "k",
    "azimuth_deg",
    "ellipticity_deg",
    "path_length_cm",
    "signal",
    "response_component",
    "show_decomposition",
    "show_amplitude",
    "show_in_phase",
    "add_pi_in_phase",
    "show_quadrature",
    "add_pi_quadrature",
)


def _probe_condition_keys(label):
    return tuple(f"probe_{field}_{label}" for field in PROBE_FIELDS)


PROBE_CONDITION_KEYS = (
    *_probe_condition_keys("A"),
    *_probe_condition_keys("B"),
)


def _pump_condition_keys(prefix):
    return tuple(
        f"{field}_{prefix}"
        for field in ("line", "transition", "det_rel", "intensity", "k", "pol")
    )


PUMP_CONDITION_KEYS = (
    *_pump_condition_keys("A1"),
    *_pump_condition_keys("A2"),
    *_pump_condition_keys("A3"),
    *_pump_condition_keys("B1"),
    *_pump_condition_keys("B2"),
    *_pump_condition_keys("B3"),
)

CONDITION_KEYS = (
    "condition_name",
    "atom_A_name",
    "atom_B_name",
    "density_mode",
    "density_ratio_B_to_A",
    "gamma_ER_A",
    "gamma_ER_B",
    "static_field_axis",
    "static_field_nT",
    "q_axis_A",
    "q_axis_B",
    "temperature_C_for_table",
    "n2_pressure_torr",
    "D1_width_A",
    "D2_width_A",
    "D1_shift_A",
    "D2_shift_A",
    "D1_width_B",
    "D2_width_B",
    "D1_shift_B",
    "D2_shift_B",
    *PUMP_CONDITION_KEYS,
    *RF_CONDITION_KEYS,
    *PROBE_CONDITION_KEYS,
    "show_allowed_only",
    "show_rate_matrices",
)

DEFAULT_STARTUP_CONDITION = {
    "condition_name": "default-dual-alkali",
    "atom_A_name": "Rb87",
    "atom_B_name": "None",
    "density_mode": "Independent saturated-vapor curves",
    "density_ratio_B_to_A": 1.0,
    "gamma_ER_A": 10.0,
    "gamma_ER_B": 10.0,
    "static_field_axis": "z",
    "static_field_nT": 0.0,
    "q_axis_A": "z",
    "q_axis_B": "z",
    "temperature_C_for_table": 23.0,
    "n2_pressure_torr": 0.0,
    "D1_width_A": 17.8,
    "D2_width_A": 18.1,
    "D1_shift_A": -8.25,
    "D2_shift_A": -5.9,
    "D1_width_B": 17.8,
    "D2_width_B": 18.1,
    "D1_shift_B": -8.25,
    "D2_shift_B": -5.9,
    "line_A1": "D1",
    "transition_A1": "2→2",
    "det_rel_A1": 450.0,
    "intensity_A1": 10.0,
    "k_A1": "x",
    "pol_A1": "linear z",
    "line_A2": "D1",
    "transition_A2": "1→2",
    "det_rel_A2": 0.0,
    "intensity_A2": 10.0,
    "k_A2": "x",
    "pol_A2": "linear z",
    "line_A3": "D1",
    "transition_A3": "2→2",
    "det_rel_A3": 0.0,
    "intensity_A3": 0.0,
    "k_A3": "x",
    "pol_A3": "linear z",
    "line_B1": "D1",
    "transition_B1": "2→2",
    "det_rel_B1": 450.0,
    "intensity_B1": 10.0,
    "k_B1": "x",
    "pol_B1": "linear z",
    "line_B2": "D1",
    "transition_B2": "1→2",
    "det_rel_B2": 0.0,
    "intensity_B2": 10.0,
    "k_B2": "x",
    "pol_B2": "linear z",
    "line_B3": "D1",
    "transition_B3": "2→2",
    "det_rel_B3": 0.0,
    "intensity_B3": 0.0,
    "k_B3": "x",
    "pol_B3": "linear z",
    "rf_axis_A": "x",
    "rf_observable_A": "Fx",
    "rf_frequency_lower_hz_A": 0.0,
    "rf_frequency_upper_hz_A": 100.0,
    "rf_show_amplitude_A": True,
    "rf_show_in_phase_A": False,
    "rf_add_pi_in_phase_A": False,
    "rf_show_quadrature_A": False,
    "rf_add_pi_quadrature_A": False,
    "rf_relaxation_normalized_A": False,
    "rf_density_factor_A": False,
    "rf_axis_B": "x",
    "rf_observable_B": "Fx",
    "rf_frequency_lower_hz_B": 0.0,
    "rf_frequency_upper_hz_B": 100.0,
    "rf_show_amplitude_B": True,
    "rf_show_in_phase_B": False,
    "rf_add_pi_in_phase_B": False,
    "rf_show_quadrature_B": False,
    "rf_add_pi_quadrature_B": False,
    "rf_relaxation_normalized_B": False,
    "rf_density_factor_B": False,
    "probe_source_A": "PumpA1 weak",
    "probe_line_A": "D1",
    "probe_transition_A": "2→2",
    "probe_det_rel_A": 1000.0,
    "probe_k_A": "x",
    "probe_azimuth_deg_A": 45.0,
    "probe_ellipticity_deg_A": 0.0,
    "probe_path_length_cm_A": 1.0,
    "probe_signal_A": "rotation",
    "probe_response_component_A": "Total",
    "probe_show_decomposition_A": False,
    "probe_show_amplitude_A": True,
    "probe_show_in_phase_A": False,
    "probe_add_pi_in_phase_A": False,
    "probe_show_quadrature_A": False,
    "probe_add_pi_quadrature_A": False,
    "probe_source_B": "PumpB1 weak",
    "probe_line_B": "D1",
    "probe_transition_B": "2→2",
    "probe_det_rel_B": 1000.0,
    "probe_k_B": "x",
    "probe_azimuth_deg_B": 45.0,
    "probe_ellipticity_deg_B": 0.0,
    "probe_path_length_cm_B": 1.0,
    "probe_signal_B": "rotation",
    "probe_response_component_B": "Total",
    "probe_show_decomposition_B": False,
    "probe_show_amplitude_B": True,
    "probe_show_in_phase_B": False,
    "probe_add_pi_in_phase_B": False,
    "probe_show_quadrature_B": False,
    "probe_add_pi_quadrature_B": False,
    "show_allowed_only": True,
    "show_rate_matrices": False,
}


def clean_condition_name(value):
    name = str(value or "").strip()
    if name.lower().endswith(".json"):
        name = name[:-5].rstrip()
    return name or "default"


def sync_condition_save_name():
    """Synchronize the download filename with the visible condition field."""
    save_name = clean_condition_name(st.session_state.get("condition_name"))
    st.session_state["_condition_save_name"] = save_name
    return save_name


def build_condition_payload(values):
    conditions = {key: values.get(key) for key in CONDITION_KEYS}
    conditions["condition_name"] = clean_condition_name(conditions.get("condition_name"))
    return {
        "app": "alkali_pumping",
        "format": "alkali_pumping_conditions",
        "version": CONDITION_SCHEMA_VERSION,
        "saved_at": datetime.now().isoformat(timespec="seconds"),
        "conditions": conditions,
    }


def current_condition_values(condition_name=None):
    values = {key: st.session_state.get(key) for key in CONDITION_KEYS}
    if condition_name is not None:
        values["condition_name"] = condition_name
    return values


def normalize_rf_frequency_bounds(label="A", prefer="lower"):
    lower_key = f"rf_frequency_lower_hz_{label}"
    upper_key = f"rf_frequency_upper_hz_{label}"
    lower = max(0.0, float(st.session_state.get(lower_key, 0.0)))
    upper = max(0.0, float(st.session_state.get(upper_key, lower)))
    if lower > upper:
        if prefer == "upper":
            lower = upper
        else:
            upper = lower
    st.session_state[lower_key] = lower
    st.session_state[upper_key] = upper


def apply_loaded_condition_dict(payload):
    if not isinstance(payload, dict):
        raise ValueError("The loaded file is not a JSON object.")
    if payload.get("app") != "alkali_pumping":
        raise ValueError("This is not an alkali_pumping condition file.")
    if payload.get("format") != "alkali_pumping_conditions":
        raise ValueError("The JSON file is not an alkali_pumping condition file.")
    conditions = payload.get("conditions")
    if not isinstance(conditions, dict):
        raise ValueError("The JSON file does not contain a conditions object.")

    version = payload.get("version")
    if version != CONDITION_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported condition-file version; expected {CONDITION_SCHEMA_VERSION}."
        )
    loaded_conditions = dict(conditions)
    missing = [key for key in CONDITION_KEYS if key not in loaded_conditions]
    if missing:
        raise ValueError("The condition file is missing required fields: " + ", ".join(missing))

    loaded_name = clean_condition_name(loaded_conditions["condition_name"])
    for key in CONDITION_KEYS:
        value = loaded_conditions.get(key)
        if value is not None:
            st.session_state[key] = value
    st.session_state["_condition_save_name"] = loaded_name
    normalize_rf_frequency_bounds("A")
    normalize_rf_frequency_bounds("B")
    st.session_state["_last_atom_names_for_defaults"] = {
        "A": loaded_conditions["atom_A_name"],
        "B": loaded_conditions["atom_B_name"],
    }
    return loaded_name


def load_condition_callback():
    uploaded = st.session_state.get("condition_file_uploader")
    if uploaded is None:
        return
    try:
        payload = json.loads(uploaded.getvalue().decode("utf-8"))
        loaded_name = apply_loaded_condition_dict(payload)
        st.session_state["_condition_load_message"] = f"Loaded condition: {loaded_name}"
        st.session_state.pop("_condition_load_error", None)
    except Exception as exc:
        st.session_state["_condition_load_error"] = str(exc)
        st.session_state.pop("_condition_load_message", None)
