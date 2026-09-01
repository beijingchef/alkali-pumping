"""Save/load support for the independent light-shift analysis page."""

import json
from datetime import datetime

import streamlit as st

from ..physics.constants import DEFAULT_N2_COEFFS


LIGHT_SHIFT_CONDITION_VERSION = "1.1"
LEGACY_LIGHT_SHIFT_CONDITION_VERSION = "1.0"
LIGHT_SHIFT_PREFIX = "ls_"

LIGHT_SHIFT_DEFAULTS = {
    "condition_name": "default-light-shift",
    "atom_name": "Rb87",
    "temperature_C": 23.0,
    "n2_pressure_torr": 0.0,
    "static_field_axis": "z",
    "static_field_nT": 0.0,
    "D1_width": DEFAULT_N2_COEFFS["Rb87"]["D1"]["width"],
    "D1_shift": DEFAULT_N2_COEFFS["Rb87"]["D1"]["shift"],
    "D2_width": DEFAULT_N2_COEFFS["Rb87"]["D2"]["width"],
    "D2_shift": DEFAULT_N2_COEFFS["Rb87"]["D2"]["shift"],
    "line": "D1",
    "k_axis": "z",
    "polarization_mode": "Preset",
    "preset": "sigma+",
    "azimuth_deg": 0.0,
    "ellipticity_deg": 0.0,
    "intensity_uW_cm2": 1.0,
    "reference": "Zero-pressure line center",
    "lower_MHz": -10000.0,
    "upper_MHz": 10000.0,
    "points": 401,
    "normalization": "Per intensity",
    "view": "Components",
    "transition_quantity": "Frequency shift",
    "y_scale": "Linear",
    "show_scalar": True,
    "show_scattering": False,
    "state_manifolds": None,
    "state_components": ["Total diagonal"],
}


def state_key(field):
    return f"{LIGHT_SHIFT_PREFIX}{field}"


def clean_light_shift_condition_name(value):
    name = str(value or "").strip()
    if name.lower().endswith(".json"):
        name = name[:-5].rstrip()
    return name or "light-shift"


def initialize_light_shift_conditions():
    for field, value in LIGHT_SHIFT_DEFAULTS.items():
        st.session_state.setdefault(state_key(field), value)
    st.session_state.setdefault(
        "_ls_condition_save_name",
        clean_light_shift_condition_name(st.session_state[state_key("condition_name")]),
    )
    st.session_state.setdefault("_ls_previous_atom", st.session_state[state_key("atom_name")])


def current_light_shift_values(condition_name=None):
    values = {
        field: st.session_state.get(state_key(field), default)
        for field, default in LIGHT_SHIFT_DEFAULTS.items()
    }
    if condition_name is not None:
        values["condition_name"] = condition_name
    return values


def build_light_shift_payload(values):
    conditions = {
        field: values.get(field, default)
        for field, default in LIGHT_SHIFT_DEFAULTS.items()
    }
    conditions["condition_name"] = clean_light_shift_condition_name(
        conditions["condition_name"]
    )
    return {
        "app": "alkali_pumping",
        "format": "alkali_pumping_light_shift_conditions",
        "version": LIGHT_SHIFT_CONDITION_VERSION,
        "saved_at": datetime.now().isoformat(timespec="seconds"),
        "conditions": conditions,
    }


def apply_light_shift_payload(payload):
    if not isinstance(payload, dict):
        raise ValueError("The loaded file is not a JSON object.")
    if payload.get("app") != "alkali_pumping":
        raise ValueError("This is not an alkali_pumping condition file.")
    if payload.get("format") != "alkali_pumping_light_shift_conditions":
        raise ValueError("This is not a light-shift condition file.")
    payload_version = payload.get("version")
    if payload_version not in (
        LIGHT_SHIFT_CONDITION_VERSION,
        LEGACY_LIGHT_SHIFT_CONDITION_VERSION,
    ):
        raise ValueError(
            f"Unsupported light-shift condition version; expected {LIGHT_SHIFT_CONDITION_VERSION}."
        )
    conditions = payload.get("conditions")
    if not isinstance(conditions, dict):
        raise ValueError("The JSON file does not contain a conditions object.")
    conditions = dict(conditions)
    if payload_version == LEGACY_LIGHT_SHIFT_CONDITION_VERSION:
        conditions.setdefault("show_scalar", True)
    missing = [field for field in LIGHT_SHIFT_DEFAULTS if field not in conditions]
    if missing:
        raise ValueError("The condition file is missing required fields: " + ", ".join(missing))

    for field in LIGHT_SHIFT_DEFAULTS:
        st.session_state[state_key(field)] = conditions[field]
    loaded_name = clean_light_shift_condition_name(conditions["condition_name"])
    st.session_state[state_key("condition_name")] = loaded_name
    st.session_state["_ls_condition_save_name"] = loaded_name
    st.session_state["_ls_previous_atom"] = conditions["atom_name"]
    st.session_state["_ls_loaded_preserve_range"] = True
    return loaded_name


def load_light_shift_callback(upload_key="light_shift_condition_upload"):
    uploaded = st.session_state.get(upload_key)
    if uploaded is None:
        return
    try:
        payload = json.loads(uploaded.getvalue().decode("utf-8"))
        loaded_name = apply_light_shift_payload(payload)
        st.session_state["_ls_condition_load_message"] = f"Loaded condition: {loaded_name}"
        st.session_state.pop("_ls_condition_load_error", None)
    except Exception as exc:
        st.session_state["_ls_condition_load_error"] = str(exc)
        st.session_state.pop("_ls_condition_load_message", None)
