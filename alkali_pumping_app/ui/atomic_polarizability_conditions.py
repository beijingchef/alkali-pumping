"""Save/load support for the atomic-polarizability analysis page."""

import json
from datetime import datetime

import streamlit as st

from ..physics.constants import DEFAULT_N2_COEFFS


ATOMIC_POLARIZABILITY_CONDITION_VERSION = "1.0"
ATOMIC_POLARIZABILITY_PREFIX = "ap_"

ATOMIC_POLARIZABILITY_DEFAULTS = {
    "condition_name": "default-atomic-polarizability",
    "atom_name": "Rb87",
    "temperature_C": 23.0,
    "n2_pressure_torr": 0.0,
    "line": "D1",
    "lower_MHz": -4600.0,
    "upper_MHz": 6100.0,
    "points": 401,
    "plot_alpha_eq": False,
    "plot_alpha_hfs": False,
    "plot_alpha_gt": True,
    "plot_alpha_br": True,
}


def state_key(field):
    return f"{ATOMIC_POLARIZABILITY_PREFIX}{field}"


def clean_atomic_polarizability_condition_name(value):
    name = str(value or "").strip()
    if name.lower().endswith(".json"):
        name = name[:-5].rstrip()
    return name or "atomic-polarizability"


def initialize_atomic_polarizability_conditions():
    for field, value in ATOMIC_POLARIZABILITY_DEFAULTS.items():
        st.session_state.setdefault(state_key(field), value)
    st.session_state.setdefault(
        "_ap_condition_save_name",
        clean_atomic_polarizability_condition_name(
            st.session_state[state_key("condition_name")]
        ),
    )


def current_atomic_polarizability_values(condition_name=None):
    values = {
        field: st.session_state.get(state_key(field), default)
        for field, default in ATOMIC_POLARIZABILITY_DEFAULTS.items()
    }
    if condition_name is not None:
        values["condition_name"] = condition_name
    return values


def build_atomic_polarizability_payload(values):
    conditions = {
        field: values.get(field, default)
        for field, default in ATOMIC_POLARIZABILITY_DEFAULTS.items()
    }
    conditions["condition_name"] = clean_atomic_polarizability_condition_name(
        conditions["condition_name"]
    )
    return {
        "app": "alkali_pumping",
        "format": "alkali_pumping_atomic_polarizability_conditions",
        "version": ATOMIC_POLARIZABILITY_CONDITION_VERSION,
        "saved_at": datetime.now().isoformat(timespec="seconds"),
        "conditions": conditions,
    }


def apply_atomic_polarizability_payload(payload):
    if not isinstance(payload, dict):
        raise ValueError("The loaded file is not a JSON object.")
    if payload.get("app") != "alkali_pumping":
        raise ValueError("This is not an alkali_pumping condition file.")
    if payload.get("format") != "alkali_pumping_atomic_polarizability_conditions":
        raise ValueError("This is not an atomic-polarizability condition file.")
    if payload.get("version") != ATOMIC_POLARIZABILITY_CONDITION_VERSION:
        raise ValueError(
            "Unsupported atomic-polarizability condition version; expected "
            f"{ATOMIC_POLARIZABILITY_CONDITION_VERSION}."
        )
    conditions = payload.get("conditions")
    if not isinstance(conditions, dict):
        raise ValueError("The JSON file does not contain a conditions object.")
    missing = [
        field for field in ATOMIC_POLARIZABILITY_DEFAULTS if field not in conditions
    ]
    if missing:
        raise ValueError("The condition file is missing required fields: " + ", ".join(missing))

    atom_name = conditions["atom_name"]
    if atom_name not in DEFAULT_N2_COEFFS:
        raise ValueError(f"Unsupported atom: {atom_name}")
    for field in ATOMIC_POLARIZABILITY_DEFAULTS:
        st.session_state[state_key(field)] = conditions[field]
    for line in ("D1", "D2"):
        for coefficient in ("width", "shift"):
            st.session_state[state_key(f"{line}_{coefficient}")] = (
                DEFAULT_N2_COEFFS[atom_name][line][coefficient]
            )
    loaded_name = clean_atomic_polarizability_condition_name(
        conditions["condition_name"]
    )
    st.session_state[state_key("condition_name")] = loaded_name
    st.session_state["_ap_condition_save_name"] = loaded_name
    st.session_state["_ap_loaded_preserve_range"] = True
    return loaded_name


def load_atomic_polarizability_callback(
    upload_key="atomic_polarizability_condition_upload",
):
    uploaded = st.session_state.get(upload_key)
    if uploaded is None:
        return
    try:
        payload = json.loads(uploaded.getvalue().decode("utf-8"))
        loaded_name = apply_atomic_polarizability_payload(payload)
        st.session_state["_ap_condition_load_message"] = (
            f"Loaded condition: {loaded_name}"
        )
        st.session_state.pop("_ap_condition_load_error", None)
    except Exception as exc:
        st.session_state["_ap_condition_load_error"] = str(exc)
        st.session_state.pop("_ap_condition_load_message", None)
