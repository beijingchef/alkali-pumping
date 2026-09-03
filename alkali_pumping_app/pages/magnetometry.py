"""Magnetometry analysis page for the Alkali Pumping application."""

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from alkali_pumping_app.physics import *
from alkali_pumping_app.physics.optical_pumping import build_optical_L
from alkali_pumping_app.ui.conditions import *
from alkali_pumping_app.ui.downloads import save_button_with_immediate_download
from alkali_pumping_app.ui.exports import (
    dataframe_to_csv_bytes,
    weak_probe_export_dataframe,
    weak_rf_export_dataframe,
)
from alkali_pumping_app.ui.plot_data import aligned_population_bar_data
from alkali_pumping_app.ui.page_state import register_persistent_page_settings
from alkali_pumping_app.ui.rf_display import (
    add_probe_decomposition_legend,
    prepare_weak_rf_plot_values,
    rf_component_legend_label,
)
from alkali_pumping_app.ui.tables import (
    ZEEMAN_COLUMN_LABELS,
    ZEEMAN_DEFAULT_VISIBLE_COLUMN_KEYS,
    ZEEMAN_OPTIONAL_COLUMN_KEYS,
    render_transition_table_html,
    render_zeeman_properties_table_html,
)
from alkali_pumping_app.ui.text import (
    accent_caption,
    inactive_aware_label,
    line_center_detuning_caption,
)
from alkali_pumping_app.ui.uploads import open_file_button


register_persistent_page_settings(
    (
        *CONDITION_KEYS,
        "pump_configuration_tab",
        "probe_configuration_tab",
        "result_species_tab",
        "zeeman_visible_columns_A",
        "zeeman_visible_columns_B",
    )
)

st.html(
    """
    <style>
    [data-testid="stMainBlockContainer"] {
        max-width: none;
        padding-left: 2rem;
        padding-right: 2rem;
    }
    @media (max-width: 700px) {
        [data-testid="stMainBlockContainer"] {
            padding-left: 1rem;
            padding-right: 1rem;
        }
    }
    div[data-testid="stPopoverBody"]:has(.st-key-zeeman_column_options_A),
    div[data-testid="stPopoverBody"]:has(.st-key-zeeman_column_options_B) {
        max-height: calc(100vh - 2rem);
    }
    div[role="listbox"] {
        max-height: calc(100vh - 8rem);
    }
    </style>
    """
)


for _key, _value in DEFAULT_STARTUP_CONDITION.items():
    st.session_state.setdefault(_key, _value)
# A browser session can survive an in-place app upgrade without loading a
# condition file. Preserve the v6.8 meaning of its existing plain pump links
# exactly once; subsequent plain-pump selections are intentional v6.9 modes.
if st.session_state.get("_probe_source_semantics_version") != "6.9":
    for _label in ("A", "B"):
        _source_key = f"probe_source_{_label}"
        _source = st.session_state.get(_source_key)
        if _source in {f"Pump{_label}{_number}" for _number in (1, 2, 3)}:
            st.session_state[_source_key] = f"{_source} weak"
    st.session_state["_probe_source_semantics_version"] = "6.9"
st.session_state.setdefault(
    "_condition_save_name",
    clean_condition_name(st.session_state["condition_name"]),
)


def _translate_transitions_for_alkali_change(label):
    """Keep each pump on the same relative hyperfine branch."""
    previous_atom_name = st.session_state["_last_atom_names_for_defaults"].get(label)
    new_atom_name = st.session_state[f"atom_{label}_name"]
    old_atom = ATOMS.get(previous_atom_name)
    new_atom = ATOMS.get(new_atom_name)
    for pump_number in (1, 2, 3):
        prefix = f"{label}{pump_number}"
        transition_key = f"transition_{prefix}"
        st.session_state[transition_key] = relative_hyperfine_transition_label(
            old_atom,
            new_atom,
            st.session_state.get(f"line_{prefix}", "D1"),
            st.session_state.get(transition_key),
            allowed_only=bool(st.session_state.get("show_allowed_only", True)),
        )
    probe_transition_key = f"probe_transition_{label}"
    st.session_state[probe_transition_key] = relative_hyperfine_transition_label(
        old_atom,
        new_atom,
        st.session_state.get(f"probe_line_{label}", "D1"),
        st.session_state.get(probe_transition_key),
        allowed_only=bool(st.session_state.get("show_allowed_only", True)),
    )
for _key in ("q_axis_A", "q_axis_B", *RF_CONDITION_KEYS):
    # These condition keys were widget-bound before v6.1.6. Reassigning them
    # protects existing sessions from Streamlit's lazy-widget cleanup.
    st.session_state[_key] = st.session_state[_key]
st.session_state.setdefault(
    "_last_atom_names_for_defaults",
    {
        "A": DEFAULT_STARTUP_CONDITION["atom_A_name"],
        "B": DEFAULT_STARTUP_CONDITION["atom_B_name"],
    },
)

def _n2_coefficients(label):
    return {
        "D1": {
            "width": float(st.session_state[f"D1_width_{label}"]),
            "shift": float(st.session_state[f"D1_shift_{label}"]),
        },
        "D2": {
            "width": float(st.session_state[f"D2_width_{label}"]),
            "shift": float(st.session_state[f"D2_shift_{label}"]),
        },
    }


def _initialize_atom_coefficients(label, atom_name):
    """Apply atom defaults only when a real selector changes manually."""
    previous = st.session_state["_last_atom_names_for_defaults"].get(label)
    if previous == atom_name:
        return
    if atom_name == "None":
        st.session_state["_last_atom_names_for_defaults"][label] = atom_name
        return
    for line in ("D1", "D2"):
        st.session_state[f"{line}_width_{label}"] = DEFAULT_N2_COEFFS[atom_name][line]["width"]
        st.session_state[f"{line}_shift_{label}"] = DEFAULT_N2_COEFFS[atom_name][line]["shift"]
    st.session_state["_last_atom_names_for_defaults"][label] = atom_name


def _migrate_legacy_beam_intensity(
    prefix, atom, n2_coeffs, line, transition, det_rel, q_axis_value
):
    legacy_inputs = st.session_state.get("_legacy_pump_inputs", {})
    legacy = legacy_inputs.get(prefix)
    if legacy is None:
        return
    relative_detuning = 0.0 if legacy.get("rate_reference") == "At resonance" else det_rel
    detuning, selected = absolute_detuning_from_transition_choice(
        atom=atom,
        line=line,
        transition_label=transition,
        relative_detuning_MHz=relative_detuning,
        n2_pressure_torr=n2_pressure_torr,
        n2_coeffs=n2_coeffs,
        allowed_only=show_allowed_only,
    )
    k_axis = st.session_state.get(f"k_{prefix}", "x")
    pol_options = allowed_polarizations(k_axis)
    pol = st.session_state.get(f"pol_{prefix}", pol_options[0])
    if pol not in pol_options:
        pol = pol_options[0]
    states = build_ground_states(atom)
    rate_scale = optical_rate_scale_from_intensity(
        atom=atom,
        line=line,
        intensity_uW_cm2=1.0,
        n2_pressure_torr=n2_pressure_torr,
        temperature_C=temperature_C,
        n2_width_MHz_per_torr=n2_coeffs[line]["width"],
    )
    _, info = build_optical_L(
        atom=atom,
        line=line,
        ground_states=states,
        detuning_MHz=detuning,
        pump_rate_s=rate_scale,
        selected_transition=selected,
        k_axis=k_axis,
        pol=pol,
        q_axis=q_axis_value,
        n2_pressure_torr=n2_pressure_torr,
        temperature_C=temperature_C,
        n2_width_MHz_per_torr=n2_coeffs[line]["width"],
        n2_shift_MHz_per_torr=n2_coeffs[line]["shift"],
        normalize_to_selected_total=False,
    )
    indices = np.ix_(info["reference_ground_indices"], info["reference_excited_indices"])
    rate_per_uW = float(info["R_ge"][indices].sum())
    st.session_state[f"intensity_{prefix}"] = (
        max(0.0, float(legacy.get("rate", 0.0))) / rate_per_uW
        if rate_per_uW > 0.0
        else 0.0
    )
    del legacy_inputs[prefix]
    if legacy_inputs:
        st.session_state["_legacy_pump_inputs"] = legacy_inputs
    else:
        st.session_state.pop("_legacy_pump_inputs", None)


def _prepare_beam_state(prefix, atom_name, n2_coeffs, default_Fg, q_axis_value):
    """Normalize a pump's state before any tab-local widgets are instantiated."""
    atom = ATOMS[atom_name]
    line_key = f"line_{prefix}"
    if st.session_state.get(line_key) not in ("D1", "D2"):
        st.session_state[line_key] = "D1"
    line = st.session_state[line_key]

    k_key = f"k_{prefix}"
    if st.session_state.get(k_key) not in ("z", "x", "y"):
        st.session_state[k_key] = "x"
    pol_key = f"pol_{prefix}"
    pol_options = allowed_polarizations(st.session_state[k_key])
    if st.session_state.get(pol_key) not in pol_options:
        st.session_state[pol_key] = pol_options[0]

    transition_options = transition_choice_labels(
        atom, line, n2_pressure_torr, n2_coeffs, allowed_only=show_allowed_only
    )
    transition_key = f"transition_{prefix}"
    if st.session_state.get(transition_key) not in transition_options:
        st.session_state[transition_key] = default_transition_label(
            atom,
            line,
            n2_pressure_torr,
            n2_coeffs,
            default_Fg,
            default_Fg,
            allowed_only=show_allowed_only,
        )
    _migrate_legacy_beam_intensity(
        prefix,
        atom,
        n2_coeffs,
        line,
        st.session_state[transition_key],
        float(st.session_state[f"det_rel_{prefix}"]),
        q_axis_value,
    )
    # These keys were widget-bound before v6.1.2. Reassigning them detaches
    # existing sessions from Streamlit's stale-widget cleanup during upgrade.
    for field in ("line", "transition", "det_rel", "intensity", "k", "pol"):
        state_key = f"{field}_{prefix}"
        st.session_state[state_key] = st.session_state[state_key]


def _beam_from_state(
    prefix, target_label, atom_name, n2_coeffs, active, q_axis_value, placeholder=None
):
    atom = ATOMS[atom_name]
    line = st.session_state[f"line_{prefix}"]
    transition = st.session_state[f"transition_{prefix}"]
    det_rel = float(st.session_state[f"det_rel_{prefix}"])
    detuning, selected = absolute_detuning_from_transition_choice(
        atom=atom,
        line=line,
        transition_label=transition,
        relative_detuning_MHz=det_rel,
        n2_pressure_torr=n2_pressure_torr,
        n2_coeffs=n2_coeffs,
        allowed_only=show_allowed_only,
    )
    return {
        "name": f"Pump{prefix}",
        "target_label": target_label,
        "target_atom": atom_name,
        "line": line,
        "transition_label": transition,
        "selected_transition": selected,
        "detuning_relative": det_rel,
        "detuning": float(detuning),
        "absolute_frequency_MHz": line_center_frequency_MHz(atom, line) + float(detuning),
        "intensity": float(st.session_state[f"intensity_{prefix}"]),
        "k_axis": st.session_state[f"k_{prefix}"],
        "pol": st.session_state[f"pol_{prefix}"],
        "q_axis": q_axis_value,
        "rate_placeholder": placeholder,
        "active": active,
    }


def _pump_widget_key(prefix, field):
    return f"_pump_widget_{field}_{prefix}"


def _store_pump_widget_value(prefix, field):
    """Copy a visible pump widget value into its persistent condition state."""
    previous_intensity = float(st.session_state[f"intensity_{prefix}"])
    st.session_state[f"{field}_{prefix}"] = st.session_state[
        _pump_widget_key(prefix, field)
    ]
    if field == "intensity" or previous_intensity > 0.0:
        st.session_state["_pump_requires_app_rerun"] = True


def _prime_pump_widget(prefix, field):
    """Restore a lazy tab widget from persistent state before it is rendered."""
    widget_key = _pump_widget_key(prefix, field)
    st.session_state[widget_key] = st.session_state[f"{field}_{prefix}"]
    return widget_key


def _beam_config_ui(
    prefix, target_label, atom_name, n2_coeffs, active, default_Fg, q_axis_value
):
    st.markdown(f"#### Pump{prefix}")
    k_axis = st.selectbox(
        "Beam direction",
        ["z", "x", "y"],
        key=_prime_pump_widget(prefix, "k"),
        on_change=_store_pump_widget_value,
        args=(prefix, "k"),
    )
    pol_options = allowed_polarizations(k_axis)
    st.selectbox(
        "Polarization",
        pol_options,
        key=_prime_pump_widget(prefix, "pol"),
        on_change=_store_pump_widget_value,
        args=(prefix, "pol"),
    )
    line = st.selectbox(
        "Reference line",
        ["D1", "D2"],
        key=_prime_pump_widget(prefix, "line"),
        on_change=_store_pump_widget_value,
        args=(prefix, "line"),
    )
    transition_options = transition_choice_labels(
        ATOMS[atom_name], line, n2_pressure_torr, n2_coeffs,
        allowed_only=show_allowed_only,
    )
    transition = st.selectbox(
        "Hyperfine transition",
        transition_options,
        key=_prime_pump_widget(prefix, "transition"),
        on_change=_store_pump_widget_value,
        args=(prefix, "transition"),
    )
    relative_detuning_MHz = st.number_input(
        "Detuning (MHz)",
        step=10.0,
        format="%g",
        key=_prime_pump_widget(prefix, "det_rel"),
        on_change=_store_pump_widget_value,
        args=(prefix, "det_rel"),
        help="Detuning relative to the selected pressure-shifted hyperfine transition.",
    )
    line_center_detuning_MHz, _selected_transition = (
        absolute_detuning_from_transition_choice(
            atom=ATOMS[atom_name],
            line=line,
            transition_label=transition,
            relative_detuning_MHz=relative_detuning_MHz,
            n2_pressure_torr=n2_pressure_torr,
            n2_coeffs=n2_coeffs,
            allowed_only=show_allowed_only,
        )
    )
    st.caption(line_center_detuning_caption(line, line_center_detuning_MHz))
    st.number_input(
        "Intensity (µW/cm²)",
        min_value=0.0,
        step=1.0,
        format="%.1f",
        key=_prime_pump_widget(prefix, "intensity"),
        on_change=_store_pump_widget_value,
        args=(prefix, "intensity"),
    )
    rate_placeholder = st.empty()
    if active:
        saved_caption = st.session_state.get(f"_pump_rate_caption_{prefix}")
        if saved_caption:
            rate_placeholder.caption(accent_caption(saved_caption))
    return _beam_from_state(
        prefix, target_label, atom_name, n2_coeffs, active, q_axis_value,
        rate_placeholder,
    )


@st.fragment
def _pump_configuration_ui(
    atom_A_name,
    atom_B_name,
    active_B,
    n2_coeffs_A,
    n2_coeffs_B,
    q_axis_A,
    q_axis_B,
):
    """Render pump controls independently until a physical beam changes."""
    pump_B_atom_name = atom_B_name if atom_B_name != "None" else atom_A_name
    pump_B_coeffs = n2_coeffs_B if atom_B_name != "None" else n2_coeffs_A
    for prep_args in (
        ("A1", atom_A_name, n2_coeffs_A, 1, q_axis_A),
        ("A2", atom_A_name, n2_coeffs_A, 2, q_axis_A),
        ("A3", atom_A_name, n2_coeffs_A, 2, q_axis_A),
        ("B1", pump_B_atom_name, pump_B_coeffs, 1, q_axis_B),
        ("B2", pump_B_atom_name, pump_B_coeffs, 2, q_axis_B),
        ("B3", pump_B_atom_name, pump_B_coeffs, 2, q_axis_B),
    ):
        _prepare_beam_state(*prep_args)

    st.header("Pump configuration")
    pump_tab_A, pump_tab_B = st.tabs(
        ["Alkali A", "Alkali B"],
        key="pump_configuration_tab",
        on_change="rerun",
    )
    beam_A1 = _beam_from_state("A1", "A", atom_A_name, n2_coeffs_A, True, q_axis_A)
    beam_A2 = _beam_from_state("A2", "A", atom_A_name, n2_coeffs_A, True, q_axis_A)
    beam_A3 = _beam_from_state("A3", "A", atom_A_name, n2_coeffs_A, True, q_axis_A)
    beam_B1 = _beam_from_state(
        "B1", "B", pump_B_atom_name, pump_B_coeffs, active_B, q_axis_B
    )
    beam_B2 = _beam_from_state(
        "B2", "B", pump_B_atom_name, pump_B_coeffs, active_B, q_axis_B
    )
    beam_B3 = _beam_from_state(
        "B3", "B", pump_B_atom_name, pump_B_coeffs, active_B, q_axis_B
    )
    if pump_tab_A.open:
        with pump_tab_A:
            pump_A_col_1, pump_A_col_2, pump_A_col_3 = st.columns(3, gap="xsmall")
            with pump_A_col_1:
                beam_A1 = _beam_config_ui(
                    "A1", "A", atom_A_name, n2_coeffs_A, True, 1, q_axis_A
                )
            with pump_A_col_2:
                beam_A2 = _beam_config_ui(
                    "A2", "A", atom_A_name, n2_coeffs_A, True, 2, q_axis_A
                )
            with pump_A_col_3:
                beam_A3 = _beam_config_ui(
                    "A3", "A", atom_A_name, n2_coeffs_A, True, 2, q_axis_A
                )
    elif pump_tab_B.open:
        with pump_tab_B:
            pump_B_col_1, pump_B_col_2, pump_B_col_3 = st.columns(3, gap="xsmall")
            with pump_B_col_1:
                beam_B1 = _beam_config_ui(
                    "B1", "B", pump_B_atom_name, pump_B_coeffs, active_B, 1, q_axis_B
                )
            with pump_B_col_2:
                beam_B2 = _beam_config_ui(
                    "B2", "B", pump_B_atom_name, pump_B_coeffs, active_B, 2, q_axis_B
                )
            with pump_B_col_3:
                beam_B3 = _beam_config_ui(
                    "B3", "B", pump_B_atom_name, pump_B_coeffs, active_B, 2, q_axis_B
                )
    if st.session_state.pop("_pump_requires_app_rerun", False):
        st.rerun()
    return beam_A1, beam_A2, beam_A3, beam_B1, beam_B2, beam_B3


def _probe_source_options(label):
    options = ["Custom"]
    for number in (1, 2, 3):
        pump_name = f"Pump{label}{number}"
        options.extend((f"{pump_name} weak", pump_name))
    return options


def _probe_source_spec(label, source):
    """Return the linked pump name and readout mode for one source value."""
    if source == "Custom":
        return {"pump_name": None, "mode": "weak"}
    pump_name = source.removesuffix(" weak")
    allowed = {f"Pump{label}{number}" for number in (1, 2, 3)}
    if pump_name not in allowed:
        return {"pump_name": None, "mode": "weak"}
    return {
        "pump_name": pump_name,
        "mode": "weak" if source.endswith(" weak") else "nonlinear",
    }


def _probe_widget_key(state_key):
    return f"_probe_widget_{state_key}"


def _store_probe_widget_value(state_key):
    """Copy a tab-local probe widget into persistent condition state."""
    st.session_state[state_key] = st.session_state[_probe_widget_key(state_key)]


def _prime_probe_widget(state_key):
    """Restore a tab-local probe widget from persistent condition state."""
    widget_key = _probe_widget_key(state_key)
    st.session_state[widget_key] = st.session_state[state_key]
    return widget_key


def _copy_selected_pump_to_probe(label):
    """Synchronize probe beam geometry/spectrum from its selected pump."""
    source_key = f"probe_source_{label}"
    source = st.session_state.get(source_key, "Custom")
    options = _probe_source_options(label)
    if source not in options:
        source = "Custom"
        st.session_state[source_key] = source
    source_spec = _probe_source_spec(label, source)
    if source_spec["pump_name"] is None:
        return
    prefix = source_spec["pump_name"].removeprefix("Pump")
    for probe_field, pump_field in (
        ("line", "line"),
        ("transition", "transition"),
        ("det_rel", "det_rel"),
        ("k", "k"),
    ):
        st.session_state[f"probe_{probe_field}_{label}"] = st.session_state[
            f"{pump_field}_{prefix}"
        ]
    azimuth, ellipticity = preset_ellipse_parameters(
        st.session_state[f"k_{prefix}"], st.session_state[f"pol_{prefix}"]
    )
    st.session_state[f"probe_azimuth_deg_{label}"] = float(azimuth)
    st.session_state[f"probe_ellipticity_deg_{label}"] = float(ellipticity)


def _store_probe_source(label):
    source_key = f"probe_source_{label}"
    _store_probe_widget_value(source_key)
    _copy_selected_pump_to_probe(label)


def _prepare_probe_state(label, atom_name, n2_coeffs):
    _copy_selected_pump_to_probe(label)
    atom = ATOMS[atom_name]
    line_key = f"probe_line_{label}"
    if st.session_state.get(line_key) not in ("D1", "D2"):
        st.session_state[line_key] = "D1"
    line = st.session_state[line_key]
    transition_options = transition_choice_labels(
        atom, line, n2_pressure_torr, n2_coeffs, allowed_only=show_allowed_only
    )
    transition_key = f"probe_transition_{label}"
    if st.session_state.get(transition_key) not in transition_options:
        upper_F = int(round(float(atom["I"]) + 0.5))
        st.session_state[transition_key] = default_transition_label(
            atom,
            line,
            n2_pressure_torr,
            n2_coeffs,
            upper_F,
            upper_F,
            allowed_only=show_allowed_only,
        )
    if st.session_state.get(f"probe_k_{label}") not in ("z", "x", "y"):
        st.session_state[f"probe_k_{label}"] = "x"


def _probe_from_state(label, atom_name, n2_coeffs):
    line = st.session_state[f"probe_line_{label}"]
    transition = st.session_state[f"probe_transition_{label}"]
    det_rel = float(st.session_state[f"probe_det_rel_{label}"])
    detuning, selected = absolute_detuning_from_transition_choice(
        atom=ATOMS[atom_name],
        line=line,
        transition_label=transition,
        relative_detuning_MHz=det_rel,
        n2_pressure_torr=n2_pressure_torr,
        n2_coeffs=n2_coeffs,
        allowed_only=show_allowed_only,
    )
    source = st.session_state[f"probe_source_{label}"]
    source_spec = _probe_source_spec(label, source)
    pump_name = source_spec["pump_name"]
    pump_prefix = pump_name.removeprefix("Pump") if pump_name is not None else None
    return {
        "source": source,
        "mode": source_spec["mode"],
        "pump_name": pump_name,
        "pump_intensity_uW_cm2": (
            float(st.session_state[f"intensity_{pump_prefix}"])
            if pump_prefix is not None
            else 0.0
        ),
        "line": line,
        "transition_label": transition,
        "selected_transition": selected,
        "detuning_MHz": float(detuning),
        "detuning_relative_MHz": det_rel,
        "k_axis": st.session_state[f"probe_k_{label}"],
        "azimuth_deg": float(st.session_state[f"probe_azimuth_deg_{label}"]),
        "ellipticity_deg": float(
            st.session_state[f"probe_ellipticity_deg_{label}"]
        ),
        "path_length_cm": float(st.session_state[f"probe_path_length_cm_{label}"]),
        # Both dynamic ranks are computed once. The result-panel selector
        # chooses either contribution or their coherent complex total.
        "include_scalar": False,
        "include_orientation": True,
        "include_alignment": True,
    }


def _render_probe_config(label, atom_name, n2_coeffs, disabled=False):
    source_key = f"probe_source_{label}"
    source_col, direction_col, path_col = st.columns([0.45, 0.20, 0.35], gap="xsmall")
    with source_col:
        source = st.selectbox(
            "Probe source",
            _probe_source_options(label),
            key=_prime_probe_widget(source_key),
            on_change=_store_probe_source,
            args=(label,),
            disabled=disabled,
            help=(
                "Custom enables an independent weak probe. A 'weak' pump source "
                "copies that pump's optical settings but ignores intensity. A "
                "plain pump source uses the physical pump intensity in a "
                "self-consistent Stokes calculation."
            ),
        )
    source_spec = _probe_source_spec(label, source)
    linked = source_spec["pump_name"] is not None
    with direction_col:
        direction_key = f"probe_k_{label}"
        st.selectbox(
            "Direction", ["z", "x", "y"],
            key=_prime_probe_widget(direction_key),
            on_change=_store_probe_widget_value, args=(direction_key,),
            disabled=disabled or linked,
        )
    with path_col:
        path_length_key = f"probe_path_length_cm_{label}"
        st.number_input(
            "Path length (cm)",
            min_value=0.0,
            step=0.1,
            format="%.3g",
            key=_prime_probe_widget(path_length_key),
            on_change=_store_probe_widget_value, args=(path_length_key,),
            disabled=disabled,
        )
    spectral_line_col, transition_col = st.columns(2, gap="xsmall")
    with spectral_line_col:
        line_key = f"probe_line_{label}"
        line = st.selectbox(
            "Reference line", ["D1", "D2"], key=_prime_probe_widget(line_key),
            on_change=_store_probe_widget_value, args=(line_key,),
            disabled=disabled or linked,
        )
    transition_options = transition_choice_labels(
        ATOMS[atom_name], line, n2_pressure_torr, n2_coeffs,
        allowed_only=show_allowed_only,
    )
    with transition_col:
        transition_key = f"probe_transition_{label}"
        transition = st.selectbox(
            "Hyperfine transition",
            transition_options,
            key=_prime_probe_widget(transition_key),
            on_change=_store_probe_widget_value,
            args=(transition_key,),
            disabled=disabled or linked,
        )
    detuning_key = f"probe_det_rel_{label}"
    detuning_col, center_caption_col = st.columns(
        2, gap="xsmall", vertical_alignment="bottom"
    )
    with detuning_col:
        det_rel = st.number_input(
            "Detuning (MHz)",
            step=10.0,
            format="%g",
            key=_prime_probe_widget(detuning_key),
            on_change=_store_probe_widget_value,
            args=(detuning_key,),
            disabled=disabled or linked,
            help="Detuning relative to the selected pressure-shifted transition.",
        )
    detuning, _ = absolute_detuning_from_transition_choice(
        atom=ATOMS[atom_name],
        line=line,
        transition_label=transition,
        relative_detuning_MHz=det_rel,
        n2_pressure_torr=n2_pressure_torr,
        n2_coeffs=n2_coeffs,
        allowed_only=show_allowed_only,
    )
    with center_caption_col:
        st.caption(line_center_detuning_caption(line, detuning))
    polarization_col, ellipticity_col = st.columns(2, gap="xsmall")
    with polarization_col:
        azimuth_key = f"probe_azimuth_deg_{label}"
        st.number_input(
            "Azimuth (°)",
            min_value=-90.0,
            max_value=90.0,
            step=5.0,
            key=_prime_probe_widget(azimuth_key),
            on_change=_store_probe_widget_value, args=(azimuth_key,),
            disabled=disabled or linked,
        )
    with ellipticity_col:
        ellipticity_key = f"probe_ellipticity_deg_{label}"
        st.number_input(
            "Ellipticity angle (°)",
            min_value=-45.0,
            max_value=45.0,
            step=5.0,
            key=_prime_probe_widget(ellipticity_key),
            on_change=_store_probe_widget_value, args=(ellipticity_key,),
            disabled=disabled or linked,
        )
    if linked:
        if source_spec["mode"] == "weak":
            message = (
                f"Using {source_spec['pump_name']} optical settings as a weak "
                "probe; pump intensity is ignored and path length remains independent."
            )
        else:
            prefix = source_spec["pump_name"].removeprefix("Pump")
            intensity = float(st.session_state[f"intensity_{prefix}"])
            message = (
                f"Using physical {source_spec['pump_name']} with intensity "
                f"{intensity:g} µW/cm² and self-consistent CBOR+LDOR propagation."
            )
        st.caption(accent_caption(message))
    else:
        st.caption(
            accent_caption(
                f"Probe-{label} is a weak, non-perturbing detector for Alkali {label}."
            )
        )


def _probe_configuration_ui(
    atom_A_name, atom_B_name, active_B, n2_coeffs_A, n2_coeffs_B
):
    probe_B_atom = atom_B_name if active_B else atom_A_name
    probe_B_coeffs = n2_coeffs_B if active_B else n2_coeffs_A
    _prepare_probe_state("A", atom_A_name, n2_coeffs_A)
    _prepare_probe_state("B", probe_B_atom, probe_B_coeffs)
    st.header("Probe configuration")
    tab_A, tab_B = st.tabs(
        ["Probe-A", "Probe-B"], key="probe_configuration_tab", on_change="rerun"
    )
    if tab_A.open:
        with tab_A:
            _render_probe_config("A", atom_A_name, n2_coeffs_A)
    elif tab_B.open:
        with tab_B:
            _render_probe_config(
                "B", probe_B_atom, probe_B_coeffs, disabled=not active_B
            )
    return (
        _probe_from_state("A", atom_A_name, n2_coeffs_A),
        _probe_from_state("B", probe_B_atom, probe_B_coeffs),
    )


with st.sidebar:
    st.header("Settings")
    condition_controls_placeholder = st.empty()

    st.header("Atom / cell")
    atom_col_A, atom_col_B = st.columns(2, gap="xsmall")
    with atom_col_A:
        atom_A_name = st.selectbox(
            "Alkali A",
            list(ATOMS),
            key="atom_A_name",
            on_change=_translate_transitions_for_alkali_change,
            args=("A",),
        )
    with atom_col_B:
        atom_B_name = st.selectbox(
            inactive_aware_label(
                "Alkali B",
                st.session_state.get("atom_B_name", "None") == "None",
            ),
            ["None", *list(ATOMS)],
            key="atom_B_name",
            on_change=_translate_transitions_for_alkali_change,
            args=("B",),
        )
    active_B = atom_B_name != "None" and atom_B_name != atom_A_name
    if atom_B_name == atom_A_name:
        st.caption(
            accent_caption(
                "Alkali B matches Alkali A, so its physical effects are ignored.",
                color="orange",
            )
        )

    _initialize_atom_coefficients("A", atom_A_name)
    _initialize_atom_coefficients("B", atom_B_name)
    n2_coeffs_A = _n2_coefficients("A")
    n2_coeffs_B = _n2_coefficients("B")

    density_model_col, density_ratio_col = st.columns(2, gap="xsmall")
    with density_model_col:
        density_mode = st.selectbox(
            "Mixture density model",
            ["Independent saturated-vapor curves", "Relative concentration"],
            format_func=lambda mode: (
                "Independent"
                if mode == "Independent saturated-vapor curves"
                else "Mixed"
            ),
            key="density_mode",
        )
    with density_ratio_col:
        density_ratio = st.number_input(
            "Molar ratio B/A",
            min_value=0.0,
            step=0.1,
            format="%.3g",
            key="density_ratio_B_to_A",
            disabled=(not active_B or density_mode != "Relative concentration"),
            help=(
                "Condensed-phase mole ratio B/A. Under Raoult's law, each "
                "pure saturated-vapor density is multiplied by that alkali's "
                "liquid mole fraction."
            ),
        )

    cell_col_1, cell_col_2 = st.columns(2, gap="xsmall")
    with cell_col_1:
        n2_pressure_torr = st.number_input(
            "N₂ pressure (Torr)", min_value=0.0, step=10.0, format="%.1f", key="n2_pressure_torr"
        )
    with cell_col_2:
        temperature_C = st.number_input(
            "Temperature (°C)", step=1.0, format="%.1f", key="temperature_C_for_table"
        )

    density_A, density_B = resolve_alkali_densities(
        atom_A_name, atom_B_name, temperature_C, density_mode, density_ratio
    )
    er_col_A, er_col_B = st.columns(2, gap="xsmall")
    with er_col_A:
        R_ER_A = st.number_input(
            r"$R_{\mathrm{ER},A}$ (s⁻¹)", min_value=0.0, step=1.0, format="%.1f", key="gamma_ER_A"
        )
    with er_col_B:
        R_ER_B = st.number_input(
            r"$R_{\mathrm{ER},B}$ (s⁻¹)", min_value=0.0, step=1.0, format="%.1f", key="gamma_ER_B"
        )
    field_col_1, field_col_2 = st.columns(2, gap="xsmall")
    with field_col_1:
        static_field_axis = st.selectbox(
            "Static field direction", ["z", "x", "y"], key="static_field_axis"
        )
    with field_col_2:
        static_field_nT = st.number_input(
            "Strength (nT)",
            step=1.0,
            format="%g",
            key="static_field_nT",
            help="A negative value reverses the selected field direction.",
        )

    q_axis_A = st.session_state["q_axis_A"]
    q_axis_B = st.session_state["q_axis_B"]

    show_allowed_only = bool(st.session_state["show_allowed_only"])
    beam_A1, beam_A2, beam_A3, beam_B1, beam_B2, beam_B3 = _pump_configuration_ui(
        atom_A_name,
        atom_B_name,
        active_B,
        n2_coeffs_A,
        n2_coeffs_B,
        q_axis_A,
        q_axis_B,
    )

    probe_A, probe_B = _probe_configuration_ui(
        atom_A_name, atom_B_name, active_B, n2_coeffs_A, n2_coeffs_B
    )

    st.header("Display")
    show_allowed_only = st.checkbox("Only show allowed hyperfine transitions", key="show_allowed_only")
    show_rate_matrices = st.checkbox("Show rate matrices", key="show_rate_matrices")

    with condition_controls_placeholder.container():
        load_col, save_col, name_col = st.columns([0.20, 0.20, 0.60], gap="xsmall")
        with load_col:
            open_file_button(
                type=["json"],
                key="condition_file_uploader",
                help="Load a v6.9 condition or migrate a v6.8/v6.7/v6.6/v6.5/v6.4/v6.3/v6.2/v6.1/v6.0/v5.0 condition.",
                on_change=load_condition_callback,
            )
        with save_col:
            save_placeholder = st.empty()
        with name_col:
            st.text_input(
                "Setting name",
                key="condition_name",
                on_change=sync_condition_save_name,
            )
        # Streamlit prepares the file before a download-button callback runs.
        # Synchronize first so the initial click uses the visible setting name.
        condition_save_name = sync_condition_save_name()
        payload = build_condition_payload(current_condition_values(condition_save_name))
        save_button_with_immediate_download(
            save_placeholder,
            data=json.dumps(payload, indent=2),
            file_name=f"{condition_save_name}.json",
            mime="application/json",
            key="save_condition_button",
        )
        if st.session_state.get("_condition_load_message"):
            st.success(st.session_state.pop("_condition_load_message"))
        if st.session_state.get("_condition_load_error"):
            st.error("Could not load condition file: " + st.session_state.pop("_condition_load_error"))


def _rf_frequency_samples(label):
    normalize_rf_frequency_bounds(label)
    lower = float(st.session_state[f"rf_frequency_lower_hz_{label}"])
    upper = float(st.session_state[f"rf_frequency_upper_hz_{label}"])
    if np.isclose(lower, upper):
        return np.array([lower], dtype=float)
    return np.linspace(lower, upper, 1201)


rf_frequencies_A = _rf_frequency_samples("A")
rf_frequencies_B = _rf_frequency_samples("B")

all_beams = [beam_A1, beam_A2, beam_A3]
if active_B:
    all_beams.extend([beam_B1, beam_B2, beam_B3])


def _physical_beam_config(beam):
    """Strip UI-only objects from an active beam before solver caching."""
    return {
        key: beam[key]
        for key in (
            "name",
            "target_label",
            "absolute_frequency_MHz",
            "intensity",
            "k_axis",
            "pol",
            "selected_transition",
            "transition_label",
        )
    }


@st.cache_data(max_entries=24, show_spinner=False)
def _compute_alkali_system_cached(
    species_A_config, species_B_config, physical_beams, common
):
    return compute_alkali_system(
        species_A_config, species_B_config, physical_beams, common
    )

species_A_config = {
    "label": "A",
    "atom_name": atom_A_name,
    "density_cm3": density_A,
    "R_ER": R_ER_A,
    "n2_coeffs": n2_coeffs_A,
    "q_axis": q_axis_A,
    "rf_axis": st.session_state["rf_axis_A"],
    "rf_observable": st.session_state["rf_observable_A"],
    "rf_frequencies_hz": rf_frequencies_A,
    "probe": probe_A,
}
species_B_config = None
if active_B:
    species_B_config = {
        "label": "B",
        "atom_name": atom_B_name,
        "density_cm3": density_B,
        "R_ER": R_ER_B,
        "n2_coeffs": n2_coeffs_B,
        "q_axis": q_axis_B,
        "rf_axis": st.session_state["rf_axis_B"],
        "rf_observable": st.session_state["rf_observable_B"],
        "rf_frequencies_hz": rf_frequencies_B,
        "probe": probe_B,
    }
common = {
    "temperature_C": temperature_C,
    "n2_pressure_torr": n2_pressure_torr,
    "static_field_axis": static_field_axis,
    "static_field_nT": static_field_nT,
}
physical_beams = [
    _physical_beam_config(beam)
    for beam in all_beams
    if float(beam["intensity"]) > 0.0
]
system = _compute_alkali_system_cached(
    species_A_config, species_B_config, physical_beams, common
)


def _populate_rate_captions():
    for beam in (beam_A1, beam_A2, beam_A3, beam_B1, beam_B2, beam_B3):
        if not beam["active"]:
            continue
        result = system[beam["target_label"]]
        diagnostic = next(
            ((candidate, info) for candidate, info in result["diagnostics"] if candidate["name"] == beam["name"]),
            None,
        )
        if diagnostic is None:
            caption = "pump rate: 0 s⁻¹"
            st.session_state[f"_pump_rate_caption_{beam['name'][4:]}"] = caption
            if beam["rate_placeholder"] is not None:
                beam["rate_placeholder"].caption(accent_caption(caption))
            continue
        _, info = diagnostic
        selected_rows = info["reference_ground_indices"]
        selected_rate = float(info["R_ge"][selected_rows, :].sum())
        caption = f"pump rate: {selected_rate:.3g} s⁻¹"
        st.session_state[f"_pump_rate_caption_{beam['name'][4:]}"] = caption
        if beam["rate_placeholder"] is not None:
            beam["rate_placeholder"].caption(accent_caption(caption))


_populate_rate_captions()


def _zeeman_display_dataframe(result):
    renamed = result["df_pop"].rename(columns={
        "hyperfine_population": "P_F",
        "population": "Pₘ",
        "population_difference": "Dₘ",
        "nu_VS": "ν^{VS} (Hz)",
        "nu_TS": "ν^{TS} (Hz)",
        "nu_LS": "ν^{LS} (Hz)",
        "nu_B": "ν^{B} (Hz)",
        "nu_m": "ν_m (Hz)",
        "Lambda": "Λ (s⁻¹)",
        "G_OP": "G^{OP} (s^-1)",
        "Gamma_OP": "Γ^{OP} (s^-1)",
        "G_ER": "G^{ER} (s^-1)",
        "Gamma_ER": "Γ^{ER} (s^-1)",
        "G_SE_self": "G^{SE,self} (s^-1)",
        "Gamma_SE_self": "Γ^{SE,self} (s^-1)",
        "G_SE_cross": "G^{SE,cross} (s^-1)",
        "Gamma_SE_cross": "Γ^{SE,cross} (s^-1)",
        "G_SE": "G^{SE} (s^-1)",
        "Gamma_SE": "Γ^{SE} (s^-1)",
        "G_total": "G (s^-1)",
        "Gamma_total": "Γ (s^-1)",
        "Gamma_total_over_2pi": "Γ/2π (Hz)",
    })
    columns = [
        "F", "m", "P_F", "Pₘ", "Dₘ", "ν^{VS} (Hz)", "ν^{TS} (Hz)",
        "ν^{LS} (Hz)", "ν^{B} (Hz)", "ν_m (Hz)", "Λ (s⁻¹)",
        "G^{OP} (s^-1)", "Γ^{OP} (s^-1)",
        "G^{ER} (s^-1)", "Γ^{ER} (s^-1)",
        "G^{SE,self} (s^-1)", "Γ^{SE,self} (s^-1)",
        "G^{SE,cross} (s^-1)", "Γ^{SE,cross} (s^-1)",
        "G^{SE} (s^-1)", "Γ^{SE} (s^-1)", "G (s^-1)",
        "Γ (s^-1)", "Γ/2π (Hz)",
    ]
    return renamed[columns].sort_values(["F", "m"], ascending=[False, False], kind="stable").reset_index(drop=True)


def _compact_title(text):
    st.markdown(
        f"<div style='text-align:center;font-size:1.25rem;font-weight:600;margin:.25rem 0 .45rem'>{text}</div>",
        unsafe_allow_html=True,
    )


def _render_population_plot(result):
    states = result["ground_states"]
    population = result["population"]
    common_m_values, groups = aligned_population_bar_data(states, population)
    fig, axes = plt.subplots(
        len(groups), 1, figsize=(4.6, 4.5), sharex=True, sharey=True
    )
    axes = np.atleast_1d(axes)
    for axis, (F, m_values, manifold_population) in zip(axes, groups):
        axis.bar(m_values, manifold_population, width=0.72)
        axis.set_title(f"F={F:g}", fontsize=11, pad=3)
        axis.set_ylabel("Population")
        axis.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.22)
    axes[-1].set_xticks(common_m_values)
    axes[-1].set_xticklabels([f"{m:g}" for m in common_m_values])
    axes[-1].set_xlabel(f"m along {result['q_axis']}")
    if len(common_m_values):
        axes[-1].set_xlim(common_m_values[0] - 0.6, common_m_values[-1] + 0.6)
    axis_max = max(0.01, 1.08 * float(np.max(population)))
    for axis in axes:
        axis.set_ylim(0.0, axis_max)
    fig.tight_layout()
    st.pyplot(fig, width="stretch")
    plt.close(fig)


def _compact_scientific_notation(formatted):
    """Convert e notation to compact LaTeX scientific notation."""
    if "e" not in formatted.lower():
        return formatted
    mantissa, exponent = formatted.lower().split("e", maxsplit=1)
    return rf"{mantissa}\!\times\!10^{{{int(exponent)}}}"


def _format_two_significant(value):
    """Format a finite value with two significant figures and useful zeros."""
    value = float(value)
    if not np.isfinite(value):
        return f"{value:g}"
    if np.isclose(value, 0.0, atol=5e-13):
        return "0"
    exponent = int(np.floor(np.log10(abs(value))))
    if -3 <= exponent < 2:
        decimal_places = max(0, 1 - exponent)
        return f"{value:.{decimal_places}f}"
    return _compact_scientific_notation(f"{value:.1e}")


def _render_population_and_rate_summary(result):
    """Show species moments and the full A/B density and collision summary."""
    states = result["ground_states"]
    population = result["population"]
    result_B = system["B"]
    R_SE_B = result_B["R_SE_self"] if result_B is not None else 0.0
    R_B_from_A = result_B["R_SE_cross"] if result_B is not None else 0.0
    mean_m = _format_two_significant(expectation_m(states, population))
    mean_m2 = _format_two_significant(expectation_m2(states, population))
    formatted_density_A = _format_two_significant(density_A)
    formatted_density_B = _format_two_significant(density_B)
    formatted_R_ER_A = _compact_scientific_notation(f"{R_ER_A:.3g}")
    formatted_R_ER_B = _compact_scientific_notation(f"{R_ER_B:.3g}")
    formatted_R_SE_A = _compact_scientific_notation(
        f"{system['A']['R_SE_self']:.3g}"
    )
    formatted_R_SE_B = _compact_scientific_notation(f"{R_SE_B:.3g}")
    formatted_R_A_from_B = _compact_scientific_notation(
        f"{system['A']['R_SE_cross']:.3g}"
    )
    formatted_R_B_from_A = _compact_scientific_notation(f"{R_B_from_A:.3g}")
    density_A_label = rf"$n_A={formatted_density_A}\ \mathrm{{cm^{{-3}}}}$"
    density_B_label = rf"$n_B={formatted_density_B}\ \mathrm{{cm^{{-3}}}}$"
    if result_B is not None and density_mode == "Relative concentration":
        saturated_A = alkali_vapor_density_cm3(atom_A_name, temperature_C)
        saturated_B = alkali_vapor_density_cm3(atom_B_name, temperature_C)
        fraction_A = density_A / saturated_A if saturated_A > 0.0 else 0.0
        fraction_B = density_B / saturated_B if saturated_B > 0.0 else 0.0
        density_A_label = (
            rf"$n_A={formatted_density_A}\ \mathrm{{cm^{{-3}}}}"
            rf"={_format_two_significant(fraction_A)}\,n_A^{{\mathrm{{sat}}}}$"
        )
        density_B_label = (
            rf"$n_B={formatted_density_B}\ \mathrm{{cm^{{-3}}}}"
            rf"={_format_two_significant(fraction_B)}\,n_B^{{\mathrm{{sat}}}}$"
        )
    summary_caption = (
        rf"$\langle m\rangle={mean_m}$, "
        rf"$\langle m^2\rangle={mean_m2}$, "
        f"{density_A_label}, {density_B_label}, "
        rf"$R_{{\mathrm{{ER}},A}}={formatted_R_ER_A}\ \mathrm{{s^{{-1}}}}$, "
        rf"$R_{{\mathrm{{ER}},B}}={formatted_R_ER_B}\ \mathrm{{s^{{-1}}}}$, "
        rf"$R_{{\mathrm{{SE}},A}}={formatted_R_SE_A}\ \mathrm{{s^{{-1}}}}$, "
        rf"$R_{{\mathrm{{SE}},B}}={formatted_R_SE_B}\ \mathrm{{s^{{-1}}}}$, "
        rf"$R_{{A\leftarrow B}}={formatted_R_A_from_B}\ \mathrm{{s^{{-1}}}}$, "
        rf"$R_{{B\leftarrow A}}={formatted_R_B_from_A}\ \mathrm{{s^{{-1}}}}$"
    )
    st.caption(accent_caption(summary_caption))


def _result_widget_key(state_key):
    return f"_result_widget_{state_key}"


def _store_result_widget_value(state_key):
    """Copy a lazy result-tab widget into persistent condition state."""
    st.session_state[state_key] = st.session_state[_result_widget_key(state_key)]


def _prime_result_widget(state_key):
    """Restore a result-tab widget from persistent state before rendering."""
    widget_key = _result_widget_key(state_key)
    st.session_state[widget_key] = st.session_state[state_key]
    return widget_key


def _store_zeeman_column_visibility(label, column_key, widget_key):
    """Update a species' ordered visible-column list from one checkbox."""
    state_key = f"zeeman_visible_columns_{label}"
    selected = set(st.session_state[state_key])
    if st.session_state[widget_key]:
        selected.add(column_key)
    else:
        selected.discard(column_key)
    st.session_state[state_key] = [
        key for key in ZEEMAN_OPTIONAL_COLUMN_KEYS if key in selected
    ]


def _rf_resonances_in_sweep(result, label):
    lower = float(st.session_state[f"rf_frequency_lower_hz_{label}"])
    upper = float(st.session_state[f"rf_frequency_upper_hz_{label}"])
    frame = result["df_pop"]
    mask = (
        np.isclose(frame["F"].to_numpy(dtype=float), result["rf_upper_F"])
        & np.isfinite(frame["nu_m"].to_numpy(dtype=float))
    )
    selected = frame.loc[mask, ["F", "m", "nu_m", "population_difference", "Gamma_total"]].copy()
    selected["resonance_Hz"] = np.abs(selected["nu_m"].to_numpy(dtype=float))
    selected = selected.loc[
        (selected["resonance_Hz"] >= lower) & (selected["resonance_Hz"] <= upper)
    ]
    if selected.empty:
        return selected
    selected["Transition"] = selected.apply(
        lambda row: f"F={row['F']:g}: m={row['m']:g}→{row['m'] - 1:g}", axis=1
    )
    selected["linewidth_Hz"] = selected["Gamma_total"] / (2.0 * np.pi)
    return selected


def _add_rf_resonance_markers(axis, result, label):
    resonances = _rf_resonances_in_sweep(result, label)
    for frequency in resonances.get("resonance_Hz", []):
        axis.axvline(float(frequency), color="0.45", linewidth=0.8, alpha=0.28)


def _render_shared_rf_controls(_result, label):
    _compact_title(f"Alkali {label} field response")
    axis_col, lower_col, upper_col = st.columns(3, gap="xsmall")
    with axis_col:
        key = f"rf_axis_{label}"
        st.selectbox(
            "RF axis", ["x", "y", "z"], key=_prime_result_widget(key),
            on_change=_store_result_widget_value, args=(key,),
        )
    with lower_col:
        key = f"rf_frequency_lower_hz_{label}"
        st.number_input(
            "Lower RF (Hz)", min_value=0.0, step=1.0, format="%g",
            key=_prime_result_widget(key), on_change=_store_result_widget_value,
            args=(key,),
        )
    with upper_col:
        key = f"rf_frequency_upper_hz_{label}"
        st.number_input(
            "Upper RF (Hz)", min_value=0.0, step=1.0, format="%g",
            key=_prime_result_widget(key), on_change=_store_result_widget_value,
            args=(key,),
        )
    st.caption(
        accent_caption(
            f"These RF-{label} settings are shared by the atomic-moment and Probe-{label} plots."
        )
    )


def _render_rf(result, label):
    rf_axis = st.session_state[f"rf_axis_{label}"]
    rf_observable = st.session_state[f"rf_observable_{label}"]
    observable_symbol = rf_observable_display_label(rf_observable, rf_axis)
    if rf_observable == "Q_ij":
        observable_indices = observable_symbol.removeprefix("Q_")
        observable_title = f"<i>Q</i><sub>{observable_indices}</sub>"
        observable_math = rf"Q_{{{observable_indices}}}"
    else:
        observable_title = observable_symbol
        observable_math = observable_symbol
    rf_show_amplitude = bool(st.session_state[f"rf_show_amplitude_{label}"])
    rf_show_in_phase = bool(st.session_state[f"rf_show_in_phase_{label}"])
    rf_add_pi_in_phase = bool(st.session_state[f"rf_add_pi_in_phase_{label}"])
    rf_show_quadrature = bool(st.session_state[f"rf_show_quadrature_{label}"])
    rf_add_pi_quadrature = bool(
        st.session_state[f"rf_add_pi_quadrature_{label}"]
    )
    rf_relaxation_normalized = bool(
        st.session_state[f"rf_relaxation_normalized_{label}"]
    )
    rf_density_factor = bool(st.session_state[f"rf_density_factor_{label}"])
    reference = result["rf_relaxation_reference"]
    gamma = reference.get("Gamma_m") if rf_relaxation_normalized and reference.get("available", False) else None
    density = result["density_cm3"] if rf_density_factor else None
    plotted = prepare_weak_rf_plot_values(
        result["rf_amplitude"], result["rf_in_phase"], result["rf_quadrature"],
        flip_in_phase=rf_add_pi_in_phase,
        flip_quadrature=rf_add_pi_quadrature,
        relaxation_gamma_s_inv=gamma, density_cm3=density,
    )
    export = weak_rf_export_dataframe(
        frequencies_hz=result["rf_frequencies_hz"],
        susceptibility_amplitude=result["rf_amplitude"],
        susceptibility_in_phase=result["rf_in_phase"],
        susceptibility_quadrature=result["rf_quadrature"],
        plotted_amplitude=plotted[0], plotted_in_phase=plotted[1], plotted_quadrature=plotted[2],
        in_phase_plot_factor=(-1.0 if rf_add_pi_in_phase else 1.0),
        quadrature_plot_factor=(-1.0 if rf_add_pi_quadrature else 1.0),
        relaxation_normalized=gamma is not None,
        normalization_gamma_s_inv=gamma,
        density_factored=rf_density_factor,
        density_cm3=density,
    )
    title_col, download_col, help_col = st.columns(
        [0.74, 0.18, 0.08], gap="small"
    )
    with title_col:
        _compact_title(
            f"Alkali {label} upper-hyperfine ⟨{observable_title}⟩ weak-RF susceptibility "
            f"(F={result['rf_upper_F']:g})"
        )
    with download_col:
        with st.container(horizontal_alignment="right"):
            st.download_button(
                "Download CSV", dataframe_to_csv_bytes(export),
                file_name=f"{condition_save_name}_alkali-{label}_weak-rf.csv",
                mime="text/csv; charset=utf-8", key=f"download_rf_{label}", width="content",
            )
    with help_col:
        with st.popover("❓"):
            st.markdown(
                rf"""
The plot is the linear weak-RF susceptibility of the selected upper-hyperfine
observable,

$$
\chi_{{{observable_math}}}(\omega)
=\frac{{d\langle {observable_math}\rangle}}{{d\Omega_{{\mathrm{{rf}}}}}},
$$

versus RF angular-drive frequency expressed in hertz. The RF field is applied
along the selected RF axis, and the other alkali's RF drive is set to zero.

**X (in phase)** is the response in phase with the applied RF drive.
**Y (quadrature)** is the component shifted by $90^\circ$.
**Amplitude** is $|\chi|=\sqrt{{X^2+Y^2}}$.

For $Q_{{ij}}$, the tensor component is chosen perpendicular to the RF axis:
$Q_{{yz}}$ for RF-$x$, $Q_{{zx}}$ for RF-$y$, and $Q_{{xy}}$ for RF-$z$.

Selecting **Add $\pi$** changes the displayed sign of that component and its
legend from $X$ or $Y$ to $-X$ or $-Y$; it does not change the underlying
calculated susceptibility.

With no optional display factors, the susceptibility is per atom and has units
of $\hbar\,\mathrm{{s}}/\mathrm{{atom}}$. **Relaxation normalized** multiplies
the curves by the selected adjacent-coherence decay rate $\Gamma_m$, removing
the time factor. **Density factor** additionally multiplies by the alkali number
density, converting the result from per atom to per volume.

The CSV contains both the raw phase-convention values and the values actually
plotted after the selected $\pi$ shifts, relaxation normalization, and density
factor.
"""
            )
    control_col, plot_col = st.columns([0.23, 0.77], gap="small")
    with control_col:
        st.caption(
            accent_caption(f"RF-{label} applied; the other RF drive is zero.")
        )
        rf_observable_key = f"rf_observable_{label}"
        st.selectbox(
            "Observable",
            ["Fx", "Fy", "Fz", "Q_ij"],
            key=_prime_result_widget(rf_observable_key),
            on_change=_store_result_widget_value,
            args=(rf_observable_key,),
            format_func=lambda value: rf_observable_display_label(value, rf_axis),
        )
        amplitude_key = f"rf_show_amplitude_{label}"
        st.checkbox(
            "Amplitude",
            key=_prime_result_widget(amplitude_key),
            on_change=_store_result_widget_value,
            args=(amplitude_key,),
        )
        for curve_label, show_field, add_pi_field in (
            ("In phase", "show_in_phase", "add_pi_in_phase"),
            ("Quadrature", "show_quadrature", "add_pi_quadrature"),
        ):
            curve_column, phase_column = st.columns([0.62, 0.38], gap="xxsmall")
            show_key = f"rf_{show_field}_{label}"
            add_pi_key = f"rf_{add_pi_field}_{label}"
            with curve_column:
                st.checkbox(
                    curve_label,
                    key=_prime_result_widget(show_key),
                    on_change=_store_result_widget_value,
                    args=(show_key,),
                )
            with phase_column:
                st.checkbox(
                    "Add π",
                    key=_prime_result_widget(add_pi_key),
                    on_change=_store_result_widget_value,
                    args=(add_pi_key,),
                )
        for checkbox_label, state_key in (
            ("Relaxation normalized", f"rf_relaxation_normalized_{label}"),
            ("Density factor", f"rf_density_factor_{label}"),
        ):
            st.checkbox(
                checkbox_label,
                key=_prime_result_widget(state_key),
                on_change=_store_result_widget_value,
                args=(state_key,),
            )

    with plot_col:
        if not result["static_field_aligned"] and abs(static_field_nT) > 0.0:
            st.warning(
                "The static field is transverse to this quantization axis. "
                "The current population model includes only the field component "
                "parallel to the quantization axis, so transverse static-field "
                "mixing is omitted."
            )
        elif not result["light_shift_available"]:
            st.warning(
                "RF response is unavailable because an active optical field has "
                "a non-diagonal light shift."
            )
        elif result["rf_info"].get("used_transitions", 0) == 0:
            st.info("No driven adjacent Zeeman transitions are available for this RF geometry.")
        else:
            selected_curves = []
            if rf_show_amplitude:
                selected_curves.append((plotted[0], "Amplitude", "-"))
            if rf_show_in_phase:
                selected_curves.append(
                    (
                        plotted[1],
                        rf_component_legend_label("X", rf_add_pi_in_phase),
                        "--",
                    )
                )
            if rf_show_quadrature:
                selected_curves.append(
                    (
                        plotted[2],
                        rf_component_legend_label("Y", rf_add_pi_quadrature),
                        ":",
                    )
                )
            if not selected_curves:
                st.info("Select at least one RF curve.")
            else:
                fig, axis = plt.subplots(figsize=(8.6, 4.2))
                for values, curve_label, linestyle in selected_curves:
                    axis.plot(
                        result["rf_frequencies_hz"], values,
                        linestyle=linestyle, label=curve_label,
                    )
                _add_rf_resonance_markers(axis, result, label)
                axis.axhline(0.0, linewidth=0.8, alpha=0.45)
                axis.set_xlabel(f"RF-{label} frequency (Hz)")
                axis.set_ylabel(
                    f"Alkali {label} "
                    + rf"$d\langle {observable_math}\rangle/d\Omega_{{\mathrm{{rf}}}}$"
                )
                axis.grid(True, alpha=0.25)
                axis.legend(frameon=False)
                fig.tight_layout()
                st.pyplot(fig, width="stretch")
                plt.close(fig)


def _render_probe_response(result, label):
    signal_labels = {
        "rotation": "Optical rotation ψ",
        "ellipticity": "Ellipticity χ",
        "s1": "Normalized Stokes s₁",
        "s2": "Normalized Stokes s₂",
        "s3": "Normalized Stokes s₃",
        "transmission": "Fractional transmission",
    }
    signal = st.session_state[f"probe_signal_{label}"]
    nonlinear = result["probe_info"].get("mode") == "nonlinear physical pump"
    pump_name = result["probe_info"].get("pump_name")
    response = result["probe_response"]
    response_choice = st.session_state[f"probe_response_component_{label}"]
    response_key = {
        "Orientation induced": "orientation",
        "Alignment induced": "alignment",
        "Total": "total",
    }[response_choice]
    selected_response = response[response_key][signal]
    rf_rabi_rad_s_per_nT = 2.0 * np.pi * upper_larmor_frequency_from_field_nT(
        result["atom_name"], 1.0
    )
    add_pi_x = bool(st.session_state[f"probe_add_pi_in_phase_{label}"])
    add_pi_y = bool(st.session_state[f"probe_add_pi_quadrature_{label}"])
    plotted_total = prepare_weak_rf_plot_values(
        selected_response["amplitude"],
        selected_response["in_phase"],
        selected_response["quadrature"],
        flip_in_phase=add_pi_x, flip_quadrature=add_pi_y,
    )
    plotted_total = tuple(
        values * rf_rabi_rad_s_per_nT for values in plotted_total
    )
    export = weak_probe_export_dataframe(
        result["rf_frequencies_hz"], response, signal,
        rf_rabi_rad_s_per_nT=rf_rabi_rad_s_per_nT,
    )
    export["readout_mode"] = "nonlinear physical pump" if nonlinear else "weak probe"
    export["physical_pump"] = pump_name or ""
    title_col, download_col, help_col = st.columns([0.74, 0.18, 0.08], gap="small")
    with title_col:
        readout_title = (
            f"{pump_name} nonlinear optical readout"
            if nonlinear
            else f"Probe-{label} weak optical readout"
        )
        _compact_title(f"{readout_title} — {signal_labels[signal]}")
    with download_col:
        with st.container(horizontal_alignment="right"):
            st.download_button(
                "Download CSV",
                dataframe_to_csv_bytes(export),
                file_name=f"{condition_save_name}_probe-{label}-{signal}.csv",
                mime="text/csv; charset=utf-8",
                key=f"download_probe_{label}",
                width="content",
            )
    with help_col:
        with st.popover("❓"):
            if nonlinear:
                st.markdown(
                    r"""
The selected physical pump is both the perturbing beam and its own detector.
Its RF-induced normalized Stokes vector is propagated self-consistently through
the cell. Rank-1 circular birefringence/dichroism and rank-2 linear
birefringence/dichroism both feed the pump's vector/tensor light shifts.

Orientation induced and Alignment induced are counterfactual nonlinear
solutions with only that optical feedback channel present. Total is the
physical coupled solution; because of feedback, it is not the sum of the two
counterfactual curves. Optical rotation, ellipticity, normalized Stokes
signals, and fractional transmission all come from the same propagated
response.
"""
                )
            else:
                st.markdown(
                    r"""
Probe-A detects Alkali A and Probe-B detects Alkali B. Each probe is treated as
a weak, non-perturbing detector: it does not pump the atoms, add light shift,
or broaden the magnetic resonance. The result is linear in both the probe
polarizability and the RF-induced density-matrix response.

The response selector shows the orientation-induced rank-1 signal, the
alignment-induced rank-2 signal, or their coherent total. For Total, the two
complex responses add before amplitude is calculated. The Mathur rank-2
coefficient is converted to the raw Cartesian Qᵢⱼ convention used by the
atomic-moment solver.

Optical rotation and ellipticity are derived from the normalized transmitted
Stokes vector. X, Y, and amplitude use the same RF phase convention and the
same RF-axis/frequency sweep as the atomic-moment plot above. Responses are
reported per RF magnetic-field amplitude using
$d\Omega_{\mathrm{rf}}/dB_{\mathrm{rf}}=2\pi\nu_L(1\,\mathrm{nT})$.
"""
                )

    control_col, plot_col = st.columns([0.23, 0.77], gap="small")
    with control_col:
        key = f"probe_signal_{label}"
        st.selectbox(
            "Probe signal",
            list(signal_labels),
            key=_prime_result_widget(key),
            on_change=_store_result_widget_value,
            args=(key,),
            format_func=signal_labels.get,
        )
        key = f"probe_response_component_{label}"
        st.selectbox(
            "Atomic response",
            ["Orientation induced", "Alignment induced", "Total"],
            key=_prime_result_widget(key),
            on_change=_store_result_widget_value,
            args=(key,),
        )
        if response_choice == "Total":
            key = f"probe_show_decomposition_{label}"
            st.checkbox(
                "Show rank decomposition",
                key=_prime_result_widget(key),
                on_change=_store_result_widget_value,
                args=(key,),
            )
        key = f"probe_show_amplitude_{label}"
        st.checkbox(
            "Amplitude", key=_prime_result_widget(key),
            on_change=_store_result_widget_value, args=(key,),
        )
        for curve_label, field, pi_field in (
            ("In phase", "show_in_phase", "add_pi_in_phase"),
            ("Quadrature", "show_quadrature", "add_pi_quadrature"),
        ):
            curve_col, phase_col = st.columns([0.62, 0.38], gap="xxsmall")
            key = f"probe_{field}_{label}"
            pi_key = f"probe_{pi_field}_{label}"
            with curve_col:
                st.checkbox(
                    curve_label, key=_prime_result_widget(key),
                    on_change=_store_result_widget_value, args=(key,),
                )
            with phase_col:
                st.checkbox(
                    "Add π", key=_prime_result_widget(pi_key),
                    on_change=_store_result_widget_value, args=(pi_key,),
                )
    with plot_col:
        availability = result["probe_info"]["signal_availability"]
        if signal == "rotation" and not availability["rotation_available"]:
            st.info("Optical rotation is undefined for a purely circular input probe.")
            return
        if signal == "ellipticity" and not availability["ellipticity_available"]:
            st.info("Ellipticity is undefined at a purely circular input state.")
            return
        if not result["light_shift_available"]:
            st.warning("Probe response is unavailable because the RF atomic response is unavailable.")
            return
        if nonlinear and not result["probe_info"].get("nonlinear_available", False):
            st.warning(result["probe_info"].get(
                "nonlinear_reason", "The nonlinear pump response is unavailable."
            ))
            return
        show_amplitude = bool(st.session_state[f"probe_show_amplitude_{label}"])
        show_x = bool(st.session_state[f"probe_show_in_phase_{label}"])
        show_y = bool(st.session_state[f"probe_show_quadrature_{label}"])
        selected = []
        if show_amplitude:
            selected.append((0, "Amplitude", "-"))
        if show_x:
            selected.append((
                1,
                f"In phase ({rf_component_legend_label('X', add_pi_x)})",
                "--",
            ))
        if show_y:
            selected.append((
                2,
                f"Quadrature ({rf_component_legend_label('Y', add_pi_y)})",
                ":",
            ))
        if not selected:
            st.info("Select at least one probe curve.")
            return
        fig, axis = plt.subplots(figsize=(8.6, 4.2))
        legend_entries = {"Total": [], "Orientation": [], "Alignment": []}
        for index, component_label, linestyle in selected:
            line, = axis.plot(
                result["rf_frequencies_hz"], plotted_total[index],
                linestyle=linestyle, label=f"{response_choice} {component_label}",
            )
            legend_entries["Total"].append((line, component_label))
        if (
            response_choice == "Total"
            and st.session_state[f"probe_show_decomposition_{label}"]
        ):
            colors = {"orientation": "C2", "alignment": "C4"}
            for rank in ("orientation", "alignment"):
                rank_result = response[rank][signal]
                plotted_rank = prepare_weak_rf_plot_values(
                    rank_result["amplitude"], rank_result["in_phase"],
                    rank_result["quadrature"], flip_in_phase=add_pi_x,
                    flip_quadrature=add_pi_y,
                )
                plotted_rank = tuple(
                    values * rf_rabi_rad_s_per_nT for values in plotted_rank
                )
                for index, component_label, linestyle in selected:
                    line, = axis.plot(
                        result["rf_frequencies_hz"], plotted_rank[index],
                        color=colors[rank], linestyle=linestyle, alpha=0.75,
                        label=f"{rank.title()} {component_label}",
                    )
                    legend_entries[rank.title()].append((line, component_label))
        _add_rf_resonance_markers(axis, result, label)
        axis.axhline(0.0, linewidth=0.8, alpha=0.45)
        axis.set_xlabel(f"RF-{label} frequency (Hz)")
        signal_axis_labels = {
            "rotation": rf"$d\psi/dB_{{\mathrm{{rf}}}}$ (rad/nT)",
            "ellipticity": rf"$d\chi/dB_{{\mathrm{{rf}}}}$ (rad/nT)",
            "s1": rf"$ds_1/dB_{{\mathrm{{rf}}}}$ (nT$^{{-1}}$)",
            "s2": rf"$ds_2/dB_{{\mathrm{{rf}}}}$ (nT$^{{-1}}$)",
            "s3": rf"$ds_3/dB_{{\mathrm{{rf}}}}$ (nT$^{{-1}}$)",
            "transmission": rf"$dT/dB_{{\mathrm{{rf}}}}$ (nT$^{{-1}}$)",
        }
        axis.set_ylabel(signal_axis_labels[signal])
        axis.grid(True, alpha=0.25)
        if response_choice == "Total" and all(legend_entries.values()):
            add_probe_decomposition_legend(axis, legend_entries)
        else:
            axis.legend(frameon=False)
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)


def _matrix_dataframe(matrix, result, prefix=""):
    labels = [f"{prefix}{state['label']}" for state in result["ground_states"]]
    return pd.DataFrame(matrix, index=labels, columns=labels)


def _render_species_result(result, label):
    q_control, q_caption = st.columns([0.22, 0.78], gap="small")
    with q_control:
        q_axis_key = f"q_axis_{label}"
        st.selectbox(
            f"Alkali {label} quantization axis",
            ["z", "x", "y"],
            key=_prime_result_widget(q_axis_key),
            on_change=_store_result_widget_value,
            args=(q_axis_key,),
        )
    with q_caption:
        st.caption(
            accent_caption(
                f"Shared static field: {static_field_nT:g} nT along {static_field_axis}; "
                f"Alkali {label} upper-manifold Larmor frequency: "
                f"{result['bias_larmor_hz']:.6g} Hz."
            )
        )
    display_df = _zeeman_display_dataframe(result)
    visible_columns_key = f"zeeman_visible_columns_{label}"
    existing_visible_columns = st.session_state.get(
        visible_columns_key, list(ZEEMAN_DEFAULT_VISIBLE_COLUMN_KEYS)
    )
    st.session_state[visible_columns_key] = [
        key for key in ZEEMAN_OPTIONAL_COLUMN_KEYS if key in existing_visible_columns
    ]
    with st.expander(
        f"{result['atom_name']} ground-state Zeeman properties and populations",
        expanded=True,
    ):
        left, right = st.columns([1, 3], gap="small")
        with left:
            _compact_title("Populations")
            _render_population_plot(result)
        with right:
            columns_col, _spacer_col, download_col, help_col = st.columns(
                [0.17, 0.55, 0.20, 0.08], gap="small"
            )
            with columns_col:
                with st.popover("Columns"):
                    st.caption("F and m are always visible.")
                    with st.container(key=f"zeeman_column_options_{label}", gap=None):
                        option_columns = st.columns(2, gap="medium")
                        for index, column_key in enumerate(ZEEMAN_OPTIONAL_COLUMN_KEYS):
                            widget_key = f"_zeeman_column_{label}_{index}"
                            st.session_state[widget_key] = (
                                column_key in st.session_state[visible_columns_key]
                            )
                            with option_columns[index % 2]:
                                st.checkbox(
                                    ZEEMAN_COLUMN_LABELS[column_key],
                                    key=widget_key,
                                    on_change=_store_zeeman_column_visibility,
                                    args=(label, column_key, widget_key),
                                )
            with download_col:
                with st.container(horizontal_alignment="right"):
                    st.download_button(
                        "Download CSV", dataframe_to_csv_bytes(display_df),
                        file_name=f"{condition_save_name}_alkali-{label}_zeeman.csv",
                        mime="text/csv; charset=utf-8", key=f"download_zeeman_{label}", width="content",
                    )
            with help_col:
                with st.popover("❓"):
                    st.markdown(
                        r"""
$F$ and $m$ label the ground-state hyperfine manifold and its Zeeman sublevel.

$P_F=\sum_{m=-F}^{F}P_m$ is the total population of manifold $F$; $P_m$ is the population of $|F,m\rangle$; and $D_m=P_m-P_{m-1}$ is the adjacent-sublevel population difference.

$\nu_{F,m}^{\mathrm{VS}}$ and $\nu_{F,m}^{\mathrm{TS}}$ are the vector and tensor contributions to the light shift of the state $|F,m\rangle$. The total state light shift is
$\nu_{F,m}^{\mathrm{LS}}=\nu_F^{\mathrm{SS}}+\nu_{F,m}^{\mathrm{VS}}+\nu_{F,m}^{\mathrm{TS}}+\nu_{F,m}^{\mathrm{res}}$, where $\nu_F^{\mathrm{SS}}$ is the scalar contribution and $\nu_{F,m}^{\mathrm{res}}$ is any residual left by the scalar-vector-tensor fit. The scalar contribution contains both a common-mode shift shared by the ground hyperfine manifolds and a manifold-dependent hyperfine shift. A common scalar shift cancels from every transition frequency; the scalar shift within one manifold also cancels from its adjacent-Zeeman transition frequencies.

$\nu_{F,m}^{B}=m(g_F/g_{F_+})\nu_{B,+}$ is the static-field Zeeman shift, and
$\nu_m=[\nu_{F,m}^{\mathrm{LS}}+\nu_{F,m}^B]-[\nu_{F,m-1}^{\mathrm{LS}}+\nu_{F,m-1}^B]$
is the total adjacent-sublevel resonance frequency.

$\Lambda_m$ is the total repopulation rate into $|F,m\rangle$, divided by its steady-state population.

$G_m^{\mathrm{OP}}$, $G_m^{\mathrm{ER}}$, and $G_m^{\mathrm{SE}}$ are signed fractional population-loss rates due to optical pumping, electron randomization, and spin exchange. Positive values mean net loss. $\mathrm{SE,self}$ means same-species spin exchange, while $\mathrm{SE,cross}$ means exchange with the other active alkali.

$\Gamma_m^{\mathrm{OP}}$, $\Gamma_m^{\mathrm{ER}}$, and $\Gamma_m^{\mathrm{SE}}$ are the corresponding self-decay rates of the adjacent coherence $\rho_{m,m-1}$. In particular,
$\Gamma_m^{\mathrm{OP}}=(G_m^{\mathrm{OP}}+G_{m-1}^{\mathrm{OP}})/2$.

$G_m$ and $\Gamma_m$ are the sums of all displayed population-loss and adjacent-coherence relaxation mechanisms. $\Gamma_m/(2\pi)$ expresses the total angular decay rate as an ordinary linewidth in hertz.
"""
                    )
            st.markdown(
                render_zeeman_properties_table_html(
                    display_df, st.session_state[visible_columns_key]
                ),
                unsafe_allow_html=True,
            )
        _render_population_and_rate_summary(result)
    _render_shared_rf_controls(result, label)
    with st.expander(f"{result['atom_name']} atomic response", expanded=True):
        _render_rf(result, label)
    with st.expander(f"{result['atom_name']} optical response", expanded=True):
        _render_probe_response(result, label)
    with st.expander("Optical transition frequencies"):
        transition_df = hyperfine_transition_table(
            atom=result["atom"], n2_pressure_torr=n2_pressure_torr,
            n2_coeffs=result["n2_coeffs"], allowed_only=show_allowed_only,
            pump_beams=[
                beam for beam in all_beams if beam["target_label"] == label
            ],
            temperature_C=temperature_C,
        )
        st.markdown(render_transition_table_html(transition_df), unsafe_allow_html=True)
    if show_rate_matrices:
        with st.expander(f"Alkali {label} rate matrices", expanded=False):
            st.markdown("**Optical population generator**")
            st.dataframe(_matrix_dataframe(result["L_optical"], result), width="stretch")
            st.markdown("**Electron-randomization map**")
            st.dataframe(_matrix_dataframe(result["M_ER"], result), width="stretch")
            st.markdown("**Self spin-exchange map**")
            st.dataframe(_matrix_dataframe(result["self_map"], result), width="stretch")
            st.markdown("**Cross spin-exchange map**")
            st.dataframe(_matrix_dataframe(result["cross_map"], result), width="stretch")
            st.markdown("**Local small-signal population Jacobian block**")
            st.dataframe(_matrix_dataframe(result["J_population"], result), width="stretch")


if active_B:
    result_labels = [f"Alkali A — {atom_A_name}", f"Alkali B — {atom_B_name}"]
    if st.session_state.get("result_species_tab") not in result_labels:
        st.session_state["result_species_tab"] = result_labels[0]
    result_tab_A, result_tab_B = st.tabs(
        result_labels,
        key="result_species_tab",
        on_change="rerun",
    )
    if result_tab_A.open:
        with result_tab_A:
            _render_species_result(system["A"], "A")
    elif result_tab_B.open:
        with result_tab_B:
            _render_species_result(system["B"], "B")
else:
    _render_species_result(system["A"], "A")

if show_rate_matrices and active_B:
    with st.expander("Coupled A/B population Jacobian", expanded=False):
        labels_A = [f"A:{state['label']}" for state in system["A"]["ground_states"]]
        labels_B = [f"B:{state['label']}" for state in system["B"]["ground_states"]]
        labels = labels_A + labels_B
        st.dataframe(pd.DataFrame(system["J_coupled"], index=labels, columns=labels), width="stretch")
        st.caption(
            accent_caption(
                "Block order: J_AA, J_AB; J_BA, J_BB. "
                "Off-diagonal blocks are cross-species spin-exchange feedback."
            )
        )

solve = system["solve"]
if not solve.get("converged", True):
    st.warning(
        f"The nonlinear spin-exchange fixed-point iteration did not fully converge "
        f"(residual={solve.get('residual', float('nan')):.3g})."
    )
