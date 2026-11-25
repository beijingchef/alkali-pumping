"""Time-evolution sidebar: ACC optical/cell controls with isolated state keys."""

from copy import deepcopy
import json

import streamlit as st

from ..physics.constants import ATOMS, DEFAULT_N2_COEFFS
from ..physics.angular_momentum import build_ground_states
from ..physics.optical_pumping import build_optical_L, optical_rate_scale_from_intensity
from ..physics.polarization import allowed_polarizations, preset_ellipse_parameters
from ..physics.spectroscopy import (
    absolute_detuning_from_transition_choice, default_transition_label,
    line_center_frequency_MHz, relative_hyperfine_transition_label, transition_choice_labels,
)
from ..physics.system_model import resolve_alkali_densities
from .downloads import save_button_with_immediate_download, update_settings_download_snapshot
from .page_state import register_persistent_page_settings
from .magnetic_field import render_field_controls
from .probe_presets import PROBE_PRESET_SOURCES, probe_preset_settings
from .text import line_center_detuning_caption, pump_rate_caption
from .time_evolution_conditions import (
    DEFAULTS, PREFIX, build_settings_payload, parse_settings_payload,
)
from .uploads import open_file_button


def key(name):
    return PREFIX + name


def initialize_settings():
    for name, value in DEFAULTS.items():
        st.session_state.setdefault(key(name), deepcopy(value))
    register_persistent_page_settings([key(name) for name in DEFAULTS])
    st.session_state.setdefault(key("last_atoms"), {"A": "Rb87", "B": "None"})


def settings_values():
    return {name: st.session_state[key(name)] for name in DEFAULTS}


def _prepare_settings_download():
    payload = build_settings_payload(settings_values())
    return json.dumps(payload, indent=2), payload["settings"]["condition_name"] + ".json"


def refresh_settings_download():
    """Keep browser-only Save current after a display-section fragment event."""
    snapshot, _ = _prepare_settings_download()
    update_settings_download_snapshot(key("save"), snapshot)


def _clamp_event_times():
    total = st.session_state[key("t2_s")]
    for prefix in (f"{label}{number}" for label in ("A", "B") for number in (1, 2, 3)):
        start_key, duration_key = key(f"pump_start_s_{prefix}"), key(f"pump_duration_s_{prefix}")
        st.session_state[start_key] = min(st.session_state[start_key], total)
        st.session_state[duration_key] = min(st.session_state[duration_key], total - st.session_state[start_key])
    start_key, duration_key = key("probe_start_s"), key("probe_duration_s")
    st.session_state[start_key] = min(st.session_state[start_key], total)
    st.session_state[duration_key] = min(st.session_state[duration_key], total - st.session_state[start_key])
    start_key, duration_key = key("field_start_s"), key("field_duration_s")
    st.session_state[start_key] = min(st.session_state[start_key], total)
    st.session_state[duration_key] = min(st.session_state[duration_key], total - st.session_state[start_key])


def _store_time_widget(name, suffix=""):
    st.session_state[key(name)] = st.session_state["_time_widget_" + key(name) + suffix]
    _clamp_event_times()


def _time_controls():
    st.session_state["_time_widget_" + key("t2_s")] = st.session_state[key("t2_s")]
    st.slider("Total time (s)", min_value=0.001, max_value=30.0, step=0.001,
              format="%.3f", key="_time_widget_" + key("t2_s"),
              on_change=_store_time_widget, args=("t2_s",))


def _probe_time_controls(label, disabled):
    # Separate durable settings from widget identity so moving the controls or
    # changing their bounds cannot reset loaded values to the slider minimum.
    # Both probe tabs edit the same existing timing settings.
    suffix = "" if label == "A" else "_B"
    for name in ("probe_start_s", "probe_duration_s"):
        st.session_state["_time_widget_" + key(name) + suffix] = st.session_state[key(name)]
    start_col, duration_col = st.columns(2, gap="xsmall")
    with start_col:
        st.number_input("Probe start (s)", min_value=0.0,
                        max_value=float(st.session_state[key("t2_s")]), step=0.001,
                        format="%.3f", key="_time_widget_" + key("probe_start_s") + suffix,
                        on_change=_store_time_widget, args=("probe_start_s", suffix), disabled=disabled)
    with duration_col:
        st.number_input("Probe duration (s)", min_value=0.0,
                        max_value=float(st.session_state[key("t2_s")] - st.session_state[key("probe_start_s")]),
                        step=0.001, format="%.3f", key="_time_widget_" + key("probe_duration_s") + suffix,
                        on_change=_store_time_widget, args=("probe_duration_s", suffix), disabled=disabled)


def _file_widget_key():
    revision = st.session_state.get(key("file_revision"), 0)
    return key(f"file_{revision}")


def _load_settings():
    uploaded = st.session_state.get(_file_widget_key())
    if uploaded is None:
        return
    try:
        values = parse_settings_payload(json.loads(uploaded.getvalue().decode("utf-8")))
        for name, value in values.items():
            st.session_state[key(name)] = value
        st.session_state[key("last_atoms")] = {label: values[f"atom_{label}_name"] for label in ("A", "B")}
        st.session_state[key("load_error")] = None
        st.session_state[key("load_success")] = "Loaded Time evolution settings."
    except (ValueError, UnicodeDecodeError) as exc:
        st.session_state[key("load_error")] = str(exc)
    finally:
        # A fresh uploader keeps Open available after success or failure,
        # instead of replacing the button with Streamlit's uploaded-file chip.
        st.session_state[key("file_revision")] = st.session_state.get(key("file_revision"), 0) + 1


def _atom_changed(label):
    previous_name = st.session_state[key("last_atoms")][label]
    # An inactive B displays A's spectrum; use that previous reference when
    # activating B so the lower/upper manifold branches are retained.
    previous_atom = ATOMS.get(previous_name, ATOMS[st.session_state[key("atom_A_name")]])
    name = st.session_state[key(f"atom_{label}_name")]
    if name == "None":
        st.session_state[key("last_atoms")][label] = name
        return
    affected_labels = [label]
    if label == "A" and st.session_state[key("atom_B_name")] == "None":
        affected_labels.append("B")
    for affected in affected_labels:
        for number in (1, 2, 3):
            prefix = f"{affected}{number}"
            transition_key = key(f"transition_{prefix}")
            st.session_state[transition_key] = relative_hyperfine_transition_label(
                previous_atom, ATOMS[name], st.session_state[key(f"line_{prefix}")],
                st.session_state[transition_key], allowed_only=True,
            )
        probe_transition = key(f"probe_transition_{affected}")
        st.session_state[probe_transition] = relative_hyperfine_transition_label(
            previous_atom, ATOMS[name], st.session_state[key(f"probe_line_{affected}")],
            st.session_state[probe_transition], allowed_only=True,
        )
    for line in ("D1", "D2"):
        for field in ("width", "shift"):
            st.session_state[key(f"{line}_{field}_{label}")] = DEFAULT_N2_COEFFS[name][line][field]
    st.session_state[key("last_atoms")][label] = name


def _coefficients(label):
    return {line: {field: float(st.session_state[key(f"{line}_{field}_{label}")])
                   for field in ("width", "shift")} for line in ("D1", "D2")}


def _transition(name, line, transition_key, pressure, coefficients, upper=True, excited_upper=True):
    options = transition_choice_labels(
        ATOMS[name], line, pressure, coefficients,
        allowed_only=True,
    )
    if st.session_state[transition_key] not in options:
        F = int(round(float(ATOMS[name]["I"]) + 0.5))
        st.session_state[transition_key] = default_transition_label(
            ATOMS[name], line, pressure, coefficients, F if upper else F - 1,
            F if excited_upper else F - 1,
            allowed_only=True,
        )
    return options


def _detuning(name, line, transition, relative, pressure, coefficients):
    return absolute_detuning_from_transition_choice(
        atom=ATOMS[name], line=line, transition_label=transition,
        relative_detuning_MHz=relative, n2_pressure_torr=pressure, n2_coeffs=coefficients,
        allowed_only=True,
    )


def _pump(label, number, name, coefficients, pressure, disabled):
    prefix = f"{label}{number}"
    st.markdown(f"#### Pump{prefix}")
    k_axis = st.selectbox("Beam direction", ["z", "x", "y"], key=key(f"k_{prefix}"), disabled=disabled)
    polarizations = allowed_polarizations(k_axis)
    if st.session_state[key(f"pol_{prefix}")] not in polarizations:
        st.session_state[key(f"pol_{prefix}")] = polarizations[0]
    pol = st.selectbox("Polarization", polarizations, key=key(f"pol_{prefix}"), disabled=disabled)
    line = st.selectbox("Reference line", ["D1", "D2"], key=key(f"line_{prefix}"), disabled=disabled)
    transition_key = key(f"transition_{prefix}")
    options = _transition(name, line, transition_key, pressure, coefficients, upper=number != 1)
    transition = st.selectbox("F→F′", options, key=transition_key, disabled=disabled)
    relative = st.number_input("Detuning (MHz)", step=10.0, format="%g", key=key(f"det_rel_{prefix}"), disabled=disabled)
    detuning, selected = _detuning(name, line, transition, relative, pressure, coefficients)
    st.caption(line_center_detuning_caption(line, detuning, pump=True))
    intensity = st.number_input("I (µW/cm²)", min_value=0.0, step=1.0, format="%.1f",
                                key=key(f"intensity_{prefix}"), disabled=disabled)
    if not disabled and intensity > 0:
        atom = ATOMS[name]
        rate_scale = optical_rate_scale_from_intensity(
            atom, line, intensity, pressure, st.session_state[key("temperature_C_for_table")],
            coefficients[line]["width"],
        )
        _, info = build_optical_L(
            atom=atom, line=line, ground_states=build_ground_states(atom),
            detuning_MHz=detuning, pump_rate_s=rate_scale,
            selected_transition=selected, k_axis=k_axis, pol=pol, q_axis="z",
            n2_pressure_torr=pressure, temperature_C=st.session_state[key("temperature_C_for_table")],
            n2_width_MHz_per_torr=coefficients[line]["width"],
            n2_shift_MHz_per_torr=coefficients[line]["shift"],
            normalize_to_selected_total=False,
        )
        rate = float(info["R_ge"][info["reference_ground_indices"], :].sum())
    else:
        rate = 0.0
    st.caption(pump_rate_caption(rate))
    start = st.number_input("Start time (s)", min_value=0.0,
                            max_value=float(st.session_state[key("t2_s")]), step=0.001,
                            format="%.3f", key=key(f"pump_start_s_{prefix}"), disabled=disabled)
    duration = st.number_input("Duration (s)", min_value=0.0,
                               max_value=float(st.session_state[key("t2_s")] - start), step=0.001,
                               format="%.3f", key=key(f"pump_duration_s_{prefix}"), disabled=disabled)
    return {"name": f"Pump{prefix}", "target_label": label,
            "absolute_frequency_MHz": line_center_frequency_MHz(ATOMS[name], line) + detuning,
            "intensity": intensity, "k_axis": k_axis, "pol": pol,
            "start_s": start, "duration_s": duration,
            "selected_transition": selected, "transition_label": transition}


def _apply_probe_source(label):
    """Populate a preset once on selection, retaining subsequent user edits."""
    source = st.session_state[key(f"probe_source_{label}")]
    name = st.session_state[key(f"atom_{label}_name")]
    if name == "None":
        name = st.session_state[key("atom_A_name")]
    for field, value in probe_preset_settings(source, ATOMS[name]).items():
        st.session_state[key(f"probe_{field}_{label}")] = value


def _probe(label, name, coefficients, pressure, disabled):
    source_col, direction_col, path_col = st.columns([0.45, 0.20, 0.35], gap="xsmall")
    with source_col:
        source = st.selectbox(
            "Probe source", ["Custom", *PROBE_PRESET_SOURCES, *[f"Pump{label}{n}" for n in (1, 2, 3)]],
            key=key(f"probe_source_{label}"), disabled=disabled,
            on_change=_apply_probe_source, args=(label,),
        )
    preset = source in PROBE_PRESET_SOURCES
    linked = source != "Custom" and not preset
    if linked:
        prefix = source.removeprefix("Pump")
        for probe_field, pump_field in (("line", "line"), ("transition", "transition"),
                                        ("det_rel", "det_rel"), ("k", "k")):
            st.session_state[key(f"probe_{probe_field}_{label}")] = st.session_state[key(f"{pump_field}_{prefix}")]
        azimuth, ellipticity = preset_ellipse_parameters(
            st.session_state[key(f"k_{prefix}")], st.session_state[key(f"pol_{prefix}")]
        )
        st.session_state[key(f"probe_azimuth_deg_{label}")] = azimuth
        st.session_state[key(f"probe_ellipticity_deg_{label}")] = ellipticity
    with direction_col:
        k_axis = st.selectbox("Direction", ["z", "x", "y"], key=key(f"probe_k_{label}"), disabled=disabled or linked)
    with path_col:
        path = st.number_input("Path length (cm)", min_value=0.0, step=0.1, format="%.3g",
                               key=key(f"probe_path_length_cm_{label}"), disabled=disabled)
    spectral_line_col, transition_col, detuning_col = st.columns(3, gap="xsmall")
    with spectral_line_col:
        line = st.selectbox("Reference line", ["D1", "D2"], key=key(f"probe_line_{label}"), disabled=disabled or linked)
    transition_key = key(f"probe_transition_{label}")
    options = _transition(name, line, transition_key, pressure, coefficients, excited_upper=False)
    with transition_col:
        transition = st.selectbox("F→F′", options, key=transition_key, disabled=disabled or linked)
    with detuning_col:
        relative = st.number_input("Detuning (MHz)", step=10.0, format="%g", key=key(f"probe_det_rel_{label}"), disabled=disabled or linked)
    detuning, selected = _detuning(name, line, transition, relative, pressure, coefficients)
    st.caption(line_center_detuning_caption(line, detuning))
    cols = st.columns(2, gap="xsmall")
    with cols[0]:
        azimuth = st.number_input("Azimuth (°)", min_value=-90.0, max_value=90.0, step=5.0,
                                  key=key(f"probe_azimuth_deg_{label}"), disabled=disabled or linked)
    with cols[1]:
        ellipticity = st.number_input("Ellipticity angle (°)", min_value=-45.0, max_value=45.0, step=5.0,
                                      key=key(f"probe_ellipticity_deg_{label}"), disabled=disabled or linked)
    weak_col, intensity_col = st.columns(2, gap="xsmall", vertical_alignment="bottom")
    with weak_col:
        weak = st.checkbox(
            "Weak", key=key(f"probe_weak_{label}"), disabled=disabled,
            help="A weak probe is a nonperturbing readout and does not alter the atomic state.",
        )
    with intensity_col:
        intensity = st.number_input(
            "I (µW/cm²)", min_value=0.0, step=1.0, format="%.3g",
            key=key(f"probe_intensity_{label}"), disabled=disabled,
            help="Used only when Weak is unchecked.",
        )
    _probe_time_controls(label, disabled)
    if preset:
        st.caption("Preset values are applied when selected. All probe fields remain editable.")
    if weak:
        st.caption("Nonperturbing readout during the configured probe window.")
    else:
        st.caption("Physical probe: absorption, repopulation and light shift act during the probe window.")
    return {"line": line, "detuning_MHz": detuning, "k_axis": k_axis,
            "azimuth_deg": azimuth, "ellipticity_deg": ellipticity, "path_length_cm": path,
            "weak": weak, "intensity": intensity, "source": source,
            "absolute_frequency_MHz": line_center_frequency_MHz(ATOMS[name], line) + detuning,
            "selected_transition": selected, "transition_label": transition}


def perturbing_probe_beam(label, probe):
    """Return the continuous physical beam for a non-weak fixed probe."""
    if probe["weak"] or float(probe["intensity"]) <= 0.0:
        return None
    return {
        "name": f"Probe{label}", "target_label": label,
        "absolute_frequency_MHz": probe["absolute_frequency_MHz"],
        "intensity": probe["intensity"], "k_axis": probe["k_axis"],
        "azimuth_deg": probe["azimuth_deg"], "ellipticity_deg": probe["ellipticity_deg"],
        "selected_transition": probe["selected_transition"],
        "transition_label": probe["transition_label"], "continuous": True,
    }


def render_sidebar():
    """Return independent species, pump, common-cell, and probe configs."""
    initialize_settings()
    _clamp_event_times()
    with st.sidebar:
        controls = st.empty()
        st.header("Atom / cell")
        cols = st.columns(2, gap="xsmall")
        with cols[0]:
            atom_A = st.selectbox("Alkali A", list(ATOMS), key=key("atom_A_name"), on_change=_atom_changed, args=("A",))
        with cols[1]:
            atom_B = st.selectbox("Alkali B", ["None", *ATOMS], key=key("atom_B_name"), on_change=_atom_changed, args=("B",))
        active_B = atom_B not in ("None", atom_A)
        if atom_B == atom_A:
            st.caption("Alkali B matches Alkali A; its pumps and physical effects are inactive.")
        cols = st.columns(2, gap="xsmall")
        with cols[0]:
            mode = st.selectbox("Mixture density model", ["Independent saturated-vapor curves", "Relative concentration"],
                                format_func=lambda value: "Independent" if value.startswith("Independent") else "Mixed",
                                key=key("density_mode"))
        with cols[1]:
            ratio = st.number_input("Molar ratio B/A", min_value=0.0, step=0.1, format="%.3g",
                                    key=key("density_ratio_B_to_A"), disabled=not active_B or mode != "Relative concentration")
        cols = st.columns(2, gap="xsmall")
        with cols[0]:
            pressure = st.number_input("N₂ pressure (Torr)", min_value=0.0, step=10.0, format="%.1f", key=key("n2_pressure_torr"))
        with cols[1]:
            temperature = st.number_input("Temperature (°C)", min_value=-100.0, max_value=500.0, step=1.0,
                                          format="%.1f", key=key("temperature_C_for_table"))
        densities = resolve_alkali_densities(atom_A, atom_B, temperature, mode, ratio)
        cols = st.columns(2, gap="xsmall")
        for label, column in zip(("A", "B"), cols):
            with column:
                st.number_input(rf"$R_{{\mathrm{{ER}},{label}}}$ (s⁻¹)", min_value=0.0, step=1.0, format="%.1f", key=key(f"gamma_ER_{label}"),
                                help="Electron spin randomization rate: relaxes the electron spin toward an unpolarized state while preserving nuclear polarization. Remains active after pump turn-off.")
        cols = st.columns(2, gap="xsmall")
        for label, column in zip(("A", "B"), cols):
            with column:
                st.number_input(rf"$R_{{u,{label}}}$ (s⁻¹)", min_value=0.0, step=1.0, format="%.1f",
                                key=key(f"gamma_u_{label}"),
                                help="Uniform relaxation toward equal ground-state populations; remains active after pump turn-off.")
        theta, phi, field = render_field_controls(PREFIX)
        field_start_col, field_duration_col = st.columns(2, gap="xsmall")
        with field_start_col:
            field_start = st.number_input(
                "Field start time (s)", min_value=0.0,
                max_value=float(st.session_state[key("t2_s")]), step=0.001, format="%.3f",
                key=key("field_start_s"), on_change=_clamp_event_times,
            )
        with field_duration_col:
            st.number_input(
                "Field duration (s)", min_value=0.0,
                max_value=float(st.session_state[key("t2_s")] - field_start), step=0.001,
                format="%.3f", key=key("field_duration_s"),
            )
        common = {"temperature_C": temperature, "n2_pressure_torr": pressure,
                  "static_field_theta_deg": theta, "static_field_phi_deg": phi, "static_field_nT": field}
        atoms = {"A": atom_A, "B": atom_B if active_B else atom_A}
        coefficients = {label: _coefficients(label if label == "A" or active_B else "A") for label in ("A", "B")}
        st.header("Pump configuration")
        beams = []
        for label, tab in zip(("A", "B"), st.tabs(["Alkali A", "Alkali B"])):
            with tab:
                for number, column in zip((1, 2, 3), st.columns(3, gap="xsmall")):
                    with column:
                        beam = _pump(label, number, atoms[label], coefficients[label], pressure,
                                     label == "B" and not active_B)
                        if label == "A" or active_B:
                            beams.append(beam)
        st.header("Probe configuration")
        probes = {}
        for label, tab in zip(("A", "B"), st.tabs(["Probe-A", "Probe-B"])):
            with tab:
                probes[label] = _probe(label, atoms[label], coefficients[label], pressure, label == "B" and not active_B)
        _time_controls()
        with controls.container():
            cols = st.columns([0.2, 0.2, 0.6], gap="xsmall", vertical_alignment="bottom")
            with cols[0]:
                open_file_button(type=["json"], key=_file_widget_key(), on_change=_load_settings,
                                 help="Load independent Time evolution settings.")
            with cols[1]:
                save_placeholder = st.empty()
            with cols[2]:
                st.text_input("Setting name", key=key("condition_name"))
            save_button_with_immediate_download(save_placeholder, prepare=_prepare_settings_download,
                                                mime="application/json", key=key("save"))
            if st.session_state.get(key("load_error")):
                st.error(st.session_state.pop(key("load_error")))
            if st.session_state.get(key("load_success")):
                st.success(st.session_state.pop(key("load_success")))
    species = {}
    for label, density in zip(("A", "B"), densities):
        species[label] = {"label": label, "atom_name": atoms[label], "density_cm3": density,
                          "R_ER": st.session_state[key(f"gamma_ER_{label}")],
                          "R_u": st.session_state[key(f"gamma_u_{label}")],
                          "n2_coeffs": coefficients[label], "q_axis": "z"}
    return species["A"], species["B"] if active_B else None, beams, common, probes
