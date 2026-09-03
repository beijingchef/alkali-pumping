"""Atomic-polarizability analysis page based on Mathur et al. (1970)."""

import json

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from alkali_pumping_app.physics.atomic_polarizability import (
    calculate_atomic_polarizability_sweep,
    default_polarizability_sweep_range_MHz,
)
from alkali_pumping_app.physics.constants import ATOMS, DEFAULT_N2_COEFFS
from alkali_pumping_app.ui.atomic_polarizability_conditions import (
    build_atomic_polarizability_payload,
    clean_atomic_polarizability_condition_name,
    current_atomic_polarizability_values,
    initialize_atomic_polarizability_conditions,
    load_atomic_polarizability_callback,
)
from alkali_pumping_app.ui.downloads import save_button_with_immediate_download
from alkali_pumping_app.ui.exports import (
    atomic_polarizability_export_dataframe,
    dataframe_to_csv_bytes,
)
from alkali_pumping_app.ui.page_state import register_persistent_page_settings
from alkali_pumping_app.ui.uploads import open_file_button


POLARIZABILITY_COMPONENTS = (
    ("alpha_eq", r"$\alpha_{\mathrm{eq}}$"),
    ("alpha_hfs", r"$\alpha_{\mathrm{hfs}}$"),
    ("alpha_gt", r"$\alpha_{\mathrm{gt}}$"),
    ("alpha_br", r"$\alpha_{\mathrm{br}}$"),
)
POLARIZABILITY_Y_AXIS_TITLES = {
    "alpha_eq": "Isotropic polarizability (10⁻¹⁸ cm³)",
    "alpha_hfs": "Hyperfine polarizability (10⁻¹⁸ cm³)",
    "alpha_gt": "Gyrotropic polarizability (10⁻¹⁸ cm³)",
    "alpha_br": "Mathur birefringent coefficient (10⁻¹⁸ cm³)",
}
TRANSITION_MARKER_COLOR = "#606060"
TRANSITION_MARKER_OPACITY = 0.58
LOWER_MANIFOLD_COLOR = "#4C78A8"
UPPER_MANIFOLD_COLOR = "#F58518"
HORIZONTAL_DOMAIN_PADDING_FRACTION = 0.02
ATOMIC_POLARIZABILITY_SETTING_KEYS = (
    "ap_condition_name",
    "ap_atom_name",
    "ap_temperature_C",
    "ap_n2_pressure_torr",
    "ap_D1_width",
    "ap_D1_shift",
    "ap_D2_width",
    "ap_D2_shift",
    "ap_line",
    "ap_lower_MHz",
    "ap_upper_MHz",
    "ap_points",
    "ap_plot_alpha_eq",
    "ap_plot_alpha_hfs",
    "ap_plot_alpha_gt",
    "ap_plot_alpha_br",
)


register_persistent_page_settings(ATOMIC_POLARIZABILITY_SETTING_KEYS)


def _initialize_state():
    initialize_atomic_polarizability_conditions()
    atom_name = st.session_state["ap_atom_name"]
    for line in ("D1", "D2"):
        for coefficient in ("width", "shift"):
            st.session_state.setdefault(
                f"ap_{line}_{coefficient}",
                DEFAULT_N2_COEFFS[atom_name][line][coefficient],
            )
    line = st.session_state["ap_line"]
    pressure_torr = float(st.session_state["ap_n2_pressure_torr"])
    n2_coeffs = {
        coefficient_line: {
            coefficient: float(
                st.session_state[f"ap_{coefficient_line}_{coefficient}"]
            )
            for coefficient in ("width", "shift")
        }
        for coefficient_line in ("D1", "D2")
    }
    range_signature = (
        atom_name,
        line,
        pressure_torr,
        n2_coeffs[line]["shift"],
    )
    if st.session_state.pop("_ap_loaded_preserve_range", False):
        st.session_state["_ap_range_signature"] = range_signature
    elif st.session_state.get("_ap_range_signature") != range_signature:
        lower, upper = default_polarizability_sweep_range_MHz(
            ATOMS[atom_name], line, pressure_torr, n2_coeffs
        )
        st.session_state["ap_lower_MHz"] = lower
        st.session_state["ap_upper_MHz"] = upper
        st.session_state["_ap_range_signature"] = range_signature


def _reset_pressure_coefficients_for_atom():
    atom_name = st.session_state["ap_atom_name"]
    for line in ("D1", "D2"):
        st.session_state[f"ap_{line}_width"] = DEFAULT_N2_COEFFS[atom_name][line]["width"]
        st.session_state[f"ap_{line}_shift"] = DEFAULT_N2_COEFFS[atom_name][line]["shift"]


def _store_condition_name():
    st.session_state["_ap_condition_save_name"] = (
        clean_atomic_polarizability_condition_name(
            st.session_state.get("ap_condition_name")
        )
    )


@st.cache_data(max_entries=32, show_spinner=False)
def _cached_polarizability_sweep(
    atom_name,
    line,
    detunings,
    pressure_torr,
    temperature_C,
    width_coefficient,
    shift_coefficient,
):
    return calculate_atomic_polarizability_sweep(
        atom=ATOMS[atom_name],
        line=line,
        detunings_MHz=np.asarray(detunings, dtype=float),
        n2_pressure_torr=pressure_torr,
        temperature_C=temperature_C,
        n2_width_MHz_per_torr=width_coefficient,
        n2_shift_MHz_per_torr=shift_coefficient,
    )


def polarizability_dataframe(sweep, component):
    """Return one complex polarizability component in long plotting form."""
    detunings = sweep["detunings_MHz"]
    values_by_series = (
        {component.replace("alpha_", "α_"): sweep[component]}
        if component in ("alpha_eq", "alpha_hfs")
        else {f"F={F:g}": values for F, values in sweep[component].items()}
    )
    rows = []
    for series, values in values_by_series.items():
        for part, part_values in (
            ("Real (phase)", np.real(values)),
            ("Imaginary (attenuation)", np.imag(values)),
        ):
            for detuning, value in zip(detunings, part_values):
                rows.append(
                    {
                        "Detuning (MHz)": float(detuning),
                        "Series": series,
                        "Part": part,
                        "Polarizability (10^-18 cm^3)": float(value * 1e18),
                        "Component": component,
                    }
                )
    return pd.DataFrame(rows)


def _render_atomic_polarizability_help(sweep):
    """Explain the four complex polarizability response functions."""
    st.markdown(
        rf"""
$\alpha_{{\mathrm{{eq}}}}$ — **isotropic polarizability.** The scalar,
orientation-independent response of an equilibrium ensemble. It gives the
common refractive and absorptive response that does not distinguish atomic
orientation or alignment.

$\alpha_{{\mathrm{{hfs}}}}$ — **hyperfine polarizability.** The scalar
hyperfine-dependent part of the response. It measures how the two ground
hyperfine manifolds depart from the equilibrium isotropic contribution.

$\alpha_{{\mathrm{{gt}}}}(F,F)$ — **gyrotropic polarizability.** The diagonal
rank-1 response of manifold $F$, associated with atomic orientation and
circular birefringence/dichroism. Its dispersive part produces optical
rotation for an oriented ensemble.

$\alpha_{{\mathrm{{br}}}}(F,F)$ — **birefringent polarizability.** The diagonal
rank-2 response of manifold $F$, associated with atomic alignment and linear
birefringence/dichroism. This page plots Mathur's coefficient itself, which
multiplies Mathur's normalized rank-2 operator. It is therefore not the
coefficient of the app's raw Cartesian
$Q_{{ij}}=(F_iF_j+F_jF_i)/2-\delta_{{ij}}F(F+1)/3$. For the raw-$Q_{{ij}}$
convention the coefficient is $c_F\alpha_{{\mathrm{{br}}}}$, where
$c_F=-\sqrt{{30/[F(F+1)(2F-1)(2F+1)(2F+3)]}}$. The probe-readout calculation
applies this conversion automatically. Off-diagonal hyperfine-coherence terms
are not shown.

For every $\alpha$, $\operatorname{{Re}}\alpha$ is the dispersive phase
response (solid curve) and $\operatorname{{Im}}\alpha$ is the absorptive or
attenuating response (dashed curve). Values are reported in
$10^{{-18}}\,\mathrm{{cm}}^3$.

For the current calculation,
$\Gamma_D={sweep['doppler_fwhm_MHz']:.4g}\,\mathrm{{MHz}}$ is the Doppler FWHM
and $\Gamma_L={sweep['lorentz_fwhm_MHz']:.4g}\,\mathrm{{MHz}}$ is the
Lorentzian FWHM. The response functions follow Mathur, Tang, and Happer,
*Phys. Rev. A* **2**, 648 (1970), Eqs. (41)–(43), (53), and (60).
"""
    )


def _transition_markers(sweep):
    rows = []
    for (Fg, Fe), center in sweep["transition_centers_MHz"].items():
        rows.append(
            {
                "Detuning (MHz)": float(center),
                "Transition": f"F={Fg:g} to F'={Fe:g}",
            }
        )
    return pd.DataFrame(rows)


def _polarizability_chart(dataframe, markers, upper_F, y_axis_title):
    detuning_min = float(dataframe["Detuning (MHz)"].min())
    detuning_max = float(dataframe["Detuning (MHz)"].max())
    detuning_padding = (
        detuning_max - detuning_min
    ) * HORIZONTAL_DOMAIN_PADDING_FRACTION
    detuning_domain = [
        detuning_min - detuning_padding,
        detuning_max + detuning_padding,
    ]
    series = list(dict.fromkeys(dataframe["Series"].tolist()))
    upper_manifold_series = f"F={float(upper_F):g}"
    color = alt.Color("Series:N", title=None)
    if upper_manifold_series in series:
        lower_series = [value for value in series if value != upper_manifold_series]
        color = alt.Color(
            "Series:N",
            title=None,
            scale=alt.Scale(
                domain=[*lower_series, upper_manifold_series],
                range=[
                    *([LOWER_MANIFOLD_COLOR] * len(lower_series)),
                    UPPER_MANIFOLD_COLOR,
                ],
            ),
        )

    line = alt.Chart(dataframe).mark_line(strokeWidth=2, clip=True).encode(
        x=alt.X(
            "Detuning (MHz):Q",
            title="Laser detuning (MHz)",
            scale=alt.Scale(domain=detuning_domain, nice=False, padding=0),
        ),
        y=alt.Y(
            "Polarizability (10^-18 cm^3):Q",
            title=y_axis_title,
            scale=alt.Scale(zero=False),
            axis=alt.Axis(titlePadding=12, labelPadding=4, titleLimit=230),
        ),
        color=color,
        strokeDash=alt.StrokeDash(
            "Part:N",
            title=None,
            scale=alt.Scale(
                domain=["Real (phase)", "Imaginary (attenuation)"],
                range=[[1, 0], [5, 3]],
            ),
        ),
        tooltip=[
            alt.Tooltip("Detuning (MHz):Q", format=".5g"),
            alt.Tooltip("Series:N"),
            alt.Tooltip("Part:N"),
            alt.Tooltip("Polarizability (10^-18 cm^3):Q", format=".6g"),
        ],
    )
    zero = alt.Chart(pd.DataFrame({"zero": [0.0]})).mark_rule(
        color="#808080", opacity=0.45, strokeWidth=1
    ).encode(y="zero:Q")
    transition_rules = alt.Chart(markers).mark_rule(
        color=TRANSITION_MARKER_COLOR,
        opacity=TRANSITION_MARKER_OPACITY,
        strokeDash=[3, 3],
        strokeWidth=1,
        clip=True,
    ).encode(
        x="Detuning (MHz):Q",
        tooltip=["Transition:N", alt.Tooltip("Detuning (MHz):Q", format=".5g")],
    )
    return alt.layer(zero, transition_rules, line).properties(height=390)


_initialize_state()

with st.sidebar:
    st.header("Atomic-polarizability settings")

    load_column, save_column, name_column = st.columns(
        [20, 20, 60], gap="xsmall", vertical_alignment="center"
    )
    with load_column:
        open_file_button(
            type=["json"],
            key="atomic_polarizability_condition_upload",
            on_change=load_atomic_polarizability_callback,
            help="Load atomic-polarizability settings from a JSON file.",
        )
    with save_column:
        save_placeholder = st.empty()
    with name_column:
        st.text_input(
            "Condition name",
            key="ap_condition_name",
            on_change=_store_condition_name,
            help="Filename used when saving these atomic-polarizability settings.",
        )
    if st.session_state.get("_ap_condition_load_message"):
        st.success(st.session_state.pop("_ap_condition_load_message"))
    if st.session_state.get("_ap_condition_load_error"):
        st.error(
            "Could not load condition: "
            + st.session_state.pop("_ap_condition_load_error")
        )

    st.header("Atom / cell")
    atom_name = st.selectbox(
        "Atom",
        list(ATOMS),
        key="ap_atom_name",
        on_change=_reset_pressure_coefficients_for_atom,
    )
    cell_columns = st.columns(2, gap="xsmall")
    with cell_columns[0]:
        temperature_C = st.number_input(
            "Temperature (°C)",
            step=1.0,
            format="%.1f",
            key="ap_temperature_C",
        )
    with cell_columns[1]:
        pressure_torr = st.number_input(
            "N₂ pressure (Torr)",
            min_value=0.0,
            step=10.0,
            format="%.1f",
            key="ap_n2_pressure_torr",
        )

    line = st.segmented_control(
        "Reference line",
        ["D1", "D2"],
        key="ap_line",
    )

    st.header("Sweep")
    range_columns = st.columns(2, gap="xsmall")
    with range_columns[0]:
        lower_MHz = st.number_input(
            "Lower (MHz)",
            step=100.0,
            format="%.6g",
            key="ap_lower_MHz",
        )
    with range_columns[1]:
        upper_MHz = st.number_input(
            "Upper (MHz)",
            step=100.0,
            format="%.6g",
            key="ap_upper_MHz",
        )
    points = st.segmented_control(
        "Points",
        [201, 401, 801],
        key="ap_points",
    )

    st.header("Plots")
    selected_components = []
    plot_columns = st.columns(4, gap="xsmall")
    for column, (component, label) in zip(plot_columns, POLARIZABILITY_COMPONENTS):
        with column:
            if st.checkbox(label, key=f"ap_plot_{component}"):
                selected_components.append(component)

    # Synchronize before preparing the download so its first click uses the
    # currently visible condition name.
    _store_condition_name()
    condition_name = st.session_state["_ap_condition_save_name"]
    payload = build_atomic_polarizability_payload(
        current_atomic_polarizability_values(condition_name)
    )
    save_button_with_immediate_download(
        save_placeholder,
        data=json.dumps(payload, indent=2),
        file_name=f"{condition_name}.json",
        mime="application/json",
        key="save_atomic_polarizability_condition",
    )

title_column, action_column = st.columns([0.68, 0.32], gap="small")
with title_column:
    st.title("Atomic polarizability")
action_placeholder = action_column.empty()

if lower_MHz >= upper_MHz:
    st.error("The lower detuning must be smaller than the upper detuning.")
    st.stop()

detunings = np.linspace(float(lower_MHz), float(upper_MHz), int(points))
width_coefficient = float(st.session_state[f"ap_{line}_width"])
shift_coefficient = float(st.session_state[f"ap_{line}_shift"])
sweep = _cached_polarizability_sweep(
    atom_name,
    line,
    tuple(detunings),
    float(pressure_torr),
    float(temperature_C),
    width_coefficient,
    shift_coefficient,
)
markers = _transition_markers(sweep)

if not selected_components:
    st.info("Select at least one polarizability component in the sidebar.")
    plotted = pd.DataFrame()
else:
    plotted_frames = []
    for component in selected_components:
        frame = polarizability_dataframe(sweep, component)
        plotted_frames.append(frame)
        with st.container(border=True):
            st.altair_chart(
                _polarizability_chart(
                    frame,
                    markers,
                    float(ATOMS[atom_name]["I"]) + 0.5,
                    POLARIZABILITY_Y_AXIS_TITLES[component],
                ),
                width="stretch",
            )

    plotted = pd.concat(plotted_frames, ignore_index=True)

export = (
    atomic_polarizability_export_dataframe(plotted)
    if not plotted.empty
    else pd.DataFrame()
)
with action_placeholder.container():
    download_column, help_column = st.columns([0.78, 0.22], gap="small")
    with download_column:
        with st.container(horizontal_alignment="right"):
            st.download_button(
                "Download CSV",
                dataframe_to_csv_bytes(export),
                file_name=f"{atom_name}_{line}_atomic-polarizability.csv",
                mime="text/csv; charset=utf-8",
                key="download_atomic_polarizability_csv",
                disabled=plotted.empty,
                width="content",
            )
    with help_column:
        with st.popover("❓"):
            _render_atomic_polarizability_help(sweep)
