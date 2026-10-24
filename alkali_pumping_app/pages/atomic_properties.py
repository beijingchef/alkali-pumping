"""Atomic-properties reference page."""

from math import isfinite

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from matplotlib.lines import Line2D

from alkali_pumping_app.physics.atomic_properties import (
    ATOMIC_PROPERTY_DATA,
    BUFFER_GAS_FILL_TEMPERATURE_C,
    BUFFER_GASES,
    ELECTRON_RANDOMIZATION_CROSS_SECTION_CM2,
    SELF_ELECTRON_RANDOMIZATION_CROSS_SECTION_CM2,
    alkali_thermal_properties,
    allowed_hyperfine_F,
    atomic_property_record,
    buffer_gas_properties,
    format_transition_strength_vertical_fraction,
    grotrian_transitions,
    magnetic_sublevels,
)
from alkali_pumping_app.ui.page_state import register_persistent_page_settings


ATOMIC_PROPERTIES_SETTING_KEYS = (
    "atomic_properties_isotope",
    "atomic_properties_tab",
    "atomic_properties_temperature_C",
    "atomic_properties_line",
    *(f"atomic_self_er_sigma_{isotope}" for isotope in ATOMIC_PROPERTY_DATA),
    *(f"atomic_pressure_{gas_name}_torr" for gas_name in BUFFER_GASES),
    *(
        f"atomic_er_sigma_{isotope}_{gas_name}"
        for isotope in ATOMIC_PROPERTY_DATA
        for gas_name in BUFFER_GASES
    ),
    "grotrian_line",
    "grotrian_Fg",
    "grotrian_Fe",
    "grotrian_sigma_minus_v522",
    "grotrian_pi_v522",
    "grotrian_sigma_plus_v522",
    "grotrian_scale_width_v522",
    "grotrian_labels_v522",
)


register_persistent_page_settings(ATOMIC_PROPERTIES_SETTING_KEYS)
for _key, _default in {
    "atomic_properties_line": "D1",
    "grotrian_line": "D1",
    "grotrian_sigma_minus_v522": True,
    "grotrian_pi_v522": True,
    "grotrian_sigma_plus_v522": True,
    "grotrian_scale_width_v522": True,
    "grotrian_labels_v522": True,
}.items():
    st.session_state.setdefault(_key, _default)


def _format_F(value):
    rounded = round(value)
    if abs(value - rounded) < 1e-9:
        return str(int(rounded))
    return f"{int(round(2 * value))}/2"


def _transition_label_position(mg, me, q, line_bottom, line_top):
    """Return a label point that lies exactly on its transition line."""
    label_y = 0.62 if q == 0 else 0.30
    line_fraction = (label_y - line_bottom) / (line_top - line_bottom)
    label_x = mg + line_fraction * (me - mg)
    return label_x, label_y


def _grotrian_figure(isotope, line, Fg, Fe, polarizations, scale_width, labels):
    atom = atomic_property_record(isotope)
    transitions = grotrian_transitions(atom["I"], line, Fg, Fe, polarizations)
    ground_m = magnetic_sublevels(Fg)
    excited_m = magnetic_sublevels(Fe)
    x_min = min(ground_m + excited_m) - 0.65
    x_max = max(ground_m + excited_m) + 0.65

    fig, axis = plt.subplots(figsize=(9.0, 5.4))
    quantum_number_fontsize = 10
    fraction_fontsize = 1.5 * quantum_number_fontsize
    transition_line_bottom = 0.055
    transition_line_top = 0.945
    level_half_width = 0.17
    level_label_gap = 0.10
    for m in ground_m:
        axis.plot([m - level_half_width, m + level_half_width], [0, 0], color="black", lw=2)
    for m in excited_m:
        axis.plot([m - level_half_width, m + level_half_width], [1, 1], color="black", lw=2)

    colors = {-1: "#2575d8", 0: "#1f9d55", 1: "#d64545"}
    max_strength = max((row["strength"] for row in transitions), default=1.0)
    for row in transitions:
        normalized = row["strength"] / max_strength
        linewidth = 0.8 + (5.5 * normalized if scale_width else 1.4)
        axis.plot(
            [row["mg"], row["me"]],
            [transition_line_bottom, transition_line_top],
            color=colors[row["q"]],
            lw=linewidth,
            alpha=0.82,
            solid_capstyle="round",
        )
        if labels:
            label_x, label_y = _transition_label_position(
                row["mg"],
                row["me"],
                row["q"],
                transition_line_bottom,
                transition_line_top,
            )
            axis.text(
                label_x,
                label_y,
                format_transition_strength_vertical_fraction(row["strength"]),
                fontsize=fraction_fontsize,
                ha="center",
                va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1},
            )

    ground_label_x = min(ground_m) - level_half_width - level_label_gap
    excited_label_x = min(excited_m) - level_half_width - level_label_gap
    axis.text(
        ground_label_x,
        0,
        f"F={_format_F(Fg)}",
        ha="right",
        va="center",
        fontsize=10,
    )
    axis.text(
        excited_label_x,
        1,
        f"F′={_format_F(Fe)}",
        ha="right",
        va="center",
        fontsize=10,
    )
    axis.set_xlim(x_min, x_max)
    axis.set_ylim(-0.16, 1.16)
    axis.set_xlabel("Magnetic quantum number m", fontsize=quantum_number_fontsize)
    axis.tick_params(axis="x", labelsize=quantum_number_fontsize)
    axis.set_yticks([])
    axis.set_title(f'{atom["label"]} {line} hyperfine–Zeeman transition strengths')
    axis.spines[["left", "right", "top"]].set_visible(False)
    axis.grid(axis="x", alpha=0.16)
    legend_handles = [
        Line2D([0], [0], color=colors[-1], lw=3, label="σ−  (Δm = −1)"),
        Line2D([0], [0], color=colors[0], lw=3, label="π  (Δm = 0)"),
        Line2D([0], [0], color=colors[1], lw=3, label="σ+  (Δm = +1)"),
    ]
    axis.legend(handles=legend_handles, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout()
    return fig, transitions


def _thermal_tab(isotope, temperature_C):
    atom = atomic_property_record(isotope)
    reference_sigma = SELF_ELECTRON_RANDOMIZATION_CROSS_SECTION_CM2[atom["element"]]
    sigma_key = f"atomic_self_er_sigma_{isotope}"
    st.session_state.setdefault(sigma_key, reference_sigma or 0.0)
    cross_sections = st.columns(2)
    cross_sections[0].metric(
        "Self spin-exchange cross section", f'{atom["spin_exchange_cross_section_cm2"]:.2e} cm²'
    )
    sigma = cross_sections[1].number_input(
        "Self electron-randomization cross section (cm²)",
        min_value=0.0,
        step=1e-18,
        format="%.4e",
        key=sigma_key,
    )
    # An unknown reference with no user-supplied value is unavailable, not zero.
    values = alkali_thermal_properties(
        isotope, temperature_C, None if reference_sigma is None and sigma == 0.0 else sigma
    )
    row1 = st.columns(3)
    row1[0].metric("Saturated-vapor number density", f'{values["density_cm3"]:.4e} cm⁻³')
    row1[1].metric("Self spin-exchange rate", f'{values["spin_exchange_rate_s"]:.4e} s⁻¹')
    self_er_rate = values["self_er_rate_s"]
    row1[2].metric(
        "Self electron-randomization rate",
        "Unavailable" if self_er_rate is None else f"{self_er_rate:.4e} s⁻¹",
    )
    if self_er_rate is None:
        st.info("No built-in self electron-randomization cross section is available for this element. "
                "Enter a measured cross section to calculate its rate.")
    row2 = st.columns(2)
    row2[0].metric("RMS atomic velocity", f'{values["rms_velocity_m_s"]:.3f} m/s')
    row2[1].metric("Mean relative velocity", f'{values["mean_relative_velocity_m_s"]:.3f} m/s')
    st.caption(
        f'Saturated vapor is treated as the {values["phase"]} phase; '
        f'calculated vapor pressure = {values["pressure_torr"]:.4e} Torr. '
        "Both self-collision rates are nσv̄rel, without nuclear slowing down. "
        "Self electron randomization uses a binary spin-destruction cross section, "
        "held constant at the entered value."
    )


_BUFFER_PRESSURE_COLUMN = "Fill pressure at 25 °C (Torr)"
_BUFFER_SIGMA_COLUMN = "Electron-randomization cross section (cm²)"


def _store_buffer_table_edits(isotope, editor_key):
    """Store edited inputs before rebuilding the table's calculated columns."""
    gas_names = list(BUFFER_GASES)
    for row_index, changes in st.session_state[editor_key]["edited_rows"].items():
        gas_name = gas_names[int(row_index)]
        input_keys = {
            _BUFFER_PRESSURE_COLUMN: f"atomic_pressure_{gas_name}_torr",
            _BUFFER_SIGMA_COLUMN: f"atomic_er_sigma_{isotope}_{gas_name}",
        }
        for column, key in input_keys.items():
            if column not in changes:
                continue
            try:
                value = float(changes[column])
            except (TypeError, ValueError):
                value = float("nan")
            if isfinite(value) and value >= 0.0:
                st.session_state[key] = value
            else:
                st.warning("Pressure and cross section must be nonnegative, finite numbers.")
    # Start a fresh edit delta after committing inputs. Old deltas must not
    # overwrite later changes to the pressures shared by different isotopes.
    st.session_state["_atomic_buffer_editor_revision"] += 1


def _buffer_tab(isotope, temperature_C):
    atom = atomic_property_record(isotope)
    element = atom["element"]
    line = st.radio("Optical transition", ["D1", "D2"], horizontal=True, key="atomic_properties_line")
    st.caption(
        f"Fill pressures are specified at {BUFFER_GAS_FILL_TEMPERATURE_C:g} °C in a sealed cell. "
        f"Results use the sidebar temperature ({temperature_C:g} °C). "
        "Operating pressure follows P = Pfill T / Tfill. Broadening is Lorentzian FWHM. "
        "Optical coefficients (MHz/Torr) and entered cross sections are held constant; "
        "their additional temperature dependence is not modeled."
    )

    rows = []
    missing_er = []
    for gas_name, gas in BUFFER_GASES.items():
        pressure_key = f"atomic_pressure_{gas_name}_torr"
        sigma_key = f"atomic_er_sigma_{isotope}_{gas_name}"
        st.session_state.setdefault(pressure_key, 0.0)
        reference_sigma = ELECTRON_RANDOMIZATION_CROSS_SECTION_CM2[element][gas_name]
        st.session_state.setdefault(sigma_key, reference_sigma or 0.0)
        pressure = float(st.session_state[pressure_key])
        sigma = float(st.session_state[sigma_key])
        values = buffer_gas_properties(
            isotope, gas_name, line, temperature_C, pressure, sigma
        )
        unknown_er = reference_sigma is None and sigma == 0.0
        if unknown_er:
            missing_er.append(gas["label"])
        rows.append(
            {
                "Buffer gas": gas["label"],
                _BUFFER_PRESSURE_COLUMN: pressure,
                "Operating pressure (Torr)": values["operating_pressure_torr"],
                "FWHM coeff. (MHz/Torr)": values["fwhm_coefficient_MHz_torr"],
                "Broadening (MHz)": values["broadening_MHz"],
                "Shift coeff. (MHz/Torr)": values["shift_coefficient_MHz_torr"],
                "Shift (MHz)": values["shift_MHz"],
                _BUFFER_SIGMA_COLUMN: sigma,
                "Electron randomization (s⁻¹)": (
                    None if unknown_er else values["electron_randomization_rate_s"]
                ),
            }
        )
    table = pd.DataFrame(rows)
    st.session_state.setdefault("_atomic_buffer_editor_revision", 0)
    editor_key = f"atomic_buffer_editor_{isotope}_{st.session_state['_atomic_buffer_editor_revision']}"
    st.data_editor(
        table,
        hide_index=True,
        width="stretch",
        num_rows="fixed",
        disabled=[column for column in table if column not in (
            _BUFFER_PRESSURE_COLUMN, _BUFFER_SIGMA_COLUMN
        )],
        column_config={
            _BUFFER_PRESSURE_COLUMN: st.column_config.NumberColumn(
                label="P_fill at 25 ℃ (Torr)",
                min_value=0.0,
                required=True,
                format="%.6g",
                help="Editable fill pressure at 25 °C in Torr.",
            ),
            "Operating pressure (Torr)": st.column_config.NumberColumn(
                label="P_T (Torr)",
                help="Operating pressure at the sidebar temperature.",
            ),
            _BUFFER_SIGMA_COLUMN: st.column_config.NumberColumn(
                label="σ_ER (cm²)",
                min_value=0.0,
                required=True,
                format="%.4e",
                help="Editable cross section in cm²; scientific notation is accepted. "
                     "For a gas with no reference value, enter a measured cross section.",
            ),
            "Electron randomization (s⁻¹)": st.column_config.NumberColumn(
                label="Γ_ER (s⁻¹)",
                help="Electron-randomization rate at the sidebar temperature.",
            ),
        },
        key=editor_key,
        on_change=_store_buffer_table_edits,
        args=(isotope, editor_key),
    )
    if missing_er:
        st.info(
            "No built-in ground-state electron-randomization cross section is available "
            f'for {atom["label"]} with {", ".join(missing_er)}. Enter a measured value '
            "in the cross-section column to calculate its rate."
        )
    st.caption(
        "Edit fill pressure and electron-randomization cross section directly in the table. "
        "Calculated columns are read-only. "
        "A blank result means a pressure coefficient or electron-randomization cross section "
        "is unavailable. Enter a measured cross section where needed. "
        "Changing pressure, temperature, or cross section recalculates the results automatically."
    )


def _transition_tab(isotope):
    atom = atomic_property_record(isotope)
    control_panel, graph_panel = st.columns([0.28, 0.72], gap="large")
    with control_panel:
        st.markdown("#### Grotrian controls")
        line = st.selectbox("Fine-structure line", ["D1", "D2"], key="grotrian_line")
        Je = 0.5 if line == "D1" else 1.5
        ground_F = allowed_hyperfine_F(atom["I"], 0.5)
        excited_F = allowed_hyperfine_F(atom["I"], Je)
        if st.session_state.get("grotrian_Fg") not in ground_F:
            st.session_state["grotrian_Fg"] = ground_F[0]
        if st.session_state.get("grotrian_Fe") not in excited_F:
            st.session_state["grotrian_Fe"] = excited_F[0]
        Fe = st.selectbox("Excited F′", excited_F, format_func=_format_F, key="grotrian_Fe")
        Fg = st.selectbox("Ground F", ground_F, format_func=_format_F, key="grotrian_Fg")
        st.markdown("**Polarization**")
        with st.container(horizontal=True, gap="small", wrap=False):
            sigma_minus = st.checkbox("σ−", key="grotrian_sigma_minus_v522")
            pi_enabled = st.checkbox("π", key="grotrian_pi_v522")
            sigma_plus = st.checkbox("σ+", key="grotrian_sigma_plus_v522")
        st.markdown("**Display**")
        scale_width = st.checkbox(
            "Scale line width by strength",
            key="grotrian_scale_width_v522",
        )
        labels = st.checkbox(
            "Show fractional strengths",
            key="grotrian_labels_v522",
        )
    selected = []
    if sigma_minus:
        selected.append(-1)
    if pi_enabled:
        selected.append(0)
    if sigma_plus:
        selected.append(1)
    with graph_panel:
        figure, transitions = _grotrian_figure(
            isotope, line, Fg, Fe, selected, scale_width, labels
        )
        st.pyplot(figure, width="stretch")
        plt.close(figure)
        summed_strength = sum(row["strength"] for row in transitions)
        st.caption(
            f'{len(transitions)} allowed transition(s); summed displayed strength = '
            f'{format_transition_strength_vertical_fraction(summed_strength)}. '
            "Colors and strength normalization follow the supplied Mathematica notebook."
        )


def render_atomic_properties():
    """Render the sidebar isotope selector and three atomic-property tabs."""
    isotope_names = list(ATOMIC_PROPERTY_DATA)
    st.session_state.setdefault("atomic_properties_isotope", "Rb87")
    st.session_state.setdefault(
        "atomic_properties_temperature_C",
        float(st.session_state.get("temperature_C_for_table", 23.0)),
    )
    with st.sidebar:
        st.header("Atomic-properties settings")
        isotope = st.selectbox(
            "Alkali isotope",
            isotope_names,
            format_func=lambda name: ATOMIC_PROPERTY_DATA[name]["label"],
            key="atomic_properties_isotope",
        )
        temperature_C = st.slider(
            "Temperature (°C)",
            min_value=-20.0,
            max_value=300.0,
            step=0.5,
            key="atomic_properties_temperature_C",
        )
    atom = atomic_property_record(isotope)
    st.caption(
        f'Atomic mass {atom["mass_amu"]:.10g} u · nuclear spin I = {_format_F(atom["I"])}'
    )
    thermal_tab, buffer_tab, transition_tab = st.tabs(
        ["Thermal & spin exchange", "Buffer-gas collisions", "Transition strengths"],
        key="atomic_properties_tab",
        on_change="rerun",
    )
    if thermal_tab.open:
        with thermal_tab:
            _thermal_tab(isotope, temperature_C)
    elif buffer_tab.open:
        with buffer_tab:
            _buffer_tab(isotope, temperature_C)
    elif transition_tab.open:
        with transition_tab:
            _transition_tab(isotope)
    st.caption(
        "Reference data: [D. A. Steck alkali data](https://steck.us/alkalidata/); "
        "[Happer, Jau & Walker, *Optically Pumped Atoms*]"
        "(https://doi.org/10.1002/9783527629646); "
        "[Pitz et al. pressure coefficients]"
        "(https://www.sciencedirect.com/science/article/pii/S0022407311004195); "
        "[Liu et al. self spin-destruction cross sections, Table 1]"
        "(https://doi.org/10.1038/s41598-017-06434-2); "
        "and the supplied transition-strength notebook."
    )


def render_atomic_properties_page():
    """Render the complete atomic-properties page."""
    st.title("Atomic properties")
    render_atomic_properties()


if __name__ == "__main__":
    render_atomic_properties_page()
