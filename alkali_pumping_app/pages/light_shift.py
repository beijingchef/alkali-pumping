"""Light-shift analysis page, controls, plots, and plot-data helpers."""

import json

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from alkali_pumping_app.physics.angular_momentum import build_ground_states
from alkali_pumping_app.physics.constants import (
    ATOMS,
    DEFAULT_N2_COEFFS,
    ground_zeeman_shifts_hz,
    upper_larmor_frequency_from_field_nT,
)
from alkali_pumping_app.physics.light_shift import (
    calculate_light_shift_sweep,
    manifold_light_shift_eigenvalues,
)
from alkali_pumping_app.physics.polarization import (
    allowed_polarizations,
    lab_e_field,
    polarization_ellipse_vector,
    preset_ellipse_parameters,
)
from alkali_pumping_app.physics.spectroscopy import (
    doppler_fwhm_MHz,
    hyperfine_transition_choices,
)
from alkali_pumping_app.ui.downloads import save_button_with_immediate_download
from alkali_pumping_app.ui.uploads import open_file_button
from alkali_pumping_app.ui.exports import dataframe_to_csv_bytes
from alkali_pumping_app.ui.light_shift_conditions import (
    LIGHT_SHIFT_DEFAULTS,
    build_light_shift_payload,
    clean_light_shift_condition_name,
    current_light_shift_values,
    initialize_light_shift_conditions,
    load_light_shift_callback,
    state_key,
)
from alkali_pumping_app.ui.text import accent_caption


COMPONENT_LABELS = {
    "scalar": "Scalar shift",
    "vector_per_m": "Vector coefficient",
    "tensor_m0": "Tensor m=0 shift",
}
FICTITIOUS_FIELD_LABEL = "Fictitious magnetic field"
SCALAR_DIFFERENCE_SERIES = "Δν_F"
STATE_COMPONENT_LABELS = {
    "scalar": "Scalar",
    "vector": "Vector",
    "tensor": "Tensor",
    "total": "Total diagonal",
}
TRANSITION_MARKER_COLOR = "#606060"
TRANSITION_MARKER_OPACITY = 0.58


def coefficient_dataframe(sweep, display_detunings_MHz, scale=1.0):
    """Return manifold light-shift coefficients in long plotting form."""
    rows = []
    manifolds = sweep["components_hz_per_uW_cm2"]["coefficients"]
    for manifold in manifolds:
        for key, label in COMPONENT_LABELS.items():
            for detuning, value in zip(display_detunings_MHz, manifold[key]):
                rows.append(
                    {
                        "Detuning (MHz)": float(detuning),
                        "F": f"F={manifold['F']:g}",
                        "Component": label,
                        "Shift": float(scale * value),
                    }
                )
    if len(manifolds) >= 2:
        lower = min(manifolds, key=lambda row: float(row["F"]))
        upper = max(manifolds, key=lambda row: float(row["F"]))
        for detuning, value in zip(
            display_detunings_MHz,
            np.asarray(upper["scalar"]) - np.asarray(lower["scalar"]),
        ):
            rows.append(
                {
                    "Detuning (MHz)": float(detuning),
                    "F": SCALAR_DIFFERENCE_SERIES,
                    "Component": "Scalar shift",
                    "Shift": float(scale * value),
                }
            )
    return pd.DataFrame(rows)


def fictitious_field_dataframe(
    sweep,
    display_detunings_MHz,
    gamma_hz_per_nT_by_F,
    scale=1.0,
):
    """Convert each manifold vector coefficient into a signed field in µG."""
    rows = []
    for manifold in sweep["components_hz_per_uW_cm2"]["coefficients"]:
        F = float(manifold["F"])
        gamma_hz_per_nT = float(gamma_hz_per_nT_by_F[F])
        if np.isclose(gamma_hz_per_nT, 0.0):
            field_values = np.full_like(manifold["vector_per_m"], np.nan, dtype=float)
        else:
            # 1 nT = 10 µG. V_F/gamma_F is the fictitious field in nT.
            field_values = 10.0 * np.asarray(manifold["vector_per_m"]) / gamma_hz_per_nT
        for detuning, value in zip(display_detunings_MHz, field_values):
            rows.append(
                {
                    "Detuning (MHz)": float(detuning),
                    "F": f"F={F:g}",
                    "Component": FICTITIOUS_FIELD_LABEL,
                    "Shift": float(scale * value),
                }
            )
    return pd.DataFrame(rows)


def state_shift_dataframe(sweep, display_detunings_MHz, scale=1.0):
    """Return every scalar/vector/tensor/total diagonal state contribution."""
    components = sweep["components_hz_per_uW_cm2"]
    values_by_component = {
        "scalar": components["scalar"],
        "vector": components["vector"],
        "tensor": components["tensor"],
        "total": sweep["diagonal_hz_per_uW_cm2"],
    }
    rows = []
    for state_index, state in enumerate(sweep["ground_states"]):
        for key, label in STATE_COMPONENT_LABELS.items():
            for point_index, detuning in enumerate(display_detunings_MHz):
                rows.append(
                    {
                        "Detuning (MHz)": float(detuning),
                        "F": float(state["F"]),
                        "m": float(state["m"]),
                        "Component": label,
                        "Shift": float(scale * values_by_component[key][point_index, state_index]),
                    }
                )
    return pd.DataFrame(rows)


def eigenvalue_dataframe(sweep, display_detunings_MHz, intensity_uW_cm2=1.0):
    """Return true within-manifold eigenvalue shifts for arbitrary polarization."""
    eigenvalues = manifold_light_shift_eigenvalues(
        sweep["ground_states"],
        sweep["hamiltonian_hz_per_uW_cm2"],
        sweep["bare_zeeman_hz"],
        intensity_uW_cm2=intensity_uW_cm2,
    )
    rows = []
    for manifold in eigenvalues:
        for branch in range(manifold["shifts"].shape[1]):
            for point_index, detuning in enumerate(display_detunings_MHz):
                rows.append(
                    {
                        "Detuning (MHz)": float(detuning),
                        "F": float(manifold["F"]),
                        "Branch": f"EV{branch + 1}",
                        "Shift": float(manifold["shifts"][point_index, branch]),
                    }
                )
    return pd.DataFrame(rows)


def adjacent_transition_dataframe(sweep, display_detunings_MHz, scale=1.0):
    """Return adjacent-m diagonal transition shifts within every F manifold."""
    diagonal = sweep["diagonal_hz_per_uW_cm2"]
    rows = []
    for F in sorted({float(state["F"]) for state in sweep["ground_states"]}):
        indices = sorted(
            [
                index
                for index, state in enumerate(sweep["ground_states"])
                if np.isclose(float(state["F"]), F)
            ],
            key=lambda index: float(sweep["ground_states"][index]["m"]),
        )
        for lower, upper in zip(indices[:-1], indices[1:]):
            lower_m = float(sweep["ground_states"][lower]["m"])
            upper_m = float(sweep["ground_states"][upper]["m"])
            values = diagonal[:, upper] - diagonal[:, lower]
            for detuning, value in zip(display_detunings_MHz, values):
                rows.append(
                    {
                        "Detuning (MHz)": float(detuning),
                        "F": float(F),
                        "Transition": f"m={lower_m:g} to {upper_m:g}",
                        "Shift": float(scale * value),
                    }
                )
    return pd.DataFrame(rows)


def scattering_dataframe(sweep, display_detunings_MHz, scale=1.0):
    """Return the manifold-average photon-scattering rate."""
    rates = sweep["scattering_s_inv_per_uW_cm2"]
    rows = []
    for F in sorted({float(state["F"]) for state in sweep["ground_states"]}):
        indices = [
            index
            for index, state in enumerate(sweep["ground_states"])
            if np.isclose(float(state["F"]), F)
        ]
        mean_rate = np.mean(rates[:, indices], axis=1)
        for detuning, value in zip(display_detunings_MHz, mean_rate):
            rows.append(
                {
                    "Detuning (MHz)": float(detuning),
                    "F": f"F={F:g}",
                    "Scattering rate": float(scale * value),
                }
            )
    return pd.DataFrame(rows)


def light_shift_export_dataframe(
    plotted,
    *,
    normalization,
    view,
    transition_quantity="Frequency shift",
    scattering=None,
):
    """Return a wide plot-data export with a units row below the headers."""
    detuning_column = "Laser detuning"
    detunings = np.sort(plotted["Detuning (MHz)"].unique())
    numeric = pd.DataFrame({detuning_column: detunings})
    units = {detuning_column: "MHz"}

    if transition_quantity == "Equivalent field" and view == "Transitions":
        shift_unit = (
            "nT/(µW/cm²)" if normalization == "Per intensity" else "nT"
        )
    else:
        shift_unit = (
            "Hz/(µW/cm²)" if normalization == "Per intensity" else "Hz"
        )

    if "Component" in plotted and "m" in plotted:
        series_fields = ["F", "m", "Component"]
        label_for = lambda key: f"F={key[0]:g}, m={key[1]:g} | {key[2]}"
    elif "Component" in plotted:
        series_fields = ["Component", "F"]
        label_for = lambda key: f"{key[0]} | {key[1]}"
    elif "Branch" in plotted:
        series_fields = ["F", "Branch"]
        label_for = lambda key: f"F={key[0]:g} | {key[1]}"
    else:
        series_fields = ["F", "Transition"]
        label_for = lambda key: f"F={key[0]:g} | {key[1]}"

    for key, group in plotted.groupby(series_fields, sort=False, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        label = label_for(key)
        values = group.set_index("Detuning (MHz)")["Shift"].reindex(detunings)
        numeric[label] = values.to_numpy(dtype=float)
        units[label] = (
            str(group["Unit"].iloc[0]) if "Unit" in group else shift_unit
        )

    if scattering is not None and not scattering.empty:
        scatter_unit = (
            "s⁻¹/(µW/cm²)" if normalization == "Per intensity" else "s⁻¹"
        )
        for F, group in scattering.groupby("F", sort=False):
            label = f"Mean scattering rate | {F}"
            values = group.set_index("Detuning (MHz)")["Scattering rate"].reindex(
                detunings
            )
            numeric[label] = values.to_numpy(dtype=float)
            units[label] = scatter_unit

    return pd.concat([pd.DataFrame([units]), numeric], ignore_index=True)


def _render_light_shift_help(
    atom_name,
    line,
    static_field_nT,
    q_axis,
    weights,
    stokes,
    sweep,
    residual_max,
):
    """Explain the quantities reported above the light-shift plots."""
    st.markdown(
        rf"""
**System.** {atom_name} {line}, with static field
$B={static_field_nT:g}\,\mathrm{{nT}}$ along the quantization axis
$\hat{{q}}=\hat{{{q_axis}}}$.

**Spherical polarization weights.** The light polarization is decomposed into
$q=-1,0,+1$ spherical components $\epsilon_q$, normalized by
$\sum_q|\epsilon_q|^2=1$:
$|\epsilon_{{-1}}|^2={weights[-1]:.3f}$,
$|\epsilon_0|^2={weights[0]:.3f}$, and
$|\epsilon_{{+1}}|^2={weights[1]:.3f}$.

**Stokes parameters.** $(s_1,s_2,s_3)=({stokes[0]:+.3f},
{stokes[1]:+.3f},{stokes[2]:+.3f})$ describe horizontal/vertical linear,
diagonal linear, and circular polarization. For normalized transverse fields,
$s_1^2+s_2^2+s_3^2=1$.

**Rank-2 polarization.** $E_{{20}}={sweep['E20']:+.4f}$ is the dimensionless
rank-2 polarization tensor component that weights the tensor light shift.

**Line widths.** The Doppler and Lorentzian full widths at half maximum are
$\Gamma_D={sweep['doppler_fwhm_MHz']:.1f}\,\mathrm{{MHz}}$ and
$\Gamma_L={sweep['lorentz_fwhm_MHz']:.1f}\,\mathrm{{MHz}}$. The latter includes
natural and pressure broadening.

**Numerical check.** The largest residual after the scalar, vector, and tensor
decomposition is
$r_{{\max}}={residual_max:.2e}\,\mathrm{{Hz}}/(\mu\mathrm{{W}}/\mathrm{{cm}}^2)$.
"""
    )


@st.cache_data(max_entries=32, show_spinner=False)
def _cached_sweep(
    atom_name,
    line,
    detunings,
    E_lab,
    k_axis,
    q_axis,
    pressure_torr,
    temperature_C,
    width_coefficient,
    shift_coefficient,
    bare_zeeman_hz,
):
    vector_components = np.array(E_lab, dtype=float)
    vector = vector_components[:3] + 1j * vector_components[3:]
    return calculate_light_shift_sweep(
        atom=ATOMS[atom_name],
        line=line,
        detunings_MHz=np.array(detunings, dtype=float),
        E_lab=vector,
        k_axis=k_axis,
        q_axis=q_axis,
        n2_pressure_torr=pressure_torr,
        temperature_C=temperature_C,
        n2_width_MHz_per_torr=width_coefficient,
        n2_shift_MHz_per_torr=shift_coefficient,
        bare_zeeman_hz=np.array(bare_zeeman_hz, dtype=float),
    )


def _layered_line_chart(
    dataframe,
    markers,
    y_title,
    color_field,
    color_title,
    dash_field=None,
    symlog=False,
    height=325,
    color_scale=None,
    show_color_legend=True,
    color_legend=None,
):
    y_scale = alt.Scale(type="symlog") if symlog else alt.Scale(zero=False)
    color_options = {"title": color_title}
    if color_scale is not None:
        color_options["scale"] = color_scale
    if color_legend is not None:
        color_options["legend"] = color_legend
    elif not show_color_legend:
        color_options["legend"] = None
    encoding = {
        "x": alt.X("Detuning (MHz):Q", title="Laser detuning (MHz)"),
        "y": alt.Y(
            "Shift:Q",
            title=y_title,
            scale=y_scale,
            axis=alt.Axis(
                titlePadding=12,
                labelPadding=4,
                titleLimit=240,
            ),
        ),
        "color": alt.Color(f"{color_field}:N", **color_options),
        "tooltip": [
            alt.Tooltip("Detuning (MHz):Q", format=".4g"),
            alt.Tooltip(f"{color_field}:N"),
            alt.Tooltip("Shift:Q", format=".6g"),
        ],
    }
    if dash_field is not None:
        encoding["strokeDash"] = alt.StrokeDash(f"{dash_field}:N", title=dash_field)
        encoding["tooltip"].insert(2, alt.Tooltip(f"{dash_field}:N"))
    line = alt.Chart(dataframe).mark_line(strokeWidth=2).encode(**encoding)
    zero = alt.Chart(pd.DataFrame({"zero": [0.0]})).mark_rule(
        color="#808080", opacity=0.45, strokeWidth=1
    ).encode(y="zero:Q")
    layers = [zero, line]
    if len(markers):
        layers.insert(1, _transition_marker_rules(markers))
    return alt.layer(*layers).properties(height=height)


def _transition_marker_rules(markers):
    """Return consistently visible transition-center rules for light-shift charts."""
    return alt.Chart(markers).mark_rule(
        color=TRANSITION_MARKER_COLOR,
        opacity=TRANSITION_MARKER_OPACITY,
        strokeDash=[3, 3],
        strokeWidth=1,
    ).encode(
        x="Detuning (MHz):Q",
        tooltip=["Transition:N", alt.Tooltip("Detuning (MHz):Q", format=".4g")],
    )


def _scattering_chart(dataframe, markers, y_title, height=325):
    """Return a scattering-rate chart with untitled F legend and center rules."""
    line = alt.Chart(dataframe).mark_line().encode(
        x=alt.X("Detuning (MHz):Q", title="Laser detuning (MHz)"),
        y=alt.Y("Scattering rate:Q", title=y_title),
        color=alt.Color("F:N", title=None),
        tooltip=[
            alt.Tooltip("Detuning (MHz):Q", format=".4g"),
            "F:N",
            alt.Tooltip("Scattering rate:Q", format=".6g"),
        ],
    )
    layers = [line]
    if len(markers):
        layers.insert(0, _transition_marker_rules(markers))
    return alt.layer(*layers).properties(height=height)


def _scalar_chart_with_difference_legend(
    dataframe,
    markers,
    y_title,
    symlog,
):
    """Render the scalar chart with a simple ΔνF difference label."""
    series = list(dict.fromkeys(dataframe["F"].tolist()))
    palette = ["#4C78A8", "#F58518", "#E45756"]
    color_scale = alt.Scale(domain=series, range=palette[: len(series)])
    legend = alt.Legend(
        title=None,
        labelExpr=(
            f"datum.label === '{SCALAR_DIFFERENCE_SERIES}' "
            "? 'ΔνF' : datum.label"
        ),
    )
    return _layered_line_chart(
        dataframe,
        markers,
        y_title,
        color_field="F",
        color_title=None,
        symlog=symlog,
        height=325,
        color_scale=color_scale,
        color_legend=legend,
    )


def _render_framed_chart(chart):
    """Render a Light shift chart in the same bordered frame as polarizability plots."""
    with st.container(border=True):
        st.altair_chart(chart, width="stretch")


def _render_component_plots(
    dataframe,
    markers,
    frequency_y_title,
    field_y_title,
    symlog,
    show_scalar=True,
):
    for component, y_title in (
        (FICTITIOUS_FIELD_LABEL, field_y_title),
        ("Tensor m=0 shift", frequency_y_title),
    ):
        st.caption(component)
        subset = dataframe[dataframe["Component"] == component]
        chart = _layered_line_chart(
            subset,
            markers,
            y_title,
            color_field="F",
            color_title=None,
            symlog=symlog,
            # At 155 px the axis labels and legend left too little plotting
            # area for narrow dispersive features; fullscreen only appeared
            # to fix the curves because it supplied many more display pixels.
            height=325,
        )
        _render_framed_chart(chart)
    if show_scalar:
        _render_scalar_component_plot(
            dataframe, markers, frequency_y_title, symlog
        )


def _render_scalar_component_plot(dataframe, markers, y_title, symlog):
    """Render the optional scalar component panel."""
    st.caption("Scalar shift")
    subset = dataframe[dataframe["Component"] == "Scalar shift"]
    _render_framed_chart(
        _scalar_chart_with_difference_legend(subset, markers, y_title, symlog)
    )


def _render_state_components_separately(dataframe, markers, y_title, symlog):
    """Render vector and tensor state shifts as separate charts."""
    dataframe = dataframe.copy()
    dataframe["State"] = dataframe.apply(
        lambda row: f"F={row['F']:g}, m={row['m']:g}", axis=1
    )
    for component in ("Vector", "Tensor"):
        st.caption(f"{component} contribution")
        subset = dataframe[dataframe["Component"] == component]
        _render_framed_chart(
            _layered_line_chart(
                subset,
                markers,
                y_title,
                color_field="State",
                color_title="Zeeman state",
                symlog=symlog,
                height=390,
            )
        )


def _automatic_range(atom, line, pressure_torr, n2_coeffs, temperature_C):
    transitions = hyperfine_transition_choices(
        atom, line, pressure_torr, n2_coeffs, allowed_only=True
    )
    centers = np.array([row["detP"] for row in transitions], dtype=float)
    lorentz = (
        float(atom[line]["gamma_nat_MHz"])
        + float(n2_coeffs[line]["width"]) * float(pressure_torr)
    )
    doppler = doppler_fwhm_MHz(atom, line, temperature_C)
    margin = max(250.0, 2.0 * (lorentz + doppler))
    return float(np.floor((centers.min() - margin) / 100.0) * 100.0), float(
        np.ceil((centers.max() + margin) / 100.0) * 100.0
    )


def render_light_shift_explorer(
    result,
    label,
    beams,
    pressure_torr,
    temperature_C,
    condition_name,
):
    """Render the self-contained light-shift controls and plots for one species.

    This renderer intentionally is not a fragment.  It is mounted conditionally
    inside a dynamic expander, and an isolated fragment rerun would lose that
    parent container and clear the plots.  The expensive sweep itself is cached.
    """
    atom_name = result["atom_name"]
    atom = ATOMS[atom_name]
    prefix = f"light_shift_{label}"
    source_options = ["Custom light field", *[beam["name"] for beam in beams]]
    control_col, plot_col = st.columns([0.24, 0.76], gap="small")

    with control_col:
        source = st.selectbox("Light field", source_options, key=f"{prefix}_source")
        source_beam = next((beam for beam in beams if beam["name"] == source), None)
        if source_beam is not None:
            line = source_beam["line"]
            k_axis = source_beam["k_axis"]
            E_lab = lab_e_field(k_axis, source_beam["pol"])
            intensity = float(source_beam["intensity"])
            st.caption(
                accent_caption(
                    f"{line}, k={k_axis}, {source_beam['pol']}, {intensity:g} µW/cm²"
                )
            )
        else:
            line = st.segmented_control(
                "Reference line", ["D1", "D2"], default="D1", key=f"{prefix}_line"
            )
            k_axis = st.selectbox(
                "Beam direction", ["z", "x", "y"], key=f"{prefix}_k_axis"
            )
            polarization_mode = st.segmented_control(
                "Polarization input",
                ["Preset", "Ellipse"],
                default="Preset",
                key=f"{prefix}_polarization_mode",
            )
            if polarization_mode == "Ellipse":
                azimuth = st.slider(
                    "Azimuth psi (degrees)",
                    0.0,
                    180.0,
                    0.0,
                    1.0,
                    key=f"{prefix}_azimuth",
                )
                ellipticity = st.slider(
                    "Ellipticity chi (degrees)",
                    -45.0,
                    45.0,
                    0.0,
                    1.0,
                    key=f"{prefix}_ellipticity",
                    help="-45 degrees is sigma-, 0 is linear, and +45 degrees is sigma+.",
                )
                E_lab = polarization_ellipse_vector(k_axis, azimuth, ellipticity)
            else:
                preset = st.selectbox(
                    "Polarization",
                    allowed_polarizations(k_axis),
                    key=f"{prefix}_preset",
                )
                E_lab = lab_e_field(k_axis, preset)
            intensity = st.number_input(
                "Intensity (µW/cm²)",
                min_value=0.0,
                value=1.0,
                step=1.0,
                format="%.3g",
                key=f"{prefix}_intensity",
            )

        transition_rows = hyperfine_transition_choices(
            atom, line, pressure_torr, result["n2_coeffs"], allowed_only=True
        )
        reference_labels = ["Zero-pressure line center"] + [
            f"F={row['Fg']:g} to F'={row['Fe']:g}" for row in transition_rows
        ]
        reference = st.selectbox(
            "Detuning reference", reference_labels, key=f"{prefix}_reference"
        )
        reference_offset = 0.0
        if reference != reference_labels[0]:
            reference_offset = float(
                transition_rows[reference_labels.index(reference) - 1]["detP"]
            )

        auto_lower, auto_upper = _automatic_range(
            atom, line, pressure_torr, result["n2_coeffs"], temperature_C
        )
        signature = (line, reference)
        if st.session_state.get(f"{prefix}_range_signature") != signature:
            st.session_state[f"{prefix}_lower"] = auto_lower - reference_offset
            st.session_state[f"{prefix}_upper"] = auto_upper - reference_offset
            st.session_state[f"{prefix}_range_signature"] = signature
        range_cols = st.columns(2, gap="xxsmall")
        with range_cols[0]:
            lower = st.number_input(
                "Lower (MHz)", step=100.0, format="%g", key=f"{prefix}_lower"
            )
        with range_cols[1]:
            upper = st.number_input(
                "Upper (MHz)", step=100.0, format="%g", key=f"{prefix}_upper"
            )
        points = st.segmented_control(
            "Sweep points", [201, 401, 801], default=401, key=f"{prefix}_points"
        )
        normalization = st.segmented_control(
            "Shift units",
            ["Per intensity", "Absolute"],
            default="Per intensity",
            key=f"{prefix}_normalization",
        )
        view = st.segmented_control(
            "View",
            ["Components", "Zeeman states", "Eigenvalues", "Transitions"],
            default="Components",
            key=f"{prefix}_view",
        )
        transition_quantity = "Frequency shift"
        if view == "Transitions":
            transition_quantity = st.segmented_control(
                "Transition quantity",
                ["Frequency shift", "Equivalent field"],
                default="Frequency shift",
                key=f"{prefix}_transition_quantity",
            )
        y_scale = st.selectbox(
            "Y scale", ["Linear", "Symmetric log"], key=f"{prefix}_y_scale"
        )
        show_scattering = st.toggle(
            "Show scattering rate", key=f"{prefix}_show_scattering"
        )

    with plot_col:
        if lower >= upper:
            st.error("The upper detuning must be greater than the lower detuning.")
            return
        display_detunings = np.linspace(float(lower), float(upper), int(points))
        absolute_detunings = display_detunings + reference_offset
        vector_array = np.asarray(E_lab, dtype=complex)
        vector_key = tuple(float(value) for value in np.concatenate(
            [np.real(vector_array), np.imag(vector_array)]
        ))
        bare_zeeman = tuple(result["df_pop"]["nu_B"].to_numpy(dtype=float))
        with st.spinner("Calculating light-shift sweep..."):
            sweep = _cached_sweep(
                atom_name,
                line,
                tuple(float(value) for value in absolute_detunings),
                vector_key,
                k_axis,
                result["q_axis"],
                float(pressure_torr),
                float(temperature_C),
                float(result["n2_coeffs"][line]["width"]),
                float(result["n2_coeffs"][line]["shift"]),
                bare_zeeman,
            )

        weights = sweep["spherical_weights"]
        stokes = sweep["stokes"]
        residual_max = max(
            float(np.max(row["residual_max"]))
            for row in sweep["components_hz_per_uW_cm2"]["coefficients"]
        )
        st.caption(
            accent_caption(
                f"{atom_name} {line}; q axis={result['q_axis']}; "
                f"|epsilon_-1|^2={weights[-1]:.3f}, |epsilon_0|^2={weights[0]:.3f}, "
                f"|epsilon_+1|^2={weights[1]:.3f}; "
                f"(s1,s2,s3)=({stokes[0]:+.3f},{stokes[1]:+.3f},{stokes[2]:+.3f}); "
                f"E20={sweep['E20']:+.4f}; "
                f"Doppler FWHM={sweep['doppler_fwhm_MHz']:.1f} MHz; "
                f"Lorentz FWHM={sweep['lorentz_fwhm_MHz']:.1f} MHz; "
                f"max rank-0/1/2 residual={residual_max:.2e} Hz/(µW/cm²)."
            )
        )
        if not result["static_field_aligned"] and abs(result["field_parallel_nT"]) == 0.0:
            st.warning(
                "The shared static field is transverse to this quantization axis. "
                "As elsewhere in the app, transverse magnetic mixing is omitted."
            )
        if not sweep["diagonal_in_selected_basis"] and view != "Eigenvalues":
            st.warning(
                "This polarization contains multiple spherical components. These curves are "
                "first-order diagonal matrix elements in the selected |F,m> basis, not generally "
                "the light-shift eigenvalues. Use the Eigenvalues view when the light defines or mixes the axis."
            )

        scale = 1.0 if normalization == "Per intensity" else float(intensity)
        y_title = (
            "Light shift / intensity (Hz/(µW/cm²))"
            if normalization == "Per intensity"
            else "Light shift (Hz)"
        )
        field_y_title = (
            "B_fic / intensity (µG/(µW/cm²))"
            if normalization == "Per intensity"
            else "B_fic (µG)"
        )
        markers = pd.DataFrame(
            {
                "Detuning (MHz)": [row["detP"] - reference_offset for row in transition_rows],
                "Transition": [
                    f"F={row['Fg']:g} to F'={row['Fe']:g}" for row in transition_rows
                ],
            }
        )
        symlog = y_scale == "Symmetric log"
        manifolds = sorted({float(state["F"]) for state in sweep["ground_states"]})

        if view == "Components":
            coefficients = coefficient_dataframe(sweep, display_detunings, scale)
            upper_gamma_hz_per_nT = upper_larmor_frequency_from_field_nT(atom_name, 1.0)
            gamma_by_F = {
                F: result["bias_info"]["ratio_by_F"][F] * upper_gamma_hz_per_nT
                for F in manifolds
            }
            fields = fictitious_field_dataframe(
                sweep, display_detunings, gamma_by_F, scale
            )
            plotted = pd.concat(
                [
                    coefficients[coefficients["Component"] != "Vector coefficient"],
                    fields,
                ],
                ignore_index=True,
            )
            plotted["Unit"] = np.where(
                plotted["Component"] == FICTITIOUS_FIELD_LABEL,
                "µG/(µW/cm²)" if normalization == "Per intensity" else "µG",
                "Hz/(µW/cm²)" if normalization == "Per intensity" else "Hz",
            )
            _render_component_plots(
                plotted, markers, y_title, field_y_title, symlog
            )
            st.caption(
                "The scalar panel also shows the upper-manifold shift minus the lower-manifold shift. "
                "The vector panel shows B_fic=V_F/γ_F using each manifold's signed gyromagnetic ratio. "
                "The tensor curve is <F,m=0|delta E2|F,m=0>/h. "
                "E20 is dimensionless and is reported separately above; it is not part of the y-axis unit."
            )
        elif view == "Zeeman states":
            selected_F = st.pills(
                "Hyperfine manifolds",
                manifolds,
                default=manifolds,
                selection_mode="multi",
                format_func=lambda value: f"F={value:g}",
                key=f"{prefix}_state_F",
            )
            selected_components = st.pills(
                "State components",
                list(STATE_COMPONENT_LABELS.values()),
                default=["Total diagonal"],
                selection_mode="multi",
                key=f"{prefix}_state_components",
            )
            plotted = state_shift_dataframe(sweep, display_detunings, scale)
            plotted = plotted[
                plotted["F"].isin(selected_F)
                & plotted["Component"].isin(selected_components)
            ]
            for F in selected_F:
                subset = plotted[np.isclose(plotted["F"], float(F))]
                st.caption(f"F={F:g}")
                _render_framed_chart(
                    _layered_line_chart(
                        subset,
                        markers,
                        y_title,
                        color_field="m",
                        color_title="m",
                        dash_field="Component",
                        symlog=symlog,
                    )
                )
        elif view == "Eigenvalues":
            plotted = eigenvalue_dataframe(
                sweep, display_detunings, intensity_uW_cm2=scale
            )
            st.info(
                "These are eigenvalues of each F-manifold light-shift block, including the "
                "parallel Zeeman splitting. Different-F Raman coupling is omitted under the hyperfine secular approximation."
            )
            for F in manifolds:
                subset = plotted[np.isclose(plotted["F"], F)]
                st.caption(f"F={F:g}")
                _render_framed_chart(
                    _layered_line_chart(
                        subset,
                        markers,
                        y_title,
                        color_field="Branch",
                        color_title="Eigenvalues",
                        symlog=symlog,
                    )
                )
        else:
            plotted = adjacent_transition_dataframe(sweep, display_detunings, scale)
            transition_y_title = y_title
            if transition_quantity == "Equivalent field":
                upper_gamma_hz_per_nT = upper_larmor_frequency_from_field_nT(
                    atom_name, 1.0
                )
                for F in manifolds:
                    gamma_hz_per_nT = (
                        result["bias_info"]["ratio_by_F"][F]
                        * upper_gamma_hz_per_nT
                    )
                    mask = np.isclose(plotted["F"], F)
                    plotted.loc[mask, "Shift"] /= gamma_hz_per_nT
                transition_y_title = (
                    "Equivalent light-shift field / intensity (nT/(µW/cm²))"
                    if normalization == "Per intensity"
                    else "Equivalent light-shift field (nT)"
                )
            for F in manifolds:
                subset = plotted[np.isclose(plotted["F"], F)]
                st.caption(f"F={F:g} adjacent-m transition shifts")
                _render_framed_chart(
                    _layered_line_chart(
                        subset,
                        markers,
                        transition_y_title,
                        color_field="Transition",
                        color_title="Transition",
                        symlog=symlog,
                    )
                )

        if show_scattering:
            scatter = scattering_dataframe(sweep, display_detunings, scale)
            scatter_title = (
                "Mean scattering rate / intensity (s⁻¹/(µW/cm²))"
                if normalization == "Per intensity"
                else "Mean scattering rate (s⁻¹)"
            )
            chart = _scattering_chart(scatter, markers, scatter_title, height=273)
            _render_framed_chart(chart)

        st.download_button(
            "Download light-shift CSV",
            dataframe_to_csv_bytes(plotted),
            file_name=f"{condition_name}_alkali-{label}-light-shifts.csv",
            mime="text/csv; charset=utf-8",
            key=f"{prefix}_download",
            width="stretch",
        )


def _light_shift_widget_key(field):
    return f"_ls_widget_{field}"


def _prime_light_shift_control(field, options=None):
    """Restore a page widget from state that survives navigation cleanup."""
    persistent_key = state_key(field)
    value = st.session_state.get(persistent_key, LIGHT_SHIFT_DEFAULTS[field])
    if options is not None and value not in options:
        value = options[0]
        st.session_state[persistent_key] = value
    widget_key = _light_shift_widget_key(field)
    st.session_state[widget_key] = value
    return widget_key


def _prime_light_shift_multiselect(field, options):
    persistent_key = state_key(field)
    value = st.session_state.get(persistent_key)
    if value is None:
        value = list(options)
    else:
        value = [item for item in value if item in options]
    st.session_state[persistent_key] = value
    widget_key = _light_shift_widget_key(field)
    st.session_state[widget_key] = value
    return widget_key


def _store_light_shift_control(field):
    st.session_state[state_key(field)] = st.session_state[_light_shift_widget_key(field)]


def _store_light_shift_name():
    _store_light_shift_control("condition_name")
    st.session_state["_ls_condition_save_name"] = clean_light_shift_condition_name(
        st.session_state[state_key("condition_name")]
    )


def _reset_light_shift_auto_range():
    st.session_state.pop("_ls_range_signature", None)


def _store_light_shift_range_input(field):
    _store_light_shift_control(field)
    _reset_light_shift_auto_range()


def _store_light_shift_atom():
    _store_light_shift_control("atom_name")
    atom_name = st.session_state[state_key("atom_name")]
    for line in ("D1", "D2"):
        st.session_state[state_key(f"{line}_width")] = DEFAULT_N2_COEFFS[atom_name][line]["width"]
        st.session_state[state_key(f"{line}_shift")] = DEFAULT_N2_COEFFS[atom_name][line]["shift"]
    st.session_state[state_key("reference")] = "Zero-pressure line center"
    st.session_state[state_key("state_manifolds")] = None
    st.session_state["_ls_previous_atom"] = atom_name
    _reset_light_shift_auto_range()


def _store_light_shift_beam_axis():
    _store_light_shift_control("k_axis")
    k_axis = st.session_state[state_key("k_axis")]
    if st.session_state.get(state_key("preset")) not in allowed_polarizations(k_axis):
        st.session_state[state_key("preset")] = allowed_polarizations(k_axis)[0]


def _store_light_shift_polarization_mode():
    """Store the mode and preserve the preset polarization when entering Ellipse."""
    previous_mode = st.session_state.get(
        state_key("polarization_mode"), LIGHT_SHIFT_DEFAULTS["polarization_mode"]
    )
    new_mode = st.session_state[_light_shift_widget_key("polarization_mode")]
    if previous_mode == "Preset" and new_mode == "Ellipse":
        k_axis = st.session_state.get(
            state_key("k_axis"), LIGHT_SHIFT_DEFAULTS["k_axis"]
        )
        preset_options = allowed_polarizations(k_axis)
        preset = st.session_state.get(
            state_key("preset"), LIGHT_SHIFT_DEFAULTS["preset"]
        )
        if preset not in preset_options:
            preset = preset_options[0]
            st.session_state[state_key("preset")] = preset
        azimuth, ellipticity = preset_ellipse_parameters(k_axis, preset)
        st.session_state[state_key("azimuth_deg")] = azimuth
        st.session_state[state_key("ellipticity_deg")] = ellipticity
    _store_light_shift_control("polarization_mode")


def _light_shift_n2_coefficients():
    return {
        line: {
            "width": float(st.session_state[state_key(f"{line}_width")]),
            "shift": float(st.session_state[state_key(f"{line}_shift")]),
        }
        for line in ("D1", "D2")
    }


def render_light_shift_page():
    """Render an independent light-shift analysis with all controls in the sidebar."""
    initialize_light_shift_conditions()

    with st.sidebar:
        st.header("Light-shift settings")
        load_column, save_column, name_column = st.columns(
            [20, 20, 60], gap="xsmall", vertical_alignment="center"
        )
        with load_column:
            open_file_button(
                type=["json"],
                key="light_shift_condition_upload",
                on_change=load_light_shift_callback,
                help="Load a light-shift settings JSON file.",
            )
        with save_column:
            save_placeholder = st.empty()
        with name_column:
            st.text_input(
                "Setting name",
                key=_prime_light_shift_control("condition_name"),
                on_change=_store_light_shift_name,
                # label_visibility="collapsed",
                help="Filename used when saving these light-shift settings.",
            )
        if st.session_state.get("_ls_condition_load_message"):
            st.success(st.session_state.pop("_ls_condition_load_message"))
        if st.session_state.get("_ls_condition_load_error"):
            st.error("Could not load settings: " + st.session_state.pop("_ls_condition_load_error"))

        st.header("Atom / cell")
        atom_temperature_columns = st.columns(2, gap="xsmall")
        with atom_temperature_columns[0]:
            atom_name = st.selectbox(
                "Atom",
                list(ATOMS),
                key=_prime_light_shift_control("atom_name", list(ATOMS)),
                on_change=_store_light_shift_atom,
            )
        with atom_temperature_columns[1]:
            temperature_C = st.number_input(
                "Temperature (°C)",
                step=1.0,
                format="%.1f",
                key=_prime_light_shift_control("temperature_C"),
                on_change=_store_light_shift_range_input,
                args=("temperature_C",),
            )
        pressure_torr = st.number_input(
            "N₂ pressure (Torr)",
            min_value=0.0,
            step=10.0,
            format="%.1f",
            key=_prime_light_shift_control("n2_pressure_torr"),
            on_change=_store_light_shift_range_input,
            args=("n2_pressure_torr",),
        )
        field_columns = st.columns(2, gap="xsmall")
        with field_columns[0]:
            q_axis = st.selectbox(
                "Static field direction",
                ["z", "x", "y"],
                key=_prime_light_shift_control("static_field_axis", ["z", "x", "y"]),
                on_change=_store_light_shift_control,
                args=("static_field_axis",),
                help="This direction also defines the quantization axis, including at zero field.",
            )
        with field_columns[1]:
            static_field_nT = st.number_input(
                "Strength (nT)",
                step=1.0,
                format="%g",
                key=_prime_light_shift_control("static_field_nT"),
                on_change=_store_light_shift_control,
                args=("static_field_nT",),
                help="A negative value reverses the selected field direction.",
            )

        with st.expander("N₂ pressure coefficients"):
            for coefficient_line in ("D1", "D2"):
                st.markdown(f"**{coefficient_line}**")
                coefficient_columns = st.columns(2, gap="xsmall")
                with coefficient_columns[0]:
                    st.number_input(
                        "Width (MHz/Torr)",
                        step=0.1,
                        format="%.4g",
                        key=_prime_light_shift_control(f"{coefficient_line}_width"),
                        on_change=_store_light_shift_range_input,
                        args=(f"{coefficient_line}_width",),
                    )
                with coefficient_columns[1]:
                    st.number_input(
                        "Shift (MHz/Torr)",
                        step=0.1,
                        format="%.4g",
                        key=_prime_light_shift_control(f"{coefficient_line}_shift"),
                        on_change=_store_light_shift_range_input,
                        args=(f"{coefficient_line}_shift",),
                    )

        st.header("Light")
        beam_columns = st.columns(2, gap="xsmall")
        with beam_columns[0]:
            k_axis = st.selectbox(
                "Beam direction",
                ["z", "x", "y"],
                key=_prime_light_shift_control("k_axis", ["z", "x", "y"]),
                on_change=_store_light_shift_beam_axis,
            )
        with beam_columns[1]:
            intensity = st.number_input(
                "Intensity (µW/cm²)",
                min_value=0.0,
                step=1.0,
                format="%.3g",
                key=_prime_light_shift_control("intensity_uW_cm2"),
                on_change=_store_light_shift_control,
                args=("intensity_uW_cm2",),
            )
        polarization_mode = st.segmented_control(
            "Polarization input",
            ["Preset", "Ellipse"],
            key=_prime_light_shift_control("polarization_mode", ["Preset", "Ellipse"]),
            on_change=_store_light_shift_polarization_mode,
        )
        if polarization_mode == "Ellipse":
            azimuth = st.slider(
                "Azimuth ψ (degrees)",
                0.0,
                180.0,
                step=1.0,
                key=_prime_light_shift_control("azimuth_deg"),
                on_change=_store_light_shift_control,
                args=("azimuth_deg",),
            )
            ellipticity = st.slider(
                "Ellipticity χ (degrees)",
                -45.0,
                45.0,
                step=1.0,
                key=_prime_light_shift_control("ellipticity_deg"),
                on_change=_store_light_shift_control,
                args=("ellipticity_deg",),
                help="−45° is σ−, 0° is linear, and +45° is σ+.",
            )
            E_lab = polarization_ellipse_vector(k_axis, azimuth, ellipticity)
        else:
            polarization_options = allowed_polarizations(k_axis)
            preset = st.selectbox(
                "Polarization",
                polarization_options,
                key=_prime_light_shift_control("preset", polarization_options),
                on_change=_store_light_shift_control,
                args=("preset",),
            )
            E_lab = lab_e_field(k_axis, preset)
        atom = ATOMS[atom_name]
        n2_coeffs = _light_shift_n2_coefficients()
        reference_columns = st.columns(2, gap="xsmall")
        with reference_columns[0]:
            line = st.segmented_control(
                "Reference line",
                ["D1", "D2"],
                key=_prime_light_shift_control("line", ["D1", "D2"]),
                on_change=_store_light_shift_range_input,
                args=("line",),
            )
        transition_rows = hyperfine_transition_choices(
            atom, line, pressure_torr, n2_coeffs, allowed_only=True
        )
        reference_labels = ["Zero-pressure line center"] + [
            f"F={row['Fg']:g} to F'={row['Fe']:g}" for row in transition_rows
        ]
        with reference_columns[1]:
            reference = st.selectbox(
                "Detuning reference",
                reference_labels,
                key=_prime_light_shift_control("reference", reference_labels),
                on_change=_store_light_shift_range_input,
                args=("reference",),
            )
        reference_offset = 0.0
        if reference != reference_labels[0]:
            reference_offset = float(transition_rows[reference_labels.index(reference) - 1]["detP"])

        auto_lower, auto_upper = _automatic_range(
            atom, line, pressure_torr, n2_coeffs, temperature_C
        )
        range_signature = (
            atom_name,
            line,
            reference,
            float(pressure_torr),
            float(temperature_C),
            float(n2_coeffs[line]["width"]),
            float(n2_coeffs[line]["shift"]),
        )
        if st.session_state.pop("_ls_loaded_preserve_range", False):
            st.session_state["_ls_range_signature"] = range_signature
        elif st.session_state.get("_ls_range_signature") != range_signature:
            st.session_state[state_key("lower_MHz")] = auto_lower - reference_offset
            st.session_state[state_key("upper_MHz")] = auto_upper - reference_offset
            st.session_state["_ls_range_signature"] = range_signature

        detuning_columns = st.columns(2, gap="xsmall")
        with detuning_columns[0]:
            lower = st.number_input(
                "Lower detuning (MHz)",
                step=100.0,
                format="%g",
                key=_prime_light_shift_control("lower_MHz"),
                on_change=_store_light_shift_control,
                args=("lower_MHz",),
            )
        with detuning_columns[1]:
            upper = st.number_input(
                "Upper detuning (MHz)",
                step=100.0,
                format="%g",
                key=_prime_light_shift_control("upper_MHz"),
                on_change=_store_light_shift_control,
                args=("upper_MHz",),
            )
        view_options = [
            "Components",
            "Zeeman states",
            "Zeeman states by component",
            "Eigenvalues",
            "Transitions",
        ]
        view = st.selectbox(
            "View",
            view_options,
            key=_prime_light_shift_control("view", view_options),
            on_change=_store_light_shift_control,
            args=("view",),
        )
        transition_quantity = "Frequency shift"
        if view == "Transitions":
            quantity_options = ["Frequency shift", "Equivalent field"]
            transition_quantity = st.segmented_control(
                "Transition quantity",
                quantity_options,
                key=_prime_light_shift_control("transition_quantity", quantity_options),
                on_change=_store_light_shift_control,
                args=("transition_quantity",),
            )
        y_scale = "Linear"
        st.session_state[state_key("y_scale")] = y_scale
        show_scalar = bool(
            st.session_state.get(
                state_key("show_scalar"), LIGHT_SHIFT_DEFAULTS["show_scalar"]
            )
        )
        if view == "Components":
            show_scalar = st.toggle(
                "Show scalar shift",
                key=_prime_light_shift_control("show_scalar"),
                on_change=_store_light_shift_control,
                args=("show_scalar",),
            )
        show_scattering = st.toggle(
            "Show scattering rate",
            key=_prime_light_shift_control("show_scattering"),
            on_change=_store_light_shift_control,
            args=("show_scattering",),
        )

        ground_states = build_ground_states(atom)
        manifolds = sorted({float(state["F"]) for state in ground_states})
        selected_F = manifolds
        selected_components = ["Total diagonal"]
        if view in ("Zeeman states", "Zeeman states by component"):
            selected_F = st.multiselect(
                "Hyperfine manifolds",
                manifolds,
                format_func=lambda value: f"F={value:g}",
                key=_prime_light_shift_multiselect("state_manifolds", manifolds),
                on_change=_store_light_shift_control,
                args=("state_manifolds",),
            )
            if view == "Zeeman states":
                component_options = list(STATE_COMPONENT_LABELS.values())
                selected_components = st.multiselect(
                    "State components",
                    component_options,
                    key=_prime_light_shift_multiselect("state_components", component_options),
                    on_change=_store_light_shift_control,
                    args=("state_components",),
                )

        sweep_columns = st.columns(2, gap="xsmall")
        with sweep_columns[0]:
            points = st.segmented_control(
                "Sweep points",
                [201, 401, 801],
                key=_prime_light_shift_control("points", [201, 401, 801]),
                on_change=_store_light_shift_control,
                args=("points",),
            )
        with sweep_columns[1]:
            normalization = st.selectbox(
                "Shift units",
                ["Per intensity", "Absolute"],
                key=_prime_light_shift_control(
                    "normalization", ["Per intensity", "Absolute"]
                ),
                on_change=_store_light_shift_control,
                args=("normalization",),
            )

        # Synchronize before preparing the download so its first click uses the
        # currently visible setting name.
        _store_light_shift_name()
        condition_name = st.session_state["_ls_condition_save_name"]
        payload = build_light_shift_payload(current_light_shift_values(condition_name))
        save_button_with_immediate_download(
            save_placeholder,
            data=json.dumps(payload, indent=2),
            file_name=f"{condition_name}.json",
            mime="application/json",
            key="save_light_shift_condition",
        )

    title_column, action_column = st.columns([0.68, 0.32], gap="small")
    with title_column:
        st.title("Light shift")
    action_placeholder = action_column.empty()
    if lower >= upper:
        st.error("The upper detuning must be greater than the lower detuning.")
        return

    display_detunings = np.linspace(float(lower), float(upper), int(points))
    absolute_detunings = display_detunings + reference_offset
    vector_array = np.asarray(E_lab, dtype=complex)
    vector_key = tuple(
        float(value)
        for value in np.concatenate([np.real(vector_array), np.imag(vector_array)])
    )
    upper_larmor_hz = upper_larmor_frequency_from_field_nT(atom_name, static_field_nT)
    bare_zeeman, bias_info = ground_zeeman_shifts_hz(
        atom_name, atom, ground_states, upper_larmor_hz
    )
    with st.spinner("Calculating light-shift sweep..."):
        sweep = _cached_sweep(
            atom_name,
            line,
            tuple(float(value) for value in absolute_detunings),
            vector_key,
            k_axis,
            q_axis,
            float(pressure_torr),
            float(temperature_C),
            float(n2_coeffs[line]["width"]),
            float(n2_coeffs[line]["shift"]),
            tuple(float(value) for value in bare_zeeman),
        )

    weights = sweep["spherical_weights"]
    stokes = sweep["stokes"]
    residual_max = max(
        float(np.max(row["residual_max"]))
        for row in sweep["components_hz_per_uW_cm2"]["coefficients"]
    )
    if not sweep["diagonal_in_selected_basis"] and view != "Eigenvalues":
        st.warning(
            "This polarization contains multiple spherical components. These curves are "
            "first-order diagonal matrix elements in the selected |F,m⟩ basis. Use "
            "Eigenvalues to display the light-shift eigenvalues when the light mixes that basis."
        )

    scale = 1.0 if normalization == "Per intensity" else float(intensity)
    y_title = (
        "Light shift / intensity (Hz/(µW/cm²))"
        if normalization == "Per intensity"
        else "Light shift (Hz)"
    )
    field_y_title = (
        "B_fic / intensity (µG/(µW/cm²))"
        if normalization == "Per intensity"
        else "B_fic (µG)"
    )
    markers = pd.DataFrame(
        {
            "Detuning (MHz)": [row["detP"] - reference_offset for row in transition_rows],
            "Transition": [
                f"F={row['Fg']:g} to F'={row['Fe']:g}" for row in transition_rows
            ],
        }
    )
    symlog = y_scale == "Symmetric log"

    if view == "Components":
        coefficients = coefficient_dataframe(sweep, display_detunings, scale)
        upper_gamma_hz_per_nT = upper_larmor_frequency_from_field_nT(atom_name, 1.0)
        gamma_by_F = {
            F: bias_info["ratio_by_F"][F] * upper_gamma_hz_per_nT
            for F in manifolds
        }
        fields = fictitious_field_dataframe(
            sweep, display_detunings, gamma_by_F, scale
        )
        plotted = pd.concat(
            [
                coefficients[coefficients["Component"] != "Vector coefficient"],
                fields,
            ],
            ignore_index=True,
        )
        plotted["Unit"] = np.where(
            plotted["Component"] == FICTITIOUS_FIELD_LABEL,
            "µG/(µW/cm²)" if normalization == "Per intensity" else "µG",
            "Hz/(µW/cm²)" if normalization == "Per intensity" else "Hz",
        )
        _render_component_plots(
            plotted,
            markers,
            y_title,
            field_y_title,
            symlog,
            show_scalar=False,
        )
        st.caption(
            "The scalar panel also shows the upper-manifold shift minus the lower-manifold shift. "
            "The vector panel shows B_fic=V_F/γ_F using each manifold's signed gyromagnetic ratio. "
            "The tensor curve is ⟨F,m=0|δE²|F,m=0⟩/h. E₂₀ is dimensionless and is reported separately above."
        )
    elif view == "Zeeman states":
        plotted = state_shift_dataframe(sweep, display_detunings, scale)
        plotted = plotted[
            plotted["F"].isin(selected_F)
            & plotted["Component"].isin(selected_components)
        ]
        if not selected_F or not selected_components:
            st.info("Select at least one hyperfine manifold and one state component in the sidebar.")
        for F in selected_F:
            subset = plotted[np.isclose(plotted["F"], float(F))]
            st.caption(f"F={F:g}")
            _render_framed_chart(
                _layered_line_chart(
                    subset,
                    markers,
                    y_title,
                    color_field="m",
                    color_title="m",
                    dash_field="Component",
                    symlog=symlog,
                )
            )
    elif view == "Zeeman states by component":
        plotted = state_shift_dataframe(sweep, display_detunings, scale)
        plotted = plotted[
            plotted["F"].isin(selected_F)
            & plotted["Component"].isin(["Vector", "Tensor"])
        ]
        if not selected_F:
            st.info("Select at least one hyperfine manifold in the sidebar.")
        else:
            _render_state_components_separately(
                plotted, markers, y_title, symlog
            )
    elif view == "Eigenvalues":
        plotted = eigenvalue_dataframe(sweep, display_detunings, intensity_uW_cm2=scale)
        st.info(
            "Each curve is an eigenvalue of one F-manifold light-shift block, including "
            "parallel Zeeman splitting. Different-F Raman coupling is omitted."
        )
        for F in manifolds:
            subset = plotted[np.isclose(plotted["F"], F)]
            st.caption(f"F={F:g}")
            _render_framed_chart(
                _layered_line_chart(
                    subset,
                    markers,
                    y_title,
                    color_field="Branch",
                    color_title="Eigenvalues",
                    symlog=symlog,
                )
            )
    else:
        plotted = adjacent_transition_dataframe(sweep, display_detunings, scale)
        transition_y_title = y_title
        if transition_quantity == "Equivalent field":
            upper_gamma_hz_per_nT = upper_larmor_frequency_from_field_nT(atom_name, 1.0)
            for F in manifolds:
                gamma_hz_per_nT = bias_info["ratio_by_F"][F] * upper_gamma_hz_per_nT
                mask = np.isclose(plotted["F"], F)
                plotted.loc[mask, "Shift"] /= gamma_hz_per_nT
            transition_y_title = (
                "Equivalent light-shift field / intensity (nT/(µW/cm²))"
                if normalization == "Per intensity"
                else "Equivalent light-shift field (nT)"
            )
        for F in manifolds:
            subset = plotted[np.isclose(plotted["F"], F)]
            st.caption(f"F={F:g} adjacent-m transition shifts")
            _render_framed_chart(
                _layered_line_chart(
                    subset,
                    markers,
                    transition_y_title,
                    color_field="Transition",
                    color_title="Transition",
                    symlog=symlog,
                )
            )

    scatter = None
    if show_scattering:
        scatter = scattering_dataframe(sweep, display_detunings, scale)
        scatter_title = (
            "Mean scattering rate / intensity (s⁻¹/(µW/cm²))"
            if normalization == "Per intensity"
            else "Mean scattering rate (s⁻¹)"
        )
        chart = _scattering_chart(scatter, markers, scatter_title)
        _render_framed_chart(chart)

    if view == "Components" and show_scalar:
        _render_scalar_component_plot(plotted, markers, y_title, symlog)

    export = light_shift_export_dataframe(
        plotted,
        normalization=normalization,
        view=view,
        transition_quantity=transition_quantity,
        scattering=scatter,
    )
    with action_placeholder.container():
        download_column, help_column = st.columns([0.78, 0.22], gap="small")
        with download_column:
            st.download_button(
                "Download CSV",
                dataframe_to_csv_bytes(export),
                file_name=f"{condition_name}_{atom_name}_{line}_light-shifts.csv",
                mime="text/csv; charset=utf-8",
                key="download_light_shift_page_csv",
                width="stretch",
            )
        with help_column:
            with st.popover("❓"):
                _render_light_shift_help(
                    atom_name,
                    line,
                    static_field_nT,
                    q_axis,
                    weights,
                    stokes,
                    sweep,
                    residual_max,
                )


if __name__ == "__main__":
    render_light_shift_page()
