"""Coherent pumping transients and pump-on steady-state probe spectra."""

import altair as alt
from hashlib import sha256
import pickle
import numpy as np
import streamlit as st

from alkali_pumping_app.physics.constants import ATOMS
from alkali_pumping_app.ui.time_evolution_cache import (
    evolution_cached as _evolution_cached, display_cached as _display_cached,
    steady_cached as _steady_cached, spectrum_cached as _spectrum_cached, probe_cached as _probe_cached,
)
from alkali_pumping_app.ui.page_state import register_persistent_page_settings
from alkali_pumping_app.ui.exports import dataframe_to_csv_bytes
from alkali_pumping_app.ui.time_evolution_conditions import PLOTTED_MOMENTS, PLOTTED_PROBE_SIGNALS
from alkali_pumping_app.ui.time_evolution_data import (
    probe_scan_export_dataframe, section_export_dataframe, series_frame, zeeman_plot_frame,
)
from alkali_pumping_app.ui.time_evolution_settings import (
    key, render_sidebar,
)


st.set_page_config(initial_sidebar_state=420)
species_A, species_B, beams, common, probes = render_sidebar()
t1 = st.session_state[key("t1_s")]
t2 = st.session_state[key("t2_s")]
SECTION_CONTROL_WIDTH = 56

# Native bordered containers have symmetric padding. Limit this exception to
# plot panels so their lower margin can be reduced independently.
st.html("""
<style>
[class*="st-key-te_plot_panel_"] { padding-bottom: 0.25rem !important; }
[class*="st-key-te_axis_label_"] {
    position: relative;
    height: 270px;
    min-height: 270px;
    overflow: visible !important;
}
[class*="st-key-te_axis_label_"] [data-testid="stElementContainer"],
[class*="st-key-te_axis_label_"] [data-testid="stMarkdown"] {
    position: static;
}
[class*="st-key-te_axis_label_"] [data-testid="stMarkdownContainer"] {
    position: absolute;
    left: 50%;
    top: 50%;
    width: max-content;
    white-space: nowrap;
    transform: translate(-50%, -50%) rotate(-90deg);
    font-size: 12px;
}
[class*="st-key-te_axis_label_"] p { margin: 0; }
</style>
""")

with st.container(horizontal=True, vertical_alignment="center"):
    st.title("Time evolution", width="stretch")
    with st.popover("❓", help="About Time evolution", key=key("page_help")):
        st.markdown(r"""
Initially unpolarized atoms have equal probability in every ground Zeeman
state. All configured pumps are on from $t=0$ to the pump turn-off time,
then off through the total time. Electron randomization, uniform relaxation,
and self/cross spin exchange remain active throughout. RF is off.

Time-domain graphs span zero to the total time. The dashed vertical line marks pump
turn-off; populations are continuous there, although their derivatives change.

The density-matrix model retains all Zeeman coherences within each ground
hyperfine level, with magnetic precession and pump AC-Stark evolution. Each
species has independent population and moment axes. The Zeeman selector is
in its section's third column; the moment selector is before its checkboxes.
Neither needs to align with the magnetic field.
The Zeeman axis defines populations and exported coherences. Electron
polarization and $Q_{qq}$ use the atomic-moment axis.
$F_x,F_y,F_z$ use laboratory axes. The alignment frame
$(\hat{\mathbf u},\hat{\mathbf v},\hat{\mathbf q})$ is right-handed:
$(y,z,x)$ for $q=x$, $(z,x,y)$ for $q=y$, and $(x,y,z)$ for $q=z$.
Inter-hyperfine coherences are averaged out. Magnetic-field shifts are omitted
from optical detunings. The continuous weak probe reads the full density matrix
without changing the evolution.

Each time-domain plot section has a narrow third column with one help popup
above one download icon. The scan has help to the left of download in its top row.
Section CSVs contain all displayed curves in the plotted units,
with timing and species metadata.

Probe signal scanning separately converges the pump-on steady state and scans
the selected probe D line. Its detuning range slider is at the top left. The
three frequency plots show rotation, ellipticity and absorption one per row;
their light dashed lines mark allowed hyperfine transitions. Scan results do
not depend on the pump turn-off or total-time sliders.

Zeeman populations, Atomic moments and Probe signal scanning start collapsed.
Only open sections create plots and controls or refresh them after sidebar
changes. Closing all time-domain sections defers transient integration until
one is opened. Turn off Auto update to edit several physical conditions, then
press Apply once. Pending curves and downloads are hidden until applied.
Atomic moments default to hyperfine polarization and $\langle F_z\rangle$.
Probe signal checkboxes select absorption, optical rotation and ellipticity;
absorption and optical rotation start checked. Each species has its own probe
selection. The scan always shows all three signals and calculates its steady
state only while open.
Selected plots start in default checkbox order. Newly checked plots append at
the end; unchecking removes a plot, and checking it again appends it at the end.
Plot order is retained in saved settings and section CSV columns.
""")


def _store_scan_range(label):
    widget_key = "_scan_widget_" + key(f"probe_scan_range_MHz_{label}")
    lower, upper = st.session_state[widget_key]
    low, high = key(f"probe_scan_min_MHz_{label}"), key(f"probe_scan_max_MHz_{label}")
    # Native range handles may meet; keep a nonzero interval for the spectrum.
    if upper - lower < 10.:
        if lower != st.session_state[low]:
            lower = max(-15000., upper - 10.)
            upper = max(upper, lower + 10.)
        else:
            upper = min(15000., lower + 10.)
            lower = min(lower, upper - 10.)
    st.session_state[low], st.session_state[high] = float(lower), float(upper)
    st.session_state[widget_key] = (float(lower), float(upper))


def _scan_range_slider(label):
    widget_key = "_scan_widget_" + key(f"probe_scan_range_MHz_{label}")
    st.session_state[widget_key] = (
        st.session_state[key(f"probe_scan_min_MHz_{label}")],
        st.session_state[key(f"probe_scan_max_MHz_{label}")],
    )
    st.slider("Detuning range (MHz)", min_value=-15000., max_value=15000.,
              step=10., key=widget_key, on_change=_store_scan_range, args=(label,),
              label_visibility="collapsed")


def _moment_expression(name, q):
    u, v = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}[q]
    if name == "Q_uu-Q_vv":
        return rf"Q_{{{u}{u}}}-Q_{{{v}{v}}}"
    if name == "Q_qq":
        return rf"Q_{{{q}{q}}}"
    if name == "Q_uv":
        return rf"Q_{{{u}{v}}}"
    if name.startswith(("F_", "Q_")):
        symbol, indices = name.split("_", 1)
        return f"{symbol}_{{{indices}}}"
    return None


def _store_plot_choice(selection_name, name, option_key):
    # Update only the toggled choice so the remaining plots keep their order.
    selected = list(st.session_state[key(selection_name)])
    if st.session_state[option_key]:
        if name not in selected:
            selected.append(name)
    else:
        selected = [item for item in selected if item != name]
    st.session_state[key(selection_name)] = selected


def _section_controls(plots, *, species, label, section_id, title, export_frame=None,
                      help_text=None, horizontal=False):
    export_species = dict(species, q_axis=species["moment_axis"]) if section_id == "atomic_moments" else species
    export = export_frame if export_frame is not None else section_export_dataframe(
        plots, species=export_species, label=label, section=title, t1_s=t1, t2_s=t2,
    )
    # Streamlit permits disabling wrapping only for horizontal containers.
    with st.container(horizontal=horizontal, wrap=not horizontal, gap="xsmall",
                      horizontal_alignment="right" if horizontal else "center"):
        with st.popover("❓", help=f"About {title}", key=key(f"help_section_{label}_{section_id}")):
            for plot in plots:
                st.markdown(f"${plot['math_title']}$" if plot["math_title"] else f"**{plot['title']}**")
                st.markdown(plot["explanation"])
            if not plots and section_id == "atomic_moments":
                st.markdown("Select moment checkboxes to display and download atomic-moment curves.")
            if not plots and section_id == "probe_signal":
                st.markdown("Select probe checkboxes to display and download probe-signal curves.")
            st.markdown(help_text or (
                "The dashed vertical line marks pump turn-off. The CSV contains all displayed "
                "curves in this section, one column per curve, with plot names and plotted units "
                "in the column labels, time in seconds, and timing/species/axis metadata."
            ))
        st.download_button(
            ":material/download:", dataframe_to_csv_bytes(export),
            file_name=f"time_evolution_{label}_{species['atom_name']}_{section_id}.csv",
            mime="text/csv; charset=utf-8", key=key(f"download_section_{label}_{section_id}"),
            help=f"Download {title} CSV", disabled=not plots, on_click="ignore", width=SECTION_CONTROL_WIDTH,
        )


def _plot(frame, y_label, title, *, species, label, plot_id, value_unit, explanation, section_plots,
          show_legend=True, math_title=None, x_label="Time (s)", x_domain=None,
          show_switch=True, reference_lines=None):
    section_plots.append(dict(frame=frame, y_label=y_label, title=title,
                              value_unit=value_unit, explanation=explanation, math_title=math_title))
    with st.container(border=True, key=key(f"plot_panel_{label}_{plot_id}")):
        st.markdown(f"${math_title}$" if math_title else f"**{title}**")
        if not np.isfinite(frame[y_label]).any():
            st.info(f"{title} is undefined for this input polarization or extinguished output.")
            return
        # Choose precision once for the graph's complete vertical extent,
        # including all displayed curves and zero, rather than per sample.
        values = frame[y_label].to_numpy(dtype=float)
        small_range = np.max(np.abs(values[np.isfinite(values)])) < 0.01
        quantity_axis = alt.Axis(
            format=".2e", labelExpr="datum.value === 0 ? '0' : format(datum.value, '.2e')",
        ) if small_range else alt.Axis()
        quantity_format = ".2e" if small_range else ".6g"
        # Vega axis titles are plain text. Reserve a narrow native Streamlit
        # label beside the chart so KaTeX can typeset subscripts and powers.
        axis_math = {
            "⟨S·I⟩/ℏ²": r"\langle\mathbf S\cdot\mathbf I\rangle/\hbar^2",
            "2⟨S_q⟩/ℏ": r"2\langle S_q\rangle/\hbar",
            **{f"⟨F_{axis}⟩/ℏ": rf"\langle F_{{{axis}}}\rangle/\hbar" for axis in "xyz"},
            "Absorption 1 − I_out/I_in": r"1-I_{\mathrm{out}}/I_{\mathrm{in}}",
        }.get(y_label)
        if y_label.startswith("⟨Q_") and math_title:
            axis_math = math_title + r"/\hbar^2"
        curve = alt.Chart(frame).mark_line().encode(
            x=alt.X(f"{x_label}:Q", scale=alt.Scale(domain=x_domain if x_domain is not None else [0, t2])),
            y=alt.Y(f"{y_label}:Q", title=None if axis_math else y_label,
                    scale=alt.Scale(zero=True), axis=quantity_axis),
            color=alt.Color("Series:N", sort=list(frame["Series"].drop_duplicates()), legend=alt.Legend(
                orient="top-right", direction="vertical", title=None, padding=5, offset=5,
            ) if show_legend else None),
            tooltip=[alt.Tooltip(f"{x_label}:Q", format=".6g"), "Series:N",
                     alt.Tooltip(f"{y_label}:Q", format=quantity_format)],
        )
        chart = alt.layer(curve)
        if show_switch:
            switch = alt.Chart({"values": [{"switch": t1}]}).mark_rule(
                color="#666666", strokeDash=[5, 4]
            ).encode(x=alt.X("switch:Q", title=x_label))
            chart += switch
        if reference_lines:
            markers = alt.Chart({"values": reference_lines}).mark_rule(
                color="#777777", opacity=.3, strokeDash=[2, 3],
            ).encode(x=alt.X("center:Q", title=x_label),
                     tooltip=["transition:N", alt.Tooltip("center:Q", format=".6g", title="Detuning (MHz)")])
            chart += markers
        chart = chart.properties(
            description=title, height=260,
            padding={"left": 5, "right": 5, "top": 5, "bottom": 5},
        )
        if axis_math:
            with st.container(horizontal=True, vertical_alignment="top", gap=None, wrap=False):
                with st.container(width=36, border=False,
                                  key=key(f"axis_label_{label}_{plot_id}")):
                    st.markdown(f"${axis_math}$")
                st.altair_chart(chart, width="stretch")
        else:
            st.altair_chart(chart, width="stretch")


# Execution controls are session UI preferences, independent of saved conditions.
st.session_state.setdefault(key("auto_update"), True)
register_persistent_page_settings([key("auto_update")])
physical_beams = [beam for beam in beams if beam["intensity"] > 0]
physical_signature = sha256(pickle.dumps((species_A, species_B, physical_beams, common))).hexdigest()
transient_signature = sha256(pickle.dumps((physical_signature, t1, t2))).hexdigest()
with st.bottom:
    with st.container(horizontal=True, wrap=False, vertical_alignment="center"):
        auto_update = st.checkbox("Auto update", key=key("auto_update"), width="content")
        apply_changes = st.button("Apply", key=key("apply_conditions"), disabled=auto_update, width="content")
if auto_update or apply_changes:
    st.session_state[key("applied_transient")] = transient_signature
    st.session_state[key("applied_steady")] = physical_signature
time_ready = st.session_state.get(key("applied_transient")) == transient_signature
scan_ready = st.session_state.get(key("applied_steady")) == physical_signature
if not time_ready or not scan_ready:
    st.info("Pending condition changes. Press Apply to update the opened sections.")

# Reserve all section headers before any expensive calculation. One trajectory
# is retrieved per rerun and shared by both species, only when a time section is open.
evolution = None
labels = ["A", "B"] if species_B is not None else ["A"]
configs = {"A": species_A, "B": species_B}
tabs = st.tabs([f"Alkali {label} · {configs[label]['atom_name']}" for label in labels])
for label, tab in zip(labels, tabs):
    with tab:
        config = configs[label]
        species = dict(atom_name=config["atom_name"],
                       q_axis=st.session_state[key(f"zeeman_axis_{label}")],
                       zeeman_axis=st.session_state[key(f"zeeman_axis_{label}")],
                       moment_axis=st.session_state[key(f"moment_axis_{label}")])
        upper_F = ATOMS[config["atom_name"]]["I"] + .5
        Fs = [upper_F - 1, upper_F]
        manifold_labels = {"Fa": f"F = {Fs[-1]:g}", "Fb": f"F = {Fs[0]:g}"}
        zeeman_section = st.expander(
            "Zeeman populations", expanded=False, key=key(f"zeeman_populations_open_{label}"),
            on_change="rerun",
        )
        moment_section = st.expander(
            "Atomic moments", expanded=False, key=key(f"atomic_moments_open_{label}"), on_change="rerun",
        )
        probe_section = st.expander(
            "Probe signal", expanded=True, key=key(f"probe_signal_open_{label}"), on_change="rerun",
        )
        scan_section = st.expander(
            "Probe signal scanning", expanded=False, key=key(f"probe_scan_open_{label}"), on_change="rerun",
        )
        if time_ready and any(section.open for section in (zeeman_section, moment_section, probe_section)):
            if evolution is None:
                with st.spinner("Calculating pumping and relaxation…"):
                    try:
                        evolution = _evolution_cached(species_A, species_B, physical_beams, common, t1, t2)
                    except (ValueError, RuntimeError) as exc:
                        st.error(str(exc))
                        st.stop()
                if evolution["sampling_limited"]:
                    st.warning("The output is limited to 20,001 samples. Fast Larmor oscillations may be "
                               "undersampled; reduce total time or magnetic field to resolve them.")
            physical_species = evolution[label]
            times = evolution["time_s"]
            if zeeman_section.open or moment_section.open:
                species = _display_cached(physical_species, species["zeeman_axis"], species["moment_axis"])
        if not time_ready:
            for section in (zeeman_section, moment_section, probe_section):
                if section.open:
                    with section:
                        st.info("Press Apply to calculate these curves for the current conditions.")
        if not scan_ready and scan_section.open:
            with scan_section:
                st.info("Press Apply to calculate the scan for the current conditions.")
        if zeeman_section.open and time_ready:
            with zeeman_section:
                zeeman_columns = st.columns([1, 1, 0.15], gap="small")
                zeeman_plots = []
                for name, F, column in zip(("Fa", "Fb"), reversed(Fs), zeeman_columns[:2]):
                    with column:
                        _plot(zeeman_plot_frame(times, species, name), "Population", f"F = {F:g}",
                              species=species, label=label, plot_id=f"zeeman_{name}", value_unit="1",
                              section_plots=zeeman_plots, math_title=rf"F={F:g}",
                              explanation=rf"Individual Zeeman populations in $F={F:g}$ along the "
                              "selected quantization axis. Values are probabilities per total atom, "
                              "without renormalizing within the manifold.")
                with zeeman_columns[2]:
                    _section_controls(zeeman_plots, species=species, label=label,
                                      section_id="zeeman_populations", title="Zeeman populations")
                    with st.container(horizontal_alignment="center", gap="xsmall"):
                        st.markdown("Axis", text_alignment="center", width=SECTION_CONTROL_WIDTH)
                        st.selectbox("Axis", ["z", "x", "y"], width=SECTION_CONTROL_WIDTH,
                                     label_visibility="collapsed", key=key(f"zeeman_axis_{label}"))
        if moment_section.open and time_ready:
            with moment_section:
                moments = species["moments"]
                moment_plots = []
                q = species["moment_axis"]
                u, v = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}[q]
                # Shared A/B display choices live only on the dynamics page.
                with st.container(horizontal=True, vertical_alignment="center", wrap=False, gap="xsmall"):
                    st.markdown("Axis", width="content")
                    st.selectbox("Axis", ["z", "x", "y"], width=SECTION_CONTROL_WIDTH, label_visibility="collapsed",
                                 key=key(f"moment_axis_{label}"))
                    st.space(20)
                    with st.container(horizontal=True, width="stretch", gap="small"):
                        if label == "A":
                            for name in PLOTTED_MOMENTS:
                                option_key = key(f"moment_option_{name}")
                                st.session_state[option_key] = name in st.session_state[key("plotted_moments")]
                                expression = _moment_expression(name, q)
                                if name == "Hyperfine polarization":
                                    caption = r"$\langle\mathbf S\cdot\mathbf I\rangle/\hbar^2$"
                                elif name == "Electron polarization":
                                    caption = r"$2\langle S_q\rangle/\hbar$"
                                else:
                                    caption = rf"$\langle {expression}\rangle$"
                                st.checkbox(caption, key=option_key, on_change=_store_plot_choice,
                                            args=("plotted_moments", name, option_key), width="content")
                selected_moments = st.session_state[key("plotted_moments")]
                frame_explanation = (
                    rf"The right-handed frame is $(\hat{{\mathbf u}},\hat{{\mathbf v}},\hat{{\mathbf q}})"
                    rf"=(\hat{{\mathbf {u}}},\hat{{\mathbf {v}}},\hat{{\mathbf {q}}})$. "
                )
                plot_specs = {
                    "Hyperfine polarization": ("hyperfine_polarization", "⟨S·I⟩/ℏ²", "Hyperfine polarization",
                        "Hyperfine polarization is the electron–nuclear spin correlation "
                        r"$\langle\mathbf S\cdot\mathbf I\rangle/\hbar^2"
                        r"=[\langle\mathbf F^2\rangle/\hbar^2-I(I+1)-3/4]/2$. "
                        r"It equals $(I/2)P(F_a)-[(I+1)/2]P(F_b)$, per total atom. "
                        "It is initially zero and independent of quantization axis."),
                    "Electron polarization": ("electron_polarization", "2⟨S_q⟩/ℏ", "Electron polarization",
                        rf"Electron polarization is $2\langle S_q\rangle/\hbar$ along $q={q}$, including both ground "
                        "hyperfine manifolds. It ranges from −1 to +1."),
                    **{f"F_{axis}": (f"orientation_{axis}", f"⟨F_{axis}⟩/ℏ", f"⟨F_{axis}⟩",
                        rf"Orientation $\langle F_{axis}\rangle/\hbar$ along the laboratory ${axis}$ axis, "
                        "including density-matrix coherences. ")
                       for axis in "xyz"},
                    "Q_qq": ("alignment", f"⟨Q_{q}{q}⟩/ℏ²", f"⟨Q_{q}{q}⟩",
                        rf"Longitudinal alignment is $Q_{{{q}{q}}}=F_{q}^2-F(F+1)\hbar^2/3$. " + frame_explanation),
                    "Q_uv": ("alignment_uv", f"⟨Q_{u}{v}⟩/ℏ²", f"⟨Q_{u}{v}⟩",
                        rf"Transverse mixed alignment is $Q_{{{u}{v}}}=(F_{u}F_{v}+F_{v}F_{u})/2$. " + frame_explanation),
                    "Q_uu-Q_vv": ("alignment_uu_minus_vv", f"⟨Q_{u}{u} − Q_{v}{v}⟩/ℏ²", f"⟨Q_{u}{u} − Q_{v}{v}⟩",
                        rf"Transverse alignment difference is $Q_{{{u}{u}}}-Q_{{{v}{v}}}=F_{u}^2-F_{v}^2$. " + frame_explanation),
                }
                moment_columns = st.columns([1, 1, 0.15], gap="small")
                for index, choice in enumerate(selected_moments):
                    with moment_columns[index % 2]:
                        observable, y_label, title, explanation = plot_specs[choice]
                        if choice in ("Hyperfine polarization", "Electron polarization"):
                            series = {y_label: moments[observable]}
                        else:
                            series = {manifold_labels[name]: moments[f"{observable}_{name}"]
                                      for name in ("Fa", "Fb")}
                            explanation += ("Moments are per total atom, without dividing by manifold population. "
                                            "Both hyperfine manifolds are always shown.")
                        expression = _moment_expression(choice, q)
                        _plot(series_frame(times, series, y_label), y_label, title,
                              species=species, label=label, plot_id=observable, value_unit="1",
                              explanation=explanation, section_plots=moment_plots,
                              show_legend=choice not in ("Hyperfine polarization", "Electron polarization"),
                              math_title=rf"\langle {expression}\rangle" if expression else None)
                with moment_columns[2]:
                    _section_controls(moment_plots, species=species, label=label,
                                      section_id="atomic_moments", title="Atomic moments")
        # The frequency scan always uses all signals, independently of the
        # time-domain section's open state and selected checkboxes.
        plots = (
            ("rotation_rad", 1000, "Rotation (mrad)", "Optical rotation", "mrad",
             "Rotation is relative to the input azimuth, in mrad. It is undefined for circular input "
             "polarization or effectively extinguished output."),
            ("ellipticity_rad", 1000, "Ellipticity (mrad)", "Ellipticity", "mrad",
             "Ellipticity is the transmitted ellipse angle in mrad, rather than its change from the "
             "input angle. It is undefined at effectively extinguished output."),
            ("absorption", 1, "Absorption 1 − I_out/I_in", "Absorption", "1",
             r"Absorption is the lost intensity fraction $1-I_{\mathrm{out}}/I_{\mathrm{in}}$, including scalar absorption. "
             "Zero means no loss and one means complete attenuation. It is not logarithmic absorbance."),
        )
        if probe_section.open and time_ready:
            with probe_section:
                with st.container(horizontal=True, gap="small"):
                    for name in PLOTTED_PROBE_SIGNALS:
                        option_key = key(f"probe_option_{label}_{name}")
                        st.session_state[option_key] = name in st.session_state[key(f"plotted_probe_signals_{label}")]
                        st.checkbox(name, key=option_key, on_change=_store_plot_choice,
                                    args=(f"plotted_probe_signals_{label}", name, option_key), width="content")
                probe_plot_specs = {plot[3]: plot for plot in plots}
                selected_plots = [probe_plot_specs[name]
                                  for name in st.session_state[key(f"plotted_probe_signals_{label}")]]
                probe_plots = []
                signals = None
                if selected_plots:
                    try:
                        signals = _probe_cached(physical_species, probes[label], common, config["n2_coeffs"])
                    except (ValueError, RuntimeError) as exc:
                        st.error("Weak-probe readout failed: " + str(exc))
                probe_columns = st.columns([1, 1, 0.15], gap="small")
                if signals is not None:
                    for index, (signal, scale, y_label, title, value_unit, explanation) in enumerate(selected_plots):
                        with probe_columns[index % 2]:
                            _plot(series_frame(times, {title: signals[signal] * scale}, y_label), y_label, title,
                                  species=species, label=label, plot_id=f"probe_{signal}", value_unit=value_unit,
                                  section_plots=probe_plots,
                                  show_legend=False,
                                  explanation=explanation + " The continuous weak probe uses uniform-path Jones "
                                  "propagation and does not alter the atomic evolution.")
                with probe_columns[2]:
                    _section_controls(probe_plots, species=species, label=label,
                                      section_id="probe_signal", title="Probe signal")
        if scan_section.open and scan_ready:
            with scan_section:
                scan_columns = st.columns([.35, 1, .8], gap="small", vertical_alignment="center")
                with scan_columns[0]:
                    st.markdown("Detuning range (MHz)")
                with scan_columns[1]:
                    _scan_range_slider(label)
                lower = st.session_state[key(f"probe_scan_min_MHz_{label}")]
                upper = st.session_state[key(f"probe_scan_max_MHz_{label}")]
                scan_plots, scan_export = [], None
                try:
                    with st.spinner("Converging the pump-on steady state and scanning the weak probe…"):
                        steady = _steady_cached(species_A, species_B, physical_beams, common)
                        spectrum = _spectrum_cached(steady[label], probes[label], common,
                                                    config["n2_coeffs"], lower, upper)
                    scan_export = probe_scan_export_dataframe(
                        spectrum, species=steady[label], label=label, probe=probes[label],
                        config=config, common=common, beams=beams, steady=steady,
                    )
                    references = [dict(center=float(center), transition=f"F={Fg:g} → F′={Fe:g}")
                                  for (Fg, Fe), center in spectrum["transition_centers_MHz"].items()
                                  if lower <= center <= upper]
                    for signal, scale, y_label, title, value_unit, explanation in plots:
                        frame = series_frame(spectrum["detuning_MHz"], {title: scale * spectrum[signal]}, y_label)
                        frame = frame.rename(columns={"Time (s)": "Detuning (MHz)"})
                        _plot(frame, y_label, title, species=species, label=label,
                              plot_id=f"probe_scan_{signal}", value_unit=value_unit,
                              explanation=explanation, section_plots=scan_plots, show_legend=False,
                              x_label="Detuning (MHz)", x_domain=[lower, upper], show_switch=False,
                              reference_lines=references)
                except (ValueError, RuntimeError) as exc:
                    st.error("Probe scan failed: " + str(exc))
                with scan_columns[2]:
                    _section_controls(
                        scan_plots, species=species, label=label, section_id="probe_scan",
                        title="Probe signal scanning", export_frame=scan_export, horizontal=True,
                        help_text=rf"""
The atoms are independently converged to a **pump-on steady state** using the
same coherent master equation, including all configured pumps, magnetic
precession, light shifts, electron randomization, uniform relaxation and spin
exchange. The total-time and pump turn-off sliders do not set this state.

The weak probe scans the sidebar's **{probes[label]['line']}** reference line,
using its direction, input polarization and path length. All allowed ground
and excited hyperfine transitions contribute. The fixed sidebar detuning and
its F→F′ reference do not limit this scan. Detuning is relative to the
**zero-pressure fine-structure D-line center**; Doppler/pressure broadening and
pressure shifts remain included. Light dashed lines identify allowed transition
centers; hover for their F→F′ labels. The probe does not change the steady state.

The shorter range slider at the top left has its label on the left and two handles
for minimum and maximum detuning within
$\pm 15000\,\mathrm{{MHz}}$. They default to $-7000$ and $+7000\,\mathrm{{MHz}}$,
and must remain at least $10\,\mathrm{{MHz}}$ apart. There are 2001 scan samples.
Rotation is relative to the input azimuth, ellipticity is the output angle,
and absorption is $1-I_{{\mathrm{{out}}}}/I_{{\mathrm{{in}}}}$.

The CSV contains all three signals in plotted units, the detuning origin,
probe/cell/pump configuration and steady-state convergence information.
Undefined optical angles remain empty. Scan limits are independent per species.
""",
                    )
