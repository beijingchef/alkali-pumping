"""Exercise the real dynamics page, slider bounds, and page-state isolation."""

import json
from pathlib import Path
import unittest
from unittest.mock import patch

from streamlit.testing.v1 import AppTest


ROOT = Path(__file__).resolve().parents[1]
PAGE = "alkali_pumping_app/pages/time_evolution.py"


class TimeEvolutionUiTests(unittest.TestCase):
    def set_moments(self, app, choices):
        from alkali_pumping_app.ui.time_evolution_conditions import PLOTTED_MOMENTS
        for name in PLOTTED_MOMENTS:
            app.checkbox(key=f"te_moment_option_{name}").set_value(name in choices)
        return app.run()

    def app(self, *, open_moments=True):
        app = AppTest.from_file(str(ROOT / "alkali_pumping.py"), default_timeout=60).run()
        self.assertFalse(app.exception)
        with patch("alkali_pumping_app.physics.system_model._evaluate_system_responses",
                   side_effect=AssertionError("Time evolution must not run RF")):
            app.switch_page(PAGE).run()
        self.assertFalse(app.exception)
        if open_moments:
            app.session_state["te_atomic_moments_open_A"] = True
            app.run()
            self.assertFalse(app.exception)
        return app

    def test_page_renders_required_plots_and_two_sliders(self):
        app = self.app(open_moments=False)
        for title in ("Atomic moments", "Probe signal scanning"):
            section = next(item for item in app.expander if item.label == title)
            self.assertFalse(section.proto.expanded)
            self.assertFalse(section.get("vega_lite_chart"))
            self.assertFalse(section.get("download_button"))
        self.assertEqual(len(app.get("vega_lite_chart")), 5)
        self.assertEqual(len(app.main.get("download_button")), 2)
        app.session_state["te_atomic_moments_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="te_gamma_u_A").value, 1.)
        self.assertEqual(app.number_input(key="te_gamma_u_B").value, 1.)
        relaxation_keys = [widget.key for widget in app.sidebar.number_input
                           if widget.key.startswith(("te_gamma_ER_", "te_gamma_u_"))]
        self.assertEqual(relaxation_keys, ["te_gamma_ER_A", "te_gamma_ER_B", "te_gamma_u_A", "te_gamma_u_B"])
        self.assertEqual(len(app.slider), 2)
        self.assertEqual(len(app.sidebar.slider), 2)
        self.assertEqual(len(app.main.slider), 0)
        self.assertIn("Probe signal scanning", [expander.label for expander in app.expander])
        self.assertEqual(app.sidebar.header[-1].value, "Event Sequence")
        self.assertNotIn("Display", [header.value for header in app.sidebar.header])
        self.assertFalse(any(widget.label == "Only show allowed hyperfine transitions" for widget in app.checkbox))
        slider_columns = [column for column in app.sidebar.get("column") if len(column.slider)]
        self.assertEqual([len(column.slider) for column in slider_columns], [1, 1])
        self.assertEqual(app.slider(key="_time_widget_te_t2_s").proto.max, 30.)
        self.assertEqual(app.slider(key="_time_widget_te_t2_s").label, "Total time (s)")
        self.assertEqual(app.slider(key="_time_widget_te_t1_s").label, "Pump turn-off time (s)")
        self.assertEqual(len(app.get("vega_lite_chart")), 9)
        section_downloads = [button for button in app.main.get("download_button")
                             if "te_download_section_" in button.proto.id]
        self.assertEqual(len(section_downloads), 3)
        self.assertTrue(all(button.proto.ignore_rerun for button in section_downloads))
        self.assertTrue(all(button.proto.label == ":material/download:" for button in section_downloads))
        self.assertEqual(len(app.main.get("popover")), 4)
        self.assertEqual(len(app.main.caption), 0)
        for chart in app.get("vega_lite_chart"):
            spec = json.loads(chart.proto.spec)
            self.assertEqual(spec["height"], 260)
            legend = spec["layer"][0]["encoding"]["color"]["legend"]
            if spec["description"] in ("Optical rotation", "Ellipticity", "Absorption",
                                       "Hyperfine polarization", "Electron polarization"):
                self.assertIsNone(legend)
            else:
                self.assertEqual(legend["orient"], "top-right")
            if spec["description"] == "Absorption":
                encoding = spec["layer"][0]["encoding"]
                self.assertEqual(encoding["y"]["axis"]["format"], ".2e")
                self.assertEqual(encoding["tooltip"][-1]["format"], ".2e")
        titles = [json.loads(chart.proto.spec)["description"] for chart in app.get("vega_lite_chart")]
        self.assertNotIn("Total population in Fa and Fb", titles)
        self.assertNotIn("Hyperfine populations", [expander.label for expander in app.expander])
        self.assertIn("Hyperfine polarization", titles)
        self.assertIn("Electron polarization", titles)
        self.assertIn("F = 2", titles)
        self.assertIn("F = 1", titles)
        self.assertIn("Ellipticity", titles)
        self.assertIn("Absorption", titles)
        self.assertNotIn("Transmission", titles)
        self.assertFalse(any(key.startswith("rf_") for key in [widget.key for widget in app.selectbox]))

    def test_probe_scan_opens_on_demand_and_retains_limits_when_closed(self):
        app = self.app()
        app.session_state["te_probe_scan_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 12)
        self.assertEqual(len(app.main.slider), 1)
        self.assertEqual(app.slider(key="_scan_widget_te_probe_scan_range_MHz_A").value,
                         (-7000., 7000.))
        app.slider(key="_scan_widget_te_probe_scan_range_MHz_A").set_value((-5000., 7000.)).run()
        app.session_state["te_probe_scan_open_A"] = False
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 9)
        self.assertEqual(len(app.main.slider), 0)
        self.assertEqual(app.session_state["te_probe_scan_min_MHz_A"], -5000.)
        app.session_state["te_probe_scan_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.slider(key="_scan_widget_te_probe_scan_range_MHz_A").value,
                         (-5000., 7000.))

    def test_sliders_clamp_switch_and_support_zero_and_end_time(self):
        app = self.app()
        app.slider(key="_time_widget_te_t2_s").set_value(30.).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.slider(key="_time_widget_te_t2_s").proto.max, 30.)
        self.assertEqual(app.slider(key="_time_widget_te_t1_s").proto.max, 30.)
        app.slider(key="_time_widget_te_t1_s").set_value(30.).run()
        self.assertFalse(app.exception)
        app.slider(key="_time_widget_te_t2_s").set_value(.1).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.slider(key="_time_widget_te_t1_s").value, .1)
        app.slider(key="_time_widget_te_t1_s").set_value(0.).run()
        self.assertFalse(app.exception)
        for chart in app.get("vega_lite_chart"):
            spec = json.loads(chart.proto.spec)
            if spec["layer"][0]["encoding"]["x"]["field"] != "Time (s)":
                continue
            self.assertEqual(spec["layer"][0]["encoding"]["x"]["scale"]["domain"], [0, .1])
        app.slider(key="_time_widget_te_t2_s").set_value(.001).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.slider(key="_time_widget_te_t1_s").value, 0.)
        self.assertEqual(app.slider(key="_time_widget_te_t2_s").proto.max, 30.)

    def test_sidebar_and_moment_choices_persist_without_changing_acc(self):
        app = self.app()
        app.number_input(key="te_intensity_A1").set_value(35.).run()
        app.number_input(key="te_gamma_ER_A").set_value(3.).run()
        app.number_input(key="te_gamma_u_A").set_value(4.).run()
        self.set_moments(app, ["Q_qq"])
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["intensity_A1"], 10.)
        self.assertEqual(app.session_state["gamma_ER_A"], 10.)
        self.assertEqual(app.session_state["gamma_u_A"], 1.)
        app.switch_page("alkali_pumping_app/pages/ACC-Magnetometer.py").run()
        self.assertFalse(app.exception)
        app.number_input(key="gamma_ER_A").set_value(16.).run()
        app.number_input(key="gamma_u_A").set_value(7.).run()
        app.switch_page(PAGE)
        app.session_state["te_atomic_moments_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="te_intensity_A1").value, 35.)
        self.assertEqual(app.number_input(key="te_gamma_ER_A").value, 3.)
        self.assertEqual(app.number_input(key="te_gamma_u_A").value, 4.)
        self.assertEqual(app.session_state["gamma_u_A"], 7.)
        self.assertEqual(app.session_state["te_plotted_moments"], ["Q_qq"])
        self.assertTrue(app.checkbox(key="te_moment_option_Q_qq").value)
        self.assertEqual(app.session_state["gamma_ER_A"], 16.)

    def test_all_moments_and_empty_selection_control_only_moment_graphs(self):
        from alkali_pumping_app.ui.time_evolution_conditions import PLOTTED_MOMENTS
        app = self.app()
        self.assertNotIn("Plotted moments", [item.value for item in app.markdown])
        self.assertEqual(app.checkbox(key="te_moment_option_Hyperfine polarization").label,
                         r"$\langle\mathbf S\cdot\mathbf I\rangle/\hbar^2$")
        self.assertEqual(app.checkbox(key="te_moment_option_Electron polarization").label,
                         r"$2\langle S_q\rangle/\hbar$")
        self.assertEqual(app.checkbox(key="te_moment_option_Q_uv").label, r"$\langle Q_{xy}\rangle$")
        self.set_moments(app, list(PLOTTED_MOMENTS))
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 13)
        specs = [json.loads(chart.proto.spec) for chart in app.get("vega_lite_chart")]
        moment_titles = {"⟨F_x⟩", "⟨F_y⟩", "⟨F_z⟩", "⟨Q_zz⟩", "⟨Q_xy⟩", "⟨Q_xx − Q_yy⟩"}
        self.assertTrue(moment_titles <= {spec["description"] for spec in specs})
        for spec in specs:
            if spec["description"] in moment_titles:
                self.assertEqual(spec["layer"][0]["encoding"]["color"]["sort"], ["F = 2", "F = 1"])
        self.set_moments(app, [])
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 5)
        self.assertEqual(len(app.main.get("download_button")), 3)
        moment_download = next(button for button in app.main.get("download_button")
                               if "te_download_section_A_atomic_moments" in button.proto.id)
        self.assertTrue(moment_download.proto.disabled)

    def test_dual_species_and_custom_probe_controls(self):
        app = self.app()
        app.selectbox(key="te_atom_B_name").select("Rb85").run()
        app.session_state["te_atomic_moments_open_B"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="te_transition_B1").value, "2→3")
        self.assertEqual(app.selectbox(key="te_transition_B2").value, "3→3")
        self.assertEqual(app.selectbox(key="te_probe_transition_B").value, "3→2")
        self.assertEqual(len(app.get("vega_lite_chart")), 18)
        self.assertEqual(len(app.main.get("download_button")), 6)
        self.assertEqual(len(app.main.get("popover")), 7)
        app.selectbox(key="te_probe_source_A").select("Custom").run()
        app.selectbox(key="te_probe_k_A").select("z").run()
        app.number_input(key="te_probe_ellipticity_deg_A").set_value(45.).run()
        self.assertFalse(app.exception)
        self.assertTrue(any("undefined" in info.value for info in app.info))
        # Returning to an inactive B must recover the one-species page.
        app.selectbox(key="te_atom_B_name").select("Rb87").run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 8)

    def test_transverse_field_and_independent_quantization_axes_render(self):
        app = self.app()
        app.slider(key="_time_widget_te_t2_s").set_value(.04).run()
        app.selectbox(key="te_atom_B_name").select("Rb85").run()
        app.session_state["te_atomic_moments_open_B"] = True
        app.run()
        app.selectbox(key="te_zeeman_axis_A").select("x").run()
        app.selectbox(key="te_moment_axis_A").select("y").run()
        app.selectbox(key="te_zeeman_axis_B").select("y").run()
        app.selectbox(key="te_moment_axis_B").select("x").run()
        app.number_input(key="te_static_field_nT").set_value(.1).run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.assertFalse(any("parallel" in warning.value for warning in app.warning))
        self.assertEqual(len(app.get("vega_lite_chart")), 18)
        self.assertEqual(app.session_state["te_zeeman_axis_A"], "x")
        self.assertEqual(app.session_state["te_moment_axis_A"], "y")
        self.assertEqual(app.checkbox(key="te_moment_option_Q_qq").label, r"$\langle Q_{yy}\rangle$")
        self.assertEqual(app.checkbox(key="te_moment_option_Q_uv").label, r"$\langle Q_{zx}\rangle$")
        self.assertEqual(app.checkbox(key="te_moment_option_Q_uu-Q_vv").label,
                         r"$\langle Q_{zz}-Q_{xx}\rangle$")
        self.assertEqual(app.session_state["te_zeeman_axis_B"], "y")
        self.assertEqual(app.session_state["te_moment_axis_B"], "x")
        self.assertFalse(any("Quantization axis" in widget.label or "quantization axis" in widget.label
                             for widget in app.sidebar.selectbox))
        self.assertEqual(app.session_state["static_field_nT"], 0.)


if __name__ == "__main__":
    unittest.main()
