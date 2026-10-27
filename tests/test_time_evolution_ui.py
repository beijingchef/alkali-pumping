"""Exercise the real dynamics page, slider bounds, and page-state isolation."""

import base64
from contextlib import ExitStack
import json
from pathlib import Path
import re
import unittest
from unittest.mock import patch

from streamlit.proto.WidgetStates_pb2 import WidgetState
from streamlit.testing.v1 import AppTest
from streamlit.testing.v1 import element_tree


ROOT = Path(__file__).resolve().parents[1]
PAGE = "alkali_pumping_app/pages/time_evolution.py"


class TimeEvolutionUiTests(unittest.TestCase):
    def setUp(self):
        # AppTest currently serializes only Widget nodes. Stateful expanders
        # are Block nodes, so their bool values are omitted from reruns and
        # reset to collapsed. Submit their real registered state as a browser
        # does; do not force sections open or bypass the app's render gates.
        original_get_widget_state = element_tree.get_widget_state

        def widget_state_with_expanders(node):
            state = original_get_widget_state(node)
            if state is None and isinstance(node, element_tree.Expander) and node.proto.id:
                return WidgetState(
                    id=node.proto.id,
                    bool_value=node.root._runner.session_state[node.proto.id],
                )
            return state

        adapter = patch.object(element_tree, "get_widget_state", side_effect=widget_state_with_expanders)
        adapter.start()
        self.addCleanup(adapter.stop)

    def set_moments(self, app, choices):
        from alkali_pumping_app.ui.time_evolution_conditions import PLOTTED_MOMENTS
        for name in PLOTTED_MOMENTS:
            app.checkbox(key=f"te_moment_option_{name}").set_value(name in choices)
        return app.run()

    def app(self, *, open_moments=True, open_zeeman=True):
        app = AppTest.from_file(str(ROOT / "alkali_pumping.py"), default_timeout=60).run()
        self.assertFalse(app.exception)
        with patch("alkali_pumping_app.physics.system_model._evaluate_system_responses",
                   side_effect=AssertionError("Time evolution must not run RF")):
            app.switch_page(PAGE).run()
        self.assertFalse(app.exception)
        if open_moments or open_zeeman:
            app.session_state["te_atomic_moments_open_A"] = open_moments
            app.session_state["te_zeeman_populations_open_A"] = open_zeeman
            app.run()
            self.assertFalse(app.exception)
        return app

    def test_page_renders_required_plots_and_two_sliders(self):
        app = self.app(open_moments=False, open_zeeman=False)
        for title in ("Zeeman populations", "Atomic moments", "Probe signal scanning"):
            section = next(item for item in app.expander if item.label == title)
            self.assertFalse(section.proto.expanded)
            self.assertFalse(section.get("vega_lite_chart"))
            self.assertFalse(section.get("download_button"))
        self.assertEqual(len(app.get("vega_lite_chart")), 2)
        self.assertEqual(len(app.main.get("download_button")), 1)
        self.assertEqual(app.session_state["te_plotted_moments"], ["Hyperfine polarization", "F_z"])
        self.assertTrue(app.checkbox(key="te_probe_option_A_Absorption").value)
        self.assertTrue(app.checkbox(key="te_probe_option_A_Optical rotation").value)
        self.assertFalse(app.checkbox(key="te_probe_option_A_Ellipticity").value)
        app.session_state["te_atomic_moments_open_A"] = True
        app.session_state["te_zeeman_populations_open_A"] = True
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
        self.assertEqual(len(app.get("vega_lite_chart")), 6)
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
        self.assertNotIn("Electron polarization", titles)
        self.assertIn("⟨F_z⟩", titles)
        self.assertNotIn("⟨Q_zz⟩", titles)
        self.assertIn("F = 2", titles)
        self.assertIn("F = 1", titles)
        self.assertNotIn("Ellipticity", titles)
        self.assertIn("Absorption", titles)
        self.assertNotIn("Transmission", titles)
        self.assertFalse(any(key.startswith("rf_") for key in [widget.key for widget in app.selectbox]))

    def test_probe_scan_opens_on_demand_and_retains_limits_when_closed(self):
        app = self.app()
        app.session_state["te_probe_scan_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 9)
        self.assertEqual(len(app.main.slider), 1)
        self.assertEqual(app.slider(key="_scan_widget_te_probe_scan_range_MHz_A").value,
                         (-7000., 7000.))
        app.slider(key="_scan_widget_te_probe_scan_range_MHz_A").set_value((-5000., 7000.)).run()
        app.session_state["te_probe_scan_open_A"] = False
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 6)
        self.assertEqual(len(app.main.slider), 0)
        self.assertEqual(app.session_state["te_probe_scan_min_MHz_A"], -5000.)
        app.session_state["te_probe_scan_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.slider(key="_scan_widget_te_probe_scan_range_MHz_A").value,
                         (-5000., 7000.))

    def test_sidebar_changes_only_render_open_sections(self):
        app = self.app(open_moments=False, open_zeeman=False)
        sections = {
            "Zeeman populations": ("zeeman_populations", 2),
            "Atomic moments": ("atomic_moments", 2),
            "Probe signal": ("probe_signal", 2),
            "Probe signal scanning": ("probe_scan", 3),
        }
        for index, opened in enumerate([None, *sections]):
            for title, (state_key, _) in sections.items():
                app.session_state[f"te_{state_key}_open_A"] = title == opened
            # A changed physical condition must not render or read out closed sections.
            with ExitStack() as stack:
                if opened in (None, "Probe signal scanning"):
                    stack.enter_context(patch("alkali_pumping_app.physics.time_evolution.compute_density_evolution",
                                              side_effect=AssertionError("Closed transient sections integrated")))
                for title, target in (
                    ("Zeeman populations", "alkali_pumping_app.ui.time_evolution_data.zeeman_plot_frame"),
                    ("Probe signal", "alkali_pumping_app.physics.optical_readout.density_probe_signals"),
                    ("Probe signal scanning", "alkali_pumping_app.physics.probe_scan.pump_on_steady_state"),
                ):
                    if title != opened:
                        stack.enter_context(patch(target, side_effect=AssertionError(f"Closed {title} ran")))
                app.number_input(key="te_gamma_ER_A").set_value(10. + index).run()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            self.assertEqual(len(app.get("vega_lite_chart")), 0 if opened is None else sections[opened][1])
            for title, (_, count) in sections.items():
                section = next(item for item in app.expander if item.label == title)
                self.assertEqual(len(section.get("vega_lite_chart")), count if title == opened else 0)
                self.assertEqual(len(section.get("download_button")), int(title == opened))

    def test_probe_checkboxes_filter_exports_persist_and_leave_scan_independent(self):
        from alkali_pumping_app.ui import exports

        app = self.app(open_moments=False, open_zeeman=False)
        app.checkbox(key="te_probe_option_A_Optical rotation").uncheck()
        app.checkbox(key="te_probe_option_A_Ellipticity").check()
        with patch.object(exports, "dataframe_to_csv_bytes", wraps=exports.dataframe_to_csv_bytes) as export:
            app.run()
        self.assertFalse(app.exception)
        frame = export.call_args.args[0]
        self.assertIn("Absorption | Absorption [1]", frame)
        self.assertIn("Ellipticity | Ellipticity [mrad]", frame)
        self.assertNotIn("Optical rotation | Optical rotation [mrad]", frame)
        self.assertEqual(app.session_state["te_plotted_probe_signals_A"], ["Absorption", "Ellipticity"])
        app.session_state["te_probe_signal_open_A"] = False
        app.number_input(key="te_probe_det_rel_A").set_value(-2100.).run()
        self.assertFalse(app.exception)
        self.assertFalse([w for w in app.main.checkbox if w.key.startswith("te_probe_option_")])
        app.session_state["te_probe_signal_open_A"] = True
        app.run()
        self.assertTrue(app.checkbox(key="te_probe_option_A_Ellipticity").value)
        self.assertFalse(app.checkbox(key="te_probe_option_A_Optical rotation").value)

        app.checkbox(key="te_probe_option_A_Absorption").uncheck()
        app.checkbox(key="te_probe_option_A_Ellipticity").uncheck().run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 0)
        self.assertTrue(app.main.get("download_button")[0].proto.disabled)
        app.session_state["te_probe_scan_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.assertEqual(len(app.get("vega_lite_chart")), 3)

        app.session_state["te_probe_scan_open_A"] = False
        app.selectbox(key="te_atom_B_name").select("Rb85").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["te_plotted_probe_signals_A"], [])
        self.assertEqual(app.session_state["te_plotted_probe_signals_B"], ["Absorption", "Optical rotation"])
        app.switch_page("alkali_pumping_app/pages/ACC-Magnetometer.py").run()
        app.switch_page(PAGE).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["te_plotted_probe_signals_A"], [])
        self.assertEqual(app.session_state["te_plotted_probe_signals_B"], ["Absorption", "Optical rotation"])

    def test_plot_order_follows_checks_and_survives_hidden_sections_and_navigation(self):
        from alkali_pumping_app.ui import exports

        app = self.app(open_zeeman=False)

        def assert_order(title, expected, *, tab_index=0):
            section = next(item for item in app.main.tabs[tab_index].expander if item.label == title)
            # Check visual row order: left, right, then the next row.
            for column, titles in zip(section.get("column")[:2], (expected[::2], expected[1::2])):
                self.assertEqual([json.loads(chart.proto.spec)["description"]
                                  for chart in column.get("vega_lite_chart")], titles)

        moment_order = ["Hyperfine polarization", "⟨F_z⟩"]
        probe_order = ["Absorption", "Optical rotation"]
        assert_order("Atomic moments", moment_order)
        assert_order("Probe signal", probe_order)
        for section, option, checked, expected in (
            ("Atomic moments", "moment_option_F_y", True, [*moment_order, "⟨F_y⟩"]),
            ("Atomic moments", "moment_option_Electron polarization", True,
             [*moment_order, "⟨F_y⟩", "Electron polarization"]),
            ("Atomic moments", "moment_option_Hyperfine polarization", False,
             ["⟨F_z⟩", "⟨F_y⟩", "Electron polarization"]),
            ("Atomic moments", "moment_option_Hyperfine polarization", True,
             ["⟨F_z⟩", "⟨F_y⟩", "Electron polarization", "Hyperfine polarization"]),
            ("Probe signal", "probe_option_A_Ellipticity", True, [*probe_order, "Ellipticity"]),
            ("Probe signal", "probe_option_A_Absorption", False, ["Optical rotation", "Ellipticity"]),
            ("Probe signal", "probe_option_A_Absorption", True, ["Optical rotation", "Ellipticity", "Absorption"]),
        ):
            with self.subTest(section=section, option=option, checked=checked):
                with patch.object(exports, "dataframe_to_csv_bytes", wraps=exports.dataframe_to_csv_bytes) as export:
                    app.checkbox(key="te_" + option).set_value(checked).run()
                self.assertFalse(app.exception)
                assert_order(section, expected)
                frame = next(call.args[0] for call in export.call_args_list
                             if not call.args[0].empty and call.args[0]["section"].iloc[0] == section)
                self.assertEqual(list(dict.fromkeys(column.split(" | ")[0] for column in frame if " | " in column)),
                                 expected)

        moment_order = ["⟨F_z⟩", "⟨F_y⟩", "Electron polarization", "Hyperfine polarization"]
        probe_order = ["Optical rotation", "Ellipticity", "Absorption"]
        for name in ("atomic_moments", "probe_signal"):
            app.session_state[f"te_{name}_open_A"] = False
        app.number_input(key="te_probe_det_rel_A").set_value(-2200.).run()
        self.assertFalse(app.exception)
        for name in ("atomic_moments", "probe_signal"):
            app.session_state[f"te_{name}_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        assert_order("Atomic moments", moment_order)
        assert_order("Probe signal", probe_order)

        app.selectbox(key="te_atom_B_name").select("Rb85").run()
        app.session_state["te_atomic_moments_open_B"] = True
        app.run()
        self.assertFalse(app.exception)
        assert_order("Atomic moments", moment_order, tab_index=1)
        assert_order("Probe signal", ["Absorption", "Optical rotation"], tab_index=1)
        app.checkbox(key="te_probe_option_B_Absorption").uncheck().run()
        app.checkbox(key="te_probe_option_B_Absorption").check().run()
        self.assertFalse(app.exception)
        assert_order("Probe signal", ["Optical rotation", "Absorption"], tab_index=1)
        assert_order("Probe signal", probe_order)

        app.switch_page("alkali_pumping_app/pages/ACC-Magnetometer.py").run()
        app.switch_page(PAGE)
        app.session_state["te_atomic_moments_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        assert_order("Atomic moments", moment_order)
        assert_order("Probe signal", probe_order)

    def test_manual_apply_batches_physical_edits_and_never_displays_stale_curves(self):
        from alkali_pumping_app.physics import time_evolution as physics
        app = self.app(open_moments=False, open_zeeman=False)
        self.assertTrue(app.checkbox(key="te_auto_update").value)
        app.checkbox(key="te_auto_update").uncheck().run()
        with patch.object(physics, "compute_density_evolution", side_effect=AssertionError("Pending edits integrated")):
            app.number_input(key="te_intensity_A1").set_value(13.).run()
            app.number_input(key="te_gamma_ER_A").set_value(12.).run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 0)
        self.assertEqual(len(app.main.get("download_button")), 0)
        self.assertTrue(any("Pending condition changes" in item.value for item in app.info))
        with patch.object(physics, "compute_density_evolution", wraps=physics.compute_density_evolution) as integrate:
            app.button(key="te_apply_conditions").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(integrate.call_count, 1)
        self.assertEqual(len(app.get("vega_lite_chart")), 2)
        self.assertFalse(any("Pending condition changes" in item.value for item in app.info))
        # Applying with every section closed records approval but still defers work.
        app.session_state["te_probe_signal_open_A"] = False
        app.run()
        with patch.object(physics, "compute_density_evolution", side_effect=AssertionError("Closed sections integrated")):
            app.number_input(key="te_intensity_A1").set_value(14.).run()
            app.button(key="te_apply_conditions").click().run()
        self.assertFalse(app.exception)
        app.session_state["te_probe_signal_open_A"] = True
        with patch.object(physics, "compute_density_evolution", wraps=physics.compute_density_evolution) as integrate:
            app.run()
        self.assertFalse(app.exception)
        self.assertEqual(integrate.call_count, 1)
        self.assertEqual(len(app.get("vega_lite_chart")), 2)
        # Returning to automatic mode calculates pending conditions immediately.
        app.number_input(key="te_intensity_A1").set_value(15.).run()
        app.checkbox(key="te_auto_update").check().run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 2)

    def test_save_includes_current_edits_and_current_format(self):
        from alkali_pumping_app.ui import downloads
        from alkali_pumping_app.ui.time_evolution_conditions import (
            DEFAULTS, SCHEMA_VERSION, parse_settings_payload,
        )

        app = self.app(open_moments=False)
        app.selectbox(key="te_atom_B_name").select("Rb85").run()
        self.assertFalse(app.exception)
        app.checkbox(key="te_probe_option_A_Absorption").uncheck().run()
        app.checkbox(key="te_probe_option_A_Absorption").check().run()
        app.session_state["te_plotted_moments"] = ["F_z", "Hyperfine polarization"]
        edits = {
            "condition_name": "current-settings-save-check",
            "intensity_A1": 35., "gamma_u_B": 7.5,
            "probe_det_rel_A": -1234., "probe_path_length_cm_B": 2.5,
            "t1_s": .6, "t2_s": 3.5,
            "plotted_moments": ["F_z", "Hyperfine polarization"],
            "plotted_probe_signals_A": ["Optical rotation", "Absorption"],
        }
        # Submit pending edits and Save together, as the actual button does.
        app.text_input(key="te_condition_name").set_value(edits["condition_name"])
        for name in ("intensity_A1", "gamma_u_B", "probe_det_rel_A", "probe_path_length_cm_B"):
            app.number_input(key="te_" + name).set_value(edits[name])
        for name in ("t1_s", "t2_s"):
            app.slider(key="_time_widget_te_" + name).set_value(edits[name])
        with patch.object(downloads.st, "html", wraps=downloads.st.html) as render_html:
            app.button(key="te_save").click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        scripts = [call.args[0] for call in render_html.call_args_list
                   if call.args and "const encoded =" in str(call.args[0])]
        self.assertEqual(len(scripts), 1)
        script = scripts[0]
        encoded = re.search(r'const encoded = ("[^"]+");', script)
        self.assertIsNotNone(encoded)
        payload = json.loads(base64.b64decode(json.loads(encoded.group(1))))
        self.assertEqual(payload["version"], SCHEMA_VERSION)
        self.assertIn('link.download = "current-settings-save-check.json";', script)
        self.assertIn('type: "application/json"', script)
        saved = parse_settings_payload(payload)
        self.assertEqual(saved, {name: app.session_state["te_" + name] for name in DEFAULTS})
        for name, value in edits.items():
            self.assertEqual(saved[name], value)
        self.assertEqual(saved["atom_B_name"], "Rb85")

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
        self.assertEqual(len(app.get("vega_lite_chart")), 12)
        specs = [json.loads(chart.proto.spec) for chart in app.get("vega_lite_chart")]
        moment_titles = {"⟨F_x⟩", "⟨F_y⟩", "⟨F_z⟩", "⟨Q_zz⟩", "⟨Q_xy⟩", "⟨Q_xx − Q_yy⟩"}
        self.assertTrue(moment_titles <= {spec["description"] for spec in specs})
        for spec in specs:
            if spec["description"] in moment_titles:
                self.assertEqual(spec["layer"][0]["encoding"]["color"]["sort"], ["F = 2", "F = 1"])
        self.set_moments(app, [])
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("vega_lite_chart")), 4)
        self.assertEqual(len(app.main.get("download_button")), 3)
        moment_download = next(button for button in app.main.get("download_button")
                               if "te_download_section_A_atomic_moments" in button.proto.id)
        self.assertTrue(moment_download.proto.disabled)

    def test_dual_species_and_custom_probe_controls(self):
        app = self.app()
        app.selectbox(key="te_atom_B_name").select("Rb85").run()
        app.session_state["te_atomic_moments_open_B"] = True
        app.session_state["te_zeeman_populations_open_B"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="te_transition_B1").value, "2→3")
        self.assertEqual(app.selectbox(key="te_transition_B2").value, "3→3")
        self.assertEqual(app.selectbox(key="te_probe_transition_B").value, "3→2")
        self.assertEqual(len(app.get("vega_lite_chart")), 12)
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
        self.assertEqual(len(app.get("vega_lite_chart")), 5)

    def test_transverse_field_and_independent_quantization_axes_render(self):
        app = self.app()
        app.slider(key="_time_widget_te_t2_s").set_value(.04).run()
        app.selectbox(key="te_atom_B_name").select("Rb85").run()
        app.session_state["te_atomic_moments_open_B"] = True
        app.session_state["te_zeeman_populations_open_B"] = True
        app.run()
        app.selectbox(key="te_zeeman_axis_A").select("x").run()
        app.selectbox(key="te_moment_axis_A").select("y").run()
        app.selectbox(key="te_zeeman_axis_B").select("y").run()
        app.selectbox(key="te_moment_axis_B").select("x").run()
        app.number_input(key="te_static_field_nT").set_value(.1).run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.assertFalse(any("parallel" in warning.value for warning in app.warning))
        self.assertEqual(len(app.get("vega_lite_chart")), 12)
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
