"""Exercise editable Time evolution presets and the original ACC source menus."""

from pathlib import Path
import unittest

from streamlit.testing.v1 import AppTest


ROOT = Path(__file__).resolve().parents[1]
ACC = "alkali_pumping_app/pages/ACC-Magnetometer.py"
TIME = "alkali_pumping_app/pages/time_evolution.py"
PRESETS = ["Preset D1a", "Preset D1b", "Preset D2a", "Preset D2b"]


class ProbePresetUiTests(unittest.TestCase):
    def assert_probe(self, app, label, expected):
        for field, value in expected.items():
            state_key = f"te_probe_{field}_{label}"
            self.assertEqual(app.session_state[state_key], value)
            kind = app.selectbox if field in ("line", "transition", "k") else app.number_input
            widget = kind(key=state_key)
            self.assertEqual(widget.value, value)
            self.assertFalse(widget.proto.disabled)

    def test_acc_source_menus_do_not_include_time_evolution_presets(self):
        app = AppTest.from_file(str(ROOT / "alkali_pumping.py"), default_timeout=60).run()
        self.assertFalse(app.exception)
        for label in ("A", "B"):
            if label == "B":
                app.session_state["probe_configuration_tab"] = "Probe-B"
                app.run()
            options = app.selectbox(key=f"_probe_widget_probe_source_{label}").options
            expected = ["Custom Weak"]
            for number in (1, 2, 3):
                expected.extend((f"Pump{label}{number} weak", f"Pump{label}{number}"))
            self.assertEqual(options, expected)
        self.assertFalse(app.exception)

    def test_time_evolution_presets_populate_once_and_keep_edits(self):
        app = AppTest.from_file(str(ROOT / "alkali_pumping.py"), default_timeout=60).run()
        acc_source = app.session_state["probe_source_A"]
        app.switch_page(TIME).run()
        self.assertFalse(app.exception)
        original_B_source = app.session_state["te_probe_source_B"]
        pump_values = {field: app.session_state["te_" + field] for field in (
            "intensity_A1", "intensity_A2", "line_A1", "transition_A1", "pol_A1", "k_A1")}
        for source, line, transition, detuning in (
            ("Preset D1a", "D1", "2→1", -2000.),
            ("Preset D1b", "D1", "1→2", 2000.),
            ("Preset D2a", "D2", "2→1", -2000.),
            ("Preset D2b", "D2", "1→2", 2000.),
        ):
            app.selectbox(key="te_probe_source_A").select(source).run()
            self.assertFalse(app.exception)
            self.assert_probe(app, "A", {
                "line": line, "transition": transition, "det_rel": detuning,
                "k": "z", "path_length_cm": 1., "azimuth_deg": 0., "ellipticity_deg": 0.,
            })
            self.assertEqual(app.session_state["te_probe_source_B"], original_B_source)
        for field, expected in pump_values.items():
            self.assertEqual(app.session_state["te_" + field], expected)
        # Every optical field can be edited while the preset stays selected.
        edits = {"line": "D1", "transition": "2→2", "det_rel": 1750., "k": "x",
                 "path_length_cm": 2.5, "azimuth_deg": 20., "ellipticity_deg": 10.}
        for field, value in edits.items():
            if field in ("line", "transition", "k"):
                app.selectbox(key=f"te_probe_{field}_A").select(value)
            else:
                app.number_input(key=f"te_probe_{field}_A").set_value(value)
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["te_probe_source_A"], "Preset D2b")
        self.assert_probe(app, "A", edits)
        app.slider(key="_time_widget_te_t2_s").set_value(5.).run()
        self.assert_probe(app, "A", edits)
        # A different preset repopulates the fields using the current isotope.
        app.selectbox(key="te_atom_A_name").select("Cs133").run()
        app.selectbox(key="te_probe_source_A").select("Preset D2a").run()
        self.assert_probe(app, "A", {
            "line": "D2", "transition": "4→3", "det_rel": -2000., "k": "z",
            "path_length_cm": 1., "azimuth_deg": 0., "ellipticity_deg": 0.,
        })
        app.number_input(key="te_probe_det_rel_A").set_value(-1234.).run()
        app.selectbox(key="te_atom_B_name").select("Rb85").run()
        app.selectbox(key="te_probe_source_B").select("Preset D1b").run()
        self.assert_probe(app, "B", {
            "line": "D1", "transition": "2→3", "det_rel": 2000., "k": "z",
            "path_length_cm": 1., "azimuth_deg": 0., "ellipticity_deg": 0.,
        })
        app.number_input(key="te_probe_path_length_cm_B").set_value(3.).run()
        app.switch_page(ACC).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["probe_source_A"], acc_source)
        self.assertFalse(set(PRESETS) & set(app.selectbox(key="_probe_widget_probe_source_A").options))
        app.switch_page(TIME).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["te_probe_source_A"], "Preset D2a")
        self.assertEqual(app.number_input(key="te_probe_det_rel_A").value, -1234.)
        self.assertEqual(app.session_state["te_probe_source_B"], "Preset D1b")
        self.assertEqual(app.number_input(key="te_probe_path_length_cm_B").value, 3.)


if __name__ == "__main__":
    unittest.main()
