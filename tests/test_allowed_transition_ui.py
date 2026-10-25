"""Both analysis sidebars enforce allowed optical references without a toggle."""

from pathlib import Path
import unittest

from streamlit.testing.v1 import AppTest


class AllowedTransitionUiTests(unittest.TestCase):
    def test_both_pages_always_filter_pump_and_probe_reference_menus(self):
        app = AppTest.from_file(
            str(Path(__file__).resolve().parents[1] / "alkali_pumping.py"), default_timeout=60,
        ).run()
        self.assertFalse(app.exception)
        # Stale in-memory values from a previous running release cannot disable filtering.
        app.session_state["show_allowed_only"] = False
        app.selectbox(key="_pump_widget_line_A1").select("D2").run()
        app.selectbox(key="_probe_widget_probe_source_A").select("Custom").run()
        app.selectbox(key="_probe_widget_probe_line_A").select("D2").run()
        self.assertFalse(app.exception)
        self.assertFalse(any(widget.label == "Only show allowed hyperfine transitions" for widget in app.checkbox))
        expected = {"1→0", "1→1", "1→2", "2→1", "2→2", "2→3"}
        for key in ("_pump_widget_transition_A1", "_probe_widget_probe_transition_A"):
            self.assertEqual(set(app.selectbox(key=key).options), expected)
        app.session_state["te_show_allowed_only"] = False
        app.switch_page("alkali_pumping_app/pages/time_evolution.py").run()
        app.selectbox(key="te_line_A1").select("D2").run()
        app.selectbox(key="te_probe_line_A").select("D2").run()
        self.assertFalse(app.exception)
        self.assertFalse(any(widget.label == "Only show allowed hyperfine transitions" for widget in app.checkbox))
        for key in ("te_transition_A1", "te_probe_transition_A"):
            self.assertEqual(set(app.selectbox(key=key).options), expected)


if __name__ == "__main__":
    unittest.main()
