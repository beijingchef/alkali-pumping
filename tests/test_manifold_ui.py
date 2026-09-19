"""Exercise real manifold controls, independent species state, and navigation."""

from pathlib import Path
import unittest

from streamlit.testing.v1 import AppTest


class ManifoldUiTests(unittest.TestCase):
    def test_manifold_controls_rerun_and_survive_navigation(self):
        app = AppTest.from_file(
            str(Path(__file__).resolve().parents[1] / "alkali_pumping.py"),
            default_timeout=60,
        ).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["rf_manifold_A"], "Both")
        app.selectbox(key="_result_widget_rf_manifold_A").select("Fb").run()
        self.assertFalse(app.exception)
        app.selectbox(key="_result_widget_rf_observable_A").select("Q_ij").run()
        app.selectbox(key="_result_widget_probe_manifold_A").select("Fb").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["rf_manifold_A"], "Fb")
        self.assertEqual(app.session_state["probe_manifold_A"], "Fb")
        app.selectbox(key="_result_widget_probe_manifold_A").select("Total").run()
        app.checkbox(key="_result_widget_probe_show_manifold_decomposition_A").check().run()
        app.checkbox(key="_result_widget_probe_show_decomposition_A").check().run()
        self.assertFalse(app.exception)
        app.switch_page("alkali_pumping_app/pages/atomic_properties.py").run()
        self.assertFalse(app.exception)
        app.switch_page("alkali_pumping_app/pages/ACC-Magnetometer.py").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["rf_manifold_A"], "Fb")
        self.assertTrue(app.session_state["probe_show_manifold_decomposition_A"])
        app.selectbox(key="atom_B_name").select("Rb85").run()
        self.assertFalse(app.exception)
        app.session_state["result_species_tab"] = "Alkali B — Rb85"
        app.run()
        self.assertFalse(app.exception)
        app.selectbox(key="_result_widget_rf_manifold_B").select("Fa").run()
        app.selectbox(key="_result_widget_probe_manifold_B").select("Fb").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["rf_manifold_A"], "Fb")
        self.assertEqual(app.session_state["rf_manifold_B"], "Fa")
        self.assertEqual(app.session_state["probe_manifold_A"], "Total")
        self.assertEqual(app.session_state["probe_manifold_B"], "Fb")
        app.session_state["result_species_tab"] = "Alkali A — Rb87"
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="_result_widget_rf_manifold_A").value, "Fb")


if __name__ == "__main__":
    unittest.main()
