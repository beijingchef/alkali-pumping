"""Check shared temperature, editable collision inputs, and navigation state."""

from math import sqrt
from pathlib import Path
import json
import unittest

from streamlit.testing.v1 import AppTest


class AtomicPropertiesUiTests(unittest.TestCase):
    @staticmethod
    def edit_buffer_table(app, edited_rows):
        """Send a table-cell edit through the real data-editor callback.

        AppTest exposes editors as read-only Dataframe elements, so send the
        same widget-state JSON used by the browser until it has a cell-edit API.
        """
        table = app.dataframe[0]
        states = app._tree.get_widget_states()
        state = next((item for item in states.widgets if item.id == table.proto.id), None)
        if state is None:
            state = states.widgets.add(id=table.proto.id)
        state.string_value = json.dumps({
            "edited_rows": edited_rows, "added_rows": [], "deleted_rows": []
        })
        app._run(states)

    def test_temperature_pressure_and_self_collision_updates_survive_navigation(self):
        app = AppTest.from_file(
            str(Path(__file__).resolve().parents[1] / "alkali_pumping.py"),
            default_timeout=60,
        ).run()
        app.switch_page("alkali_pumping_app/pages/atomic_properties.py").run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.main.slider), 0)
        app.sidebar.slider(key="atomic_properties_temperature_C").set_value(25.0).run()

        def self_rate():
            return next(metric.value for metric in app.metric
                        if metric.label == "Self electron-randomization rate")

        cold_self_rate = self_rate()
        app.number_input(key="atomic_self_er_sigma_Rb87").set_value(3.2e-17).run()
        self.assertFalse(app.exception)
        self.assertAlmostEqual(float(self_rate().split()[0]) /
                               float(cold_self_rate.split()[0]), 2, places=3)

        app.session_state["atomic_properties_tab"] = "Buffer-gas collisions"
        app.run()
        self.assertEqual(len(app.number_input), 0)
        config = json.loads(app.dataframe[0].proto.columns)
        editable = {"Fill pressure at 25 °C (Torr)", "Electron-randomization cross section (cm²)"}
        for column in app.dataframe[0].value.columns:
            if column in editable:
                self.assertEqual(config[column]["type_config"]["min_value"], 0)
            else:
                self.assertTrue(config[column]["disabled"])
        self.edit_buffer_table(app, {0: {"Fill pressure at 25 °C (Torr)": 10.0}})
        self.assertFalse(app.exception)
        cold = app.dataframe[0].value.iloc[0]
        self.assertAlmostEqual(cold["Operating pressure (Torr)"], 10)
        self.assertAlmostEqual(cold["Broadening (MHz)"], 178)
        self.assertAlmostEqual(cold["Shift (MHz)"], -82.5)
        self.assertGreater(cold["Electron randomization (s⁻¹)"], 0)

        app.sidebar.slider(key="atomic_properties_temperature_C").set_value(100.0).run()
        self.assertFalse(app.exception)
        hot = app.dataframe[0].value.iloc[0]
        ratio = 373.15 / 298.15
        self.assertAlmostEqual(hot["Operating pressure (Torr)"], 10 * ratio)
        self.assertAlmostEqual(hot["Broadening (MHz)"], 178 * ratio)
        self.assertAlmostEqual(hot["Shift (MHz)"], -82.5 * ratio)
        self.assertAlmostEqual(hot["Electron randomization (s⁻¹)"] /
                               cold["Electron randomization (s⁻¹)"], sqrt(ratio))
        self.edit_buffer_table(app, {0: {"Fill pressure at 25 °C (Torr)": 20.0}})
        doubled = app.dataframe[0].value.iloc[0]
        for column in ("Broadening (MHz)", "Shift (MHz)", "Electron randomization (s⁻¹)"):
            self.assertAlmostEqual(doubled[column] / hot[column], 2)
        self.edit_buffer_table(app, {2: {
            "Electron-randomization cross section (cm²)": 1e-22,
            "Fill pressure at 25 °C (Torr)": 5.0,
        }})
        self.assertGreater(app.dataframe[0].value.iloc[2]["Electron randomization (s⁻¹)"], 0)

        # Pressure is shared; isotope-specific cross sections must stay separate.
        app.sidebar.selectbox(key="atomic_properties_isotope").select("Cs133").run()
        self.assertEqual(app.dataframe[0].value.iloc[0]["Fill pressure at 25 °C (Torr)"], 20)
        self.edit_buffer_table(app, {0: {"Fill pressure at 25 °C (Torr)": 30.0}})
        app.sidebar.selectbox(key="atomic_properties_isotope").select("Rb87").run()
        self.assertEqual(app.dataframe[0].value.iloc[0]["Fill pressure at 25 °C (Torr)"], 30)
        self.assertEqual(app.dataframe[0].value.iloc[2]["Electron-randomization cross section (cm²)"], 1e-22)
        self.edit_buffer_table(app, {0: {"Electron-randomization cross section (cm²)": 2e-22}})
        self.assertEqual(app.dataframe[0].value.iloc[0]["Fill pressure at 25 °C (Torr)"], 30)
        self.assertAlmostEqual(app.dataframe[0].value.iloc[0]["Electron randomization (s⁻¹)"] /
                               hot["Electron randomization (s⁻¹)"], 6)

        app.session_state["atomic_properties_tab"] = "Thermal & spin exchange"
        app.run()
        self.assertNotEqual(self_rate(), cold_self_rate)
        self.assertEqual(app.number_input(key="atomic_self_er_sigma_Rb87").value, 3.2e-17)
        app.switch_page("alkali_pumping_app/pages/ACC-Magnetometer.py").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="atom_A_name").value, "Rb87")
        app.switch_page("alkali_pumping_app/pages/atomic_properties.py").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.sidebar.slider(key="atomic_properties_temperature_C").value, 100)
        self.assertEqual(app.number_input(key="atomic_self_er_sigma_Rb87").value, 3.2e-17)
        app.session_state["atomic_properties_tab"] = "Buffer-gas collisions"
        app.run()
        self.assertEqual(app.dataframe[0].value.iloc[0]["Fill pressure at 25 °C (Torr)"], 30)
        self.assertEqual(app.dataframe[0].value.iloc[2]["Electron-randomization cross section (cm²)"], 1e-22)
        self.assertFalse(app.exception)


if __name__ == "__main__":
    unittest.main()
