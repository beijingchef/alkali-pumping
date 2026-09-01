import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from alkali_pumping_app.pages import light_shift as light_shift_page

from alkali_pumping_app.physics import (
    ATOMS,
    calculate_light_shift_sweep,
    polarization_ellipse_vector,
)
from alkali_pumping_app.pages.light_shift import (
    TRANSITION_MARKER_COLOR,
    TRANSITION_MARKER_OPACITY,
    _layered_line_chart,
    _scattering_chart,
    adjacent_transition_dataframe,
    coefficient_dataframe,
    eigenvalue_dataframe,
    fictitious_field_dataframe,
    light_shift_export_dataframe,
    scattering_dataframe,
    state_shift_dataframe,
)


class LightShiftPlotDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.detunings = np.array([-500.0, 500.0])
        cls.sweep = calculate_light_shift_sweep(
            atom=ATOMS["Rb87"],
            line="D1",
            detunings_MHz=cls.detunings,
            E_lab=polarization_ellipse_vector("z", 17.0, 12.0),
            k_axis="z",
            q_axis="z",
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
        )

    def test_component_table_contains_three_manifold_quantities(self):
        frame = coefficient_dataframe(self.sweep, self.detunings)
        self.assertEqual(
            set(frame["Component"]),
            {"Scalar shift", "Vector coefficient", "Tensor m=0 shift"},
        )
        self.assertEqual(set(frame["F"]), {"F=1", "F=2", "Δν_F"})
        self.assertEqual(len(frame), 2 * 2 * 3 + 2)

    def test_scalar_difference_is_upper_minus_lower_manifold(self):
        frame = coefficient_dataframe(self.sweep, self.detunings)
        scalar = frame[frame["Component"] == "Scalar shift"]
        upper = scalar[scalar["F"] == "F=2"]["Shift"].to_numpy()
        lower = scalar[scalar["F"] == "F=1"]["Shift"].to_numpy()
        difference = scalar[scalar["F"] == "Δν_F"]["Shift"].to_numpy()
        np.testing.assert_allclose(difference, upper - lower)

    def test_vector_coefficient_converts_to_signed_microgauss(self):
        coefficients = coefficient_dataframe(self.sweep, self.detunings)
        fields = fictitious_field_dataframe(
            self.sweep,
            self.detunings,
            gamma_hz_per_nT_by_F={1.0: -2.0, 2.0: 4.0},
        )
        for F, gamma in ((1.0, -2.0), (2.0, 4.0)):
            vector = coefficients[
                (coefficients["Component"] == "Vector coefficient")
                & (coefficients["F"] == f"F={F:g}")
            ]["Shift"].to_numpy()
            field = fields[fields["F"] == f"F={F:g}"]["Shift"].to_numpy()
            np.testing.assert_allclose(field, 10.0 * vector / gamma)

    def test_state_table_contains_all_components_and_states(self):
        frame = state_shift_dataframe(self.sweep, self.detunings)
        self.assertEqual(
            set(frame["Component"]),
            {"Scalar", "Vector", "Tensor", "Total diagonal"},
        )
        self.assertEqual(len(frame), 8 * 2 * 4)

    def test_eigen_transition_and_scattering_tables_have_expected_groups(self):
        eigen = eigenvalue_dataframe(self.sweep, self.detunings)
        transitions = adjacent_transition_dataframe(self.sweep, self.detunings)
        scattering = scattering_dataframe(self.sweep, self.detunings)
        self.assertEqual(set(eigen["F"]), {1.0, 2.0})
        self.assertTrue(set(eigen["Branch"]).issubset({f"EV{i}" for i in range(1, 6)}))
        self.assertIn("EV1", set(eigen["Branch"]))
        self.assertEqual(len(eigen), 8 * 2)
        self.assertEqual(len(transitions), (2 + 4) * 2)
        self.assertEqual(set(scattering["F"]), {"F=1", "F=2"})

    def test_plot_export_is_wide_and_places_units_below_headers(self):
        plotted = pd.DataFrame(
            {
                "Detuning (MHz)": [-1.0, 1.0, -1.0, 1.0],
                "F": ["F=1", "F=1", "F=2", "F=2"],
                "Component": ["Scalar shift"] * 4,
                "Shift": [1.0, 2.0, 3.0, 4.0],
                "Unit": ["Hz/(µW/cm²)"] * 4,
            }
        )
        scattering = pd.DataFrame(
            {
                "Detuning (MHz)": [-1.0, 1.0],
                "F": ["F=1", "F=1"],
                "Scattering rate": [5.0, 6.0],
            }
        )
        exported = light_shift_export_dataframe(
            plotted,
            normalization="Per intensity",
            view="Components",
            scattering=scattering,
        )
        self.assertEqual(exported.iloc[0, 0], "MHz")
        self.assertEqual(exported.iloc[0, 1], "Hz/(µW/cm²)")
        self.assertEqual(exported.iloc[0, 3], "s⁻¹/(µW/cm²)")
        self.assertEqual(exported.columns[0], "Laser detuning")
        self.assertEqual(exported.shape, (3, 4))
        np.testing.assert_allclose(exported.iloc[1:, 1].astype(float), [1.0, 2.0])

    def test_transition_center_rules_are_dark_and_scattering_legend_is_untitled(self):
        markers = pd.DataFrame(
            {"Detuning (MHz)": [0.0], "Transition": ["F=1 to F'=2"]}
        )
        shifts = pd.DataFrame(
            {"Detuning (MHz)": [-1.0, 1.0], "F": ["F=1", "F=1"], "Shift": [0.0, 1.0]}
        )
        line_spec = _layered_line_chart(
            shifts, markers, "Shift", color_field="F", color_title=None
        ).to_dict()
        transition_rule = next(
            layer for layer in line_spec["layer"] if layer.get("mark", {}).get("strokeDash")
        )
        self.assertEqual(transition_rule["mark"]["color"], TRANSITION_MARKER_COLOR)
        self.assertEqual(
            transition_rule["mark"]["opacity"], TRANSITION_MARKER_OPACITY
        )

        scattering = pd.DataFrame(
            {
                "Detuning (MHz)": [-1.0, 1.0],
                "F": ["F=1", "F=1"],
                "Scattering rate": [1.0, 2.0],
            }
        )
        scattering_spec = _scattering_chart(
            scattering, markers, "Scattering rate"
        ).to_dict()
        scattering_line = next(
            layer
            for layer in scattering_spec["layer"]
            if layer.get("mark", {}).get("type") == "line"
        )
        self.assertIsNone(scattering_line["encoding"]["color"]["title"])
        self.assertTrue(
            any(layer.get("mark", {}).get("strokeDash") for layer in scattering_spec["layer"])
        )

    def test_entering_ellipse_mode_copies_the_current_preset(self):
        state = {
            "ls_polarization_mode": "Preset",
            "_ls_widget_polarization_mode": "Ellipse",
            "ls_k_axis": "x",
            "ls_preset": "linear z",
            "ls_azimuth_deg": 12.0,
            "ls_ellipticity_deg": 8.0,
        }
        with patch.object(light_shift_page.st, "session_state", state):
            light_shift_page._store_light_shift_polarization_mode()
        self.assertEqual(state["ls_polarization_mode"], "Ellipse")
        self.assertEqual(state["ls_azimuth_deg"], 90.0)
        self.assertEqual(state["ls_ellipticity_deg"], 0.0)


if __name__ == "__main__":
    unittest.main()
