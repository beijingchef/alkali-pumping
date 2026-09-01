import unittest
from unittest.mock import patch

from alkali_pumping_app.physics.constants import DEFAULT_N2_COEFFS
from alkali_pumping_app.ui import light_shift_conditions as conditions


class LightShiftConditionTests(unittest.TestCase):
    def test_default_temperature_is_23_C(self):
        self.assertEqual(conditions.LIGHT_SHIFT_DEFAULTS["temperature_C"], 23.0)

    def test_payload_has_independent_versioned_shape(self):
        payload = conditions.build_light_shift_payload(
            conditions.LIGHT_SHIFT_DEFAULTS
        )
        self.assertEqual(payload["app"], "alkali_pumping")
        self.assertEqual(
            payload["format"], "alkali_pumping_light_shift_conditions"
        )
        self.assertEqual(payload["version"], "1.1")
        self.assertEqual(
            set(payload["conditions"]), set(conditions.LIGHT_SHIFT_DEFAULTS)
        )

    def test_round_trip_restores_every_light_shift_setting(self):
        values = dict(conditions.LIGHT_SHIFT_DEFAULTS)
        values.update(
            {
                "condition_name": "cesium-ellipse.json",
                "atom_name": "Cs133",
                "temperature_C": 72.0,
                "n2_pressure_torr": 90.0,
                "static_field_axis": "x",
                "static_field_nT": -35.0,
                "D1_width": DEFAULT_N2_COEFFS["Cs133"]["D1"]["width"],
                "polarization_mode": "Ellipse",
                "azimuth_deg": 27.0,
                "ellipticity_deg": -13.0,
                "view": "Eigenvalues",
                "show_scalar": False,
            }
        )
        payload = conditions.build_light_shift_payload(values)
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            loaded_name = conditions.apply_light_shift_payload(payload)
        self.assertEqual(loaded_name, "cesium-ellipse")
        self.assertEqual(session_state["ls_atom_name"], "Cs133")
        self.assertEqual(session_state["ls_static_field_nT"], -35.0)
        self.assertEqual(session_state["ls_polarization_mode"], "Ellipse")
        self.assertEqual(session_state["ls_ellipticity_deg"], -13.0)
        self.assertEqual(session_state["ls_view"], "Eigenvalues")
        self.assertFalse(session_state["ls_show_scalar"])
        self.assertTrue(session_state["_ls_loaded_preserve_range"])

    def test_v10_payload_defaults_to_showing_scalar_plot(self):
        payload = conditions.build_light_shift_payload(
            conditions.LIGHT_SHIFT_DEFAULTS
        )
        payload["version"] = "1.0"
        payload["conditions"].pop("show_scalar")
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_light_shift_payload(payload)
        self.assertTrue(session_state["ls_show_scalar"])

    def test_wrong_condition_format_is_rejected(self):
        payload = conditions.build_light_shift_payload(
            conditions.LIGHT_SHIFT_DEFAULTS
        )
        payload["format"] = "alkali_pumping_conditions"
        with patch.object(conditions.st, "session_state", {}):
            with self.assertRaisesRegex(ValueError, "light-shift condition"):
                conditions.apply_light_shift_payload(payload)

    def test_missing_field_is_rejected(self):
        payload = conditions.build_light_shift_payload(
            conditions.LIGHT_SHIFT_DEFAULTS
        )
        payload["conditions"].pop("ellipticity_deg")
        with patch.object(conditions.st, "session_state", {}):
            with self.assertRaisesRegex(ValueError, "ellipticity_deg"):
                conditions.apply_light_shift_payload(payload)


if __name__ == "__main__":
    unittest.main()
