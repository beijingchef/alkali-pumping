import unittest
from unittest.mock import patch

from alkali_pumping_app.physics.constants import DEFAULT_N2_COEFFS
from alkali_pumping_app.ui import atomic_polarizability_conditions as conditions


class AtomicPolarizabilityConditionTests(unittest.TestCase):
    def test_default_temperature_is_23_C(self):
        self.assertEqual(
            conditions.ATOMIC_POLARIZABILITY_DEFAULTS["temperature_C"], 23.0
        )

    def test_payload_has_independent_versioned_shape(self):
        payload = conditions.build_atomic_polarizability_payload(
            conditions.ATOMIC_POLARIZABILITY_DEFAULTS
        )
        self.assertEqual(payload["app"], "alkali_pumping")
        self.assertEqual(
            payload["format"],
            "alkali_pumping_atomic_polarizability_conditions",
        )
        self.assertEqual(payload["version"], "1.0")
        self.assertEqual(
            set(payload["conditions"]),
            set(conditions.ATOMIC_POLARIZABILITY_DEFAULTS),
        )

    def test_round_trip_restores_every_visible_setting(self):
        values = dict(conditions.ATOMIC_POLARIZABILITY_DEFAULTS)
        values.update(
            {
                "condition_name": "cesium-polarizability.json",
                "atom_name": "Cs133",
                "temperature_C": 68.0,
                "n2_pressure_torr": 120.0,
                "line": "D2",
                "lower_MHz": -7200.0,
                "upper_MHz": 8800.0,
                "points": 801,
                "plot_alpha_eq": True,
                "plot_alpha_gt": False,
            }
        )
        payload = conditions.build_atomic_polarizability_payload(values)
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            loaded_name = conditions.apply_atomic_polarizability_payload(payload)
        self.assertEqual(loaded_name, "cesium-polarizability")
        self.assertEqual(session_state["ap_atom_name"], "Cs133")
        self.assertEqual(session_state["ap_n2_pressure_torr"], 120.0)
        self.assertEqual(session_state["ap_line"], "D2")
        self.assertEqual(session_state["ap_lower_MHz"], -7200.0)
        self.assertTrue(session_state["ap_plot_alpha_eq"])
        self.assertFalse(session_state["ap_plot_alpha_gt"])
        self.assertEqual(
            session_state["ap_D2_width"],
            DEFAULT_N2_COEFFS["Cs133"]["D2"]["width"],
        )
        self.assertTrue(session_state["_ap_loaded_preserve_range"])

    def test_wrong_condition_format_is_rejected(self):
        payload = conditions.build_atomic_polarizability_payload(
            conditions.ATOMIC_POLARIZABILITY_DEFAULTS
        )
        payload["format"] = "alkali_pumping_light_shift_conditions"
        with patch.object(conditions.st, "session_state", {}):
            with self.assertRaisesRegex(ValueError, "atomic-polarizability condition"):
                conditions.apply_atomic_polarizability_payload(payload)

    def test_missing_field_is_rejected(self):
        payload = conditions.build_atomic_polarizability_payload(
            conditions.ATOMIC_POLARIZABILITY_DEFAULTS
        )
        payload["conditions"].pop("line")
        with patch.object(conditions.st, "session_state", {}):
            with self.assertRaisesRegex(ValueError, "line"):
                conditions.apply_atomic_polarizability_payload(payload)


if __name__ == "__main__":
    unittest.main()
