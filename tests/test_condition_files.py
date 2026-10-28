import unittest
from unittest.mock import patch

from alkali_pumping_app.ui import conditions
from alkali_pumping_app.version import CONDITION_SCHEMA_VERSION


class ConditionFileTests(unittest.TestCase):
    def test_manifold_settings_round_trip_independently(self):
        values = dict(conditions.DEFAULT_STARTUP_CONDITION)
        values.update(rf_manifold_A="Fb", rf_manifold_B="Both",
                      probe_manifold_A="Fa", probe_manifold_B="Fb",
                      probe_show_manifold_decomposition_A=True)
        payload = conditions.build_condition_payload(values)
        self.assertIn("two-manifold", payload["physics_model"])
        restored = {}
        with patch.object(conditions.st, "session_state", restored):
            conditions.apply_loaded_condition_dict(payload)
        for key in conditions.CONDITION_KEYS:
            self.assertEqual(restored[key], values[key])

    def test_missing_new_field_and_invalid_manifold_are_rejected_without_partial_load(self):
        for field, value in (("rf_manifold_A", "invalid"), ("probe_manifold_B", "Both"),
                             ("rf_manifold_A", None)):
            payload = conditions.build_condition_payload(conditions.DEFAULT_STARTUP_CONDITION)
            if value is None:
                payload["conditions"].pop(field)
            else:
                payload["conditions"][field] = value
            restored = {"condition_name": "unchanged"}
            with patch.object(conditions.st, "session_state", restored):
                with self.assertRaises(ValueError):
                    conditions.apply_loaded_condition_dict(payload)
            self.assertEqual(restored, {"condition_name": "unchanged"})

    def test_incomplete_current_file_is_not_filled_with_new_session_defaults(self):
        payload = conditions.build_condition_payload(conditions.DEFAULT_STARTUP_CONDITION)
        payload["conditions"].pop("static_field_nT")
        with patch.object(conditions.st, "session_state", {}):
            with self.assertRaisesRegex(ValueError, "static_field_nT"):
                conditions.apply_loaded_condition_dict(payload)

    def test_condition_schema_is_v615(self):
        self.assertEqual(CONDITION_SCHEMA_VERSION, "6.15")
        self.assertEqual(conditions.CONDITION_SCHEMA_VERSION, "6.15")

    def test_payload_has_current_shape(self):
        payload = conditions.build_condition_payload(conditions.DEFAULT_STARTUP_CONDITION)
        self.assertEqual(payload["app"], "alkali_pumping")
        self.assertEqual(payload["format"], "alkali_pumping_conditions")
        self.assertEqual(payload["version"], "6.15")
        self.assertEqual(set(payload["conditions"]), set(conditions.CONDITION_KEYS))
        self.assertNotIn("include_spin_exchange", payload["conditions"])
        self.assertNotIn("show_allowed_only", payload["conditions"])

    def test_default_atoms_and_pump_intensities(self):
        defaults = conditions.DEFAULT_STARTUP_CONDITION
        self.assertEqual(defaults["atom_A_name"], "Rb87")
        self.assertEqual(defaults["atom_B_name"], "None")
        self.assertEqual(defaults["gamma_ER_A"], 10.0)
        self.assertEqual(defaults["gamma_ER_B"], 10.0)
        self.assertEqual(defaults["gamma_u_A"], 1.0)
        self.assertEqual(defaults["gamma_u_B"], 1.0)
        self.assertEqual(defaults["temperature_C_for_table"], 23.0)
        intensities = [
            defaults[f"intensity_{prefix}"]
            for prefix in ("A1", "A2", "A3", "B1", "B2", "B3")
        ]
        self.assertEqual(intensities, [10.0, 10.0, 0.0, 10.0, 10.0, 0.0])
        self.assertTrue(all(value >= 0.0 for value in intensities))
        for label in ("A", "B"):
            self.assertEqual(defaults[f"transition_{label}1"], "2→2")
            self.assertEqual(defaults[f"det_rel_{label}1"], 450.0)
            self.assertEqual(defaults[f"transition_{label}2"], "1→2")
            self.assertEqual(defaults[f"det_rel_{label}2"], 0.0)
        self.assertEqual(defaults["probe_source_A"], "PumpA1 weak")
        self.assertEqual(defaults["probe_source_B"], "PumpB1 weak")
        self.assertEqual(defaults["rf_frequency_upper_hz_A"], 100.0)
        self.assertEqual(defaults["rf_frequency_upper_hz_B"], 100.0)

    def test_save_click_synchronizes_clean_filename_with_visible_name(self):
        session_state = {"condition_name": "  measurement-42.json  "}
        with patch.object(conditions.st, "session_state", session_state):
            save_name = conditions.sync_condition_save_name()
        self.assertEqual(save_name, "measurement-42")
        self.assertEqual(session_state["_condition_save_name"], "measurement-42")

    def test_unknown_version_is_rejected(self):
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "2.23",
            "conditions": conditions.DEFAULT_STARTUP_CONDITION,
        }
        with patch.object(conditions.st, "session_state", {}):
            with self.assertRaisesRegex(
                ValueError, "expected 6.15"
            ):
                conditions.apply_loaded_condition_dict(payload)

    def test_missing_current_field_is_rejected(self):
        values = dict(conditions.DEFAULT_STARTUP_CONDITION)
        values.pop("rf_axis_A")
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": CONDITION_SCHEMA_VERSION,
            "conditions": values,
        }
        with patch.object(conditions.st, "session_state", {}):
            with self.assertRaisesRegex(ValueError, "rf_axis_A"):
                conditions.apply_loaded_condition_dict(payload)

if __name__ == "__main__":
    unittest.main()
