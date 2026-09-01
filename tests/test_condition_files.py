import unittest
from unittest.mock import patch

from alkali_pumping_app.ui import conditions
from alkali_pumping_app.version import CONDITION_SCHEMA_VERSION


def legacy_v5_conditions():
    return {
        "condition_name": "legacy",
        "atom_name": "Rb85",
        "gamma_ER": 7.0,
        "q_axis": "y",
        "bias_larmor_hz": 12.0,
        "temperature_C_for_table": 40.0,
        "n2_pressure_torr": 120.0,
        "include_spin_exchange": True,
        "D1_width": 17.8,
        "D2_width": 18.1,
        "D1_shift": -8.25,
        "D2_shift": -5.9,
        "line1": "D1",
        "transition1": "2→3",
        "det_rel1": 10.0,
        "intensity1": 3.0,
        "k1": "x",
        "pol1": "linear z",
        "line2": "D2",
        "transition2": "3→4",
        "det_rel2": -20.0,
        "intensity2": 4.0,
        "k2": "z",
        "pol2": "sigma+",
        "line3": "D1",
        "transition3": "3→3",
        "det_rel3": 0.0,
        "intensity3": 99.0,
        "k3": "z",
        "pol3": "sigma-",
        "rf_axis": "x",
        "rf_observable": "Fx",
        "rf_frequency_lower_hz": 0.0,
        "rf_frequency_upper_hz": 50.0,
        "rf_show_amplitude": True,
        "rf_show_in_phase": False,
        "rf_show_quadrature": False,
        "rf_relaxation_normalized": False,
        "show_allowed_only": True,
        "show_rate_matrices": False,
    }


def legacy_v60_conditions():
    excluded = {
        "static_field_axis", "static_field_nT", "q_axis_A", "q_axis_B",
        *conditions._pump_condition_keys("A3"),
        *conditions._pump_condition_keys("B3"),
    }
    values = {
        key: value
        for key, value in conditions.DEFAULT_STARTUP_CONDITION.items()
        if key not in excluded and not key.startswith("rf_")
    }
    values.update({
        "q_axis": "x",
        "bias_larmor_hz_A": 10.0,
        "rf_axis": "y",
        "rf_observable": "Fz",
        "rf_frequency_lower_hz": 2.0,
        "rf_frequency_upper_hz": 80.0,
        "rf_show_amplitude": False,
        "rf_show_in_phase": True,
        "rf_show_quadrature": False,
        "rf_relaxation_normalized": True,
        "rf_density_factor": True,
    })
    return values


class ConditionFileTests(unittest.TestCase):
    def test_condition_schema_is_v68(self):
        self.assertEqual(CONDITION_SCHEMA_VERSION, "6.8")
        self.assertEqual(conditions.CONDITION_SCHEMA_VERSION, "6.8")

    def test_payload_has_current_shape(self):
        payload = conditions.build_condition_payload(conditions.DEFAULT_STARTUP_CONDITION)
        self.assertEqual(payload["app"], "alkali_pumping")
        self.assertEqual(payload["format"], "alkali_pumping_conditions")
        self.assertEqual(payload["version"], "6.8")
        self.assertEqual(set(payload["conditions"]), set(conditions.CONDITION_KEYS))
        self.assertNotIn("include_spin_exchange", payload["conditions"])

    def test_default_atoms_and_pump_intensities(self):
        defaults = conditions.DEFAULT_STARTUP_CONDITION
        self.assertEqual(defaults["atom_A_name"], "Rb87")
        self.assertEqual(defaults["atom_B_name"], "None")
        self.assertEqual(defaults["gamma_ER_A"], 10.0)
        self.assertEqual(defaults["gamma_ER_B"], 10.0)
        self.assertEqual(defaults["temperature_C_for_table"], 23.0)
        intensities = [
            defaults[f"intensity_{prefix}"]
            for prefix in ("A1", "A2", "A3", "B1", "B2", "B3")
        ]
        self.assertEqual(intensities, [5.0, 5.0, 0.0, 0.0, 0.0, 0.0])
        self.assertTrue(all(value >= 0.0 for value in intensities))
        self.assertEqual(defaults["probe_source_A"], "PumpA2")
        self.assertEqual(defaults["probe_source_B"], "PumpB2")
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
                ValueError, "Expected 6.8, legacy 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1, 6.0, or 5.0"
            ):
                conditions.apply_loaded_condition_dict(payload)

    def test_missing_current_field_is_rejected(self):
        values = dict(conditions.DEFAULT_STARTUP_CONDITION)
        values.pop("rf_axis_A")
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.8",
            "conditions": values,
        }
        with patch.object(conditions.st, "session_state", {}):
            with self.assertRaisesRegex(ValueError, "rf_axis_A"):
                conditions.apply_loaded_condition_dict(payload)

    def test_v60_condition_migrates_field_and_independent_rf_controls(self):
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.0",
            "conditions": legacy_v60_conditions(),
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertEqual(session_state["static_field_axis"], "x")
        self.assertEqual(session_state["q_axis_A"], "x")
        self.assertEqual(session_state["q_axis_B"], "x")
        self.assertGreater(session_state["static_field_nT"], 0.0)
        self.assertEqual(session_state["rf_axis_A"], "y")
        self.assertEqual(session_state["rf_axis_B"], "y")
        self.assertEqual(session_state["rf_frequency_upper_hz_A"], 80.0)
        self.assertEqual(session_state["rf_frequency_upper_hz_B"], 80.0)
        self.assertEqual(session_state["intensity_A3"], 0.0)
        self.assertEqual(session_state["intensity_B3"], 0.0)

    def test_v61_condition_adds_zero_intensity_third_pumps(self):
        values = {
            key: value
            for key, value in conditions.DEFAULT_STARTUP_CONDITION.items()
            if key not in {
                *conditions._pump_condition_keys("A3"),
                *conditions._pump_condition_keys("B3"),
            }
        }
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.1",
            "conditions": values,
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertEqual(session_state["intensity_A3"], 0.0)
        self.assertEqual(session_state["intensity_B3"], 0.0)

    def test_v62_condition_drops_obsolete_spin_exchange_toggle(self):
        values = dict(conditions.DEFAULT_STARTUP_CONDITION)
        values["include_spin_exchange"] = False
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.2",
            "conditions": values,
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertNotIn("include_spin_exchange", session_state)

    def test_v63_condition_preserves_relative_concentration_value(self):
        values = dict(conditions.DEFAULT_STARTUP_CONDITION)
        values["density_mode"] = "Relative concentration"
        values["density_ratio_B_to_A"] = 0.375
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.3",
            "conditions": values,
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertEqual(session_state["density_mode"], "Relative concentration")
        self.assertEqual(session_state["density_ratio_B_to_A"], 0.375)

    def test_v64_condition_adds_disabled_pi_shift_controls(self):
        values = {
            key: value
            for key, value in conditions.DEFAULT_STARTUP_CONDITION.items()
            if "add_pi" not in key
        }
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.4",
            "conditions": values,
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertFalse(session_state["rf_add_pi_in_phase_A"])
        self.assertFalse(session_state["rf_add_pi_quadrature_A"])
        self.assertFalse(session_state["rf_add_pi_in_phase_B"])
        self.assertFalse(session_state["rf_add_pi_quadrature_B"])

    def test_v65_condition_adds_independent_probe_defaults(self):
        values = {
            key: value
            for key, value in conditions.DEFAULT_STARTUP_CONDITION.items()
            if key not in conditions.PROBE_CONDITION_KEYS
        }
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.5",
            "conditions": values,
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertEqual(session_state["probe_line_A"], "D1")
        self.assertEqual(session_state["probe_line_B"], "D1")
        self.assertEqual(session_state["probe_response_component_A"], "Total")
        self.assertEqual(session_state["probe_response_component_B"], "Total")

    def test_v66_condition_adds_custom_probe_source_defaults(self):
        values = {
            key: value
            for key, value in conditions.DEFAULT_STARTUP_CONDITION.items()
            if key not in ("probe_source_A", "probe_source_B")
        }
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.6",
            "conditions": values,
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertEqual(session_state["probe_source_A"], "Custom")
        self.assertEqual(session_state["probe_source_B"], "Custom")

    def test_v67_condition_migrates_probe_rank_switches(self):
        values = dict(conditions.DEFAULT_STARTUP_CONDITION)
        values.pop("probe_response_component_A")
        values.pop("probe_response_component_B")
        values.update({
            "probe_include_scalar_A": False,
            "probe_include_orientation_A": True,
            "probe_include_alignment_A": False,
            "probe_include_scalar_B": False,
            "probe_include_orientation_B": False,
            "probe_include_alignment_B": True,
        })
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "6.7",
            "conditions": values,
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertEqual(
            session_state["probe_response_component_A"], "Orientation induced"
        )
        self.assertEqual(
            session_state["probe_response_component_B"], "Alignment induced"
        )

    def test_v5_condition_migrates_legacy_third_pump_to_A3(self):
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "5.0",
            "conditions": legacy_v5_conditions(),
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertEqual(session_state["atom_A_name"], "Rb85")
        self.assertEqual(session_state["atom_B_name"], "None")
        self.assertEqual(session_state["intensity_A1"], 3.0)
        self.assertEqual(session_state["intensity_A2"], 4.0)
        self.assertEqual(session_state["intensity_A3"], 99.0)
        self.assertEqual(session_state["pol_A3"], "sigma-")
        self.assertEqual(session_state["intensity_B1"], 0.0)
        self.assertNotIn("intensity3", session_state)
        self.assertFalse(session_state["rf_density_factor_A"])
        self.assertFalse(session_state["rf_density_factor_B"])
        self.assertEqual(session_state["static_field_axis"], "y")
        self.assertEqual(session_state["q_axis_A"], "y")

    def test_legacy_rate_fields_are_staged_for_intensity_migration(self):
        values = legacy_v5_conditions()
        for number, rate in ((1, 1200.0), (2, 400.0)):
            values.pop(f"intensity{number}")
            values[f"rate{number}"] = rate
            values[f"rate_reference{number}"] = "At resonance"
        payload = {
            "app": "alkali_pumping",
            "format": "alkali_pumping_conditions",
            "version": "5.0",
            "conditions": values,
        }
        session_state = {}
        with patch.object(conditions.st, "session_state", session_state):
            conditions.apply_loaded_condition_dict(payload)
        self.assertEqual(session_state["_legacy_pump_inputs"]["A1"]["rate"], 1200.0)
        self.assertEqual(
            session_state["_legacy_pump_inputs"]["A2"]["rate_reference"],
            "At resonance",
        )


if __name__ == "__main__":
    unittest.main()
