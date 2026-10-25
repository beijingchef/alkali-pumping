"""Check preset branches and compatibility of saved probe selections."""

from copy import deepcopy
import unittest
from unittest.mock import patch

from alkali_pumping_app.physics.constants import ATOMS, DEFAULT_N2_COEFFS
from alkali_pumping_app.physics.spectroscopy import transition_choice_labels
from alkali_pumping_app.ui import conditions
from alkali_pumping_app.ui.probe_presets import PROBE_PRESET_SOURCES, probe_preset_settings
from alkali_pumping_app.ui.time_evolution_conditions import (
    DEFAULTS, SCHEMA_VERSION, build_settings_payload, parse_settings_payload,
)


class ProbePresetTests(unittest.TestCase):
    def test_requested_branches_are_allowed_for_each_isotope_and_line(self):
        branches = {"Rb87": ("2→1", "1→2"), "Rb85": ("3→2", "2→3"),
                    "Cs133": ("4→3", "3→4"), "K39": ("2→1", "1→2")}
        self.assertEqual(set(branches), set(ATOMS))
        for name, (a, b) in branches.items():
            for source, line, transition, detuning in (
                ("Preset D1a", "D1", a, -2000.), ("Preset D1b", "D1", b, 2000.),
                ("Preset D2a", "D2", a, -2000.), ("Preset D2b", "D2", b, 2000.),
            ):
                with self.subTest(atom=name, source=source):
                    values = probe_preset_settings(source, ATOMS[name])
                    self.assertEqual(values, {
                        "line": line, "transition": transition, "det_rel": detuning,
                        "k": "z", "path_length_cm": 1., "azimuth_deg": 0.,
                        "ellipticity_deg": 0.,
                    })
                    self.assertIn(transition, transition_choice_labels(
                        ATOMS[name], line, 0., DEFAULT_N2_COEFFS[name], allowed_only=True))

    def test_other_sources_do_not_override_custom_or_pump_settings(self):
        for source in ("Custom", "PumpA1 weak", "PumpA1", "PumpB3 weak"):
            self.assertEqual(probe_preset_settings(source, ATOMS["Rb87"]), {})

    def test_edited_presets_round_trip_in_time_evolution_settings(self):
        for source in PROBE_PRESET_SOURCES:
            with self.subTest(source=source):
                values = deepcopy(DEFAULTS)
                values["atom_B_name"] = "Rb85"
                for label, name in (("A", "Rb87"), ("B", "Rb85")):
                    values[f"probe_source_{label}"] = source
                    for field, value in probe_preset_settings(source, ATOMS[name]).items():
                        values[f"probe_{field}_{label}"] = value
                values.update(probe_det_rel_A=-1234., probe_path_length_cm_B=2.5,
                              probe_azimuth_deg_A=20., probe_ellipticity_deg_B=-10.)
                payload = build_settings_payload(values)
                self.assertEqual(payload["version"], SCHEMA_VERSION)
                self.assertEqual(parse_settings_payload(payload), values)

    def test_current_formats_preserve_custom_optics_and_uniform_rates(self):
        for page in ("acc", "time"):
            values = deepcopy(conditions.DEFAULT_STARTUP_CONDITION if page == "acc" else DEFAULTS)
            values.update(probe_source_A="Custom", probe_path_length_cm_A=2.7,
                          probe_line_A="D2", probe_transition_A="2→3",
                          probe_det_rel_A=123., probe_azimuth_deg_A=25.,
                          probe_ellipticity_deg_A=-10., gamma_u_A=7.5)
            if page == "acc":
                payload = conditions.build_condition_payload(values)
                state = {}
                with patch.object(conditions.st, "session_state", state):
                    conditions.apply_loaded_condition_dict(payload)
                restored = {key: state[key] for key in conditions.CONDITION_KEYS}
            else:
                payload = build_settings_payload(values)
                restored = parse_settings_payload(payload)
            self.assertEqual(restored, values)


if __name__ == "__main__":
    unittest.main()
