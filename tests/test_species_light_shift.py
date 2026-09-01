import unittest

import numpy as np

from alkali_pumping_app.physics import ATOMS, compute_alkali_system
from alkali_pumping_app.physics.spectroscopy import line_center_frequency_MHz


class SpeciesLightShiftTests(unittest.TestCase):
    @staticmethod
    def _species_config(label, atom_name, q_axis):
        return {
            "label": label,
            "atom_name": atom_name,
            "density_cm3": 1.0e10,
            "R_ER": 4.0,
            "n2_coeffs": {
                "D1": {"width": 17.8, "shift": -4.52},
                "D2": {"width": 18.0, "shift": -4.52},
            },
            "q_axis": q_axis,
            "rf_axis": "y",
            "rf_observable": "Fz",
            "rf_frequencies_hz": np.array([10.0]),
        }

    @staticmethod
    def _beam(name, target_label, atom_name, axis):
        return {
            "name": name,
            "target_label": target_label,
            "absolute_frequency_MHz": line_center_frequency_MHz(
                ATOMS[atom_name], "D1"
            ),
            "intensity": 5.0,
            "k_axis": axis,
            "pol": "sigma+",
            "selected_transition": None,
            "transition_label": "D1 test pump",
        }

    def test_each_species_ignores_other_species_pump_light_shift(self):
        config_A = self._species_config("A", "Rb87", "z")
        config_B = self._species_config("B", "Cs133", "x")
        beams = [
            self._beam("PumpA1", "A", "Rb87", "z"),
            self._beam("PumpB1", "B", "Cs133", "x"),
        ]
        common = {
            "temperature_C": 23.5,
            "n2_pressure_torr": 0.0,
            # A stale caller flag must not disable the now-unconditional model.
            "include_spin_exchange": False,
            "static_field_axis": "z",
            "static_field_nT": 0.0,
        }

        system = compute_alkali_system(config_A, config_B, beams, common)

        for label in ("A", "B"):
            result = system[label]
            self.assertTrue(result["light_shift_available"])
            self.assertTrue(np.isfinite(result["nu_LS"]).all())
            self.assertEqual(
                [beam["target_label"] for beam, _info in result["light_shift_diagnostics"]],
                [label],
            )
            # Both beams remain in the general optical diagnostics; filtering
            # is intentionally limited to the reported AC-Stark shift.
            self.assertEqual(len(result["diagnostics"]), 2)
            self.assertGreater(result["R_SE_self"], 0.0)
            self.assertGreater(result["R_SE_cross"], 0.0)


if __name__ == "__main__":
    unittest.main()
