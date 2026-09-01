import unittest

import numpy as np

from alkali_pumping_app.physics.atomic_polarizability import (
    calculate_atomic_polarizability_sweep,
    default_polarizability_sweep_range_MHz,
    oscillator_strength_from_natural_width,
)
from alkali_pumping_app.physics.constants import ATOMS, DEFAULT_N2_COEFFS
from alkali_pumping_app.physics.spectroscopy import hyperfine_transition_choices


class AtomicPolarizabilityTests(unittest.TestCase):
    @staticmethod
    def _sweep(atom=None, line="D1", detunings=None):
        atom = atom or ATOMS["Rb87"]
        detunings = (
            np.linspace(-10000.0, 10000.0, 101)
            if detunings is None
            else np.asarray(detunings, dtype=float)
        )
        coefficients = DEFAULT_N2_COEFFS["Rb87"][line]
        return calculate_atomic_polarizability_sweep(
            atom=atom,
            line=line,
            detunings_MHz=detunings,
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=coefficients["width"],
            n2_shift_MHz_per_torr=coefficients["shift"],
        )

    def test_rb_oscillator_strengths_match_known_d_line_scale(self):
        atom = ATOMS["Rb87"]
        self.assertAlmostEqual(
            oscillator_strength_from_natural_width(atom, "D1"), 0.3423, places=3
        )
        self.assertAlmostEqual(
            oscillator_strength_from_natural_width(atom, "D2"), 0.6962, places=3
        )

    def test_sweep_returns_all_four_finite_complex_responses(self):
        sweep = self._sweep()
        self.assertEqual(sweep["alpha_eq"].shape, (101,))
        self.assertEqual(sweep["alpha_hfs"].shape, (101,))
        self.assertEqual(set(sweep["alpha_gt"]), {1.0, 2.0})
        self.assertEqual(set(sweep["alpha_br"]), {1.0, 2.0})
        for values in (
            sweep["alpha_eq"],
            sweep["alpha_hfs"],
            *sweep["alpha_gt"].values(),
            *sweep["alpha_br"].values(),
        ):
            self.assertTrue(np.iscomplexobj(values))
            self.assertTrue(np.isfinite(values).all())

    def test_absorptive_wing_decays_faster_than_dispersive_wing(self):
        sweep = self._sweep(detunings=[50000.0, 100000.0])
        real_ratio = abs(sweep["alpha_eq"].real[1] / sweep["alpha_eq"].real[0])
        imaginary_ratio = abs(
            sweep["alpha_eq"].imag[1] / sweep["alpha_eq"].imag[0]
        )
        self.assertLess(imaginary_ratio, real_ratio)

    def test_transition_centers_include_pressure_shift(self):
        atom = ATOMS["Rb87"]
        common = dict(
            atom=atom,
            line="D1",
            detunings_MHz=np.array([-1.0, 1.0]),
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
        )
        zero = calculate_atomic_polarizability_sweep(
            n2_pressure_torr=0.0, **common
        )
        shifted = calculate_atomic_polarizability_sweep(
            n2_pressure_torr=10.0, **common
        )
        for transition, center in zero["transition_centers_MHz"].items():
            self.assertAlmostEqual(
                shifted["transition_centers_MHz"][transition] - center,
                -82.5,
            )

    def test_default_range_has_1500_MHz_outward_rounded_margins(self):
        atom = ATOMS["Rb87"]
        coefficients = DEFAULT_N2_COEFFS["Rb87"]
        lower, upper = default_polarizability_sweep_range_MHz(
            atom,
            "D1",
            25.0,
            coefficients,
        )
        transitions = hyperfine_transition_choices(
            atom,
            "D1",
            25.0,
            coefficients,
            allowed_only=True,
        )
        minimum = min(row["detP"] for row in transitions)
        maximum = max(row["detP"] for row in transitions)
        self.assertEqual(lower % 100.0, 0.0)
        self.assertEqual(upper % 100.0, 0.0)
        self.assertGreaterEqual(minimum - lower, 1500.0)
        self.assertLess(minimum - lower, 1600.0)
        self.assertGreaterEqual(upper - maximum, 1500.0)
        self.assertLess(upper - maximum, 1600.0)

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            calculate_atomic_polarizability_sweep(
                atom=ATOMS["Rb87"],
                line="D3",
                detunings_MHz=np.array([-1.0, 1.0]),
                n2_pressure_torr=0.0,
                temperature_C=23.5,
                n2_width_MHz_per_torr=17.8,
                n2_shift_MHz_per_torr=-8.25,
            )
        with self.assertRaises(ValueError):
            calculate_atomic_polarizability_sweep(
                atom=ATOMS["Rb87"],
                line="D1",
                detunings_MHz=np.array([0.0]),
                n2_pressure_torr=0.0,
                temperature_C=23.5,
                n2_width_MHz_per_torr=17.8,
                n2_shift_MHz_per_torr=-8.25,
            )


if __name__ == "__main__":
    unittest.main()
