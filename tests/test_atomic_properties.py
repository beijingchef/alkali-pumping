import unittest
from math import pi, sqrt

from alkali_pumping_app.physics.atomic_properties import (
    ATOMIC_PROPERTY_DATA,
    alkali_thermal_properties,
    allowed_hyperfine_F,
    buffer_gas_collision_rate_s,
    buffer_gas_properties,
    format_transition_strength_fraction,
    format_transition_strength_vertical_fraction,
    grotrian_transitions,
    notebook_transition_strength,
)


class AtomicPropertiesTests(unittest.TestCase):
    def test_requested_isotopes_are_supported(self):
        self.assertEqual(
            list(ATOMIC_PROPERTY_DATA),
            ["Na23", "K39", "K41", "Rb85", "Rb87", "Cs133"],
        )

    def test_vapor_density_and_self_exchange_rate_increase_with_temperature(self):
        cold = alkali_thermal_properties("Rb87", 20.0)
        hot = alkali_thermal_properties("Rb87", 100.0)
        self.assertGreater(hot["density_cm3"], cold["density_cm3"])
        self.assertGreater(hot["spin_exchange_rate_s"], cold["spin_exchange_rate_s"])
        self.assertGreater(cold["rms_velocity_m_s"], 0.0)
        self.assertGreater(cold["mean_relative_velocity_m_s"], 0.0)

    def test_buffer_gas_rate_is_linear_in_pressure_and_cross_section(self):
        first = buffer_gas_collision_rate_s("K39", "N2", 50.0, 10.0, 1e-22)
        second = buffer_gas_collision_rate_s("K39", "N2", 50.0, 20.0, 2e-22)
        self.assertGreater(first, 0.0)
        self.assertAlmostEqual(second / first, 4.0, places=12)

    def test_self_randomization_uses_binary_relative_speed_and_saturated_density(self):
        for isotope in ATOMIC_PROPERTY_DATA:
            with self.subTest(isotope=isotope):
                values = alkali_thermal_properties(isotope, 100.0, 2e-17)
                mass_kg = ATOMIC_PROPERTY_DATA[isotope]["mass_amu"] * 1.66053906660e-27
                speed_cm_s = 100 * sqrt(16 * 1.380649e-23 * 373.15 / (pi * mass_kg))
                expected = values["density_cm3"] * 2e-17 * speed_cm_s
                self.assertAlmostEqual(values["self_er_rate_s"] / expected, 1.0, places=12)
                doubled = alkali_thermal_properties(isotope, 100.0, 4e-17)
                self.assertAlmostEqual(doubled["self_er_rate_s"] / expected, 2.0, places=12)
                self.assertEqual(alkali_thermal_properties(isotope, 100.0, 0)["self_er_rate_s"], 0)

    def test_self_randomization_defaults_and_missing_data(self):
        self.assertEqual(alkali_thermal_properties("Rb87", 25)["self_er_cross_section_cm2"], 1.6e-17)
        self.assertIsNone(alkali_thermal_properties("Na23", 25)["self_er_rate_s"])
        self.assertGreater(alkali_thermal_properties("Na23", 25, 1e-18)["self_er_rate_s"], 0)
        self.assertGreater(alkali_thermal_properties("Rb87", 100)["self_er_rate_s"],
                           alkali_thermal_properties("Rb87", 25)["self_er_rate_s"])

    def test_sealed_buffer_gas_conserves_density_and_has_sqrt_temperature_rate(self):
        cold = buffer_gas_properties("Rb87", "N2", "D1", 25, 10, 1e-22)
        hot = buffer_gas_properties("Rb87", "N2", "D1", 323.15, 10, 1e-22)
        self.assertAlmostEqual(cold["operating_pressure_torr"], 10)
        self.assertAlmostEqual(hot["operating_pressure_torr"], 20)
        self.assertEqual(hot["density_cm3"], cold["density_cm3"])
        self.assertAlmostEqual(hot["electron_randomization_rate_s"] /
                               cold["electron_randomization_rate_s"], sqrt(2), places=12)
        self.assertAlmostEqual(cold["broadening_MHz"], 178)
        self.assertAlmostEqual(cold["shift_MHz"], -82.5)
        self.assertAlmostEqual(hot["broadening_MHz"], 356)
        self.assertAlmostEqual(hot["shift_MHz"], -165)

    def test_buffer_pressure_and_cross_section_changes_and_zero_pressure(self):
        first = buffer_gas_properties("K39", "He4", "D2", 100, 10, 1e-22)
        doubled = buffer_gas_properties("K39", "He4", "D2", 100, 20, 2e-22)
        for key in ("broadening_MHz", "shift_MHz"):
            self.assertAlmostEqual(doubled[key] / first[key], 2)
        self.assertAlmostEqual(doubled["electron_randomization_rate_s"] /
                               first["electron_randomization_rate_s"], 4)
        zero = buffer_gas_properties("K39", "He4", "D2", 100, 0, 1e-22)
        for key in ("operating_pressure_torr", "density_cm3", "broadening_MHz",
                    "shift_MHz", "electron_randomization_rate_s"):
            self.assertEqual(zero[key], 0)
        missing = buffer_gas_properties("Na23", "CH4", "D1", 100, 10, 1e-22)
        self.assertIsNone(missing["broadening_MHz"])
        self.assertIsNone(missing["shift_MHz"])
        self.assertGreater(missing["electron_randomization_rate_s"], 0)

    def test_invalid_collision_inputs_are_rejected(self):
        for temperature in (-273.15, float("nan"), float("inf")):
            with self.subTest(temperature=temperature), self.assertRaises(ValueError):
                buffer_gas_properties("Rb87", "N2", "D1", temperature, 10, 1e-22)
        for value in (-1, float("nan"), float("inf")):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    buffer_gas_properties("Rb87", "N2", "D1", 25, value, 1e-22)
                with self.assertRaises(ValueError):
                    buffer_gas_properties("Rb87", "N2", "D1", 25, 10, value)
                with self.assertRaises(ValueError):
                    alkali_thermal_properties("Rb87", 25, value)

    def test_hyperfine_manifolds_match_angular_momentum_addition(self):
        self.assertEqual(allowed_hyperfine_F(1.5, 0.5), [1.0, 2.0])
        self.assertEqual(allowed_hyperfine_F(1.5, 1.5), [0.0, 1.0, 2.0, 3.0])
        self.assertEqual(allowed_hyperfine_F(3.5, 1.5), [2.0, 3.0, 4.0, 5.0])

    def test_notebook_strength_obeys_zeeman_selection_rule(self):
        allowed = notebook_transition_strength(1.5, 0.5, 1.0, -1.0, 2.0, 0.0)
        forbidden = notebook_transition_strength(1.5, 0.5, 1.0, -1.0, 2.0, 1.0)
        self.assertGreater(allowed, 0.0)
        self.assertEqual(forbidden, 0.0)

    def test_grotrian_transition_list_filters_polarization(self):
        transitions = grotrian_transitions(1.5, "D1", 1.0, 2.0, polarizations=(1,))
        self.assertTrue(transitions)
        self.assertTrue(all(row["q"] == 1 for row in transitions))
        self.assertTrue(all(row["strength"] > 0.0 for row in transitions))

    def test_transition_strengths_are_formatted_as_reduced_fractions(self):
        self.assertEqual(format_transition_strength_fraction(1.0 / 12.0), "1/12")
        self.assertEqual(format_transition_strength_fraction(2.0 / 6.0), "1/3")
        self.assertEqual(format_transition_strength_fraction(1.0), "1")

    def test_transition_strengths_support_vertical_fraction_mathtext(self):
        self.assertEqual(
            format_transition_strength_vertical_fraction(1.0 / 12.0),
            r"$\frac{1}{12}$",
        )
        self.assertEqual(format_transition_strength_vertical_fraction(1.0), "$1$")


if __name__ == "__main__":
    unittest.main()
