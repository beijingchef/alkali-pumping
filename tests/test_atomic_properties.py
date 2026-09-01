import unittest

from alkali_pumping_app.physics.atomic_properties import (
    ATOMIC_PROPERTY_DATA,
    alkali_thermal_properties,
    allowed_hyperfine_F,
    buffer_gas_collision_rate_s,
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
