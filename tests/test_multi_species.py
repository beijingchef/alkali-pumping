import unittest

import numpy as np

from alkali_pumping_app.physics import (
    ATOMS,
    alkali_vapor_density_cm3,
    build_ER_matrix,
    build_cross_spin_exchange_matrix,
    build_ground_states,
    coupled_population_jacobian,
    electron_marginal_from_population,
    resolve_alkali_densities,
    steady_state_two_species,
)


class MultiSpeciesTests(unittest.TestCase):
    def setUp(self):
        self.atom_A = ATOMS["Rb87"]
        self.atom_B = ATOMS["Rb85"]
        self.states_A = build_ground_states(self.atom_A)
        self.states_B = build_ground_states(self.atom_B)
        self.uniform_A = np.ones(len(self.states_A)) / len(self.states_A)
        self.uniform_B = np.ones(len(self.states_B)) / len(self.states_B)

    def test_cross_exchange_map_preserves_target_population(self):
        electron_B = electron_marginal_from_population(
            self.atom_B, self.states_B, self.uniform_B
        )
        matrix = build_cross_spin_exchange_matrix(
            self.atom_A, self.states_A, electron_B
        )
        np.testing.assert_allclose(matrix.sum(axis=0), 1.0, atol=1e-12)
        self.assertTrue(np.all(matrix >= -1e-14))

    def test_coupled_solver_normalizes_both_species(self):
        linear_A = 2.0 * (build_ER_matrix(self.atom_A, self.states_A) - np.eye(len(self.states_A)))
        linear_B = 3.0 * (build_ER_matrix(self.atom_B, self.states_B) - np.eye(len(self.states_B)))
        p_A, p_B, info = steady_state_two_species(
            linear_A, self.atom_A, self.states_A, 4.0, 7.0,
            linear_B, self.atom_B, self.states_B, 5.0, 6.0,
        )
        self.assertAlmostEqual(float(p_A.sum()), 1.0)
        self.assertAlmostEqual(float(p_B.sum()), 1.0)
        self.assertTrue(np.all(p_A >= 0.0))
        self.assertTrue(np.all(p_B >= 0.0))
        self.assertTrue(info["converged"])

    def test_full_jacobian_has_conserving_blocks(self):
        linear_A = build_ER_matrix(self.atom_A, self.states_A) - np.eye(len(self.states_A))
        linear_B = build_ER_matrix(self.atom_B, self.states_B) - np.eye(len(self.states_B))
        full, blocks = coupled_population_jacobian(
            linear_A, self.atom_A, self.states_A, self.uniform_A, 2.0, 3.0,
            linear_B, self.atom_B, self.states_B, self.uniform_B, 4.0, 5.0,
        )
        expected_size = len(self.states_A) + len(self.states_B)
        self.assertEqual(full.shape, (expected_size, expected_size))
        for block in blocks.values():
            np.testing.assert_allclose(block.sum(axis=0), 0.0, atol=1e-11)

    def test_density_modes_and_inactive_B(self):
        temperature_C = 60.0
        liquid_ratio = 0.25
        density_A, density_B = resolve_alkali_densities(
            "Rb87", "Cs133", temperature_C, "Relative concentration", liquid_ratio
        )
        saturated_A = alkali_vapor_density_cm3("Rb87", temperature_C)
        saturated_B = alkali_vapor_density_cm3("Cs133", temperature_C)
        self.assertAlmostEqual(density_A / saturated_A, 0.8)
        self.assertAlmostEqual(density_B / saturated_B, 0.2)
        self.assertAlmostEqual(
            density_B / density_A,
            liquid_ratio * saturated_B / saturated_A,
        )

        independent_A, independent_B = resolve_alkali_densities(
            "Rb87", "Cs133", temperature_C,
            "Independent saturated-vapor curves", liquid_ratio,
        )
        self.assertAlmostEqual(independent_A, saturated_A)
        self.assertAlmostEqual(independent_B, saturated_B)

        _, same_density_B = resolve_alkali_densities(
            "Rb87", "Rb87", temperature_C, "Relative concentration", 1.0
        )
        _, none_density_B = resolve_alkali_densities(
            "Rb87", "None", temperature_C, "Relative concentration", 1.0
        )
        self.assertEqual(same_density_B, 0.0)
        self.assertEqual(none_density_B, 0.0)

    def test_zero_liquid_B_ratio_recovers_pure_A_vapor(self):
        temperature_C = 80.0
        density_A, density_B = resolve_alkali_densities(
            "Rb87", "Cs133", temperature_C, "Relative concentration", 0.0
        )
        self.assertAlmostEqual(
            density_A, alkali_vapor_density_cm3("Rb87", temperature_C)
        )
        self.assertEqual(density_B, 0.0)


if __name__ == "__main__":
    unittest.main()
