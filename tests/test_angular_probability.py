"""Angular probability normalization and the periodic weak RF state."""

import unittest
from copy import deepcopy

import numpy as np

from alkali_pumping_app.physics.rf_response import (
    _spin_operator, build_rf_response_contexts, weak_rf_density_phasor,
    weak_rf_matrix_readouts,
)
from alkali_pumping_app.physics.spin_coherent import angular_probability_surface
from alkali_pumping_app.physics.system_model import compute_alkali_equilibrium
from tests.test_two_manifold_rf import system_inputs


class AngularProbabilityTests(unittest.TestCase):
    def test_coherent_probability_integrates_to_manifold_population(self):
        states = [{"F": 0.5, "m": -0.5}, {"F": 0.5, "m": 0.5}]
        density = np.array([[0.2, 0.04j], [-0.04j, 0.5]])
        theta, phi = np.meshgrid(np.linspace(0, np.pi, 181),
                                 np.linspace(0, 2 * np.pi, 361), indexing="ij")
        probability = angular_probability_surface(density, states, 0.5, theta, phi)
        integral = np.trapezoid(np.trapezoid(probability * np.sin(theta), phi[0], axis=1),
                                theta[:, 0])
        self.assertAlmostEqual(integral, 0.7, places=4)
        self.assertGreaterEqual(float(probability.min()), 0.0)

    def test_density_phasor_matches_existing_atomic_response(self):
        config, beams, common = system_inputs()
        physical = {key: value for key, value in config.items()
                    if key not in ("probe", "rf_axis", "rf_observable", "rf_frequencies_hz")}
        equilibrium = compute_alkali_equilibrium(physical, None, beams, common)
        result = dict(equilibrium["A"], rf_axis="x", rf_frequencies_hz=np.array([10.0]))
        context = build_rf_response_contexts(equilibrium["A"])["A"]
        F = max(float(state["F"]) for state in result["ground_states"])
        operator = _spin_operator(result["ground_states"], result["q_axis"], "x", F)
        response = weak_rf_matrix_readouts(result, {"Fx": operator}, context=context)["Fx"]
        phasor = weak_rf_density_phasor(result, 10.0, context=context)
        moment = np.trace(phasor @ operator)
        self.assertAlmostEqual(2 * moment.real, response[1][0], places=9)
        self.assertAlmostEqual(2 * moment.imag, response[2][0], places=9)
        population = np.asarray(result["population"])
        for phase in (0.0, 0.7, np.pi, 2 * np.pi):
            oscillation = phasor * np.exp(-1j * phase)
            density = np.diag(population) + 0.01 * (oscillation + oscillation.conj().T)
            np.testing.assert_allclose(density, density.conj().T, atol=1e-12)
            self.assertAlmostEqual(np.trace(density).real, 1.0, places=12)

    def test_dual_species_density_phasor_uses_coupled_modes(self):
        config, beams, common = system_inputs()
        other = deepcopy(config)
        other.update(label="B", atom_name="Rb85", density_cm3=2e10)
        physical = lambda item: {key: value for key, value in item.items()
                                 if key not in ("probe", "rf_axis", "rf_observable", "rf_frequencies_hz")}
        equilibrium = compute_alkali_equilibrium(physical(config), physical(other), beams, common)
        contexts = build_rf_response_contexts(
            equilibrium["A"], equilibrium["B"], equilibrium["J_coupled"]
        )
        for label in ("A", "B"):
            with self.subTest(label=label):
                result = dict(equilibrium[label], rf_axis="x", rf_frequencies_hz=np.array([10.0]))
                F = max(float(state["F"]) for state in result["ground_states"])
                operator = _spin_operator(result["ground_states"], result["q_axis"], "x", F)
                response = weak_rf_matrix_readouts(
                    result, {"Fx": operator}, context=contexts[label]
                )["Fx"]
                phasor = weak_rf_density_phasor(result, 10.0, context=contexts[label])
                moment = np.trace(phasor @ operator)
                self.assertAlmostEqual(2 * moment.real, response[1][0], places=9)
                self.assertAlmostEqual(2 * moment.imag, response[2][0], places=9)


if __name__ == "__main__":
    unittest.main()
