import unittest

import numpy as np
import pandas as pd

from alkali_pumping_app.physics import (
    ATOMS,
    build_ground_states,
    coupled_weak_rf_matrix_readouts,
    coupled_weak_rf_observable_susceptibilities,
    quadrupole_axes_for_rf_axis,
    quadrupole_operator,
    rf_observable_display_label,
    weak_rf_observable_susceptibility,
)


def response_fixture(atom_name, frequencies, cross_rate):
    atom = ATOMS[atom_name]
    states = build_ground_states(atom)
    population = np.arange(1.0, len(states) + 1.0)
    population /= population.sum()
    zeeman = np.array([5.0 * float(state["m"]) for state in states])
    return {
        "atom": atom,
        "ground_states": states,
        "population": population,
        "nu_LS": np.zeros(len(states)),
        "df_pop": pd.DataFrame({
            "nu_B": zeeman,
            "G_OP": np.zeros(len(states)),
        }),
        "R_ER": 4.0,
        "R_SE_self": 2.0,
        "R_SE_cross": cross_rate,
        "rf_upper_F": max(float(state["F"]) for state in states),
        "q_axis": "z",
        "rf_axis": "x",
        "rf_observable": "Fx",
        "rf_frequencies_hz": np.asarray(frequencies, dtype=float),
        "light_shift_available": True,
    }


class CoupledRfTests(unittest.TestCase):
    def test_qij_axes_follow_the_rf_axis(self):
        self.assertEqual(quadrupole_axes_for_rf_axis("x"), ("y", "z"))
        self.assertEqual(quadrupole_axes_for_rf_axis("y"), ("z", "x"))
        self.assertEqual(quadrupole_axes_for_rf_axis("z"), ("x", "y"))
        self.assertEqual(rf_observable_display_label("Q_ij", "x"), "Q_yz")
        self.assertEqual(rf_observable_display_label("Q_ij", "y"), "Q_zx")
        self.assertEqual(rf_observable_display_label("Q_ij", "z"), "Q_xy")

    def test_off_diagonal_quadrupole_operator_is_hermitian(self):
        states = build_ground_states(ATOMS["Rb87"])
        upper_F = max(float(state["F"]) for state in states)
        operator = quadrupole_operator(
            states, "z", "y", "z", upper_F
        )
        np.testing.assert_allclose(operator, operator.conj().T, atol=1e-14)
        self.assertGreater(float(np.max(np.abs(operator))), 0.0)
        lower_indices = [
            index
            for index, state in enumerate(states)
            if not np.isclose(float(state["F"]), upper_F)
        ]
        np.testing.assert_allclose(operator[lower_indices, :], 0.0, atol=1e-14)
        np.testing.assert_allclose(operator[:, lower_indices], 0.0, atol=1e-14)

    def test_single_species_qij_response_is_finite_and_nonzero(self):
        states = build_ground_states(ATOMS["Rb87"])
        population = np.arange(1.0, len(states) + 1.0) ** 2
        population /= population.sum()
        state_index = {
            (float(state["F"]), float(state["m"])): index
            for index, state in enumerate(states)
        }
        transitions = np.full(len(states), np.nan)
        for index, state in enumerate(states):
            F = float(state["F"])
            m = float(state["m"])
            if (F, m - 1.0) in state_index:
                transitions[index] = 8.0 + 0.35 * m
        response = weak_rf_observable_susceptibility(
            frequencies_hz=np.linspace(0.0, 16.0, 33),
            ground_states=states,
            populations=population,
            adjacent_transition_hz=transitions,
            gamma_op=np.full(len(states), 1.5),
            gamma_er=np.full(len(states), 2.0),
            gamma_se=np.full(len(states), 3.0),
            q_axis="z",
            rf_axis="x",
            observable="Q_ij",
        )
        self.assertTrue(np.isfinite(response[0]).all())
        self.assertGreater(float(np.max(response[0])), 1e-12)
        self.assertEqual(response[3]["observable_label"], "Q_yz")

    def test_A_and_B_drives_use_independent_frequency_arrays(self):
        result_A = response_fixture("Rb87", [1.0, 2.0, 3.0], 3.0)
        result_B = response_fixture("Rb85", [5.0, 7.0], 4.0)
        response = coupled_weak_rf_observable_susceptibilities(result_A, result_B)
        self.assertEqual(response["A"][0].shape, (3,))
        self.assertEqual(response["B"][0].shape, (2,))
        self.assertTrue(response["A"][3]["coupled"])
        self.assertTrue(response["B"][3]["coupled"])
        self.assertTrue(np.all(np.isfinite(response["A"][0])))
        self.assertTrue(np.all(np.isfinite(response["B"][0])))

    def test_coupled_matrix_readouts_reuse_each_species_frequency_grid(self):
        result_A = response_fixture("Rb87", [1.0, 2.0, 3.0], 3.0)
        result_B = response_fixture("Rb85", [5.0, 7.0], 4.0)
        operator_A = quadrupole_operator(
            result_A["ground_states"], "z", "y", "z", result_A["rf_upper_F"]
        )
        operator_B = quadrupole_operator(
            result_B["ground_states"], "z", "y", "z", result_B["rf_upper_F"]
        )
        response = coupled_weak_rf_matrix_readouts(
            result_A, result_B, {"alignment": operator_A}, {"alignment": operator_B}
        )
        self.assertEqual(response["A"]["alignment"][0].shape, (3,))
        self.assertEqual(response["B"]["alignment"][0].shape, (2,))
        self.assertTrue(response["A"]["alignment"][3]["coupled"])
        self.assertGreater(float(np.max(response["A"]["alignment"][0])), 0.0)

    def test_cross_exchange_changes_species_local_response(self):
        uncoupled_A = response_fixture("Rb87", [1.0, 2.0], 0.0)
        uncoupled_B = response_fixture("Rb85", [1.0, 2.0], 0.0)
        coupled_A = response_fixture("Rb87", [1.0, 2.0], 8.0)
        coupled_B = response_fixture("Rb85", [1.0, 2.0], 9.0)
        response_uncoupled = coupled_weak_rf_observable_susceptibilities(
            uncoupled_A, uncoupled_B
        )
        response_coupled = coupled_weak_rf_observable_susceptibilities(
            coupled_A, coupled_B
        )
        self.assertGreater(
            float(np.max(np.abs(response_uncoupled["A"][0] - response_coupled["A"][0]))),
            1e-18,
        )
        self.assertGreater(
            float(np.max(np.abs(response_uncoupled["B"][0] - response_coupled["B"][0]))),
            1e-18,
        )

    def test_coupled_qij_response_uses_the_dynamic_tensor_readout(self):
        result_A = response_fixture("Rb87", [1.0, 5.0, 9.0], 3.0)
        result_B = response_fixture("Rb85", [1.0, 5.0, 9.0], 4.0)
        result_A["rf_observable"] = "Q_ij"
        result_B["rf_observable"] = "Q_ij"
        response = coupled_weak_rf_observable_susceptibilities(
            result_A, result_B
        )
        self.assertTrue(np.isfinite(response["A"][0]).all())
        self.assertTrue(np.isfinite(response["B"][0]).all())
        self.assertGreater(float(np.max(response["A"][0])), 1e-18)
        self.assertGreater(float(np.max(response["B"][0])), 1e-18)


if __name__ == "__main__":
    unittest.main()
