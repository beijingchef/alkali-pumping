"""Analytic limits and end-to-end propagation of uniform ground relaxation."""

from copy import deepcopy
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import expm

from alkali_pumping_app.physics import (
    atomic_rf_readout_operators, build_rf_response_contexts, compute_alkali_system,
    largest_abs_Dm_relaxation_reference, magnetic_rf_operator, weak_rf_matrix_readouts,
    weak_rf_observable_susceptibility,
)
from alkali_pumping_app.physics.observables import add_total_relaxation_columns
from alkali_pumping_app.physics.rf_response import _single_species_coherence_generator
from alkali_pumping_app.physics.time_evolution import compute_population_evolution
from alkali_pumping_app.physics.uniform_relaxation import uniform_relaxation_generator
from alkali_pumping_app.ui import conditions
from alkali_pumping_app.ui.time_evolution_conditions import (
    DEFAULTS, build_settings_payload, parse_settings_payload,
)
from tests.test_two_manifold_rf import phasor, synthetic_result, system_inputs


class UniformRelaxationPhysicsTests(unittest.TestCase):
    def test_generator_preserves_trace_and_has_exact_uniform_decay(self):
        for size in (8, 12, 16):
            rate = 2.3
            generator = uniform_relaxation_generator(size, rate)
            p = np.arange(1., size + 1)
            p /= p.sum()
            np.testing.assert_allclose(generator.sum(axis=0), 0, atol=1e-15)
            np.testing.assert_allclose(generator @ np.ones(size), 0, atol=1e-15)
            for time in (0., .1, 1.):
                expected = 1 / size + (p - 1 / size) * np.exp(-rate * time)
                np.testing.assert_allclose(expm(generator * time) @ p, expected, atol=1e-14)
            np.testing.assert_array_equal(uniform_relaxation_generator(size, 0.), np.zeros((size, size)))
        for rate in (-1., np.nan, np.inf):
            with self.subTest(rate=rate), self.assertRaises(ValueError):
                uniform_relaxation_generator(8, rate)

    def test_dark_population_and_all_longitudinal_moments_decay_at_Ru(self):
        config, beams, common = system_inputs()
        config.update(density_cm3=0., R_ER=0., R_u=3.)
        beams[0].update(k_axis="z", pol="sigma+")
        result = compute_population_evolution(config, None, beams, common, .05, 1., 61)
        dark = result["time_s"] >= .05
        times = result["time_s"][dark] - .05
        populations = result["A"]["population"][dark]
        expected = 1 / 8 + (populations[0] - 1 / 8)[None, :] * np.exp(-3 * times[:, None])
        np.testing.assert_allclose(populations, expected, atol=3e-8)
        for name, values in result["A"]["moments"].items():
            if name.startswith("population_"):
                continue
            np.testing.assert_allclose(values[dark], values[dark][0] * np.exp(-3 * times), atol=3e-8)

    def test_population_steady_state_and_transient_long_time_limit_agree(self):
        config, beams, common = system_inputs()
        config["R_u"] = 4.
        beams[0].update(k_axis="z", pol="sigma+")
        steady = compute_alkali_system(config, None, beams, common)["A"]
        evolution = compute_population_evolution(config, None, beams, common, 8., 8., 31)
        np.testing.assert_allclose(evolution["A"]["population"][-1], steady["population"], atol=3e-8)
        zero = deepcopy(config)
        zero["R_u"] = 0.
        previous = compute_alkali_system(zero, None, beams, common)["A"]
        self.assertGreater(np.max(np.abs(steady["population"] - previous["population"])), 1e-5)
        np.testing.assert_allclose(steady["L_uniform"].sum(axis=0), 0, atol=1e-15)

    def test_uniform_relaxation_damps_every_secular_coherence(self):
        result = synthetic_result()
        result.update(R_ER=2., R_SE_self=3.)
        pairs, old = _single_species_coherence_generator(result)
        result["R_u"] = 2.7
        _, new = _single_species_coherence_generator(result)
        np.testing.assert_allclose(new - old, -2.7 * np.eye(len(pairs)), atol=1e-13)

    def test_two_manifold_rf_matches_analytic_uniform_linewidth(self):
        result = synthetic_result()
        result["R_u"] = 3.2
        actual = weak_rf_matrix_readouts(result, atomic_rf_readout_operators(result))
        _, info = magnetic_rf_operator(result)
        states = result["ground_states"]
        shifts = result["df_pop"]["nu_B"].to_numpy()
        transition = np.full(len(states), np.nan)
        index = {(state["F"], state["m"]): i for i, state in enumerate(states)}
        for i, state in enumerate(states):
            lower = index.get((state["F"], state["m"] - 1))
            if lower is not None:
                transition[i] = shifts[i] - shifts[lower]
        for manifold, F in (("Fa", 2.), ("Fb", 1.)):
            ratio = info["drive_ratios"][manifold]
            expected = weak_rf_observable_susceptibility(
                result["rf_frequencies_hz"], states, result["population"], transition,
                np.full(8, 4.), np.zeros(8), np.zeros(8), "z", "x", "Fx", target_F=F, gamma_u=3.2,
            )
            np.testing.assert_allclose(phasor(actual[manifold]), ratio * phasor(expected), atol=1e-12)

    def test_coupled_coherence_blocks_use_independent_uniform_rates(self):
        a, b = synthetic_result(), synthetic_result("Rb85")
        a.update(R_ER=2., R_SE_self=3., R_SE_cross=1.)
        b.update(R_ER=4., R_SE_self=5., R_SE_cross=2.)
        old_context = build_rf_response_contexts(a, b)["A"]
        a["R_u"], b["R_u"] = 2., 5.
        contexts = build_rf_response_contexts(a, b)
        def generator(context):
            eigenvalues, vectors, inverse = context["coherence_modes"]
            return (vectors * eigenvalues) @ inverse
        first = len(contexts["A"]["coherence_pairs"])
        second = len(contexts["B"]["coherence_pairs"])
        expected = -np.diag(np.r_[np.full(first, 2.), np.full(second, 5.)])
        np.testing.assert_allclose(generator(contexts["A"]) - generator(old_context), expected, atol=3e-11)

    def test_diagnostics_and_linewidth_normalization_include_uniform_rate(self):
        config, beams, common = system_inputs()
        config["R_u"] = 6.
        result = compute_alkali_system(config, None, beams, common)["A"]
        frame = result["df_pop"]
        populated = frame["population"] > 1e-15
        np.testing.assert_allclose(frame.loc[populated, "G_u"],
            6 * (1 - 1 / (len(frame) * frame.loc[populated, "population"])))
        finite = frame["Gamma_u"].notna()
        np.testing.assert_allclose(frame.loc[finite, "Gamma_u"], 6.)
        totals = add_total_relaxation_columns(frame)
        np.testing.assert_allclose(totals.loc[finite, "Gamma_total"],
            frame.loc[finite, ["Gamma_OP", "Gamma_ER", "Gamma_SE", "Gamma_u"]].sum(axis=1))
        old = frame.drop(columns=["G_u", "Gamma_u"])
        reference = largest_abs_Dm_relaxation_reference(frame, target_F="all")
        old_reference = largest_abs_Dm_relaxation_reference(old, target_F="all")
        self.assertAlmostEqual(reference["Gamma_m"] - old_reference["Gamma_m"], 6.)

    def test_physical_pump_readout_inherits_uniform_relaxation(self):
        config, beams, common = system_inputs(physical=True)
        old = compute_alkali_system(config, None, beams, common)["A"]
        config["R_u"] = 7.
        new = compute_alkali_system(config, None, beams, common)["A"]
        old_signal = phasor(old["probe_response"]["total"]["rotation"])
        new_signal = phasor(new["probe_response"]["total"]["rotation"])
        self.assertGreater(np.max(np.abs(new_signal - old_signal)), 1e-14)
        for manifold in ("Fa", "Fb"):
            self.assertTrue(np.isfinite(new["rf_response"][manifold]["amplitude"]).all())
        config["R_u"] = 0.
        zero = compute_alkali_system(config, None, beams, common)["A"]
        np.testing.assert_array_equal(zero["population"], old["population"])
        np.testing.assert_array_equal(phasor(zero["probe_response"]["total"]["rotation"]), old_signal)

    def test_dark_dual_species_retain_independent_uniform_rates(self):
        config, beams, common = system_inputs()
        config.update(density_cm3=0., R_ER=0., R_u=2.)
        beams[0].update(k_axis="z", pol="sigma+")
        other = deepcopy(config)
        other.update(label="B", atom_name="Rb85", R_u=5.)
        result = compute_population_evolution(config, other, beams, common, .05, 1., 41)
        dark = result["time_s"] >= .05
        elapsed = result["time_s"][dark] - .05
        for label, rate in (("A", 2.), ("B", 5.)):
            populations = result[label]["population"][dark]
            uniform = 1 / populations.shape[1]
            expected = uniform + (populations[0] - uniform)[None, :] * np.exp(-rate * elapsed[:, None])
            np.testing.assert_allclose(populations, expected, atol=3e-8)


class UniformRelaxationSettingsTests(unittest.TestCase):
    def test_new_defaults_and_independent_rate_round_trips(self):
        for label in ("A", "B"):
            self.assertEqual(conditions.DEFAULT_STARTUP_CONDITION[f"gamma_u_{label}"], 1.)
            self.assertEqual(DEFAULTS[f"gamma_u_{label}"], 1.)
        values = dict(conditions.DEFAULT_STARTUP_CONDITION, gamma_u_A=2.5, gamma_u_B=7.)
        state = {}
        with patch.object(conditions.st, "session_state", state):
            conditions.apply_loaded_condition_dict(conditions.build_condition_payload(values))
        self.assertEqual(state["gamma_u_A"], 2.5)
        self.assertEqual(state["gamma_u_B"], 7.)
        values = dict(DEFAULTS, gamma_u_A=0., gamma_u_B=4.)
        self.assertEqual(parse_settings_payload(build_settings_payload(values)), values)

    def test_current_files_require_finite_nonnegative_uniform_rates(self):
        for value in (-1., None, True, float("nan"), float("inf")):
            payload = conditions.build_condition_payload(conditions.DEFAULT_STARTUP_CONDITION)
            payload["conditions"]["gamma_u_A"] = value
            state = {"condition_name": "unchanged"}
            with patch.object(conditions.st, "session_state", state), self.assertRaises(ValueError):
                conditions.apply_loaded_condition_dict(payload)
            self.assertEqual(state, {"condition_name": "unchanged"})
        for label in ("A", "B"):
            payload = conditions.build_condition_payload(conditions.DEFAULT_STARTUP_CONDITION)
            payload["conditions"].pop(f"gamma_u_{label}")
            with patch.object(conditions.st, "session_state", {}), self.assertRaises(ValueError):
                conditions.apply_loaded_condition_dict(payload)
            payload = build_settings_payload(DEFAULTS)
            payload["settings"].pop(f"gamma_u_{label}")
            with self.assertRaises(ValueError):
                parse_settings_payload(payload)


if __name__ == "__main__":
    unittest.main()
