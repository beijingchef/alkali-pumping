"""Physical limits and integration checks for the common-field Fa/Fb model."""

from copy import deepcopy
import unittest

import numpy as np
import pandas as pd

from alkali_pumping_app.physics import (
    ATOMS, build_ground_states, ground_zeeman_shifts_hz,
    ground_hyperfine_lande_g, line_center_frequency_MHz,
    compute_alkali_system, compute_alkali_responses,
    atomic_rf_readout_operators, build_rf_response_contexts,
    magnetic_rf_operator, weak_rf_matrix_readouts,
    weak_rf_observable_susceptibility, weak_liouvillian_source_matrix_readouts,
)
from alkali_pumping_app.physics.rf_response import (
    _spin_operator, _single_species_coherence_generator, _drive_vector,
)


def phasor(response):
    if isinstance(response, dict):
        return response["in_phase"] + 1j * response["quadrature"]
    return response[1] + 1j * response[2]


def synthetic_result(atom_name="Rb87", frequencies=(-51.0, -17.0, 0.0, 17.0, 51.0)):
    atom = ATOMS[atom_name]
    states = build_ground_states(atom)
    p = np.arange(1.0, len(states) + 1.0) ** 2
    p /= p.sum()
    zeeman, _ = ground_zeeman_shifts_hz(atom_name, atom, states, 50.0)
    size = len(states)
    return {
        "atom_name": atom_name, "atom": atom, "ground_states": states,
        "population": p, "nu_LS": np.zeros(size),
        "df_pop": pd.DataFrame({"nu_B": zeeman, "G_OP": np.full(size, 4.0)}),
        "R_ER": 0.0, "R_SE_self": 0.0, "R_SE_cross": 0.0,
        "rf_upper_F": max(state["F"] for state in states),
        "q_axis": "z", "rf_axis": "x", "rf_observable": "Fx",
        "rf_frequencies_hz": np.asarray(frequencies), "light_shift_available": True,
        "J_population": -4.0 * (np.eye(size) - np.ones((size, size)) / size),
    }


def system_inputs(physical=False):
    probe = {
        "source": "PumpA1" if physical else "Custom", "mode": "physical" if physical else "weak",
        "pump_name": "PumpA1", "pump_intensity_uW_cm2": 20.0,
        "line": "D1", "detuning_MHz": 500.0, "k_axis": "x",
        "azimuth_deg": 90.0, "ellipticity_deg": 0.0, "path_length_cm": 2.0,
        "include_scalar": False, "include_orientation": True, "include_alignment": True,
    }
    species = {
        "label": "A", "atom_name": "Rb87", "density_cm3": 1e10, "R_ER": 10.0,
        "n2_coeffs": {"D1": {"width": 17.8, "shift": -8.25},
                      "D2": {"width": 18.1, "shift": -5.9}},
        "q_axis": "z", "rf_axis": "x", "rf_observable": "Q_ij",
        "rf_frequencies_hz": np.array([0.0, 10.0, 30.0, 50.0]), "probe": probe,
    }
    beam = {
        "name": "PumpA1", "target_label": "A", "target_atom": "Rb87",
        "absolute_frequency_MHz": line_center_frequency_MHz(ATOMS["Rb87"], "D1") + 500.0,
        "intensity": 20.0, "k_axis": "x", "pol": "linear z",
        "selected_transition": None, "transition_label": "D1 test pump",
    }
    common = {"temperature_C": 23.0, "n2_pressure_torr": 0.0,
              "static_field_axis": "z", "static_field_nT": 3.0}
    return species, [beam], common


class TwoManifoldPhysicsTests(unittest.TestCase):
    def test_signed_g_factors_include_negative_nuclear_zeeman_energy(self):
        ga = ground_hyperfine_lande_g("Rb87", ATOMS["Rb87"], 2.0)
        gb = ground_hyperfine_lande_g("Rb87", ATOMS["Rb87"], 1.0)
        # Steck values, allowing for the tabulated-moment approximation.
        self.assertAlmostEqual(ga, 0.4998364, delta=5e-6)
        self.assertAlmostEqual(gb, -0.5018267, delta=5e-6)
        for atom_name in ATOMS:
            result = synthetic_result(atom_name)
            drive, info = magnetic_rf_operator(result)
            np.testing.assert_allclose(drive, drive.conj().T)
            self.assertEqual(info["drive_ratios"]["Fa"], 1.0)
            self.assertLess(info["drive_ratios"]["Fb"], 0.0)
            for name, F in info["manifolds"].items():
                indices = [i for i, state in enumerate(result["ground_states"]) if state["F"] == F]
                spin = _spin_operator(result["ground_states"], "z", "x", F)
                np.testing.assert_allclose(drive[np.ix_(indices, indices)],
                    info["drive_ratios"][name] * spin[np.ix_(indices, indices)])

    def test_independent_transition_limit_matches_analytic_both_frequency_terms(self):
        for atom_name in ATOMS:
            result = synthetic_result(atom_name)
            states = result["ground_states"]
            operators = atomic_rf_readout_operators(result)
            responses = weak_rf_matrix_readouts(result, operators)
            _, info = magnetic_rf_operator(result)
            transitions = np.full(len(states), np.nan)
            by_state = {(state["F"], state["m"]): i for i, state in enumerate(states)}
            zeeman = result["df_pop"]["nu_B"].to_numpy()
            for i, state in enumerate(states):
                j = by_state.get((state["F"], state["m"] - 1.0))
                if j is not None:
                    transitions[i] = zeeman[i] - zeeman[j]
            for name, F in info["manifolds"].items():
                expected = weak_rf_observable_susceptibility(
                    result["rf_frequencies_hz"], states, result["population"], transitions,
                    np.full(len(states), 4.0), np.zeros(len(states)), np.zeros(len(states)),
                    "z", "x", "Fx", target_F=F,
                )
                np.testing.assert_allclose(phasor(responses[name]),
                    info["drive_ratios"][name] * phasor(expected), rtol=2e-10, atol=1e-13)
                values = phasor(responses[name])
                np.testing.assert_allclose(values, values[::-1].conj(), atol=1e-13)
                self.assertGreater(np.max(np.abs(values)), 1e-8)

    def test_lower_only_population_gradient_drives_lower_manifold(self):
        result = synthetic_result()
        for i, state in enumerate(result["ground_states"]):
            if state["F"] == 2.0:
                result["population"][i] = 0.1
        result["population"] /= result["population"].sum()
        responses = weak_rf_matrix_readouts(result, atomic_rf_readout_operators(result))
        np.testing.assert_allclose(responses["Fa"][0], 0.0, atol=1e-14)
        self.assertGreater(np.max(responses["Fb"][0]), 1e-6)
        # Responses are per total atom; halving a manifold's source halves its response.
        result["population"] *= 0.5
        half = weak_rf_matrix_readouts(result, atomic_rf_readout_operators(result))
        np.testing.assert_allclose(phasor(half["Fb"]), 0.5 * phasor(responses["Fb"]), atol=1e-14)

    def test_collision_transfer_generates_upper_response_from_lower_source(self):
        result = synthetic_result()
        for i, state in enumerate(result["ground_states"]):
            if state["F"] == 2.0:
                result["population"][i] = 0.1
        result["population"] /= result["population"].sum()
        result["R_ER"] = 5.0
        responses = weak_rf_matrix_readouts(result, atomic_rf_readout_operators(result))
        self.assertGreater(np.max(responses["Fa"][0]), 1e-7)

    def test_shared_resolvent_matches_direct_linear_solve_and_secular_basis(self):
        result = synthetic_result()
        result.update(R_ER=3.0, R_SE_self=6.0)
        pairs, generator = _single_species_coherence_generator(result)
        self.assertTrue(all(result["ground_states"][a]["F"] == result["ground_states"][b]["F"]
                            for a, b in pairs))
        source, _ = _drive_vector(result, pairs, len(pairs), 0)
        operators = atomic_rf_readout_operators(result)
        responses = weak_rf_matrix_readouts(result, operators)
        for name, operator in operators.items():
            row = np.array([operator[b, a] for a, b in pairs])
            reference = np.array([2 * row @ np.linalg.solve(
                -2j * np.pi * frequency * np.eye(len(pairs)) - generator, source)
                for frequency in result["rf_frequencies_hz"]])
            np.testing.assert_allclose(phasor(responses[name]), reference, rtol=2e-9, atol=1e-12)

    def test_zero_cross_exchange_matches_single_species_for_rf_and_optical_sources(self):
        a, b = synthetic_result(), synthetic_result("Rb85", (10.0, 20.0))
        a.update(R_ER=3.0, R_SE_self=5.0)
        b.update(R_ER=4.0, R_SE_self=7.0)
        full_population = np.block([
            [a["J_population"], np.zeros((8, 12))],
            [np.zeros((12, 8)), b["J_population"]],
        ])
        contexts = build_rf_response_contexts(a, b, full_population)
        for label, result in (("A", a), ("B", b)):
            operators = atomic_rf_readout_operators(result)
            local = weak_rf_matrix_readouts(result, operators)
            coupled = weak_rf_matrix_readouts(result, operators, context=contexts[label])
            for name in operators:
                np.testing.assert_allclose(phasor(coupled[name]), phasor(local[name]), atol=2e-12)
            drive, _ = magnetic_rf_operator(result)
            rho = np.diag(result["population"])
            source = -1j * (drive @ rho - rho @ drive)
            np.testing.assert_allclose(source, source.conj().T)
            self.assertAlmostEqual(np.trace(source), 0.0)
            local = weak_liouvillian_source_matrix_readouts(result, source, operators)
            coupled = weak_liouvillian_source_matrix_readouts(result, source, operators, context=contexts[label])
            for name in operators:
                np.testing.assert_allclose(phasor(coupled[name]), phasor(local[name]), atol=2e-12)

    def test_complex_readouts_add_and_can_cancel(self):
        result = synthetic_result()
        operators = atomic_rf_readout_operators(result)
        operators["sum"] = operators["Fa"] + operators["Fb"]
        operators["negative"] = -operators["Fb"]
        operators["cancel"] = operators["Fb"] + operators["negative"]
        responses = weak_rf_matrix_readouts(result, operators)
        np.testing.assert_allclose(phasor(responses["sum"]),
            phasor(responses["Fa"]) + phasor(responses["Fb"]), atol=1e-14)
        np.testing.assert_allclose(responses["cancel"][0], 0.0)
        self.assertTrue(np.any(responses["sum"][0] < responses["Fa"][0] + responses["Fb"][0] - 1e-10))


class TwoManifoldSystemTests(unittest.TestCase):
    def check_manifold_sum(self, result):
        for rank, signals in result["probe_response"].items():
            for signal, total in signals.items():
                parts = result["probe_manifold_response"]
                expected = phasor(parts["Fa"][rank][signal]) + phasor(parts["Fb"][rank][signal])
                np.testing.assert_allclose(phasor(total), expected, rtol=1e-8, atol=1e-15)
                np.testing.assert_allclose(total["amplitude"], np.abs(expected), rtol=1e-8, atol=1e-15)

    def test_weak_probe_and_cached_response_use_identical_two_manifold_model(self):
        species, beams, common = system_inputs()
        system = compute_alkali_system(species, None, beams, common)
        result = system["A"]
        self.check_manifold_sum(result)
        for response in result["rf_response"].values():
            self.assertTrue(np.isfinite(response["amplitude"]).all())
            self.assertGreater(np.max(response["amplitude"]), 0.0)
        updated = compute_alkali_responses(system, species, None, common)["A"]
        for manifold in ("Fa", "Fb"):
            np.testing.assert_allclose(phasor(updated["rf_response"][manifold]),
                                       phasor(result["rf_response"][manifold]), atol=1e-14)
        other_observable = deepcopy(species)
        other_observable["rf_observable"] = "Fy"
        alternate = compute_alkali_responses(system, other_observable, None, common)["A"]
        np.testing.assert_allclose(phasor(alternate["probe_response"]["total"]["rotation"]),
                                   phasor(result["probe_response"]["total"]["rotation"]), atol=1e-15)
        changed = deepcopy(species)
        changed["probe"]["detuning_MHz"] += 1000.0
        changed["rf_observable"] = "Fy"
        new = compute_alkali_responses(system, changed, None, common)["A"]
        np.testing.assert_array_equal(new["population"], result["population"])
        np.testing.assert_array_equal(new["nu_LS"], result["nu_LS"])
        np.testing.assert_array_equal(new["L_optical"], result["L_optical"])
        self.check_manifold_sum(new)

    def test_physical_manifold_signals_share_the_same_feedback_closure(self):
        species, beams, common = system_inputs(physical=True)
        result = compute_alkali_system(species, None, beams, common)["A"]
        self.check_manifold_sum(result)
        for manifold in ("Fa", "Fb"):
            self.assertTrue(result["probe_manifold_response"][manifold]["total"]["rotation"]["info"]["shared_manifold_feedback"])
        species["probe"]["path_length_cm"] = 0.0
        zero = compute_alkali_system(species, None, beams, common)["A"]
        self.check_manifold_sum(zero)
        np.testing.assert_allclose(zero["probe_response"]["total"]["rotation"]["amplitude"], 0.0)

    def test_dual_species_independent_grids_and_shared_manifold_readout(self):
        species, beams, common = system_inputs(physical=True)
        other = deepcopy(species)
        other.update(label="B", atom_name="Rb85", density_cm3=2e10,
                     rf_frequencies_hz=np.array([5.0, 15.0]))
        other["probe"]["mode"] = "weak"
        system = compute_alkali_system(species, other, beams, common)
        for label, count in (("A", 4), ("B", 2)):
            result = system[label]
            self.check_manifold_sum(result)
            for response in result["rf_response"].values():
                self.assertEqual(response["amplitude"].shape, (count,))
                self.assertTrue(np.isfinite(response["amplitude"]).all())


if __name__ == "__main__":
    unittest.main()
