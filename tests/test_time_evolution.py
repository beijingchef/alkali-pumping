"""Physical and numerical checks for pump transients independent of RF."""

from copy import deepcopy
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import expm

from alkali_pumping_app.physics.constants import ATOMS
from alkali_pumping_app.physics.multi_species import build_cross_spin_exchange_matrix
from alkali_pumping_app.physics.optical_readout import (
    population_probe_signals, weak_probe_readout_operators,
)
from alkali_pumping_app.physics.spin_exchange import build_spin_exchange_matrix
from alkali_pumping_app.physics.time_evolution import (
    build_population_dynamics, compute_population_evolution, population_derivative, population_moments,
)
from alkali_pumping_app.ui.time_evolution_data import (
    evolution_export_dataframe, plot_export_dataframe, section_export_dataframe, series_frame, zeeman_plot_frame,
)
from tests.test_two_manifold_rf import system_inputs


def inputs():
    config, beams, common = system_inputs()
    beams[0].update(k_axis="z", pol="sigma+")
    return config, beams, common


def probe_signals(species, config, common, **overrides):
    probe = dict(config["probe"])
    probe = {key: probe[key] for key in ("line", "detuning_MHz", "k_axis", "azimuth_deg",
                                        "ellipticity_deg", "path_length_cm")}
    probe.update(overrides)
    line = probe["line"]
    return population_probe_signals(
        atom=ATOMS[config["atom_name"]], ground_states=species["ground_states"],
        population=species["population"], q_axis=config["q_axis"], **probe,
        density_cm3=config["density_cm3"], temperature_C=common["temperature_C"],
        n2_pressure_torr=common["n2_pressure_torr"],
        n2_width_MHz_per_torr=config["n2_coeffs"][line]["width"],
        n2_shift_MHz_per_torr=config["n2_coeffs"][line]["shift"],
    )


class PopulationEvolutionTests(unittest.TestCase):
    def test_hyperfine_polarization_matches_spin_dot_product_for_all_isotopes(self):
        from alkali_pumping_app.physics.angular_momentum import build_ground_states, coupled_basis_amplitudes
        from alkali_pumping_app.physics.coherent_evolution import spin_matrices
        from alkali_pumping_app.physics.spin_exchange import hyperfine_uncoupled_probabilities
        rng = np.random.default_rng(510)
        for atom_name, atom in ATOMS.items():
            states = build_ground_states(atom)
            n = len(states)
            probabilities, _, _ = hyperfine_uncoupled_probabilities(atom, states)
            model = {"atom": atom, "ground_states": states, "electron_by_state": probabilities.sum(axis=1)}
            I = atom["I"]
            nuclear = spin_matrices([{"F": I, "m": m} for m in np.arange(-I, I+1)])
            electron = spin_matrices([{"F": .5, "m": m} for m in (-.5, .5)])
            uncoupled_operator = sum(np.kron(i, s) for i, s in zip(nuclear, electron))
            amplitudes = coupled_basis_amplitudes(atom, states).reshape(n, n)
            operator = amplitudes @ uncoupled_operator @ amplitudes.T
            mask = np.array([[a["F"] == b["F"] for b in states] for a in states])
            random = (rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))) * mask
            rho = random @ random.conj().T
            rho /= np.trace(rho)
            actual = population_moments(model, np.diag(rho).real[None])["hyperfine_polarization"][0]
            with self.subTest(atom=atom_name):
                self.assertAlmostEqual(actual, np.trace(rho @ operator).real, places=13)
                initial = population_moments(model, np.full((1, n), 1/n))["hyperfine_polarization"][0]
                self.assertAlmostEqual(initial, 0., places=13)
                for F, expected in ((I+.5, I/2), (I-.5, -(I+1)/2)):
                    p = np.array([s["F"] == F for s in states], dtype=float)
                    p /= p.sum()
                    result = population_moments(model, p[None])["hyperfine_polarization"][0]
                    self.assertAlmostEqual(result, expected, places=13)

    def test_nonlinear_derivative_matches_existing_self_and_cross_maps(self):
        config, beams, common = inputs()
        other = deepcopy(config)
        other.update(label="B", atom_name="Rb85", density_cm3=3e10)
        models = build_population_dynamics(config, other, beams, common)
        populations = [np.arange(1., model["size"] + 1)**2 for model in models]
        populations = [p / p.sum() for p in populations]
        expected = []
        electrons = []
        for model, p in zip(models, populations):
            matrix, electron = build_spin_exchange_matrix(model["atom"], model["ground_states"], p)
            expected.append(model["L_linear"] @ p + model["R_SE_self"] * (matrix @ p - p))
            electrons.append(electron)
        for index, (model, p) in enumerate(zip(models, populations)):
            cross = build_cross_spin_exchange_matrix(model["atom"], model["ground_states"], electrons[1-index])
            expected[index] += model["R_SE_cross"] * (cross @ p - p)
        actual = population_derivative(np.concatenate(populations), models, True)
        np.testing.assert_allclose(actual, np.concatenate(expected), atol=1e-12)

    def test_linear_limit_matches_matrix_exponentials_through_switch(self):
        config, beams, common = inputs()
        config["density_cm3"] = 0.0  # no self or unlike-alkali exchange
        evolution = compute_population_evolution(config, None, beams, common, .03, .1, 51)
        model = build_population_dynamics(config, None, beams, common)[0]
        p0 = np.full(model["size"], 1 / model["size"])
        switch = expm(model["L_linear"] * .03) @ p0
        expected = np.array([
            expm(model["L_linear"] * t) @ p0 if t <= .03
            else expm(model["L_dark"] * (t - .03)) @ switch
            for t in evolution["time_s"]
        ])
        np.testing.assert_allclose(evolution["A"]["population"], expected, atol=2e-8)
        self.assertEqual(np.count_nonzero(evolution["time_s"] == .03), 1)

    def test_dual_species_normalization_initial_state_and_exchange_transfer(self):
        config, beams, common = inputs()
        other = deepcopy(config)
        other.update(label="B", atom_name="Rb85", density_cm3=2e10)
        result = compute_population_evolution(config, other, beams, common, .05, .2, 61)
        for label in ("A", "B"):
            population = result[label]["population"]
            np.testing.assert_allclose(population.sum(axis=1), 1, atol=1e-12)
            self.assertGreaterEqual(population.min(), 0)
            np.testing.assert_allclose(population[0], 1 / population.shape[1])
            moments = result[label]["moments"]
            np.testing.assert_allclose(moments["population_Fa"] + moments["population_Fb"], 1)
        self.assertGreater(np.max(result["A"]["moments"]["electron_polarization"]), .001)
        self.assertFalse(result["rf_enabled"])

    def test_no_pumping_and_endpoint_switches(self):
        config, beams, common = inputs()
        for switch in (0., .1):
            result = compute_population_evolution(config, None, beams, common, switch, .1, 21)
            self.assertEqual(result["time_s"][0], 0)
            self.assertEqual(result["time_s"][-1], .1)
            if switch == 0:
                np.testing.assert_allclose(result["A"]["population"], 1 / 8, atol=1e-10)
        for beam in beams:
            beam["intensity"] = 0.
        result = compute_population_evolution(config, None, beams, common, .05, .1, 21)
        np.testing.assert_allclose(result["A"]["population"], 1 / 8, atol=1e-10)

    def test_dynamics_never_calls_rf_solver_and_ignores_probe_and_rf_settings(self):
        config, beams, common = inputs()
        with patch("alkali_pumping_app.physics.system_model._evaluate_system_responses",
                   side_effect=AssertionError("RF solver must not run")):
            first = compute_population_evolution(config, None, beams, common, .02, .05, 21)
            config.update(rf_axis="z", rf_frequencies_hz=np.array([1e6]))
            config["probe"].update(detuning_MHz=-4000., path_length_cm=100.)
            second = compute_population_evolution(config, None, beams, common, .02, .05, 21)
        np.testing.assert_array_equal(first["A"]["population"], second["A"]["population"])

    def test_dark_no_relaxation_preserves_total_angular_momentum(self):
        config, beams, common = inputs()
        config["R_ER"] = 0.
        result = compute_population_evolution(config, None, beams, common, .04, .4, 81)
        moments = result["A"]["moments"]
        total = moments["orientation_Fa"] + moments["orientation_Fb"]
        dark = result["time_s"] >= .04
        np.testing.assert_allclose(total[dark], total[dark][0], atol=1e-8)

    def test_inactive_B_and_zero_cross_exchange_recover_single_species(self):
        config, beams, common = inputs()
        single = compute_population_evolution(config, None, beams, common, .02, .05, 21)
        same = deepcopy(config)
        same["label"] = "B"
        inactive = compute_population_evolution(config, same, beams, common, .02, .05, 21)
        self.assertIsNone(inactive["B"])
        np.testing.assert_array_equal(single["A"]["population"], inactive["A"]["population"])
        other = deepcopy(config)
        other.update(label="B", atom_name="Rb85", density_cm3=0.)
        dual = compute_population_evolution(config, other, beams, common, .02, .05, 21)
        np.testing.assert_allclose(single["A"]["population"], dual["A"]["population"], atol=2e-8)

    def test_invalid_times_rejected(self):
        config, beams, common = inputs()
        for t1, t2 in ((-1, 1), (2, 1), (0, 0), (0, np.nan)):
            with self.subTest(t1=t1, t2=t2), self.assertRaises(ValueError):
                compute_population_evolution(config, None, beams, common, t1, t2)


class PopulationProbeTests(unittest.TestCase):
    def setUp(self):
        self.config, self.beams, self.common = inputs()
        self.evolution = compute_population_evolution(self.config, None, self.beams, self.common, .02, .05, 21)

    def test_zero_path_and_unpolarized_isotropic_state(self):
        result = probe_signals(self.evolution["A"], self.config, self.common, path_length_cm=0., ellipticity_deg=12.)
        np.testing.assert_allclose(result["rotation_rad"], 0, atol=1e-14)
        np.testing.assert_allclose(result["ellipticity_rad"], np.deg2rad(12.), atol=1e-14)
        np.testing.assert_allclose(result["transmission"], 1, atol=1e-14)
        np.testing.assert_allclose(result["absorption"], 0, atol=1e-14)
        result = probe_signals(self.evolution["A"], self.config, self.common, k_axis="z", azimuth_deg=0.)
        self.assertAlmostEqual(result["rotation_rad"][0], 0, places=12)
        self.assertAlmostEqual(result["ellipticity_rad"][0], 0, places=12)
        self.assertTrue(np.all(result["transmission"] > 0))
        self.assertTrue(np.all(result["transmission"] <= 1 + 1e-10))
        np.testing.assert_allclose(result["absorption"] + result["transmission"], 1, atol=1e-14)
        self.assertTrue(np.all(result["absorption"] >= -1e-10))
        self.assertTrue(np.all(result["absorption"] <= 1))
        np.testing.assert_allclose(result["s1"]**2 + result["s2"]**2 + result["s3"]**2, 1, atol=1e-12)

    def test_small_optical_depth_recovers_existing_linear_readout(self):
        config = deepcopy(self.config)
        config["density_cm3"] = 1e5
        probe = dict(config["probe"])
        probe.update(k_axis="z", azimuth_deg=0.)
        config["probe"] = probe
        species = self.evolution["A"]
        actual = probe_signals(species, config, self.common)
        operators, _ = weak_probe_readout_operators(
            atom=ATOMS[config["atom_name"]], ground_states=species["ground_states"], q_axis=config["q_axis"],
            **{key: probe[key] for key in ("line", "detuning_MHz", "k_axis", "azimuth_deg", "ellipticity_deg", "path_length_cm")},
            density_cm3=config["density_cm3"], temperature_C=self.common["temperature_C"],
            n2_pressure_torr=self.common["n2_pressure_torr"], n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25, include_scalar=True,
        )
        for signal in ("rotation", "ellipticity", "transmission"):
            expected = species["population"] @ np.diag(operators["total"][signal]).real
            key = signal if signal == "transmission" else signal + "_rad"
            if signal == "transmission":
                expected += 1
            np.testing.assert_allclose(actual[key], expected, atol=1e-12)
        expected_absorption = -(species["population"] @ np.diag(operators["total"]["transmission"]).real)
        np.testing.assert_allclose(actual["absorption"], expected_absorption, atol=1e-12)

    def test_probe_readout_is_pure_and_circular_input_has_no_rotation_angle(self):
        before = self.evolution["A"]["population"].copy()
        result = probe_signals(self.evolution["A"], self.config, self.common, ellipticity_deg=45.)
        self.assertTrue(np.isnan(result["rotation_rad"]).all())
        self.assertTrue(np.isfinite(result["transmission"]).all())
        np.testing.assert_array_equal(before, self.evolution["A"]["population"])

    def test_plots_and_csv_preserve_units_without_zeeman_totals(self):
        signals = probe_signals(self.evolution["A"], self.config, self.common)
        table = evolution_export_dataframe(self.evolution, "A", signals, self.config,
                                            self.config["probe"], self.common, self.beams)
        np.testing.assert_array_equal(table["rotation_mrad"], 1000 * signals["rotation_rad"])
        np.testing.assert_array_equal(table["absorption_1_minus_Iout_over_Iin"], signals["absorption"])
        np.testing.assert_allclose(table["absorption_1_minus_Iout_over_Iin"] + table["transmission_Iout_over_Iin"], 1)
        np.testing.assert_array_equal(table["population_Fa"], self.evolution["A"]["moments"]["population_Fa"])
        np.testing.assert_array_equal(table["hyperfine_polarization"], self.evolution["A"]["moments"]["hyperfine_polarization"])
        self.assertFalse(table["RF_enabled"].any())
        for manifold, count in (("Fa", 5), ("Fb", 3)):
            plot = zeeman_plot_frame(self.evolution["time_s"], self.evolution["A"], manifold)
            self.assertEqual(plot["Series"].nunique(), count)
            self.assertTrue(all(name.startswith("m =") for name in plot["Series"].unique()))

    def test_graph_csv_contains_only_shown_curves_with_units_and_event_metadata(self):
        species = self.evolution["A"]
        times = self.evolution["time_s"]
        displayed = 1000 * probe_signals(species, self.config, self.common)["rotation_rad"]
        frames = (
            (series_frame(times, {"Optical rotation": displayed}, "Rotation (mrad)"),
             "Rotation (mrad)", "Optical rotation", "mrad", {"Optical rotation": displayed}),
            (series_frame(times, {"Fa": species["moments"]["orientation_Fa"]}, "Moment"),
             "Moment", "Orientation", "hbar", {"Fa": species["moments"]["orientation_Fa"]}),
        )
        for frame, y_label, title, unit, expected in frames:
            table = plot_export_dataframe(frame, y_label, species=species, label="A", title=title,
                                          value_unit=unit, t1_s=.02, t2_s=.05)
            self.assertEqual(list(table.columns[:1 + len(expected)]), ["time_s", *expected])
            np.testing.assert_array_equal(table["time_s"], times)
            for name, values in expected.items():
                np.testing.assert_array_equal(table[name], values)
            self.assertNotIn("Fb", table)
            np.testing.assert_array_equal(table["pumps_on"], times < .02)
            self.assertEqual(set(table["pump_off_time_s"]), {.02})
            self.assertEqual(set(table["total_time_s"]), {.05})
            self.assertEqual(set(table["value_units"]), {unit})
            self.assertEqual(set(table["atom"]), {species["atom_name"]})
        frame = zeeman_plot_frame(times, species, "Fa")
        table = plot_export_dataframe(frame, "Population", species=species, label="A", title="Fa",
                                      value_unit="1", t1_s=.02, t2_s=.05)
        self.assertEqual(len([name for name in table if name.startswith("m =")]), 5)
        self.assertNotIn("Fa", table)
        self.assertNotIn("Fb", table)

    def test_graph_csv_keeps_undefined_signal_samples(self):
        times = self.evolution["time_s"]
        frame = series_frame(times, {"Optical rotation": np.full(len(times), np.nan)}, "Rotation (mrad)")
        table = plot_export_dataframe(frame, "Rotation (mrad)", species=self.evolution["A"], label="A",
                                      title="Optical rotation", value_unit="mrad", t1_s=.02, t2_s=.05)
        self.assertEqual(len(table), len(times))
        self.assertTrue(table["Optical rotation"].isna().all())

    def test_section_csv_disambiguates_curves_and_keeps_undefined_samples(self):
        species, times = self.evolution["A"], self.evolution["time_s"]
        plots = [dict(frame=zeeman_plot_frame(times, species, manifold), y_label="Population",
                      title=title, value_unit="1")
                 for manifold, title in (("Fa", "F = 2"), ("Fb", "F = 1"))]
        plots.append(dict(frame=series_frame(times, {"Optical rotation": np.full(len(times), np.nan)}, "Angle"),
                          y_label="Angle", title="Optical rotation", value_unit="mrad"))
        table = section_export_dataframe(plots, species=species, label="A", section="Example",
                                         t1_s=.02, t2_s=.05)
        self.assertEqual(len(table), len(times))
        self.assertEqual(len([name for name in table if " | m = " in name]), 8)
        for F in (1, 2):
            index = next(i for i, state in enumerate(species["ground_states"]) if state["F"] == F and state["m"] == 0)
            np.testing.assert_array_equal(table[f"F = {F} | m = 0 [1]"], species["population"][:, index])
        self.assertTrue(table["Optical rotation | Optical rotation [mrad]"].isna().all())
        np.testing.assert_array_equal(table["time_s"], times)
        self.assertTrue((table["section"] == "Example").all())


if __name__ == "__main__":
    unittest.main()
