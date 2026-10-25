"""Stationary coherent pumping, spectral propagation and frequency exports."""

from copy import deepcopy
import unittest
from unittest.mock import patch

import numpy as np

from alkali_pumping_app.physics.coherent_evolution import density_derivative, prepare_density_model
from alkali_pumping_app.physics.constants import ATOMS
from alkali_pumping_app.physics.optical_readout import density_probe_signals, density_probe_spectrum
from alkali_pumping_app.physics.probe_scan import pump_on_steady_state
from alkali_pumping_app.physics.time_evolution import build_population_dynamics
from alkali_pumping_app.ui.time_evolution_data import probe_scan_export_dataframe
from tests.test_time_evolution import inputs


class ProbeScanTests(unittest.TestCase):
    def test_stationary_state_matches_linear_master_equation_without_spin_exchange(self):
        config, beams, common = inputs()
        config.update(density_cm3=0., R_ER=2., R_u=1.)
        common.update(static_field_axis="x", static_field_nT=.05)
        model = prepare_density_model(build_population_dynamics(config, None, beams, common)[0], common, "z")
        size, indices = model["size"], model["density_indices"]
        identity = np.eye(len(indices), dtype=complex)
        generator = np.column_stack([density_derivative(column, [model], True) for column in identity.T])
        generator[0] = np.eye(size).ravel()[indices]
        target = np.zeros(len(indices), complex)
        target[0] = 1
        exact = np.linalg.solve(generator, target)
        with patch("alkali_pumping_app.physics.system_model._evaluate_system_responses",
                   side_effect=AssertionError("A probe scan must not solve RF")):
            steady = pump_on_steady_state(config, None, beams, common)
        actual = steady["A"]["density_matrix"][0]
        np.testing.assert_allclose(actual.ravel()[indices], exact, atol=2e-8)
        self.assertAlmostEqual(np.trace(actual).real, 1.)
        self.assertGreaterEqual(np.linalg.eigvalsh(actual).min(), -1e-10)
        self.assertLess(steady["stationary_residual_s-1"], 1e-7)

    def test_zero_pump_is_unpolarized_including_dual_species(self):
        config, _, common = inputs()
        config.update(density_cm3=1e10, R_u=1.)
        other = deepcopy(config)
        other.update(label="B", atom_name="Rb85")
        common.update(static_field_axis="y", static_field_nT=.1)
        steady = pump_on_steady_state(config, other, [], common)
        for label in ("A", "B"):
            rho = steady[label]["density_matrix"][0]
            np.testing.assert_allclose(rho, np.eye(len(rho)) / len(rho), atol=1e-12)
        self.assertEqual(steady["convergence_time_s"], 0.)

    def test_spectrum_matches_single_frequency_readout_for_coherent_states(self):
        config, _, common = inputs()
        model = build_population_dynamics(config, None, [], common)[0]
        states = model["ground_states"]
        vector = np.zeros(len(states), complex)
        vector[-1], vector[-2] = 1 / np.sqrt(2), 1j / np.sqrt(2)
        rho = np.outer(vector, vector.conj())
        original = rho.copy()
        frequencies = np.array([-14000., -2000., 1000., 14000.])
        for line in ("D1", "D2"):
            arguments = dict(atom=ATOMS[config["atom_name"]], ground_states=states, q_axis="z", line=line,
                k_axis="z", azimuth_deg=17., ellipticity_deg=8., path_length_cm=.1, density_cm3=1e10,
                n2_pressure_torr=common["n2_pressure_torr"], temperature_C=common["temperature_C"],
                n2_width_MHz_per_torr=17.8, n2_shift_MHz_per_torr=-8.25)
            spectrum = density_probe_spectrum(density_matrix=rho, detunings_MHz=frequencies, **arguments)
            for index, frequency in enumerate(frequencies):
                single = density_probe_signals(density_matrix=rho[None], detuning_MHz=frequency, **arguments)
                for signal in ("transmission", "absorption", "s1", "s2", "s3", "ellipticity_rad"):
                    np.testing.assert_allclose(spectrum[signal][index], single[signal][0], atol=1e-12)
                np.testing.assert_allclose(np.exp(2j * spectrum["rotation_rad"][index]),
                                           np.exp(2j * single["rotation_rad"][0]), atol=1e-12)
            self.assertEqual(len(spectrum["transition_centers_MHz"]), 4 if line == "D1" else 6)
        np.testing.assert_array_equal(rho, original)

    def test_scan_export_has_frequency_and_plotted_units_without_event_times(self):
        config, _, common = inputs()
        config.update(R_u=1.)
        steady = pump_on_steady_state(config, None, [], common)
        species = steady["A"]
        probe = dict(config["probe"], line="D1")
        arguments = {name: probe[name] for name in ("line", "k_axis", "azimuth_deg", "ellipticity_deg", "path_length_cm")}
        spectrum = density_probe_spectrum(atom=ATOMS[config["atom_name"]], ground_states=species["ground_states"],
            density_matrix=species["density_matrix"][0], q_axis="z", detunings_MHz=np.array([-100., 100.]),
            density_cm3=config["density_cm3"], n2_pressure_torr=common["n2_pressure_torr"],
            temperature_C=common["temperature_C"], n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25, **arguments)
        table = probe_scan_export_dataframe(spectrum, species=species, label="A", probe=probe,
                                            config=config, common=common, beams=[], steady=steady)
        np.testing.assert_array_equal(table["detuning_MHz"], [-100., 100.])
        np.testing.assert_array_equal(table["rotation_mrad"], 1000 * spectrum["rotation_rad"])
        np.testing.assert_array_equal(table["ellipticity_mrad"], 1000 * spectrum["ellipticity_rad"])
        np.testing.assert_array_equal(table["absorption_1_minus_Iout_over_Iin"], spectrum["absorption"])
        self.assertTrue(table["pumps_on"].all())
        self.assertNotIn("time_s", table)
        self.assertNotIn("pump_off_time_s", table)


if __name__ == "__main__":
    unittest.main()
