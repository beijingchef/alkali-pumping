import unittest

import numpy as np

from alkali_pumping_app.physics import ATOMS, build_ground_states, compute_alkali_system
from alkali_pumping_app.physics.nonlinear_readout import (
    _SIGNALS,
    propagate_stokes_feedback,
)
from alkali_pumping_app.physics.rf_response import (
    _spin_operator,
    weak_drive_matrix_readouts,
    weak_rf_matrix_susceptibility,
)
from alkali_pumping_app.physics.spectroscopy import line_center_frequency_MHz


class StokesFeedbackPropagationTests(unittest.TestCase):
    def test_zero_feedback_reproduces_weak_full_cell_response(self):
        base = {
            signal: np.array([index + 1j * (index + 0.5)])
            for index, signal in enumerate(_SIGNALS)
        }
        feedback = {
            signal: {
                stokes: np.zeros(1, dtype=complex)
                for stokes in ("s1", "s2", "s3")
            }
            for signal in _SIGNALS
        }

        result, radius = propagate_stokes_feedback(base, feedback)

        for signal in _SIGNALS:
            np.testing.assert_allclose(result[signal], base[signal], atol=1e-14)
        np.testing.assert_allclose(radius, 0.0, atol=1e-14)

    def test_single_stokes_loop_matches_analytic_exponential(self):
        loop_gain = 0.3
        base = {
            signal: np.zeros(1, dtype=complex)
            for signal in _SIGNALS
        }
        base["s2"][:] = 2.0 - 0.5j
        feedback = {
            signal: {
                stokes: np.zeros(1, dtype=complex)
                for stokes in ("s1", "s2", "s3")
            }
            for signal in _SIGNALS
        }
        feedback["s2"]["s2"][:] = loop_gain

        result, radius = propagate_stokes_feedback(base, feedback)

        expected = base["s2"] * np.expm1(loop_gain) / loop_gain
        np.testing.assert_allclose(result["s2"], expected, rtol=1e-13)
        np.testing.assert_allclose(radius, loop_gain, rtol=1e-13)


class GeneralizedDriveTests(unittest.TestCase):
    def test_spin_drive_matches_existing_weak_rf_solver(self):
        states = build_ground_states(ATOMS["Rb87"])
        size = len(states)
        upper_F = max(float(state["F"]) for state in states)
        frequencies = np.array([7.0, 10.0, 13.0])
        population = np.linspace(1.0, 2.0, size)
        population /= population.sum()
        transition_hz = np.full(size, 10.0)
        rates = np.full(size, 3.0)
        operator = _spin_operator(states, "z", "x", upper_F)

        expected = weak_rf_matrix_susceptibility(
            frequencies,
            states,
            population,
            transition_hz,
            rates,
            rates,
            rates,
            "z",
            "x",
            operator,
            upper_F,
        )
        actual = weak_drive_matrix_readouts(
            frequencies,
            states,
            population,
            transition_hz,
            rates,
            rates,
            rates,
            operator,
            {"readout": operator},
            upper_F,
        )["readout"]

        for expected_values, actual_values in zip(expected[:3], actual[:3]):
            np.testing.assert_allclose(actual_values, expected_values, atol=1e-14)


class PhysicalPumpReadoutIntegrationTests(unittest.TestCase):
    def test_physical_pump_returns_every_signal_and_component(self):
        atom = ATOMS["Rb87"]
        absolute_frequency = line_center_frequency_MHz(atom, "D1") + 500.0
        probe = {
            "source": "PumpA1",
            "mode": "nonlinear",
            "pump_name": "PumpA1",
            "pump_intensity_uW_cm2": 20.0,
            "line": "D1",
            "detuning_MHz": 500.0,
            "k_axis": "x",
            "azimuth_deg": 90.0,
            "ellipticity_deg": 0.0,
            "path_length_cm": 2.0,
            "include_scalar": False,
            "include_orientation": True,
            "include_alignment": True,
        }
        species = {
            "label": "A",
            "atom_name": "Rb87",
            "density_cm3": 1.0e10,
            "R_ER": 10.0,
            "n2_coeffs": {
                "D1": {"width": 17.8, "shift": -8.25},
                "D2": {"width": 18.1, "shift": -5.9},
            },
            "q_axis": "z",
            "rf_axis": "x",
            "rf_observable": "Fx",
            "rf_frequencies_hz": np.array([10.0, 20.0]),
            "probe": probe,
        }
        beam = {
            "name": "PumpA1",
            "target_label": "A",
            "target_atom": "Rb87",
            "absolute_frequency_MHz": absolute_frequency,
            "intensity": 20.0,
            "k_axis": "x",
            "pol": "linear z",
            "selected_transition": None,
            "transition_label": "D1 test pump",
        }
        common = {
            "temperature_C": 23.0,
            "n2_pressure_torr": 0.0,
            "static_field_axis": "z",
            "static_field_nT": 0.0,
        }

        result = compute_alkali_system(species, None, [beam], common)["A"]

        self.assertTrue(result["probe_info"]["nonlinear_available"])
        self.assertEqual(result["probe_info"]["pump_name"], "PumpA1")
        self.assertIn("probe_weak_response", result)
        self.assertEqual(
            set(result["probe_response"]),
            {"total", "scalar", "orientation", "alignment"},
        )
        for component in result["probe_response"].values():
            self.assertEqual(set(component), set(_SIGNALS))
            for response in component.values():
                self.assertEqual(response["amplitude"].shape, (2,))
                self.assertTrue(np.isfinite(response["amplitude"]).all())
        def phasor(component, signal):
            response = result["probe_response"][component][signal]
            return response["in_phase"] + 1j * response["quadrature"]

        total = phasor("total", "rotation")
        rank_sum = phasor("orientation", "rotation") + phasor(
            "alignment", "rotation"
        )
        self.assertFalse(
            np.allclose(total, rank_sum, rtol=1e-8, atol=1e-20),
            "The coupled nonlinear total must not be reconstructed by adding "
            "the two counterfactual rank responses.",
        )


if __name__ == "__main__":
    unittest.main()
