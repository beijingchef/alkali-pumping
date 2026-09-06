import unittest

import numpy as np

from alkali_pumping_app.physics import ATOMS, build_ground_states, compute_alkali_system
from alkali_pumping_app.physics.nonlinear_readout import (
    _SIGNALS,
    propagate_uniform_atomic_feedback,
)
from alkali_pumping_app.physics.optical_pumping import (
    build_optical_L,
    optical_pumping_stokes_liouvillian_sources,
    optical_rate_scale_from_intensity,
)
from alkali_pumping_app.physics.rf_response import (
    _spin_operator,
    weak_drive_matrix_readouts,
    weak_rf_matrix_susceptibility,
)
from alkali_pumping_app.physics.spectroscopy import line_center_frequency_MHz


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


class DissipativeStokesSourceTests(unittest.TestCase):
    def test_sources_are_hermitian_and_trace_preserving(self):
        atom = ATOMS["Rb87"]
        states = build_ground_states(atom)
        population = np.linspace(1.0, 2.0, len(states))
        population /= population.sum()
        sources = optical_pumping_stokes_liouvillian_sources(
            atom=atom,
            line="D1",
            ground_states=states,
            detuning_MHz=500.0,
            intensity_uW_cm2=20.0,
            k_axis="x",
            q_axis="z",
            input_stokes=np.array([-1.0, 0.0, 0.0]),
            density_matrix=np.diag(population),
            n2_pressure_torr=0.0,
            temperature_C=23.0,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
        )

        self.assertEqual(
            set(sources), {"transmission", "s1", "s2", "s3"}
        )
        for source in sources.values():
            np.testing.assert_allclose(source, source.conj().T, atol=1e-13)
            self.assertAlmostEqual(float(np.real(np.trace(source))), 0.0, places=12)
            self.assertAlmostEqual(float(np.imag(np.trace(source))), 0.0, places=12)
        self.assertGreater(np.linalg.norm(sources["s2"]), 0.0)
        self.assertGreater(np.linalg.norm(sources["s3"]), 0.0)

    def test_fractional_intensity_source_matches_population_generator(self):
        atom = ATOMS["Rb87"]
        states = build_ground_states(atom)
        population = np.linspace(1.0, 2.0, len(states))
        population /= population.sum()
        rate_scale = optical_rate_scale_from_intensity(
            atom=atom,
            line="D1",
            intensity_uW_cm2=20.0,
            n2_pressure_torr=0.0,
            temperature_C=23.0,
            n2_width_MHz_per_torr=17.8,
        )
        generator, _info = build_optical_L(
            atom=atom,
            line="D1",
            ground_states=states,
            detuning_MHz=500.0,
            pump_rate_s=rate_scale,
            selected_transition=None,
            k_axis="x",
            pol="linear z",
            q_axis="z",
            n2_pressure_torr=0.0,
            temperature_C=23.0,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
            normalize_to_selected_total=False,
        )
        sources = optical_pumping_stokes_liouvillian_sources(
            atom=atom,
            line="D1",
            ground_states=states,
            detuning_MHz=500.0,
            intensity_uW_cm2=20.0,
            k_axis="x",
            q_axis="z",
            input_stokes=np.array([-1.0, 0.0, 0.0]),
            density_matrix=np.diag(population),
            n2_pressure_torr=0.0,
            temperature_C=23.0,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
        )

        np.testing.assert_allclose(
            np.real(np.diag(sources["transmission"])),
            generator @ population,
            rtol=1e-12,
            atol=1e-13,
        )


class UniformAtomicFeedbackTests(unittest.TestCase):
    @staticmethod
    def _base(**overrides):
        values = {
            signal: np.array([0.0j, 0.0j]) for signal in _SIGNALS
        }
        values.update(overrides)
        return values

    @staticmethod
    def _feedback():
        return {
            signal: {
                coordinate: np.array([0.0j, 0.0j])
                for coordinate in ("transmission", "s1", "s2", "s3")
            }
            for signal in _SIGNALS
        }

    def test_zero_feedback_reproduces_direct_response(self):
        base = self._base(rotation=np.array([1.0 + 2.0j, -3.0j]))
        propagated, spectral_radius = propagate_uniform_atomic_feedback(
            base, self._feedback()
        )

        for signal in _SIGNALS:
            np.testing.assert_allclose(propagated[signal], base[signal])
        np.testing.assert_allclose(spectral_radius, 0.0)

    def test_path_average_halves_the_closed_loop_gain(self):
        base = self._base(s1=np.array([1.0 + 0.0j, 2.0 + 0.0j]))
        feedback = self._feedback()
        feedback["s1"]["s1"] = np.array([0.4 + 0.0j, 0.4 + 0.0j])

        propagated, spectral_radius = propagate_uniform_atomic_feedback(
            base, feedback
        )

        np.testing.assert_allclose(propagated["s1"], base["s1"] / 0.8)
        np.testing.assert_allclose(spectral_radius, 0.2)

    def test_path_average_drives_cross_signal_response(self):
        base = self._base(s1=np.array([2.0 + 0.0j, 4.0 + 0.0j]))
        feedback = self._feedback()
        feedback["rotation"]["s1"] = np.array([3.0 + 0.0j, 3.0 + 0.0j])

        propagated, spectral_radius = propagate_uniform_atomic_feedback(
            base, feedback
        )

        np.testing.assert_allclose(propagated["rotation"], np.array([3.0, 6.0]))
        np.testing.assert_allclose(spectral_radius, 0.0)


class PhysicalPumpReadoutIntegrationTests(unittest.TestCase):
    def test_physical_pump_uses_one_path_averaged_atomic_response(self):
        atom = ATOMS["Rb87"]
        absolute_frequency = line_center_frequency_MHz(atom, "D1") + 500.0
        probe = {
            "source": "PumpA1",
            "mode": "physical",
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

        self.assertTrue(result["probe_info"]["physical_pump_available"])
        self.assertEqual(result["probe_info"]["pump_name"], "PumpA1")
        self.assertIn("probe_weak_response", result)
        self.assertEqual(
            result["probe_info"]["atomic_state_spatial_model"],
            "one density matrix shared by the full cell",
        )
        self.assertEqual(
            result["probe_info"]["stokes_atomic_feedback"],
            "optical-path average",
        )
        self.assertFalse(result["probe_info"]["input_intensity_adjusted"])
        self.assertEqual(
            set(result["probe_info"]["feedback_coordinates"]),
            {"transmission", "s1", "s2", "s3"},
        )
        self.assertIn(
            "dissipative optical pumping",
            result["probe_info"]["feedback_contributions"],
        )
        self.assertEqual(
            set(result["probe_response"]),
            {"total", "scalar", "orientation", "alignment"},
        )
        for component in result["probe_response"].values():
            self.assertEqual(set(component), set(_SIGNALS))
            for response in component.values():
                self.assertEqual(response["amplitude"].shape, (2,))
                self.assertTrue(np.isfinite(response["amplitude"]).all())
        self.assertGreater(
            result["probe_info"]["liouvillian_source_diagnostics"]["s2"][
                "dissipative_frobenius_norm_s_inv"
            ],
            0.0,
        )

        def phasor(component, signal):
            response = result["probe_response"][component][signal]
            return response["in_phase"] + 1j * response["quadrature"]

        total = phasor("total", "rotation")
        rank_sum = phasor("orientation", "rotation") + phasor(
            "alignment", "rotation"
        )
        self.assertGreater(
            np.max(np.abs(total - rank_sum)),
            1e-20,
        )


if __name__ == "__main__":
    unittest.main()
