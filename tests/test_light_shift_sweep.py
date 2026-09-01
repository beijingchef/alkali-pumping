import unittest

import numpy as np

from alkali_pumping_app.physics import (
    ATOMS,
    build_ground_states,
    build_optical_L,
    calculate_light_shift_sweep,
    dipole_amplitude,
    dipole_strength,
    lab_e_field,
    optical_rate_scale_from_intensity,
    polarization_ellipse_vector,
    preset_ellipse_parameters,
    spherical_weights_from_lab,
    stokes_from_ellipse,
    tensor_geometry_E20,
)


class ArbitraryPolarizationTests(unittest.TestCase):
    def test_every_preset_has_an_equivalent_ellipse_setting(self):
        for k_axis in ("x", "y", "z"):
            for preset in ("sigma+", "sigma-", *[f"linear {axis}" for axis in ("x", "y", "z") if axis != k_axis]):
                azimuth, ellipticity = preset_ellipse_parameters(k_axis, preset)
                ellipse = polarization_ellipse_vector(k_axis, azimuth, ellipticity)
                expected = lab_e_field(k_axis, preset)
                self.assertAlmostEqual(abs(np.vdot(ellipse, expected)), 1.0, places=12)

    def test_ellipse_special_cases_match_existing_presets(self):
        for k_axis in ("x", "y", "z"):
            sigma_plus = polarization_ellipse_vector(k_axis, 0.0, 45.0)
            sigma_minus = polarization_ellipse_vector(k_axis, 0.0, -45.0)
            np.testing.assert_allclose(
                np.abs(np.vdot(sigma_plus, lab_e_field(k_axis, "sigma+"))),
                1.0,
                atol=1e-12,
            )
            np.testing.assert_allclose(
                np.abs(np.vdot(sigma_minus, lab_e_field(k_axis, "sigma-"))),
                1.0,
                atol=1e-12,
            )

    def test_stokes_and_tensor_geometry_are_normalized(self):
        np.testing.assert_allclose(stokes_from_ellipse(19.0, 45.0), [0.0, 0.0, 1.0], atol=1e-12)
        linear_z = polarization_ellipse_vector("x", 90.0, 0.0)
        weights = spherical_weights_from_lab(linear_z, "z")
        self.assertAlmostEqual(weights[0], 1.0, places=12)
        self.assertAlmostEqual(tensor_geometry_E20(linear_z, "z"), -2.0 / np.sqrt(6.0), places=12)
        circular_x = polarization_ellipse_vector("x", 0.0, 45.0)
        self.assertAlmostEqual(
            tensor_geometry_E20(circular_x, "z"),
            -1.0 / (2.0 * np.sqrt(6.0)),
            places=12,
        )

    def test_dipole_amplitude_square_matches_strength(self):
        atom = ATOMS["Rb87"]
        amplitude = dipole_amplitude(atom["I"], 0.5, 0.5, 2.0, 1.0, 2.0, 2.0, 1)
        strength = dipole_strength(atom["I"], 0.5, 0.5, 2.0, 1.0, 2.0, 2.0, 1)
        self.assertAlmostEqual(amplitude**2, strength, places=14)


class LightShiftSweepTests(unittest.TestCase):
    def _sweep(self, E_lab, q_axis="z"):
        return calculate_light_shift_sweep(
            atom=ATOMS["Rb87"],
            line="D1",
            detunings_MHz=np.array([-1000.0, -250.0, 500.0]),
            E_lab=E_lab,
            k_axis="z",
            q_axis=q_axis,
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
        )

    def test_sweep_hamiltonian_is_hermitian_and_components_reconstruct_diagonal(self):
        result = self._sweep(polarization_ellipse_vector("z", 23.0, 17.0))
        matrix = result["hamiltonian_hz_per_uW_cm2"]
        np.testing.assert_allclose(matrix, np.swapaxes(matrix.conj(), 1, 2), atol=1e-12)
        components = result["components_hz_per_uW_cm2"]
        reconstructed = components["scalar"] + components["vector"] + components["tensor"] + components["residual"]
        np.testing.assert_allclose(reconstructed, result["diagonal_hz_per_uW_cm2"], atol=1e-11)
        self.assertFalse(result["diagonal_in_selected_basis"])

    def test_preset_diagonal_matches_existing_optical_diagnostic(self):
        atom = ATOMS["Rb87"]
        states = build_ground_states(atom)
        detuning = 500.0
        rate_scale = optical_rate_scale_from_intensity(
            atom=atom,
            line="D1",
            intensity_uW_cm2=1.0,
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
        )
        _, info = build_optical_L(
            atom=atom,
            line="D1",
            ground_states=states,
            detuning_MHz=detuning,
            pump_rate_s=rate_scale,
            k_axis="z",
            pol="sigma+",
            q_axis="z",
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
            normalize_to_selected_total=False,
        )
        expected = np.sum(info["light_shift_ge_angular"], axis=1) / (2.0 * np.pi)
        sweep = calculate_light_shift_sweep(
            atom=atom,
            line="D1",
            detunings_MHz=np.array([detuning, detuning + 1.0]),
            E_lab=lab_e_field("z", "sigma+"),
            k_axis="z",
            q_axis="z",
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
        )
        np.testing.assert_allclose(sweep["diagonal_hz_per_uW_cm2"][0], expected, rtol=1e-11, atol=1e-12)
        self.assertTrue(sweep["diagonal_in_selected_basis"])

    def test_eigenvalues_are_independent_of_basis_axis_at_zero_field(self):
        common = dict(
            atom=ATOMS["Rb87"],
            line="D1",
            detunings_MHz=np.array([500.0, 501.0]),
            E_lab=lab_e_field("x", "linear z"),
            k_axis="x",
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
        )
        along_polarization = calculate_light_shift_sweep(q_axis="z", **common)
        transverse_basis = calculate_light_shift_sweep(q_axis="x", **common)
        for expected, actual in zip(
            along_polarization["eigenvalues_hz_per_uW_cm2"],
            transverse_basis["eigenvalues_hz_per_uW_cm2"],
        ):
            np.testing.assert_allclose(actual["shifts"], expected["shifts"], atol=1e-12)

    def test_linear_aligned_tensor_shift_has_mathur_m_pattern(self):
        sweep = calculate_light_shift_sweep(
            atom=ATOMS["Rb87"],
            line="D1",
            detunings_MHz=np.array([500.0, 501.0]),
            E_lab=lab_e_field("x", "linear z"),
            k_axis="x",
            q_axis="z",
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
        )
        states = sweep["ground_states"]
        tensor = sweep["components_hz_per_uW_cm2"]["tensor"][0]
        f2 = {
            float(state["m"]): tensor[index]
            for index, state in enumerate(states)
            if np.isclose(float(state["F"]), 2.0)
        }
        K2 = f2[0.0]
        np.testing.assert_allclose(
            [f2[-2.0], f2[-1.0], f2[0.0], f2[1.0], f2[2.0]],
            [-K2, 0.5 * K2, K2, 0.5 * K2, -K2],
            atol=1e-12,
        )


if __name__ == "__main__":
    unittest.main()
