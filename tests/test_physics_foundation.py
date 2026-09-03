import unittest

import numpy as np

from alkali_pumping_app.physics import (
    ATOMS,
    allowed_F,
    build_ER_matrix,
    build_excited_states,
    build_ground_states,
    build_optical_L,
    decompose_nu_LS_components,
    add_total_relaxation_columns,
    optical_rate_scale_from_intensity,
    steady_state_from_L,
    transition_shift_MHz,
    field_nT_from_upper_larmor_frequency,
    upper_larmor_frequency_from_field_nT,
)
from alkali_pumping_app.version import DISPLAY_VERSION, __version__


class PhysicsFoundationTests(unittest.TestCase):
    def test_version_metadata(self):
        self.assertEqual(__version__, "6.9.4")
        self.assertEqual(DISPLAY_VERSION, "6.9.4")

    def test_static_field_and_larmor_frequency_conversion_round_trip(self):
        for atom_name in ATOMS:
            frequency = upper_larmor_frequency_from_field_nT(atom_name, 123.4)
            field = field_nT_from_upper_larmor_frequency(atom_name, frequency)
            self.assertAlmostEqual(field, 123.4, places=10)

    def test_rb87_ground_basis_is_preserved(self):
        atom = ATOMS["Rb87"]
        self.assertEqual(allowed_F(atom["I"], atom["ground"]["J"]), [1.0, 2.0])
        states = build_ground_states(atom)
        self.assertEqual(len(states), 8)
        self.assertEqual([state["label"] for state in states[:3]], ["F=1, m=-1", "F=1, m=0", "F=1, m=1"])

    def test_er_map_preserves_total_population(self):
        atom = ATOMS["Rb87"]
        states = build_ground_states(atom)
        matrix = build_ER_matrix(atom, states)
        np.testing.assert_allclose(matrix.sum(axis=0), 1.0, atol=1e-12)

    def test_steady_state_solver_normalizes_population(self):
        generator = np.array([[-2.0, 1.0], [2.0, -1.0]])
        population = steady_state_from_L(generator)
        np.testing.assert_allclose(population, [1.0 / 3.0, 2.0 / 3.0], atol=1e-12)
        self.assertAlmostEqual(float(population.sum()), 1.0)

    def test_pump_rate_reference_can_be_resonance_center(self):
        atom = ATOMS["Rb87"]
        states = build_ground_states(atom)
        selected = {"Fg": 1.0, "Fe": 2.0}
        selected_ground = next(state for state in states if state["F"] == 1.0)
        selected_excited = next(
            state for state in build_excited_states(atom, "D1")
            if state["F"] == 2.0
        )
        center = transition_shift_MHz(selected_ground, selected_excited)
        common = dict(
            atom=atom,
            line="D1",
            ground_states=states,
            detuning_MHz=center + 300.0,
            pump_rate_s=100.0,
            selected_transition=selected,
            k_axis="z",
            pol="linear z",
            q_axis="z",
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=0.0,
            n2_shift_MHz_per_torr=0.0,
        )
        _, detuned_info = build_optical_L(**common)
        _, center_info = build_optical_L(
            **common, reference_at_resonance_center=True
        )

        indices = np.ix_(
            detuned_info["reference_ground_indices"],
            detuned_info["reference_excited_indices"],
        )
        # The rate-reference choice changes only the intensity scale. Both
        # modes must use the same laser frequency and transition detunings.
        np.testing.assert_allclose(
            detuned_info["delta_ge_MHz"],
            center_info["delta_ge_MHz"],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            detuned_info["delta_ge_MHz"][indices],
            300.0,
            rtol=0.0,
            atol=1e-12,
        )
        detuned_total = detuned_info["R_ge"][indices].sum()
        center_referenced_total = center_info["R_ge"][indices].sum()
        self.assertAlmostEqual(float(detuned_total), 100.0, places=9)
        self.assertLess(float(center_referenced_total), 100.0)

    def test_rb87_d1_total_rate_sums_over_ground_zeeman_sublevels(self):
        atom = ATOMS["Rb87"]
        states = build_ground_states(atom)
        results = {}

        for ground_F in (1.0, 2.0):
            selected = {"Fg": ground_F, "Fe": 2.0}
            selected_ground = next(
                state for state in states if state["F"] == ground_F
            )
            selected_excited = next(
                state for state in build_excited_states(atom, "D1")
                if state["F"] == 2.0
            )
            center = transition_shift_MHz(selected_ground, selected_excited)
            _, info = build_optical_L(
                atom=atom,
                line="D1",
                ground_states=states,
                detuning_MHz=center,
                pump_rate_s=100.0,
                selected_transition=selected,
                k_axis="z",
                pol="sigma+",
                q_axis="z",
                n2_pressure_torr=0.0,
                temperature_C=23.5,
                n2_width_MHz_per_torr=0.0,
                n2_shift_MHz_per_torr=0.0,
            )
            indices = np.ix_(
                info["reference_ground_indices"],
                info["reference_excited_indices"],
            )
            selected_rates = info["R_ge"][indices]
            self.assertAlmostEqual(float(selected_rates.sum()), 100.0, places=9)
            self.assertAlmostEqual(
                float(selected_rates.sum(axis=1).mean()),
                100.0 / (2.0 * ground_F + 1.0),
                places=9,
            )
            results[ground_F] = info["normalization_scale"]

        # The two Rb87 D1 F=1,2 -> F'=2 transitions have equal total
        # Zeeman-summed strength, so equal total rates represent equal intensity.
        self.assertAlmostEqual(results[1.0], results[2.0], places=12)

    def test_optical_rate_scale_is_linear_in_intensity(self):
        atom = ATOMS["Rb87"]
        common = dict(
            atom=atom,
            line="D1",
            n2_pressure_torr=0.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
        )
        one_uW = optical_rate_scale_from_intensity(
            intensity_uW_cm2=1.0,
            **common,
        )
        seven_uW = optical_rate_scale_from_intensity(
            intensity_uW_cm2=7.0,
            **common,
        )
        self.assertGreater(one_uW, 0.0)
        self.assertAlmostEqual(seven_uW, 7.0 * one_uW, places=10)

    def test_pressure_broadening_reduces_rate_at_fixed_intensity(self):
        atom = ATOMS["Rb87"]
        common = dict(
            atom=atom,
            line="D1",
            intensity_uW_cm2=1.0,
            temperature_C=23.5,
            n2_width_MHz_per_torr=17.8,
        )
        unbroadened = optical_rate_scale_from_intensity(
            n2_pressure_torr=0.0,
            **common,
        )
        broadened = optical_rate_scale_from_intensity(
            n2_pressure_torr=400.0,
            **common,
        )
        self.assertLess(broadened, unbroadened)

    def test_light_shift_decomposition_returns_state_contributions(self):
        states = [
            {"F": 1.0, "m": -1.0},
            {"F": 1.0, "m": 0.0},
            {"F": 1.0, "m": 1.0},
        ]
        # scalar=10, vector=2*m, tensor=3*[3m^2-F(F+1)]
        total = np.array([11.0, 4.0, 15.0])
        components = decompose_nu_LS_components(states, total)

        np.testing.assert_allclose(components["scalar"], [10.0, 10.0, 10.0])
        np.testing.assert_allclose(components["vector"], [-2.0, 0.0, 2.0])
        np.testing.assert_allclose(components["tensor"], [3.0, -6.0, 3.0])
        np.testing.assert_allclose(
            components["scalar"] + components["vector"] + components["tensor"],
            total,
            atol=1e-12,
        )

    def test_unavailable_light_shift_components_remain_unavailable(self):
        states = [{"F": 1.0, "m": -1.0}, {"F": 1.0, "m": 0.0}]
        components = decompose_nu_LS_components(states, [np.nan, np.nan])
        self.assertTrue(np.isnan(components["vector"]).all())
        self.assertTrue(np.isnan(components["tensor"]).all())

    def test_total_relaxation_columns_sum_all_components(self):
        import pandas as pd

        component_rates = pd.DataFrame({
            "G_OP": [10.0, 20.0],
            "G_ER": [2.0, 3.0],
            "G_SE": [-1.0, 4.0],
            "Gamma_OP": [8.0, np.nan],
            "Gamma_ER": [1.5, np.nan],
            "Gamma_SE": [0.5, np.nan],
        })

        totals = add_total_relaxation_columns(component_rates)

        np.testing.assert_allclose(totals["G_total"], [11.0, 27.0])
        self.assertAlmostEqual(totals.loc[0, "Gamma_total"], 10.0)
        self.assertAlmostEqual(
            totals.loc[0, "Gamma_total_over_2pi"],
            10.0 / (2.0 * np.pi),
        )
        self.assertTrue(np.isnan(totals.loc[1, "Gamma_total"]))
        self.assertTrue(np.isnan(totals.loc[1, "Gamma_total_over_2pi"]))


if __name__ == "__main__":
    unittest.main()
