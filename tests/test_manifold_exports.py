import unittest

import numpy as np

from alkali_pumping_app.physics import compute_alkali_system, upper_larmor_frequency_from_field_nT
from alkali_pumping_app.ui.exports import manifold_rf_export_dataframe, weak_probe_export_dataframe
from alkali_pumping_app.ui.rf_display import prepare_weak_rf_plot_values
from tests.test_two_manifold_rf import system_inputs


class ManifoldExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        species, beams, common = system_inputs()
        cls.result = compute_alkali_system(species, None, beams, common)["A"]

    def test_atomic_export_uses_one_scale_and_correct_quadrupole_units(self):
        result = self.result
        gamma = result["rf_relaxation_reference"]["Gamma_m"]
        plotted = {name: prepare_weak_rf_plot_values(
            values["amplitude"], values["in_phase"], values["quadrature"],
            flip_in_phase=True, relaxation_gamma_s_inv=gamma,
        ) for name, values in result["rf_response"].items()}
        frame = manifold_rf_export_dataframe(
            result, plotted, ("Fb",), in_phase_plot_factor=-1.0,
            relaxation_normalized=True, normalization_gamma_s_inv=gamma,
        )
        self.assertEqual(set(frame["manifold"]), {"Fa", "Fb"})
        self.assertEqual(set(frame["plotted_units"]), {"hbar^2/atom"})
        self.assertEqual(set(frame["normalization_gamma_s_inv"]), {gamma})
        for manifold, group in frame.groupby("manifold"):
            np.testing.assert_allclose(group["in_phase_plotted"],
                -gamma * group["in_phase_raw_hbar2_s_per_atom"])
            self.assertTrue((group["selected_for_plot"] == (manifold == "Fb")).all())
            self.assertEqual(set(group["F"]), {result["rf_manifolds"][manifold]})

    def test_probe_export_preserves_coherent_manifold_sum_in_field_units(self):
        result = self.result
        scale = 2.0 * np.pi * upper_larmor_frequency_from_field_nT("Rb87", 1.0)
        frame = weak_probe_export_dataframe(
            result["rf_frequencies_hz"], result["probe_response"], "rotation",
            rf_rabi_rad_s_per_nT=scale, manifold_response=result["probe_manifold_response"],
            manifold_F=result["rf_manifolds"], displayed_manifold="Fb",
            quadrature_plot_factor=-1.0,
        )
        for rank in ("total", "scalar", "orientation", "alignment"):
            for field in ("in_phase", "quadrature"):
                np.testing.assert_allclose(frame[f"{rank}_{field}_per_nT"],
                    frame[f"Fa_{rank}_{field}_per_nT"] + frame[f"Fb_{rank}_{field}_per_nT"],
                    atol=1e-14)
        self.assertEqual(set(frame["units"]), {"rad/nT"})
        self.assertEqual(set(frame["displayed_manifold"]), {"Fb"})
        self.assertEqual(set(frame["quadrature_plot_factor"]), {-1.0})


if __name__ == "__main__":
    unittest.main()
