import unittest

import numpy as np
import pandas as pd

from alkali_pumping_app.ui.exports import (
    atomic_polarizability_export_dataframe,
    dataframe_to_csv_bytes,
    weak_rf_export_dataframe,
)


class DataExportTests(unittest.TestCase):
    def test_steady_density_csv_preserves_all_complex_elements_and_basis(self):
        from io import BytesIO
        from alkali_pumping_app.ui.time_evolution_data import steady_density_matrix_export_dataframe
        rho = np.array([[.5, .1j, 0], [-.1j, .3, 0], [0, 0, .2]])
        states = [dict(F=1, m=-1), dict(F=1, m=0), dict(F=0, m=0)]
        species = dict(density_matrix=rho[None], ground_states=states, q_axis='z', atom_name='test')
        frame = steady_density_matrix_export_dataframe(species, label='B')
        restored = pd.read_csv(BytesIO(dataframe_to_csv_bytes(frame)))
        self.assertEqual(len(restored), rho.size)
        rebuilt = np.zeros_like(rho)
        rebuilt[restored.row_index, restored.column_index] = restored.rho_real + 1j * restored.rho_imag
        np.testing.assert_array_equal(rebuilt, rho)
        for row in restored.itertuples():
            self.assertEqual((row.row_F, row.row_m), (states[row.row_index]['F'], states[row.row_index]['m']))
            self.assertEqual((row.column_F, row.column_m), (states[row.column_index]['F'], states[row.column_index]['m']))
        self.assertEqual(set(restored.quantization_axis), {'z'})
        self.assertEqual(set(restored.alkali), {'B'})
        self.assertEqual(set(restored.normalization), {'per total atom'})

    def test_atomic_polarizability_export_has_one_curve_per_column_and_units_row(self):
        plotted = pd.DataFrame(
            {
                "Detuning (MHz)": [-1.0, 1.0, -1.0, 1.0],
                "Series": ["F=1"] * 4,
                "Part": ["Real (phase)"] * 2 + ["Imaginary (attenuation)"] * 2,
                "Polarizability (10^-18 cm^3)": [1.0, 2.0, 3.0, 4.0],
                "Component": ["alpha_gt"] * 4,
            }
        )

        exported = atomic_polarizability_export_dataframe(plotted)

        self.assertEqual(
            list(exported.columns),
            [
                "Laser detuning",
                "α_gt | F=1 | Real (phase)",
                "α_gt | F=1 | Imaginary (attenuation)",
            ],
        )
        self.assertEqual(list(exported.iloc[0]), ["MHz", "10⁻¹⁸ cm³", "10⁻¹⁸ cm³"])
        np.testing.assert_allclose(exported.iloc[1:, 1].astype(float), [1.0, 2.0])

    def test_csv_download_uses_utf8_bom_and_preserves_unicode_headers(self):
        dataframe = pd.DataFrame({"ν_m (Hz)": [1.25], "Pₘ": [0.375]})

        payload = dataframe_to_csv_bytes(dataframe)

        self.assertTrue(payload.startswith(b"\xef\xbb\xbf"))
        self.assertIn("ν_m (Hz),Pₘ", payload.decode("utf-8-sig"))

    def test_weak_rf_export_contains_raw_and_plotted_values(self):
        dataframe = weak_rf_export_dataframe(
            frequencies_hz=[10.0, 20.0],
            susceptibility_amplitude=[1.0, 2.0],
            susceptibility_in_phase=[0.25, -0.5],
            susceptibility_quadrature=[-0.75, 1.0],
            plotted_amplitude=[3.0, 6.0],
            plotted_in_phase=[-0.75, 1.5],
            plotted_quadrature=[2.25, -3.0],
            in_phase_plot_factor=-1.0,
            quadrature_plot_factor=-1.0,
            relaxation_normalized=True,
            normalization_gamma_s_inv=3.0,
        )

        np.testing.assert_allclose(dataframe["frequency_Hz"], [10.0, 20.0])
        np.testing.assert_allclose(
            dataframe["in_phase_raw_hbar_s_per_atom"], [0.25, -0.5]
        )
        np.testing.assert_allclose(dataframe["in_phase_plotted"], [-0.75, 1.5])
        self.assertEqual(set(dataframe["plotted_units"]), {"hbar/atom"})
        self.assertEqual(set(dataframe["in_phase_plot_factor"]), {-1.0})
        self.assertEqual(set(dataframe["quadrature_plot_factor"]), {-1.0})
        self.assertEqual(set(dataframe["relaxation_normalized"]), {True})
        self.assertEqual(set(dataframe["normalization_gamma_s_inv"]), {3.0})
        self.assertEqual(set(dataframe["density_factored"]), {False})
        self.assertTrue(np.isnan(dataframe.loc[0, "density_cm3"]))

    def test_weak_rf_export_marks_unnormalized_plot_units(self):
        dataframe = weak_rf_export_dataframe(
            frequencies_hz=[10.0],
            susceptibility_amplitude=[1.0],
            susceptibility_in_phase=[0.25],
            susceptibility_quadrature=[-0.75],
            plotted_amplitude=[1.0],
            plotted_in_phase=[-0.25],
            plotted_quadrature=[0.75],
            relaxation_normalized=False,
        )

        self.assertEqual(dataframe.loc[0, "plotted_units"], "hbar s/atom")
        self.assertFalse(dataframe.loc[0, "relaxation_normalized"])
        self.assertTrue(np.isnan(dataframe.loc[0, "normalization_gamma_s_inv"]))

    def test_weak_rf_export_records_density_factor_and_units(self):
        dataframe = weak_rf_export_dataframe(
            frequencies_hz=[10.0],
            susceptibility_amplitude=[1.0],
            susceptibility_in_phase=[0.25],
            susceptibility_quadrature=[-0.75],
            plotted_amplitude=[4.2e9],
            plotted_in_phase=[-1.05e9],
            plotted_quadrature=[3.15e9],
            relaxation_normalized=False,
            density_factored=True,
            density_cm3=4.2e9,
        )

        self.assertEqual(dataframe.loc[0, "plotted_units"], "hbar s/cm^3")
        self.assertTrue(dataframe.loc[0, "density_factored"])
        self.assertEqual(dataframe.loc[0, "density_cm3"], 4.2e9)

    def test_weak_rf_export_rejects_mismatched_sample_counts(self):
        with self.assertRaisesRegex(ValueError, "same length"):
            weak_rf_export_dataframe(
                frequencies_hz=[10.0, 20.0],
                susceptibility_amplitude=[1.0],
                susceptibility_in_phase=[0.25, 0.5],
                susceptibility_quadrature=[-0.75, -1.0],
                plotted_amplitude=[1.0, 2.0],
                plotted_in_phase=[-0.25, -0.5],
                plotted_quadrature=[0.75, 1.0],
                relaxation_normalized=False,
            )


if __name__ == "__main__":
    unittest.main()
