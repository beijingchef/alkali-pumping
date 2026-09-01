import unittest

import matplotlib.pyplot as plt
import numpy as np

from alkali_pumping_app.ui.rf_display import (
    add_probe_decomposition_legend,
    prepare_weak_rf_plot_values,
    rf_component_legend_label,
)


class RfDisplayTests(unittest.TestCase):
    def test_density_and_relaxation_factors_multiply_all_plotted_components(self):
        amplitude, in_phase, quadrature = prepare_weak_rf_plot_values(
            amplitude=[1.0, 2.0],
            in_phase=[0.25, -0.5],
            quadrature=[-0.75, 1.0],
            relaxation_gamma_s_inv=3.0,
            density_cm3=10.0,
        )

        np.testing.assert_allclose(amplitude, [30.0, 60.0])
        np.testing.assert_allclose(in_phase, [7.5, -15.0])
        np.testing.assert_allclose(quadrature, [-22.5, 30.0])

    def test_no_optional_factor_preserves_per_atom_scale(self):
        amplitude, in_phase, quadrature = prepare_weak_rf_plot_values(
            amplitude=[1.0],
            in_phase=[0.25],
            quadrature=[-0.75],
        )

        np.testing.assert_allclose(amplitude, [1.0])
        np.testing.assert_allclose(in_phase, [0.25])
        np.testing.assert_allclose(quadrature, [-0.75])

    def test_pi_shifts_flip_only_the_selected_components(self):
        _amplitude, in_phase, quadrature = prepare_weak_rf_plot_values(
            amplitude=[1.0],
            in_phase=[0.25],
            quadrature=[-0.75],
            flip_in_phase=True,
            flip_quadrature=False,
        )
        np.testing.assert_allclose(in_phase, [-0.25])
        np.testing.assert_allclose(quadrature, [-0.75])

        _amplitude, in_phase, quadrature = prepare_weak_rf_plot_values(
            amplitude=[1.0],
            in_phase=[0.25],
            quadrature=[-0.75],
            flip_in_phase=False,
            flip_quadrature=True,
        )
        np.testing.assert_allclose(in_phase, [0.25])
        np.testing.assert_allclose(quadrature, [0.75])

    def test_pi_shift_updates_each_legend_label(self):
        self.assertEqual(rf_component_legend_label("X", False), "X")
        self.assertEqual(rf_component_legend_label("X", True), "−X")
        self.assertEqual(rf_component_legend_label("Y", False), "Y")
        self.assertEqual(rf_component_legend_label("Y", True), "−Y")

    def test_probe_decomposition_legend_has_contribution_columns(self):
        figure, axis = plt.subplots()
        entries = {}
        for column in ("Total", "Orientation", "Alignment"):
            entries[column] = []
            for row, linestyle in (
                ("Amplitude", "-"),
                ("In phase (X)", "--"),
                ("Quadrature (Y)", ":"),
            ):
                line, = axis.plot([0.0, 1.0], [0.0, 1.0], linestyle=linestyle)
                entries[column].append((line, row))

        legend = add_probe_decomposition_legend(axis, entries)

        self.assertEqual(legend._ncols, 3)
        self.assertEqual(
            [text.get_text() for text in legend.get_texts()],
            [
                "Total", "Amplitude", "In phase (X)", "Quadrature (Y)",
                "Orientation", "Amplitude", "In phase (X)", "Quadrature (Y)",
                "Alignment", "Amplitude", "In phase (X)", "Quadrature (Y)",
            ],
        )
        for index in (0, 4, 8):
            self.assertEqual(legend.get_texts()[index].get_weight(), "bold")
        plt.close(figure)


if __name__ == "__main__":
    unittest.main()
