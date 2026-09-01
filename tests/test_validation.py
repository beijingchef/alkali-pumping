import unittest

import numpy as np

from alkali_pumping_app.physics.validation import (
    density_matrix_diagnostics,
    population_diagnostics,
)


class DensityMatrixValidationTests(unittest.TestCase):
    def test_valid_density_matrix(self):
        rho = np.array([[0.6, 0.1j], [-0.1j, 0.4]], dtype=complex)
        result = density_matrix_diagnostics(rho)
        self.assertTrue(result.valid)
        self.assertAlmostEqual(result.trace.real, 1.0)
        self.assertGreaterEqual(result.minimum_eigenvalue, 0.0)

    def test_nonhermitian_matrix_is_rejected_by_diagnostic(self):
        rho = np.array([[0.5, 0.2], [0.0, 0.5]], dtype=complex)
        result = density_matrix_diagnostics(rho)
        self.assertFalse(result.hermitian_ok)
        self.assertFalse(result.valid)

    def test_negative_population_is_not_positive(self):
        result = population_diagnostics([1.1, -0.1])
        self.assertFalse(result.positive_semidefinite_ok)
        self.assertFalse(result.valid)

    def test_population_trace(self):
        result = population_diagnostics([0.25, 0.75])
        self.assertTrue(result.trace_ok)
        self.assertTrue(result.hermitian_ok)
        self.assertTrue(result.positive_semidefinite_ok)

    def test_non_square_input_raises(self):
        with self.assertRaises(ValueError):
            density_matrix_diagnostics(np.ones((2, 3)))


if __name__ == "__main__":
    unittest.main()

