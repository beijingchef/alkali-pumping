import unittest

import numpy as np

from alkali_pumping_app.ui.plot_data import aligned_population_bar_data


class PopulationBarLayoutTests(unittest.TestCase):
    def test_hyperfine_manifolds_share_the_same_numerical_m_axis(self):
        states = [
            {"F": 2.0, "m": -2.0, "E": 1.0},
            {"F": 2.0, "m": -1.0, "E": 1.0},
            {"F": 2.0, "m": 0.0, "E": 1.0},
            {"F": 2.0, "m": 1.0, "E": 1.0},
            {"F": 2.0, "m": 2.0, "E": 1.0},
            {"F": 1.0, "m": -1.0, "E": 0.0},
            {"F": 1.0, "m": 0.0, "E": 0.0},
            {"F": 1.0, "m": 1.0, "E": 0.0},
        ]
        population = np.arange(1.0, 9.0)

        ticks, groups = aligned_population_bar_data(states, population)

        np.testing.assert_allclose(ticks, [-2.0, -1.0, 0.0, 1.0, 2.0])
        self.assertEqual([group[0] for group in groups], [2.0, 1.0])
        np.testing.assert_allclose(groups[0][1], [-2.0, -1.0, 0.0, 1.0, 2.0])
        np.testing.assert_allclose(groups[1][1], [-1.0, 0.0, 1.0])

    def test_population_length_must_match_states(self):
        with self.assertRaisesRegex(ValueError, "Population length"):
            aligned_population_bar_data(
                [{"F": 1.0, "m": 0.0, "E": 0.0}],
                [],
            )


if __name__ == "__main__":
    unittest.main()
