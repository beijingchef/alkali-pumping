import unittest

from alkali_pumping_app.physics.constants import ATOMS
from alkali_pumping_app.physics.spectroscopy import (
    relative_hyperfine_transition_label,
)


class TransitionTranslationTests(unittest.TestCase):
    def test_d1_transition_keeps_ground_and_excited_manifold_ranks(self):
        self.assertEqual(
            relative_hyperfine_transition_label(
                ATOMS["Rb87"], ATOMS["Rb85"], "D1", "1→2"
            ),
            "2→3",
        )
        self.assertEqual(
            relative_hyperfine_transition_label(
                ATOMS["Rb87"], ATOMS["Cs133"], "D1", "2→1"
            ),
            "4→3",
        )

    def test_d2_excited_manifold_rank_is_preserved(self):
        self.assertEqual(
            relative_hyperfine_transition_label(
                ATOMS["Rb87"], ATOMS["Rb85"], "D2", "1→0"
            ),
            "2→1",
        )
        self.assertEqual(
            relative_hyperfine_transition_label(
                ATOMS["Cs133"], ATOMS["K39"], "D2", "4→5"
            ),
            "2→3",
        )

    def test_invalid_or_filtered_transition_requests_normal_default(self):
        self.assertIsNone(
            relative_hyperfine_transition_label(
                ATOMS["Rb87"], ATOMS["Rb85"], "D1", "not-a-transition"
            )
        )
        self.assertIsNone(
            relative_hyperfine_transition_label(
                ATOMS["Rb87"], ATOMS["Rb85"], "D2", "1→3", allowed_only=True
            )
        )


if __name__ == "__main__":
    unittest.main()
