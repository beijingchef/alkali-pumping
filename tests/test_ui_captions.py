import unittest

from alkali_pumping_app.ui.captions import input_conditions_caption


class InputConditionsCaptionTests(unittest.TestCase):
    def test_caption_omits_zero_rate_beams(self):
        beams = [
            {
                "name": "Beam 1",
                "line": "D1",
                "transition_label": "1→2",
                "detuning_relative": 0.0,
                "intensity": 7.0,
                "rate": 1200.0,
                "rate_at_resonance": 1200.0,
                "k_axis": "x",
                "pol": "linear z",
            },
            {
                "name": "Beam 2",
                "line": "D1",
                "transition_label": "2→2",
                "detuning_relative": 400.0,
                "intensity": 0.0,
                "rate": 0.0,
                "rate_at_resonance": 400.0,
                "k_axis": "x",
                "pol": "linear z",
            },
        ]

        caption = input_conditions_caption(
            atom_name="Rb87",
            temperature_C=23.5,
            n2_pressure_torr=400.0,
            R_ER=4.0,
            R_SE=25.0,
            q_axis="z",
            bias_larmor_hz=10.0,
            beam_inputs=beams,
        )

        self.assertIn("Sidebar inputs — Rb87", caption)
        self.assertIn("Beam 1", caption)
        self.assertIn("I=7 µW/cm²", caption)
        self.assertIn("R_F,res=1200 s⁻¹", caption)
        self.assertNotIn("Beam 2", caption)

    def test_caption_reports_no_active_pumps(self):
        caption = input_conditions_caption(
            atom_name="Rb87",
            temperature_C=23.5,
            n2_pressure_torr=0.0,
            R_ER=4.0,
            R_SE=0.0,
            q_axis="z",
            bias_larmor_hz=0.0,
            beam_inputs=[],
        )

        self.assertIn("SE=on (R_SE=0", caption)
        self.assertTrue(caption.endswith("Active pumps — none."))
