import unittest

from alkali_pumping_app.ui.text import (
    accent_caption,
    inactive_aware_label,
    line_center_detuning_caption,
    pump_rate_caption,
)


class UiTextTests(unittest.TestCase):
    def test_pump_rate_caption_uses_fixed_notation_at_small_and_large_rates(self):
        self.assertEqual(pump_rate_caption(0), r"$R_{\mathrm{OP}}$: 0 s$^{-1}$")
        self.assertEqual(pump_rate_caption(0.000012345), r"$R_{\mathrm{OP}}$: 0.0000123 s$^{-1}$")
        self.assertEqual(pump_rate_caption(123456), r"$R_{\mathrm{OP}}$: 123000 s$^{-1}$")

    def test_only_inactive_widget_label_is_gray(self):
        self.assertEqual(inactive_aware_label("Alkali B", True), ":gray[Alkali B]")
        self.assertEqual(inactive_aware_label("Alkali B", False), "Alkali B")

    def test_active_captions_use_an_obvious_accent(self):
        self.assertEqual(accent_caption("Target rate"), ":blue[Target rate]")
        self.assertEqual(
            accent_caption("Ignored duplicate", color="orange"),
            ":orange[Ignored duplicate]",
        )

    def test_line_center_detuning_caption_shows_reference_line_and_sign(self):
        self.assertEqual(
            line_center_detuning_caption("D1", 1234.5),
            ":blue[δν from D1 center: +1234 MHz]",
        )
        self.assertEqual(
            line_center_detuning_caption("D2", -25.0),
            ":blue[δν from D2 center: -25 MHz]",
        )
        self.assertEqual(
            line_center_detuning_caption("D1", 1234.5, pump=True),
            r":blue[$\delta\nu_{\mathrm{D1}}$: +1234 MHz]",
        )


if __name__ == "__main__":
    unittest.main()
