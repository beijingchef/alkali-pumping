import unittest

import pandas as pd

from alkali_pumping_app.ui.tables import (
    ZEEMAN_COLUMN_LABELS,
    ZEEMAN_DEFAULT_VISIBLE_COLUMN_KEYS,
    render_zeeman_properties_table_html,
)


class ZeemanTableTests(unittest.TestCase):
    def test_default_visible_columns_omit_manifold_population_and_repopulation_rate(self):
        self.assertNotIn("P_F", ZEEMAN_DEFAULT_VISIBLE_COLUMN_KEYS)
        self.assertNotIn("Λ (s⁻¹)", ZEEMAN_DEFAULT_VISIBLE_COLUMN_KEYS)
        self.assertIn("Pₘ", ZEEMAN_DEFAULT_VISIBLE_COLUMN_KEYS)

    def test_column_selector_labels_are_symbols_without_descriptions(self):
        self.assertEqual(ZEEMAN_COLUMN_LABELS["P_F"], r"$P_F$")
        self.assertEqual(ZEEMAN_COLUMN_LABELS["ν^{VS} (Hz)"], r"$\nu^{\mathrm{VS}}$")
        self.assertEqual(
            ZEEMAN_COLUMN_LABELS["G^{SE,self} (s^-1)"],
            r"$G^{\mathrm{SE}}_{\mathrm{self}}$",
        )
        self.assertEqual(ZEEMAN_COLUMN_LABELS["Γ/2π (Hz)"], r"$\Gamma/(2\pi)$")

    def test_relaxation_headers_have_requested_order_and_totals(self):
        row = {
            "F": 2.0,
            "m": 2.0,
            "P_F": 0.5,
            "Pₘ": 0.2,
            "Dₘ": 0.1,
            "ν^{VS} (Hz)": 1.0,
            "ν^{TS} (Hz)": 2.0,
            "ν^{LS} (Hz)": 3.0,
            "ν^{B} (Hz)": 4.0,
            "ν_m (Hz)": 5.0,
            "Λ (s⁻¹)": 6.0,
            "G^{OP} (s^-1)": 7.0,
            "Γ^{OP} (s^-1)": 8.0,
            "G^{ER} (s^-1)": 9.0,
            "Γ^{ER} (s^-1)": 10.0,
            "G^{SE,self} (s^-1)": 4.0,
            "Γ^{SE,self} (s^-1)": 5.0,
            "G^{SE,cross} (s^-1)": 7.0,
            "Γ^{SE,cross} (s^-1)": 7.0,
            "G^{SE} (s^-1)": 11.0,
            "Γ^{SE} (s^-1)": 12.0,
            "G (s^-1)": 27.0,
            "Γ (s^-1)": 30.0,
            "Γ/2π (Hz)": 30.0 / (2.0 * 3.141592653589793),
        }

        html = render_zeeman_properties_table_html(pd.DataFrame([row]))

        ordered_headers = [
            "G<sup>OP</sup>",
            "Γ<sup>OP</sup>",
            "G<sup>ER</sup>",
            "Γ<sup>ER</sup>",
            "zeeman-scripted-symbol'>G<span class='zeeman-script-stack'><span>SE</span><span>self</span>",
            "zeeman-scripted-symbol'>Γ<span class='zeeman-script-stack'><span>SE</span><span>self</span>",
            "zeeman-scripted-symbol'>G<span class='zeeman-script-stack'><span>SE</span><span>cross</span>",
            "zeeman-scripted-symbol'>Γ<span class='zeeman-script-stack'><span>SE</span><span>cross</span>",
            "G<sup>SE</sup>",
            "Γ<sup>SE</sup>",
            "<div class='zeeman-header-quantity'>G</div>",
            "<div class='zeeman-header-quantity'>Γ</div>",
            "Γ/2π",
        ]
        positions = [html.index(header) for header in ordered_headers]

        self.assertEqual(positions, sorted(positions))
        for old_two_decimal_value in (
            "8.00", "9.00", "10.00", "4.00", "5.00", "7.00", "11.00", "12.00"
        ):
            self.assertNotIn(f"<td>{old_two_decimal_value}</td>", html)
        for one_decimal_value in (
            "8.0", "9.0", "10.0", "4.0", "5.0", "7.0", "11.0", "12.0"
        ):
            self.assertIn(f"<td>{one_decimal_value}</td>", html)

    def test_visible_columns_filter_keeps_f_and_m_in_canonical_order(self):
        frame = pd.DataFrame([{
            "F": 2.0, "m": 1.0, "P_F": 0.5, "Pₘ": 0.2,
            "G (s^-1)": 27.0,
        }])

        html = render_zeeman_properties_table_html(
            frame, ["G (s^-1)", "Pₘ"]
        )

        self.assertIn("<div class='zeeman-header-quantity'>F</div>", html)
        self.assertIn("<div class='zeeman-header-quantity'>m</div>", html)
        self.assertIn("P<sub>m</sub>", html)
        self.assertIn("<div class='zeeman-header-quantity'>G</div>", html)
        self.assertNotIn("P<sub>F</sub>", html)
        self.assertLess(html.index("P<sub>m</sub>"), html.index("<div class='zeeman-header-quantity'>G</div>"))


if __name__ == "__main__":
    unittest.main()
