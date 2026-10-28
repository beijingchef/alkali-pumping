"""Independent dynamics settings use only the current versioned format."""

from copy import deepcopy
import unittest

from alkali_pumping_app.ui.time_evolution_conditions import (
    DEFAULTS, PLOTTED_MOMENTS, PLOTTED_PROBE_SIGNALS, build_settings_payload, parse_settings_payload,
)


class TimeEvolutionConditionTests(unittest.TestCase):
    def test_complete_current_settings_round_trip(self):
        settings = deepcopy(DEFAULTS)
        settings.update(t1_s=0., t2_s=2., intensity_A1=0., probe_source_A="Custom")
        self.assertEqual(parse_settings_payload(build_settings_payload(settings)), settings)
        settings.update(t1_s=30., t2_s=30.)
        self.assertEqual(parse_settings_payload(build_settings_payload(settings)), settings)
        payload = build_settings_payload(DEFAULTS)
        self.assertEqual(payload["version"], "2.7")
        self.assertNotIn("q_axis_A", payload["settings"])
        self.assertNotIn("q_axis_B", payload["settings"])
        self.assertNotIn("moment_types", payload["settings"])
        self.assertNotIn("moment_manifolds", payload["settings"])
        self.assertNotIn("show_allowed_only", payload["settings"])
        self.assertFalse(any(key.startswith("rf_") for key in payload["settings"]))

    def test_missing_or_obsolete_fields_are_rejected_without_default_filling(self):
        for name in DEFAULTS:
            payload = build_settings_payload(DEFAULTS)
            payload["settings"].pop(name)
            with self.subTest(missing=name), self.assertRaisesRegex(ValueError, name):
                parse_settings_payload(payload)
        payload = build_settings_payload(DEFAULTS)
        payload["settings"]["show_allowed_only"] = False
        with self.assertRaisesRegex(ValueError, "Unknown.*show_allowed_only"):
            parse_settings_payload(payload)

    def test_wrong_format_and_future_version_rejected(self):
        for field, value in (("format", "alkali_pumping_conditions"), ("version", "2.4"), ("version", "3.0")):
            payload = build_settings_payload(DEFAULTS)
            payload[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                parse_settings_payload(payload)

    def test_probe_scan_limits_round_trip_and_validate_per_species(self):
        settings = dict(DEFAULTS, probe_scan_min_MHz_A=-10000., probe_scan_max_MHz_A=10000.,
                        probe_scan_min_MHz_B=-100., probe_scan_max_MHz_B=100.)
        self.assertEqual(parse_settings_payload(build_settings_payload(settings)), settings)
        for updates in ({"probe_scan_min_MHz_A": 10000.}, {"probe_scan_max_MHz_A": -10001.},
                        {"probe_scan_min_MHz_B": 100.}, {"probe_scan_max_MHz_B": -95.}):
            with self.subTest(updates=updates), self.assertRaisesRegex(ValueError, "scan limits"):
                build_settings_payload(dict(settings, **updates))

    def test_independent_section_axes_round_trip_and_reject_invalid_axes(self):
        settings = dict(DEFAULTS, zeeman_axis_A="x", moment_axis_A="y", zeeman_axis_B="y", moment_axis_B="z")
        self.assertEqual(parse_settings_payload(build_settings_payload(settings)), settings)
        for field in ("zeeman_axis_A", "moment_axis_A", "zeeman_axis_B", "moment_axis_B"):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, field):
                build_settings_payload(dict(settings, **{field: "invalid"}))
        with self.assertRaisesRegex(ValueError, "Unknown.*q_axis_A"):
            build_settings_payload(dict(settings, q_axis_A="z"))

    def test_moment_selection_round_trips_including_empty_and_all(self):
        for selection in ([], list(PLOTTED_MOMENTS), ["Q_uv", "F_y"]):
            settings = dict(DEFAULTS, plotted_moments=selection)
            self.assertEqual(parse_settings_payload(build_settings_payload(settings)), settings)
        for selection in (["Orientation"], ["Alignment"], ["F_x", "F_x"], "F_z"):
            with self.subTest(selection=selection), self.assertRaisesRegex(ValueError, "plotted_moments"):
                build_settings_payload(dict(DEFAULTS, plotted_moments=selection))

    def test_missing_or_obsolete_physics_model_is_rejected(self):
        for model in (None, "v2-nonlinear-populations-with-uniform-relaxation"):
            payload = build_settings_payload(DEFAULTS)
            payload["physics_model"] = model
            with self.assertRaisesRegex(ValueError, "physics model"):
                parse_settings_payload(payload)

    def test_probe_signal_selections_round_trip_and_validate_independently(self):
        for label in ("A", "B"):
            field = f"plotted_probe_signals_{label}"
            self.assertEqual(DEFAULTS[field], ["Absorption", "Optical rotation"])
            for selection in ([], list(PLOTTED_PROBE_SIGNALS), ["Ellipticity"],
                              ["Optical rotation", "Ellipticity", "Absorption"]):
                settings = dict(DEFAULTS, **{field: selection})
                self.assertEqual(parse_settings_payload(build_settings_payload(settings)), settings)
            for selection in (None, "Absorption", ["Transmission"], ["Absorption", "Absorption"]):
                with self.subTest(label=label, selection=selection), self.assertRaisesRegex(ValueError, field):
                    build_settings_payload(dict(DEFAULTS, **{field: selection}))

    def test_invalid_fields_and_physical_probe_sources_rejected(self):
        for key, value in (("t1_s", 5.), ("t2_s", 0.), ("t2_s", 30.001), ("intensity_A1", -1.),
                           ("probe_source_A", "PumpA1"), ("probe_path_length_cm_A", -1.),
                           ("gamma_ER_A", float("nan")), ("gamma_u_A", -1.),
                           ("gamma_u_B", float("inf")), ("atom_A_name", "invalid"),
                           ("pol_A1", "linear z"), ("plotted_moments", ["RF amplitude"])):
            settings = deepcopy(DEFAULTS)
            settings[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                build_settings_payload(settings)

    def test_scan_signal_selections_are_saved_and_validated_independently(self):
        for label in ('A', 'B'):
            field = f'plotted_scan_signals_{label}'
            self.assertEqual(DEFAULTS[field], ['Optical rotation', 'Absorption'])
            for selected in ([], ['Ellipticity'], ['Absorption', 'Ellipticity', 'Optical rotation']):
                settings = dict(DEFAULTS, **{field: selected})
                self.assertEqual(parse_settings_payload(build_settings_payload(settings)), settings)
            for invalid in (None, 'Absorption', ['Transmission'], ['Absorption', 'Absorption']):
                with self.subTest(label=label, invalid=invalid), self.assertRaisesRegex(ValueError, field):
                    build_settings_payload(dict(DEFAULTS, **{field: invalid}))
        old = build_settings_payload(DEFAULTS)
        old['version'] = '2.6'
        with self.assertRaisesRegex(ValueError, 'version'):
            parse_settings_payload(old)


if __name__ == "__main__":
    unittest.main()
