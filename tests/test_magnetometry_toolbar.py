"""ACC manual execution, deferred sections, and independent page preferences."""
from contextlib import ExitStack
import unittest
from unittest.mock import patch

import streamlit as st

from alkali_pumping_app.ui import magnetometry_cache as cache
from tests.streamlit_fragments import FragmentSession
from tests import test_magnetometry_responsiveness as ui_helpers


class MagnetometryToolbarTests(unittest.TestCase):
    setUp = ui_helpers.MagnetometryUiTests.setUp
    app = ui_helpers.MagnetometryUiTests.app

    def test_toolbar_and_manual_apply_defer_all_physics_and_exports(self):
        app = self.app()
        self.assertTrue(app.checkbox(key="acc_auto_update").value)
        self.assertTrue(app.button(key="acc_apply_conditions").disabled)
        self.assertFalse(any(widget.label == "Auto update" for widget in app.sidebar.checkbox))
        app.checkbox(key="acc_auto_update").uncheck().run()
        self.assertFalse(app.button(key="acc_apply_conditions").disabled)
        with ExitStack() as stack:
            for name in ("build_species_linear_model", "solve_alkali_equilibrium", "compute_alkali_response"):
                stack.enter_context(patch.object(cache, name, side_effect=AssertionError("Physics before Apply")))
            stack.enter_context(patch.object(st, "download_button", side_effect=AssertionError("Pending CSV")))
            stack.enter_context(patch.object(st, "pyplot", side_effect=AssertionError("Pending plot")))
            for value in (17., 21.):
                app.number_input(key="gamma_ER_A").set_value(value).run()
                self.assertFalse(app.exception)
                self.assertTrue(any("Pending condition changes" in item.value for item in app.info))
                self.assertTrue(any("pump rate: pending Apply" in item.value for item in app.caption))
        cache.equilibrium_cached.clear()
        with patch.object(cache, "solve_alkali_equilibrium", wraps=cache.solve_alkali_equilibrium) as solve:
            app.button(key="acc_apply_conditions").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(solve.call_count, 1)
        self.assertFalse(any("Pending condition changes" in item.value for item in app.info))
        self.assertEqual(len(app.get("download_button")), 3)

    def test_pending_fragment_toggles_and_closed_apply_do_not_solve(self):
        with FragmentSession() as session:
            app = self.app()
            app.checkbox(key="acc_auto_update").uncheck().run()
            app.number_input(key="gamma_ER_A").set_value(24.).run()
            app.session_state["acc_optical_open_A"] = False
            session.run(app, "acc_optical_open_A")
            app.session_state["acc_optical_open_A"] = True
            with patch.object(cache, "compute_alkali_response", side_effect=AssertionError("Pending fragment computed")):
                session.run(app, "acc_optical_open_A")
            self.assertFalse(app.exception)
            self.assertTrue(any("Press Apply" in item.value for item in app.info))
            for section in ("zeeman", "atomic", "optical", "rf_transitions", "transitions"):
                app.session_state[f"acc_{section}_open_A"] = False
            with patch.object(cache, "solve_alkali_equilibrium", side_effect=AssertionError("Closed Apply solved")):
                app.button(key="acc_apply_conditions").click().run()
            self.assertFalse(app.exception)
            app.session_state["acc_zeeman_open_A"] = True
            with patch.object(cache, "compute_alkali_response", side_effect=AssertionError("Closed response computed")):
                session.run(app, "acc_zeeman_open_A")
            self.assertFalse(app.exception)
            self.assertEqual(len(app.get("download_button")), 1)

    def test_readout_edits_remain_immediate_and_auto_update_applies_pending_changes(self):
        app = self.app()
        app.checkbox(key="acc_auto_update").uncheck().run()
        applied = app.session_state["acc_applied_physics"]
        with ExitStack() as stack:
            for name in ("build_species_linear_model", "solve_alkali_equilibrium"):
                stack.enter_context(patch.object(cache, name, side_effect=AssertionError("Readout rebuilt physics")))
            app.number_input(key="_probe_widget_probe_path_length_cm_A").set_value(2.).run()
            app.selectbox(key="_result_widget_rf_observable_A").select("Fy").run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["acc_applied_physics"], applied)
        self.assertFalse(any("Pending condition changes" in item.value for item in app.info))
        app.number_input(key="gamma_ER_A").set_value(28.).run()
        self.assertEqual(app.session_state["acc_applied_physics"], applied)
        app.checkbox(key="acc_auto_update").check().run()
        self.assertFalse(app.exception)
        self.assertNotEqual(app.session_state["acc_applied_physics"], applied)
        self.assertTrue(app.button(key="acc_apply_conditions").disabled)
        self.assertEqual(len(app.get("download_button")), 3)

    def test_preferences_and_pending_edits_survive_navigation_independently(self):
        app = self.app()
        app.checkbox(key="acc_auto_update").uncheck().run()
        applied = app.session_state["acc_applied_physics"]
        app.number_input(key="gamma_ER_A").set_value(31.).run()
        app.switch_page("alkali_pumping_app/pages/time_evolution.py").run()
        self.assertTrue(app.checkbox(key="te_auto_update").value)
        app.checkbox(key="te_auto_update").uncheck().run()
        app.switch_page("alkali_pumping_app/pages/ACC-Magnetometer.py").run()
        self.assertFalse(app.exception)
        self.assertFalse(app.checkbox(key="acc_auto_update").value)
        self.assertEqual(app.session_state["acc_applied_physics"], applied)
        self.assertTrue(any("Pending condition changes" in item.value for item in app.info))
        self.assertEqual(len(app.get("download_button")), 0)
        app.checkbox(key="acc_auto_update").check().run()
        app.switch_page("alkali_pumping_app/pages/time_evolution.py").run()
        self.assertFalse(app.checkbox(key="te_auto_update").value)


if __name__ == "__main__":
    unittest.main()
