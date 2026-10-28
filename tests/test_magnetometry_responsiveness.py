"""Verify unchanged ACC physics, deferred work, and native fragment isolation."""
import base64
from contextlib import ExitStack
from copy import deepcopy
import json
from pathlib import Path
import re
import unittest
from unittest.mock import patch

import numpy as np
import streamlit as st
from streamlit.proto.WidgetStates_pb2 import WidgetState
from streamlit.testing.v1 import AppTest, element_tree

from alkali_pumping_app.physics import system_model as physics
from alkali_pumping_app.ui import magnetometry_cache as cache
from tests.streamlit_fragments import FragmentSession
from tests.test_two_manifold_rf import system_inputs


ROOT = Path(__file__).resolve().parents[1]
PAGE = "alkali_pumping_app/pages/ACC-Magnetometer.py"


class MagnetometryCalculationTests(unittest.TestCase):
    def tearDown(self):
        for function in (cache.models_cached, cache.equilibrium_cached, cache.contexts_cached,
                         cache.atomic_cached, cache.optical_cached):
            function.clear()

    def test_separate_readouts_match_combined_single_dual_and_physical_results(self):
        for dual in (False, True):
            for physical_pump in (False, True):
                with self.subTest(dual=dual, physical=physical_pump):
                    config, beams, common = system_inputs(physical=physical_pump)
                    other = deepcopy(config) if dual else None
                    if other is not None:
                        other.update(label="B", atom_name="Rb85", density_cm3=2e10)
                    combined = physics.compute_alkali_system(config, other, beams, common)
                    physical = [{key: value for key, value in item.items()
                                 if key not in ("probe", "rf_axis", "rf_observable", "rf_frequencies_hz")}
                                if item is not None else None for item in (config, other)]
                    with patch.object(physics, "weak_probe_readout_operators", side_effect=AssertionError("Readout during equilibrium")):
                        equilibrium = physics.compute_alkali_equilibrium(*physical, beams, common)
                    contexts = physics.prepare_alkali_response_contexts(equilibrium, common)
                    for label, item in (("A", config), ("B", other)):
                        if item is None:
                            continue
                        rf = {name: item[name] for name in ("rf_axis", "rf_observable", "rf_frequencies_hz")}
                        with patch.object(physics, "weak_probe_readout_operators", side_effect=AssertionError("Probe during atomic readout")):
                            atomic = physics.compute_alkali_response(equilibrium, label, rf, common, contexts)
                        with patch.object(physics, "atomic_rf_readout_operators", side_effect=AssertionError("Atomic readout during optical")):
                            optical = physics.compute_alkali_response(
                                equilibrium, label, dict(rf, probe=item["probe"]), common, contexts, optical=True)
                        np.testing.assert_allclose(equilibrium[label]["population"], combined[label]["population"], atol=1e-13)
                        np.testing.assert_allclose(equilibrium["J_coupled"], combined["J_coupled"], atol=1e-13)
                        for manifold in ("Fa", "Fb"):
                            for field in ("amplitude", "in_phase", "quadrature"):
                                np.testing.assert_allclose(atomic["rf_response"][manifold][field],
                                                           combined[label]["rf_response"][manifold][field], atol=1e-13)
                        for rank in ("total", "orientation", "alignment"):
                            for signal, values in optical["probe_response"][rank].items():
                                for field in ("amplitude", "in_phase", "quadrature"):
                                    np.testing.assert_allclose(values[field], combined[label]["probe_response"][rank][signal][field], atol=1e-13)

    def test_probe_and_observable_changes_reuse_equilibrium_and_modes(self):
        config, beams, common = system_inputs()
        physical = {name: value for name, value in config.items()
                    if name not in ("probe", "rf_axis", "rf_observable", "rf_frequencies_hz")}
        models = cache.models_cached(physical, None, beams, common)
        equilibrium = cache.equilibrium_cached(models, common)
        rf = {name: config[name] for name in ("rf_axis", "rf_observable", "rf_frequencies_hz")}
        cache.atomic_cached(equilibrium, "A", rf, common)
        with ExitStack() as stack:
            for name in ("build_species_linear_model", "solve_alkali_equilibrium", "prepare_alkali_response_contexts"):
                stack.enter_context(patch.object(cache, name, side_effect=AssertionError("Shared physics recomputed")))
            rf["rf_observable"] = "Fy"
            cache.atomic_cached(cache.equilibrium_cached(models, common), "A", rf, common)
            for detuning in (800., 1400.):
                probe = dict(config["probe"], detuning_MHz=detuning)
                cache.optical_cached(equilibrium, "A", {k: v for k, v in rf.items() if k != "rf_observable"}, probe, common)


class MagnetometryUiTests(unittest.TestCase):
    def setUp(self):
        original = element_tree.get_widget_state

        def widget_state(node):
            state = original(node)
            if state is None and isinstance(node, element_tree.Expander) and node.proto.id:
                return WidgetState(id=node.proto.id, bool_value=node.root._runner.session_state[node.proto.id])
            return state

        adapter = patch.object(element_tree, "get_widget_state", side_effect=widget_state)
        adapter.start()
        self.addCleanup(adapter.stop)

    def app(self):
        app = AppTest.from_file(str(ROOT / "alkali_pumping.py"), default_timeout=60).run()
        self.assertFalse(app.exception)
        return app

    def test_matrices_are_nested_closed_and_defer_dataframe_work(self):
        app = self.app()
        self.assertNotIn("Display", [item.value for item in app.sidebar.header])
        self.assertNotIn("Show rate matrices", [item.label for item in app.checkbox])
        self.assertFalse(app.exception)
        zeeman = next(item for item in app.expander if "ground-state Zeeman" in item.label)
        matrices = next(item for item in zeeman.expander if item.label == "Alkali A rate matrices")
        self.assertFalse(matrices.proto.expanded)
        self.assertFalse(matrices.dataframe)
        app.session_state["acc_rate_matrices_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        zeeman = next(item for item in app.expander if "ground-state Zeeman" in item.label)
        self.assertEqual(len(zeeman.expander[0].dataframe), 6)
        app.session_state["acc_rate_matrices_open_A"] = False
        with patch.object(st, "dataframe", side_effect=AssertionError("Closed matrix rendering")):
            app.run()
        self.assertFalse(app.exception)

    def test_closed_sections_defer_solve_readouts_exports_and_refresh_when_opened(self):
        app = self.app()
        for section in ("zeeman", "atomic", "optical", "transitions", "rf_transitions"):
            app.session_state[f"acc_{section}_open_A"] = False
        cache.equilibrium_cached.clear()
        cache.contexts_cached.clear()
        with ExitStack() as stack:
            for name in ("solve_alkali_equilibrium", "prepare_alkali_response_contexts", "compute_alkali_response"):
                stack.enter_context(patch.object(cache, name, side_effect=AssertionError("Closed physics executed")))
            stack.enter_context(patch.object(st, "download_button", side_effect=AssertionError("Closed export executed")))
            app.number_input(key="gamma_ER_A").set_value(13.5).run()
        self.assertFalse(app.exception)
        self.assertFalse(app.get("imgs"))
        app.session_state["acc_zeeman_open_A"] = True
        with (patch.object(cache, "compute_alkali_response", side_effect=AssertionError("Closed response executed")),
              patch.object(st, "pyplot", wraps=st.pyplot) as plots):
            app.run()
        self.assertFalse(app.exception)
        self.assertEqual(plots.call_count, 1)
        app.session_state["acc_optical_open_A"] = True
        app.run()
        self.assertFalse(app.exception)
        self.assertTrue(app.selectbox(key="_result_widget_probe_signal_A"))

    def test_fragment_controls_are_local_and_save_hidden_choices(self):
        from alkali_pumping_app.ui import downloads
        with FragmentSession() as session:
            app = self.app()
            with patch.object(st, "pyplot", wraps=st.pyplot) as plot:
                app.selectbox(key="_result_widget_rf_observable_A").select("Fy")
                session.run(app, "_result_widget_rf_observable_A")
            self.assertFalse(app.exception)
            self.assertEqual(plot.call_count, 1)
            app.selectbox(key="_result_widget_probe_signal_A").select("ellipticity")
            with patch.object(downloads.st, "html", wraps=downloads.st.html) as html:
                session.run(app, "_result_widget_probe_signal_A")
            self.assertFalse(app.exception)
            updates = [call.args[0] for call in html.call_args_list if call.args and 'button.dataset.settingsSnapshot =' in str(call.args[0])]
            self.assertEqual(len(updates), 1)
            encoded = re.search(r'button.dataset.settingsSnapshot = "([^"]+)"', updates[0]).group(1)
            saved = json.loads(base64.b64decode(encoded))["conditions"]
            self.assertEqual(saved["rf_observable_A"], "Fy")
            self.assertEqual(saved["probe_signal_A"], "ellipticity")
            app.session_state["acc_optical_open_A"] = False
            with patch.object(st, "pyplot", side_effect=AssertionError("Closing section plotted")):
                session.run(app, "acc_optical_open_A")
            self.assertFalse(app.exception)
            app.switch_page("alkali_pumping_app/pages/atomic_properties.py").run()
            app.switch_page(PAGE).run()
            self.assertFalse(app.exception)
            self.assertFalse(app.session_state["acc_optical_open_A"])
            self.assertEqual(app.selectbox(key="_result_widget_rf_observable_A").value, "Fy")
            app.session_state["acc_optical_open_A"] = True
            app.run()
            self.assertEqual(app.selectbox(key="_result_widget_probe_signal_A").value, "ellipticity")


if __name__ == "__main__":
    unittest.main()
