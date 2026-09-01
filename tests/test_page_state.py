import unittest
from unittest.mock import patch

from alkali_pumping_app.ui import page_state


class _RecordingState(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.assignments = []

    def __setitem__(self, key, value):
        self.assignments.append((key, value))
        super().__setitem__(key, value)


class PageStateTests(unittest.TestCase):
    def test_registered_setting_is_detached_on_later_entrypoint_run(self):
        state = _RecordingState({"ap_line": "D2", "transient_button": False})
        with patch.object(page_state.st, "session_state", state):
            page_state.register_persistent_page_settings(["ap_line"])
            state.assignments.clear()
            page_state.preserve_persistent_page_settings()

        self.assertIn(("ap_line", "D2"), state.assignments)
        self.assertNotIn(("transient_button", False), state.assignments)

    def test_registration_accumulates_settings_from_multiple_pages(self):
        state = _RecordingState()
        with patch.object(page_state.st, "session_state", state):
            page_state.register_persistent_page_settings(["ap_line"])
            page_state.register_persistent_page_settings(["atomic_properties_tab"])

        self.assertEqual(
            state[page_state._REGISTRY_KEY],
            ("ap_line", "atomic_properties_tab"),
        )


if __name__ == "__main__":
    unittest.main()
