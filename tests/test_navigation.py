import unittest
from pathlib import Path
import runpy
from types import SimpleNamespace
from unittest.mock import patch

from alkali_pumping_app.ui import page_state


class _Navigation:
    def __init__(self, title):
        self.title = title
        self.was_run = False

    def run(self):
        self.was_run = True


class _Streamlit:
    def __init__(self, selected_title):
        self.selected_title = selected_title
        self.page_config_calls = []
        self.active_navigation = None
        self.navigation_groups = None
        self.navigation_kwargs = None
        self.session_state = {}

    def set_page_config(self, **kwargs):
        self.page_config_calls.append(kwargs)

    @staticmethod
    def Page(path, **kwargs):
        return SimpleNamespace(path=path, **kwargs)

    def navigation(self, groups, **kwargs):
        self.navigation_groups = groups
        self.navigation_kwargs = kwargs
        pages = [page for group in groups.values() for page in group]
        selected = next(page for page in pages if page.title == self.selected_title)
        self.active_navigation = _Navigation(selected.title)
        return self.active_navigation


class NavigationTests(unittest.TestCase):
    def test_browser_title_uses_active_navigation_page(self):
        entrypoint = Path(__file__).resolve().parents[1] / "alkali_pumping.py"

        for page_title in (
            "ACC Magnetometer",
            "Light shift",
            "Atomic polarizability",
            "Atomic properties",
        ):
            with self.subTest(page_title=page_title):
                fake_streamlit = _Streamlit(page_title)
                with (
                    patch.dict("sys.modules", {"streamlit": fake_streamlit}),
                    patch.object(
                        page_state.st,
                        "session_state",
                        fake_streamlit.session_state,
                    ),
                ):
                    runpy.run_path(str(entrypoint), run_name="__navigation_test__")

                self.assertEqual(
                    fake_streamlit.page_config_calls[-1]["page_title"],
                    f"Optical pumping: {page_title}",
                )
                self.assertTrue(fake_streamlit.active_navigation.was_run)

    def test_navigation_groups_and_acc_page_route(self):
        entrypoint = Path(__file__).resolve().parents[1] / "alkali_pumping.py"
        fake_streamlit = _Streamlit("ACC Magnetometer")
        with (
            patch.dict("sys.modules", {"streamlit": fake_streamlit}),
            patch.object(
                page_state.st,
                "session_state",
                fake_streamlit.session_state,
            ),
        ):
            runpy.run_path(str(entrypoint), run_name="__navigation_test__")

        self.assertEqual(
            list(fake_streamlit.navigation_groups),
            ["Analysis", "Reference"],
        )
        self.assertEqual(
            [page.title for page in fake_streamlit.navigation_groups["Analysis"]],
            ["ACC Magnetometer"],
        )
        self.assertEqual(
            [page.title for page in fake_streamlit.navigation_groups["Reference"]],
            ["Light shift", "Atomic polarizability", "Atomic properties"],
        )
        acc_page = fake_streamlit.navigation_groups["Analysis"][0]
        self.assertEqual(
            acc_page.path,
            "alkali_pumping_app/pages/ACC-Magnetometer.py",
        )
        self.assertEqual(acc_page.url_path, "acc-magnetometer")
        self.assertTrue(acc_page.default)
        self.assertEqual(fake_streamlit.navigation_kwargs, {"position": "top"})


if __name__ == "__main__":
    unittest.main()
