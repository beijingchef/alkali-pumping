"""Top-level navigation for the Alkali Pumping application."""

import streamlit as st

from alkali_pumping_app.ui.page_state import preserve_persistent_page_settings

st.set_page_config(
    layout="wide",
    initial_sidebar_state="expanded",
)

preserve_persistent_page_settings()

navigation = st.navigation(
    {
        "Analysis": [
            st.Page(
                "alkali_pumping_app/pages/light_shift.py",
                title="Light shift",
                icon=":material/waves:",
                url_path="light-shift",
                default=True,
            ),
            st.Page(
                "alkali_pumping_app/pages/atomic_polarizability.py",
                title="Atomic polarizability",
                icon=":material/ssid_chart:",
                url_path="atomic-polarizability",
            ),
            st.Page(
                "alkali_pumping_app/pages/magnetometry.py",
                title="Magnetometry",
                icon=":material/explore:",
                url_path="magnetometry",
            ),
        ],
        "Reference": [
            st.Page(
                "alkali_pumping_app/pages/atomic_properties.py",
                title="Atomic properties",
                icon=":material/science:",
                url_path="atomic-properties",
            )
        ],
    },
    position="top",
)

st.set_page_config(page_title=f"Optical pumping: {navigation.title}")
navigation.run()
