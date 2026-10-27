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
                "alkali_pumping_app/pages/ACC-Magnetometer.py",
                title="ACC Magnetometer",
                icon=":material/explore:",
                url_path="acc-magnetometer",
                default=True,
            ),
            st.Page(
                "alkali_pumping_app/pages/time_evolution.py",
                title="Time evolution",
                icon=":material/timeline:",
                url_path="time-evolution",
            ),
        ],
        "Reference": [
            st.Page(
                "alkali_pumping_app/pages/light_shift.py",
                title="Light shift",
                icon=":material/waves:",
                url_path="light-shift",
            ),
            st.Page(
                "alkali_pumping_app/pages/atomic_polarizability.py",
                title="Atomic polarizability",
                icon=":material/ssid_chart:",
                url_path="atomic-polarizability",
            ),
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

# Keep the 3.75rem header clear and halve the remaining native content gap:
# (6rem + 3.75rem)/2 without top navigation, (8rem + 3.75rem)/2 with it.
# The native bottom container reserves 3.5rem below its controls; halve that too.
st.html("""
<style>
@media screen {
    [data-testid="stMainBlockContainer"] {
        padding-top: 4.875rem;
    }
    [data-testid="stAppViewContainer"]:has([data-testid="stTopNavSection"])
    [data-testid="stMainBlockContainer"] {
        padding-top: 5.875rem;
    }
    [data-testid="stBottomBlockContainer"] {
        padding-bottom: 1.75rem;
    }
    [data-testid="stSidebarHeader"] {
        height: 1.875rem;
        margin-bottom: 0.5rem;
    }
}
</style>
""")
navigation.run()
