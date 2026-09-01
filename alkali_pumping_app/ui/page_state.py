"""Keep registered page-setting widget values across Streamlit navigation."""

import streamlit as st


_REGISTRY_KEY = "_persistent_page_setting_keys"


def register_persistent_page_settings(keys):
    """Register setting keys that should survive while their page is hidden."""
    registered = set(st.session_state.get(_REGISTRY_KEY, ()))
    registered.update(str(key) for key in keys)
    st.session_state[_REGISTRY_KEY] = tuple(sorted(registered))


def preserve_persistent_page_settings():
    """Detach registered settings from Streamlit's unrendered-widget cleanup."""
    for key in st.session_state.get(_REGISTRY_KEY, ()):
        if key in st.session_state:
            st.session_state[key] = st.session_state[key]
