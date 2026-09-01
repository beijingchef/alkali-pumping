"""Compact native file-upload controls."""

import streamlit as st


_COMPACT_UPLOADER_STYLE = """
<style>
[data-testid="stFileUploaderDropzone"] {
  min-height: 0;
  padding: 0;
  border: 0;
  background: transparent;
}
[data-testid="stFileUploaderDropzoneInstructions"],
[data-testid="stFileUploaderFile"] {
  display: none;
}
[data-testid="stFileUploaderDropzone"] button {
  width: 100%;
  margin: 0;
  gap: 0;
}
[data-testid="stFileUploaderDropzone"] button svg,
[data-testid="stFileUploaderDropzone"] button [data-testid="stIconMaterial"],
[data-testid="stFileUploaderDropzone"] button span[translate="no"] {
  display: none !important;
}
[data-testid="stFileUploaderDropzone"] button p {
  font-size: 0;
}
[data-testid="stFileUploaderDropzone"] button p::after {
  content: "Open";
  font-size: 0.875rem;
}
</style>
"""


def open_file_button(*, type, key, on_change, help=None):
    """Render the native uploader as a compact Open button."""
    st.html(_COMPACT_UPLOADER_STYLE)
    return st.file_uploader(
        "Open",
        type=type,
        key=key,
        on_change=on_change,
        label_visibility="collapsed",
        help=help,
        width="stretch",
    )
