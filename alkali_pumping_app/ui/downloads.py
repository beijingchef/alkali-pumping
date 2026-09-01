"""Download controls whose files depend on newly committed widget values."""

import base64
import json

import streamlit as st


def save_button_with_immediate_download(
    container,
    *,
    data,
    file_name,
    mime="application/octet-stream",
    key,
    label="Save",
    width="stretch",
):
    """Commit pending widgets, then download the freshly constructed file.

    Unlike ``st.download_button``, a regular button submits an uncommitted text
    input before Python constructs the file. The short browser-side script then
    starts the download produced by that same rerun.
    """
    if not container.button(label, key=key, width=width):
        return False

    raw_data = data.encode("utf-8") if isinstance(data, str) else bytes(data)
    encoded_data = base64.b64encode(raw_data).decode("ascii")
    script = f"""
    <script>
    (() => {{
      const encoded = {json.dumps(encoded_data)};
      const binary = atob(encoded);
      const bytes = new Uint8Array(binary.length);
      for (let index = 0; index < binary.length; index += 1) {{
        bytes[index] = binary.charCodeAt(index);
      }}
      const blob = new Blob([bytes], {{type: {json.dumps(mime)}}});
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = {json.dumps(file_name)};
      document.body.appendChild(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    }})();
    </script>
    """
    st.html(script, unsafe_allow_javascript=True)
    return True
