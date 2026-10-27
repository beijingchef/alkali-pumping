"""Browser-only settings downloads; Save never submits widget state."""
import base64
from hashlib import sha256
from html import escape
import json
import streamlit as st


def _button_id(key):
    return "settings-save-" + sha256(key.encode("utf-8")).hexdigest()[:16]


def _encode(data):
    raw = data.encode("utf-8") if isinstance(data, str) else bytes(data)
    return base64.b64encode(raw).decode("ascii")


def update_settings_download_snapshot(key, data):
    """Refresh committed settings after an ACC pump-only fragment rerun."""
    button_id = _button_id(key)
    marker = button_id + "-update"
    st.html(f"""
    <span id="{marker}" hidden></span>
    <style>[data-testid="stElementContainer"]:has(#{marker}) {{ display: none; }}</style>
    <script>
    (() => {{
      const button = document.getElementById({json.dumps(button_id)});
      if (button) button.dataset.settingsSnapshot = {json.dumps(_encode(data))};
    }})();
    </script>
    """, unsafe_allow_javascript=True)


def save_button_with_immediate_download(
    container, *, data=None, file_name=None, prepare=None,
    mime="application/json", key, label="Save", width="stretch", fallback_name="default",
):
    """Publish committed settings, then download entirely in the browser.

    Only the live Setting name is read on Save. Preventing mouse-down focus
    changes leaves other dirty inputs uncommitted. No Streamlit widget event,
    callback, fragment rerun or calculation is requested by this button.
    """
    if prepare is not None:
        data, file_name = prepare()
    button_id = _button_id(key)
    body = f"""
    <style>
    #{button_id} {{
      display: inline-flex; align-items: center; justify-content: center;
      width: 100%; min-height: 2.5rem; padding: .25rem .75rem;
      border: 1px solid var(--st-border-color, rgba(49,51,63,.2));
      border-radius: var(--st-base-radius, .5rem);
      color: var(--st-text-color, #31333f);
      background: var(--st-background-color, white);
      font: inherit; font-size: .875rem; cursor: pointer;
    }}
    #{button_id}:hover {{ border-color: var(--st-primary-color, #ff4b4b); color: var(--st-primary-color, #ff4b4b); }}
    #{button_id}:focus-visible {{ outline: 2px solid var(--st-primary-color, #ff4b4b); outline-offset: 2px; }}
    </style>
    <button type="button" id="{button_id}" class="alkali-settings-save"
      title="Save committed settings using the current setting name"
      data-settings-snapshot="{_encode(data)}">{escape(label)}</button>
    <script>
    (() => {{
      const button = document.getElementById({json.dumps(button_id)});
      if (!button || button.dataset.downloadBound) return;
      button.dataset.downloadBound = "true";
      // Do not blur an unfinished input: Streamlit commits inputs on blur.
      for (const eventName of ["pointerdown", "mousedown"]) {{
        button.addEventListener(eventName, event => event.preventDefault());
      }}
      button.addEventListener("click", event => {{
        event.preventDefault();
        event.stopPropagation();
        const bytes = Uint8Array.from(atob(button.dataset.settingsSnapshot), c => c.charCodeAt(0));
        const payload = JSON.parse(new TextDecoder().decode(bytes));
        const settings = payload.conditions || payload.settings;
        const input = document.querySelector('[data-testid="stSidebar"] input[aria-label="Setting name"]');
        let name = (input ? input.value : settings.condition_name).trim().replace(/\\.json$/i, "").trim();
        name = name || {json.dumps(fallback_name)};
        settings.condition_name = name;
        if ("saved_at" in payload) {{
          const now = new Date();
          payload.saved_at = new Date(now.getTime() - now.getTimezoneOffset()*60000).toISOString().slice(0,19);
        }}
        const blob = new Blob([JSON.stringify(payload, null, 2)], {{type: {json.dumps(mime)}}});
        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");
        link.href = url;
        link.download = name + ".json";
        document.body.appendChild(link);
        link.click();
        link.remove();
        setTimeout(() => URL.revokeObjectURL(url), 1000);
      }});
    }})();
    </script>
    """
    with container.container():
        st.html(body, unsafe_allow_javascript=True, width=width)
