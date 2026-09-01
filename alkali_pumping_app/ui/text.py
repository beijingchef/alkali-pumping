"""Consistent native Streamlit rich-text treatments for UI labels."""


def inactive_aware_label(label, inactive):
    """Gray a widget label only while the corresponding input is inactive."""
    return f":gray[{label}]" if inactive else label


def accent_caption(text, color="blue"):
    """Return native Streamlit colored text suitable for ``st.caption``."""
    return f":{color}[{text}]"


def line_center_detuning_caption(line, detuning_MHz):
    """Format the physical detuning from a zero-pressure D-line center."""
    return accent_caption(
        f"δν from {line} center: {float(detuning_MHz):+.0f} MHz"
    )
