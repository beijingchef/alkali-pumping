"""Consistent native Streamlit rich-text treatments for UI labels."""

from decimal import Decimal


def inactive_aware_label(label, inactive):
    """Gray a widget label only while the corresponding input is inactive."""
    return f":gray[{label}]" if inactive else label


def accent_caption(text, color="blue"):
    """Return native Streamlit colored text suitable for ``st.caption``."""
    return f":{color}[{text}]"


def pump_rate_caption(rate):
    """Show three significant digits without exponent notation."""
    number = format(Decimal(f"{float(rate):.3g}"), "f")
    if "." in number:
        number = number.rstrip("0").rstrip(".")
    return rf"$R_{{\mathrm{{OP}}}}$: {number} s$^{{-1}}$"


def line_center_detuning_caption(line, detuning_MHz, *, pump=False):
    """Format the physical detuning from a zero-pressure D-line center."""
    if pump:
        return accent_caption(
            rf"$\delta\nu_{{\mathrm{{{line}}}}}$: {float(detuning_MHz):+.0f} MHz"
        )
    return accent_caption(
        f"δν from {line} center: {float(detuning_MHz):+.0f} MHz"
    )
