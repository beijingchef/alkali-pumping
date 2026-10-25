"""Isotope-relative optical settings for nonperturbing probe presets."""


PROBE_PRESET_SOURCES = ("Preset D1a", "Preset D1b", "Preset D2a", "Preset D2b")


def probe_preset_settings(source, atom):
    """Return probe field values, or no overrides for a non-preset source.

    Both D lines use F = I ± 1/2 and F′ = I ∓ 1/2. In particular, the
    D2 presets do not use the outermost excited-state hyperfine manifolds.
    Detuning is relative to the chosen pressure-shifted transition, in MHz.
    """
    if source not in PROBE_PRESET_SOURCES:
        return {}
    branch = source[-1]
    upper = float(atom["I"]) + 0.5
    lower = float(atom["I"]) - 0.5
    ground, excited = (upper, lower) if branch == "a" else (lower, upper)
    return {
        "line": source[-3:-1],
        "transition": f"{ground:g}→{excited:g}",
        "det_rel": -2000.0 if branch == "a" else 2000.0,
        "k": "z",
        "path_length_cm": 1.0,
        "azimuth_deg": 0.0,
        "ellipticity_deg": 0.0,
    }
