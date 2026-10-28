"""Independent, versioned settings for coherent time evolution."""

from copy import deepcopy
from math import isfinite

from .conditions import DEFAULT_STARTUP_CONDITION, PUMP_CONDITION_KEYS, clean_condition_name
from .magnetic_field import validate_field_conditions
from .probe_presets import PROBE_PRESET_SOURCES
from ..physics.constants import ATOMS
from ..physics.polarization import allowed_polarizations
from ..physics.time_evolution import TIME_EVOLUTION_MODEL_VERSION


FORMAT = "alkali_pumping_time_evolution"
SCHEMA_VERSION = "2.7"
PREFIX = "te_"
PLOTTED_MOMENTS = (
    "Hyperfine polarization", "Electron polarization", "F_x", "F_y", "F_z",
    "Q_qq", "Q_uv", "Q_uu-Q_vv",
)
PLOTTED_PROBE_SIGNALS = ("Absorption", "Optical rotation", "Ellipticity")
PLOTTED_SCAN_SIGNALS = ("Optical rotation", "Ellipticity", "Absorption")
SHARED_KEYS = (
    "condition_name", "atom_A_name", "atom_B_name", "density_mode", "density_ratio_B_to_A",
    "gamma_ER_A", "gamma_ER_B", "static_field_theta_deg", "static_field_phi_deg", "static_field_nT",
    "gamma_u_A", "gamma_u_B",
    "temperature_C_for_table", "n2_pressure_torr",
    *[f"{line}_{field}_{label}" for label in ("A", "B")
      for line in ("D1", "D2") for field in ("width", "shift")],
    *PUMP_CONDITION_KEYS,
    *[f"probe_{field}_{label}" for label in ("A", "B") for field in (
        "source", "line", "transition", "det_rel", "k", "azimuth_deg", "ellipticity_deg", "path_length_cm")],
)
DEFAULTS = {key: deepcopy(DEFAULT_STARTUP_CONDITION[key]) for key in SHARED_KEYS}
DEFAULTS.update(condition_name="time-evolution", t1_s=1.0, t2_s=4.0,
                plotted_moments=["Hyperfine polarization", "F_z"])
for _label in ("A", "B"):
    _atom_name = DEFAULTS[f"atom_{_label}_name"]
    _atom = ATOMS[DEFAULTS["atom_A_name"] if _atom_name == "None" else _atom_name]
    _upper_F = float(_atom["I"]) + 0.5
    _lower_F = _upper_F - 1
    DEFAULTS.update({
        f"zeeman_axis_{_label}": "z", f"moment_axis_{_label}": "z",
        f"plotted_probe_signals_{_label}": ["Absorption", "Optical rotation"],
        f"plotted_scan_signals_{_label}": ["Optical rotation", "Absorption"],
        f"probe_scan_min_MHz_{_label}": -7000.0, f"probe_scan_max_MHz_{_label}": 7000.0,
        f"transition_{_label}1": f"{_lower_F:g}→{_upper_F:g}",
        f"det_rel_{_label}1": 0.0,
        f"pol_{_label}1": "sigma+", f"k_{_label}1": "z",
        f"transition_{_label}2": f"{_upper_F:g}→{_upper_F:g}",
        f"intensity_{_label}2": 0.0,
        f"pol_{_label}2": "sigma+", f"k_{_label}2": "z",
        f"probe_source_{_label}": "Custom", f"probe_k_{_label}": "z",
        f"probe_path_length_cm_{_label}": 1.0,
        f"probe_transition_{_label}": f"{_upper_F:g}→{_lower_F:g}",
        f"probe_det_rel_{_label}": -2000.0,
        f"probe_ellipticity_deg_{_label}": 0.0,
    })


def validate_settings(values):
    """Validate a complete current-format configuration without migrations."""
    if not isinstance(values, dict):
        raise ValueError("Settings must be a JSON object.")
    missing = [name for name in DEFAULTS if name not in values]
    if missing:
        raise ValueError("Missing Time evolution settings: " + ", ".join(missing))
    unknown = [name for name in values if name not in DEFAULTS]
    if unknown:
        raise ValueError("Unknown Time evolution settings: " + ", ".join(unknown))
    result = deepcopy(values)
    for key, default in DEFAULTS.items():
        value = result[key]
        if isinstance(default, bool):
            if not isinstance(value, bool):
                raise ValueError(f"{key} must be true or false.")
        elif isinstance(default, (float, int)):
            if isinstance(value, bool) or not isinstance(value, (float, int)) or not isfinite(value):
                raise ValueError(f"{key} must be a finite number.")
            result[key] = float(value)
        elif isinstance(default, str) and not isinstance(value, str):
            raise ValueError(f"{key} must be text.")
    choices = {
        "atom_A_name": list(ATOMS), "atom_B_name": ["None", *ATOMS],
        "density_mode": ["Independent saturated-vapor curves", "Relative concentration"],
    }
    validate_field_conditions(result)
    for field in ('static_field_theta_deg', 'static_field_phi_deg'):
        result[field] = int(result[field])
    for label in ("A", "B"):
        choices[f"zeeman_axis_{label}"] = ["x", "y", "z"]
        choices[f"moment_axis_{label}"] = ["x", "y", "z"]
        for number in (1, 2, 3):
            prefix = f"{label}{number}"
            choices[f"line_{prefix}"] = ["D1", "D2"]
            choices[f"k_{prefix}"] = ["x", "y", "z"]
        choices[f"probe_line_{label}"] = ["D1", "D2"]
        choices[f"probe_k_{label}"] = ["x", "y", "z"]
        choices[f"probe_source_{label}"] = ["Custom", *PROBE_PRESET_SOURCES, *[f"Pump{label}{n} weak" for n in (1, 2, 3)]]
    for key, allowed in choices.items():
        if result[key] not in allowed:
            raise ValueError(f"Invalid {key}.")
    for label in ("A", "B"):
        for number in (1, 2, 3):
            prefix = f"{label}{number}"
            if result[f"pol_{prefix}"] not in allowed_polarizations(result[f"k_{prefix}"]):
                raise ValueError(f"Invalid pump polarization for {prefix}.")
    for key, value in result.items():
        if (key.startswith(("intensity_", "gamma_ER_", "gamma_u_", "probe_path_length"))
                or key in ("n2_pressure_torr", "density_ratio_B_to_A")
                or "_width_" in key):
            if value < 0:
                raise ValueError(f"{key} must be nonnegative.")
    if not -100 <= result["temperature_C_for_table"] <= 500:
        raise ValueError("Temperature must be between -100 and 500 °C.")
    for label in ("A", "B"):
        if not -90 <= result[f"probe_azimuth_deg_{label}"] <= 90:
            raise ValueError("Probe azimuth must be between -90 and 90 degrees.")
        if not -45 <= result[f"probe_ellipticity_deg_{label}"] <= 45:
            raise ValueError("Probe ellipticity must be between -45 and 45 degrees.")
        lower, upper = result[f"probe_scan_min_MHz_{label}"], result[f"probe_scan_max_MHz_{label}"]
        if not (-10000 <= lower <= 9990 and -9990 <= upper <= 10000 and upper - lower >= 10):
            raise ValueError("Probe scan limits must lie within ±10000 MHz and be at least 10 MHz apart.")
    if not 0.001 <= result["t2_s"] <= 30 or not 0 <= result["t1_s"] <= result["t2_s"]:
        raise ValueError("Times must satisfy 0 ≤ t1 ≤ t2 and 0.001 ≤ t2 ≤ 30 seconds.")
    for field, allowed in (
        ("plotted_moments", PLOTTED_MOMENTS),
        *[(f"plotted_probe_signals_{label}", PLOTTED_PROBE_SIGNALS) for label in ("A", "B")],
        *[(f"plotted_scan_signals_{label}", PLOTTED_SCAN_SIGNALS) for label in ("A", "B")],
    ):
        selection = result[field]
        if (not isinstance(selection, list)
                or any(item not in allowed for item in selection)
                or len(selection) != len(set(selection))):
            raise ValueError(f"Invalid {field} selection.")
    result["condition_name"] = clean_condition_name(result["condition_name"])
    return result


def build_settings_payload(values):
    return {"app": "alkali_pumping", "format": FORMAT, "version": SCHEMA_VERSION,
            "physics_model": TIME_EVOLUTION_MODEL_VERSION, "settings": validate_settings(values)}


def parse_settings_payload(payload):
    if not isinstance(payload, dict) or payload.get("app") != "alkali_pumping" or payload.get("format") != FORMAT:
        raise ValueError("This is not a Time evolution settings file.")
    version = payload.get("version")
    if version != SCHEMA_VERSION:
        raise ValueError(f"Unsupported Time evolution settings version; expected {SCHEMA_VERSION}.")
    if payload.get("physics_model") != TIME_EVOLUTION_MODEL_VERSION:
        raise ValueError("Unsupported or missing Time evolution physics model.")
    return validate_settings(payload.get("settings"))
