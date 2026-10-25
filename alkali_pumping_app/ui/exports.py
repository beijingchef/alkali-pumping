"""CSV export helpers for calculated app results."""

import numpy as np
import pandas as pd

from ..version import PHYSICS_MODEL_VERSION


def dataframe_to_csv_bytes(dataframe):
    """Serialize a dataframe as an Excel-friendly UTF-8 CSV download."""
    return dataframe.to_csv(index=False).encode("utf-8-sig")


def atomic_polarizability_export_dataframe(plotted):
    """Return one column per polarizability curve with a units row."""
    component_labels = {
        "alpha_eq": "α_eq",
        "alpha_hfs": "α_hfs",
        "alpha_gt": "α_gt",
        "alpha_br": "α_br",
    }
    detuning_column = "Laser detuning"
    detunings = np.sort(plotted["Detuning (MHz)"].unique())
    numeric = pd.DataFrame({detuning_column: detunings})
    units = {detuning_column: "MHz"}

    for (component, series, part), group in plotted.groupby(
        ["Component", "Series", "Part"], sort=False
    ):
        symbol = component_labels[component]
        label = (
            f"{symbol} | {part}"
            if component in ("alpha_eq", "alpha_hfs")
            else f"{symbol} | {series} | {part}"
        )
        values = group.set_index("Detuning (MHz)")[
            "Polarizability (10^-18 cm^3)"
        ].reindex(detunings)
        numeric[label] = values.to_numpy(dtype=float)
        units[label] = "10⁻¹⁸ cm³"

    return pd.concat([pd.DataFrame([units]), numeric], ignore_index=True)


def weak_rf_export_dataframe(
    frequencies_hz,
    susceptibility_amplitude,
    susceptibility_in_phase,
    susceptibility_quadrature,
    plotted_amplitude,
    plotted_in_phase,
    plotted_quadrature,
    *,
    in_phase_plot_factor=1.0,
    quadrature_plot_factor=1.0,
    relaxation_normalized,
    normalization_gamma_s_inv=None,
    density_factored=False,
    density_cm3=None,
    observable="Fx",
):
    """Return raw and plotted weak-RF susceptibility samples for export.

    Raw susceptibility values retain the calculation's phase convention. The
    plotted signed components include their independently selected phase
    factors, relaxation normalization, and the alkali-density factor when
    those display options are active.
    """
    spin_unit = "hbar2" if observable == "Q_ij" else "hbar"
    arrays = {
        "frequency_Hz": np.asarray(frequencies_hz, dtype=float),
        f"amplitude_raw_{spin_unit}_s_per_atom": np.asarray(
            susceptibility_amplitude, dtype=float
        ),
        f"in_phase_raw_{spin_unit}_s_per_atom": np.asarray(
            susceptibility_in_phase, dtype=float
        ),
        f"quadrature_raw_{spin_unit}_s_per_atom": np.asarray(
            susceptibility_quadrature, dtype=float
        ),
        "amplitude_plotted": np.asarray(plotted_amplitude, dtype=float),
        "in_phase_plotted": np.asarray(plotted_in_phase, dtype=float),
        "quadrature_plotted": np.asarray(plotted_quadrature, dtype=float),
    }
    sample_counts = {len(values) for values in arrays.values()}
    if len(sample_counts) != 1:
        raise ValueError("All weak-RF export arrays must have the same length.")

    if density_factored:
        plotted_units = "hbar/cm^3" if relaxation_normalized else "hbar s/cm^3"
    else:
        plotted_units = "hbar/atom" if relaxation_normalized else "hbar s/atom"
    normalization_gamma = (
        float(normalization_gamma_s_inv)
        if relaxation_normalized and normalization_gamma_s_inv is not None
        else np.nan
    )
    density_value = (
        float(density_cm3)
        if density_factored and density_cm3 is not None
        else np.nan
    )

    dataframe = pd.DataFrame(arrays)
    dataframe["plotted_units"] = plotted_units
    if spin_unit == "hbar2":
        dataframe["plotted_units"] = plotted_units.replace("hbar", "hbar^2")
    dataframe["in_phase_plot_factor"] = float(in_phase_plot_factor)
    dataframe["quadrature_plot_factor"] = float(quadrature_plot_factor)
    dataframe["relaxation_normalized"] = bool(relaxation_normalized)
    dataframe["normalization_gamma_s_inv"] = normalization_gamma
    dataframe["density_factored"] = bool(density_factored)
    dataframe["density_cm3"] = density_value
    return dataframe


def weak_probe_export_dataframe(
    frequencies_hz, response, signal, *, rf_rabi_rad_s_per_nT,
    manifold_response=None, manifold_F=None, displayed_manifold="Total",
    in_phase_plot_factor=1.0, quadrature_plot_factor=1.0,
    drive_info=None,
    uniform_relaxation_rate_s=0.0,
):
    """Return total and rank-resolved weak-probe readouts per RF field."""
    data = {"frequency_Hz": np.asarray(frequencies_hz, dtype=float)}
    for component in ("total", "scalar", "orientation", "alignment"):
        values = response[component][signal]
        for field in ("amplitude", "in_phase", "quadrature"):
            data[f"{component}_{field}_per_nT"] = (
                np.asarray(values[field], dtype=float)
                * float(rf_rabi_rad_s_per_nT)
            )
    frame = pd.DataFrame(data)
    frame["signal"] = signal
    frame["units"] = (
        "rad/nT" if signal in ("rotation", "ellipticity") else "1/nT"
    )
    frame["rf_rabi_rad_s_per_nT"] = float(rf_rabi_rad_s_per_nT)
    for manifold, parts in (manifold_response or {}).items():
        for component, signals in parts.items():
            for field in ("amplitude", "in_phase", "quadrature"):
                frame[f"{manifold}_{component}_{field}_per_nT"] = (
                    np.asarray(signals[signal][field], dtype=float) * rf_rabi_rad_s_per_nT
                )
        if manifold_F:
            frame[f"{manifold}_F"] = manifold_F[manifold]
    frame["displayed_manifold"] = displayed_manifold
    frame["in_phase_plot_factor"] = float(in_phase_plot_factor)
    frame["quadrature_plot_factor"] = float(quadrature_plot_factor)
    frame["rf_normalization"] = "Omega_a = g_a mu_B B_rf / hbar"
    frame["physics_model"] = PHYSICS_MODEL_VERSION
    frame["R_u_s-1"] = float(uniform_relaxation_rate_s)
    for manifold, ratio in (drive_info or {}).get("drive_ratios", {}).items():
        frame[f"{manifold}_rf_drive_ratio_gF_over_ga"] = ratio
        frame[f"{manifold}_g_F"] = drive_info["g_factors"][manifold]
    return frame


def manifold_rf_export_dataframe(result, plotted, selected_manifolds, **display_options):
    """Export both manifolds on the same field/linewidth/population scale."""
    frames = []
    for manifold, response in result["rf_response"].items():
        frame = weak_rf_export_dataframe(
            result["rf_frequencies_hz"], response["amplitude"],
            response["in_phase"], response["quadrature"], *plotted[manifold],
            observable=result["rf_observable"], **display_options,
        )
        frame["manifold"] = manifold
        frame["F"] = result["rf_manifolds"][manifold]
        frame["observable"] = result["rf_observable"]
        frame["selected_for_plot"] = manifold in selected_manifolds
        frame["g_F"] = result["rf_drive_info"]["g_factors"][manifold]
        frame["rf_drive_ratio_gF_over_ga"] = result["rf_drive_info"]["drive_ratios"][manifold]
        frame["rf_normalization"] = result["rf_drive_info"]["normalization"]
        frame["physics_model"] = PHYSICS_MODEL_VERSION
        frame["R_u_s-1"] = float(result.get("R_u", 0.0))
        reference = result["rf_relaxation_reference"]
        frame["normalization_reference_F"] = reference.get("F", np.nan)
        frame["normalization_reference_m"] = reference.get("m", np.nan)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)
