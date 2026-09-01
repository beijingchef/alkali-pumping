"""Display-only transformations for weak-RF susceptibility curves."""

from matplotlib.lines import Line2D
import numpy as np


def rf_component_legend_label(component, add_pi=False):
    """Return X/Y or the π-shifted −X/−Y legend label."""
    component = str(component)
    if component not in ("X", "Y"):
        raise ValueError(f"Unsupported RF component: {component}")
    return f"−{component}" if add_pi else component


def prepare_weak_rf_plot_values(
    amplitude,
    in_phase,
    quadrature,
    *,
    flip_in_phase=False,
    flip_quadrature=False,
    relaxation_gamma_s_inv=None,
    density_cm3=None,
):
    """Apply independent optional π shifts and scalar display factors."""
    factor = 1.0
    if relaxation_gamma_s_inv is not None:
        factor *= float(relaxation_gamma_s_inv)
    if density_cm3 is not None:
        factor *= float(density_cm3)

    plotted_amplitude = np.asarray(amplitude, dtype=float).copy() * factor
    in_phase_factor = -1.0 if flip_in_phase else 1.0
    quadrature_factor = -1.0 if flip_quadrature else 1.0
    plotted_in_phase = (
        in_phase_factor * np.asarray(in_phase, dtype=float).copy() * factor
    )
    plotted_quadrature = (
        quadrature_factor * np.asarray(quadrature, dtype=float).copy() * factor
    )
    return plotted_amplitude, plotted_in_phase, plotted_quadrature


def add_probe_decomposition_legend(axis, entries_by_column):
    """Add a contribution-column/response-row legend to a probe plot.

    ``entries_by_column`` maps the ordered headers Total, Orientation, and
    Alignment to equal-length lists of ``(line_handle, row_label)`` pairs.
    Matplotlib fills multi-column legends down each column, so grouping each
    header with its row entries produces the requested table layout.
    """
    headers = ("Total", "Orientation", "Alignment")
    missing = [header for header in headers if header not in entries_by_column]
    if missing:
        raise ValueError("Missing probe legend columns: " + ", ".join(missing))
    row_counts = {len(entries_by_column[header]) for header in headers}
    if len(row_counts) != 1:
        raise ValueError("Probe legend columns must contain the same rows.")

    handles = []
    labels = []
    header_indices = []
    for header in headers:
        header_indices.append(len(labels))
        handles.append(Line2D([], [], linestyle="none"))
        labels.append(header)
        for line_handle, row_label in entries_by_column[header]:
            handles.append(line_handle)
            labels.append(row_label)

    legend = axis.legend(
        handles,
        labels,
        ncol=len(headers),
        frameon=False,
        columnspacing=1.6,
        handletextpad=0.6,
    )
    texts = legend.get_texts()
    for index in header_indices:
        texts[index].set_weight("bold")
    return legend
