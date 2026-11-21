"""Shared chart text styling for application plots."""

from matplotlib.legend import Legend


PLOT_TEXT_COLOR = "#000000"


def style_altair_text(chart):
    """Return an Altair chart with black axes, tick marks, and legends."""
    return chart.configure_axis(
        labelColor=PLOT_TEXT_COLOR,
        titleColor=PLOT_TEXT_COLOR,
        tickColor=PLOT_TEXT_COLOR,
        domainColor=PLOT_TEXT_COLOR,
    ).configure_legend(
        labelColor=PLOT_TEXT_COLOR,
        titleColor=PLOT_TEXT_COLOR,
    )


def style_matplotlib_axis_text(axis):
    """Set axis labels, ticks, domains, and every attached legend to black."""
    axis.xaxis.label.set_color(PLOT_TEXT_COLOR)
    axis.yaxis.label.set_color(PLOT_TEXT_COLOR)
    axis.xaxis.offsetText.set_color(PLOT_TEXT_COLOR)
    axis.yaxis.offsetText.set_color(PLOT_TEXT_COLOR)
    axis.tick_params(axis="both", colors=PLOT_TEXT_COLOR)
    for spine in axis.spines.values():
        spine.set_color(PLOT_TEXT_COLOR)
    for legend in axis.findobj(match=Legend):
        legend.get_title().set_color(PLOT_TEXT_COLOR)
        for text in legend.get_texts():
            text.set_color(PLOT_TEXT_COLOR)
    return axis
