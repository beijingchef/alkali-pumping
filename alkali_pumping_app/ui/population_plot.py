"""Shared Zeeman population bars for stationary states."""
import matplotlib.pyplot as plt
import altair as alt
import numpy as np
import streamlit as st
from .plot_data import aligned_population_bar_data


def render_population_plot(result, *, figure_width=4.6, height=None, bar_width=.72):
    states, population = result['ground_states'], np.asarray(result['population'])
    common_m_values, groups = aligned_population_bar_data(states, population)
    if height is not None:
        # Native chart viewports keep the combined height fixed while the
        # 30% population column changes width.
        limit = max(.01, 1.08*float(np.max(population)))
        # Separate native chart viewports include their axes in the height;
        # concatenated Vega views otherwise add axes outside the stated height.
        with st.container(gap=None):
            for index, (F, m_values, manifold_population) in enumerate(groups):
                last = index == len(groups)-1
                chart = alt.Chart({'values': [dict(m=float(m), left=float(m)-bar_width/2, right=float(m)+bar_width/2,
                                                  population=float(p))
                                             for m, p in zip(m_values, manifold_population)]}).mark_bar(
                    color='#1f77b4', orient='vertical',
                ).encode(
                    x=alt.X('left:Q', scale=alt.Scale(domain=[common_m_values[0]-max(.6, bar_width/2),
                                                             common_m_values[-1]+max(.6, bar_width/2)]),
                            axis=alt.Axis(values=list(common_m_values), labels=last, ticks=last, grid=False,
                                          title=f"m along {result['q_axis']}" if last else None)),
                    x2=alt.X2('right:Q'),
                    y=alt.Y('population:Q', title=f'F={F:g}', scale=alt.Scale(domain=[0, limit]),
                            axis=alt.Axis(tickCount=3, grid=False)),
                    y2=alt.datum(0),
                    tooltip=[alt.Tooltip('m:Q'), alt.Tooltip('population:Q', format='.6g')],
                ).properties(height=height/len(groups), description=f'Steady-state populations F={F:g}')
                st.altair_chart(chart, width='stretch', height=int(height/len(groups)))
        return
    fig, axes = plt.subplots(len(groups), 1, figsize=(figure_width, 4.5), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for axis, (F, m_values, manifold_population) in zip(axes, groups):
        axis.bar(m_values, manifold_population, width=bar_width)
        axis.set_ylabel(f'F={F:g}')
        axis.grid(False)
    axes[-1].set_xticks(common_m_values)
    axes[-1].set_xticklabels([f'{m:g}' for m in common_m_values])
    axes[-1].set_xlabel(f"m along {result['q_axis']}")
    if len(common_m_values):
        axes[-1].set_xlim(common_m_values[0]-.6, common_m_values[-1]+.6)
    for axis in axes:
        axis.set_ylim(0., max(.01, 1.08*float(np.max(population))))
    fig.tight_layout()
    st.pyplot(fig, width='stretch')
    plt.close(fig)
