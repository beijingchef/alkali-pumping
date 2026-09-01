"""Pure data-layout helpers for application plots."""

import numpy as np


def aligned_population_bar_data(states, population):
    """Return hyperfine population groups on one common numerical m axis."""
    population = np.asarray(population, dtype=float)
    if len(states) != len(population):
        raise ValueError("Population length must match the number of states.")

    energies = {}
    for state in states:
        energies.setdefault(float(state["F"]), float(state["E"]))
    manifolds = sorted(energies, key=energies.get, reverse=True)
    common_m_values = np.asarray(
        sorted({float(state["m"]) for state in states}),
        dtype=float,
    )

    groups = []
    for F in manifolds:
        indices = sorted(
            (
                index
                for index, state in enumerate(states)
                if np.isclose(float(state["F"]), F)
            ),
            key=lambda index: float(states[index]["m"]),
        )
        groups.append(
            (
                F,
                np.asarray([float(states[index]["m"]) for index in indices]),
                population[indices],
            )
        )
    return common_m_values, groups
