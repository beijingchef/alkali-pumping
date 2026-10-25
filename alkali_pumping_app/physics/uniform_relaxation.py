"""Uniform ground-state relaxation toward the completely unpolarized state.

D_u[rho] = R_u (Tr(rho) I/N - rho). This trace-preserving depolarizing
channel redistributes populations over all ground Zeeman states and damps
every off-diagonal coherence at R_u, with no transfer between coherences.
Unlike electron randomization, it also erases nuclear polarization.
"""

import numpy as np


def uniform_relaxation_generator(state_count, rate_s):
    """Return the population generator, with columns summing to zero."""
    rate_s = float(rate_s)
    if not np.isfinite(rate_s) or rate_s < 0:
        raise ValueError("Uniform relaxation rate must be finite and nonnegative.")
    if state_count < 1:
        raise ValueError("Uniform relaxation needs at least one ground state.")
    return rate_s * (np.full((state_count, state_count), 1 / state_count) - np.eye(state_count))


def uniform_relaxation_diagnostics(ground_states, population, rate_s):
    """Signed fractional population loss and adjacent-coherence decay rates."""
    population = np.asarray(population, dtype=float)
    derivative = uniform_relaxation_generator(len(ground_states), rate_s) @ population
    fractional = np.divide(-derivative, population, out=np.full_like(population, np.nan),
                           where=population > 1e-15)
    states = {(float(state["F"]), float(state["m"])) for state in ground_states}
    coherence = np.array([
        float(rate_s) if (float(state["F"]), float(state["m"]) - 1) in states else np.nan
        for state in ground_states
    ])
    return fractional, coherence
