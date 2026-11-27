"""Angular momentum coherent-state probability on a sphere."""

from math import comb, pi

import numpy as np


def angular_probability_surface(density_matrix, ground_states, F, theta, phi):
    """Return Q_F(n)=(2F+1)/(4pi) <F,n|rho|F,n> per steradian.

    The density block is not renormalized: integrating Q_F over solid angle
    yields the population of manifold F, not unity. The coherent ket uses the
    local quantization axis as z and phases exp(-i m phi).
    """
    indices = [i for i, state in enumerate(ground_states)
               if np.isclose(float(state["F"]), float(F))]
    if not indices:
        raise ValueError(f"No ground states for F={F:g}.")
    theta, phi = np.broadcast_arrays(theta, phi)
    spin_twice = int(round(2.0 * float(F)))
    amplitudes = np.stack([
        np.sqrt(comb(spin_twice, int(round(float(F) + ground_states[i]["m"]))))
        * np.cos(theta / 2.0) ** (float(F) + ground_states[i]["m"])
        * np.sin(theta / 2.0) ** (float(F) - ground_states[i]["m"])
        * np.exp(-1j * float(ground_states[i]["m"]) * phi)
        for i in indices
    ], axis=-1)
    block = np.asarray(density_matrix, dtype=complex)[np.ix_(indices, indices)]
    probability = np.einsum("...i,ij,...j->...", amplitudes.conj(), block, amplitudes)
    return (2.0 * float(F) + 1.0) / (4.0 * pi) * np.real(probability)
