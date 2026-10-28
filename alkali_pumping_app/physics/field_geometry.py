"""Static-field geometry in the common laboratory Cartesian frame."""
import numpy as np


def static_field_direction(common):
    """Unit direction: theta from +z, phi from +x toward +y, in degrees."""
    theta = float(common['static_field_theta_deg'])
    phi = float(common['static_field_phi_deg'])
    if not np.isfinite([theta, phi]).all() or not 0 <= theta <= 180 or not 0 <= phi <= 360:
        raise ValueError('Field angles must be finite, with 0 ≤ theta ≤ 180 and 0 ≤ phi ≤ 360 degrees.')
    theta, phi = np.deg2rad([theta, phi])
    direction = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    # Cardinal directions must retain exact zero transverse components.
    direction[np.abs(direction) < 1e-15] = 0.
    return direction


def static_field_vector_nT(common):
    strength = float(common['static_field_nT'])
    if not np.isfinite(strength):
        raise ValueError('Static field strength must be finite.')
    return strength * static_field_direction(common)
