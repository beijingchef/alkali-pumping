"""Independent Zeeman display bases for stationary probe-scan populations."""
import numpy as np
from scipy.linalg import expm

from .coherent_evolution import basis_rotation, spin_matrices
from .field_geometry import static_field_direction
from .polarization import transverse_basis_for_k

SCAN_POPULATION_AXES = ("x", "y", "z", "B direction", "probe major", "probe minor")


def scan_axis_direction(axis, common, probe):
    """Resolve a directed display axis in laboratory Cartesian coordinates.

    A zero field uses its configured spherical direction. Ellipse axes use
    the input azimuth and a right-handed major/minor/k frame, independent of
    helicity. For circular input these are the azimuth reference axes.
    """
    if axis in ("x", "y", "z"):
        return np.eye(3)["xyz".index(axis)]
    if axis == "B direction":
        return (-1. if common["static_field_nT"] < 0 else 1.) * static_field_direction(common)
    if axis in ("probe major", "probe minor"):
        a, b = (np.asarray(v.real, dtype=float) for v in transverse_basis_for_k(probe["k_axis"]))
        angle = np.deg2rad(float(probe["azimuth_deg"]))
        major = np.cos(angle)*a + np.sin(angle)*b
        minor = -np.sin(angle)*a + np.cos(angle)*b
        return major if axis == "probe major" else minor
    raise ValueError(f"Unknown scan population axis: {axis}")


def scan_population_view(species, axis, direction):
    """Rotate the full steady matrix into |F,m> along direction; no dynamics.

    Columns of U=exp(-i phi Fz) exp(-i theta Fy) are the exported basis
    vectors in lab-z. This Euler convention fixes coherence phases as well
    as populations. Source Cartesian bases retain their existing convention.
    """
    direction = np.asarray(direction, dtype=float)
    if direction.shape != (3,) or not np.isfinite(direction).all() or np.linalg.norm(direction) == 0:
        raise ValueError("Scan quantization direction must be a finite nonzero 3-vector.")
    direction = direction / np.linalg.norm(direction)
    theta = np.arccos(np.clip(direction[2], -1., 1.))
    phi = np.arctan2(direction[1], direction[0]) if np.linalg.norm(direction[:2]) > 1e-14 else 0.
    states = species["ground_states"]
    _, Fy, Fz = spin_matrices(states)
    rotation = expm(-1j * phi * Fz) @ expm(-1j * theta * Fy)
    transform = rotation.conj().T @ basis_rotation(states, species["q_axis"])
    density = transform @ species["density_matrix"] @ transform.conj().T
    return dict(species, q_axis=axis, zeeman_axis=axis, density_matrix=density,
                population=np.diagonal(density, axis1=1, axis2=2).real.copy(),
                quantization_axis_vector=direction,
                basis_theta_deg=float(np.rad2deg(theta)), basis_phi_deg=float(np.rad2deg(phi)),
                basis_rotation_convention="U=exp(-i phi F_z) exp(-i theta F_y); dimensionless F")
