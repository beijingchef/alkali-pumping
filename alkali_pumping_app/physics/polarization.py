"""Beam geometry and polarization-basis transformations."""

from math import cos, radians, sin, sqrt

import numpy as np

# ============================================================
# 4. Light geometry and polarization
# ============================================================

def unit_vector(axis):
    if axis == "x":
        return np.array([1.0, 0.0, 0.0], dtype=complex)
    if axis == "y":
        return np.array([0.0, 1.0, 0.0], dtype=complex)
    if axis == "z":
        return np.array([0.0, 0.0, 1.0], dtype=complex)
    raise ValueError(axis)


def transverse_basis_for_k(k_axis):
    """
    Return two transverse unit vectors a,b such that a x b = k.
    This defines the local beam frame.
    """
    if k_axis == "z":
        return unit_vector("x"), unit_vector("y")
    if k_axis == "x":
        return unit_vector("y"), unit_vector("z")
    if k_axis == "y":
        return unit_vector("z"), unit_vector("x")
    raise ValueError(k_axis)


def allowed_polarizations(k_axis):
    """
    For each propagation direction, allow sigma+, sigma-, and the two lab-linear
    directions perpendicular to k.
    """
    if k_axis == "z":
        return ["sigma+", "sigma-", "linear x", "linear y"]
    if k_axis == "x":
        return ["sigma+", "sigma-", "linear y", "linear z"]
    if k_axis == "y":
        return ["sigma+", "sigma-", "linear z", "linear x"]
    raise ValueError(k_axis)


def lab_e_field(k_axis, pol):
    """
    Complex electric field in lab x,y,z coordinates.

    Convention:
      If k = z and quantization axis = z, sigma+ gives q=+1.
    """
    a, b = transverse_basis_for_k(k_axis)

    if pol == "sigma+":
        return -(a + 1j * b) / sqrt(2)
    if pol == "sigma-":
        return (a - 1j * b) / sqrt(2)
    if pol.startswith("linear"):
        ax = pol.split()[-1]
        return unit_vector(ax)

    raise ValueError(pol)


def polarization_ellipse_vector(k_axis, azimuth_deg, ellipticity_deg):
    """Return a normalized Jones vector for an arbitrary pure polarization.

    ``azimuth_deg`` rotates the ellipse's major axis from the first vector of
    the local transverse basis. ``ellipticity_deg`` is in [-45, 45] degrees;
    zero is linear, +45 is sigma+, and -45 is sigma- in this app's convention.
    The overall optical phase is intentionally fixed because it is unobservable.
    """
    azimuth = radians(float(azimuth_deg))
    ellipticity = radians(float(ellipticity_deg))
    if not -45.0 <= float(ellipticity_deg) <= 45.0:
        raise ValueError("Ellipticity angle must lie between -45 and 45 degrees.")

    e1, e2 = transverse_basis_for_k(k_axis)
    component_1 = (
        cos(ellipticity) * cos(azimuth)
        - 1j * sin(ellipticity) * sin(azimuth)
    )
    component_2 = (
        cos(ellipticity) * sin(azimuth)
        + 1j * sin(ellipticity) * cos(azimuth)
    )
    vector = component_1 * e1 + component_2 * e2
    norm = float(np.linalg.norm(vector))
    if norm <= 0.0:
        raise ValueError("Polarization vector has zero norm.")
    return np.asarray(vector / norm, dtype=complex)


def preset_ellipse_parameters(k_axis, pol):
    """Return ellipse angles that reproduce one of the preset polarizations."""
    if pol == "sigma+":
        return 0.0, 45.0
    if pol == "sigma-":
        return 0.0, -45.0
    if pol.startswith("linear"):
        linear_axis = pol.split()[-1]
        target = unit_vector(linear_axis)
        e1, e2 = transverse_basis_for_k(k_axis)
        if np.allclose(target, e1):
            return 0.0, 0.0
        if np.allclose(target, e2):
            return 90.0, 0.0
    raise ValueError(f"{pol!r} is not an allowed polarization for k={k_axis!r}.")


def stokes_from_ellipse(azimuth_deg, ellipticity_deg):
    """Return normalized (s1, s2, s3) with s3=+1 for sigma+."""
    azimuth = radians(float(azimuth_deg))
    ellipticity = radians(float(ellipticity_deg))
    return np.array(
        [
            cos(2.0 * ellipticity) * cos(2.0 * azimuth),
            cos(2.0 * ellipticity) * sin(2.0 * azimuth),
            sin(2.0 * ellipticity),
        ],
        dtype=float,
    )


def local_components(E_lab, q_axis):
    """
    Components in a local frame with local z along the chosen quantization axis.
    """
    if q_axis == "z":
        ux, uy, uz = unit_vector("x"), unit_vector("y"), unit_vector("z")
    elif q_axis == "x":
        ux, uy, uz = unit_vector("y"), unit_vector("z"), unit_vector("x")
    elif q_axis == "y":
        ux, uy, uz = unit_vector("z"), unit_vector("x"), unit_vector("y")
    else:
        raise ValueError(q_axis)

    return np.array([
        np.vdot(ux, E_lab),
        np.vdot(uy, E_lab),
        np.vdot(uz, E_lab),
    ], dtype=complex)


def spherical_components_from_lab(E_lab, q_axis):
    """Return complex q=-1,0,+1 components in the app's convention."""
    Ex, Ey, Ez = local_components(np.asarray(E_lab, dtype=complex), q_axis)
    E_plus = -(Ex - 1j * Ey) / sqrt(2)
    E_zero = Ez
    E_minus = (Ex + 1j * Ey) / sqrt(2)
    return {-1: E_minus, 0: E_zero, +1: E_plus}


def spherical_weights_from_lab(E_lab, q_axis):
    """Return normalized spherical intensities for a lab-frame Jones vector."""
    components = spherical_components_from_lab(E_lab, q_axis)
    weights = {q: float(abs(value) ** 2) for q, value in components.items()}
    total = sum(weights.values())
    if total <= 0.0:
        return {-1: 0.0, 0: 0.0, +1: 0.0}
    return {q: value / total for q, value in weights.items()}


def tensor_geometry_E20(E_lab, q_axis):
    """Return Mathur's dimensionless rank-2 q=0 polarization factor."""
    pi_fraction = spherical_weights_from_lab(E_lab, q_axis)[0]
    return float((1.0 - 3.0 * pi_fraction) / sqrt(6.0))


def spherical_weights_relative_to_quant_axis(k_axis, pol, q_axis):
    """
    Return |E_q|^2 for q=-1,0,+1 relative to the selected quantization axis.
    """
    return spherical_weights_from_lab(lab_e_field(k_axis, pol), q_axis)
