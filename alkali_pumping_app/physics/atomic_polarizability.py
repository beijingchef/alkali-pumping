"""Mathur polarizability response functions for alkali D lines.

The implementation follows Eqs. (41)-(43), (53), and (60) of
Mathur, Tang, and Happer, Phys. Rev. A 2, 648 (1970).  Only the diagonal
ground-hyperfine responses ``alpha_gt(F,F)`` and ``alpha_br(F,F)`` are
returned because those are the responses present without hyperfine coherence.
"""

from functools import lru_cache
from math import ceil, floor, log, pi, sqrt

import numpy as np
from scipy.special import wofz
from sympy import S
from sympy.physics.wigner import racah

from .angular_momentum import (
    allowed_F,
    hfs_energy_MHz,
    hyperfine_transition_allowed,
)
from .spectroscopy import hyperfine_transition_choices


CLASSICAL_ELECTRON_RADIUS_CM = 2.8179403262e-13
ATOMIC_MASS_KG = 1.66053906660e-27
BOLTZMANN_J_K = 1.380649e-23
LIGHT_SPEED_M_S = 299792458.0
ELECTRON_MASS_KG = 9.1093837139e-31
ELEMENTARY_CHARGE_C = 1.602176634e-19
VACUUM_PERMITTIVITY_F_M = 8.8541878128e-12


def mathur_quadrupole_scale(F):
    """Convert Mathur's normalized rank-2 operator to raw Cartesian Q_ij."""
    F = float(F)
    if F < 1.0:
        return 0.0
    denominator = (
        F * (F + 1.0) * (2.0 * F - 1.0)
        * (2.0 * F + 1.0) * (2.0 * F + 3.0)
    )
    return -sqrt(30.0 / denominator)


def default_polarizability_sweep_range_MHz(
    atom,
    line,
    n2_pressure_torr,
    n2_coeffs,
    margin_MHz=1500.0,
    rounding_MHz=100.0,
):
    """Return an outward-rounded range around all allowed transitions."""
    transitions = hyperfine_transition_choices(
        atom,
        line,
        n2_pressure_torr,
        n2_coeffs,
        allowed_only=True,
    )
    centers = np.asarray([row["detP"] for row in transitions], dtype=float)
    lower = rounding_MHz * floor(
        (float(centers.min()) - float(margin_MHz)) / rounding_MHz
    )
    upper = rounding_MHz * ceil(
        (float(centers.max()) + float(margin_MHz)) / rounding_MHz
    )
    return float(lower), float(upper)


def oscillator_strength_from_natural_width(atom, line):
    """Return the absorption oscillator strength inferred from the D-line width."""
    wavelength_m = float(atom[f"lambda_{line}_nm"]) * 1e-9
    angular_frequency = 2.0 * pi * LIGHT_SPEED_M_S / wavelength_m
    decay_rate_s = 2.0 * pi * float(atom[line]["gamma_nat_MHz"]) * 1e6
    ground_degeneracy = 2.0 * float(atom["ground"]["J"]) + 1.0
    excited_degeneracy = 2.0 * float(atom[line]["Jp"]) + 1.0
    return float(
        decay_rate_s
        * 2.0
        * pi
        * VACUUM_PERMITTIVITY_F_M
        * ELECTRON_MASS_KG
        * LIGHT_SPEED_M_S**3
        / (ELEMENTARY_CHARGE_C**2 * angular_frequency**2)
        * excited_degeneracy
        / ground_degeneracy
    )


def mathur_scale_G_cm3(atom, line, temperature_C):
    """Return Mathur's common polarizability scale G in cm^3."""
    temperature_K = float(temperature_C) + 273.15
    if temperature_K <= 0.0:
        raise ValueError("Temperature must be above absolute zero.")
    wavelength_cm = float(atom[f"lambda_{line}_nm"]) * 1e-7
    mass_kg = float(atom["mass_amu"]) * ATOMIC_MASS_KG
    thermal_factor = LIGHT_SPEED_M_S / sqrt(
        2.0 * BOLTZMANN_J_K * temperature_K / mass_kg
    )
    oscillator_strength = oscillator_strength_from_natural_width(atom, line)
    return float(
        wavelength_cm**2
        * CLASSICAL_ELECTRON_RADIUS_CM
        * oscillator_strength
        * thermal_factor
        / (8.0 * pi**2)
    )


@lru_cache(maxsize=512)
def _racah(a, b, c, d, e, f):
    """Evaluate a Racah W coefficient using exact half-integer arguments."""
    return float(racah(*(S(str(value)) for value in (a, b, c, d, e, f))))


def _phase(exponent):
    rounded = int(round(float(exponent)))
    if not np.isclose(float(exponent), rounded):
        raise ValueError("Angular-momentum phase exponent is not integral.")
    return -1.0 if rounded % 2 else 1.0


def _plasma_dispersion_response(
    detunings_MHz,
    transition_center_MHz,
    lorentz_fwhm_MHz,
    doppler_fwhm_MHz,
):
    """Return Z(x+iy), with Re(Z) dispersive and Im(Z) absorptive."""
    sigma_MHz = float(doppler_fwhm_MHz) / (2.0 * sqrt(2.0 * log(2.0)))
    sigma_MHz = max(sigma_MHz, 1e-15)
    denominator = sigma_MHz * sqrt(2.0)
    x = (np.asarray(detunings_MHz, dtype=float) - float(transition_center_MHz)) / denominator
    y = max(float(lorentz_fwhm_MHz), 1e-15) / (2.0 * denominator)
    return 1j * sqrt(pi) * wofz(x + 1j * y)


def calculate_atomic_polarizability_sweep(
    atom,
    line,
    detunings_MHz,
    n2_pressure_torr,
    temperature_C,
    n2_width_MHz_per_torr,
    n2_shift_MHz_per_torr,
):
    """Return the four complex Mathur polarizability responses over detuning.

    Detuning is measured from the zero-pressure fine-structure line center.
    Values are in cm^3.  ``alpha_gt`` and ``alpha_br`` are dictionaries keyed
    by the ground hyperfine quantum number F and contain the diagonal F,F
    response functions relevant when ground-hyperfine coherence is absent.
    """
    detunings = np.asarray(detunings_MHz, dtype=float)
    if detunings.ndim != 1 or len(detunings) < 2:
        raise ValueError("detunings_MHz must be a one-dimensional sweep.")
    if line not in ("D1", "D2"):
        raise ValueError("line must be D1 or D2.")
    if float(n2_pressure_torr) < 0.0:
        raise ValueError("N2 pressure must be nonnegative.")

    temperature_K = float(temperature_C) + 273.15
    if temperature_K <= 0.0:
        raise ValueError("Temperature must be above absolute zero.")

    I = float(atom["I"])
    Jg = float(atom["ground"]["J"])
    Je = float(atom[line]["Jp"])
    ground_F = [float(value) for value in allowed_F(I, Jg)]
    excited_F = [float(value) for value in allowed_F(I, Je)]

    ground_energy = {
        F: hfs_energy_MHz(
            I,
            Jg,
            F,
            float(atom["ground"]["A"]),
            float(atom["ground"]["B"]),
        )
        for F in ground_F
    }
    excited_energy = {
        F: hfs_energy_MHz(
            I,
            Je,
            F,
            float(atom[line]["A"]),
            float(atom[line]["B"]),
        )
        for F in excited_F
    }

    pressure_torr = float(n2_pressure_torr)
    pressure_shift_MHz = float(n2_shift_MHz_per_torr) * pressure_torr
    lorentz_fwhm_MHz = (
        float(atom[line]["gamma_nat_MHz"])
        + float(n2_width_MHz_per_torr) * pressure_torr
    )
    frequency_Hz = LIGHT_SPEED_M_S / (float(atom[f"lambda_{line}_nm"]) * 1e-9)
    mass_kg = float(atom["mass_amu"]) * ATOMIC_MASS_KG
    doppler_fwhm_MHz = (
        2.0
        * frequency_Hz
        * sqrt(2.0 * BOLTZMANN_J_K * temperature_K * log(2.0) / (mass_kg * LIGHT_SPEED_M_S**2))
        / 1e6
    )
    G = mathur_scale_G_cm3(atom, line, temperature_C)

    response = {}
    transition_centers = {}
    for Fg in ground_F:
        for Fe in excited_F:
            center = excited_energy[Fe] - ground_energy[Fg] + pressure_shift_MHz
            if hyperfine_transition_allowed(Fg, Fe):
                transition_centers[(Fg, Fe)] = float(center)
            response[(Fg, Fe)] = _plasma_dispersion_response(
                detunings,
                center,
                lorentz_fwhm_MHz,
                doppler_fwhm_MHz,
            )

    A0 = {}
    for Fg in ground_F:
        total = np.zeros_like(detunings, dtype=complex)
        for Fe in excited_F:
            W = _racah(Je, Fe, Jg, Fg, I, 1.0)
            total += (2.0 * Fe + 1.0) * W**2 * response[(Fg, Fe)]
        A0[Fg] = 2.0 * sqrt(3.0) * G * sqrt(2.0 * Fg + 1.0) * total

    alpha_eq = sum(
        sqrt(2.0 * Fg + 1.0) * A0[Fg] for Fg in ground_F
    ) / (2.0 * sqrt(3.0) * (2.0 * I + 1.0))
    alpha_hfs = (
        2.0
        / (sqrt(3.0) * (2.0 * I + 1.0))
        * sum(
            A0[Fg]
            / sqrt(2.0 * Fg + 1.0)
            * _phase(I + 0.5 - Fg)
            for Fg in ground_F
        )
    )

    alpha_gt = {}
    alpha_br = {}
    for Fg in ground_F:
        gt_total = np.zeros_like(detunings, dtype=complex)
        br_total = np.zeros_like(detunings, dtype=complex)
        gt_denominator = _racah(1.0, 0.5, Fg, I, 0.5, Fg)
        if np.isclose(gt_denominator, 0.0):
            raise ValueError(f"Undefined gyrotropic normalization for F={Fg:g}.")
        for Fe in excited_F:
            electronic_W = _racah(Je, Fe, Jg, Fg, I, 1.0)
            z = response[(Fg, Fe)]
            gt_total += (
                _phase(Fe - Fg)
                * z
                * (2.0 * Fe + 1.0)
                * _racah(1.0, 1.0, Fg, Fg, 1.0, Fe)
                / gt_denominator
                * electronic_W**2
            )
            br_total += (
                _phase(Fe - Fg - 1.0)
                * z
                * (2.0 * Fg + 1.0)
                * (2.0 * Fe + 1.0)
                * _racah(1.0, 1.0, Fg, Fg, 2.0, Fe)
                * electronic_W**2
            )
        alpha_gt[Fg] = -6.0 * G * gt_total
        alpha_br[Fg] = 6.0 * G * br_total

    return {
        "detunings_MHz": detunings,
        "alpha_eq": np.asarray(alpha_eq, dtype=complex),
        "alpha_hfs": np.asarray(alpha_hfs, dtype=complex),
        "alpha_gt": alpha_gt,
        "alpha_br": alpha_br,
        "ground_F": ground_F,
        "transition_centers_MHz": transition_centers,
        "G_cm3": G,
        "oscillator_strength": oscillator_strength_from_natural_width(atom, line),
        "lorentz_fwhm_MHz": lorentz_fwhm_MHz,
        "doppler_fwhm_MHz": doppler_fwhm_MHz,
        "pressure_shift_MHz": pressure_shift_MHz,
    }
