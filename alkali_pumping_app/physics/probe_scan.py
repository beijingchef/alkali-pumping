"""Pump-on stationary states for a non-perturbing probe frequency scan."""

import numpy as np
from scipy.integrate import solve_ivp

from .coherent_evolution import density_derivative, prepare_density_model
from .time_evolution import build_population_dynamics


def pump_on_steady_state(species_A, species_B, beams, common):
    """Converge the same nonlinear coherent master equation with pumps on.

    Start from the unpolarized state, preserving the transient model's choice
    of conserved quantities even when a stationary state is not unique. Use
    exponentially longer autonomous BDF intervals, independently of t1/t2.
    Accept only a stationary residual AND two consecutive small state changes.
    Refine integration accuracy when a nearly unchanged state has a stalled
    residual; longer intervals cannot remove integration error at that floor.
    A persistent oscillation or integration failure is reported, not plotted
    as a stationary spectrum.
    """
    first = dict(species_A, q_axis="z")
    second = dict(species_B, q_axis="z") if species_B is not None else None
    models = build_population_dynamics(first, second, beams, common)
    for model in models:
        prepare_density_model(model, common, "z")
    boundaries = np.cumsum([0, *[model["density_size"] for model in models]])
    state = np.concatenate([
        (np.eye(model["size"]) / model["size"]).ravel()[model["density_indices"]]
        for model in models
    ]).astype(complex)
    slow_rates, roundoff, initial_roundoff = [], [], []
    for model in models:
        optical_rate = float(np.linalg.norm(model["density_optical"], ord=np.inf))
        rates = [model["R_ER"], model["R_u"], model["R_SE_self"], model["R_SE_cross"], optical_rate]
        positive = [float(rate) for rate in rates if rate > 0]
        slow_rates.append(min(positive) if positive else 1.)
        initial_roundoff.append(100 * np.finfo(float).eps * max([*rates, 1e-30]))
        fastest = max([*rates, float(np.linalg.norm(model["density_dark"], ord=np.inf)), 1e-30])
        roundoff.append(100 * np.finfo(float).eps * fastest)
    tolerance = np.array(slow_rates) * 1e-10 + roundoff

    def residual(values):
        derivative = density_derivative(values, models, True)
        return np.array([np.max(np.abs(derivative[a:b])) for a, b in zip(boundaries[:-1], boundaries[1:])])

    elapsed, stable = 0., 0
    rtol, atol = 2e-10, 1e-12
    refinements = 0
    initial_duration = 1 / min(slow_rates)
    duration = initial_duration
    previous_residual = residual(state)
    # An exactly stationary unpolarized state needs no integration.
    converged = bool(np.all(previous_residual <= np.array(slow_rates) * 1e-14 + initial_roundoff))
    for _ in range(30):
        if converged:
            break
        solution = solve_ivp(
            lambda time, values: density_derivative(values, models, True),
            (0., duration), state, method="BDF", t_eval=[duration], rtol=rtol, atol=atol,
        )
        if not solution.success:
            raise RuntimeError("Pump-on steady state failed: " + solution.message)
        candidate = solution.y[:, -1]
        change = float(np.max(np.abs(candidate - state)))
        state = candidate
        elapsed += duration
        current_residual = residual(state)
        stationary = np.all(current_residual <= tolerance)
        stable = stable + 1 if stationary and change <= 1e-9 else 0
        converged = stable >= 2
        # The residual has units s^-1, whereas BDF controls state error.
        # Small numerical coherences can therefore pass BDF's error test but
        # fail stationarity, especially in the presence of rapid precession.
        stalled = change <= 1e-9 and np.any(
            (current_residual > tolerance) & (current_residual >= .5 * previous_residual)
        )
        previous_residual = current_residual
        if stalled:
            if refinements >= 3:
                raise RuntimeError(
                    "Pump-on steady-state residual stalled after three accuracy refinements; "
                    "no steady-state scan is shown."
                )
            refinements += 1
            atol *= .01
            rtol = max(rtol * .01, 2.5e-14)
            # The state has already settled: refine over one initial interval
            # instead of repeating an exponentially enlarged interval.
            duration = initial_duration
        elif change > 1e-9:
            duration *= 2
    if not converged:
        raise RuntimeError("The pump-on state did not converge to a stationary state; no steady-state scan is shown.")
    result = {"A": None, "B": None, "convergence_time_s": elapsed,
              "stationary_residual_s-1": float(np.max(residual(state)))}
    for model, start, stop in zip(models, boundaries[:-1], boundaries[1:]):
        flat = np.zeros(model["size"]**2, complex)
        flat[model["density_indices"]] = state[start:stop]
        rho = flat.reshape(model["size"], model["size"])
        trace = np.trace(rho)
        if (not np.isfinite(rho).all() or abs(trace - 1) > 1e-7
                or np.max(np.abs(rho - rho.conj().T)) > 1e-7
                or np.min(np.linalg.eigvalsh((rho + rho.conj().T) / 2)) < -1e-8):
            raise RuntimeError("Pump-on steady state did not preserve a physical density matrix.")
        rho = (rho + rho.conj().T) / 2
        rho /= np.trace(rho).real
        result[model["label"]] = dict(
            atom_name=model["atom_name"], ground_states=model["ground_states"], q_axis="z",
            density_matrix=rho[None], density_cm3=model["density_cm3"],
            R_ER=model["R_ER"], R_u=model["R_u"],
            R_SE_self=model["R_SE_self"], R_SE_cross=model["R_SE_cross"],
        )
    return result
