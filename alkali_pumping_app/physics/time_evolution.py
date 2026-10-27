"""Coherent pump-on / pump-off evolution, with population reference helpers.

The initial state is unpolarized: equal probability per ground Zeeman state.
The density solver retains all within-F coherences and signed Larmor dynamics.
The population helpers below provide diagonal reference limits. Collision
targets recombine instantaneous nuclear/electron marginals rather than using
an equilibrium Jacobian. Optical terms switch off exactly at t1; the complete
atomic state is continuous across that switch.
"""

import numpy as np
from scipy.integrate import solve_ivp

from .constants import ATOMS, cross_spin_exchange_rate_info
from .spin_exchange import hyperfine_uncoupled_probabilities
from .system_model import build_species_linear_model


TIME_EVOLUTION_MODEL_VERSION = "v3-coherent-ground-hfs-secular-master-equation"
TIME_EVOLUTION_SOLVER_VERSION = "v2-prepared-rhs-analytic-jacobian-and-rate-selection"


def build_population_dynamics(species_A, species_B, beams, common):
    """Build population generators without calculating any RF/probe response."""
    if species_B is not None and species_B["atom_name"] == species_A["atom_name"]:
        species_B = None
    active_beams = [beam for beam in beams
                    if beam.get("target_label") != "B" or species_B is not None]
    for beam in active_beams:
        intensity = float(beam["intensity"])
        if not np.isfinite(intensity) or intensity < 0:
            raise ValueError("Pump intensities must be finite and nonnegative.")
    models = []
    for config in (species_A, species_B):
        if config is None:
            continue
        for key in ("R_ER", "density_cm3"):
            if not np.isfinite(float(config[key])) or float(config[key]) < 0:
                raise ValueError(f"{key} must be finite and nonnegative.")
        model = build_species_linear_model(config, active_beams, common)
        probabilities, _, _ = hyperfine_uncoupled_probabilities(
            model["atom"], model["ground_states"]
        )
        model.update(
            collision_probabilities=probabilities,
            nuclear_by_state=probabilities.sum(axis=2),
            electron_by_state=probabilities.sum(axis=1),
            L_dark=model["L_relaxation"],
            R_SE_cross=0.0,
        )
        models.append(model)
    if len(models) == 2:
        for index, model in enumerate(models):
            partner = models[1 - index]
            model["R_SE_cross"] = cross_spin_exchange_rate_info(
                model["atom_name"], partner["atom_name"], common["temperature_C"],
                partner["density_cm3"],
            )["rate_s"]
    return models


def population_derivative(population, models, pumps_on):
    """Evaluate the coupled nonlinear map using both species at the same time."""
    boundaries = np.cumsum([0, *[model["size"] for model in models]])
    local = [population[boundaries[i]:boundaries[i + 1]]
             for i in range(len(models))]
    nuclear = [p @ model["nuclear_by_state"] for p, model in zip(local, models)]
    electron = [(p @ model["electron_by_state"]) / p.sum()
                for p, model in zip(local, models)]
    derivatives = []
    for index, (model, p) in enumerate(zip(models, local)):
        probabilities = model["collision_probabilities"]
        self_target = np.einsum("i,s,ais->a", nuclear[index], electron[index], probabilities)
        derivative = model["L_dark"] @ p + model["R_SE_self"] * (self_target - p)
        if pumps_on:
            derivative = derivative + model["L_optical"] @ p
        if len(models) == 2:
            cross_target = np.einsum(
                "i,s,ais->a", nuclear[index], electron[1 - index], probabilities
            )
            derivative = derivative + model["R_SE_cross"] * (cross_target - p)
        derivatives.append(derivative)
    return np.concatenate(derivatives)


def population_moments(model, population):
    """Scalar/longitudinal moments per atom, without manifold renormalization."""
    states = model["ground_states"]
    F = np.array([state["F"] for state in states], dtype=float)
    m = np.array([state["m"] for state in states], dtype=float)
    electron_weights = model["electron_by_state"][:, 1] - model["electron_by_state"][:, 0]
    # F = I + S implies S.I/hbar^2 = [F(F+1)-I(I+1)-S(S+1)]/2.
    # This scalar is constant within each F block, so within-F coherences
    # do not contribute. It is not a normalized difference of populations.
    nuclear_spin = float(model["atom"]["I"])
    electron_spin = 0.5
    hyperfine_weights = .5 * (F * (F + 1) - nuclear_spin * (nuclear_spin + 1)
                              - electron_spin * (electron_spin + 1))
    result = {"electron_polarization": population @ electron_weights,
              "hyperfine_polarization": population @ hyperfine_weights}
    for name, value in (("Fa", F.max()), ("Fb", F.min())):
        mask = np.isclose(F, value)
        result[f"population_{name}"] = population[:, mask].sum(axis=1)
        result[f"orientation_{name}"] = population @ (mask * m)
        # Q_qq / hbar^2 is traceless alignment, not the raw second moment.
        result[f"alignment_{name}"] = population @ (mask * (m**2 - value * (value + 1) / 3))
    return result


def density_moments(model, density_matrix):
    """Coherent moments: laboratory orientation and cyclic local alignment.

    The density matrix is expressed in the output-q basis. Its local Cartesian
    axes are (u,v,q); laboratory operators must therefore be rotated into that
    basis before taking Tr(rho O). All values are normalized by hbar**rank and
    are per total atom, including the population weight of each F block.
    """
    from .coherent_evolution import spin_matrices

    rho = np.asarray(density_matrix, dtype=complex)
    result = population_moments(model, np.diagonal(rho, axis1=1, axis2=2).real)
    spins = spin_matrices(model["ground_states"])
    rotation = model["basis_rotation"]
    lab_spins = [rotation.conj().T @ spin @ rotation for spin in spins]
    Fu, Fv, _ = spins
    alignment_uv = (Fu @ Fv + Fv @ Fu) / 2
    alignment_difference = Fu @ Fu - Fv @ Fv
    Fs = np.array([state["F"] for state in model["ground_states"]])
    for name, F in (("Fa", Fs.max()), ("Fb", Fs.min())):
        mask = np.isclose(Fs, F)
        operators = {f"orientation_{axis}": spin for axis, spin in zip("xyz", lab_spins)}
        operators.update(alignment_uv=alignment_uv, alignment_uu_minus_vv=alignment_difference)
        for observable, operator in operators.items():
            result[f"{observable}_{name}"] = np.einsum(
                "tij,ji->t", rho, mask[:, None] * operator,
            ).real
    return result


def density_display_view(species, zeeman_axis, moment_axis):
    """Read one physical trajectory in independent population/moment frames.

    U_axis maps the local axis basis into laboratory z. For an input basis U0,
    rho_axis = (U_axis^dagger U0) rho_input (U0^dagger U_axis).
    Populations and exported coherences use zeeman_axis; Sq and alignment use
    moment_axis. Laboratory Fx,Fy,Fz and scalar S.I remain invariant.
    """
    from .coherent_evolution import basis_rotation

    states = species["ground_states"]
    source_rotation = basis_rotation(states, species["q_axis"])
    zeeman_rotation = basis_rotation(states, zeeman_axis)
    moment_rotation = basis_rotation(states, moment_axis)

    def rotate(target):
        transform = target.conj().T @ source_rotation
        return transform @ species["density_matrix"] @ transform.conj().T

    zeeman_density = rotate(zeeman_rotation)
    moment_density = zeeman_density if moment_axis == zeeman_axis else rotate(moment_rotation)
    atom = ATOMS[species["atom_name"]]
    probabilities, _, _ = hyperfine_uncoupled_probabilities(atom, states)
    moment_model = dict(atom=atom, ground_states=states, electron_by_state=probabilities.sum(axis=1),
                        basis_rotation=moment_rotation)
    return dict(species, q_axis=zeeman_axis, zeeman_axis=zeeman_axis, moment_axis=moment_axis,
                density_matrix=zeeman_density,
                population=np.diagonal(zeeman_density, axis1=1, axis2=2).real.copy(),
                moments=density_moments(moment_model, moment_density))


def _density_integration_method(models, pumps_on, duration):
    """Use DOP853 for sustained coherent motion without faster damping.

    Infinity norms bound coherent and optical-damping generator rates in s^-1.
    Twice the collision rates conservatively accounts for the two marginal
    factors in the spin-exchange target. This is a performance heuristic, not
    a change to the master equation or integration accuracy. Apply it to the
    complete coupled system and separately on either side of pump turn-off.
    BDF remains the default for damping-dominated or short intervals.
    """
    coherent_rate, damping_rate = 0., 0.
    for model in models:
        hamiltonian = model["density_dark"]
        damping = (model["R_u"] + 2 * (model["R_ER"] + model["R_SE_self"]
                                      + model["R_SE_cross"]))
        if pumps_on:
            light_shift = model["density_light_shift"]
            hamiltonian = hamiltonian + light_shift
            damping += np.linalg.norm(model["density_optical"] - light_shift, ord=np.inf)
        coherent_rate = max(coherent_rate, float(np.linalg.norm(hamiltonian, ord=np.inf)))
        damping_rate = max(damping_rate, float(damping))
    if coherent_rate * duration >= 2 * np.pi and damping_rate <= coherent_rate:
        return "DOP853"
    return "BDF"


def compute_density_evolution(species_A, species_B, beams, common, t1_s, t2_s,
                              sample_count=801):
    """Evolve full within-F density blocks; q axes specify output bases only.

    The laboratory field drives signed-g_F Larmor precession. Pumps contribute
    coherent absorption/repopulation and AC-Stark evolution until the switch.
    Optical detunings omit magnetic shifts. No RF or probe enters this equation.
    """
    from .coherent_evolution import prepare_density_model, DensityDynamics
    t1_s, t2_s = float(t1_s), float(t2_s)
    if not np.isfinite([t1_s, t2_s]).all() or not 0 <= t1_s <= t2_s or t2_s <= 0:
        raise ValueError("Times must satisfy 0 ≤ t1 ≤ t2, with finite t2 > 0.")
    if sample_count < 3:
        raise ValueError("At least three time samples are required.")
    configs = [config for config in (species_A, species_B) if config is not None]
    lab_configs = [dict(config, q_axis="z") for config in configs]
    models = build_population_dynamics(lab_configs[0], lab_configs[1] if len(lab_configs) == 2 else None,
                                       beams, common)
    for model, config in zip(models, configs):
        prepare_density_model(model, common, config["q_axis"])
    dynamics = DensityDynamics(models)
    initial = np.concatenate([(np.eye(m["size"]) / m["size"]).ravel()[m["density_indices"]]
                              for m in models]).astype(complex)
    # Zeeman populations can contain every tensor rank through 2*F_max.
    fastest_harmonic = max(2 * max(s["F"] for s in m["ground_states"]) * m["max_larmor_hz"]
                           for m in models)
    requested_samples = max(int(sample_count), int(np.ceil(20 * fastest_harmonic * t2_s)) + 1)
    sample_count = min(requested_samples, 20001)
    if 0 < t1_s < t2_s:
        first_count = min(sample_count - 1, max(2, int(round((sample_count - 1) * t1_s / t2_s)) + 1))
        times = np.r_[np.linspace(0, t1_s, first_count),
                      np.linspace(t1_s, t2_s, sample_count - first_count + 1)[1:]]
    else:
        times = np.linspace(0, t2_s, sample_count)
    density = np.empty((len(times), len(initial)), complex)
    current = initial
    for start, stop, pumps_on in ((0.0, t1_s, True), (t1_s, t2_s, False)):
        if stop <= start:
            continue
        indices = np.flatnonzero((times >= start) & (times <= stop))
        method = _density_integration_method(models, pumps_on, stop - start)
        jacobian = lambda time, rho, active=pumps_on: dynamics.jacobian(rho, active)
        solution = solve_ivp(
            lambda time, rho, active=pumps_on: dynamics.derivative(rho, active),
            (start, stop), current, method=method, t_eval=times[indices],
            rtol=2e-10, atol=1e-12,
            **({"jac": jacobian} if method == "BDF" else {}),
        )
        if method == "DOP853" and (not solution.success or not np.isfinite(solution.y).all()):
            # Restart the complete interval from its original state and keep
            # exactly the same output grid and tolerances on the BDF retry.
            solution = solve_ivp(
                lambda time, rho, active=pumps_on: dynamics.derivative(rho, active),
                (start, stop), current, method="BDF", t_eval=times[indices],
                rtol=2e-10, atol=1e-12, jac=jacobian,
            )
        if not solution.success:
            raise RuntimeError("Density-matrix integration failed: " + solution.message)
        density[indices] = solution.y.T
        current = solution.y[:, -1]
    result = {"time_s": times, "t1_s": t1_s, "t2_s": t2_s, "A": None, "B": None,
              "model": TIME_EVOLUTION_MODEL_VERSION,
              "initial_state": "equal population per ground Zeeman state", "rf_enabled": False,
              "sampling_limited": requested_samples > sample_count}
    offset = 0
    for model in models:
        flat = np.zeros((len(times), model["size"]**2), complex)
        flat[:, model["density_indices"]] = density[:, offset:offset + model["density_size"]]
        offset += model["density_size"]
        rho = flat.reshape(-1, model["size"], model["size"])
        rotation = model["basis_rotation"]
        rho = rotation.conj().T @ rho @ rotation
        if not np.isfinite(rho).all():
            raise RuntimeError("Density integration produced a nonfinite state.")
        minimum_population = float(np.diagonal(rho, axis1=1, axis2=2).real.min())
        trace_error = float(np.max(np.abs(np.trace(rho, axis1=1, axis2=2) - 1)))
        hermiticity_error = float(np.max(np.abs(rho - rho.conj().transpose(0, 2, 1))))
        rho = .5 * (rho + rho.conj().transpose(0, 2, 1))
        eigenvalues, vectors = np.linalg.eigh(rho)
        minimum_eigenvalue = float(eigenvalues.min())
        if (trace_error > 1e-6 or hermiticity_error > 1e-6
                or minimum_eigenvalue < -1e-7):
            raise RuntimeError("Density integration did not preserve trace, Hermiticity, or positivity.")
        # Correct only verified numerical roundoff; preserve every coherence.
        eigenvalues = np.maximum(eigenvalues, 0)
        eigenvalues /= eigenvalues.sum(axis=1, keepdims=True)
        rho = (vectors * eigenvalues[:, None, :]) @ vectors.conj().transpose(0, 2, 1)
        local = np.diagonal(rho, axis1=1, axis2=2).real.copy()
        result[model["label"]] = {
            "atom_name": model["atom_name"], "ground_states": model["ground_states"],
            "q_axis": model["output_axis"], "population": local, "density_matrix": rho,
            "moments": density_moments(model, rho),
            "density_cm3": model["density_cm3"], "R_SE_self": model["R_SE_self"],
            "R_u": model["R_u"],
            "R_SE_cross": model["R_SE_cross"], "max_trace_error": trace_error,
            "minimum_population_before_roundoff": minimum_population,
            "minimum_eigenvalue_before_roundoff": minimum_eigenvalue,
            "max_hermiticity_error": hermiticity_error,
        }
    return result


# Numerical population references above remain useful for the collinear limit.
# Existing Python callers receive the coherent calculation and its populations.
compute_population_evolution = compute_density_evolution
