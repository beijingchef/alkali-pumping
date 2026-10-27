"""Ground-HFS secular master equation in a common laboratory z basis.

All m,m' coherences within each F are retained. GHz inter-F coherences and
excited-state dynamics are adiabatically eliminated. Collective spontaneous
emission jumps retain polarization transfer; resolving individual m emissions
would spuriously select a quantization axis. See docs/technical/time-evolution.qmd.
"""

import numpy as np
from scipy.linalg import expm

from .angular_momentum import build_excited_states, coupled_basis_amplitudes, dipole_amplitude
from .constants import MU_B_OVER_H_HZ_PER_T, ground_hyperfine_lande_g
from .optical_pumping import optical_rate_scale_from_intensity
from .polarization import lab_e_field, spherical_components_from_lab
from .spectroscopy import complex_voigt_response_relative, doppler_fwhm_MHz, transition_shift_MHz


def spin_matrices(states):
    """Dimensionless F_x,F_y,F_z in the standard |F,m_z> convention."""
    raising = np.zeros((len(states), len(states)), complex)
    lookup = {(s['F'], s['m']): i for i, s in enumerate(states)}
    for i, s in enumerate(states):
        j = lookup.get((s['F'], s['m'] + 1))
        if j is not None:
            raising[j, i] = np.sqrt(s['F'] * (s['F'] + 1) - s['m'] * (s['m'] + 1))
    return ((raising + raising.T) / 2, (raising - raising.T) / (2j),
            np.diag([s['m'] for s in states]).astype(complex))


def basis_rotation(states, q_axis):
    """Columns are local-q basis vectors in lab-z; cyclic transverse frames."""
    if q_axis not in ('x', 'y', 'z'):
        raise ValueError('Quantization axis must be x, y, or z.')
    angle = {'z': 0., 'x': 2 * np.pi / 3, 'y': -2 * np.pi / 3}[q_axis]
    return expm(-1j * angle * sum(spin_matrices(states)) / np.sqrt(3))


def _optical_terms(model, common):
    """Weak-excitation Lindblad jumps and pump AC-Stark H/hbar (s^-1)."""
    atom, states, n = model['atom'], model['ground_states'], model['size']
    jumps, h_light = [], np.zeros((n, n), complex)
    for beam in model['beams']:
        if beam['intensity'] == 0:
            continue
        line = beam['line']
        coeff = model['n2_coeffs'][line]
        excited = build_excited_states(atom, line)
        dipoles = np.array([[[dipole_amplitude(atom['I'], atom['ground']['J'], atom[line]['Jp'],
                                              g['F'], g['m'], e['F'], e['m'], q)
                             for g in states] for e in excited] for q in (-1, 0, 1)])
        spherical = spherical_components_from_lab(lab_e_field(beam['k_axis'], beam['pol']), 'z')
        absorption = sum(spherical[q] * dipoles[q + 1] for q in (-1, 0, 1))
        rate = optical_rate_scale_from_intensity(atom, line, beam['intensity'],
                    common['n2_pressure_torr'], common['temperature_C'], coeff['width'])
        width = atom[line]['gamma_nat_MHz'] + coeff['width'] * common['n2_pressure_torr']
        doppler = doppler_fwhm_MHz(atom, line, common['temperature_C'])
        for Fe in sorted({e['F'] for e in excited}):
            ei = np.array([i for i, e in enumerate(excited) if e['F'] == Fe])
            # Sum_q D_q^dagger D_q = I_e, including all outgoing F.
            decay_norm = np.sum(dipoles[:, ei, :] ** 2, axis=(0, 2))
            for Fi in sorted({g['F'] for g in states}):
                gi = np.array([i for i, g in enumerate(states) if g['F'] == Fi])
                a = np.zeros((len(ei), n), complex)
                a[:, gi] = absorption[np.ix_(ei, gi)]
                delta = beam['detuning'] - transition_shift_MHz(states[gi[0]], excited[ei[0]])
                response = complex_voigt_response_relative(
                    delta - coeff['shift'] * common['n2_pressure_torr'], width, doppler)
                # Same Voigt scattering and dispersive normalization as Light shift.
                h_light += .5 * rate * response.imag * (a.conj().T @ a)
                for Fo in sorted({g['F'] for g in states}):
                    mask = np.array([g['F'] == Fo for g in states])
                    for dipole in dipoles:
                        emission = dipole[ei].T / np.sqrt(decay_norm)[None, :]
                        emission = emission * mask[:, None]
                        jump = np.sqrt(rate * max(0., response.real)) * emission @ a
                        if np.any(jump):
                            jumps.append(jump)
    return jumps, h_light


def prepare_density_model(model, common, output_axis):
    n, states = model['size'], model['ground_states']
    mask = np.array([[a['F'] == b['F'] for b in states] for a in states])
    indices = np.flatnonzero(mask.ravel())
    spins = spin_matrices(states)
    field = float(common['static_field_nT'])
    if not np.isfinite(field) or common['static_field_axis'] not in ('x', 'y', 'z'):
        raise ValueError('Static field must be finite with axis x, y, or z.')
    g = np.array([ground_hyperfine_lande_g(model['atom_name'], model['atom'], s['F']) for s in states])
    h_b = 2 * np.pi * MU_B_OVER_H_HZ_PER_T * 1e-9 * field * g[:, None] * spins['xyz'.index(common['static_field_axis'])]
    jumps, h_light = _optical_terms(model, common)
    identity = np.eye(n)
    def commutator(h):
        return -1j * (np.kron(h, identity) - np.kron(identity, h.T))
    light_shift = commutator(h_light)
    optical = light_shift.copy()
    for jump in jumps:
        loss = jump.conj().T @ jump
        optical += np.kron(jump, jump.conj()) - .5 * (np.kron(loss, identity) + np.kron(identity, loss.T))
    amplitudes = coupled_basis_amplitudes(model['atom'], states)
    # Common lab frame makes unlike-species electron transfer basis independent.
    nuclear_map = np.einsum('ais,bjs->ijab', amplitudes, amplitudes).reshape(-1, n*n)[:, indices]
    electron_map = np.einsum('ais,bit->stab', amplitudes, amplitudes).reshape(4, n*n)[:, indices]
    recombine = np.einsum('ais,bjt->abijst', amplitudes, amplitudes).reshape(n*n, -1)[indices]
    model.update(density_indices=indices, density_size=len(indices), nuclear_map=nuclear_map,
                 electron_map=electron_map, recombine=recombine,
                 nuclear_size=amplitudes.shape[1], output_axis=output_axis,
                 basis_rotation=basis_rotation(states, output_axis),
                 density_dark=commutator(h_b)[np.ix_(indices, indices)],
                 density_light_shift=light_shift[np.ix_(indices, indices)],
                 density_optical=optical[np.ix_(indices, indices)],
                 max_larmor_hz=float(np.max(np.abs(g))) * MU_B_OVER_H_HZ_PER_T * 1e-9 * abs(field))
    return model


class DensityDynamics:
    """Prepared master equation and exact complex Jacobian for one solve.

    For a collision target T(n,e)=R(n tensor e), differentiate both marginals.
    Electron normalization contributes d(e/trace e)=(de-e_normalized*dtrace)/trace.
    Keeping that term is essential for BDF's off-trace trial states and for
    cross-species Jacobian blocks. No Hermitian conjugates enter this derivative.
    """

    def __init__(self, models):
        self.models = tuple(models)
        self.boundaries = np.cumsum([0, *[m['density_size'] for m in models]])
        self.slices = tuple(slice(a, b) for a, b in zip(self.boundaries[:-1], self.boundaries[1:]))
        self.size = int(self.boundaries[-1])
        self.generators = {False: [], True: []}
        self.nuclear_identities = []
        self.electron_identity = np.eye(4)
        self.electron_traces = []
        unpolarized_electron = np.array([.5, 0, 0, .5])
        for m in models:
            identity = np.eye(m['density_size'])
            trace = np.eye(m['size']).ravel()[m['density_indices']]
            nuclear_identity = np.eye(m['nuclear_map'].shape[0])
            randomization = m['recombine'] @ np.kron(nuclear_identity, unpolarized_electron[:, None]) @ m['nuclear_map']
            uniform = np.outer(trace / m['size'], trace)
            cross_rate = m['R_SE_cross'] if len(models) == 2 else 0.
            dark = (m['density_dark'] + m['R_ER'] * (randomization - identity)
                    + m['R_u'] * (uniform - identity) - (m['R_SE_self'] + cross_rate) * identity)
            self.generators[False].append(dark)
            self.generators[True].append(dark + m['density_optical'])
            self.nuclear_identities.append(nuclear_identity)
            self.electron_traces.append(m['electron_map'][0] + m['electron_map'][3])

    def _marginals(self, state):
        local = [state[s] for s in self.slices]
        nuclear = [m['nuclear_map'] @ rho for m, rho in zip(self.models, local)]
        electron = [m['electron_map'] @ rho for m, rho in zip(self.models, local)]
        return local, nuclear, electron

    def derivative(self, state, pumps_on):
        local, nuclear, electron = self._marginals(state)
        result = np.empty(self.size, complex)
        for i, (m, rho, generator) in enumerate(zip(self.models, local, self.generators[pumps_on])):
            value = generator @ rho
            if m['R_SE_self']:
                e = electron[i] / (electron[i][0] + electron[i][3])
                value += m['R_SE_self'] * (m['recombine'] @ np.outer(nuclear[i], e).ravel())
            if len(self.models) == 2 and m['R_SE_cross']:
                partner = electron[1-i]
                e = partner / (partner[0] + partner[3])
                value += m['R_SE_cross'] * (m['recombine'] @ np.outer(nuclear[i], e).ravel())
            result[self.slices[i]] = value
        return result

    def jacobian(self, state, pumps_on):
        _, nuclear, electron = self._marginals(state)
        result = np.zeros((self.size, self.size), complex)
        for i, (m, generator) in enumerate(zip(self.models, self.generators[pumps_on])):
            local = generator.copy()
            # Partial derivative of the target with respect to electron state.
            to_electron = m['recombine'] @ np.kron(nuclear[i][:, None], self.electron_identity)
            if m['R_SE_self']:
                trace = electron[i][0] + electron[i][3]
                e = electron[i] / trace
                de = (m['electron_map'] - np.outer(e, self.electron_traces[i])) / trace
                to_nuclear = m['recombine'] @ np.kron(self.nuclear_identities[i], e[:, None])
                local += m['R_SE_self'] * (to_nuclear @ m['nuclear_map'] + to_electron @ de)
            if len(self.models) == 2 and m['R_SE_cross']:
                j = 1-i
                trace = electron[j][0] + electron[j][3]
                e = electron[j] / trace
                de = (self.models[j]['electron_map'] - np.outer(e, self.electron_traces[j])) / trace
                to_nuclear = m['recombine'] @ np.kron(self.nuclear_identities[i], e[:, None])
                local += m['R_SE_cross'] * (to_nuclear @ m['nuclear_map'])
                result[self.slices[i], self.slices[j]] = m['R_SE_cross'] * (to_electron @ de)
            result[self.slices[i], self.slices[i]] = local
        return result


def density_derivative(state, models, pumps_on):
    """Reference RHS for standalone callers; integrators prepare DensityDynamics once."""
    boundaries = np.cumsum([0, *[m['density_size'] for m in models]])
    local = [state[boundaries[i]:boundaries[i+1]] for i in range(len(models))]
    nuclear = [m['nuclear_map'] @ rho for m, rho in zip(models, local)]
    electron = [m['electron_map'] @ rho for m, rho in zip(models, local)]
    derivatives = []
    for i, (m, rho) in enumerate(zip(models, local)):
        def target(e):
            return m['recombine'] @ np.outer(nuclear[i], e).ravel()
        # trace(rho) is conserved; normalize the partner electron marginal only.
        d = m['density_dark'] @ rho
        if pumps_on:
            d = d + m['density_optical'] @ rho
        d += m['R_ER'] * (target(np.array([.5, 0, 0, .5])) - rho)
        if m['R_SE_self']:
            own_e = electron[i] / (electron[i][0] + electron[i][3])
            d += m['R_SE_self'] * (target(own_e) - rho)
        if len(models) == 2 and m['R_SE_cross']:
            partner = electron[1-i]
            d += m['R_SE_cross'] * (target(partner / (partner[0] + partner[3])) - rho)
        full = np.zeros(m['size']**2, complex)
        full[m['density_indices']] = rho
        trace = np.trace(full.reshape(m['size'], m['size']))
        equilibrium = (np.eye(m['size']) * trace / m['size']).ravel()[m['density_indices']]
        d += m['R_u'] * (equilibrium - rho)
        derivatives.append(d)
    return np.concatenate(derivatives)
