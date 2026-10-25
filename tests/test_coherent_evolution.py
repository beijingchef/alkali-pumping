"""Coherent precession, covariance, collisions, optical transfer and detection."""

from copy import deepcopy
import unittest
import numpy as np
from scipy.linalg import expm

from alkali_pumping_app.physics.coherent_evolution import (
    basis_rotation, density_derivative, prepare_density_model, spin_matrices,
)
from alkali_pumping_app.physics.time_evolution import (
    build_population_dynamics, compute_density_evolution, density_display_view, density_moments,
)
from alkali_pumping_app.physics.optical_readout import density_probe_signals
from alkali_pumping_app.physics.constants import ATOMS, ground_hyperfine_lande_g, MU_B_OVER_H_HZ_PER_T
from tests.test_time_evolution import inputs


def models_for(config, other, beams, common):
    return [prepare_density_model(m, common, 'z')
            for m in build_population_dynamics(config, other, beams, common)]


def read(species, config, common):
    return density_probe_signals(atom=ATOMS[config['atom_name']],
        ground_states=species['ground_states'], density_matrix=species['density_matrix'],
        q_axis=species['q_axis'], **{k: config['probe'][k] for k in (
            'line', 'detuning_MHz', 'k_axis', 'azimuth_deg', 'ellipticity_deg', 'path_length_cm')},
        density_cm3=config['density_cm3'], temperature_C=common['temperature_C'],
        n2_pressure_torr=common['n2_pressure_torr'], n2_width_MHz_per_torr=17.8,
        n2_shift_MHz_per_torr=-8.25)


class CoherentEvolutionTests(unittest.TestCase):
    def test_independent_display_axes_preserve_state_and_probe(self):
        config, beams, common = inputs()
        common.update(static_field_axis='y', static_field_nT=.1)
        physical = compute_density_evolution(config, None, beams, common, .02, .04, 11)['A']
        original_density = physical['density_matrix'].copy()
        original_signal = read(physical, config, common)
        for moment_axis in 'xyz':
            reference = density_display_view(physical, 'z', moment_axis)
            for zeeman_axis in 'xyz':
                view = density_display_view(physical, zeeman_axis, moment_axis)
                with self.subTest(zeeman=zeeman_axis, moment=moment_axis):
                    rotation = basis_rotation(view['ground_states'], zeeman_axis)
                    np.testing.assert_allclose(rotation @ view['density_matrix'] @ rotation.conj().T,
                                               original_density, atol=1e-12)
                    np.testing.assert_allclose(view['population'].sum(axis=1), 1, atol=1e-12)
                    for observable in view['moments']:
                        np.testing.assert_allclose(view['moments'][observable], reference['moments'][observable], atol=1e-12)
                    signal = read(view, config, common)
                    for name in ('rotation_rad', 'ellipticity_rad', 'absorption'):
                        np.testing.assert_allclose(signal[name], original_signal[name], atol=1e-12)
        np.testing.assert_array_equal(physical['density_matrix'], original_density)

    def test_cartesian_and_local_alignment_for_known_coherent_spin_states(self):
        config, _, common = inputs()
        model = models_for(config, None, [], common)[0]
        states, size = model['ground_states'], model['size']
        _, Fy, Fz = spin_matrices(states)
        direction = np.array([1., 2., 3.]) / np.sqrt(14.)
        theta, phi = np.arccos(direction[2]), np.arctan2(direction[1], direction[0])
        turn = expm(-1j * phi * Fz) @ expm(-1j * theta * Fy)
        rho = np.zeros((size, size), complex)
        Fs = sorted({state['F'] for state in states})
        for F, weight in ((Fs[-1], .7), (Fs[0], .3)):
            index = next(i for i, state in enumerate(states) if state['F'] == F and state['m'] == F)
            vector = turn[:, index]
            rho += weight * np.outer(vector, vector.conj())
        for q, (u, v) in {'x': ('y', 'z'), 'y': ('z', 'x'), 'z': ('x', 'y')}.items():
            rotation = basis_rotation(states, q)
            local_rho = rotation.conj().T @ rho @ rotation
            moments = density_moments(dict(model, basis_rotation=rotation), local_rho[None])
            nq, nu, nv = (direction['xyz'.index(axis)] for axis in (q, u, v))
            for name, F, weight in (('Fa', Fs[-1], .7), ('Fb', Fs[0], .3)):
                with self.subTest(q=q, F=F):
                    for axis, component in zip('xyz', direction):
                        self.assertAlmostEqual(moments[f'orientation_{axis}_{name}'][0], weight * F * component)
                    # A spin coherent state has Qij/hbar^2 = F(F-1/2)(ni*nj-deltaij/3).
                    scale = weight * F * (F - .5)
                    self.assertAlmostEqual(moments[f'alignment_{name}'][0], scale * (nq**2 - 1/3))
                    self.assertAlmostEqual(moments[f'alignment_uv_{name}'][0], scale * nu * nv)
                    self.assertAlmostEqual(moments[f'alignment_uu_minus_vv_{name}'][0], scale * (nu**2 - nv**2))
                    self.assertAlmostEqual(moments[f'orientation_{name}'][0], weight * F * nq)

    def test_signed_larmor_rotation_and_uniform_coherence_decay(self):
        config, beams, common = inputs()
        config.update(density_cm3=0., R_ER=0., R_u=2.3)
        common.update(static_field_axis='x', static_field_nT=1.)
        model = models_for(config, None, [], common)[0]
        n, states = model['size'], model['ground_states']
        Fx, Fy, Fz = spin_matrices(states)
        for F in sorted({s['F'] for s in states}):
            rho = np.zeros((n, n), complex)
            i = next(i for i, s in enumerate(states) if s['F'] == F and s['m'] == F)
            rho[i, i] = 1
            y = rho.ravel()[model['density_indices']]
            derivative = density_derivative(y, [model], False)
            d = np.zeros(n*n, complex)
            d[model['density_indices']] = derivative
            d = d.reshape(n, n)
            omega = 2*np.pi*MU_B_OVER_H_HZ_PER_T*1e-9*ground_hyperfine_lande_g(config['atom_name'], model['atom'], F)
            self.assertAlmostEqual(np.trace(d @ Fy).real, -omega*F, places=10)
            self.assertAlmostEqual(np.trace(d @ Fz).real, -2.3*F, places=10)
            np.testing.assert_allclose(d, d.conj().T, atol=1e-13)
            self.assertAlmostEqual(np.trace(d).real, 0, places=13)
            # Exact generator evolution preserves coherences and positive eigenvalues.
            generator = np.column_stack([density_derivative(np.eye(len(y), dtype=complex)[:, j], [model], False)
                                         for j in range(len(y))])
            evolved = np.zeros(n*n, complex)
            evolved[model['density_indices']] = expm(.013 * generator) @ y
            evolved = evolved.reshape(n, n)
            self.assertAlmostEqual(np.trace(evolved @ Fz).real, F*np.cos(omega*.013)*np.exp(-2.3*.013), places=11)
            self.assertGreaterEqual(np.linalg.eigvalsh(evolved).min(), -1e-12)

    def test_pumping_creates_coherences_transverse_field_precesses_after_switch(self):
        config, beams, common = inputs()
        config.update(density_cm3=0., R_ER=0., R_u=0.)
        common.update(static_field_axis='x', static_field_nT=.2)
        result = compute_density_evolution(config, None, beams, common, .03, .3, 101)
        rho, times = result['A']['density_matrix'], result['time_s']
        model = models_for(config, None, [], common)[0]
        switch = np.flatnonzero(times == .03)[0]
        self.assertGreater(np.max(np.abs(rho - np.eye(8)[None] * rho)), 1e-4)
        for index in (switch+1, len(times)-1):
            # Build exact dark unitary from the signed Zeeman Hamiltonian.
            g = np.array([ground_hyperfine_lande_g(config['atom_name'], model['atom'], s['F']) for s in model['ground_states']])
            h = 2*np.pi*MU_B_OVER_H_HZ_PER_T*1e-9*.2*g[:, None]*spin_matrices(model['ground_states'])[0]
            u = expm(-1j*h*(times[index]-.03))
            np.testing.assert_allclose(rho[index], u @ rho[switch] @ u.conj().T, atol=2e-9)
        np.testing.assert_allclose(np.trace(rho, axis1=1, axis2=2), 1, atol=1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(), -1e-12)

    def test_quantization_basis_changes_no_lab_observables_in_dual_species(self):
        config, beams, common = inputs()
        other = deepcopy(config)
        other.update(label='B', atom_name='Rb85', density_cm3=2e10)
        common.update(static_field_axis='y', static_field_nT=.1)
        beams[0].update(k_axis='x', pol='sigma+')
        reference = compute_density_evolution(config, other, beams, common, .02, .05, 21)
        original_signal = read(reference['A'], config, common)
        for axis in ('x', 'y'):
            config['q_axis'] = axis
            other['q_axis'] = 'y' if axis == 'x' else 'x'
            result = compute_density_evolution(config, other, beams, common, .02, .05, 21)
            for label in ('A', 'B'):
                species = result[label]
                u = basis_rotation(species['ground_states'], species['q_axis'])
                np.testing.assert_allclose(u @ species['density_matrix'] @ u.conj().T,
                                           reference[label]['density_matrix'], atol=1e-11)
                for name in ('population_Fa', 'population_Fb', 'hyperfine_polarization',
                             *[f'orientation_{axis}_{manifold}' for axis in 'xyz' for manifold in ('Fa', 'Fb')]):
                    np.testing.assert_allclose(species['moments'][name], reference[label]['moments'][name], atol=1e-12)
            signal = read(result['A'], config, common)
            for name in ('rotation_rad', 'ellipticity_rad', 'absorption'):
                np.testing.assert_allclose(signal[name], original_signal[name], atol=1e-12)

    def test_optical_jump_population_limit_and_collective_covariance(self):
        config, beams, common = inputs()
        common['static_field_nT'] = 0.
        config.update(density_cm3=0., R_ER=0.)
        model = models_for(config, None, beams, common)[0]
        p = np.arange(1., 9); p /= p.sum()
        y = np.diag(p).ravel()[model['density_indices']]
        flat = np.zeros(64, complex)
        flat[model['density_indices']] = model['density_optical'] @ y
        np.testing.assert_allclose(np.diag(flat.reshape(8, 8)).real, model['L_optical'] @ p, atol=1e-11)
        # Rotating both the state and pump must rotate the optical derivative.
        u = basis_rotation(model['ground_states'], 'x')
        rotated_rho = u @ np.diag(p) @ u.conj().T
        beams[0].update(k_axis='x', pol='sigma+')
        rotated = models_for(config, None, beams, common)[0]
        flat2 = np.zeros(64, complex)
        flat2[rotated['density_indices']] = rotated['density_optical'] @ rotated_rho.ravel()[rotated['density_indices']]
        np.testing.assert_allclose(flat2.reshape(8, 8), u @ flat.reshape(8, 8) @ u.conj().T, atol=1e-10)

    def test_coherent_collisions_preserve_total_angular_momentum(self):
        config, beams, common = inputs()
        config.update(R_ER=0., R_u=0.)
        common['static_field_nT'] = 0.
        models = models_for(config, None, [], common)
        model = models[0]
        u = basis_rotation(model['ground_states'], 'y')
        p = np.arange(1., 9)**2; p /= p.sum()
        rho = u @ np.diag(p) @ u.conj().T
        flat = np.zeros(64, complex)
        flat[model['density_indices']] = density_derivative(rho.ravel()[model['density_indices']], models, False)
        for spin in spin_matrices(model['ground_states']):
            self.assertAlmostEqual(np.trace(flat.reshape(8, 8) @ spin).real, 0, places=11)

    def test_probe_detects_offdiagonal_state_and_does_not_modify_it(self):
        config, beams, common = inputs()
        result = compute_density_evolution(config, None, beams, common, .03, .03, 11)
        species = result['A']
        u = basis_rotation(species['ground_states'], 'x')
        species['density_matrix'] = u @ species['density_matrix'] @ u.conj().T
        before = species['density_matrix'].copy()
        config['probe'].update(k_axis='x', azimuth_deg=0.)
        coherent = read(species, config, common)
        np.testing.assert_array_equal(before, species['density_matrix'])
        species['density_matrix'] = before * np.eye(8)[None]
        diagonal = read(species, config, common)
        self.assertGreater(np.max(np.abs(coherent['rotation_rad'] - diagonal['rotation_rad'])), 1e-8)

    def test_cross_exchange_conserves_density_weighted_angular_momentum(self):
        config, beams, common = inputs()
        config.update(R_ER=0., R_u=0.)
        other = deepcopy(config)
        other.update(label='B', atom_name='Rb85', density_cm3=3e10)
        common['static_field_nT'] = 0.
        models = models_for(config, other, [], common)
        local = []
        for m, axis in zip(models, ('x', 'y')):
            m['R_SE_self'] = 0.
            p = np.arange(1., m['size']+1)**2; p /= p.sum()
            u = basis_rotation(m['ground_states'], axis)
            local.append((u @ np.diag(p) @ u.conj().T).ravel()[m['density_indices']])
        derivative = density_derivative(np.concatenate(local), models, False)
        offset, angular = 0, np.zeros(3)
        for m in models:
            flat = np.zeros(m['size']**2, complex)
            flat[m['density_indices']] = derivative[offset:offset+m['density_size']]
            offset += m['density_size']
            d = flat.reshape(m['size'], m['size'])
            self.assertAlmostEqual(np.trace(d).real, 0, places=12)
            angular += m['density_cm3'] * np.array([np.trace(d @ f).real for f in spin_matrices(m['ground_states'])])
        np.testing.assert_allclose(angular / sum(m['density_cm3'] for m in models), 0, atol=1e-12)

    def test_all_isotopes_and_D2_preserve_density_invariants(self):
        from alkali_pumping_app.physics.spectroscopy import line_center_frequency_MHz
        for atom_name in ATOMS:
            config, beams, common = inputs()
            config.update(atom_name=atom_name, density_cm3=0., R_ER=2., R_u=1., q_axis='y')
            common.update(static_field_axis='x', static_field_nT=.1)
            beams[0].update(absolute_frequency_MHz=line_center_frequency_MHz(ATOMS[atom_name], 'D2'),
                            k_axis='z', pol='linear x')
            result = compute_density_evolution(config, None, beams, common, .004, .01, 11)
            rho = result['A']['density_matrix']
            with self.subTest(atom=atom_name):
                np.testing.assert_allclose(np.trace(rho, axis1=1, axis2=2), 1, atol=1e-12)
                np.testing.assert_allclose(rho, rho.conj().transpose(0, 2, 1), atol=1e-12)
                self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(), -1e-12)

    def test_complete_export_contains_independent_complex_coherences(self):
        from alkali_pumping_app.ui.time_evolution_data import evolution_export_dataframe
        config, beams, common = inputs()
        common.update(static_field_axis='x', static_field_nT=.1)
        result = compute_density_evolution(config, None, beams, common, .02, .04, 11)
        table = evolution_export_dataframe(result, 'A', read(result['A'], config, common),
                                           config, config['probe'], common, beams)
        names = [name for name in table if name.startswith('rho_')]
        self.assertEqual(len(names), 26)  # 2*(5 choose 2 + 3 choose 2)
        self.assertGreater(np.max(np.abs(table[names].to_numpy())), 1e-6)
        self.assertTrue((table['alignment_u_axis'] == 'x').all())
        self.assertTrue((table['alignment_v_axis'] == 'y').all())
        for name in ('orientation_x_Fa', 'orientation_y_Fb', 'orientation_z_Fa',
                     'alignment_uv_Fa', 'alignment_uu_minus_vv_Fb'):
            suffix = '_hbar' if name.startswith('orientation') else '_hbar2'
            np.testing.assert_array_equal(table[name + suffix], result['A']['moments'][name])
        states = result['A']['ground_states']
        for i, s in enumerate(states):
            for j, t in enumerate(states):
                if j > i and s['F'] == t['F']:
                    name = f"rho_F{s['F']:g}_m{s['m']:g}_m{t['m']:g}"
                    np.testing.assert_array_equal(table[name+'_real'] + 1j*table[name+'_imag'],
                                                   result['A']['density_matrix'][:, i, j])


if __name__ == '__main__':
    unittest.main()
