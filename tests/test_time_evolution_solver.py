"""Automatic transient solvers preserve the coherent master equation."""

from copy import deepcopy
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from alkali_pumping_app.physics import time_evolution as te
from alkali_pumping_app.physics.coherent_evolution import density_derivative, prepare_density_model
from alkali_pumping_app.physics.constants import ATOMS, DEFAULT_N2_COEFFS
from alkali_pumping_app.physics.spectroscopy import (
    absolute_detuning_from_transition_choice, line_center_frequency_MHz,
)
from tests.test_time_evolution import inputs


def prepared(config, other, beams, common):
    return [prepare_density_model(model, common, 'z')
            for model in te.build_population_dynamics(config, other, beams, common)]


class TimeEvolutionSolverTests(unittest.TestCase):
    def assert_physical(self, result):
        for label in ('A', 'B'):
            if result[label] is None:
                continue
            rho = result[label]['density_matrix']
            np.testing.assert_allclose(np.trace(rho, axis1=1, axis2=2), 1., atol=2e-12)
            np.testing.assert_allclose(rho, rho.conj().transpose(0, 2, 1), atol=2e-12)
            self.assertGreaterEqual(float(np.linalg.eigvalsh(rho).min()), -1e-12)

    def test_method_choice_tracks_stiffness_and_pump_turn_off(self):
        config, beams, common = inputs()
        config['R_u'] = 1.
        common.update(static_field_nT=0.)
        no_field = prepared(config, None, [], common)
        self.assertEqual(te._density_integration_method(no_field, True, 1.), 'BDF')
        self.assertEqual(te._density_integration_method(no_field, False, 3.), 'BDF')
        common.update(static_field_nT=5., static_field_theta_deg=90, static_field_phi_deg=0)
        oscillatory = prepared(config, None, beams, common)
        self.assertEqual(te._density_integration_method(oscillatory, True, 1.), 'DOP853')
        self.assertEqual(te._density_integration_method(oscillatory, False, 3.), 'DOP853')
        # Fast relaxation in either species limits the entire coupled solve.
        other = deepcopy(config)
        other.update(label='B', atom_name='Rb85', R_ER=1e5)
        coupled = prepared(config, other, beams, common)
        for pumps_on in (False, True):
            self.assertEqual(te._density_integration_method(coupled, pumps_on, 1.), 'BDF')
        # Strong optical damping disappears at turn-off; the method may change.
        strong = deepcopy(beams)
        resonant_detuning, _ = absolute_detuning_from_transition_choice(
            atom=ATOMS['Rb87'], line='D1', transition_label='1→2', relative_detuning_MHz=0.,
            n2_pressure_torr=common['n2_pressure_torr'], n2_coeffs=config['n2_coeffs'], allowed_only=True,
        )
        strong[0]['absolute_frequency_MHz'] = line_center_frequency_MHz(ATOMS['Rb87'], 'D1') + resonant_detuning
        for beam in strong:
            beam['intensity'] *= 10
        pumped = prepared(config, None, strong, common)
        self.assertEqual(te._density_integration_method(pumped, True, 1.), 'BDF')
        self.assertEqual(te._density_integration_method(pumped, False, 3.), 'DOP853')

    def test_short_intervals_and_fast_spin_exchange_retain_bdf(self):
        config, beams, common = inputs()
        common.update(static_field_nT=5., static_field_theta_deg=90, static_field_phi_deg=90)
        short = prepared(config, None, beams, common)
        self.assertEqual(te._density_integration_method(short, True, 1e-6), 'BDF')
        config['density_cm3'] = 1e15
        fast_exchange = prepared(config, None, beams, common)
        for pumps_on in (False, True):
            self.assertEqual(te._density_integration_method(fast_exchange, pumps_on, 1.), 'BDF')

    def test_bdf_uses_prepared_analytic_jacobian(self):
        from alkali_pumping_app.physics.coherent_evolution import DensityDynamics
        config, beams, common = inputs()
        common.update(static_field_nT=0.)
        with patch.object(te, 'solve_ivp', wraps=solve_ivp) as integrate:
            result = te.compute_density_evolution(config, None, beams, common, .01, .03, 21)
        dynamics = DensityDynamics(prepared(config, None, beams, common))
        for call, pumps_on in zip(integrate.call_args_list, (True, False)):
            self.assertEqual(call.kwargs['method'], 'BDF')
            np.testing.assert_allclose(call.kwargs['jac'](0., call.args[2]),
                                       dynamics.jacobian(call.args[2], pumps_on))
        self.assert_physical(result)

    def test_oscillatory_linear_solution_matches_exact_propagation(self):
        config, beams, common = inputs()
        config.update(density_cm3=0., R_ER=2., R_u=1.)
        common.update(static_field_theta_deg=90, static_field_phi_deg=0, static_field_nT=5.)
        with patch.object(te, 'solve_ivp', wraps=solve_ivp) as integrate:
            result = te.compute_density_evolution(config, None, beams, common, .02, .06, 21)
        self.assertEqual([call.kwargs['method'] for call in integrate.call_args_list], ['DOP853', 'DOP853'])
        model = prepared(config, None, beams, common)[0]
        indices = model['density_indices']
        initial = (np.eye(model['size']) / model['size']).ravel()[indices]
        identity = np.eye(len(initial), dtype=complex)
        generators = [np.column_stack([density_derivative(column, [model], pumps_on)
                                       for column in identity.T]) for pumps_on in (True, False)]
        at_switch = expm(.02 * generators[0]) @ initial
        for index in np.unique([0, len(result['time_s']) // 6,
                                np.flatnonzero(result['time_s'] == .02)[0],
                                len(result['time_s']) // 2, len(result['time_s']) - 1]):
            time = result['time_s'][index]
            exact = (expm(time * generators[0]) @ initial if time <= .02
                     else expm((time - .02) * generators[1]) @ at_switch)
            np.testing.assert_allclose(result['A']['density_matrix'][index].ravel()[indices],
                                       exact, rtol=2e-8, atol=2e-10)
        self.assertEqual(np.count_nonzero(result['time_s'] == .02), 1)
        self.assert_physical(result)

    def test_all_isotopes_agree_with_tighter_bdf_reference(self):
        for atom_name, atom in ATOMS.items():
            with self.subTest(atom=atom_name):
                config, beams, common = inputs()
                config.update(atom_name=atom_name, density_cm3=1e10, R_ER=2., R_u=1.,
                              n2_coeffs=deepcopy(DEFAULT_N2_COEFFS[atom_name]))
                beams[0].update(intensity=2., absolute_frequency_MHz=line_center_frequency_MHz(atom, 'D1') + 500.)
                common.update(static_field_theta_deg=90, static_field_phi_deg=0, static_field_nT=5.)
                with patch.object(te, 'solve_ivp', wraps=solve_ivp) as integrate:
                    actual = te.compute_density_evolution(config, None, beams, common, .012, .04, 31)
                self.assertEqual([c.kwargs['method'] for c in integrate.call_args_list], ['DOP853', 'DOP853'])
                def reference_solve(*args, **kwargs):
                    kwargs.update(method='BDF', rtol=2e-12, atol=1e-14)
                    return solve_ivp(*args, **kwargs)
                with patch.object(te, 'solve_ivp', side_effect=reference_solve):
                    reference = te.compute_density_evolution(config, None, beams, common, .012, .04, 31)
                np.testing.assert_array_equal(actual['time_s'], reference['time_s'])
                np.testing.assert_allclose(actual['A']['density_matrix'], reference['A']['density_matrix'],
                                           rtol=2e-7, atol=3e-10)
                for name in actual['A']['moments']:
                    np.testing.assert_allclose(actual['A']['moments'][name], reference['A']['moments'][name],
                                               rtol=2e-7, atol=3e-9)
                self.assert_physical(actual)

    def test_dual_species_matches_tighter_bdf_reference(self):
        config, beams, common = inputs()
        config['R_u'] = 1.
        other = deepcopy(config)
        other.update(label='B', atom_name='Rb85', density_cm3=2e10,
                     n2_coeffs=deepcopy(DEFAULT_N2_COEFFS['Rb85']))
        common.update(static_field_theta_deg=90, static_field_phi_deg=90, static_field_nT=5.)
        actual = te.compute_density_evolution(config, other, beams, common, .012, .04, 31)
        def reference_solve(*args, **kwargs):
            kwargs.update(method='BDF', rtol=2e-12, atol=1e-14)
            return solve_ivp(*args, **kwargs)
        with patch.object(te, 'solve_ivp', side_effect=reference_solve):
            reference = te.compute_density_evolution(config, other, beams, common, .012, .04, 31)
        for label in ('A', 'B'):
            np.testing.assert_allclose(actual[label]['density_matrix'], reference[label]['density_matrix'],
                                       rtol=2e-7, atol=3e-10)
        self.assert_physical(actual)

    def test_failed_or_nonfinite_explicit_solve_retries_original_interval(self):
        config, beams, common = inputs()
        common.update(static_field_theta_deg=90, static_field_phi_deg=0, static_field_nT=5.)
        for success in (False, True):
            calls = []
            def fail_first(*args, **kwargs):
                calls.append((args, kwargs))
                if len(calls) == 1:
                    return SimpleNamespace(success=success, y=np.full((len(args[2]), 1), np.nan), message='test failure')
                return solve_ivp(*args, **kwargs)
            with self.subTest(success=success), patch.object(te, 'solve_ivp', side_effect=fail_first):
                result = te.compute_density_evolution(config, None, beams, common, .02, .06, 21)
            self.assertEqual([kw['method'] for _, kw in calls], ['DOP853', 'BDF', 'DOP853'])
            self.assertEqual(calls[0][0][1], calls[1][0][1])
            np.testing.assert_array_equal(calls[0][0][2], calls[1][0][2])
            np.testing.assert_array_equal(calls[0][1]['t_eval'], calls[1][1]['t_eval'])
            for _, kwargs in calls:
                self.assertEqual((kwargs['rtol'], kwargs['atol']), (2e-10, 1e-12))
            self.assert_physical(result)

    def test_failure_of_both_solvers_is_reported(self):
        config, beams, common = inputs()
        common.update(static_field_theta_deg=90, static_field_phi_deg=0, static_field_nT=5.)
        failure = SimpleNamespace(success=False, y=np.zeros((34, 1)), message='integration failed')
        with patch.object(te, 'solve_ivp', return_value=failure) as integrate:
            with self.assertRaisesRegex(RuntimeError, 'integration failed'):
                te.compute_density_evolution(config, None, beams, common, .02, .06, 21)
        self.assertEqual([c.kwargs['method'] for c in integrate.call_args_list], ['DOP853', 'BDF'])


if __name__ == '__main__':
    unittest.main()
