"""Prepared derivatives and Jacobians retain the nonlinear collision equation."""
from copy import deepcopy
import unittest
import numpy as np
from alkali_pumping_app.physics.coherent_evolution import (
    DensityDynamics, density_derivative, prepare_density_model,
)
from alkali_pumping_app.physics.constants import ATOMS
from alkali_pumping_app.physics.time_evolution import build_population_dynamics
from tests.test_time_evolution import inputs


class DensityJacobianTests(unittest.TestCase):
    def models(self, atom='Rb87', dual=False):
        config, beams, common = inputs()
        config.update(atom_name=atom, R_u=3., density_cm3=1e10)
        common.update(static_field_axis='y', static_field_nT=.4)
        other = deepcopy(config) if dual else None
        if other is not None:
            other.update(label='B', atom_name='Rb85', density_cm3=2e10, R_ER=17., R_u=7.)
        return [prepare_density_model(m, common, 'z')
                for m in build_population_dynamics(config, other, beams, common)]

    def state(self, models):
        rng = np.random.default_rng(410)
        values = []
        for i, model in enumerate(models):
            n = model['size']
            F = np.array([s['F'] for s in model['ground_states']])
            a = (rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))) * (F[:, None] == F)
            rho = a @ a.conj().T
            rho *= (1.2 + i*.3) / np.trace(rho)
            values.append(rho.ravel()[model['density_indices']])
        return np.concatenate(values)

    def test_prepared_rhs_matches_reference_for_all_isotopes_and_coupled_states(self):
        for atom, dual in [*[(atom, False) for atom in ATOMS], ('Rb87', True)]:
            models = self.models(atom, dual)
            dynamics = DensityDynamics(models)
            state = self.state(models)
            for pumps_on in (False, True):
                with self.subTest(atom=atom, dual=dual, pumps_on=pumps_on):
                    np.testing.assert_allclose(dynamics.derivative(state, pumps_on),
                                               density_derivative(state, models, pumps_on),
                                               rtol=3e-13, atol=3e-12)

    def test_full_complex_jacobian_matches_finite_differences_including_trace_changes(self):
        for dual in (False, True):
            models = self.models(dual=dual)
            dynamics = DensityDynamics(models)
            state = self.state(models)
            eps = 1e-6
            identity = np.eye(len(state), dtype=complex)
            for pumps_on in (False, True):
                exact = dynamics.jacobian(state, pumps_on)
                for phase in (1., 1j):
                    numerical = np.column_stack([
                        (dynamics.derivative(state + eps*phase*e, pumps_on)
                         - dynamics.derivative(state - eps*phase*e, pumps_on)) / (2*eps*phase)
                        for e in identity
                    ])
                    with self.subTest(dual=dual, pumps_on=pumps_on, phase=phase):
                        np.testing.assert_allclose(exact, numerical, rtol=3e-7, atol=3e-7)
                traces = np.concatenate([np.eye(m['size']).ravel()[m['density_indices']] for m in models])
                np.testing.assert_allclose(traces @ exact, 0., atol=2e-11)
                if dual:
                    a, b = dynamics.slices
                    self.assertGreater(np.linalg.norm(exact[a, b]), 0.)
                    self.assertGreater(np.linalg.norm(exact[b, a]), 0.)

    def test_zero_exchange_is_exactly_linear_and_state_independent(self):
        models = self.models()
        models[0].update(R_SE_self=0., R_SE_cross=0.)
        dynamics = DensityDynamics(models)
        state = self.state(models)
        for pumps_on in (False, True):
            jacobian = dynamics.jacobian(state, pumps_on)
            np.testing.assert_allclose(jacobian @ state, dynamics.derivative(state, pumps_on), atol=2e-12)
            np.testing.assert_array_equal(jacobian, dynamics.jacobian(state*1.3, pumps_on))


if __name__ == '__main__':
    unittest.main()
