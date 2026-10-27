"""Density cache limits and display values avoid duplicated trajectories."""
from copy import deepcopy
import unittest
from unittest.mock import patch
import numpy as np
from alkali_pumping_app.physics import time_evolution as physics
from alkali_pumping_app.ui import time_evolution_cache as cache
from tests.test_time_evolution import inputs


class TimeEvolutionCacheTests(unittest.TestCase):
    def setUp(self):
        for fn in (cache._small_evolution, cache._large_evolution, cache.display_cached, cache._steady_state):
            fn.clear()

    tearDown = setUp

    def test_small_and_large_cache_routing_and_bounded_eviction(self):
        config, beams, common = inputs()
        common.update(static_field_nT=0.)
        self.assertLess(cache.estimated_density_bytes(config, None, common, 4.), cache.LARGE_DENSITY_BYTES)
        large = dict(common, static_field_nT=50.)
        cs = dict(config, atom_name='Cs133')
        self.assertEqual(cache.estimated_density_bytes(cs, None, large, 30.), 20001*16*16*16)
        with patch.object(physics, 'compute_density_evolution', return_value={'marker': 1}) as compute:
            for rate in (10., 11., 10.):
                cache.evolution_cached(dict(config, R_ER=rate), None, beams, common, 1., 4.)
            self.assertEqual(compute.call_count, 2)
            for rate in (20., 21., 20.):
                cache.evolution_cached(dict(config, R_ER=rate), None, beams, large, 1., 4.)
            self.assertEqual(compute.call_count, 5)  # one large entry: each new input evicts it
            cache.evolution_cached(dict(config, R_ER=10.), None, beams, common, 1., 4.)
            self.assertEqual(compute.call_count, 5)  # separate small cache survives

    def test_display_cache_contains_only_observables_and_preserves_curves(self):
        config, beams, common = inputs()
        result = physics.compute_density_evolution(config, None, beams, common, .01, .02, 21)
        species = result['A']
        original = deepcopy(species['density_matrix'])
        expected = physics.density_display_view(species, 'x', 'y')
        compact = cache.display_cached(species, 'x', 'y')
        self.assertNotIn('density_matrix', compact)
        self.assertFalse(any(isinstance(value, np.ndarray) and value.ndim == 3 for value in compact.values()))
        np.testing.assert_array_equal(compact['population'], expected['population'])
        for name in expected['moments']:
            np.testing.assert_array_equal(compact['moments'][name], expected['moments'][name])
        np.testing.assert_array_equal(species['density_matrix'], original)

    def test_solver_version_invalidates_transient_and_stationary_caches(self):
        config, beams, common = inputs()
        with patch.object(physics, 'compute_density_evolution', return_value={'value': 1}) as transient, \
                patch.object(cache.probe_scan, 'pump_on_steady_state', return_value={'value': 2}) as steady:
            for version in ('first', 'first', 'second', 'second'):
                with patch.object(physics, 'TIME_EVOLUTION_SOLVER_VERSION', version):
                    cache.evolution_cached(config, None, beams, common, 1., 4.)
                    cache.steady_cached(config, None, beams, common)
            self.assertEqual(transient.call_count, 2)
            self.assertEqual(steady.call_count, 2)


if __name__ == '__main__':
    unittest.main()
