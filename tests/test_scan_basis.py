"""Independent scan axes preserve the full density matrix and physical moments."""
import unittest
import numpy as np

from alkali_pumping_app.physics.angular_momentum import build_ground_states
from alkali_pumping_app.physics.constants import ATOMS
from alkali_pumping_app.physics.coherent_evolution import spin_matrices
from alkali_pumping_app.physics.field_geometry import static_field_direction
from alkali_pumping_app.physics.polarization import polarization_ellipse_vector
from alkali_pumping_app.physics.scan_basis import scan_axis_direction, scan_population_view, SCAN_POPULATION_AXES
from alkali_pumping_app.physics.time_evolution import density_display_view
from alkali_pumping_app.ui.time_evolution_data import steady_density_matrix_export_dataframe


class ScanBasisTests(unittest.TestCase):
    def species(self):
        states = build_ground_states(ATOMS['Rb87'])
        rng = np.random.default_rng(91)
        a = rng.normal(size=(len(states), len(states))) + 1j*rng.normal(size=(len(states), len(states)))
        rho = a @ a.conj().T
        rho /= np.trace(rho)
        return dict(atom_name='Rb87', ground_states=states, q_axis='z', density_matrix=rho[None])

    def test_field_sign_zero_and_probe_ellipse_axes(self):
        common = dict(static_field_theta_deg=53, static_field_phi_deg=127, static_field_nT=3.)
        probe = dict(k_axis='z', azimuth_deg=31., ellipticity_deg=19.)
        for strength, sign in ((3., 1.), (-3., -1.), (0., 1.)):
            common['static_field_nT'] = strength
            np.testing.assert_allclose(scan_axis_direction('B direction', common, probe), sign*static_field_direction(common))
        for k in 'xyz':
            for azimuth in (-90., 0., 31., 90.):
                probe.update(k_axis=k, azimuth_deg=azimuth)
                major = scan_axis_direction('probe major', common, probe)
                minor = scan_axis_direction('probe minor', common, probe)
                np.testing.assert_allclose(np.cross(major, minor), np.eye(3)['xyz'.index(k)], atol=1e-14)
                self.assertAlmostEqual(np.dot(major, minor), 0.)
                for ellipticity in (0., 19., 45., -45.):
                    field = polarization_ellipse_vector(k, azimuth, ellipticity)
                    np.testing.assert_allclose(np.outer(field.real, field.real) + np.outer(field.imag, field.imag),
                        np.cos(np.deg2rad(ellipticity))**2*np.outer(major, major) +
                        np.sin(np.deg2rad(ellipticity))**2*np.outer(minor, minor), atol=1e-14)
                    probe['ellipticity_deg'] = ellipticity
                    np.testing.assert_array_equal(scan_axis_direction('probe minor', common, probe), minor)

    def test_all_axes_preserve_matrix_invariants_and_projection(self):
        species = self.species()
        original = species['density_matrix'].copy()
        common = dict(static_field_theta_deg=53, static_field_phi_deg=127, static_field_nT=-3.)
        probe = dict(k_axis='x', azimuth_deg=31.)
        spins = spin_matrices(species['ground_states'])
        m = np.array([state['m'] for state in species['ground_states']])
        for axis in SCAN_POPULATION_AXES:
            with self.subTest(axis=axis):
                direction = scan_axis_direction(axis, common, probe)
                view = scan_population_view(species, axis, direction)
                rho = view['density_matrix'][0]
                np.testing.assert_allclose(rho, rho.conj().T, atol=1e-14)
                self.assertAlmostEqual(np.trace(rho).real, 1.)
                np.testing.assert_allclose(np.linalg.eigvalsh(rho), np.linalg.eigvalsh(original[0]), atol=1e-14)
                projected = sum(value*spin for value, spin in zip(direction, spins))
                self.assertAlmostEqual(float(view['population'][0] @ m), np.trace(original[0] @ projected).real)
                if axis in 'xyz':
                    legacy = density_display_view(species, axis, 'z')
                    np.testing.assert_allclose(view['population'], legacy['population'], atol=1e-14)
                frame = steady_density_matrix_export_dataframe(view, label='A')
                exported = (frame.rho_real.to_numpy() + 1j*frame.rho_imag.to_numpy()).reshape(rho.shape)
                np.testing.assert_array_equal(exported, rho)
                self.assertEqual(set(frame.quantization_axis), {axis})
                np.testing.assert_allclose(frame[['quantization_axis_x','quantization_axis_y','quantization_axis_z']].iloc[0], direction)
        np.testing.assert_array_equal(species['density_matrix'], original)

    def test_invalid_axis_and_direction_are_rejected(self):
        with self.assertRaises(ValueError): scan_axis_direction('invalid', {}, {})
        for vector in ((0.,0.,0.), (1.,2.), (float('nan'),0.,1.)):
            with self.assertRaises(ValueError): scan_population_view(self.species(), 'invalid', vector)
