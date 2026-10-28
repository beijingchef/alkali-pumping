"""Spherical fields retain Cartesian signs, Hamiltonians and saved meaning."""
from copy import deepcopy
import unittest
import numpy as np
from alkali_pumping_app.physics.field_geometry import static_field_direction, static_field_vector_nT
from alkali_pumping_app.physics.coherent_evolution import prepare_density_model, spin_matrices
from alkali_pumping_app.physics.constants import ATOMS, MU_B_OVER_H_HZ_PER_T, ground_hyperfine_lande_g
from alkali_pumping_app.physics.time_evolution import build_population_dynamics
from alkali_pumping_app.physics import compute_alkali_system
from alkali_pumping_app.ui import conditions, time_evolution_conditions as te_conditions
from tests.test_two_manifold_rf import system_inputs


class FieldGeometryTests(unittest.TestCase):
    def test_cardinal_directions_poles_periodicity_and_signed_strength(self):
        for theta, phi, expected in [(0, 73, [0,0,1]), (180, 17, [0,0,-1]),
                                      (90, 0, [1,0,0]), (90,90,[0,1,0]),
                                      (90,180,[-1,0,0]), (90,270,[0,-1,0]), (90,360,[1,0,0])]:
            common = dict(static_field_theta_deg=theta, static_field_phi_deg=phi, static_field_nT=-3.)
            np.testing.assert_allclose(static_field_direction(common), expected, atol=1e-14)
            np.testing.assert_allclose(static_field_vector_nT(common), -3*np.array(expected), atol=1e-14)
        direction = static_field_direction(dict(static_field_theta_deg=53, static_field_phi_deg=127))
        self.assertAlmostEqual(float(direction@direction), 1.)

    def test_oblique_hamiltonian_is_full_lab_cartesian_dot_product_for_all_isotopes(self):
        config, _, common = system_inputs()
        common.update(static_field_theta_deg=61, static_field_phi_deg=43, static_field_nT=2.)
        vector = static_field_vector_nT(common)
        for atom_name, atom in ATOMS.items():
            with self.subTest(atom=atom_name):
                local = dict(config, atom_name=atom_name)
                model = build_population_dynamics(local, None, [], common)[0]
                prepare_density_model(model, common, 'z')
                n, states = model['size'], model['ground_states']
                g = np.array([ground_hyperfine_lande_g(atom_name, atom, s['F']) for s in states])
                H = 2*np.pi*MU_B_OVER_H_HZ_PER_T*1e-9*g[:,None]*sum(B*F for B,F in zip(vector,spin_matrices(states)))
                expected = -1j*(np.kron(H,np.eye(n))-np.kron(np.eye(n),H.T))
                indices = model['density_indices']
                np.testing.assert_allclose(model['density_dark'], expected[np.ix_(indices,indices)], atol=1e-12)
                np.testing.assert_allclose(H,H.conj().T,atol=1e-12)

    def test_acc_projection_and_transverse_guard(self):
        config, _, common = system_inputs()
        for theta, phi, axis, projection, aligned in [(0,0,'z',3.,True), (180,52,'z',-3.,True),
                                                      (90,180,'x',-3.,True), (60,30,'z',1.5,False)]:
            with self.subTest(theta=theta,phi=phi,axis=axis):
                common.update(static_field_theta_deg=theta,static_field_phi_deg=phi)
                result=compute_alkali_system(dict(config,q_axis=axis),None,[],common)['A']
                self.assertAlmostEqual(result['field_parallel_nT'],projection)
                self.assertEqual(result['static_field_aligned'],aligned)
                self.assertEqual(result['rf_available'],aligned)

    def test_oblique_dark_evolution_matches_exact_unitary_precession(self):
        from scipy.linalg import expm
        from alkali_pumping_app.physics.time_evolution import compute_density_evolution
        from tests.test_time_evolution import inputs
        config, beams, common = inputs()
        config.update(density_cm3=0., R_ER=0., R_u=0.)
        common.update(static_field_theta_deg=57,static_field_phi_deg=123,static_field_nT=2.)
        result=compute_density_evolution(config,None,beams,common,.01,.04,31)
        rho=result['A']['density_matrix']
        states=result['A']['ground_states']
        g=np.array([ground_hyperfine_lande_g(config['atom_name'],ATOMS[config['atom_name']],s['F']) for s in states])
        H=2*np.pi*MU_B_OVER_H_HZ_PER_T*1e-9*g[:,None]*sum(B*F for B,F in zip(static_field_vector_nT(common),spin_matrices(states)))
        switch=np.flatnonzero(result['time_s']==.01)[0]
        self.assertGreater(np.max(np.abs(rho[switch]-np.eye(len(states))/len(states))),1e-8)
        for index in (switch+1,len(rho)-1):
            U=expm(-1j*H*(result['time_s'][index]-.01))
            np.testing.assert_allclose(rho[index],U@rho[switch]@U.conj().T,atol=2e-10)
        np.testing.assert_allclose(np.trace(rho,axis1=1,axis2=2),1.,atol=1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(rho).min(),-1e-12)

    def test_new_condition_formats_round_trip_and_reject_invalid_or_old_angles(self):
        for defaults, build, parse, field, old_version in [
            (conditions.DEFAULT_STARTUP_CONDITION,conditions.build_condition_payload,conditions.apply_loaded_condition_dict,'conditions','6.14'),
            (te_conditions.DEFAULTS,te_conditions.build_settings_payload,te_conditions.parse_settings_payload,'settings','2.5')]:
            from unittest.mock import patch
            values=deepcopy(defaults)
            values.update(static_field_theta_deg=147,static_field_phi_deg=317)
            payload=build(values)
            with patch.object(conditions.st,'session_state',{}):
                parsed=parse(payload)
                if field=='conditions':
                    self.assertEqual(conditions.st.session_state['static_field_theta_deg'],147)
                else:
                    self.assertEqual(parsed['static_field_phi_deg'],317)
                for key, bad in [('static_field_theta_deg',181),('static_field_phi_deg',361),
                                  ('static_field_theta_deg',-1),('static_field_phi_deg',30.5),
                                  ('static_field_theta_deg',True),('static_field_phi_deg',float('nan'))]:
                    invalid=deepcopy(payload)
                    invalid[field][key]=bad
                    with self.assertRaises(ValueError): parse(invalid)
                incomplete=deepcopy(payload)
                incomplete[field].pop('static_field_theta_deg')
                with self.assertRaises(ValueError): parse(incomplete)
                obsolete=deepcopy(payload)
                obsolete[field]['static_field_axis']='z'
                with self.assertRaises(ValueError): parse(obsolete)
                old=deepcopy(payload)
                old['version']=old_version
                with self.assertRaises(ValueError): parse(old)
