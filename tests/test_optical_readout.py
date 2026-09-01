import unittest

import numpy as np

from alkali_pumping_app.physics import (
    ATOMS,
    build_ground_states,
    mathur_quadrupole_scale,
    weak_probe_readout_operators,
)
from alkali_pumping_app.physics.optical_pumping import build_optical_L
from alkali_pumping_app.physics.optical_readout import (
    _single_detuning_polarizability,
    polarizability_tensor_operators,
)


def probe_operators(**overrides):
    atom = ATOMS["Rb87"]
    values = {
        "atom": atom,
        "ground_states": build_ground_states(atom),
        "q_axis": "z",
        "line": "D1",
        "detuning_MHz": 1200.0,
        "k_axis": "x",
        "azimuth_deg": 35.0,
        "ellipticity_deg": 5.0,
        "path_length_cm": 1.0,
        "density_cm3": 1e11,
        "n2_pressure_torr": 50.0,
        "temperature_C": 40.0,
        "n2_width_MHz_per_torr": 17.8,
        "n2_shift_MHz_per_torr": -8.25,
        "include_scalar": False,
        "include_orientation": True,
        "include_alignment": True,
    }
    values.update(overrides)
    return weak_probe_readout_operators(**values)


class OpticalReadoutTests(unittest.TestCase):
    def test_mathur_quadrupole_conversion_for_supported_manifolds(self):
        expected = {
            1.0: -1.0,
            2.0: -1.0 / np.sqrt(21.0),
            3.0: -1.0 / np.sqrt(126.0),
            4.0: -1.0 / np.sqrt(462.0),
        }
        for F, scale in expected.items():
            self.assertAlmostEqual(mathur_quadrupole_scale(F), scale, places=14)

    def test_alignment_absorption_matches_direct_zeeman_sum(self):
        """Protect the sign and normalization of the LDOR tensor response."""
        atom = ATOMS["Rb87"]
        ground_states = build_ground_states(atom)
        detuning_MHz = -1808.355979089
        response = _single_detuning_polarizability(
            atom, "D1", detuning_MHz, 0.0, 23.0, 17.8, -8.25
        )
        tensors = polarizability_tensor_operators(
            atom, ground_states, "z", response
        )
        _, info = build_optical_L(
            atom=atom,
            line="D1",
            ground_states=ground_states,
            detuning_MHz=detuning_MHz,
            pump_rate_s=1.0,
            selected_transition=None,
            k_axis="x",
            pol="linear z",
            q_axis="z",
            n2_pressure_torr=0.0,
            temperature_C=23.0,
            n2_width_MHz_per_torr=17.8,
            n2_shift_MHz_per_torr=-8.25,
            normalize_to_selected_total=False,
        )
        absorption = info["R_ge"].sum(axis=1)
        F = 2.0
        indices = [
            index for index, state in enumerate(ground_states)
            if np.isclose(float(state["F"]), F)
        ]
        m = np.asarray([ground_states[index]["m"] for index in indices])
        raw_qzz = m**2 - F * (F + 1.0) / 3.0
        direct_scalar, direct_tensor = np.linalg.lstsq(
            np.column_stack([np.ones_like(raw_qzz), raw_qzz]),
            absorption[indices],
            rcond=None,
        )[0]
        reference = next(
            index for index in indices
            if np.isclose(float(ground_states[index]["m"]), 0.0)
        )
        model_scalar = tensors["scalar"][2, 2, reference, reference].imag
        model_tensor = (
            response["alpha_br"][F].imag * mathur_quadrupole_scale(F)
        )
        self.assertGreater(direct_scalar, 0.0)
        self.assertGreater(direct_tensor, 0.0)
        self.assertAlmostEqual(
            direct_tensor / direct_scalar,
            model_tensor / model_scalar,
            places=12,
        )

    def test_all_signal_operators_are_hermitian(self):
        operators, _ = probe_operators()
        for component in operators.values():
            for operator in component.values():
                np.testing.assert_allclose(operator, operator.conj().T, atol=1e-18)

    def test_total_is_sum_of_enabled_rank_contributions(self):
        operators, _ = probe_operators(
            include_scalar=True, include_orientation=True, include_alignment=True
        )
        for signal, total in operators["total"].items():
            expected = sum(operators[rank][signal] for rank in (
                "scalar", "orientation", "alignment"
            ))
            np.testing.assert_allclose(total, expected, rtol=1e-12, atol=1e-30)

    def test_zero_density_or_path_gives_zero_readout(self):
        for override in ({"density_cm3": 0.0}, {"path_length_cm": 0.0}):
            operators, _ = probe_operators(**override)
            for component in operators.values():
                for operator in component.values():
                    np.testing.assert_allclose(operator, 0.0, atol=0.0)

    def test_orientation_and_alignment_channels_exist(self):
        operators, _ = probe_operators()
        self.assertGreater(
            max(np.max(np.abs(value)) for value in operators["orientation"].values()),
            0.0,
        )
        self.assertGreater(
            max(np.max(np.abs(value)) for value in operators["alignment"].values()),
            0.0,
        )

    def test_rotation_is_reported_unavailable_for_circular_input(self):
        operators, info = probe_operators(ellipticity_deg=45.0)
        self.assertFalse(info["signal_availability"]["rotation_available"])
        np.testing.assert_allclose(operators["total"]["rotation"], 0.0, atol=0.0)


if __name__ == "__main__":
    unittest.main()
