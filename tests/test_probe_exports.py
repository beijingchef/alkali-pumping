import unittest

import numpy as np

from alkali_pumping_app.ui.exports import weak_probe_export_dataframe


class WeakProbeExportTests(unittest.TestCase):
    @staticmethod
    def _response(signal):
        return {
            rank: {
                signal: {
                    "amplitude": np.array([1.0, 2.0]),
                    "in_phase": np.array([-3.0, 4.0]),
                    "quadrature": np.array([5.0, -6.0]),
                }
            }
            for rank in ("total", "scalar", "orientation", "alignment")
        }

    def test_rotation_is_converted_from_per_rabi_frequency_to_rad_per_nT(self):
        frame = weak_probe_export_dataframe(
            [10.0, 20.0],
            self._response("rotation"),
            "rotation",
            rf_rabi_rad_s_per_nT=2.5,
        )

        np.testing.assert_allclose(
            frame["total_in_phase_per_nT"].to_numpy(), [-7.5, 10.0]
        )
        self.assertTrue((frame["units"] == "rad/nT").all())
        self.assertTrue((frame["rf_rabi_rad_s_per_nT"] == 2.5).all())

    def test_normalized_stokes_signal_uses_inverse_nanotesla(self):
        frame = weak_probe_export_dataframe(
            [10.0, 20.0],
            self._response("s3"),
            "s3",
            rf_rabi_rad_s_per_nT=4.0,
        )

        np.testing.assert_allclose(
            frame["alignment_quadrature_per_nT"].to_numpy(), [20.0, -24.0]
        )
        self.assertTrue((frame["units"] == "1/nT").all())


if __name__ == "__main__":
    unittest.main()
