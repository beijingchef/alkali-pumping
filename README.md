# Alkali Pumping v6.9.0

This Streamlit application models optical pumping, electron randomization,
self spin exchange, and unlike-alkali spin exchange in a one- or two-alkali
vapor.

The app starts on **Analysis → Light shift**, and the top navigation separates
four workflows:

- **Analysis → Light shift** has independent sidebar controls and settings
  files for atom, cell, static field, arbitrary beam polarization, detuning,
  and plot selection. Its main panel plots scalar, vector, tensor,
  state-resolved, eigenvalue, adjacent-transition, and scattering results.
- **Analysis → Atomic polarizability** plots the complex Mathur equilibrium,
  hyperfine, gyrotropic, and birefringent response functions for D1 or D2
  light. The real curves describe phase response and the imaginary curves
  describe attenuation; diagonal gyrotropic and birefringent responses are
  shown for each ground hyperfine manifold.
- **Analysis → Magnetometry** provides dual-alkali population and Zeeman
  diagnostics, independent RF-A/RF-B drives, atomic-moment susceptibilities,
  and independent weak Probe-A/Probe-B Stokes readouts.
- **Reference → Atomic properties** presents the former Settings dialog as a
  full page with thermal, buffer-gas, and transition-strength tabs.

Light-shift settings persist while navigating between pages and can be saved
or loaded with their own JSON format. Custom polarization is specified by the
beam-frame ellipse azimuth and signed ellipticity, with spherical fractions
and the dimensionless `E^(2)_0` geometry factor shown separately from shift
units. Switching the polarization input from Preset to Ellipse initializes the
ellipse angles to the equivalent current preset. The Components view includes the upper-minus-lower scalar manifold
shift (labeled `ΔνF`) as an optional bottom panel and reports the vector component as the signed
fictitious magnetic field `B_fic` in µG. The Zeeman states by component view
separates vector and tensor state contributions into two plots. Eigenvalue
curves are labeled EV1, EV2, and so on.

Run from this directory:

```powershell
streamlit run alkali_pumping.py
```

## Dual-alkali behavior

- **Alkali A** is always active. **Alkali B** defaults to `None`.
- Selecting the same isotope for A and B leaves B inactive and ignores its
  pump settings.
- Selecting different isotopes enables a coupled A/B steady-state solve and
  separate result tabs.
- Density can use independent saturated-vapor curves or Raoult's law. In the
  relative-concentration mode, the input is the condensed-phase mole ratio
  `B/A`; each pure saturated-vapor density is multiplied by its corresponding
  liquid mole fraction.
- Self and unlike-alkali spin exchange are always included. The population
  result row reports both species' densities, electron-randomization rates,
  self-exchange rates, and directional cross-exchange rates together.
- PumpA1–PumpA3 and PumpB1–PumpB3 are configured in persistent sidebar tabs.
  Every active laser frequency is evaluated against both species.
- The optional rate-matrix display includes each local map and the full block
  population Jacobian `[[J_AA, J_AB], [J_BA, J_BB]]`.
- The sidebar contains the shared static-field direction and strength in nT.
- Each result tab has its own quantization axis, RF axis, and frequency range.
  Within that tab, the RF axis and lower/upper frequency bounds are shared by
  the atomic-moment and probe-readout plots; all other plot controls are
  independent.
- RF-A and RF-B have independent frequency ranges. In the A tab only RF-A is
  applied and the atomic observable and Probe-A readout belong to A; the B tab
  is defined conversely. New sessions default both RF upper sweep bounds to
  100 Hz.
- Probe-A and Probe-B have independent line, hyperfine-transition detuning,
  propagation direction, input ellipse, and path length. A Custom probe is a
  weak detector only and does not pump, shift, or broaden the atoms.
- Each probe source lists both `PumpN weak` and `PumpN`. The weak item inherits
  that pump's spectrum, direction, and polarization while ignoring its
  intensity. The plain pump item represents the physical pump beam and uses
  its intensity in a distributed, self-consistent Stokes calculation. New
  sessions default Probe-A to `PumpA2 weak` and Probe-B to `PumpB2 weak`.
- Physical-pump propagation includes rank-1 circular birefringence/dichroism
  and rank-2 linear birefringence/dichroism feedback into the vector and tensor
  light shifts. It retains the weak-RF approximation for the modulated atomic
  response. Orientation-induced and alignment-induced curves are nonlinear
  counterfactual solutions; the Total curve is the coupled physical solution
  and is not their sum.
- Probe signals include optical rotation, ellipticity, normalized Stokes
  s1/s2/s3, and fractional transmission. A single selector displays the
  orientation-induced, alignment-induced, or coherent total response; the
  total can optionally be decomposed. Its legend uses Total, Orientation, and
  Alignment columns with Amplitude, In phase, and Quadrature rows. The
  Mathur rank-2 polarizability is converted to the raw Cartesian Q_ij
  convention used by the atomic-moment plots.
- The Zeeman table is upstream of both field-response plots. A compact table
  lists upper-manifold transitions inside the shared RF sweep, and the same
  resonance markers appear on both plots.
- The dual-alkali RF solver retains coherence feedback through self and cross
  spin exchange even though the other RF drive is zero.

The condition-file schema is version 6.9. Complete v6.9 files are required;
v6.8, v6.7, v6.6, v6.5, v6.4, v6.3, v6.2, v6.1, v6.0, and v5.0 files are migrated automatically. Existing
v6.8 `PumpN` probe links migrate to `PumpN weak`, preserving their original
non-perturbing meaning. A relative-
concentration value loaded from an older file is retained numerically and is
interpreted as the liquid mole ratio `B/A` by the v6.4 Raoult-law model. Legacy Pump1–Pump3
become PumpA1–PumpA3. A legacy A Larmor frequency is converted to the
corresponding static field in nT.

The population model is diagonal in each selected quantization basis. When a
nonzero static-field direction is transverse to a tab's quantization axis, the
app warns that transverse static-field mixing is omitted and does not present
that RF curve. Choose the static-field direction as the quantization axis for
the supported secular weak-response calculation.

## Unlike-alkali spin-exchange defaults

| Pair | Cross section (cm²) |
| --- | ---: |
| Rb85–Rb87 | 1.70 × 10⁻¹⁴ |
| Rb–Cs | 2.30 × 10⁻¹⁴ |
| K–Rb | 2.00 × 10⁻¹⁴ |
| K–Cs | 2.24 × 10⁻¹⁴ |

The Rb-isotope value follows [Jarrett, *Phys. Rev.* 133, A111
(1964)](https://doi.org/10.1103/PhysRev.133.A111), and the Rb–Cs value follows
[Gibbs and Hull, *Phys. Rev.* 153, 132
(1967)](https://doi.org/10.1103/PhysRev.153.132). The K–Rb value is the
approximately 200 Å² hybrid-vapor convention; the K–Cs value converts the
approximately 800 a₀² thermal result from [Kartoshkin, *Optics and
Spectroscopy* 113, 235
(2012)](https://doi.org/10.1134/S0030400X12090081). Rates use the Maxwellian
mean relative speed and the collision partner's number density.

Native Streamlit theme settings can be added in `.streamlit/config.toml`.
See `CHANGELOG.md` for user-visible and physics-model updates.

## Layout

- `alkali_pumping.py`: top-level Streamlit navigation entry point.
- `alkali_pumping_app/pages/`: complete Light shift, Atomic polarizability, Magnetometry, and Atomic properties page implementations.
- `alkali_pumping_app/physics/`: single- and dual-species numerical models.
- `alkali_pumping_app/ui/`: condition-file and table-rendering helpers.
- `tests/`: regression and physical-consistency tests.
- `../archive/alkali_pumping_v6.7.13/`: application snapshot immediately before the independent dual-probe readout implementation.

Run the tests with:

```powershell
python -m unittest discover -s tests -v
```
