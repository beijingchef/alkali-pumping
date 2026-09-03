# Changelog

## 6.9.0 - 2026-09-02

### Added

- Added an explicit weak source beside every pump source in the Probe source
  selectors. `PumpN weak` uses the existing non-perturbing linear probe solver;
  `PumpN` treats the selected pump as the physical detected beam.
- Added distributed weak-RF, full-Stokes feedback for physical-pump readout.
  The nonlinear propagation includes CBOR and LDOR contributions, the linked
  pump's ellipticity, and every existing probe signal and RF quadrature.
- Added nonlinear orientation-only and alignment-only counterfactual responses
  alongside the coupled Total response, plus readout-mode metadata in CSV
  exports.

### Changed

- Advanced the condition schema to v6.9. Legacy v6.8 pump-linked probes migrate
  to their explicit `weak` variants so loading an old file does not silently
  change its physics.

## 6.8.9 - 2026-09-01

### Changed

- Restored native Streamlit field heights and vertical spacing in all condition
  sidebars by removing the custom compact-sizing CSS. The Probe source, Beam
  direction, and Path length row introduced in v6.8.8 remains unchanged.

## 6.8.8 - 2026-09-01

### Changed

- Renamed Probe beam source to Probe source and placed Probe source, Beam
  direction, and Path length on one row.
- Added compact field heights and vertical gaps across the condition sidebars;
  this styling was rolled back in v6.8.9 in favor of native Streamlit sizing.

## 6.8.7 - 2026-09-01

### Changed

- Made the mixture density model and molar-ratio controls equal width.
- Placed the probe reference-line and hyperfine-transition selectors on one
  row, and moved the D1/D2-center detuning caption to the right of the
  detuning input.

## 6.8.6 - 2026-09-01

### Changed

- Changed new-session Probe-A and Probe-B sources to PumpA2 and PumpB2.
- Changed the default RF-A and RF-B upper sweep frequencies from 50 Hz to
  100 Hz. Existing saved conditions retain their stored frequency bounds and
  legacy probe configurations.

## 6.8.5 - 2026-09-01

### Changed

- Arranged the decomposed optical-readout legend as a three-column grid for
  Total, Orientation, and Alignment, with Amplitude, In phase, and Quadrature
  aligned by row.

## 6.8.4 - 2026-09-01

### Changed

- Replaced the independent probe orientation/alignment switches with one
  selector for Orientation induced, Alignment induced, or their coherent
  Total. Existing v6.7 conditions migrate their previous switch combination.
- Clarified that the Atomic polarizability page displays Mathur's normalized
  birefringent coefficient directly; only coupling it to the raw Cartesian
  Q_ij observable requires the sign-and-normalization conversion.

## 6.8.3 - 2026-09-01

### Fixed

- Corrected the weak-probe alignment readout by converting Mathur's normalized
  rank-2 polarizability operator to the raw Cartesian Q_ij convention used by
  the magnetometry solver. This fixes both the LDOR sign and its normalization
  for every supported hyperfine manifold, so CBOR and LDOR combine with the
  physically correct relative phase.

## 6.8.2 - 2026-08-31

### Changed

- Changed weak-probe field-response plots and CSV exports from response per RF
  Rabi angular frequency to response per RF magnetic-field amplitude. Rotation
  and ellipticity are now shown in rad/nT; normalized Stokes parameters and
  fractional transmission are shown in nT⁻¹.
- Removed the redundant alkali-channel suffix from the probe-response vertical
  axis quantity.
- Preserved independent Probe-A and Probe-B configuration values when switching
  between their sidebar tabs.

## 6.8.1 - 2026-08-31

### Added

- Added a probe source selector. Probe-A can inherit PumpA1–PumpA3 and Probe-B
  can inherit PumpB1–PumpB3, copying spectral, direction, and polarization
  settings while ignoring intensity and retaining an independent path length.

### Changed

- Swapped the displayed order of probe path length and ellipticity angle.
- Removed the redundant upper-manifold-transition table above the RF plots;
  resonance markers remain available on both plots.
- Advanced Magnetometry condition files to schema v6.7 with automatic v6.6
  migration to the Custom probe-source setting.

## 6.8.0 - 2026-08-31

### Added

- Added independent weak Probe-A and Probe-B optical readouts for the two
  alkalis, with line/transition detuning, direction, polarization ellipse, and
  path-length settings.
- Added optical rotation, ellipticity, normalized Stokes s1/s2/s3, and
  fractional-transmission RF spectra derived from the Mathur scalar,
  gyrotropic, and birefringent polarizabilities.
- Added independent scalar, orientation, and alignment contribution switches,
  optional rank-resolved curves, X/Y/amplitude controls, and probe CSV exports.
- Added a Zeeman-transition-in-sweep summary and matching resonance markers to
  both atomic-moment and optical-readout plots.

### Changed

- Moved each alkali's RF axis and lower/upper frequency controls to a shared
  field-response header used by both plots. All other atomic and probe plot
  settings remain independent.
- Advanced the condition schema to v6.6; v6.5 files acquire independent probe
  defaults during migration.
- Advanced the physics model to the coupled independent-dual-probe Stokes
  readout model. Probes are non-perturbing and RF-A/Probe-A and RF-B/Probe-B
  remain separate diagonal channels while collision feedback stays coupled.

## Unreleased

### Added

- Added the cyclic off-diagonal alignment observable Q_ij to Magnetometry
  weak-RF response plots: RF x reads Q_yz, RF y reads Q_zx, and RF z reads
  Q_xy.
- RF susceptibility labels use Ω_rf, and off-diagonal Q components render
  their tensor indices as subscripts.
- Changing an alkali isotope now translates each pump's hyperfine transition
  by relative ground/excited manifold rank, without storing per-isotope
  transition history.
- RF in-phase X and quadrature Y are no longer sign-flipped automatically.
  Independent Add π controls now flip each displayed curve and its legend
  label when requested; condition files advance to schema v6.5.
- Added **Analysis → Atomic polarizability**, with D1/D2 complex response plots
  for Mathur's equilibrium, hyperfine, gyrotropic, and birefringent
  polarizabilities. The page includes atom/cell, sweep, and per-component plot
  controls; gyrotropic and birefringent plots are enabled by default.
- Added versioned JSON condition load/save controls to Atomic polarizability.
  Temperature and N₂ pressure share one row, the fixed pressure-coefficient
  editor and Light heading are removed, and all four subscripted alpha plot
  checkboxes share one row.
- Atomic-polarizability charts now color the upper ground-state manifold
  F=I+1/2 orange and use a controlled 2% horizontal margin around the requested
  sweep bounds instead of Altair's larger automatic axis expansion.
- Removed the individual Atomic-polarizability plot headings, replaced the
  generic y-axis title with component-specific polarizability labels, and
  increased plot height by 30%.
- Atomic-polarizability sweep defaults now extend 1500 MHz below and above the
  lowest and highest allowed transitions, rounded outward to 100 MHz. Removed
  the Y-scale selector.
- Page-specific settings now survive navigation between all analysis and
  reference pages for the lifetime of the browser session.

### Removed

- Removed the obsolete v4.23 source-snapshot regression test and remaining
  active-app metadata/documentation references to that deleted snapshot.

## 6.7.13 - 2026-08-26

### Changed

- Promoted the complete Light shift implementation from `ui/light_shift.py`
  into `pages/light_shift.py`.
- Promoted the Atomic properties implementation from `ui/atomic_settings.py`
  into `pages/atomic_properties.py` and removed the obsolete dialog entry point.
- Archived the previous application as `archive/alkali_pumping_v6.7.12`.

## 6.7.12 - 2026-08-18

### Changed

- Changed the browser tab title to `Optical pumping: <page name>`, using the
  active navigation page name.
- Archived the previous application as `archive/alkali_pumping_v6.7.11`.

## 6.7.11 - 2026-08-18

### Changed

- Placed the Light shift Upload button, Save button, and Settings name input on
  one sidebar row with relative widths 30:25:45.
- Archived the previous application as `archive/alkali_pumping_v6.7.10`.

## 6.7.10 - 2026-08-18

### Changed

- Added a saved **Show scalar shift** toggle to the Light shift Components view;
  it defaults to enabled.
- Moved the optional scalar plot below the other plots, including the optional
  scattering-rate plot.
- Updated light-shift settings files to version 1.1 and automatically migrates
  version 1.0 files with the scalar plot enabled.
- Archived the previous application as `archive/alkali_pumping_v6.7.9`.

## 6.7.9 - 2026-08-18

### Changed

- When the Light shift polarization input changes from **Preset** to
  **Ellipse**, initialize the ellipse azimuth and ellipticity to reproduce the
  currently selected preset for any beam direction.
- Archived the previous application as `archive/alkali_pumping_v6.7.8`.

## 6.7.8 - 2026-08-18

### Changed

- Shortened the Magnetometry sidebar label from `Static field strength (nT)`
  to `Strength (nT)`.
- Darkened the dashed transition-center markers across all Light shift plots.
- Added the transition-center markers to the scattering-rate plot and removed
  its `Hyperfine manifold` legend title.
- Archived the previous application as `archive/alkali_pumping_v6.7.7`.

## 6.7.7 - 2026-08-18

### Changed

- Shortened the Magnetometry sidebar detuning caption from
  `from D1/D2 line center` to `from D1/D2 center`.
- Archived the previous application as `archive/alkali_pumping_v6.7.6`.

## 6.7.6 - 2026-08-18

### Changed

- Simplified the Components scalar-difference legend entry to plain `ΔνF`.
- Archived the previous application as `archive/alkali_pumping_v6.7.5`.

## 6.7.5 - 2026-08-18

### Fixed

- Restored the Components scalar-shift legend in the normal chart view by
  using Vega-Lite's native legend instead of a clipped concatenated panel.
- Kept `F` visibly lowered beneath `Δν` in the scalar-difference entry.
- Archived the previous application as `archive/alkali_pumping_v6.7.4`.

## 6.7.4 - 2026-08-18

### Changed

- Made **Analysis → Light shift** the default startup page and root route.
- Typeset the Components scalar-difference legend as `Δν` with a visually
  lowered, smaller `F` subscript instead of displaying a literal underscore.
- Archived the previous application as `archive/alkali_pumping_v6.7.3`.

## 6.7.3 - 2026-08-18

### Changed

- Relabeled the Components scalar-difference curve from the explicit manifold
  expression to `Δν_F`; its value remains the upper-manifold scalar shift minus
  the lower-manifold scalar shift.
- Archived the previous application as `archive/alkali_pumping_v6.7.2`.

## 6.7.2 - 2026-08-18

### Changed

- Removed the `Hyperfine manifold` title from Components legends while keeping
  the individual manifold entries.
- Removed the Scalar panel from **Zeeman states by component**, leaving separate
  Vector and Tensor state plots.
- Renamed the Eigenvalues legend title from `Eigenvalue branch` to `Eigenvalues`
  and shortened branch labels from `Branch 1`, `Branch 2`, ... to `EV1`, `EV2`, ....
- Archived the previous application as `archive/alkali_pumping_v6.7.1`.

## 6.7.1 - 2026-08-18

### Fixed

- Restored the Components vector plot's y-axis title in the normal app view by
  shortening the title and explicitly reserving axis-title padding and limits.

### Changed

- Renamed the vector y-axis quantity from `Fictitious magnetic field` to
  `B_fic`, with units of µG or µG/(µW/cm²).
- Archived the previous application as `archive/alkali_pumping_v6.7.0`.

## 6.7.0 - 2026-08-18

### Added

- Added the upper-ground-manifold scalar shift minus the lower-manifold scalar
  shift to the Components scalar panel (`F=2 − F=1` for Rb87 and K39).
- Added **Zeeman states by component**, which renders separate Scalar, Vector,
  and Tensor panels with one curve per selected `|F,m⟩` state.

### Changed

- Replaced the Components vector-coefficient plot with the signed fictitious
  magnetic field, `B_fict,F = V_F/γ_F`, in µG or µG/(µW/cm²). Each manifold's
  own signed gyromagnetic ratio is used.
- Added an explicit unit column to Components CSV exports because their vector
  and scalar/tensor rows now use different physical units.
- Archived the previous application as `archive/alkali_pumping_v6.6.0`.

## 6.6.0 - 2026-08-17

### Added

- Added grouped top navigation with **Analysis → Light shift**, **Analysis →
  Magnetometry**, and **Reference → Atomic properties**.
- Added an independent Light shift page whose sidebar contains atom, cell,
  static-field, pressure-coefficient, beam, arbitrary-polarization, sweep, and
  display controls while the main panel contains the plots and CSV export.
- Added dedicated versioned JSON save/load support for Light shift settings.
- Added page-safe Light shift state so its settings survive navigation to
  another page and back.

### Changed

- Moved the existing dual-alkali analysis unchanged to the Magnetometry page
  and removed its embedded light-shift explorer.
- Promoted the Atomic properties dialog contents to a full reference page.
- Archived the previous application as `archive/alkali_pumping_v6.5.0`.

## 6.5.0 - 2026-08-16

### Fixed

- Kept the light-shift charts mounted when changing explorer controls. The
  cached explorer is rendered in its dynamic expander's normal app context
  instead of as a nested fragment, which previously cleared the Components
  chart subtree on an isolated widget rerun.
- Increased the normal Components-chart height and line weight so narrow
  dispersive structure is visible without entering fullscreen mode.

### Added

- Added a lazy light-shift detuning explorer to every active alkali result.
- Added scalar, vector, tensor-m=0, state-resolved, adjacent-transition, and
  within-manifold eigenvalue views with CSV export and optional scattering-rate
  curves. Adjacent-transition shifts can also be displayed as equivalent
  light-shift magnetic fields.
- Added arbitrary pure polarization through beam-frame ellipse azimuth and
  signed ellipticity, alongside the existing polarization presets.
- Added complex signed dipole amplitudes and a Hermitian, symmetrized
  light-shift Hamiltonian. Its diagonal retains exact agreement with the
  existing complex-Voigt light shifts, while its per-F blocks include coherent
  Raman terms for arbitrary polarization.
- Added explicit spherical-polarization fractions, Stokes parameters, and the
  dimensionless Mathur `E^(2)_0` geometry factor. Tensor plot units are kept as
  frequency per intensity rather than treating `E^(2)_0` as a unit.
- Added regression coverage for polarization special cases, Hamiltonian
  Hermiticity, component reconstruction, legacy-diagonal agreement, and every
  light-shift export table.

### Changed

- Extended the complex Voigt response to accept vectorized detuning arrays for
  efficient cached sweeps.
- Bumped the physics model to the v6.5 arbitrary-polarization light-shift model;
  the condition-file schema remains v6.4 because explorer display settings are
  session-local and existing condition files remain fully compatible.
- Archived the previous application as `archive/alkali_pumping_v6.4.12`.

## 6.4.12 - 2026-08-16

### Changed

- Formatted the pump line-center detuning caption with no decimal places.
- Renamed the state-specific pump-rate caption to the simpler `pump rate`.
- Archived the previous application as `archive/alkali_pumping_v6.4.11`.

## 6.4.11 - 2026-08-16

### Changed

- Shortened the pump line-center detuning caption from `Detuning` to `δν` and
  formatted its value with one digit after the decimal point.
- Archived the previous application as `archive/alkali_pumping_v6.4.10`.

## 6.4.10 - 2026-08-16

### Changed

- Added a live caption below every pump detuning input showing the detuning
  from the selected zero-pressure D1/D2 fine-structure line center used by
  the optical-pumping solver.
- Archived the previous application as `archive/alkali_pumping_v6.4.9`.

## 6.4.9 - 2026-08-11

### Changed

- When Alkali B is `None`, gray only its selector label and remove the
  redundant inactive caption, including inactive pump-rate captions.
- Display active informational captions in blue and the duplicate-alkali
  ignored notice in orange, reserving gray for inactive inputs.
- Archived the previous application as `archive/alkali_pumping_v6.4.8`.

## 6.4.8 - 2026-08-11

### Changed

- Plotted both hyperfine population manifolds on one shared numerical `m`
  axis, so bars with the same magnetic quantum number align vertically.
- Added aligned-axis regression coverage for the population plot.
- Continued from the restored v6.4.6 baseline; the withdrawn v6.4.7 theme
  experiment remains archived separately.

## 6.4.6 - 2026-08-11

### Changed

- Restored the adjacent-coherence optical-pumping relaxation column
  `Γ^{OP} (s^-1)` immediately after `G^{OP} (s^-1)` in the Zeeman table;
  the separate `Γ^{OP}/2π` column remains omitted.
- Rendered the `self` and `cross` labels as horizontally aligned subscripts
  beneath the `SE` superscript in the spin-exchange rate headings.
- Archived the previous application as `archive/alkali_pumping_v6.4.5`.

## 6.4.5 - 2026-08-10

### Changed

- Display all population- and coherence-relaxation rate columns (`G` and `Γ`)
  in the Zeeman table with one digit after the decimal point.
- Archived the previous application as `archive/alkali_pumping_v6.4.4`.

## 6.4.4 - 2026-08-10

### Changed

- Changed the population-graph and Zeeman-table result columns to an exact
  25%/75% width split (`1:3`).
- Archived the previous application as `archive/alkali_pumping_v6.4.3`.

## 6.4.3 - 2026-08-10

### Changed

- Rendered scientific notation throughout the population/rate caption as
  compact math such as `5.5×10^9` instead of `5.5e+09`.
- Replaced middle-dot separators between adjacent caption variables with
  commas.
- Preserved the archived v6.4.2 application.

## 6.4.2 - 2026-08-10

### Changed

- Formatted `⟨m⟩`, `⟨m²⟩`, `n_A`, and `n_B` in the result caption with two
  significant figures. Tiny numerical moment residuals are displayed as zero;
  collision and relaxation rate precision is unchanged.
- Preserved the existing `archive/alkali_pumping_v6.4.1` snapshot.

## 6.4.1 - 2026-08-10

### Changed

- In relative-concentration mode, the result caption now displays each actual
  vapor density together with its Raoult factor, for example
  `n_A = 5.51e9 cm^-3 = 0.500 n_A^sat`, making the reduction from the pure
  saturated-vapor density explicit.
- Archived the previous application as `archive/alkali_pumping_v6.4.0`.

## 6.4.0 - 2026-08-10

### Changed

- Replaced the former relative vapor-density rule with Raoult's law. For an
  entered condensed-phase mole ratio `r = B/A`, the model now uses
  `x_A = 1/(1+r)` and `x_B = r/(1+r)`, then calculates
  `n_A = x_A n_A^sat(T)` and `n_B = x_B n_B^sat(T)`.
- Relabeled the sidebar input as `Liquid mole B/A` and added help explaining
  its condensed-phase meaning.
- Added v6.3 condition migration. Existing numeric ratio values are retained
  but use the new liquid-mole-ratio interpretation.
- Archived the previous application as `archive/alkali_pumping_v6.3.0`.

## 6.3.0 - 2026-08-10

### Changed

- Spin exchange is now always included in both the one- and two-alkali
  solvers; the sidebar checkbox and its condition-file field were removed.
- Moved `nA`, `nB`, `R_A←B`, and `R_B←A` out of the sidebar. Each population
  result now shows one consolidated line below the population graph and
  Zeeman table containing `⟨m⟩`, `⟨m²⟩`, both densities, both electron-
  randomization rates, both self-exchange rates, and both directional
  cross-exchange rates.
- Added migration for v6.2 conditions, discarding the obsolete spin-exchange
  toggle, and archived the previous v6.2.0 source.

## 6.2.0 - 2026-08-10

### Added

- Added PumpA3 and PumpB3 to the isolated sidebar pump tabs. Both new pumps
  default to zero intensity and use the same fragment-only rerun behavior as
  the existing pumps.
- Added PumpA3 and PumpB3 settings to v6.2 condition files and their optical
  frequencies to the transition table.

### Changed

- The condition Save action now synchronizes its persistent sanitized filename
  with the visible condition-name field through both name-change and save-click
  callbacks.
- RF signed-component curves continue to plot `-X` and `-Y` and now use those
  exact labels in the legend.
- Added automatic migration of v6.1 condition files with zero-intensity A3/B3
  defaults and restored legacy v5 Beam 3 as PumpA3.
- Archived the previous 6.1.7 source without copying `.venv` or Python caches.

## 6.1.7 - 2026-08-08

### Fixed

- Restricted each alkali's reported light shift and light-shift diagonality
  check to pumps targeting that alkali. A nonzero Pump A can no longer blank
  Alkali B's light-shift columns, and the converse is also true.
- Preserved off-resonant optical-pumping contributions and all self/cross spin
  exchange physics; only cross-alkali AC-Stark shifts are neglected.
- Archived the previous 6.1.6 source without copying the incomplete `.venv`.

## 6.1.6 - 2026-08-08

### Fixed

- Made Alkali A/B result tabs lazy and strictly isolated. The selected tab now
  contains its own quantization-axis control, population graph, Zeeman table,
  RF settings, and RF susceptibility plot as one complete result section.
- Added persistent backing state for all lazy result-tab controls so switching
  tabs cannot reset or mix A/B quantization and RF configurations.
- Preserved the selected result tab through quantization-axis and RF-control
  reruns, including when the selected isotope changes.
- Archived the previous 6.1.5 source without copying the incomplete `.venv`.

## 6.1.5 - 2026-08-08

### Fixed

- Moved the pending full-rerun check to the end of the pump fragment so an
  intensity change reliably updates the physical solution and rate caption.
- Confirmed that, after PumpA1 or PumpA2 reaches zero intensity, direction,
  polarization, line, transition, and detuning edits remain fragment-only and
  do not rerun the full physical system.
- Archived the previous 6.1.4 source without copying the incomplete `.venv`.

## 6.1.4 - 2026-08-08

### Fixed

- Removed the `Calling st.rerun() within a callback is a no-op` warning from
  pump controls. Callbacks now set a pending full-rerun flag, and the pump
  fragment consumes that flag before requesting the supported app rerun.
- Preserved fragment-only reruns for zero-intensity beam settings and full
  physical recomputation for intensity changes or active-beam edits.
- Archived the previous 6.1.3 source without copying the incomplete `.venv`.

## 6.1.3 - 2026-08-08

### Fixed

- Kept the selected Alkali A/B result tab active when a quantization-axis or
  RF control triggers a rerun, so Alkali B controls no longer appear to alter
  or jump back to Alkali A.
- Isolated zero-intensity pump configuration edits in a sidebar fragment.
  Direction, polarization, line, transition, and detuning changes now avoid a
  full app rerun until that beam has nonzero intensity.
- Cached physical-system solutions and excluded zero-intensity beams from the
  solver input, while retaining their stored UI configuration.
- Archived the previous 6.1.2 source without copying the incomplete `.venv`.

## 6.1.2 - 2026-08-08

### Fixed

- Preserved every Alkali A and Alkali B pump setting when switching between
  the lazy pump-configuration tabs. Visible widget values are now copied into
  persistent condition state instead of being lost during hidden-widget
  cleanup.
- Kept condition-file loading and saving connected to the persistent pump
  settings while using separate temporary keys for visible tab controls.
- Archived the previous 6.1.1 application before applying this update.

## 6.1.1 - 2026-08-08

### Fixed

- Kept the Alkali A and Alkali B pump configurations in persistent keyed tabs
  when a pump transition or another tab-local setting triggers a rerun.
- Rendered only the open pump tab so the two pump configurations cannot appear
  stacked in the sidebar after a rerun.

### Changed

- Moved the `n(B) / n(A)` input onto the same sidebar row as the mixture density
  model selector.
- Archived the previous 6.1.0 application before applying this update.

## 6.1.0 - 2026-08-08

### Added

- Added a shared static-field direction and signed field strength in nT.
- Added independent Alkali A and Alkali B quantization-axis controls above
  their population and Zeeman-result regions.
- Added independent RF-A and RF-B axes, observables, frequency ranges, curve
  selections, and normalization settings beside their susceptibility plots.
- Added a coupled coherence-response generator. Each result applies only its
  own RF drive while retaining self- and cross-species spin-exchange feedback.
- Added automatic v6.0 condition migration, including conversion of the old
  A upper-manifold Larmor frequency to static-field strength.

### Changed

- Moved all RF-response controls out of the sidebar and to the left of the
  corresponding Alkali A or Alkali B susceptibility plot.
- Archived the pre-update application in `archive/alkali_pumping_v6.0.0`.

## 6.0.0 - 2026-08-07

### Added

- Added optional **Alkali B** selection beside **Alkali A**. `None` is the
  default, and a B selection identical to A is intentionally inactive.
- Added independent saturated-vapor and relative-concentration density modes.
- Added persistent sidebar pump tabs with PumpA1, PumpA2, PumpB1, and PumpB2.
- Added coupled unlike-alkali spin exchange, including two-species fixed-point
  populations and the full block small-signal population Jacobian.
- Added separate Alkali A and Alkali B result tabs when B is active.
- Added self- and cross-spin-exchange contributions to the Zeeman table.
- Added v5.0 condition migration to the v6.0 condition schema.

### Changed

- Retired the third A pump. PumpA1 and PumpA2 default to 5.0 µW/cm²; PumpB1
  and PumpB2 default to 0.0 µW/cm².
- Every active pump is evaluated at its absolute optical frequency for both
  species, including isotope cross-pumping.
- Rate-matrix output now includes the coupled A/B Jacobian and its local maps.
- Archived the pre-upgrade application in `archive/alkali_pumping_v5.2.16`.

## 5.2.16 - 2026-07-23

### Added

- Added a beam-intensity input in µW/cm² for each of the three pump beams.
- Added calculated sidebar captions for the selected ground manifold's total
  pump rate at the reference resonance and at the specified detuning.
- Added automatic conversion of legacy v5 rate-referenced condition files to
  the equivalent physical beam intensities.

### Changed

- Replaced the pump-rate reference selector and pump-rate input with an
  absolute weak-light rate calculation based on photon flux, D-line
  wavelength, natural linewidth, and pressure/Doppler broadening.

## 5.2.15 - 2026-07-21

### Added

- Added a fully commented `.streamlit/config.toml` template for future native
  Streamlit theme, font, color, border, and sidebar customization.

### Changed

- Documented the theme configuration file in the README.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.14`.

## 5.2.14 - 2026-07-21

### Changed

- Moved the calculated `R_SE` caption directly below the **Include spin
  exchange** checkbox.
- Removed alkali density and spin-exchange cross-section values from that
  sidebar caption.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.13`.

## 5.2.13 - 2026-07-21

### Changed

- Removed all injected sidebar widget CSS, including custom heights, fonts,
  file-uploader sizing, and number-input step-button sizing.
- Restored Streamlit's native appearance and dimensions for every sidebar
  input, selection, checkbox, uploader, and button.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.12`.

## 5.2.12 - 2026-07-21

### Changed

- Removed the editable N₂ pressure-broadening and shift coefficient section
  from the sidebar while retaining the stored coefficients in calculations and
  condition files.
- Set every sidebar input field and selection box to a uniform 25-pixel height.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.11`.

## 5.2.11 - 2026-07-21

### Changed

- Made Streamlit's native decrease and increase controls visible for N₂
  pressure, temperature, and every beam's detuning and pump-rate fields.
- Reduced the requested native step-button width and input height for a compact
  sidebar appearance.
- Moved N₂ pressure and temperature into a two-column row so each field is
  wide enough for Streamlit to render its native controls.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.10`.

## 5.2.10 - 2026-07-21

### Changed

- Changed the fresh-app pump-rate reference default to **At resonance** for
  all three beams.
- Updated the bundled default condition accordingly; explicitly loaded saved
  conditions continue to retain their own pump-rate-reference choices.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.9`.

## 5.2.9 - 2026-07-21

### Added

- Added vector light shift $\nu^{\mathrm{VS}}$ and tensor light shift
  $\nu^{\mathrm{TS}}$ columns immediately before the total light shift in the
  Zeeman-sublevel properties table and its CSV export.
- Decomposed each diagonal total AC-Stark shift within its hyperfine manifold
  into scalar, vector, and tensor state contributions.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.8`.

## 5.2.8 - 2026-07-21

### Changed

- Shortened the Grotrian hyperfine-level labels to `F=…` and `F′=…`.
- Right-aligned each label with a fixed gap to the left of its nearest Zeeman
  level segment, preventing overlap with the level-denoting line.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.7`.

## 5.2.7 - 2026-07-21

### Changed

- Increased the Grotrian stacked-fraction annotation font from 10 to 15
  points, exactly 50%, while leaving magnetic-quantum-number labels unchanged.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.6`.

## 5.2.6 - 2026-07-21

### Changed

- Replaced linear transition-strength labels such as `1/12` with vertically
  stacked MathText fractions.
- Applied the same vertical fraction format to the summed displayed strength.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.5`.

## 5.2.5 - 2026-07-21

### Changed

- Lowered the common σ−/σ+ transition-strength label level from 0.40 to 0.30
  while retaining exact placement on each slanted transition line.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.4`.

## 5.2.4 - 2026-07-21

### Changed

- Raised all π-transition strength labels to a common upper level.
- Placed σ− and σ+ strength labels at the same lower level and calculated each
  horizontal label coordinate from its transition line, keeping every label
  centered on the corresponding slanted line.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.3`.

## 5.2.3 - 2026-07-21

### Changed

- Increased Grotrian fractional-strength annotations to the same 10-point
  font size as the magnetic-quantum-number labels.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.2`.

## 5.2.2 - 2026-07-21

### Changed

- Moved every Grotrian-diagram control into a vertical panel to the left of
  the graph.
- Made every polarization and display option checked by default, including
  transition-strength labels.
- Changed individual and summed transition-strength labels from decimals to
  reduced fractions.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.1`.

## 5.2.1 - 2026-07-21

### Fixed

- Moved the atomic-properties temperature slider inside the first tab.
- Decoupled the buffer-gas tab from that slider by reporting its calculated
  electron-randomization rates at a stated 20 °C reference temperature.
- Made the modal tabs stateful and render only the active tab, so a
  Grotrian-control rerun preserves the third tab without replaying the first
  two tabs above it.
- Archived the pre-update application in `archive/alkali_pumping_v5.2.0`.

## 5.2.0 - 2026-07-21

### Added

- Added a **Settings** button that opens a large modal atomic-properties dialog.
- Added isotope selection for ²³Na, ³⁹K, ⁴¹K, ⁸⁵Rb, ⁸⁷Rb, and ¹³³Cs.
- Added temperature-dependent saturated-vapor density, RMS velocity, mean
  relative velocity, and self spin-exchange rate calculations.
- Added optical pressure broadening and shift tables for N₂, ⁴He, and CH₄,
  together with editable ground-state electron-randomization cross sections
  and calculated collision rates.
- Added an interactive hyperfine–Zeeman Grotrian diagram whose selection rules,
  colors, and line strengths reproduce the supplied Mathematica notebook.
- Archived the pre-update application in `archive/alkali_pumping_v5.1.0`.

- Added an optional **Density factor** for the weak-RF susceptibility plot.
  It multiplies every plotted component by the calculated saturated alkali
  vapor density in cm⁻³.
- Added density-factor status, density, and resulting plotted units to the
  weak-RF CSV export.
- Added backward-compatible loading of v5.0 condition files that predate the
  density-factor field; the new option defaults to off for those files.

### Changed

- Replaced the RF relaxation-normalization caption with an always-visible
  summary of scientific sidebar inputs and active pump beams.
- Shortened the pump-rate input label to **Pump rate** while retaining the
  total Zeeman-summed definition introduced in v5.1.

## 5.1.0 - 2026-07-20

### Changed

- Changed the pump-rate input from an average over ground Zeeman sublevels to
  the total selected-transition rate summed over ground and excited Zeeman
  sublevels.
- Updated the pump-rate documentation, normalization tests, and application
  metadata for the new definition.
