# Changelog

## 6.20.14 - 2026-10-28

- Restore Time evolution's separate Open and Save buttons beside Setting name,
  replacing the File dropdown. Remove Settings headings from all five sidebars.
- Make settings Save a browser-only download on both analysis and reference
  pages: no widget submission, app/fragment rerun or calculation. Download the
  committed configuration and read only the currently typed setting name.
  Prevent mouse-down blur so other unfinished inputs are not applied by Save.
- Keep ACC's download snapshot current when inactive-pump controls commit in
  their own fragment. Preserve existing filename defaults and condition formats.
- Validation: all 245 tests pass (169.880 seconds). Browser downloads on both
  analysis pages use an unfinished filename while leaving unfinished physical
  inputs focused and unapplied; downloaded values retain the committed settings
  and Save does not refresh its snapshot. Verify committed ACC fragment edits,
  restored buttons and removal of Settings titles across all five pages.
  Git diff whitespace checks pass.

## 6.20.13 - 2026-10-28

- Put Time evolution's Open and Save controls in a native File dropdown at
  the top of its sidebar, with the setting name beside it. Retain the existing
  pending-edit Save download, validated Open callback and load feedback.
  Opening the dropdown does not rerun calculations. Saved settings remain
  format 2.5; other pages retain their file controls.
- Reset the uploader after each load attempt so File retains its Open action
  and users can select the same file again. Keep validation feedback outside
  the menu and leave saved condition fields unchanged.
- Validation: all 244 tests pass (179.333 seconds), including pending-edit Save,
  File menu placement, and repeated valid/invalid Open callbacks. Browser checks
  verify the two menu actions, load a current-format file, and validate the actual
  Save download against its loaded values. Git diff whitespace checks pass.

## 6.20.12 - 2026-10-28

- Halve the visible gap between navigation and main content across all pages,
  retaining clearance for the native header. With the default font, desktop
  clearance changes from 68 to 34 px; main top padding changes from 128 to 94 px.
- Halve the pinned bottom container's lower padding from 3.5rem to 1.75rem
  (56 to 28 px), bringing Auto update and Apply closer to the browser edge.
  Use shared screen-only styles; native toolbar pinning and print spacing remain.
- Halve the sidebar's upper space from 76 to 38 px at the default font by
  reducing its native header height and following margin, keeping the collapse
  button accessible. Apply the same spacing rules across pages.
- Validation: all 243 tests pass (177.709 seconds). Browser measurements confirm
  34 px below navigation, 38 px above sidebar content, and 28 px below the
  toolbar controls. Verify shared spacing on Time evolution and Atomic
  properties; Git diff whitespace checks pass.

## 6.20.11 - 2026-10-28

- Move Time evolution Auto update and Apply into a native pinned bottom toolbar,
  keeping both controls visible while scrolling plots. Preserve batching,
  pending-result handling and page-session settings; page help stays beside
  the title. No saved-condition fields or numerical calculations change.
- Validation: all 243 tests pass (153.404 seconds), including the real-page
  Apply workflow, control placement and navigation persistence. Git diff
  whitespace checks pass.

## 6.20.10 - 2026-10-28

- Prepare constant density generators, state boundaries, electron-randomization
  and uniform-relaxation maps once per calculation. Supply transient and
  stationary BDF integration with an analytic complex Jacobian, including
  normalized electron marginals and both cross-species derivative blocks.
- Defer transient integration when every time-domain section is closed.
  Opening only Probe signal scanning performs its stationary calculation
  independently. Share one transient retrieval between species per rerun.
- Add Auto update (initially checked) and Apply controls for batching physical
  edits. Hide affected curves and downloads while changes are pending. Applying
  with closed sections records the inputs without integrating. Execution controls
  are session preferences outside saved-condition files.
- Bound transient caches separately to two trajectories below 8 MiB and one
  large trajectory. Cache four population/moment display views without retaining
  their full rotated density trajectories; reduce stationary/probe caches to
  four entries each and spectra to eight. Evicted inputs are recomputed.

- Select the Time evolution transient integrator independently for pump-on
  and dark intervals. Retain BDF for damping-dominated/short intervals and
  choose DOP853 for sustained coherent evolution without faster damping,
  considering magnetic precession, AC-Stark shifts, optical damping and
  both species' collision/relaxation rates.
- Preserve the master equation, output times, tolerances, and trace,
  Hermiticity and positivity checks. Retry failed/nonfinite DOP853 intervals
  from their original state with BDF; preserve errors if the retry fails.
- Include a separate numerical-solver version in transient and stationary cache keys.
  Settings format 2.5, the physical-model identifier and stationary probe-scan
  convergence criteria remain unchanged.
- Local four-second Rb87 timings (single runs, original BDF versus complete
  optimization): default zero field 0.447 → 0.337 s; transverse 5 nT
  8.165 → 1.165 s (7.0×); tenfold pump at zero field 0.474 → 0.362 s;
  tenfold pump with transverse 5 nT 6.951 → 2.912 s. Prepared BDF alone takes
  3.844 s in the transverse example (2.1×). Maximum density differences from
  original BDF are below 2.3×10⁻⁹; tolerances and physical checks are unchanged.
  These are development-machine measurements, not hosted timing guarantees.
- Validation: all 243 tests pass with no failures or skips (155.209 seconds),
  including finite-difference Jacobian comparisons, tighter BDF and exact linear
  references, physical-state invariants, bounded-cache eviction/version checks,
  and real-page closed-section/manual-Apply/navigation/export checks.
  Git diff whitespace checks pass.

## 6.20.9 - 2026-10-27

- Order Time evolution moment and probe plots by their selection lists.
  Defaults follow checkbox order: hyperfine polarization then ⟨F_z⟩, and
  Absorption then Optical rotation. Each newly checked plot appends at the end;
  unchecking removes it without reordering the others, and rechecking appends it.
- Preserve the order across sidebar edits, section closing/reopening,
  navigation and settings save/load. Section CSV columns follow plot order.
  Moment order remains shared by both species; probe order stays independent.
- Reuse the existing ordered selection lists in settings format 2.5.
  No saved fields, defaults, validation rules or physics calculations change.
- Validation: the focused plot-order and Save UI checks pass. The complete
  suite passes all 228 tests with no failures or skips (147.459 seconds).
  Git diff whitespace checks pass.

## 6.20.8 - 2026-10-27

- Render Time evolution plots, section controls and CSVs only while their
  sections are open, including after sidebar edits. Zeeman populations now
  starts collapsed; Atomic moments and Probe signal scanning remain collapsed.
- Default atomic-moment selection to hyperfine polarization and ⟨F_z⟩.
- Add a row of Absorption, Optical rotation and Ellipticity checkboxes to each
  species' Probe signal section, with the first two checked initially. Export
  only selected signals; disable the download for an empty selection. Preserve
  selections across closing/reopening and navigation. Scans retain all signals.
- Advance Time evolution settings to format 2.5 for the moment defaults and
  independent saved probe selections. Only the new format is supported, as
  requested; no migration is performed. Physics and ACC settings are unchanged.
- Fix the Time evolution UI test harness to submit stateful expander values
  on reruns, matching browser behavior. This restores moment-selection,
  dual-species, probe-scan, navigation-persistence, and display-axis checks
  without changing application behavior or weakening their assertions.
- Document the test-only expander adapter and current-format validation scope.
- Validation: all 35 focused Time evolution tests pass, including 10 real-page
  UI checks. The complete suite passes all 227 tests with no failures or skips
  (153.399 seconds). Git diff whitespace checks pass.

## 6.20.7 - 2026-10-25

- Verify that Time evolution Save includes pending widget edits in the current
  format 2.4 JSON and correct filename. Update the preset test's stale 2.0
  expectation. Remove historical-format fixtures and tests, including skipped
  migration tests, while retaining current-schema and field validation.
- Fix stalled probe-scan stationary convergence by tightening BDF tolerances
  when state changes are small but the residual stops improving. Retry from
  the current state over the initial interval instead of doubling settled
  trajectories; report an error after three unsuccessful accuracy refinements.
- Preserve the stationary residual threshold, two consecutive small-change
  checks, and density-matrix trace, Hermiticity and positivity checks. Ordinary
  converging cases keep their original integration tolerances; transient
  evolution, RF calculations and saved-condition schemas are unchanged.
- Add regression coverage for the original transverse-field stall, default
  accuracy on an easy case, bounded failure, stationary initial states, and
  pumped dual-species stationarity and physical density matrices.
- Representative timings on the development machine: the original Rb87
  transverse-field case completes in 1.26 s (previously stopped after 100 s);
  aligned and self-exchange cases remain about 0.21 s and 0.38 s, compared with
  0.21 s and 0.36 s previously. A 5 nT transverse field, tenfold pump intensity,
  and Rb87/Rb85 cross exchange converge in about 19 s, 14 s, and 14 s. Linear
  cases agree with direct stationary solutions within 7e-13; the nonlinear
  dual-species result agrees with a stricter integration within 3.1e-13.
- All seven probe-scan tests and the release-metadata check pass (8 tests).
  All 36 focused current-format settings and Save checks pass. Five existing
  Time evolution UI tests still have plot-count or missing-widget failures,
  reproduced on an isolated copy of committed HEAD. These checks do not
  load or save settings; the stale preset-schema expectation is corrected.

## 6.20.6 - 2026-10-25

- Format small-range Time evolution plots with three-significant-figure
  scientific notation on y-axis ticks and hover values; show zero ticks as 0.
  Keep the existing whole-graph threshold of absolute extent below 0.01 and
  full-precision calculations/CSV exports.
- Remove redundant legends from the single-curve hyperfine-polarization and
  electron-polarization graphs. Multi-curve moment and Zeeman legends remain.
- Update the focused UI check and documentation. Settings and physics are unchanged.
- Validation: the single page-render UI check passed in 7.348 seconds, including
  scientific tick/hover formatting and hidden single-curve polarization legends.

## 6.20.5 - 2026-10-25

- Fix the Time evolution layout exception introduced by passing `wrap=False`
  to vertical section-control containers. Disable wrapping only for the
  horizontal scan toolbar; preserve the requested help/download arrangement.
- Shorten the scan slider's visible and accessible label to Detuning range (MHz).
- Confirmed the incompatible parameter combination from the installed Streamlit
  implementation. Tests and browser checks were not run, following the user's
  instructions.

## 6.20.4 - 2026-10-25

- Halve the probe-scan range slider's width and place its label on the left.
  Align help immediately left of download on the same header row, using native
  Streamlit layout controls. Other section controls keep their current layout.
- Settings and calculations are unchanged. Tests and browser checks were not
  run, following the user's instructions.

## 6.20.3 - 2026-10-25

- Set the Time evolution probe detuning slider limits to ±15000 MHz, with
  a default selection of ±7000 MHz for each species. Retain the 10 MHz step
  and minimum span. Synchronize help and saved-setting validation.
- Advance Time evolution settings to schema 2.4 for the changed defaults and
  bounds; only the new format is supported under the current compatibility policy.
  Update documentation and existing checks. Tests and browser checks were not
  run, following the user's instructions.

## 6.20.2 - 2026-10-25

- Replace the two probe-scan detuning sliders with one native two-handle range
  slider at the top left. Preserve independent per-species bounds, the ±30000 MHz
  limits, 10 MHz minimum span, collapsed lazy rendering and existing saved format.
- Update help, documentation and existing UI expectations. Tests and browser
  checks were not run, following the user's instructions.

## 6.20.1 - 2026-10-25

- Start Atomic moments and Probe signal scanning collapsed on Time evolution.
  Use native Streamlit open-state tracking to create plots, controls and section
  CSVs only while open. Defer pump-on steady-state convergence and probe spectra
  until the scan is opened; retain cached calculations and saved display choices.
- Settings formats and physics are unchanged. Update existing UI expectations
  for collapsed defaults and reopening. Tests and browser checks were not run,
  following the user's instructions.

## 6.20.0 - 2026-10-25

- Add Probe signal scanning to each Time evolution species tab, with one
  frequency plot per row for rotation, ellipticity and absorption. Two detuning
  sliders at the top left scan the sidebar probe's D1/D2 line relative to its
  zero-pressure fine-structure center. All allowed hyperfine channels contribute,
  with transition-center markers. Help and an icon CSV download sit at top right.
- Independently converge the pump-on coherent master equation from the
  unpolarized state using adaptive BDF intervals. Check stationarity, state
  convergence, trace, Hermiticity and positivity; report nonconvergence rather
  than plotting an unconverged state. The scan is independent of transient timing
  and display axes, including self/cross spin exchange for both species.
- Share complex polarizability and Jones/Stokes propagation between transient
  and spectral readout. Contract unit-coefficient operators once per state and
  combine spectral coefficients for 2001 scan points. CSVs contain plotted units,
  frequency origin, pump/cell/probe settings and convergence information.
- Use a native 420-pixel initial sidebar width on Time evolution. Settings
  schema 2.3 adds independent per-species scan bounds, default ±12000 MHz,
  adjustable within ±30000 MHz with a minimum span of 10 MHz. Older settings
  are rejected without migration, per the user's explicit preference.
- Added physical-consistency and export coverage and updated UI expectations.
  Tests and browser checks were not run, following the user's instructions.

## 6.19.4 - 2026-10-25

- Remove the complete CSV download at the bottom of each Time evolution
  species tab and its unused page-side data preparation. The three section
  CSV downloads remain. Update page help and export documentation.
- Settings and calculations are unchanged. Tests and browser checks were not
  run, following the user's instructions.

## 6.19.3 - 2026-10-25

- Resolve atomic alignment names to the selected moment axis throughout
  checkbox labels, graph headings, y-axis labels, help and section CSV column
  names. The cyclic frames give Qzz, Qxy, Qxx−Qyy for z; Qxx, Qyz, Qyy−Qzz
  for x; and Qyy, Qzx, Qzz−Qxx for y. Stable selection keys and calculated
  values are unchanged; each species' plots follow its own moment axis.
- Tests and browser checks were not run, following the user's instructions.

## 6.19.2 - 2026-10-25

- Give Time evolution section download buttons and both Axis selectors the
  same 56-pixel width. Center the Zeeman Axis caption above its selector,
  aligned with the download button in the third column.
- Settings and calculations are unchanged. Tests and browser checks were not
  run, following the user's instructions.

## 6.19.1 - 2026-10-25

- Shorten the Zeeman quantization selector label to Axis and fit it to the
  third column below the help/download controls.
- Move the atomic-moment Axis label and compact selector before the moment
  checkboxes, with the label on the selector's left and a 20-pixel spacer
  before the checkbox group. Axis values, persistence and calculations are unchanged.
- Tests and browser checks were not run, following the user's instructions.

## 6.19.0 - 2026-10-25

- Remove the Time evolution sidebar quantization axes. Each species now has
  a Zeeman-population axis at that section's bottom right and an independent
  atomic-moment axis to the right of the moment checkboxes.
- Evolve and cache the physical state in laboratory z; changing a display
  axis re-expresses the cached density matrix without repeating integration.
  The Zeeman axis defines population/coherence readout. The moment axis defines
  electron polarization and local alignment; laboratory Fx,Fy,Fz and scalar
  hyperfine polarization are invariant. Probe signals use the physical state.
- Section and complete CSVs identify both axes. Complete CSV populations and
  coherences use the Zeeman basis; longitudinal/transverse moments use the
  moment frame. Settings schema 2.2 replaces q_axis_A/B with zeeman_axis_A/B
  and moment_axis_A/B, initially z. Older formats are rejected without migration,
  following the user's explicit compatibility preference.
- Tests and browser checks were not run, following the user's instructions.

## 6.18.2 - 2026-10-25

- Remove the Plotted moments caption. Label the hyperfine/electron checkboxes
  with their normalized mathematical observables, ⟨S·I⟩/ℏ² and 2⟨Sq⟩/ℏ.
- Give each collapsible section two plot columns and a narrow third column
  containing help above the download icon. In Probe signal, absorption sits
  below optical rotation; ellipticity remains in the second column.
- Typeset mathematical y-axis labels with native Streamlit LaTeX beside the
  interactive charts, rotated vertically with page-scoped CSS. Subscripts,
  superscripts and expectation brackets render as math; numeric axes and CSV
  column names are unchanged. No dynamics or saved-setting changes.
- Tests and browser checks were not run, following the user's instructions.

## 6.18.1 - 2026-10-25

- Replace per-graph downloads and explanation popups on Time evolution with
  one download icon and one help popup at the top right inside each collapsible
  Zeeman populations, Atomic moments and Probe signal section, per species.
  Individual graph panels retain their titles and borders.
- Each section CSV combines all displayed curves into one time-indexed table.
  Column labels identify plot, curve and plotted units; shared timing, species
  and axis metadata appears once. Help combines that section's explanations.
  Empty moment selection disables its download; undefined optical samples are
  retained. The complete per-species CSV remains available below the sections.
- Render moment choices as native checkboxes so their mathematical labels use
  LaTeX. Plot-panel headings and section-help equations also use LaTeX. The
  Plotted moments caption stays on the left; stored selections are unchanged.
  Native chart-axis labels and CSV column names remain plain text.
- Settings and calculations are unchanged. Tests and browser checks were not
  run, following the user's instruction to stop testing.

## 6.18.0 - 2026-10-25

- Replace Longitudinal moments with a Plotted moments selector, with its label
  to the left. Hyperfine polarization and electron polarization are selectable,
  alongside laboratory Fx, Fy, Fz and local Qqq, Quv, Quu−Qvv. Selected plots
  use two columns; each manifold-resolved plot always includes both F levels.
- Evaluate Cartesian and transverse moments from the complete density matrix.
  The right-handed (u,v,q) frames are (y,z,x), (z,x,y), (x,y,z) for q=x,y,z.
  Graph help and CSV metadata identify these axes. The complete CSV includes
  all new moments regardless of plot selection; normalized graph units are 1.
- Time evolution settings use schema 2.1 with a `plotted_moments` selection;
  obsolete moment-type/manifold fields are removed. Older formats are rejected
  without migration, consistent with the requested current-format-only policy.
  Defaults show hyperfine/electron polarization, Fz and Qqq. Dynamics is unchanged.
- Tests and browser refresh were not run, following the user's instructions.

## 6.17.6 - 2026-10-25

- Remove the Hyperfine manifolds selector from Time evolution Atomic moments.
  Orientation and alignment always show both manifolds for each species,
  including when a saved display choice selected fewer manifolds. Graph CSVs
  include both curves. The longitudinal moment type selector remains.
- Saved file structure and physics are unchanged. No tests were run, following
  the user's instruction to stop testing.

## 6.17.5 - 2026-10-25

- Restore native Time evolution section outlines and header dividers. Remove
  tracked section keys and refresh-on-toggle behavior; expanding or collapsing
  uses normal client-side Streamlit behavior again.
- Remove the regression test specific to tracked section refresh. No tests
  were run for this reversal, following the user's instruction to stop testing.

## 6.17.4 - 2026-10-25

- Track Time evolution section open/closed states and refresh the results
  layout on every toggle, so the page scroll range updates with visible content.
  Mixed section states are preserved through that refresh. Atomic dynamics
  stays cached; toggling does not resolve the atoms or change inputs.
- Invisible section outlines and individual plot borders remain. Settings
  formats, physics and exports are unchanged.

## 6.17.3 - 2026-10-25

- Make the outlines and header dividers of Time evolution's three collapsible
  result sections transparent, retaining headings, toggle controls, layout spacing
  and individual graph-panel borders. Styling is scoped to this page's
  section keys. Settings, calculations and exports are unchanged.
- Verification: live browser confirms transparent outlines/dividers and working
  collapse/expand controls; full suite: 226 tests run, 11 existing skips, no failures.

## 6.17.2 - 2026-10-25

- Set the Time evolution electron plot title to Electron polarization and its
  y-axis label to 2⟨S_q⟩/ℏ.
- Reduce bottom whitespace on all Time evolution plots: 5-pixel bottom chart
  padding and 0.25rem bottom panel padding, scoped to these plot panels.
  Chart height remains 260 pixels. Settings and calculations are unchanged.
- Verification: live browser confirms all nine panels use 4-pixel bottom
  padding, 260-pixel charts and clear x-axis labels; electron title/axis checked.
  Full suite: 226 tests run, 11 existing skips, no failures.

## 6.17.1 - 2026-10-25

- Align each Time evolution plot title, Download CSV button and question-mark
  popup in one native horizontal header row. Keep the title on the left and
  controls on the right, vertically centered. Remove the duplicate chart title
  below that row and retain it as the chart's accessibility description.
- Chart regions remain 260 pixels tall; calculations, exports and settings
  formats are unchanged.
- Verification: all nine native title/control header rows checked in the real
  Streamlit page; full suite: 226 tests run, 11 existing skips, no failures.

## 6.17.0 - 2026-10-25

- Remove the Time evolution hyperfine-population panel. Add the scalar
  hyperfine polarization ⟨S·I⟩/ℏ² to Atomic moments, alongside electron
  polarization on the first plot row. Selected orientation/alignment plots
  occupy a separate second row.
- Calculate the correlation as [F(F+1)−I(I+1)−3/4]/2 averaged over the full
  ground state, without manifold renormalization. Include it in graph and
  complete CSV exports; retain hyperfine totals in the complete export.
- Rename probe Output ellipticity to Ellipticity and the electron plot title
  to 2⟨S_q⟩/ℏ. Hide legends on the three single-curve probe plots.
- Atomic dynamics and settings formats are unchanged. Invalidate cached
  evolution results for the newly added observable.
- Verification: direct electron–nuclear spin-operator expectations, unpolarized
  and pure-manifold limits for all isotopes, quantization-axis independence,
  CSV values, real Streamlit row placement, optical legends and deselection
  behavior. Full suite: 226 tests run, 11 existing skips, no failures.

## 6.16.3 - 2026-10-25

- Remove the Fa/Fb prefixes from Time evolution Zeeman plot titles. Hyperfine
  population, orientation and alignment legends now use the isotope's actual
  F values for each species. Per-graph CSV curve names match the legends.
  Calculations, complete CSV field names and settings formats are unchanged.
- Verification: actual Rb87/Cs133 Streamlit titles and legends checked;
  full suite: 225 tests run, 11 existing skips, no failures.

## 6.16.2 - 2026-10-25

- Use three significant figures for Time evolution quantity labels and hover
  values only when the graph's entire vertical data extent, across all shown
  curves and including zero, lies strictly between -0.01 and +0.01 in plotted
  units. Other graphs retain automatic axis labels and their original
  six-significant-figure hover format, even for individual small samples.
  Coordinates, calculations and CSVs retain full precision; settings formats
  are unchanged.
- Verification: whole-graph small, mixed and threshold ranges checked;
  real Streamlit page renders all nine graphs. Full suite: 225 tests run,
  11 existing skips, no failures.

## 6.16.1 - 2026-10-25

- Format Time evolution plotted quantities to four significant figures in
  y-axis labels and hover tooltips. Calculations, curve coordinates and CSV
  exports retain full precision; settings formats are unchanged.
- Verification: all nine chart specifications checked through the real
  Streamlit page; 225 regression tests run, 11 existing skips, no failures.

## 6.16.0 - 2026-10-25

- Replaced Time evolution's diagonal population dynamics with a ground-HFS
  secular density-matrix master equation. All within-F Zeeman coherences,
  signed-g-factor magnetic precession, collective optical repopulation and
  pump AC-Stark dynamics are retained. Both optical terms stop at pump turn-off.
- Each alkali's quantization axis is independent of the field and partner's
  axis. Self/cross exchange uses full nuclear/electron reductions in a shared
  laboratory frame. ER and uniform relaxation act on the coherent state.
- Weak-probe Jones propagation now uses the full density matrix. Complete
  CSVs include independent real/imaginary Zeeman coherences. Display sampling
  follows Larmor harmonics, with an explicit warning at the 20,001-sample cap.
- Time evolution settings require format 2.0 and the new model identifier;
  older files are rejected without migration, as requested. ACC settings and
  calculations retain their existing behavior. Magnetic shifts in optical
  detunings and inter-ground-hyperfine coherences remain excluded.
- Verification: signed Larmor rotation, exact dark unitary evolution,
  trace/Hermiticity/positivity, collective optical covariance, collinear
  population limits, angular momentum conservation, dual-species basis and
  probe invariance, new-format validation, and real Streamlit transverse-axis
  controls are covered by regression tests.
- Full suite: 225 tests run, 11 existing skips, no failures. Live-browser check
  confirms damped probe precession with a 10 nT x field and a z quantization
  axis, preserving the existing graph layout and CSV controls.

## 6.15.2 - 2026-10-25

- Reduced all Time evolution plot regions from 520 to 260 pixels, retaining
  top-right legends, download buttons, and explanation popups.
- Verification: 214 tests run, 11 existing skips, no failures; Streamlit
  checks confirm every plotted region is 260 pixels tall.

## 6.15.1 - 2026-10-25

### Changed

- Added top-right CSV and explanation-popup controls to every Time evolution
  graph panel. Each graph exports only its displayed curves, in plotted units,
  with one column per curve and timing/species metadata. Complete species CSVs
  remain available, including all signals and configuration.
- Moved page-level descriptions into a question-mark popup beside the page
  title, and moved graph explanations into their panel popups.
- Doubled plot regions from 260 to 520 pixels and placed legends inside the
  top-right corner. Zeeman legend ordering follows the numerical m sequence.
  Event timing, calculations, and saved settings formats are unchanged.

### Verification

- Full unittest suite: 214 tests run, 11 existing skips, no failures.
- Checked per-graph CSV units, displayed-curve selection and event metadata;
  graph/page popup counts for single and dual species; 520-pixel chart regions;
  top-right legends; and undefined-angle downloads. Visually verified the live
  page header, help popup, and narrow three-column graph panels in the browser.

## 6.15.0 - 2026-10-25

### Changed

- Removed the remaining Time evolution settings compatibility behavior:
  missing-field defaults, ignored unknown fields, and extended duration ranges.
  Schema 1.5 requires a complete current configuration and enforces total time
  from 0.001 to 30 s, with turn-off from zero to total time. No older schemas
  are accepted. ACC conditions remain at schema 6.14.
- Shortened the slider labels to Total time (s) and Pump turn-off time (s).
- Replaced the Time evolution transmission plot with fractional absorption,
  1 − I_out/I_in, and exported the same absorption values alongside transmission.
  Uniform-path Jones propagation and the atomic evolution are unchanged.

### Verification

- Full unittest suite: 212 tests run, 11 existing skips, no failures.
- Checked strict complete-file loading, fixed 30 s limits, shortened slider
  labels, absorption plot titles, zero-path and optically thin absorption,
  and matching absorption/transmission CSV values.

## 6.14.0 - 2026-10-25

### Changed

- Removed Only show allowed hyperfine transitions from ACC Magnetometer and
  Time evolution. Pump/probe reference menus and ACC transition tables now
  always filter allowed transitions. Optical calculations continue to include
  all relevant allowed channels, including off-resonant excitation.
- Renamed the Time evolution sidebar's Display section to Event Sequence and
  moved the two time sliders into one row below that heading.
- Removed the saved show_allowed_only field. ACC now writes/loads schema 6.14
  and Time evolution schema 1.4. Per the user's compatibility choice, older
  analysis formats are rejected and their migration paths are removed.

### Verification

- Full unittest suite: 211 tests run, 11 existing skips, no failures.
- Streamlit checks cover allowed-only pump/probe menus on both pages,
  stale checkbox state, Event Sequence heading and slider row, timing
  persistence, and current-format round trips with older-file rejection.

## 6.13.1 - 2026-10-25

### Changed

- Restricted the four D1/D2 probe presets to Time evolution. ACC Magnetometer
  again lists only Custom Weak and its original weak/physical pump sources.
- Selecting a preset fills the probe fields once; every field remains
  editable while the preset stays selected. Edits survive reruns, navigation,
  and save/load. Pump-linked weak sources keep their existing behavior.
- Recorded the requirement to ask about older saved-condition compatibility
  before future condition-sidebar changes. This UI correction retains the
  existing saved field formats and adds no migrations.

### Verification

- Full unittest suite: 210 tests run, 11 existing skips, no failures.
- Real Streamlit checks confirm both ACC menus exclude presets and all Time
  evolution preset fields remain editable, retaining edits through reruns
  and navigation. Edited settings also pass save/load round trips.

## 6.13.0 - 2026-10-25

### Added

- Added Preset D1a, D1b, D2a, and D2b weak-probe sources independently for A
  and B on ACC Magnetometer and Time evolution. All use z, a 1 cm path,
  azimuth 0°, and ellipticity angle 0°. The a branches use
  I+1/2→I−1/2 at −2000 MHz, and b branches use I−1/2→I+1/2 at +2000 MHz,
  relative to the selected pressure-shifted D1/D2 transition.
- Preset fields follow the selected isotope and stay fixed until Custom
  Weak is selected for editing. Existing startup defaults and pump physics
  are unchanged.

### Compatibility

- ACC condition schema 6.13 and Time evolution settings schema 1.3 save the
  new source choices. ACC 6.12 and Time evolution 1.2 load unchanged; all
  previously supported migrations retain their stored-value semantics.

### Verification

- Full unittest suite: 210 tests run, 11 existing skips, no failures.
- Real Streamlit checks cover all four presets on both pages, isotope
  changes, Custom Weak editing, A/B independence, and navigation persistence.

## 6.12.1 - 2026-10-25

### Changed

- Updated new Time evolution settings for both species: Pump1 uses D1
  Fa−1→Fa at zero detuning with sigma+ light along z and retains 10 µW/cm²;
  Pump2 uses Fa→Fa with zero intensity and sigma+ along z. Pump3 is unchanged.
- New weak probes are Custom Weak along z, with a 1 cm path, the D1 upper-to-
  lower hyperfine branch (Rb87 2→1), −2000 MHz relative detuning, and zero
  ellipticity. Their input azimuth remains 45°. Total time defaults to 4 s
  and pump turn-off to 1 s. The atomic initial state remains unpolarized.
- Advanced Time evolution settings to v1.2. Versions 1.0 and 1.1 preserve all
  stored inputs and use their original defaults for omitted fields. ACC
  startup settings and conditions are unchanged. Activating B translates
  the displayed surrogate spectrum to the same relative hyperfine branches.

### Verification

- Full unittest suite: 204 tests run, 11 existing skips, no failures.
- Verified all requested startup values in the real Streamlit page, relative
  branches after Rb85/Cs133 selection, independent ACC defaults, and migration
  of older stored configurations without applying the new startup values.

## 6.12.0 - 2026-10-25

### Added

- Added independent R_u,A and R_u,B controls below R_ER on both analysis
  sidebars, defaulting to 1 s⁻¹ in new settings. Uniform relaxation uses
  `D_u[rho] = R_u (Tr(rho) I/N - rho)`, resetting the entire ground-state
  distribution and damping all coherences at R_u, including nuclear moments.
- Included the channel in single/dual-species steady states and population
  Jacobians, pump-on/dark dynamics, secular RF generators, weak optical and
  physical-pump feedback responses, and independent-transition RF helpers.
- Added uniform population/coherence rate columns to Zeeman diagnostics and
  CSVs, included them in total relaxation and linewidth normalization, and
  exposed the uniform population generator in the optional matrix display.
- Advanced ACC conditions to 6.12 and Time evolution settings to 1.1. Complete
  ACC v6.9/v6.10 and dynamics v1.0 files migrate with R_u=0 to preserve their
  previous relaxation model. Current files require valid, finite, nonnegative
  rates; A/B values and page-specific configurations remain independent.

### Verification

- Full unittest suite: 203 tests run, 11 existing skips, no failures.
- Verified trace-preserving uniform decay, independent A/B rates, agreement
  between long-time evolution and the steady state, RF linewidths, physical
  pump readout, zero-rate recovery, rate diagnostics, settings migrations, and
  real Streamlit controls and navigation persistence.

## 6.11.3 - 2026-10-25

### Changed

- Moved both Time evolution sliders into one row in the sidebar immediately
  below the probe configuration. The normal total-time range now ends at 30 s, and
  pump turn-off remains bounded by the selected total time. Existing saved
  durations above 30 s retain their values until explicitly reduced.

## 6.11.2 - 2026-10-25

### Changed

- Matched the Time evolution probe configuration layout to ACC Magnetometer:
  source/direction/path length in a 45%/20%/35% row, reference line and F→F′
  side by side, detuning beside its line-center caption, and azimuth beside
  ellipticity. Applies to both probes with their independent settings retained.

## 6.11.1 - 2026-10-25

### Changed

- Shortened the ACC Magnetometer and Time evolution sidebar labels from
  Hyperfine transition to F→F′ for pumps and probes, and from Intensity to I
  for pumps, retaining the µW/cm² units.

## 6.11.0 - 2026-10-25

### Added

- Added **Analysis → Time evolution**, with independent atom/cell, six-pump,
  quantization-axis, and weak-probe configurations. All ground Zeeman states
  start equally populated. Two sliders set pump turn-off time t₁ and total
  plotted time t₂; every curve spans zero through t₂ with a turn-off marker.
- Integrated the nonlinear single/dual-species population equations in separate
  pump-on and dark intervals. Electron randomization and self/cross spin
  exchange continue after turn-off, and exchange uses instantaneous marginals.
  The new path does not call the RF solver or physical-pump feedback closure.
- Added Fa/Fb population totals, separate individual-Zeeman population plots
  for each manifold (without manifold totals), electron polarization, and
  selectable longitudinal orientation and traceless alignment per total atom.
- Added separate weak-probe rotation, output ellipticity, and transmission
  plots. Uniform-state Jones propagation reuses the existing Mathur tensor,
  includes scalar absorption, and retains accumulated phase and attenuation.
  The non-perturbing probe stays on when pumps are off. Undefined rotation at
  circular input polarization is reported rather than displayed as zero.
- Added independent Time evolution settings format v1.0, persistent navigation
  state, and CSV exports with population, moment, optical, timing, and model
  metadata. ACC conditions remain at schema 6.10 and its solver is unchanged.

### Verification

- Full unittest suite: 190 tests run, 11 existing skips, no failures.
- Checked the nonlinear derivative against the existing collision maps, the
  linear limit against matrix exponentials, probability and angular-momentum
  conservation, timing endpoints, probe sign/units and zero-path limits, and
  settings serialization. Streamlit AppTest exercises both species, slider
  bounds, moment/probe controls, and independent settings across navigation.

## 6.10.1 - 2026-10-24

### Changed

- Added an Atomic properties sidebar and moved the alkali isotope selector
  into it, preserving its default and persistent selection across pages.
- Arranged the Transition strengths polarization checkboxes (σ−, π, σ+)
  in one horizontal row.
- Removed the description caption below the Atomic properties page title.
- Swapped the Transition strengths hyperfine selectors so Excited F′ appears
  above Ground F.
- Moved temperature to the Atomic properties sidebar and shared it between
  thermal and buffer-gas calculations. Added self electron-randomization rates
  from editable binary spin-destruction cross sections.
- Buffer-gas inputs now represent fill pressures at 25 °C in a sealed cell.
  Pressure, temperature, and cross-section edits immediately recalculate
  operating pressure, broadening, shift, and electron-randomization rates.
  Optical coefficients and entered cross sections are held constant; unknown
  cross sections are displayed as unavailable until a value is supplied.
- Combined buffer-gas inputs and results into one table: fill pressure and
  electron-randomization cross section are editable cells; gas names and
  calculated results remain read-only. Removed the separate input boxes.
- Shortened the buffer-gas cross-section and rate column headings to σ_ER
  and Γ_ER, retaining their units.
- Shortened the pressure headings to P_fill at 25 ℃ and P_T (Torr), and placed σ_ER
  immediately before Γ_ER.

## 6.10.0 - 2026-09-19

### Changed

- Renamed the probe-source option `Custom` to `Custom Weak` for both species,
  preserving its saved value and weak-probe behavior.
- ACC Magnetometer now applies one RF magnetic field to both ground manifolds,
  with signed couplings 1 and g_b/g_a relative to the upper-manifold angular
  frequency. Fa and Fb atomic moments can be displayed separately or together.
- Single- and dual-species RF, optical, and physical-pump source responses use
  a shared secular coherence model. It retains ER/SE transfer between the two
  manifolds while removing coherences across the ground hyperfine splitting.
  The single-species RF response therefore includes coherence transfer instead
  of the former independent-linewidth approximation. Dual-species pump feedback
  now reuses the coupled population and coherence dynamics.
- Probe signals include both manifold responses with their own optical
  coefficients. Fa/Fb complex contributions add before amplitude is evaluated.
  Physical-pump manifold contributions are projections of the same feedback
  closure; existing separately closed rank counterfactuals are retained.
- Added manifold selectors, probe manifold decomposition, colored transition
  markers for both manifolds, and manifold-resolved CSVs with model, phase,
  magnetic-coupling, and normalization metadata. Relaxation normalization uses
  one reference transition selected across both manifolds.
- Corrected the nuclear contribution to the hyperfine Landé factors: positive
  nuclear magnetic moments contribute negative magnetic energy in the selected
  Hamiltonian convention. Static and RF magnetic couplings share this correction.
- Corrected quadrupole CSV units to hbar² rather than hbar.
- Condition schema 6.10 adds persistent manifold display settings. Complete
  v6.9 files migrate without changing stored inputs, retain the Fa atomic view,
  and are recalculated with the new model. New sessions show both atomic
  manifolds. There is no legacy physics/comparison mode.
- Archived the complete pre-upgrade working version, including local
  documentation edits, in `../archive/alkali_pumping_v6.9.9/`, with a SHA-256
  manifest. Repository metadata and Python caches are omitted from the snapshot.

### Verification

- Full unittest suite: 165 tests run, 11 pre-existing skips, no failures.
- Added analytic-response, signed-g, coherence-transfer, direct-resolvent,
  zero-cross-exchange, shared physical-feedback, migration, and export checks.
- Exercised real Streamlit manifold controls, both species tabs, and navigation
  persistence through AppTest. Cached calculations include the physics version.

## 6.9.9 - 2026-09-06

### Changed

- Renamed the Magnetometry navigation item and page to **ACC Magnetometer**,
  short for antirelaxation coated cell magnetometer, and renamed its source
  file to `ACC-Magnetometer.py` with route `acc-magnetometer`.
- Made **Analysis → ACC Magnetometer** the default page.
- Moved **Light shift** and **Atomic polarizability** from the Analysis menu to
  the Reference menu alongside **Atomic properties**.
- Added an unnumbered, cross-linked Quarto book containing an initial user
  manual, physical theory, technical notes, glossary, bibliography, and links
  to the retained legacy derivations.

## 6.9.8 - 2026-09-06

### Changed

- Removed slice-by-slice atomic-state changes from physical-pump readout.
- Physical-pump readout now uses one spatially uniform RF-driven density-matrix
  perturbation for the full cell, closed self-consistently against the
  optical-path average of the generated fractional intensity and normalized
  Stokes changes.
- The entered pump intensity remains unchanged; the uniform-state model does
  not introduce beam/cell filling-factor intensity averaging.
- Full dispersive and dissipative CBOR/LDOR feedback remains active, but no
  independent downstream-slice atomic states are created.

## 6.9.7 - 2026-09-06

### Changed

- Replaced the physical-pump readout's light-shift-only Stokes feedback with a
  combined dispersive and dissipative ground-state Liouvillian source.
- Added fractional-intensity feedback and polarization-dependent optical
  pumping, including trace-preserving population redistribution and secular
  ground-manifold Raman coherences generated by Stokes changes.
- Evolved the new dissipative sources with the existing population and
  coherence Jacobians while preserving the weak-RF approximation and the
  counterfactual Orientation/Alignment displays.

## 6.9.6 - 2026-09-06

### Changed

- Added a Detuning reference selector beside the Atomic polarizability D-line
  selector. It offers the zero-pressure line center and every allowed
  hyperfine transition, and consistently shifts the automatic range, plot
  axis, transition markers, and CSV detunings.
- Added a Light shift Scattering rate checkbox and a Components-view
  multi-select dropdown for the scalar, vector, and tensor plots. Component
  selections also determine which component curves are exported.
- Added the F-independent common scalar shift to the Light shift scalar panel
  and CSV. It is calculated as the ground-state-degeneracy-weighted mean of
  the two hyperfine-manifold scalar shifts.
- Placed the Light shift azimuth and ellipticity sliders on one row and swapped
  the sidebar positions of Shift units and Sweep points.
- Made the Light shift component plot headings black and approximately 10%
  larger, and relabeled the Transition-view field axis as Equivalent
  field/intensity with its units.
- Advanced Atomic polarizability conditions to v1.1 and Light shift conditions
  to v1.2, with migrations for the prior formats.

## 6.9.5 - 2026-09-04

### Changed

- Renamed the Atomic polarizability y-axis labels for
  $\alpha_{\mathrm{gt}}$ and $\alpha_{\mathrm{br}}$ to Vector polarizability
  and Tensor polarizability. The help popover retains Mathur's coefficient
  names and explains the rank-1/rank-2 interpretation and normalization.

## 6.9.4 - 2026-09-03

### Changed

- Added independent column-visibility selectors to the Alkali A and B Zeeman
  tables. The `F` and `m` columns remain fixed; all optional columns except
  `P_F` and `Λ` are visible by default; and CSV downloads continue to contain
  the complete table. The
  selector is positioned at the table area's upper-left, and its menus expand
  vertically to show all choices when the viewport permits.
- Reduced the horizontal margins of the Magnetometry main panel without
  changing the layout of other pages.
- Increased the condition/control columns in both Magnetometry response
  expanders from 18% to 23% of the available width.
- Changed every `Download CSV` button to size itself to its caption instead of
  stretching across its layout column.
- Right-aligned each caption-sized CSV button so it sits immediately beside
  its corresponding help popover.
- Set the probe-configuration row to 45% source, 20% direction, and 35% path
  length, and shortened the direction label.

## 6.9.3 - 2026-09-03

### Changed

- Moved the shared RF settings outside the atomic-response expander.
- Simplified the population-plot heading to `Populations`, removed the
  redundant Zeeman-table heading, and clarified in the enclosing expander
  label that the displayed properties belong to the ground state.

## 6.9.2 - 2026-09-03

### Changed

- Made the atomic RF response and optical probe response independently
  collapsible in default-open expanders on each magnetometry result tab.

## 6.9.1 - 2026-09-03

### Changed

- Made the complete Zeeman sublevel-properties area, including the ground-state
  population plot and population/rate summary, collapsible in a default-open
  expander.
- Changed both species' Pump1 defaults to the D1 upper-to-upper hyperfine
  transition at +450 MHz and 10 µW/cm², and Pump2 defaults to the D1
  lower-to-upper transition at zero detuning and 10 µW/cm².
- Changed the default weak probe links to PumpA1 and PumpB1.

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
