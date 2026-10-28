# Alkali Pumping v6.21.1

This Streamlit application models optical pumping, electron randomization,
uniform ground-state relaxation, self spin exchange, and unlike-alkali spin
exchange in a one- or two-alkali vapor.

The main-content gap below the navigation is half Streamlit's default on all
pages. The pinned bottom toolbar's lower padding is also halved (56 to 28 px
with the default 16 px base font). The sidebar's upper space is reduced from
76 to 38 px, including its compact collapse-button header.

The app starts on **Analysis → ACC Magnetometer**, and the top navigation separates
five workflows:

- **Analysis → ACC Magnetometer** models an antirelaxation coated cell
  magnetometer with dual-alkali population and Zeeman diagnostics, independent
  RF-A/RF-B drives, atomic-moment susceptibilities, and independent weak or
  physical Pump-A/Pump-B Stokes readouts. Its Zeeman tables have per-species
  column selectors while retaining complete CSV exports.
- **Analysis → Time evolution** independently configures the same atoms, cell,
  pumps, and weak probes. Initially unpolarized ground-state populations evolve
  with all pumps on from zero to the turn-off slider `t₁`, then relax with the
  pumps off through the total-time slider `t₂`. Self and cross spin exchange
  follow the instantaneous density matrices, including all within-hyperfine
  Zeeman coherences. Magnetic precession and pump AC-Stark evolution are included.
  Transient integration automatically chooses BDF for damping-dominated or short
  intervals and DOP853 for sustained coherent motion without faster damping.
  Each pump-on/dark interval is selected separately, with unchanged error
  tolerances and physical-state checks. A failed or nonfinite DOP853 solve
  restarts that interval with BDF. The steady-state probe scan retains BDF.
  Both solvers reuse prepared constant matrices and collision maps; BDF uses
  an analytic Jacobian, including self and cross spin exchange.
  Each alkali has independent Zeeman-population and atomic-moment display axes,
  which can differ from the magnetic-field axis. The compact Zeeman **Axis**
  selector sits in its section's top row, with its label on the left. Axis
  selectors share the download buttons' width. The moment **Axis** label sits to the left
  of its selector; both precede the checkboxes with a gap before the group.
  Both replace the former sidebar axis. Changing either re-expresses the cached
  density matrix without rerunning integration or changing probe signals.
  Plots show individual Zeeman populations (one plot per manifold, without
  totals), hyperfine polarization ⟨S·I⟩/ℏ², electron polarization,
  selectable longitudinal orientation/alignment, and separate
  weak-probe rotation, output ellipticity, and absorption signals. The probe
  stays on throughout; RF and physical-pump feedback are absent. Settings and
  CSV exports are independent of ACC Magnetometer.
  The sidebar has separate **Open** and **Save** buttons beside the setting name.
  Save downloads committed settings immediately in the browser without rerunning
  the app. It uses the currently typed setting name and leaves other unfinished
  edits unapplied. Open validates the current settings format before loading it.
  Both timing sliders share one row under **Event Sequence** below the probe configuration, with
  a maximum total time of 30 s and turn-off bounded by total time. Their labels
  are Total time (s) and Pump turn-off time (s). Absorption is the fraction
  $1-I_{\mathrm{out}}/I_{\mathrm{in}}$ and is included in CSV exports.
  Each section places help and download at the right end of its top control row.
  Time-domain plots use two equal columns across the full width. The CSV contains all displayed curves in that
  section, with plot/curve names and units in the column labels. The four section
  downloads export plotted curves; the population title row separately exports the full steady density matrix. Individual graph panels retain their
  titles. Plot regions are 260 pixels tall
  with compact bottom margins; the electron plot is titled Electron polarization
  with y-axis label 2⟨S_q⟩/ℏ. Native moment checkboxes have no group caption;
  the hyperfine/electron choices show ⟨S·I⟩/ℏ² and 2⟨S_q⟩/ℏ.
  Mathematical checkbox labels, plot headings, y-axis labels and help equations
  use LaTeX; CSV names remain plain text. A row of probe checkboxes selects
  Absorption, Optical rotation and Ellipticity independently for each species;
  Absorption and Optical rotation start checked. Only selected curves are
  plotted and exported, and an empty selection disables that section's download.
  Moment and probe plots initially follow their default checkbox order, reading
  left to right then down. Newly checked plots append at the end; unchecking
  removes a plot without reordering the others, and rechecking appends it again.
  The order survives section closing, navigation and saving/loading settings;
  section CSV columns follow the same order.
  Choose hyperfine/electron polarization, laboratory
  ⟨F_x⟩, ⟨F_y⟩, ⟨F_z⟩, or local ⟨Q_qq⟩, ⟨Q_uv⟩, ⟨Q_uu−Q_vv⟩.
  Alignment labels resolve to actual Cartesian components: for moment axis z,
  they show ⟨Q_zz⟩, ⟨Q_xy⟩ and ⟨Q_xx−Q_yy⟩, including in section CSV headers.
  Selected plots use two columns and always show both hyperfine manifolds
  where applicable. The right-handed (u,v,q) frames are (y,z,x), (z,x,y),
  (x,y,z) for q=x,y,z. Defaults show hyperfine polarization and Fz.
  Collapsible result sections use their native outlines. Zeeman populations,
  Atomic moments and Steady-state probe signal scanning start collapsed; Probe signal starts
  open. All sections create and refresh plots and controls only while open,
  including after sidebar edits, using native Streamlit open-state tracking.
  Closing all time-domain sections defers transient integration until one is
  opened; opening only the scan computes its steady state independently.
  **Auto update** and **Apply** share a native toolbar pinned to the bottom of
  the main area, staying visible while plots scroll. Auto update starts checked.
  Uncheck it to batch physical-condition edits,
  then press **Apply** once. Pending curves and downloads are hidden until
  applied. Display and probe-readout edits continue to reuse the physical state.
  Transient caches retain two small trajectories and one trajectory of at least
  8 MiB. Four display-cache entries retain populations and moments without
  duplicated full density trajectories. Cached results are reused only for
  matching inputs; eviction can require recalculation.
  Selections, axes and scan limits survive closing and reopening. Individual
  graph panels retain borders.
  Legends are inside the top right. Single-curve probe, hyperfine-polarization
  and electron-polarization plots omit legends. Manifold legends show F values;
  Time-domain Zeeman titles show only F. Steady-state population bars place F
  on their vertical axes and omit gridlines.
  Graphs whose entire vertical data extent lies between -0.01 and +0.01 use
  scientific notation with three significant figures for y-axis ticks and
  hover values, with zero ticks shown as 0; other graphs retain their usual formatting.
  Calculations and CSV exports retain full precision.
  Page-level descriptions are in question-mark popups. Time evolution omits
  its page title and keeps its help beside Auto update and Apply in the bottom toolbar.
  New sessions use sigma+ D1 pumping along z: Pump1 on the lower-to-upper
  hyperfine branch at zero detuning, and Pump2 on the upper-to-upper branch
  with zero intensity. The Custom Weak probe is along z with a 1 cm path,
  the Rb87 2→1 reference transition, −2000 MHz detuning, and zero input
  ellipticity. Total time defaults to 4 s and turn-off to 1 s.
  **Steady-state probe signal scanning**, when opened, independently converges the pump-on steady state
  and scans all allowed hyperfine transitions of the probe's selected D line.
  If numerical error stalls stationary convergence, the solver tightens its
  integration accuracy with bounded retries rather than extending settled
  trajectories indefinitely. Normally converging cases retain their tolerances.
  The left 30% shows steady-state Zeeman populations using the same bars as
  ACC; the right 70% shows the selected probe signals, one plot per row.
  The two population manifolds are 338 pixels tall in total, 30% taller than
  one signal graph, with gaps about 30% of each bar's width. Three signal checkboxes in one row
  to the left of detuning selects optical rotation, ellipticity and absorption;
  optical rotation and absorption start checked. Each species' selection and
  plot order persist through navigation and saved settings. Scan CSVs include
  the selected signals; clearing all signals retains population bars and exports.
  Each section remembers its open or closed state across page navigation,
  independently for Alkali A and B.
  Scan CSVs include these populations and their quantization axis. The Zeeman
  populations title row also downloads the entire steady-state density matrix
  as a CSV with all real/imaginary elements and lab-z basis labels.
  A two-handle detuning range slider above the signals defaults to ±7000 MHz relative to the
  zero-pressure D-line center, adjustable within ±10000 MHz. Its label appears
  on one line to the left of the slider, whose length is now doubled. Help sits immediately left of
  the download icon on the same row at the top right. Transient timing and fixed
  probe detuning do not set the scan.
  The Time evolution sidebar starts at a native width of 420 pixels.
- **Reference → Light shift** has independent sidebar controls and settings
  files for atom, cell, static field, arbitrary beam polarization, detuning,
  and plot selection. Its main panel plots scalar, vector, tensor,
  state-resolved, eigenvalue, adjacent-transition, and scattering results.
- **Reference → Atomic polarizability** plots the complex Mathur equilibrium,
  hyperfine, vector $\alpha_{\mathrm{gt}}$, and tensor $\alpha_{\mathrm{br}}$
  response functions for D1 or D2
  light. The real curves describe phase response and the imaginary curves
  describe attenuation; diagonal vector and tensor responses are
  shown for each ground hyperfine manifold.
- **Reference → Atomic properties** presents the former Settings dialog as a
  full page with alkali isotope and temperature settings in the sidebar and thermal,
  buffer-gas, and transition-strength tabs in the main panel.

Atomic properties calculates self spin-exchange and self electron-randomization
rates at the sidebar temperature. The self electron-randomization cross section
is editable; missing reference values require a measured input. Buffer-gas
pressures are **fill pressures at 25 °C** in a sealed cell. Edit fill pressure
and electron-randomization cross section directly in the buffer-gas table;
calculated columns are read-only. Changing pressure,
temperature, or an electron-randomization cross section immediately recalculates
the table, including operating pressure, optical broadening, shift, and collision
rate. Optical coefficients in MHz/Torr and entered cross sections are held
constant; additional temperature dependence of those coefficients is not modeled.

Light-shift settings persist while navigating between pages and can be saved
or loaded with their own JSON format. Custom polarization is specified by the
beam-frame ellipse azimuth and signed ellipticity, with spherical fractions
and the dimensionless `E^(2)_0` geometry factor shown separately from shift
units. Switching the polarization input from Preset to Ellipse initializes the
ellipse angles to the equivalent current preset. The Components view provides
a multi-select dropdown for the scalar, vector, and tensor plots, includes the
degeneracy-weighted, F-independent common scalar shift and the
upper-minus-lower scalar manifold shift (labeled `ΔνF`) in the scalar panel,
and reports the vector component as the signed fictitious magnetic field
`B_fic` in µG. The Zeeman states by component view
separates vector and tensor state contributions into two plots. Eigenvalue
curves are labeled EV1, EV2, and so on.

Atomic-polarizability detuning can be referenced either to the zero-pressure
line center or to any allowed hyperfine transition on the selected D line.

Run from this directory:

```powershell
streamlit run alkali_pumping.py
```

## Dual-alkali behavior

- ACC Magnetometer and Time evolution each have independent `R_u,A` and
  `R_u,B` controls below the electron-randomization row, defaulting to 1 s⁻¹
  for new settings. Uniform relaxation acts as
  `dρ/dt = R_u [Tr(ρ) I/N − ρ]`: it drives all ground Zeeman populations
  toward equal probabilities and damps every coherence at `R_u`. It remains
  active in the dark and enters steady states, RF/optical responses,
  physical-pump feedback, relaxation totals, and CSV exports.
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
- Both analysis sidebars place magnetic-field polar angle θ, azimuth φ and
  strength in nT in one row. Angles use integer degrees: θ from +z (0–180°),
  φ from +x toward +y (0–360°). Negative strength reverses this direction.
- Each result tab has its own quantization axis, RF axis, and frequency range.
  Within that tab, the RF axis and lower/upper frequency bounds are shared by
  the atomic-moment and probe-readout plots; all other plot controls are
  independent.
- RF-A and RF-B have independent frequency ranges. In the A tab only RF-A is
  applied and the atomic observable and Probe-A readout belong to A; the B tab
  is defined conversely. New sessions default both RF upper sweep bounds to
  100 Hz.
- Within the driven species, the same RF field excites both ground manifolds
  Fa=I+1/2 and Fb=I−1/2 with their signed hyperfine g factors. Atomic plots offer
  Fa, Fb, or both overlaid. All moments use populations per total atom and the
  common upper-manifold angular-drive normalization Ωa=g_a μB B_rf/ℏ.
- Both single- and dual-species responses retain ER/SE coherence transfer
  between the two manifolds. Coherences spanning the ground hyperfine splitting
  are secularly removed. Static and RF couplings use the same signed magnetic
  constants, including the corrected negative nuclear magnetic-energy term.
- Probe-A and Probe-B have independent line, hyperfine-transition detuning,
  propagation direction, input ellipse, and path length. A Custom Weak probe is a
  weak detector only and does not pump, shift, or broaden the atoms.
- Time evolution alone offers `Preset D1a`, `Preset D1b`, `Preset D2a`, and
  `Preset D2b` as weak probe sources. Each uses direction z, path length 1 cm,
  azimuth 0°, and ellipticity angle 0°. The `a` presets use
  $F=I+1/2\rightarrow F'=I-1/2$ at −2000 MHz; the `b` presets use
  $F=I-1/2\rightarrow F'=I+1/2$ at +2000 MHz, relative to the selected
  pressure-shifted transition on D1 or D2. Selecting a preset fills these
  fields once for the chosen isotope. All fields remain editable while the
  preset is selected, and edits persist across reruns, navigation, and save/load.
- Each probe source lists both `PumpN weak` and `PumpN`. The weak item inherits
  that pump's spectrum, direction, and polarization while ignoring its
  intensity. The plain pump item represents the physical pump beam and uses
  one spatially uniform atomic density matrix driven by the path-averaged
  Stokes perturbation. The configured pump intensity is used unchanged; no
  beam-area or cell-area intensity scaling is applied. New
  sessions default Probe-A to `PumpA1 weak` and Probe-B to `PumpB1 weak`.
- Pump1 defaults to the D1 $F=I+1/2\rightarrow F'=I+1/2$ transition at
  +450 MHz and 10 µW/cm². Pump2 defaults to
  $F=I-1/2\rightarrow F'=I+1/2$ at zero detuning and 10 µW/cm². These
  relative hyperfine branches are translated when the selected isotope changes.
- Physical-pump readout includes rank-1 circular birefringence/dichroism and
  rank-2 linear birefringence/dichroism. Rather than assigning a different
  atomic state to every longitudinal slice, the solver closes the full
  dispersive and dissipative Liouvillian feedback through one cell-wide
  RF-driven density matrix. That state responds to the path-average of the
  generated fractional intensity and normalized Stokes changes. Total closes
  the coupled rank-1 plus rank-2 response; the Orientation and Alignment curves
  are separately closed counterfactuals and therefore need not add to Total.
- Probe signals include optical rotation, ellipticity, normalized Stokes
  s1/s2/s3, and fractional transmission. A single selector displays the
  orientation-induced, alignment-induced, or coherent total response; the
  total can optionally be decomposed. Its legend uses Total, Orientation, and
  Alignment columns with Amplitude, In phase, and Quadrature rows. The
  Mathur rank-2 polarizability is converted to the raw Cartesian Q_ij
  convention used by the atomic-moment plots.
- Each optical signal includes the coherent sum of Fa and Fb. Detected manifold
  selects Total, Fa, or Fb without changing the drive; Show manifold decomposition
  overlays the two contributions. Complex signals add before amplitude is taken.
  For physical-pump readout the manifold contributions share the same closure,
  while the Orientation/Alignment displays retain their counterfactual meaning.
- The Zeeman table is upstream of both field-response plots. A compact table
  lists transitions from both manifolds inside the shared RF sweep, and matching
  Fa/Fb-colored markers appear on both plots. These mark bare transition
  frequencies; coupled response peaks need not lie exactly on the markers.
- Atomic CSVs include both manifolds and their raw/plotted values. Probe CSVs
  include total, rank, and manifold contributions per nT. Exports identify the
  model and common drive normalization. Optional atomic linewidth normalization
  applies one shared reference rate, selected at the largest population
  difference across both manifolds.
- The ground-state Zeeman properties/populations, atomic response, and optical
  response are three independently collapsible, default-open result sections.
  The shared RF settings remain visible above the atomic-response section.
- The dual-alkali RF solver retains coherence feedback through self and cross
  spin exchange even though the other RF drive is zero.

Both analysis pages always list allowed hyperfine reference transitions;
the former checkbox and saved `show_allowed_only` field have been removed.
The selected reference defines detuning; calculations still include all
relevant allowed transitions, including off-resonant contributions.
ACC writes and loads condition schema **6.15**. Time evolution writes and
loads its independent settings schema **2.7** and coherent-model identifier. Time evolution requires every
current setting, rejects unknown fields, and enforces 0 ≤ t₁ ≤ t₂ ≤ 30 s
(with t₂ ≥ 0.001 s). It does not fill omitted fields or extend the time range
for previous configurations. Older analysis condition
formats are unsupported, as requested; no migration is performed.
New sessions show both atomic manifolds and default R_u,A/B to 1 s⁻¹.
The complete pre-upgrade v6.9.9 working
version is backed up in `../archive/alkali_pumping_v6.9.9/` with a SHA-256
manifest for reproducing previous calculations.

The ACC steady population model is diagonal in each selected quantization basis. When a
nonzero static field has any transverse component relative to a tab's quantization axis, the
app warns that transverse static-field mixing is omitted and does not present
that RF curve. Choose parallel or antiparallel field and quantization axes for
the supported secular weak-response calculation.

Time evolution instead evolves complete density blocks within Fa and Fb, with
independent display bases and coherent self/cross spin exchange in a common
laboratory frame. Its probe includes coherences; section CSVs contain the
displayed curves and axis metadata. Inter-F coherences and magnetic shifts in optical
detunings remain excluded. Adaptive integration is separate from display
sampling, which increases with the field and warns when capped at 20,001 samples.
See [coherent model details](docs/technical/time-evolution.qmd).

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

The unnumbered, cross-linked Quarto documentation starts at
[`docs/index.qmd`](docs/index.qmd). From the `docs` directory, render its HTML
book with `quarto render` after installing Quarto. Generated output is written
to `docs/_book/` and is intentionally excluded from Git.

## Layout

- `alkali_pumping.py`: top-level Streamlit navigation entry point.
- `alkali_pumping_app/pages/`: ACC Magnetometer, Time evolution, Light shift, Atomic polarizability, and Atomic properties pages.
- `alkali_pumping_app/physics/`: single- and dual-species numerical models.
- `alkali_pumping_app/ui/`: condition-file and table-rendering helpers.
- `tests/`: regression and physical-consistency tests.
- `docs/`: Quarto user manual, physical theory, technical notes, bibliography,
  and retained legacy derivations.
- `../archive/alkali_pumping_v6.7.13/`: application snapshot immediately before the independent dual-probe readout implementation.

Run the tests with:

```powershell
python -m unittest discover -s tests -v
```

Small ticks beneath the probe-scan detuning slider track mark all allowed
hyperfine transition centers within its full +/-10000 MHz domain. They include
the current pressure shift and use the same centers as the signal plots.
