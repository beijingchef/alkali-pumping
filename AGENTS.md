# Alkali Pumping project guidance

## Project orientation

- Treat this directory as the active application and Git repository.
- Read `README.md`, the newest entries in `CHANGELOG.md`, and
  `alkali_pumping_app/version.py` before making a substantial change.
- The Streamlit entry point is `alkali_pumping.py`. Page implementations live
  in `alkali_pumping_app/pages/`, reusable UI/state helpers in
  `alkali_pumping_app/ui/`, numerical models in `alkali_pumping_app/physics/`,
  and regression/physical-consistency tests in `tests/`.
- Treat `../archive/` and `../streamlit_top_nav_test/` as historical or
  experimental material. Do not edit or copy code from them unless the user
  explicitly asks or a comparison is needed to diagnose a regression.

## Working style

- Inspect the current implementation and relevant tests before proposing or
  applying a change. Preserve established behavior unless the request clearly
  changes it.
- Prefer focused, maintainable edits over broad rewrites. Reuse existing
  helpers and conventions before introducing a new abstraction.
- When a request is ambiguous, infer intent from the current UI, physics model,
  tests, README, and changelog. Ask only when the choice would materially alter
  scientific meaning, saved-data compatibility, or user-visible behavior.
- Do not add a production dependency without explaining why it is needed and
  obtaining confirmation.
- Never commit secrets, machine-specific paths, virtual environments, caches,
  or generated temporary files.

## Scientific and numerical requirements

- Treat physical correctness, sign conventions, normalization, axes, and units
  as part of the public behavior. Do not silently change them.
- Keep angular-frequency versus ordinary-frequency conversions explicit.
  Preserve the documented positive/negative-frequency, RF phase, spherical
  tensor, and Cartesian quadrupole conventions.
- Preserve population normalization, Hermiticity where required, conservation
  properties, and dimensional consistency. Add or update tests for these
  invariants whenever their implementation changes.
- Keep single-species and dual-species behavior consistent. Verify that an
  inactive Alkali B and zero-coupling limits recover the corresponding simpler
  model.
- Maintain the distinction between perturbing pump beams and non-perturbing
  probe readout. Do not let probe settings alter pumping, light shifts,
  broadening, or steady-state populations.
- When implementing a physics change, state the governing assumption or
  equation in code comments or `docs/` when the reasoning would not be obvious
  to a future maintainer.

## Streamlit UI and state

- Preserve the top-navigation structure and page-specific persistent settings.
  Verify navigation does not reset unrelated controls.
- Follow the existing compact scientific-dashboard layout, but prefer native
  Streamlit widget behavior and spacing unless custom styling is explicitly
  requested.
- Keep labels, symbols, subscripts, signs, displayed units, plot legends,
  captions, tables, downloads, and visible defaults synchronized with the
  underlying calculation.
- Changes to plotted values must also be reflected in CSV exports and their
  units/metadata where applicable.
- Preserve independent Alkali A/B, RF-A/RF-B, Pump, and Probe controls unless a
  request explicitly makes them shared.

## Saved conditions and compatibility

- Treat condition-file formats as versioned public interfaces.
- When changing a saved field, default, or structure, update the appropriate
  schema version and provide migration from every supported older format.
- Existing saved conditions must retain their stored values and meaning when
  loaded. New-session default changes must not overwrite values loaded from a
  file.
- Add round-trip, missing-field, invalid-version, and migration tests for
  condition-schema changes.

## Verification

- Run the focused tests for the modified area while iterating.
- Before handing off a completed code change, run the full suite from this
  directory:

  ```powershell
  python -m unittest discover -s tests -v
  ```

- For UI or navigation changes, also launch or exercise the Streamlit app when
  practical and check for runtime exceptions and state regressions:

  ```powershell
  streamlit run alkali_pumping.py
  ```

- Do not describe a change as complete if relevant tests fail. Report any test
  that could not be run and the reason.

## Versions and documentation

- For user-visible behavior or physics-model changes, update `CHANGELOG.md` and
  keep `alkali_pumping_app/version.py` consistent with the intended release.
- Update `README.md` when workflows, defaults, supported behavior, setup, or
  project layout changes.
- Keep durable design decisions and scientific conventions in versioned
  documentation rather than relying on chat history or local Codex memory.
- Use concise commit-ready summaries that distinguish UI changes, physics
  changes, compatibility changes, and verification performed.

## Git safety

- Preserve unrelated user changes in the working tree.
- Do not commit, push, rewrite history, delete branches, or modify archived
  snapshots unless explicitly requested.
- Before work moves to another computer, ensure all intended code,
  documentation, tests, and this `AGENTS.md` are committed and pushed.
