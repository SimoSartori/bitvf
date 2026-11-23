# Changelog

## 0.1.0 — YYYY-MM-DD

The first release.

### Added

- The void finder in box geometry: `findVoidsBox` (Python `find_voids_box`) on
  displacement arrays or on an otswap reconstruction, with optional box
  limits and a reused grid; `findVoidsBoxFromDivergence` on a divergence
  computed earlier; `findVoidsBoxFromField` on any field given at the centres
  of a box lattice; `boxGrid` and `watershedBox`.
- The box pipeline `runBox` (Python `run_box`): otswap's reconstruction,
  optionally the correction of redshift-space distortions and a second
  reconstruction in real space, and the void search, with optional voids on
  the redshift-space pass.
- Writers of the voids, overdense basins, divergence grid and per-tracer
  assignment, in ASCII or FITS, and readers of a divergence file, of otswap's
  mean displacement file and of a field at cell centres (`bitvf::io`, Python
  `bitvf.io`).
- The void finder and the pipeline in lightcone geometry: `findVoidsLightcone`
  (Python `find_voids_lightcone`) on sky positions or on an otswap
  reconstruction, `findVoidsLightconeFromCartesian`, `lightconeGrid`, and
  `runLightcone`, with the correction of redshift-space distortions and a
  second reconstruction in real space; the divergence grid of a lightcone is
  written and read back by `bitvf::io`, and the executable `bitvf` runs it.
  In Python, sky arrays and the sky fields of the results follow the
  `angle_unit` of the call.
- The executables `bitvf` and `bitvf-field`, driven by a parameter file.
- Progress reports at otswap's verbosity levels, each line printed as its
  step finishes.
- Results identical bit for bit whatever the number of threads, on every
  supported platform and compiler.
