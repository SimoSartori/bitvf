# Lightcone executable

```{include} ../shared/executable.md
:start-after: <!-- intro-start -->
:end-before: <!-- intro-end -->
```

The commented file at the end of this page lists every key that applies to
its run.

## What `bitvf` runs

The input keys choose the kind of run; at most one of them is given.

| Run | Input key | What it does |
|---|---|---|
| pipeline | none | reads the tracers, the randoms and the mask, builds the distance table from the cosmology, and runs {cpp:func}`bitvf::runLightcone`: the reconstruction; with `perform_recons = true`, the correction of redshift-space distortions and a second reconstruction; then the search |
| from a reconstruction | `displacement_file` | reads otswap's complete displacement file, the randoms and the mask, and searches it with {cpp:func}`bitvf::findVoidsLightcone` |
| from mean displacements | `displacement_mean_file` | reads otswap's mean displacement file, the randoms and the mask, and searches it as one realization with {cpp:func}`bitvf::findVoidsLightcone`, recording the file's NRECMEAN |

A displacement file holds the tracers, the redshift cut, the cosmology and
the mps(z) profile of its reconstruction, so a run from one takes none of
the tracer, cut, cosmology and reconstruction keys. It still needs the
randoms, which set the extent of the grid and the expected counts of
`delta2Mps`, and the mask: pass those of the reconstruction.

## The keys of `bitvf`

```{include} ../shared/executable.md
:start-after: <!-- legend-start -->
:end-before: <!-- legend-end -->
```

### Run and input

| Key | Value | Default | Runs |
|---|---|---|---|
| `geometry` | `lightcone` | required | all |
| `tracer_file` | the tracers | required | pipeline |
| `column_tracer` | 3 columns, right ascension, declination and observed redshift: names, or 0-based indices in ASCII | required | pipeline |
| `id_column` | an integer column, not one of `column_tracer`, with no value twice | none: the row | pipeline |
| `units` | `degrees` or `radians`, of the right ascension and declination of every catalogue read | required | all |
| `random_file` | the randoms; they describe the survey: its footprint, its redshift distribution and every systematic. At least `Nrec` times as many as the tracers; 30 times, the minimum with the default `Nrec` of 30, is a reasonable start | required | all |
| `column_random` | 3 columns, as `column_tracer` | required | all |
| `delimiter`, `comment` | one character, for the ASCII tracer, random and bias files | `" "`, `"#"` | all |
| `veto_mask_filename` | the footprint: a full-sky HEALPix map, a pixel is observed when its value is greater than 0; a `.fits.gz` file is read as it is | required | all |
| `z_min`, `z_max` | the closed range of observed redshift of the run, either bound alone; `z_max` not below `z_min` | none: no bound | pipeline |

Tracers on unobserved pixels or outside the redshift range keep their row
and take no part; randoms there are removed.

```{include} ../shared/executable.md
:start-after: <!-- keys-run-start -->
:end-before: <!-- keys-run-end -->
```

```{include} ../shared/executable.md
:start-after: <!-- keys-files-start -->
:end-before: <!-- keys-files-end -->
```

### Cosmology

The redshift-distance table of the run, which converts the redshifts into
comoving distances; coordinates and lengths are then in Mpc/h.

| Key | Value | Default | Runs |
|---|---|---|---|
| `Omega_matter` | in (0, 1] | required | pipeline |
| `little_h` | positive | required | pipeline |
| `w0`, `wa` | the dark energy equation of state | `-1`, `0` | pipeline |

### Reconstruction

Read by otswap; [otswap's executable
page](https://simosartori.github.io/otswap/executable.html) explains them.

```{include} ../shared/executable.md
:start-after: <!-- keys-reconstruction-start -->
:end-before: <!-- keys-reconstruction-end -->
```

| Key | Value | Default | Runs |
|---|---|---|---|
| `bins` | the redshift bins of the mps(z) profile, at least 1; each bin needs enough tracers, so a small sample needs fewer bins | required | pipeline |
| `reject_mask_crossings` | `true` or `false`: reject a displacement that crosses unobserved pixels | `true` | pipeline |
| `max_unobserved_pixels_crossed` | the longest run of unobserved pixels a displacement may cross in a row; a displacement crossing a longer run is rejected | `0` | pipeline, with `reject_mask_crossings` |

### Void finder

```{include} ../shared/executable.md
:start-after: <!-- keys-voids-start -->
:end-before: <!-- keys-voids-end -->
```

### Correction of redshift-space distortions

With `perform_recons = true` the tracers are corrected in redshift,
reconstructed again and searched in real space, on the grid of the first
reconstruction. The other keys of this section apply only then.

| Key | Value | Default | Runs |
|---|---|---|---|
| `perform_recons` | `true` or `false` | `false` | pipeline |
| `bias_file` | the linear bias of the tracers against redshift: two columns, redshift and bias, read with `delimiter` and `comment` | required | pipeline, with `perform_recons` |
| `random_recons_mode` | the randoms of the second reconstruction: `resample` gives the randoms given the redshift distribution of the corrected tracers; `none` keeps them as they are; `external` takes `random_file_recons` | `resample` | pipeline, with `perform_recons` |
| `random_file_recons` | the randoms of the second reconstruction, built for the real-space sample, with the columns of `column_random` and the `units` | required | pipeline, with `random_recons_mode = external` |

```{include} ../shared/executable.md
:start-after: <!-- keys-correction-start -->
:end-before: <!-- keys-correction-end -->
```

### Output

```{include} ../shared/executable.md
:start-after: <!-- keys-output-start -->
:end-before: <!-- keys-output-end -->
```

| Key | Value | Default | Runs |
|---|---|---|---|
| `printMps` | the mps(z) profile of each reconstruction; written only when `printDisplMean` and `printDisplComplete` are both `false`, since a displacement file already holds it | `false` | pipeline |

## Output files

Written after the run, in this order, each into `output_dir` with the name
`<stem>_<output_file>.<output_format>`:

1. otswap's files, pipeline only: with `perform_recons`,
   `displacement_mean_zspace_` and `displacement_complete_zspace_` of the
   first reconstruction, or `MPS_zspace_` with `printMps`; then
   `displacement_mean_` and `displacement_complete_` of the reconstruction
   the voids come from, or `MPS_`; then `reconstructed_catalogue_`; each as
   its `print` key asks. In ASCII, each displacement file comes with two
   twin files, `<stem>.mps.dat` and `<stem>.distances.dat`, holding its
   mps(z) profile and its distance table.
2. With `save_zspace_voids`, the redshift-space pass: `voids_zspace_`,
   `overdense_basins_zspace_`, `divergence_zspace_`, `galaxies_zspace_`.
3. The voids: `voids_` always, then `overdense_basins_`, `divergence_`,
   `galaxies_`, as asked.

In ASCII, each divergence file comes with its twin `<stem>.mps.dat`, the
mps(z) profile of the grid. Their columns and keywords are on the [lightcone
file formats](formats.md) page. An existing file is overwritten.

## Exit status and messages

```{include} ../shared/executable.md
:start-after: <!-- messages-start -->
:end-before: <!-- messages-end -->
```

## The commented parameter file

{repo-file}`examples/executable/lightcone.ini` runs `bitvf` on the
lightcone tracers of the examples, with the correction of redshift-space
distortions, and searches both passes. From the repository root:

```
bitvf examples/executable/lightcone.ini
```

```{literalinclude} ../../examples/executable/lightcone.ini
:language: ini
```
