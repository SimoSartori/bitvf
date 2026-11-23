# Box executables

Two programs search a box from a parameter file: `bitvf`, on a tracer
catalogue, and `bitvf-field`, on a field given at the centres of the cells
of a box lattice, with no tracers ([below](#bitvf-field)).

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
| pipeline | none | reads the tracers (and the randoms, when given) and runs {cpp:func}`bitvf::runBox`: the reconstruction; with `perform_recons = true`, the correction of redshift-space distortions and a second reconstruction; then the search |
| from a reconstruction | `displacement_file` | reads otswap's complete displacement file and searches it with {cpp:func}`bitvf::findVoidsBox` |
| from mean displacements | `displacement_mean_file` | reads otswap's mean displacement file and searches it as one realization with {cpp:func}`bitvf::findVoidsBox`, recording the file's NRECMEAN; a file without box keywords holds no box limits and is refused |
| from a divergence | `divergence_input_file` | reads a divergence file written with `printDiv` and the tracers, and searches the divergence again with the smoothing given, {cpp:func}`bitvf::findVoidsBoxFromDivergence` |

The box of a run from a file is the one the file records; only a pipeline
run takes `box_min`, `box_max` and `mps`.

## The keys of `bitvf`

```{include} ../shared/executable.md
:start-after: <!-- legend-start -->
:end-before: <!-- legend-end -->
```

### Run and input

| Key | Value | Default | Runs |
|---|---|---|---|
| `geometry` | `box` | required | all |
| `tracer_file` | the tracers | required | pipeline, `divergence_input_file` |
| `column_tracer` | 3 columns, x y z in coordinate units: names, or 0-based indices in ASCII | required | pipeline, `divergence_input_file` |
| `id_column` | an integer column, not one of `column_tracer`, with no value twice | none: the row | pipeline, `divergence_input_file` |
| `random_file` | the randoms of the reconstruction | none: drawn uniformly in the box, `Nrec` times as many as the tracers inside it | pipeline |
| `column_random` | 3 columns, x y z | required | pipeline, with `random_file` |
| `delimiter`, `comment` | one character, for the ASCII tracer and random files | `" "`, `"#"` | pipeline, `divergence_input_file` |
| `divergence_input_file` | a divergence file written with `printDiv` | — | its run |

```{include} ../shared/executable.md
:start-after: <!-- keys-run-start -->
:end-before: <!-- keys-run-end -->
```

```{include} ../shared/executable.md
:start-after: <!-- keys-files-start -->
:end-before: <!-- keys-files-end -->
```

### Reconstruction

Read by otswap; [otswap's executable
page](https://simosartori.github.io/otswap/executable.html) explains them.

```{include} ../shared/executable.md
:start-after: <!-- keys-reconstruction-start -->
:end-before: <!-- keys-reconstruction-end -->
```

| Key | Value | Default | Runs |
|---|---|---|---|
| `box_min`, `box_max` | the box, three values each, finite, max above min; both or neither | none: the bounding box of the tracers, and of the randoms when given | pipeline |
| `mps` | mean particle separation, in coordinate units, positive | none: otswap's | pipeline |

### Void finder

```{include} ../shared/executable.md
:start-after: <!-- keys-voids-start -->
:end-before: <!-- keys-voids-end -->
```

### Correction of redshift-space distortions

With `perform_recons = true` the tracers are taken as redshift-space
positions: they are corrected along the line of sight, reconstructed again
and searched in real space. The other keys of this section apply only then.

| Key | Value | Default | Runs |
|---|---|---|---|
| `perform_recons` | `true` or `false` | `false` | pipeline |
| `Omega_matter` | in (0, 1] | required | pipeline, with `perform_recons` |
| `little_h` | positive | required | pipeline, with `perform_recons` |
| `w0`, `wa` | the dark energy equation of state | `-1`, `0` | pipeline, with `perform_recons` |
| `redshift` | of the box, in [0, 10] | required | pipeline, with `perform_recons` |
| `bias` | linear bias of the tracers, positive | required | pipeline, with `perform_recons` |
| `los_axis` | the line of sight, 0, 1 or 2 for x, y or z | required | pipeline, with `perform_recons` |

```{include} ../shared/executable.md
:start-after: <!-- keys-correction-start -->
:end-before: <!-- keys-correction-end -->
```

### Output

```{include} ../shared/executable.md
:start-after: <!-- keys-output-start -->
:end-before: <!-- keys-output-end -->
```

## Output files

Written after the run, in this order, each into `output_dir` with the name
`<stem>_<output_file>.<output_format>`:

1. otswap's files, pipeline only: with `perform_recons`,
   `displacement_mean_zspace_` and `displacement_complete_zspace_` of the
   first reconstruction; then `displacement_mean_` and
   `displacement_complete_` of the reconstruction the voids come from; then
   `reconstructed_catalogue_`; each as its `print` key asks.
2. With `save_zspace_voids`, the redshift-space pass: `voids_zspace_`,
   `overdense_basins_zspace_`, `divergence_zspace_`, `galaxies_zspace_`.
3. The voids: `voids_` always, then `overdense_basins_`, `divergence_`,
   `galaxies_`, as asked.

Their columns and keywords are on the [box file formats](formats.md) page.
An existing file is overwritten.

## `bitvf-field`

`bitvf-field` searches the voids of a field given at the centres of the
cells of a box lattice, with no tracers, with
{cpp:func}`bitvf::findVoidsBoxFromField`. It is built and installed with
`bitvf` (`bin/bitvf-field`), and takes the same options and the same
parameter-file syntax:

```
bitvf-field <parameter-file>
```

| Key | Value | Default |
|---|---|---|
| `verbosity` | `silent`, `normal` or `detailed` | `normal` |
| `field_file` | the table of the field, one row per cell, any row order | required |
| `column_centre` | 3 columns, x y z of the cell centres | required |
| `column_value` | the column of the field | required |
| `column_weight` | the column of the smoothing weights | none: uniform |
| `delimiter`, `comment` | one character | `" "`, `"#"` |
| `mps` | mean particle separation, positive, the unit of a sigma in mps | none |
| `sigma` | width of the Gaussian smoothing, non-negative; 0 for none | `1` |
| `sigma_unit` | `mps` or `coordinate` | `mps` |
| `output_dir`, `output_file`, `output_format` | as `bitvf` | |
| `printOverdense` | the overdense basins file | `false` |
| `printDiv` | the divergence file: the field and the weight, also smoothed with sigma > 0 | `false` |

The defaults together are an error: sigma is in mps and mps has no value. A
run sets `mps`, or `sigma_unit = coordinate`, or `sigma = 0`. The field is
split at 0, voids below it: a density, or a contrast with voids above 0, must
be shifted or flipped first. A divergence file of `bitvf` is valid input,
with `cellX cellY cellZ`, `divergence` and `weight`: the field before
smoothing, since the run smooths it, and `divergenceSmoothed` would be
smoothed twice.

It writes `voids_` always, then `overdense_basins_` and `divergence_` as
asked, with the names of `bitvf`, and no galaxies file.

## Exit status and messages

```{include} ../shared/executable.md
:start-after: <!-- messages-start -->
:end-before: <!-- messages-end -->
```

`bitvf-field` reports in the same way, with `bitvf-field` in place of
`bitvf` in its own lines.

## The commented parameter files

{repo-file}`examples/executable/box.ini` runs `bitvf` on the halos of the
examples, as the [full reconstruction example](examples.md) does, and writes
the same voids file; {repo-file}`examples/executable/field.ini` runs
`bitvf-field` on the field of the field example. From the repository root:

```
bitvf examples/executable/box.ini
bitvf-field examples/executable/field.ini
```

::::{tab-set}

:::{tab-item} bitvf

```{literalinclude} ../../examples/executable/box.ini
:language: ini
```
:::

:::{tab-item} bitvf-field

```{literalinclude} ../../examples/executable/field.ini
:language: ini
```
:::

::::
