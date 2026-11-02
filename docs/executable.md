# The executables

Two programs run a search from a parameter file, with no code to write:

- `bitvf` searches the voids of a tracer catalogue, from the reconstruction
  on, or from files written earlier;
- `bitvf-field` searches the voids of a field given at the centres of the
  cells of a box lattice, with no tracers.

Both are built and installed with the C++ library (`bin/bitvf` and
`bin/bitvf-field` under the install prefix); they are not in the Python
package. Only box runs are implemented.

```
bitvf <parameter-file>
bitvf-field <parameter-file>
```

`--help` prints the usage, the two options and the address of this page,
`--version` prints `bitvf <version>` (`bitvf-field <version>`); both go to
stdout. Any other argument that starts with `--` is refused with the usage
line; a parameter file whose name starts with `--` is given as `./--name`.
Without an argument, or with more than one, the program prints the usage
line to stderr and exits with 1.

The parameter file has otswap's syntax: one `key = value` per line, `#` for
comments, double quotes around a value that holds spaces. A key left out
takes the library's default; a key with no default must be set. `key =`,
with no value, means "none" for the keys that accept it, the same as leaving
the key out. A key that is unknown, or that does not apply to the run, is an
error. Two commented files, one per program, are at the end of this page,
with every key that applies to their run.

## What `bitvf` runs

The input keys choose the kind of run; at most one of them is given.

| Run | Input key | What it does |
|---|---|---|
| pipeline | none | reads the tracers (and the randoms, when given), runs {cpp:func}`bitvf::runBox`: the reconstruction, with `perform_recons = true` the correction of redshift-space distortions and a second reconstruction, then the search |
| from a reconstruction | `displacement_file` | reads otswap's complete displacement file and searches it, {cpp:func}`bitvf::findVoidsBox` |
| from mean displacements | `displacement_mean_file` | reads otswap's mean displacement file, {cpp:func}`bitvf::io::readDisplacementMean`, and searches it as one realization |
| from a divergence | `divergence_input_file` | reads a divergence file written with `printDiv` and `printDivUnsmoothed` and the tracers, and searches it again with the smoothing given, {cpp:func}`bitvf::findVoidsBoxFromDivergence` |

## The keys of `bitvf`

**R** required; **o** optional, absent taking the default; **—** not a key
of that run, an error when given. The runs: **P** pipeline, **D**
`displacement_file`, **M** `displacement_mean_file`, **V**
`divergence_input_file`. "Default" is the library's, as the commented file
lists it; "none" means the library computes or skips the value.

### Run

| Key | Value | Default | P | D | M | V |
|---|---|---|---|---|---|---|
| `geometry` | `box` or `lightcone` | — | R | R | R | R |
| `verbosity` | `silent`, `normal` or `detailed`: every step of the run, the reconstruction, the correction and the search | `normal` | o | o | o | o |

### Input

| Key | Value | Default | P | D | M | V |
|---|---|---|---|---|---|---|
| `tracer_file` | path | — | R | — | — | R |
| `column_tracer` | 3 columns, x y z: names, or 0-based indices in ASCII | — | R | — | — | R |
| `id_column` | an integer column, not one of `column_tracer` | none: the row | o | — | — | o |
| `random_file` | path | none: drawn in the box | o | — | — | — |
| `column_random` | 3 columns | — | R with `random_file` | — | — | — |
| `delimiter`, `comment` | one character | `" "`, `"#"` | o | — | — | o |
| `displacement_file` | otswap's complete displacement file | — | — | R | — | — |
| `displacement_mean_file` | otswap's mean displacement file | — | — | — | R | — |
| `divergence_input_file` | a divergence file with `divergenceUnsmoothed` | — | — | — | — | R |
| `divergence_file` | a divergence file whose grid is reused | none: the search builds its own | — | o | o | — |

### Reconstruction

Read by otswap; see [otswap's executable
page](https://simosartori.github.io/otswap/executable.html) for their
meaning. Pipeline only.

| Key | Value | Default |
|---|---|---|
| `Nrec` | realizations, at least 1 | `30` |
| `threshold` | convergence, strictly between 0 and 1 | `1e-3` |
| `seed` | unsigned; 0 draws one, recorded in the files | `0` |
| `mesh_cellsize` | neighbour-search cell, in mps, positive; speed only | `4` |
| `box_min`, `box_max` | the box, three values each, max above min; both or neither | none: the bounding box of the tracers, and of the randoms when given |
| `mps` | mean particle separation, positive | none: otswap's |

### Void finder

| Key | Value | Default | P | D | M | V |
|---|---|---|---|---|---|---|
| `sigma` | smoothing width, non-negative; 0 for none | `1` | o | o | o | o |
| `sigma_unit` | `mps` or `coordinate` | `mps` | o | o | o | o |
| `cell_size` | divergence cell, in mps, positive | none: min(1, (20 / R)^(1/3)), R the realizations searched, 1 from a mean file | o | o, not with `divergence_file` | o, not with `divergence_file` | — |

### Correction of redshift-space distortions

Pipeline only. With `perform_recons = true` the tracers are taken as
redshift-space positions; the other keys apply only then.

| Key | Value | Default |
|---|---|---|
| `perform_recons` | `true` or `false` | `false` |
| `Omega_matter` | in (0, 1] | — (required) |
| `little_h` | positive | — (required) |
| `w0`, `wa` | the dark energy equation of state | `-1`, `0` |
| `redshift` | of the box, in [0, 10] | — (required) |
| `bias` | linear bias of the tracers, positive | — (required) |
| `los_axis` | line of sight, 0, 1 or 2 for x, y or z | — (required) |
| `rsd_smoothing_scale` | width of the neighbour average, in coordinate units, non-negative; for coordinates in Mpc/h, 10 is a reasonable start | — (required) |
| `rsd_weight_by_realizations` | `true` or `false` | `false` |
| `save_zspace_voids` | also search the redshift-space pass | `false` |
| `printRecCat` | write the corrected catalogue | `false` |

### Output

| Key | Value | Default | P | D | M | V |
|---|---|---|---|---|---|---|
| `output_dir` | a directory, created if missing | — | R | R | R | R |
| `output_file` | the name in every file name, with no `/` | — | R | R | R | R |
| `output_format` | `dat` (ASCII) or `fits` | `dat` | o | o | o | o |
| `printDisplMean`, `printDisplComplete` | otswap's mean and complete displacement files | `false` | o | — | — | — |
| `printDiv` | the divergence file | `false` | o | o | o | o |
| `printDivUnsmoothed` | adds the divergence before smoothing; only with `printDiv = true` | `false` | o | o | o | o |
| `printGal` | the galaxies file | `false` | o | o | o | o |
| `printOverdense` | the overdense basins file | `false` | o | o | o | o |

The keys of a lightcone run (`units`, `random_recons_mode`,
`random_file_recons`, `veto_mask_filename`, `z_min`, `z_max`, `bins`,
`reject_mask_crossings`, `max_unobserved_pixels_crossed`, `bias_file`,
`printMps`) are refused in a box run, as "does not apply to a box".

## The keys of `bitvf-field`

| Key | Value | Default |
|---|---|---|
| `verbosity` | `silent`, `normal` or `detailed` | `normal` |
| `field_file` | the table of the field, one row per cell, any row order | — (required) |
| `column_centre` | 3 columns, x y z of the cell centres | — (required) |
| `column_value` | the column of the field | — (required) |
| `column_weight` | the column of the smoothing weights | none: uniform |
| `delimiter`, `comment` | one character | `" "`, `"#"` |
| `mps` | mean particle separation, positive, the unit of a sigma in mps | none |
| `sigma` | smoothing width, non-negative; 0 for none | `1` |
| `sigma_unit` | `mps` or `coordinate` | `mps` |
| `output_dir`, `output_file`, `output_format` | as `bitvf` | |
| `printOverdense` | the overdense basins file | `false` |
| `printDiv` | the divergence file: the field after smoothing | `false` |
| `printDivUnsmoothed` | adds the field before smoothing; only with `printDiv = true` | `false` |

The defaults together are an error: sigma is in mps and mps has no value. A
run sets `mps`, or `sigma_unit = coordinate`, or `sigma = 0`. The field is
split at 0, voids below it: a density, or a contrast with voids above 0, must
be shifted or flipped first. A divergence file of `bitvf` is valid input,
with `cellX cellY cellZ`, `divergenceUnsmoothed` and `weight`: the field
before smoothing, since the run smooths it, and the smoothed `divergence`
would be smoothed twice.

## Output files

Written after the run, in this order, each into `output_dir` with the name
`<stem>_<output_file>.<output_format>`:

1. otswap's files, `bitvf` pipeline only: with `perform_recons`,
   `displacement_mean_zspace_` and `displacement_complete_zspace_` of the
   first reconstruction; then `displacement_mean_` and
   `displacement_complete_` of the reconstruction the voids come from; then
   `reconstructed_catalogue_`; each as its `print` key asks.
2. With `save_zspace_voids`, the redshift-space pass: `voids_zspace_`,
   `overdense_basins_zspace_`, `divergence_zspace_`, `galaxies_zspace_`.
3. The voids: `voids_` always, then `overdense_basins_`, `divergence_`,
   `galaxies_`, as asked. `bitvf-field` writes no galaxies file.

Their columns and keywords are on the [file formats](box/formats.md) page.
An existing file is overwritten.

## Exit status and messages

A run that succeeds exits with 0. The library's report lines start with
`bitvf: `; at `normal` and `detailed`, each file written is reported as
`bitvf: wrote <path>` (`bitvf-field: wrote <path>`).

Every error is one line on stderr, `bitvf: error: <message>` (or
`bitvf-field: error: <message>`), and exit status 1. Errors in the parameter
file name the file, and the line of the key when it is given:
`run.ini:7: key 'Nrec': '8.0' is not an integer`. Unknown keys and keys that
do not apply are reported all at once, one line each, before anything is
read. Apart from those of the parameter file, an error raised by otswap
starts `otswap: `. When a writer fails, the files already written stay.

## The commented parameter files

{repo-file}`examples/executable/box.ini` runs `bitvf` on the halos of the
examples, as the [full reconstruction example](box/examples.md) does, and
writes the same voids file;
{repo-file}`examples/executable/field.ini` runs `bitvf-field` on the field
of the field example. From the repository root:

```
bitvf examples/executable/box.ini
bitvf-field examples/executable/field.ini
```

::::{tab-set}

:::{tab-item} bitvf

```{literalinclude} ../examples/executable/box.ini
:language: ini
```
:::

:::{tab-item} bitvf-field

```{literalinclude} ../examples/executable/field.ini
:language: ini
```
:::

::::
