# Shared parts of the executable pages

Not a page of its own: the box and lightcone executable pages include the
parts between the markers.

<!-- intro-start -->
`bitvf` runs a void search from a parameter file, with no code to write:
from the tracers on (a pipeline run), or from a file written earlier. It is
built and installed with the C++ library (`bin/bitvf` under the install
prefix); it is not in the Python package.

```
bitvf <parameter-file>
```

`--help` prints the usage, the two options and the address of the
documentation of the keys, `--version` prints `bitvf <version>`; both go to
stdout. Any other argument that starts with `--` is refused with the usage
line; a parameter file whose name starts with `--` is given as `./--name`.
Without an argument, or with more than one, the program prints the usage
line to stderr and exits with 1.

The parameter file has otswap's syntax: one `key = value` per line, `#` for
comments, double quotes around a value that holds spaces. A key left out
takes the library's default; a key with no default must be set. `key =`,
with no value, means "none" for the keys that accept it, the same as leaving
the key out. A key that is unknown, or that does not apply to the run, is an
error: a key of the other geometry is reported as "does not apply to a box"
or "does not apply to a lightcone", a key of another kind of run or whose
switch is off as "does not apply to a run ...", naming the key that decides.
<!-- intro-end -->

<!-- legend-start -->
In the tables below, "Default" is the library's, as the commented file lists
it; "none" means that the library computes or skips the value, and
"required" that the key must be set in every run it applies to. "Runs" lists
the runs the key applies to: "all", "pipeline", or the input key of the run;
"with `<key>`" names the switch the key needs. In any other run the key is
an error.
<!-- legend-end -->

<!-- keys-run-start -->
| Key | Value | Default | Runs |
|---|---|---|---|
| `verbosity` | `silent`, `normal` or `detailed`: every step of the run, the reconstruction, the correction and the search | `normal` | all |
<!-- keys-run-end -->

<!-- keys-files-start -->
| Key | Value | Default | Runs |
|---|---|---|---|
| `displacement_file` | otswap's complete displacement file, as `printDisplComplete` writes it | — | its run |
| `displacement_mean_file` | otswap's mean displacement file, as `printDisplMean` writes it | — | its run |
| `divergence_file` | a divergence file, as `printDiv` writes it, whose grid the search reuses: the same cells, and their mean particle separation | none: the search builds its own grid | `displacement_file`, `displacement_mean_file` |
<!-- keys-files-end -->

<!-- keys-reconstruction-start -->
| Key | Value | Default | Runs |
|---|---|---|---|
| `Nrec` | realizations of the reconstruction, at least 1 | `30` | pipeline |
| `threshold` | convergence, strictly between 0 and 1 | `1e-3` | pipeline |
| `seed` | unsigned; 0 draws one, recorded in the files | `0` | pipeline |
| `mesh_cellsize` | cell of the neighbour search, in mps, positive; changes the speed only | `4` | pipeline |
<!-- keys-reconstruction-end -->

<!-- keys-voids-start -->
| Key | Value | Default | Runs |
|---|---|---|---|
| `sigma` | width of the Gaussian smoothing, non-negative; 0 for none | `1` | all |
| `sigma_unit` | `mps` or `coordinate` | `mps` | all |
| `cell_size` | side of the divergence cell, in mps, positive | none: min(1, (20 / R)^(1/3)), R the realizations searched, 1 from a mean displacement file | pipeline, `displacement_file`, `displacement_mean_file`; not with `divergence_file` |
<!-- keys-voids-end -->

<!-- keys-correction-start -->
| Key | Value | Default | Runs |
|---|---|---|---|
| `rsd_smoothing_scale` | width of the neighbour average of the displacements, in coordinate units, non-negative; 0 for none. For coordinates in Mpc/h, 10 is a reasonable start | required | pipeline, with `perform_recons` |
| `rsd_weight_by_realizations` | `true` or `false`: weight each neighbour of the average by its count of valid realizations as well | `false` | pipeline, with `perform_recons` |
| `save_zspace_voids` | `true` or `false`: also search the voids of the first, redshift-space, reconstruction | `false` | pipeline, with `perform_recons` |
| `printRecCat` | `true` or `false`: write the corrected catalogue | `false` | pipeline, with `perform_recons` |
<!-- keys-correction-end -->

<!-- keys-output-start -->
| Key | Value | Default | Runs |
|---|---|---|---|
| `output_dir` | a directory, created if missing | required | all |
| `output_file` | the name in every file name, with no `/` | required | all |
| `output_format` | `dat` (ASCII) or `fits` | `dat` | all |
| `printDisplMean`, `printDisplComplete` | otswap's mean and complete displacement files of each reconstruction | `false` | pipeline |
| `printDiv` | the divergence file: the divergence and the weight of each cell, also smoothed with sigma > 0 | `false` | all |
| `printGal` | the galaxies file: where each tracer ended up | `false` | all |
| `printOverdense` | the overdense basins file | `false` | all |
<!-- keys-output-end -->

<!-- messages-start -->
A run that succeeds exits with 0. The library's report lines start with
`bitvf: `; at `normal` and `detailed`, each file written is reported as
`bitvf: wrote <path>`.

Every error is one line on stderr, `bitvf: error: <message>`, and exit
status 1. Errors in the parameter file name the file, and the line of the
key when it is given: `run.ini:7: key 'Nrec': '8.0' is not an integer`. A
missing key has no line: `run.ini: key 'output_dir' is missing`. Unknown
keys and keys that do not apply are reported all at once, one line each,
before anything is read. Apart from those of the parameter file, an error
raised by otswap starts `otswap: `. When a writer fails, the files already
written stay.
<!-- messages-end -->
