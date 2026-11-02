# Box examples

Three examples, each a C++ program and a Python script that do the same
thing: the same inputs, parameters and seed, and the same output files, byte
for byte. Each folder in {repo-tree}`examples` is a CMake project of its own,
built against an installed bitvf.

## The data

The examples read the files of {repo-tree}`examples/data`, which come with a
clone of the repository; they are not in the Python package or the source
distribution. They are provided for running the examples, may be
redistributed unchanged with them, and may not be used for any other purpose:
{repo-file}`examples/data/README.md` gives the terms and the content of each
file.

## Building and running

For an example `<name>`, with bitvf installed under `<prefix>`:

```
cmake -S examples/<name> -B build/examples/<name> -DCMAKE_PREFIX_PATH=<prefix>
cmake --build build/examples/<name>
build/examples/<name>/<name> [data_dir]
python examples/<name>/<name>.py [data_dir]
```

Both read the examples' data folder unless `data_dir` is given, and write to
`examples/<name>/output/`, created if missing. Where cfitsio is not in a
system path, give its prefix with `-DCFITSIO_ROOT=...` (and `-DZLIB_ROOT=...`
for a static cfitsio). The field example reads the output of
`full_recons_box`, which is run first.

## Box

Reconstructs halos in a cube, in real space, with otswap, then finds the
voids of that reconstruction with {cpp:func}`bitvf::findVoidsBox`
({py:func}`bitvf.find_voids_box`): two separate calls, so that the
reconstruction can be searched several times.

Reads {repo-file}`examples/data/box_halos_real_space.dat`.

Writes `output/voids_box.dat`, `output/overdense_basins_box.dat` and
`output/galaxies_box.dat`; their columns and keywords are on the [file
formats](formats.md) page.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../../examples/box/box.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../../examples/box/box.cpp
:language: cpp
```
:::

::::

## Full reconstruction

Runs the whole pipeline on the same halos observed in redshift space, in one
call to {cpp:func}`bitvf::runBox` ({py:func}`bitvf.run_box`): the
reconstruction, the correction of redshift-space distortions, a second
reconstruction in real space and the void search, also on the
redshift-space pass. The randoms are drawn by bitvf.

Reads {repo-file}`examples/data/box_halos_redshift_space.dat`.

Writes `output/voids_zspace_box.dat`, `output/divergence_zspace_box.dat`,
`output/voids_box.dat` and `output/divergence_box.dat`.
{repo-file}`examples/executable/box.ini` runs the same pipeline with the
`bitvf` executable and writes the same files.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../../examples/full_recons_box/full_recons_box.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../../examples/full_recons_box/full_recons_box.cpp
:language: cpp
```
:::

::::

## Field

Finds the voids of a field given at the centres of the cells of a box
lattice, with no tracers, with {cpp:func}`bitvf::io::readField` and
{cpp:func}`bitvf::findVoidsBoxFromField` ({py:func}`bitvf.io.read_field`,
{py:func}`bitvf.find_voids_box_from_field`). The field is the divergence
before smoothing, `divergenceUnsmoothed`, of the full reconstruction's
divergence file, smoothed more widely: the search smooths the field it
reads, so the smoothed `divergence` column would be smoothed twice.

Reads `examples/full_recons_box/output/divergence_box.dat`, or the file given
as argument in place of `data_dir`.

Writes `output/voids_field.dat` and `output/overdense_basins_field.dat`.
{repo-file}`examples/executable/field.ini` runs the same search with the
`bitvf-field` executable, on the divergence file of `box.ini`, and writes the
same files.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../../examples/field/field.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../../examples/field/field.cpp
:language: cpp
```
:::

::::
