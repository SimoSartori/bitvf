# Lightcone examples

Two examples, in the lightcone geometry, on the same survey-like catalogue:
tracers between redshifts 0.885 and 1.10, a random catalogue with the same
footprint and redshift distribution, about eight randoms per tracer, and
the footprint as a HEALPix mask.

```{include} ../shared/examples.md
:start-after: <!-- intro-start -->
:end-before: <!-- intro-end -->
```

## The data

```{include} ../shared/examples.md
:start-after: <!-- data-start -->
:end-before: <!-- data-end -->
```

The mask, {repo-file}`examples/data/lightcone_mask.fits.gz`, is read as it
is: otswap reads a gzip-compressed FITS file directly.

## Building and running

```{include} ../shared/examples.md
:start-after: <!-- building-start -->
:end-before: <!-- building-end -->
```

The catalogues give the angles in degrees. The C++ programs convert them to
radians, the unit of the C++ library; the Python scripts pass them as they
are, with `angle_unit="deg"`.

## Lightcone

Reconstructs the tracers, in real space, with otswap, then finds the voids
of that reconstruction with {cpp:func}`bitvf::findVoidsLightcone`
({py:func}`bitvf.find_voids_lightcone`): two separate calls, so that the
reconstruction can be searched several times.

Reads {repo-file}`examples/data/lightcone_tracers.dat`,
{repo-file}`examples/data/lightcone_randoms.dat` and
{repo-file}`examples/data/lightcone_mask.fits.gz`.

Writes `output/voids_lightcone.dat`, `output/overdense_basins_lightcone.dat`
and `output/galaxies_lightcone.dat`; their columns and keywords are on the
[file formats](formats.md) page.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../../examples/lightcone/lightcone.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../../examples/lightcone/lightcone.cpp
:language: cpp
```
:::

::::

## Full reconstruction

Runs the whole pipeline on the same tracers taken as observed in redshift
space, in one call to {cpp:func}`bitvf::runLightcone`
({py:func}`bitvf.run_lightcone`): the reconstruction, the correction of
redshift-space distortions with the bias of
{repo-file}`examples/data/lightcone_bias.dat`, a second reconstruction in
real space with resampled randoms, and the void search, also on the
redshift-space pass.

Reads the files of the previous example and
{repo-file}`examples/data/lightcone_bias.dat`.

Writes `output/voids_zspace_lightcone.dat`,
`output/divergence_zspace_lightcone.dat`, `output/voids_lightcone.dat` and
`output/divergence_lightcone.dat`, each divergence file with its twin file
of the mps(z) profile, `.mps.dat`.
{repo-file}`examples/executable/lightcone.ini` runs the same pipeline with
the `bitvf` executable and writes the same files.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../../examples/full_recons_lightcone/full_recons_lightcone.py
:language: python
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../../examples/full_recons_lightcone/full_recons_lightcone.cpp
:language: cpp
```
:::

::::
