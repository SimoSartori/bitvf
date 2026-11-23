# Lightcone guide

The lightcone geometry is the one of a galaxy survey or of a mock catalogue
built like one: tracers observed on a part of the sky, each with a right
ascension, a declination and a redshift, seen by an observer at the origin.
This page explains how to run a search in this geometry, from the C++
library, from Python and with the executable, how to choose its parameters
and how to read what it returns. [The method](../method.md) explains what
each step computes and why; the [C++](../cpp.rst) and [Python](../python.rst)
references give every function in full.

## Which function

| Starting from | C++ | Python |
|---|---|---|
| a tracer catalogue: reconstruction, optional correction of redshift-space distortions, search | {cpp:func}`bitvf::runLightcone` | {py:func}`bitvf.run_lightcone` |
| an otswap lightcone reconstruction, kept in memory or read from otswap's complete or mean displacement file | {cpp:func}`bitvf::findVoidsLightcone`, the `otswap::Result` overload | {py:func}`bitvf.find_voids_lightcone` with the `otswap.Result` |
| sky positions and displacement arrays, from any source | {cpp:func}`bitvf::findVoidsLightcone` | {py:func}`bitvf.find_voids_lightcone` |
| Cartesian positions, the observer at the origin, and displacement arrays | {cpp:func}`bitvf::findVoidsLightconeFromCartesian` | {py:func}`bitvf.find_voids_lightcone_from_cartesian` |
| a grid to share among several searches | {cpp:func}`bitvf::lightconeGrid` | {py:func}`bitvf.lightcone_grid` |

The pipeline is the simplest start. Calling otswap and the search
separately keeps the reconstruction, so that it can be written, or searched
again with other parameters without being computed again.

```{include} ../shared/guide.md
:start-after: <!-- results-start -->
:end-before: <!-- results-end -->
```

## The inputs

A lightcone search needs five inputs: the tracers, the randoms, the mask,
the cosmology and the parameters. All of them except the parameters
describe the survey, and should be those of the analysis, not chosen for
the void search.

**Tracers.** One row per tracer: right ascension, declination and observed
redshift. The C++ library takes the angles in radians, as a flat array
`[tracer][ra, dec, z]`; `otswap::skyToRadians` converts a catalogue read in
degrees. Python takes an array of shape (N, 3) in the unit named by
`angle_unit`, `"deg"` or `"rad"`, which every lightcone call requires, and
returns the sky arrays of its result in the same unit. Any finite right
ascension is accepted; a declination must lie between −90 and 90 degrees.
Identities can be given as `ids`, one 64-bit integer per tracer; they are
copied to the result and to the files, and otherwise the row number stands
for them.

**Randoms.** A random catalogue sampling the survey volume, in the same
format as the tracers. The randoms are required: they are the only
description of the survey that bitvf has, so they must follow its
footprint, its redshift distribution and every systematic, such as
variations of completeness on the sky. bitvf applies no weight: whatever the
survey selects must be in the randoms. Each realization of the
reconstruction uses as many randoms as there are tracers, so at least the
number of realizations times the tracers are needed. 30 randoms per tracer,
the minimum with otswap's default of 30 realizations, is a reasonable
start. The search reads them too, to set the extent of the grid and the
expected counts of `delta2Mps`, which are less noisy with more randoms.

**Mask.** An `otswap::Mask` (Python `otswap.Mask`): a full-sky HEALPix map
of the footprint, in RING or NESTED order, read from a FITS file, also
gzip-compressed, or built from an array; a pixel is observed when its value
is greater than 0. Tracers and randoms on unobserved pixels are left out,
and so is every cell of the grid whose centre falls on one. Give the search
the mask of the reconstruction, or a stricter one.

**Cosmology.** An `otswap::DistanceTable` (Python `otswap.DistanceTable`)
converts redshifts into comoving distances: use the fiducial cosmology of
the analysis. Built from Ωm, h, w0 and wa of a flat cosmology, it is
tabulated from z = 0 to 10 by default, gives distances in Mpc/h, and
carries the growth rate the correction of redshift-space distortions needs.
Every length bitvf returns, radii and centres included, is in the unit of
the table. The table must cover the redshifts of the tracers and randoms
kept, and, with the correction, the corrected ones, which may fall slightly
outside the catalogue's range: the default range leaves room for that.

**Redshift cut.** The closed range of observed redshift searched, given to
the reconstruction and to the search. Tracers outside keep their rows, with
a flag, and take no part; randoms outside are removed. By default nothing is
cut: set a cut to search a slice of the survey.

## A typical run

The runs below reconstruct the tracers with otswap and search the voids of
the reconstruction, with 8 realizations and a fixed seed. They are the core
of the [lightcone example](examples.md), which also locates the data folder
`data` and summarizes and writes the result; the example page shows it
whole.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../../examples/lightcone/lightcone.py
:language: python
:start-at: output.mkdir(parents=True, exist_ok=True)
:end-at: result = bitvf.find_voids_lightcone(
:lines: 3-
:dedent: 4
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../../examples/lightcone/lightcone.cpp
:language: cpp
:start-at: std::filesystem::create_directories(output);
:end-at: const bitvf::Result result = bitvf::findVoidsLightcone(
:lines: 4-
:dedent: 4
```
:::

::::

The search reads everything else from the reconstruction: the tracers,
their displacements, the mps(z) profile, the distance table and the
redshift cut. The randoms and the mask are given again, since the
reconstruction does not keep them.

The pipeline, {cpp:func}`bitvf::runLightcone` ({py:func}`bitvf.run_lightcone`),
does both steps in one call, and adds the correction of redshift-space
distortions when asked ([below](#the-correction-of-redshift-space-distortions));
the [full reconstruction example](examples.md#full-reconstruction) runs it.

The executable runs the pipeline, or a search of otswap's displacement
files, from a parameter file, with no code to write:

```
bitvf examples/executable/lightcone.ini
```

The [lightcone executable](executable.md) page lists its keys, and the
commented file {repo-file}`examples/executable/lightcone.ini` sets each of
them.

## Choosing the main parameters

Few parameters change the voids found, and each has a reasonable starting
value. [The method](../method.md) explains each step they control.

- **Realizations** (`otswap::Config::nRealizations`, Python
  `n_realizations`; executable `Nrec`): more realizations give a less noisy
  divergence and allow smaller cells, at the cost of time. otswap's default
  of 30 is a reasonable start; a few realizations are enough to try a
  setup.
- **Redshift bins of the mps(z) profile** (`LightconeConfig::nBins`, Python
  `n_bins`; executable `bins`): no default; 30 is a reasonable start. Each
  bin must hold enough tracers for its mean separation to be measured, so
  use fewer bins for a small sample.
- **Cell size** (`Config::divergenceCellSize`, Python
  `divergence_cell_size`; executable `cell_size`): in units of mps(z), left
  empty for the default min(1, (20 / R)^(1/3)), R the realizations, that is
  1 up to 20 realizations and smaller above. Smaller cells resolve smaller
  voids and follow the footprint more closely, but each holds fewer
  displacements. The value used is in `Result::config`.
- **Smoothing** (`Config::smoothing`, Python `smoothing`; executable `sigma`
  and `sigma_unit`): a Gaussian of width 1 mps by default, a reasonable
  start. A larger width merges small basins into larger ones; 0 turns the
  smoothing off.
- **Paths across the mask** (`otswap::Config::rejectCrossings` and
  `maxUnobservedPixelsCrossed`, Python `reject_crossings` and
  `max_unobserved_pixels_crossed`; executable `reject_mask_crossings` and
  `max_unobserved_pixels_crossed`): whether the reconstruction rejects the
  displacements whose path on the sky crosses unobserved pixels, and how
  many in a row a path may cross. Rejection with no pixel allowed is a
  reasonable start: a displacement across a hole of the footprint carries
  no information on the flow there.

## The grid

Without a grid given, the search builds one that covers the tracers and
randoms kept: shells in comoving distance from the observer, each divided
into rings of declination and each ring into sectors of right ascension
over the footprint. A shell is the cell size times mps(z) wide at its
centre, so the cells grow with redshift and hold a similar number of
tracers at every distance. When the sample comes close to the observer, the
first shell is a sphere around it; near a pole, the last ring of a shell is
a cap over the whole circle.

A cell whose centre falls on an unobserved pixel is *masked*: it has no
divergence and belongs to no basin, and a displacement that starts or ends
in one is not counted. A tracer on an observed pixel inside a masked cell
is in no basin; the result counts these tracers in
`Counts::tracersInMaskedCells` (Python `tracers_in_masked_cells`), and the
report prints them. Many such tracers mean that the cells are coarse for the
footprint: a smaller cell size makes them fewer.

The search reports the grid it built at Normal verbosity: its cells, its
shells, the cell size and the ranges of distance, right ascension and
declination. The whole grid, with the divergence of every cell, is kept
with `keepDivergence`.

## The correction of redshift-space distortions

The observed redshift of a tracer includes its peculiar velocity, which
distorts the voids along the line of sight. `runLightcone` corrects it when
`LightconeConfig::rsd` is set (Python `rsd`; executable
`perform_recons = true`): after the first reconstruction, each tracer is
moved along its line of sight, in redshift only, to its estimated
real-space position, the corrected tracers are reconstructed again with the
same seed, and the voids are searched on this second reconstruction, on the
grid of the first. [The method](../method.md#6-correction-of-redshift-space-distortions)
explains the correction.

`bitvf::LightconeRsd` (Python `bitvf.LightconeRsd`) sets it:

- `bias`, no default: the linear bias of the tracers as a table in
  redshift, interpolated linearly between its nodes; read it from a file
  with columns z and b with `otswap::io::readBiasTable` (Python
  `otswap.io.read_bias_table`). It must be positive at the redshift of
  every tracer kept.
- `sigma`, no default: the width of the average of the displacements over
  the neighbours of each tracer, in the unit of the distance table. 10 Mpc/h
  is a reasonable start; the best value depends on the sample, and should
  be checked in each analysis.
- `secondPassRandoms` (Python `second_pass_randoms`): the randoms of the
  second reconstruction, which should follow the real-space sample.
  `Resample`, the default, keeps the sky positions of the first randoms and
  draws their redshifts from the corrected tracers, so that the randoms
  follow the corrected redshift distribution; it is the reasonable start.
  `None` reuses the first randoms unchanged, keeping the observed redshift
  distribution. `External` takes a random catalogue of your own built for
  the real-space sample, passed as `randomsSkySecondPass` (Python
  `randoms_sky_second_pass`); it must leave at least the realizations times
  the corrected tracers on observed pixels, inside the cut and inside the
  grid.
- `zspaceVoids` (Python `zspace_voids`): also search the first,
  redshift-space, reconstruction, on the same grid, and return it in
  `BitVFResult::zspaceVoids`. Comparing the two catalogues shows how much
  the distortions change the voids of the sample.

A tracer the correction cannot use keeps its row in the real-space result
and takes no part in the second reconstruction. Each such row carries one
flag, which says why: `uncorrected` when the correction did not move it,
because no tracer with a valid displacement lies around it or because the
mask or the cut left it out of the first pass; `outsideCut` when it was moved outside the redshift
cut; `outsideGrid` when it was moved outside the grid of the first pass.
The counts of each are in `Result::counts`, and the example prints them.

With `keepReconstructions` (Python `keep_reconstructions`), the run also
returns both otswap reconstructions and the corrected catalogue, so that
they can be written with otswap's writers.

## Reading and checking the results

The catalogues of voids and overdense basins have the same columns, one row
per basin, in the order of the seed cells, whatever the number of threads.
[The method](../method.md#5-voids-and-overdense-regions) defines each one;
in the lightcone geometry the centre is given both in Cartesian coordinates,
the observer at the origin, and as right ascension, declination and
redshift (`centreSky`, Python `centre_sky`). In Python, the columns are
NumPy arrays with snake-case names, such as `result.voids.radius`,
`radius_ext`, `delta_av`, `delta_2mps`, `edge_flag` and `contamination`;
[names side by side](../naming.md) gives each C++ name.

A few checks are worth making on every run.

- **What was left out.** `Result::counts` gives the tracers on unobserved
  pixels, outside the cut, in masked cells and without a valid
  displacement, and the randoms removed; the report prints the same counts.
  With `keepAssignment`, the flags `inMask` and `inCut` (Python `in_mask`,
  `in_cut`) mark each tracer, and `used` those whose displacements entered
  the divergence.
- **Voids at the boundary.** A basin that touches a masked cell or the
  border of the grid, at the edge of the footprint, at a hole of the mask
  or at the redshift limits, is flagged in `edgeFlag`: 1 when only the
  extended basin touches it, 2 when the core does. The flow is not known
  beyond the boundary, so the size and the shape of such a void are less
  certain. For statistics that need whole voids, such as a size function,
  keep the voids with `edgeFlag` 0, or 0 and 1; `contamination`, the
  fraction of the basin's surface on the boundary, grades the same effect.
- **`delta2Mps`.** It compares the tracers within 2 mps(z) of the centre
  with the randoms there, so its noise depends on the number of randoms per
  tracer; it is NaN when no random lies in the sphere, which happens near
  the edges of a sparse random catalogue.
- **The smoothed weight.** With `keepDivergence`, the grid holds, beside
  the divergence of each cell, its smoothed value and the smoothed weight
  (`weightSmoothed`, Python `weight_smoothed`). Where the smoothed weight is
  close to 0, the smoothed divergence comes mostly from cells that no
  displacement reached, and should be read with care.

## Reusing a grid

Searches on one grid give catalogues on the same cells, which makes them
directly comparable: the two passes of the pipeline share one this way, and
so can the mocks of a survey and the survey itself. A grid is reused by
passing it as `grid` to `findVoidsLightcone`; it comes from:

- `lightconeGrid` (Python `lightcone_grid`), which builds, without the
  search, the grid a search of the same reconstruction, or of the same
  tracers and randoms, would build;
- the result of an earlier search kept with `keepDivergence`
  (`Result::divergence`);
- a divergence file, read with `bitvf::io::readDivergence` (Python
  `bitvf.io.read_divergence`). In ASCII, the mps(z) profile of the grid is
  in a twin file, `<stem>.mps.<ext>`, which must stay beside it.

A reused grid brings its mps(z) profile, which then sets the smoothing in
mps and the densities in place of the reconstruction's. Leave the cell size
of the configuration empty, and use the same distance table that built the
grid: the cells are rebuilt from it, and a different table is refused only
when it gives a different number of cells. Tracers outside the grid keep
their rows, flagged in `outsideGrid`, and randoms outside are removed. A
grid read from a file stores its edges with 17 significant digits, so a
tracer within rounding of a cell edge can change cell.

## Writing the results

```{include} ../shared/guide.md
:start-after: <!-- writing-start -->
:end-before: <!-- writing-end -->
```

Lightcone files give angles in degrees, whatever the unit of the call, and
lengths in the unit of the distance table. In ASCII, a divergence file comes
with its twin file of the mps(z) profile. The [lightcone file
formats](formats.md) page lists every column and keyword.

## Reproducibility and threads

```{include} ../shared/guide.md
:start-after: <!-- reproducibility-start -->
:end-before: <!-- reproducibility-end -->
```

The randoms resampled for the second pass of a pipeline are drawn from a
seed derived from the reconstruction's, recorded in `Result::randomSeed`
(Python `random_seed`) of the real-space result and in the files as
`RANDSEED`.

## What it prints

```{include} ../shared/guide.md
:start-after: <!-- prints-start -->
:end-before: <!-- prints-end -->
```

At Normal verbosity a lightcone search prints the tracers and randoms kept
by the mask and the cut, the grid, the tracers in masked cells when there
are any, and a last line with the number of voids and overdense basins and
the time:

```
bitvf: kept <k> of <N> tracers: <m> on unobserved pixels, <c> outside the redshift cut
bitvf: kept <k> of <M> randoms: <m> on unobserved pixels, <c> outside the redshift cut
bitvf: grid of <C> cells in <S> shells, cell size <s> mps (<source>), distances [<d0>, <d1>], RA [<a0>, <a1>] deg, Dec [<b0>, <b1>] deg (<from>) in <t> s
bitvf: <n> tracers start in a masked cell, in no basin
bitvf: findVoidsLightcone: <V> voids and <O> overdense basins from <N> tracers in <T> s; <W> without a valid displacement
```

A pipeline adds its own lines on the second pass and a last line for the
whole run; {cpp:func}`bitvf::runLightcone` lists them.

## Errors

```{include} ../shared/guide.md
:start-after: <!-- errors-start -->
:end-before: <!-- errors-end -->
```
