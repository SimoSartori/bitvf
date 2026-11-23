# Box guide

The box geometry is the one of a simulation snapshot, or of a cube cut from
one: tracers in a rectangular volume, in comoving Cartesian coordinates,
with a single mean particle separation (mps). This page explains how to
run a search in this geometry, from the C++ library, from Python and with
the executables, how to choose its parameters and how to read what it
returns. [The method](../method.md) explains what each step computes and
why; the [C++](../cpp.rst) and [Python](../python.rst) references give
every function in full.

## Which function

| Starting from | C++ | Python |
|---|---|---|
| a tracer catalogue: reconstruction, optional correction of redshift-space distortions, search | {cpp:func}`bitvf::runBox` | {py:func}`bitvf.run_box` |
| an otswap box reconstruction, kept in memory or read from otswap's complete or mean displacement file | {cpp:func}`bitvf::findVoidsBox`, the `otswap::Result` overload | {py:func}`bitvf.find_voids_box` with the `otswap.Result` |
| displacement arrays, from any source | {cpp:func}`bitvf::findVoidsBox` | {py:func}`bitvf.find_voids_box` |
| a divergence grid kept or written earlier, to search it again with another smoothing | {cpp:func}`bitvf::findVoidsBoxFromDivergence` | {py:func}`bitvf.find_voids_box_from_divergence` |
| any field at the centres of a box lattice, with no tracers | {cpp:func}`bitvf::findVoidsBoxFromField` | {py:func}`bitvf.find_voids_box_from_field` |
| a field on a grid, split into basins only | {cpp:func}`bitvf::watershedBox` | {py:func}`bitvf.watershed_box` |
| a grid to share among several searches | {cpp:func}`bitvf::boxGrid` | {py:func}`bitvf.box_grid` |

The pipeline is the simplest start. Calling otswap and the search
separately keeps the reconstruction, so that it can be written, or searched
again with other parameters without being computed again. The [box
examples](examples.md) show both, and the [box
executables](executable.md) run them from a parameter file.

```{include} ../shared/guide.md
:start-after: <!-- results-start -->
:end-before: <!-- results-end -->
```

## A typical run

The runs below are the pipeline on the halos of a simulation box observed
in redshift space: otswap reconstructs them with 8 realizations and a fixed
seed, drawing the randoms uniformly in the box, the redshift-space
distortions along the y axis are corrected, and the voids are searched on
the second, real-space, reconstruction and on the first as well. They are
the core of the [full reconstruction example](examples.md#full-reconstruction),
which also locates the data folder `data` and summarizes and writes the
result; the example page shows it whole.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../../examples/full_recons_box/full_recons_box.py
:language: python
:start-at: "# Comoving Cartesian coordinates in Mpc/h, shape (N, 3)"
:end-at: run = bitvf.run_box(
:dedent: 4
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../../examples/full_recons_box/full_recons_box.cpp
:language: cpp
:start-at: "// Comoving Cartesian coordinates in Mpc/h, x, y, z for each halo"
:end-at: const bitvf::BitVFResult run = bitvf::runBox(
:dedent: 4
```
:::

::::

The real-space voids are in `run.voids`, the redshift-space ones in
`run.zspaceVoids` (Python `run.zspace_voids`). Without `rsd`, the pipeline
reconstructs once and searches that reconstruction. To search a
reconstruction kept from otswap instead, pass it to `findVoidsBox`
(`find_voids_box`), as the [box example](examples.md#box) does.

The executable runs the same pipeline from a parameter file, with no code
to write:

```
bitvf examples/executable/box.ini
```

The [box executables](executable.md) page lists its keys, and the commented
file {repo-file}`examples/executable/box.ini` sets each of them.

## Coordinates and units

- Cartesian coordinates and every length (the mps, the box limits, the
  displacements, the smoothing width, the radii) are in coordinate units:
  the unit of the tracers given. bitvf assumes none, and its files give
  lengths no unit.
- In C++ the arrays are flat and row-major, `[object][xyz]`, with
  displacements `[realization][tracer][xyz]`; the [array
  layout](../cpp.rst) section of the reference lists every array of a
  result.
- In Python a coordinate array has shape (N, 3), converted to float64.
  Every per-cell array of a box grid has shape (nx, ny, nz), C order, so
  that `a[i, j, k]` is cell (i, j, k). Flags are `bool`.

## The box and its limits

The box of a pipeline is the one given in `BoxConfig::limits` (Python
`limits`) or, without it, the smallest box containing the tracers and, when
given, the randoms. It is fixed once, before the reconstruction, and shared
by every later step. For a simulation snapshot, give the box of the
simulation.

Box limits, given to the search or recorded in the reconstruction, act like
a mask: a tracer outside them keeps its row in every per-tracer output but
takes no part, flagged in `TracerAssignment::inBox` (Python `in_box`). A
point on a face is inside. When searching displacement arrays, pass the
limits the reconstruction used: a valid displacement ending outside the
grid is an error.

## Mean particle separation and cell size

The mps sets the scale of the search: the cell size, the smoothing in mps,
the sphere of `delta2Mps` and the mean density of the density contrasts.
bitvf takes it from the reconstruction, `otswap::Result::mps`, or from the
caller (`BoxConfig::mps`); it never computes its own.

The cell side is `divergenceCellSize` (Python `divergence_cell_size`) times
the mps. Left empty, it is min(1, (20 / R)^(1/3)) mps, R the number of
realizations: 1 mps up to 20 realizations, smaller above. A search on a
mean displacement file counts as one realization, so its default is 1.
Smaller cells resolve smaller voids but hold fewer displacements each, so
the divergence is noisier. The value used is in `Result::config`.

## The grid

The grid is made of cubic cells covering the box limits of the search,
those of the reconstruction (`otswap::Result::box`) or, for displacement
arrays without limits, the tracers and the endpoints of their valid
displacements. It is padded by 0.1 mps on each side, and has the same number
of cells on the three axes, so it may overhang a box that is not a cube. It
is not periodic, even for a periodic simulation: a cell on a face has fewer
neighbours, and a basin touching the faces is flagged.

Each cell has 26 neighbours, the cells sharing a face, an edge or a corner
with it. A point on a face between two cells belongs to the upper one.

A grid can be reused: the `grid` argument of `findVoidsBox` takes one from
an earlier result (`Result::divergence`), from `boxGrid` (Python
`box_grid`) or from a divergence file (`bitvf::io::readDivergence`, Python
`bitvf.io.read_divergence`). The search then uses its cells and its mps,
and excludes the tracers outside it, flagged in `Result::outsideGrid`.
Several searches on one grid give catalogues on the same cells.

## Smoothing

`Config::smoothing` (Python `smoothing`) sets the width of the Gaussian
and its unit, mps (the default) or coordinate units; 1 mps, the default, is
a reasonable start, and a width of 0 turns the smoothing off. In a box, each
cell takes the mean of the cells within 3.5 widths along each axis, a cube
rather than a sphere.

`findVoidsBoxFromDivergence` (Python `find_voids_box_from_divergence`)
searches a divergence again with another width, from a grid kept with
`keepDivergence` or read from a divergence file, without the
reconstruction.

## The correction of redshift-space distortions

With `BoxConfig::rsd` (Python `rsd`), `runBox` takes the tracers as
redshift-space positions: after the first reconstruction they are moved
along the line of sight with otswap's correction, reconstructed again with
the same randoms and seed, and searched in real space on the grid of the
first pass, with its mps. [The
method](../method.md#6-correction-of-redshift-space-distortions) explains
the correction.

`bitvf::BoxRsd` (Python `bitvf.BoxRsd`) sets the axis of the line of sight,
the redshift of the box and the linear bias of the tracers, all fixed by
the simulation, and the width of the average of the displacements over
neighbouring tracers: 10 Mpc/h is a reasonable start for coordinates in
Mpc/h, to be checked in each analysis. The distance table given to the run
carries the growth rate. `zspaceVoids` (Python `zspace_voids`) also
searches the first, redshift-space, reconstruction.

Tracers the correction leaves unmoved, or moves outside the grid, keep
their rows and take no part in the second pass, flagged in
`Result::uncorrected` and `Result::outsideGrid`.

## Reading and checking the results

The catalogues of voids and overdense basins have the same columns, one row
per basin, in the order of the seed cells, whatever the number of threads.
[The method](../method.md#5-voids-and-overdense-regions) defines each one.
In Python, the columns are NumPy arrays with snake-case names, such as
`result.voids.radius`, `radius_ext`, `delta_av`, `delta_2mps`, `edge_flag`
and `contamination`; [names side by side](../naming.md) gives each C++
name.

`edgeFlag` marks the basins that touch the faces of the grid: 1 when only
the extended basin does, 2 when the core does. The flow is not known beyond
the faces, so for statistics that need whole voids keep those with
`edgeFlag` 0, or 0 and 1. With `keepDivergence`, the smoothed weight of
each cell (`weightSmoothed`, Python `weight_smoothed`) shows how much of
its smoothed divergence comes from displacements: near 0, it comes mostly
from cells that no displacement reached, and should be read with care.

## Writing the results

```{include} ../shared/guide.md
:start-after: <!-- writing-start -->
:end-before: <!-- writing-end -->
```

The [box file formats](formats.md) page lists every column and keyword.

## Reproducibility and threads

```{include} ../shared/guide.md
:start-after: <!-- reproducibility-start -->
:end-before: <!-- reproducibility-end -->
```

A box pipeline given no randoms draws them uniformly in the box, as many as
the realizations times the tracers inside it, from a seed derived from the
reconstruction's, recorded in `Result::randomSeed` (Python `random_seed`)
and in the files as `RANDSEED`.

## What it prints

```{include} ../shared/guide.md
:start-after: <!-- prints-start -->
:end-before: <!-- prints-end -->
```

At Normal verbosity a box search prints the tracers left out by the box
limits or a reused grid, when there are limits or a grid, and a last line
with the number of voids and overdense basins and the time; `bitvf::Config`
lists every line.

## Errors

```{include} ../shared/guide.md
:start-after: <!-- errors-start -->
:end-before: <!-- errors-end -->
```
