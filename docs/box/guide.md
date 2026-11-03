# Box guide

A box is a volume with Cartesian coordinates and a single mean particle
separation, such as a simulation snapshot or a cube cut from one. This page
covers what applies to every box search, in both APIs; the
[reference](../cpp.rst) gives each function in full.

## Which function

| Starting from | C++ | Python |
|---|---|---|
| a tracer catalogue: reconstruction, optional correction of redshift-space distortions, search | {cpp:func}`bitvf::runBox` | {py:func}`bitvf.run_box` |
| an otswap box reconstruction | {cpp:func}`bitvf::findVoidsBox` | {py:func}`bitvf.find_voids_box` |
| displacement arrays, from any source, or otswap's mean displacement file | {cpp:func}`bitvf::findVoidsBox`, {cpp:func}`bitvf::io::readDisplacementMean` | {py:func}`bitvf.find_voids_box`, {py:func}`bitvf.io.read_displacement_mean` |
| a divergence grid kept or written earlier, to search it again with another smoothing | {cpp:func}`bitvf::findVoidsBoxFromDivergence` | {py:func}`bitvf.find_voids_box_from_divergence` |
| any field at the centres of a box lattice, with no tracers | {cpp:func}`bitvf::findVoidsBoxFromField` | {py:func}`bitvf.find_voids_box_from_field` |
| a field on a grid, split into basins only | {cpp:func}`bitvf::watershedBox` | {py:func}`bitvf.watershed_box` |

Each search returns a result holding the voids and the overdense basins,
with what it counted; the divergence grid and where each tracer ended up are
kept only when asked for (`keepDivergence`, `keepAssignment`). Files are
written only by the writers of `bitvf::io` (Python `bitvf.io`).

## What a search computes

1. **Divergence.** On a grid of cubic cells, each valid displacement, from
   a tracer to the random matched to it, that leaves its tracer's cell adds
   its outward component through the face it crosses, in units of the cell
   side, to the cell it starts in and to the cell it ends in. A point on a
   face between two cells belongs to the upper one; a displacement that
   starts or ends on a face, an edge or a vertex, or runs along a face,
   counts as that of the same tracer moved one ulp into its cell. The
   divergence of a cell is that sum over the displacement ends counted in
   it, its weight that count over the number of realizations. A cell that
   no displacement reaches takes the divergence of its neighbours and keeps
   a weight of 0.
2. **Smoothing.** The divergence is averaged with a Gaussian of width sigma,
   each cell weighted by its tracer weight.
3. **Watershed.** Each cell flows to the neighbour, among the 26 around it,
   where the divergence drops the most; a cell with no lower neighbour is a
   seed, and a basin is a seed with every cell that flows to it. Seeds below
   0 start voids, the others overdense basins. Every cell belongs to exactly
   one basin.
4. **Catalogue.** For each basin: a centre, two radii, two density
   contrasts, its contact with the faces of the grid, and the tracers in it.

## Coordinates and units

- Cartesian coordinates and every length (the mps, box limits,
  displacements, sigma, radii) are in coordinate units: the unit of the
  tracers given. bitvf assumes none, and its files give lengths no unit.
- In C++ the arrays are flat and row-major, `[object][xyz]`, with
  displacements `[realization][tracer][xyz]`; the [array
  layout](../cpp.rst) section of the reference lists every array of a result.
- In Python a coordinate array has shape (N, 3), converted to float64. Every
  per-cell array of a box grid has shape (nx, ny, nz), C order, so that
  `a[i, j, k]` is cell (i, j, k). Flags are `bool`.

## Mean particle separation and cell size

The mean particle separation (mps) sets the scale of the search: the cell
size, the smoothing in units of mps, the sphere of `delta2Mps` and the mean
density of the density contrasts. bitvf takes it from the reconstruction,
`otswap::Result::mps`, or from the caller; it never computes its own.

The cell side is `divergenceCellSize` (Python `divergence_cell_size`) times
the mps. Left empty, it is min(1, (20 / R)^(1/3)) mps, R the number of
realizations: 1 mps up to 20 realizations, smaller above. A search on a mean
displacement file counts as one realization, so its default is 1. Smaller
cells resolve smaller voids but hold fewer displacements each, so the field
is noisier. The value used is in `Result::config`.

## The grid

The grid covers the box limits of the search; on an otswap reconstruction,
the box it used (`otswap::Result::box`); otherwise the tracers and the
endpoints of their valid displacements. It is padded by 0.1 mps on each
side. It has the same number of cells on the three axes, so it may overhang
a box that is not a cube. It is not periodic: a cell on a face has fewer
neighbours, and a basin touching the faces is flagged.

A grid can be reused: the `grid` argument of `findVoidsBox` takes one from
an earlier result (`Result::divergence`), from
{cpp:func}`bitvf::boxGrid` ({py:func}`bitvf.box_grid`) or from a divergence
file ({cpp:func}`bitvf::io::readDivergence`). The search then uses its cells
and its mps, and excludes the tracers outside it (`Result::outsideGrid`).
Several searches on one grid give catalogues on the same cells.

## Box limits

Box limits, given to the search or recorded in the reconstruction, act like
a mask: a tracer outside them keeps its row in every per-tracer output but is
skipped, flagged in `TracerAssignment::inBox`, and the grid is built on the
limits. A point on a face is inside. Pass the limits the reconstruction used:
a valid displacement ending outside the grid is an error.

## Smoothing

`Config::smoothing` sets the width sigma of the Gaussian and its unit, mps
(the default) or coordinate units. Each cell takes the mean of the cells
within 3.5 sigma along each axis, weighted by their tracer weight times
exp(-r² / (2 sigma²)). A larger sigma merges small basins into larger ones;
sigma = 0 disables the smoothing. 1 mps, the default, is a reasonable
starting value.

With sigma > 0 the grid keeps, beside the divergence and the weight, the
smoothed divergence (`divergenceSmoothed`), the field the watershed splits,
and the smoothed weight (`weightSmoothed`), the mean of the weight over the
same cells with the Gaussian alone as weight. The smoothed weight shows how
much of a cell's smoothed divergence comes from displacements: near 0 where
the cells around were mostly filled from their neighbours. The watershed
does not read it.

{cpp:func}`bitvf::findVoidsBoxFromDivergence` searches a divergence again
with another sigma, from a grid kept with `keepDivergence` or read from a
divergence file, without the reconstruction.

## The catalogues

Voids and overdense basins have the same columns. Row b is basin b, in the
order of its seed cell, whatever the number of threads.

- **Core and extended basin.** The core of a void is the part where the
  divergence is below 0; the extended void is the whole basin. For an
  overdense basin the two coincide.
- **centre**: the seed cell's centre x₀ shifted by
  −Σₙ divₙ (xₙ − x₀) / Σₙ |divₙ| over the seed's neighbouring cells, and
  kept inside the seed cell.
- **radius**, **radiusExt**: the radius of the sphere with the volume of the
  core, of the extended basin.
- **deltaAv**: the tracers starting in the core over its volume, divided by
  the mean density mps⁻³, minus 1.
- **delta2Mps**: the density contrast of the tracers within 2 mps of the
  centre.
- **edgeFlag**: 0 when the basin does not touch the faces of the grid, 1
  when only the extended basin does, 2 when the core does.
- **contamination**: the fraction of the basin's surface on the faces of the
  grid.
- **nTracers**: the tracers starting in the core.

The [file formats](formats.md) page lists the columns as written.

## The pipeline

{cpp:func}`bitvf::runBox` runs everything on a tracer catalogue:

1. checks every argument, so that an invalid run fails before any work;
2. fixes the box: the limits given, or the bounding box of the tracers (and
   of the randoms, when given);
3. without randoms, draws R × N of them uniformly in the box, N the tracers
   inside it, from the seed `config.ot.seed + 1299709` (a seed of 0 counting
   as the one drawn to replace it), recorded in `Result::randomSeed`;
4. reconstructs with otswap and searches the voids.

With `BoxConfig::rsd` the tracers are taken as redshift-space positions:
after the first reconstruction they are corrected along the line of sight
with otswap's correction, reconstructed again with the same randoms and
seed, and searched in real space on the grid of the first pass, with its
mps. Tracers the correction leaves unmoved, or moves outside that grid, keep
their row and are flagged (`Result::uncorrected`, `Result::outsideGrid`).
`BoxRsd::zspaceVoids` also searches the first, redshift-space,
reconstruction. The correction needs a distance table carrying the growth
rate, the redshift of the box, the linear bias and the width of the
neighbour average: 10 Mpc/h is a reasonable starting value for coordinates
in Mpc/h, to check in each analysis.

## Reproducibility and threads

Given the same inputs and a seed other than 0, every result is the same bit
for bit whatever the number of threads, on every supported platform and
compiler, but for the elapsed time. A seed of 0 is replaced by a drawn one,
recorded in otswap's result, so any run can be repeated. The number of
threads follows `OMP_NUM_THREADS`.

## What it prints

`Config::verbosity` sets the report, at otswap's levels: `Silent`, nothing;
`Normal`, what the box limits or a reused grid left out and one line per
call with its time; `Detailed`, also one line per step of the search. Each
line starts with `bitvf: ` and is written as soon as its step finishes, to
`Config::log` in C++ (`std::clog` by default), to `sys.stdout` in Python.
Every count printed can be read from the result. In Python, Ctrl+C stops a
search at its next report line.

## Errors

Every failure raises `bitvf::Error` (Python `bitvf.Error`), with a message
that names the quantity and the rule broken; an error raised by otswap
inside a bitvf call has a message starting `otswap: `. Arguments are checked
before any work. `bitvf::Error` is unrelated to `otswap::Error`: catching one
does not catch the other.
