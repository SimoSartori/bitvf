# Shared parts of the file format pages

Not a page of its own: the box and lightcone file format pages include the
parts between the markers.

<!-- intro-start -->
bitvf writes four products of a search, each a table. The format follows the
extension of the file name: `.fits`, `.fit` and `.fits.gz` are FITS,
anything else is ASCII, with the layouts of [otswap's
tables](https://simosartori.github.io/otswap/formats.html): in ASCII, the
keywords as comment lines, then a `###` line naming the columns, then one
row per line; in FITS, a binary table in HDU 2. The C++ functions are in the
namespace `bitvf::io`, the Python ones in the module `bitvf.io`; both run
the same code, so a file written from Python is the same, byte for byte, as
the one the C++ library writes.

An existing file is overwritten. Every check runs before the file is opened,
so a call that raises leaves no file.
<!-- intro-end -->

<!-- products-start -->
| C++ | Python | Writes | Rows |
|---|---|---|---|
| {cpp:func}`io::writeVoids <bitvf::io::writeVoids>` | {py:func}`io.write_voids <bitvf.io.write_voids>` | the voids | one per void |
| {cpp:func}`io::writeOverdense <bitvf::io::writeOverdense>` | {py:func}`io.write_overdense <bitvf.io.write_overdense>` | the overdense basins | one per basin |
| {cpp:func}`io::writeDivergence <bitvf::io::writeDivergence>` | {py:func}`io.write_divergence <bitvf.io.write_divergence>` | the divergence grid (`keepDivergence`) | one per cell |
| {cpp:func}`io::writeGalaxies <bitvf::io::writeGalaxies>` | {py:func}`io.write_galaxies <bitvf.io.write_galaxies>` | where each tracer ended up (`keepAssignment`) | one per input tracer, in input order |

The conventional file name of a product is `<directory>/<stem>_<name>`, with
the stems `voids`, `overdense_basins`, `divergence` and `galaxies`, and
`<stem>_zspace_<name>` for the redshift-space pass of a run with correction:
{cpp:func}`bitvf::io::fileName`, {py:func}`bitvf.io.file_name`. The
executable names its files this way.
<!-- products-end -->

<!-- values-start -->
- Lengths are in coordinate units, which bitvf cannot know: their columns
  carry no unit, and their descriptions say "in coordinate units".
- Counts and indices are 64-bit integers (FITS `K`); flags are 32-bit
  integers (FITS `J`) holding 0 or 1, and so is `edgeFlag`, holding 0, 1 or
  2; every other column is a double (FITS `D`).
- ASCII writes doubles with 9 significant digits, or 17 in the divergence
  file, which then reads back as the same double. FITS stores the doubles
  themselves. NaN is written as NaN.
- `index` is the row of the file itself; `voidIndex`, `overdenseIndex` and
  `cellIndex` are rows of the voids, overdense basins and divergence files of
  the same result, -1 for none.
<!-- values-end -->

<!-- keywords-start -->
| Keyword | Value |
|---|---|
| `PRODUCT` | what the file holds: `voids`, `overdense_basins`, `divergence` or `galaxies` |
| `BITVFV` | the bitvf version that wrote it |
| `GEOMETRY` | the geometry of the search, `box` or `lightcone` |
| `NBASINS`, `NCELLS` or `NOBJECTS` | the rows: basins, cells, or input tracers |
| `NREC` | the realizations of the search |
| `NRECMEAN` | the realizations averaged into each displacement given: 1, or the NREC of a mean displacement file searched |
| `SIGMA`, `SIGUNIT` | the smoothing sigma, and its unit, `mps` or `coordinate` |
| `CELLSIZE` | the divergence cell size used, in units of the mps |
<!-- keywords-end -->

<!-- numbers-start -->
A number is written as the shortest decimal of at most 17 significant digits
that reads back as the same double, an integer as an integer; otswap's
`io::readKeywords` (Python `otswap.io.read_keywords`) returns them as
written, in both formats.
<!-- numbers-end -->

<!-- voids-start -->
| Column | Type | Content |
|---|---|---|
| `index` | K | row of the basin |
| `centreX`, `centreY`, `centreZ` | D | the centre, in coordinate units |
| `radius` | D | equal-volume radius of the core, in coordinate units |
| `radiusExt` | D | equal-volume radius of the extended basin, in coordinate units |
| `deltaAv` | D | mean density contrast of the tracers in the core |
| `delta2Mps` | D | density contrast within 2 mps of the centre |
| `edgeFlag` | J | contact with the boundary: 0 none, 1 extended basin only, 2 core |
| `contamination` | D | fraction of the surface on the boundary |
| `nTracers` | K | tracers starting in the core |
<!-- voids-end -->

<!-- divergence-start -->
| Column | Type | Content |
|---|---|---|
| `index` | K | row of the cell |
| `cellX`, `cellY`, `cellZ` | D | the cell centre, in coordinate units |
| `divergence` | D | the divergence of the displacement field, before smoothing |
| `divergenceSmoothed` | D | after smoothing; only with sigma > 0 |
| `volume` | D | volume of the cell, in coordinate units cubed |
| `weight` | D | tracer weight of the cell |
| `weightSmoothed` | D | tracer weight of the cell, after smoothing; only with sigma > 0 |
| `nDispl` | K | displacement ends counted in the cell |
| `nTracers` | K | tracers starting in the cell |
| `voidIndex` | K | row of the void of the cell, -1 none |
| `overdenseIndex` | K | row of the overdense basin of the cell, -1 none |
<!-- divergence-end -->

<!-- galaxies-start -->
| Column | Type | Content |
|---|---|---|
| `id` | K | the tracer's ID, or its input row from 0 when none were given |
| `cellIndex` | K | row of the cell the tracer starts in, -1 none |
| `voidIndex` | K | row of the void of the tracer, -1 none |
| `overdenseIndex` | K | row of the overdense basin of the tracer, -1 none |
| `distance` | D | distance from the centre of its basin, in coordinate units; NaN when in none |
| `used` | J | 1 when a displacement of the tracer entered the divergence |
| `uncorrected` | J | 1 when left unmoved by the redshift-space correction, excluded; only in the real-space pass of a run with correction |
| `outsideGrid` | J | 1 when outside the reused grid, excluded; only when a grid was reused, as in the real-space pass of a run with correction |
<!-- galaxies-end -->

<!-- parameter-files-start -->
The executable reads a parameter file with otswap's syntax, described on
[otswap's file formats page](https://simosartori.github.io/otswap/formats.html):
one `key = value` per line, `#` for comments, double quotes around a value
that holds spaces.
<!-- parameter-files-end -->
