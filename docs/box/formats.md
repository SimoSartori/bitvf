# File formats

bitvf writes four products of a search, each a table, and reads three kinds
of file. The format follows the extension of the file name: `.fits`, `.fit`
and `.fits.gz` are FITS, anything else is ASCII, with the layouts of
[otswap's tables](https://simosartori.github.io/otswap/formats.html): in
ASCII, the keywords as comment lines, then a `###` line naming the columns,
then one row per line; in FITS, a binary table in HDU 2. The C++ functions
are in the namespace `bitvf::io`, the Python ones in the module `bitvf.io`;
both run the same code, so a file written from Python is the same, byte for
byte, as the one the C++ library writes.

An existing file is overwritten. Every check runs before the file is opened,
so a call that raises leaves no file. Only box results are written.

## The products

| C++ | Python | Writes | Rows |
|---|---|---|---|
| {cpp:func}`io::writeVoids <bitvf::io::writeVoids>` | {py:func}`io.write_voids <bitvf.io.write_voids>` | the voids | one per void |
| {cpp:func}`io::writeOverdense <bitvf::io::writeOverdense>` | {py:func}`io.write_overdense <bitvf.io.write_overdense>` | the overdense basins | one per basin |
| {cpp:func}`io::writeDivergence <bitvf::io::writeDivergence>` | {py:func}`io.write_divergence <bitvf.io.write_divergence>` | the divergence grid (`keepDivergence`) | one per cell |
| {cpp:func}`io::writeGalaxies <bitvf::io::writeGalaxies>` | {py:func}`io.write_galaxies <bitvf.io.write_galaxies>` | where each tracer ended up (`keepAssignment`) | one per input tracer, in input order |

The conventional file name of a product is
`<directory>/<stem>_<name>`, with the stems `voids`, `overdense_basins`,
`divergence` and `galaxies`, and `<stem>_zspace_<name>` for the
redshift-space pass of a run with correction:
{cpp:func}`bitvf::io::fileName`, {py:func}`bitvf.io.file_name`. The
executables name their files this way.

### Values

- Lengths are in coordinate units, which bitvf cannot know: their columns
  carry no unit, and their descriptions say "in coordinate units".
- Counts and indices are 64-bit integers (FITS `K`); flags are 32-bit
  integers (FITS `J`) holding 0 or 1, and so is `edgeFlag`, 0 to 2; every
  other column is a double (FITS `D`).
- ASCII writes doubles with 9 significant digits, or 17 in the divergence
  file, which then reads back as the same double. FITS stores the doubles
  themselves. NaN is written as NaN.
- `index` is the row of the file itself; `voidIndex`, `overdenseIndex` and
  `cellIndex` are rows of the voids, overdense basins and divergence files of
  the same result, -1 for none.

### Keywords

Every file starts with these keywords, in this order:

| Keyword | Value |
|---|---|
| `PRODUCT` | what the file holds: `voids`, `overdense_basins`, `divergence` or `galaxies` |
| `BITVFV` | the bitvf version that wrote it |
| `GEOMETRY` | `box` |
| `NBASINS`, `NCELLS` or `NOBJECTS` | the rows: basins, cells, or input tracers |
| `NREC` | the realizations of the search |
| `NRECMEAN` | the realizations averaged into each displacement given (`Config::nRealizationsMean`) |
| `SIGMA`, `SIGUNIT` | the smoothing sigma, and its unit, `mps` or `coordinate` |
| `CELLSIZE` | the divergence cell size used, in units of the mps |
| `MPS` | the mean particle separation, in coordinate units |
| `RANDSEED` | the seed of the randoms bitvf drew, when it drew them |
| `BOXLIM` = 1, `BOXXLO`, `BOXXHI`, `BOXYLO`, `BOXYHI`, `BOXZLO`, `BOXZHI` | the box limits of the search, in coordinate units, when it had limits |

A number is written as the shortest decimal of at most 17 significant digits
that reads back as the same double, an integer as an integer; otswap's
`io::readKeywords` (Python `otswap.io.read_keywords`) returns them as
written, in both formats.

### Voids and overdense basins

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

The [box guide](guide.md) defines each quantity.

### Divergence

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

After the common keywords, the record of the grid: `LOWERX`, `LOWERY`,
`LOWERZ`, `UPPERX`, `UPPERY`, `UPPERZ` and `CELLSIDE`, in coordinate units,
and `NCELLX`, `NCELLY`, `NCELLZ`. With `MPS` and `CELLSIZE` they hold every
value that defines the cells, so the grid read back is the same, bit for bit,
and can be reused by another search.

### Galaxies

| Column | Type | Content |
|---|---|---|
| `id` | K | the tracer's ID, or its input row from 0 when none were given |
| `cellIndex` | K | row of the cell the tracer starts in, -1 none |
| `voidIndex` | K | row of the void of the tracer, -1 none |
| `overdenseIndex` | K | row of the overdense basin of the tracer, -1 none |
| `distance` | D | distance from the centre of its basin, in coordinate units; NaN when in none |
| `inBox` | J | 0 outside the box limits, skipped; only when the search had limits |
| `used` | J | 1 when a displacement of the tracer entered the divergence; not after a search from a divergence |
| `uncorrected` | J | 1 when left unmoved by the redshift-space correction, excluded; only in the real-space pass of a run with correction |
| `outsideGrid` | J | 1 when outside the reused grid, excluded; only when a grid was reused |

### A field result

The result of a search on a field ({cpp:func}`bitvf::findVoidsBoxFromField`)
has no tracer: its files carry no `NREC` and `NRECMEAN`, `MPS` and
`CELLSIZE` only when an mps was given, and none of the columns that count
tracers or displacements (`deltaAv`, `delta2Mps`, `nTracers`, `nDispl`). Its
voids and overdense basins files add `fieldCentre`, the field at the centre,
after `contamination`. In its divergence file, `divergence` and, with
sigma > 0, `divergenceSmoothed` hold the field before and after smoothing,
`weight` the weight of the cell in the smoothing and, with sigma > 0,
`weightSmoothed` its smoothed mean. It has no galaxies file.

## Reading

| C++ | Python | Reads |
|---|---|---|
| {cpp:func}`io::readDivergence <bitvf::io::readDivergence>` | {py:func}`io.read_divergence <bitvf.io.read_divergence>` | a divergence file, back into a grid, bit for bit |
| {cpp:func}`io::readDisplacementMean <bitvf::io::readDisplacementMean>` | {py:func}`io.read_displacement_mean <bitvf.io.read_displacement_mean>` | otswap's mean displacement file of a box reconstruction |
| {cpp:func}`io::readField <bitvf::io::readField>` | {py:func}`io.read_field <bitvf.io.read_field>` | a field at cell centres, from any table |

- **Divergence file.** Read back with every value the writer wrote, in both
  formats. The grid can be reused by
  {cpp:func}`bitvf::findVoidsBox`, and searched again with another
  smoothing by {cpp:func}`bitvf::findVoidsBoxFromDivergence`. The divergence file of a
  field result is read with `readField`.
- **Mean displacement file.** The tracers, their mean displacement as one
  realization, its validity, the mps, the box limits when the file records
  them, and the IDs, ready for `findVoidsBox`. The divergence is not linear
  in the displacement, so the search on the mean of R realizations does not
  give the search on the R realizations; set `Config::nRealizationsMean` to
  the file's `NREC` so that bitvf's files record it. ASCII holds the
  positions with 9 significant digits: prefer FITS.
- **Field.** Columns are named in FITS; in ASCII, named from the last `###`
  line before the data, or given as 0-based indices in digits. Every value
  read must be finite. A divergence file is valid input, with `cellX`,
  `cellY`, `cellZ`, `divergence` and `weight`.

## Parameter files

The executables read a parameter file with otswap's syntax, described on
[otswap's file formats page](https://simosartori.github.io/otswap/formats.html):
one `key = value` per line, `#` for comments, double quotes around a value
that holds spaces. Their keys are on the [executables](../executable.md)
page.
