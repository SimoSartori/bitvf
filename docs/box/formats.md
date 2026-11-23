# Box file formats

```{include} ../shared/formats.md
:start-after: <!-- intro-start -->
:end-before: <!-- intro-end -->
```

bitvf also reads its own divergence files and fields at cell centres
([Reading](#reading)).

## The products

```{include} ../shared/formats.md
:start-after: <!-- products-start -->
:end-before: <!-- products-end -->
```

### Values

```{include} ../shared/formats.md
:start-after: <!-- values-start -->
:end-before: <!-- values-end -->
```

### Keywords

Every file starts with these keywords, in this order:

```{include} ../shared/formats.md
:start-after: <!-- keywords-start -->
:end-before: <!-- keywords-end -->
```

followed, in a box, by:

| Keyword | Value |
|---|---|
| `MPS` | the mean particle separation, in coordinate units |
| `RANDSEED` | the seed of the randoms bitvf drew, when it drew them |
| `BOXLIM` = 1, `BOXXLO`, `BOXXHI`, `BOXYLO`, `BOXYHI`, `BOXZLO`, `BOXZHI` | the box limits of the search, in coordinate units, when it had limits |

`GEOMETRY` is `box`.

```{include} ../shared/formats.md
:start-after: <!-- numbers-start -->
:end-before: <!-- numbers-end -->
```

### Voids and overdense basins

```{include} ../shared/formats.md
:start-after: <!-- voids-start -->
:end-before: <!-- voids-end -->
```

The boundary of `edgeFlag` and `contamination` is the faces of the grid;
the [method](../method.md) defines each quantity.

### Divergence

```{include} ../shared/formats.md
:start-after: <!-- divergence-start -->
:end-before: <!-- divergence-end -->
```

After the common keywords, the record of the grid: `LOWERX`, `LOWERY`,
`LOWERZ`, `UPPERX`, `UPPERY`, `UPPERZ` and `CELLSIDE`, in coordinate units,
and `NCELLX`, `NCELLY`, `NCELLZ`. With `MPS` and `CELLSIZE` they hold every
value that defines the cells, so the grid read back is the same, bit for bit,
and can be reused by another search.

### Galaxies

```{include} ../shared/formats.md
:start-after: <!-- galaxies-start -->
:end-before: <!-- galaxies-end -->
```

In a box, `inBox` follows `distance`:

| Column | Type | Content |
|---|---|---|
| `inBox` | J | 0 outside the box limits, skipped; only when the search had limits |

`used` is left out after a search from a divergence, which reads no
displacement.

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
| {cpp:func}`io::readField <bitvf::io::readField>` | {py:func}`io.read_field <bitvf.io.read_field>` | a field at cell centres, from any table |

- **Divergence file.** Read back with every value the writer wrote, in both
  formats. The grid can be reused by {cpp:func}`bitvf::findVoidsBox`, and
  searched again with another smoothing by
  {cpp:func}`bitvf::findVoidsBoxFromDivergence`. The divergence file of a
  field result is read with `readField`.
- **Mean displacement file.** bitvf has no reader for it: read it with
  `otswap::io::readDisplacementMean` (Python `otswap.io.read_displacement_mean`)
  and pass the `otswap::Result` to the `Result` overload of
  {cpp:func}`bitvf::findVoidsBox`. The box is the file's box keywords, with
  or without `BOXLIM`; a mean file without box keywords holds no box limits
  and is refused. The divergence is not linear in the displacement, so the
  search on the mean of R realizations does not give the search on the R
  realizations; set `Config::nRealizationsMean` to the `Result`'s
  `nRealizationsMean` so that bitvf's files record it. ASCII holds the
  positions with 9 significant digits: prefer FITS.
- **Field.** Columns are named in FITS; in ASCII, named from the last `###`
  line before the data, or given as 0-based indices in digits. Every value
  read must be finite. A divergence file is valid input, with `cellX`,
  `cellY`, `cellZ`, `divergence` and `weight`.

## Parameter files

```{include} ../shared/formats.md
:start-after: <!-- parameter-files-start -->
:end-before: <!-- parameter-files-end -->
```

Its keys are on the [box executables](executable.md) page.
