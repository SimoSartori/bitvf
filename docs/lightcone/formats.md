# Lightcone file formats

```{include} ../shared/formats.md
:start-after: <!-- intro-start -->
:end-before: <!-- intro-end -->
```

bitvf also reads its own divergence files back, as a grid to reuse
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

- Angles are in degrees in every file, whatever the unit of the arrays of
  the call. The redshift of a centre is that of its comoving distance in
  the distance table of the search.
- Coordinates are Cartesian, with the observer at the origin, in the
  distance unit of the distance table: Mpc/h when it comes from otswap's
  cosmology.

### Keywords

Every file starts with these keywords, in this order:

```{include} ../shared/formats.md
:start-after: <!-- keywords-start -->
:end-before: <!-- keywords-end -->
```

followed, in a lightcone, by:

| Keyword | Value |
|---|---|
| `RANDSEED` | the seed of the randoms bitvf drew, when it drew them: the redshifts of the resampled randoms of a second pass |
| `ZCUTMIN`, `ZCUTMAX` | the bounds of the redshift cut of the search, each only when finite |
| `NRANDUSE` | the randoms the search used: kept by the mask and the redshift cut, and inside the grid |
| `MASKNSID` | the NSIDE of the mask of the search |

`GEOMETRY` is `lightcone`. `CELLSIZE` is in units of mps(z), the mean
particle separation at the redshift of each shell; the mps(z) profile itself
is in the divergence file.

```{include} ../shared/formats.md
:start-after: <!-- numbers-start -->
:end-before: <!-- numbers-end -->
```

### Voids and overdense basins

```{include} ../shared/formats.md
:start-after: <!-- voids-start -->
:end-before: <!-- voids-end -->
```

In a lightcone, the sky position of the centre follows `index`:

| Column | Type | Content |
|---|---|---|
| `centreRA`, `centreDec` | D | right ascension and declination of the centre, in degrees |
| `centreRed` | D | redshift of the centre |

The 2 mps of `delta2Mps` are 2 mps(z) at the centre's redshift. The
boundary of `edgeFlag` and `contamination` is the masked cells and the
border of the grid; the [method](../method.md) defines each quantity.

### Divergence

```{include} ../shared/formats.md
:start-after: <!-- divergence-start -->
:end-before: <!-- divergence-end -->
```

In a lightcone, the sky position and the limits of each cell follow
`cellZ`, and `masked` comes last:

| Column | Type | Content |
|---|---|---|
| `cellRA`, `cellDec` | D | right ascension and declination of the cell centre, in degrees |
| `cellRed` | D | redshift of the cell centre |
| `cellDist` | D | comoving distance of the cell centre, in coordinate units |
| `cellRAMin`, `cellRAMax` | D | right-ascension limits of the cell, in degrees |
| `cellDecMin`, `cellDecMax` | D | declination limits of the cell, in degrees |
| `cellDistMin`, `cellDistMax` | D | inner and outer distance of the cell, in coordinate units |
| `masked` | J | 1 when the cell centre is on an unobserved pixel of the mask |

The limits of the cells describe them; they are rebuilt from the record of
the grid when the file is written, and not read back. Right-ascension limits
are continuous across 0 within a ring, so they may exceed 360; a ring of a
single sector spans 360 degrees. The central sphere, when the grid has one,
spans the distances from 0 and the whole sky; a polar cap reaches
declination 90 (or -90) over the whole circle.

After the common keywords, the record of the grid: `DISTMIN`, `DISTMAX`,
the range of comoving distance of the tracers and randoms, in coordinate
units, and `RAMIN`, `RAMAX`, `DECMIN`, `DECMAX`, the footprint, in degrees.
With `CELLSIZE`, the mps(z) profile and the distance table of the search,
they define the cells: a search with the same distance table rebuilds the
same grid from them.

The mps(z) profile of the grid follows in a second table: HDU `MPSPROF` of
a FITS file, or, in ASCII, the twin file `<stem>.mps.<ext>`, the divergence
file's name with `.mps` before its extension
(`divergence_lightcone.mps.dat` beside `divergence_lightcone.dat`). Its
columns, with 17 significant digits:

| Column | Type | Content |
|---|---|---|
| `redshift` | D | centre of the redshift bin |
| `MPS` | D | mean particle separation in the bin, in coordinate units |
| `nTracers` | J | tracers in the bin |

and its keywords: `PRODUCT` = `divergence mps profile`, `BITVFV`, `ZMIN` and
`ZMAX` (the lowest and highest redshift binned) and `MPSREPR` (the
count-weighted median mps, in coordinate units). It has no `GEOMETRY`.

A grid is written only if it reads back as a grid to reuse: it must be
rebuilt from its record and the distance table of the result with the same
number of cells, which costs one grid build.

### Galaxies

```{include} ../shared/formats.md
:start-after: <!-- galaxies-start -->
:end-before: <!-- galaxies-end -->
```

In a lightcone, `inMask` and `inCut` follow `distance`, and `outsideCut`
comes last:

| Column | Type | Content |
|---|---|---|
| `inMask` | J | 0 on an unobserved pixel of the mask, skipped |
| `inCut` | J | 0 outside the redshift cut, skipped; only when the cut has a finite bound |
| `outsideCut` | J | 1 when moved outside the redshift cut by the correction, excluded; only in the real-space pass of a run with correction |

A tracer both on an unobserved pixel and outside the cut has 0 in both
flags. A tracer that starts in a masked cell keeps that cell in `cellIndex`
but is in no basin.

## Reading

| C++ | Python | Reads |
|---|---|---|
| {cpp:func}`io::readDivergence <bitvf::io::readDivergence>` | {py:func}`io.read_divergence <bitvf.io.read_divergence>` | a divergence file, back into a grid |

A lightcone divergence file, with its profile, gives a lightcone grid: the
record, the mps(z) profile, the cell centres and their sky positions, the
fields, the counts, the basin indices and `masked`. The distances and
`CELLSIZE` read back bit for bit; the angles equal those of the grid
written up to the rounding of their conversion from degrees, so a tracer
within rounding of a right-ascension or declination edge can change cell.
The limits of the cells are not read. In Python, `angle_unit` (`"deg"` or
`"rad"`) is required and sets the unit of the grid's angles.

The grid can be reused by {cpp:func}`bitvf::findVoidsLightcone`, with the
distance table of the search that wrote it, to search other catalogues on
the same cells. There is no search of a lightcone from a divergence file.

## Parameter files

```{include} ../shared/formats.md
:start-after: <!-- parameter-files-start -->
:end-before: <!-- parameter-files-end -->
```

Its keys are on the [lightcone executable](executable.md) page.
