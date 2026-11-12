"""Writing the products of a search, and reading files back.

Each writer writes one product of a ``bitvf.Result``, a table with one row
per void, overdense basin, cell or input tracer. The format of a file
follows its extension: ``.fits``, ``.fit`` and ``.fits.gz`` are FITS,
anything else is ASCII. An existing file is overwritten. Every check runs
before the file is opened, so a call that raises leaves no file. A path is a
str or an ``os.PathLike``.

Lengths are in coordinate units, which bitvf cannot know, so their columns
carry no unit and their descriptions say "in coordinate units". Index
columns hold rows of the products: ``index`` is the row of the file itself,
``voidIndex``, ``overdenseIndex`` and ``cellIndex`` the rows of the voids,
overdense basins and divergence files of the same result, -1 for none.
Counts and indices are 64-bit integers (FITS K), flags 32-bit integers
(FITS J) holding 0 or 1, and every other column a double (FITS D), written
in ASCII with 9 significant digits, or 17 in the divergence file, which then
reads back as the same double. FITS stores the doubles themselves. NaN is
written as NaN.

Every file starts with these keywords, in this order: PRODUCT (voids,
overdense_basins, divergence or galaxies); BITVFV, the bitvf version that
wrote it; GEOMETRY, box; the row count, NBASINS, NCELLS or NOBJECTS; NREC,
``Result.n_realizations``; NRECMEAN, ``Config.n_realizations_mean``; SIGMA
and SIGUNIT, the smoothing sigma and its unit, mps or coordinate; CELLSIZE,
the divergence cell size used, in units of the mps; MPS, in coordinate
units; RANDSEED, ``Result.random_seed``, when set; BOXLIM = 1, then BOXXLO,
BOXXHI, BOXYLO, BOXYHI, BOXZLO, BOXZHI, when ``Result.limits`` is set. A
number is written as the shortest decimal of at most 17 significant digits
that reads back as the same double; ``otswap.io.read_keywords`` returns the
keywords as written, in both formats.

The files of a field result (``find_voids_box_from_field``) carry no NREC
and NRECMEAN, MPS and CELLSIZE only when an mps was given, and none of the
columns that count tracers or displacements; its voids and overdense basins
files carry ``fieldCentre``. It has no galaxies file.

Only box results are written.
"""

import enum
import os
from typing import Optional, Sequence, Union

import numpy as np
import otswap
from numpy.typing import NDArray

from . import DivergenceGrid, Result

Path = Union[str, os.PathLike[str]]


class Product(enum.Enum):
    """The products, each with its file-name stem."""

    VOIDS = 0
    """Stem ``voids``."""
    OVERDENSE = 1
    """Stem ``overdense_basins``."""
    DIVERGENCE = 2
    """Stem ``divergence``."""
    GALAXIES = 3
    """Stem ``galaxies``."""


def file_name(product: Product, directory: Path, name: str, *, zspace: bool) -> str:
    """The conventional name of the file of a product.

    ``<directory>/<stem>_[zspace_]<name>``, with the stem of the product. A
    ``/`` is added only when ``directory`` is not empty and does not end in
    one; an empty directory gives ``<stem>_[zspace_]<name>``.

    Parameters
    ----------
    product : Product
        The product, which sets the stem.
    directory : str or os.PathLike
        The folder of the file; empty for none.
    name : str
        The rest of the name; its extension sets the format.
    zspace : bool
        True only for the redshift-space pass of a run with correction,
        whose files then do not collide with those of the real-space pass.

    Returns
    -------
    str

    Raises
    ------
    Error
        If ``name`` is empty, or an argument is not of its type.
    """


def write_voids(file: Path, result: Result) -> None:
    """Write the voids, one row per void (``Result.voids``).

    Columns: index (K); centreX, centreY, centreZ; radius, radiusExt;
    deltaAv; delta2Mps; edgeFlag (J); contamination; nTracers (K). A field
    result has no deltaAv, delta2Mps and nTracers, and has fieldCentre after
    contamination. Keywords as in the module description, the row count
    NBASINS. ASCII: 9 significant digits.

    Parameters
    ----------
    file : str or os.PathLike
        The file to write; its extension sets the format.
    result : bitvf.Result
        A box result.

    Raises
    ------
    Error
        If the result is not a box result, or is malformed (a result made by
        bitvf never is).
    Error
        If the file cannot be written: the message starts with
        ``"otswap: "``.
    """


def write_overdense(file: Path, result: Result) -> None:
    """Write the overdense basins, one row per basin
    (``Result.overdense``), with the columns, keywords and errors of
    ``write_voids``.

    Parameters
    ----------
    file : str or os.PathLike
        The file to write; its extension sets the format.
    result : bitvf.Result
        A box result.
    """


def write_divergence(file: Path, result: Result) -> None:
    """Write the divergence grid, one row per cell (``Result.divergence``).

    Columns: index (K); cellX, cellY, cellZ, the cell centre; divergence,
    before smoothing; divergenceSmoothed, only with sigma > 0
    (``DivergenceGrid.divergence_smoothed``); volume, in coordinate units
    cubed; weight; weightSmoothed, only with sigma > 0
    (``DivergenceGrid.weight_smoothed``); nDispl (K), the displacement ends
    counted; nTracers (K); voidIndex (K); overdenseIndex (K). Rows in the
    order of the cells, cell (i, j, k) row (i * ny + j) * nz + k. A field
    result has no nDispl and nTracers. ASCII: 17 significant digits.

    Keywords as in the module description, the row count NCELLS, followed by
    the record of the grid: LOWERX, LOWERY, LOWERZ, UPPERX, UPPERY, UPPERZ,
    CELLSIDE, in coordinate units, and NCELLX, NCELLY, NCELLZ. With MPS and
    CELLSIZE they hold every field that defines the cells, so the grid read
    back by ``read_divergence`` is the same bit for bit, and the file can be
    searched again with another smoothing
    (``find_voids_box_from_divergence``, ``find_voids_box_from_field``).

    Parameters
    ----------
    file : str or os.PathLike
        The file to write; its extension sets the format.
    result : bitvf.Result
        A box result of a search with ``Config.keep_divergence``.

    Raises
    ------
    Error
        If the result is not a box result, or holds no divergence grid.
    Error
        If the grid fails the checks of a grid to reuse, the message starting
        with ``"writeDivergence: cannot write the divergence: "``; or if
        the result is malformed, the message starting with
        ``"writeDivergence: "``.
    Error
        If the file cannot be written: the message starts with
        ``"otswap: "``.
    """


def write_galaxies(file: Path, result: Result) -> None:
    """Write where each tracer ended up, one row per input tracer, in input
    order (``Result.assignment``).

    Columns: id (K), ``Result.ids``; cellIndex (K); voidIndex (K);
    overdenseIndex (K); distance, from the centre of the tracer's basin, NaN
    for a tracer in no basin; inMask (J), only for a lightcone; inBox (J), only when
    ``Result.limits`` is set; used (J), except after
    ``find_voids_box_from_divergence``; uncorrected (J) and outsideGrid (J),
    only when those flags are filled. Keywords as in the module description,
    the row count NOBJECTS. ASCII: 9 significant digits.

    Parameters
    ----------
    file : str or os.PathLike
        The file to write; its extension sets the format.
    result : bitvf.Result
        A box result of a search with ``Config.keep_assignment``.

    Raises
    ------
    Error
        If the result is not a box result, is a field result, or holds no
        tracer assignment; or if it is malformed.
    Error
        If the file cannot be written: the message starts with
        ``"otswap: "``.
    """


class Field:
    """A field at the centres of cells, for ``find_voids_box_from_field``.
    Read-only."""

    @property
    def centres(self) -> NDArray[np.float64]:
        """The centre of each cell, in coordinate units, shape (n, 3)."""
    @property
    def values(self) -> NDArray[np.float64]:
        """The field at each centre, shape (n,)."""
    @property
    def weights(self) -> Optional[NDArray[np.float64]]:
        """The weight of each cell, shape (n,); None when no weight column
        was named, so that it passes straight to ``weights=``."""


def read_divergence(file: Path, *, angle_unit: Optional[str] = None) -> DivergenceGrid:
    """Read a divergence file written by ``write_divergence``.

    Every field the writer wrote is read back bit for bit, in both formats:
    the cells (``lower``, ``upper``, ``cell_size``, ``n_cells_axis``,
    ``n_cells``, ``mps``, ``cell_size_mps``), the box limits,
    ``n_realizations``, ``n_realizations_mean`` and every per-cell column,
    ``divergence_smoothed`` and ``weight_smoothed`` when the file has them.
    The grid passes the checks of a grid to reuse, so it can be passed to
    ``find_voids_box`` as ``grid``, or searched again with
    ``find_voids_box_from_divergence``.

    Parameters
    ----------
    file : str or os.PathLike
        A divergence file, ASCII or FITS.
    angle_unit : {"deg", "rad"}, optional
        The unit of a lightcone grid's angles: required for a lightcone
        file, an error for a box file, which holds no angles.

    Returns
    -------
    DivergenceGrid

    Raises
    ------
    Error
        If ``angle_unit`` is neither ``"deg"`` nor ``"rad"``.
    Error
        If ``angle_unit`` is given for a box file, or missing for a
        lightcone file.
    Error
        If the file cannot be read, or lacks a column: the message starts
        with ``"otswap: "``.
    Error
        If it has no PRODUCT keyword, its PRODUCT is not divergence, or
        GEOMETRY is not box (only box files are read); if a keyword of the
        grid record is missing, not a number, or not a whole number in
        range; if NREC or NRECMEAN is 0 or missing (a field result's file:
        read it with ``read_field``); if BOXLIM holds anything but 1, or
        BOXLIM and the six box keywords do not come together.
    Error
        If the rows are not NCELLS, index is not 0 to NCELLS - 1 in order, a
        count is negative, or a basin index is below -1: the message names
        the data row, counted from 1 among the data rows (comment and blank
        lines do not count), as "data row k (from 1)". If the grid fails the checks of a grid to
        reuse.
    """


def read_field(file: Path, *, centres: Sequence[str], value: str, weight: Optional[str] = None,
               delimiter: str = " ", comment: str = "#") -> Field:
    """Read a field at cell centres from a table, for
    ``find_voids_box_from_field``, in the order of the rows.

    The columns are read as ``otswap.io.read`` reads them: in FITS by name,
    from the table of HDU 2; in ASCII by name, from the last line starting
    with ``###`` before the data, or by 0-based index given as digits alone.
    Every value read must be finite; ``find_voids_box_from_field`` checks
    the rest. A bitvf divergence file is valid input::

        f = bitvf.io.read_field("divergence.fits", centres=["cellX", "cellY", "cellZ"],
                                value="divergence", weight="weight")

    Parameters
    ----------
    file : str or os.PathLike
        The table, ASCII or FITS.
    centres : sequence of three str
        The columns of x, y and z of the centres.
    value : str
        The column of the field.
    weight : str, optional
        The column of the weights; if omitted, ``Field.weights`` is None.
    delimiter : str, optional
        One character, the separator of the columns of an ASCII file;
        ignored for FITS.
    comment : str, optional
        One character, the one that starts a comment line of an ASCII file;
        ignored for FITS.

    Returns
    -------
    Field

    Raises
    ------
    Error
        If ``centres`` is not three column names, or ``delimiter`` or
        ``comment`` is not a single character.
    Error
        If the file cannot be read, a column is missing or named twice: the
        message starts with ``"otswap: "``.
    Error
        If a value read is not finite: the message names the file, the data
        row, from 1, and the column as selected.
    """


__all__: list[str]
