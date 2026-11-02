"""Writing the products of a search, and reading files back.

The format of a file follows its extension: ``.fits``, ``.fit`` and
``.fits.gz`` are FITS, anything else is ASCII. An existing file is
overwritten; a call that raises leaves no file. A path is a str or an
``os.PathLike``.
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
    OVERDENSE = 1
    DIVERGENCE = 2
    GALAXIES = 3


def file_name(product: Product, directory: Path, name: str, *, zspace: bool) -> str:
    """``<directory>/<stem>_[zspace_]<name>``; the extension, part of
    ``name``, sets the format. ``zspace`` only for the redshift-space pass of
    a run with correction."""


def write_voids(file: Path, result: Result) -> None:
    """Write the voids, one row per void."""


def write_overdense(file: Path, result: Result) -> None:
    """Write the overdense basins, one row per basin."""


def write_divergence(file: Path, result: Result) -> None:
    """Write the divergence grid, one row per cell (``keep_divergence``); read
    back bit for bit by ``read_divergence``."""


def write_galaxies(file: Path, result: Result) -> None:
    """Write where each tracer ended up, one row per tracer
    (``keep_assignment``)."""


class DisplacementMean:
    """otswap's mean displacement file, shaped for ``find_voids_box`` as one
    realization::

        m = bitvf.io.read_displacement_mean("displacement_mean.fits")
        config.n_realizations_mean = m.n_realizations_averaged
        r = bitvf.find_voids_box(m.tracers, m.displacement, m.valid, mps=m.mps, box=m.limits,
                                 ids=m.ids, config=config)

    Read-only."""

    @property
    def tracers(self) -> NDArray[np.float64]:
        """Shape (N, 3)."""
    @property
    def displacement(self) -> NDArray[np.float64]:
        """Shape (1, N, 3); NaN where no realization was valid."""
    @property
    def valid(self) -> NDArray[np.bool_]:
        """Shape (1, N)."""
    @property
    def mps(self) -> float: ...
    @property
    def limits(self) -> Optional[otswap.BoxLimits]:
        """The box limits, when the file records them as limits."""
    @property
    def ids(self) -> NDArray[np.int64]:
        """A new copy on each access; 0..N-1 when the file has none."""
    @property
    def n_realizations_averaged(self) -> int: ...


class Field:
    """A field at the centres of cells, for ``find_voids_box_from_field``.
    Read-only."""

    @property
    def centres(self) -> NDArray[np.float64]:
        """Shape (n, 3)."""
    @property
    def values(self) -> NDArray[np.float64]:
        """Shape (n,)."""
    @property
    def weights(self) -> Optional[NDArray[np.float64]]:
        """Shape (n,); None when no weight column was named."""


def read_divergence(file: Path, *, angle_unit: str) -> DivergenceGrid:
    """Read a divergence file written by ``write_divergence``; ``angle_unit``,
    ``"deg"`` or ``"rad"``, is that of a lightcone grid's angles, unused by a
    box grid."""


def read_displacement_mean(file: Path) -> DisplacementMean:
    """Read otswap's mean displacement file of a box reconstruction."""


def read_field(file: Path, *, centres: Sequence[str], value: str, weight: Optional[str] = None,
               delimiter: str = " ", comment: str = "#") -> Field:
    """Read a field at cell centres from a table: ``centres`` the names of the
    x, y and z columns; in ASCII, names from the last ``###`` line before the
    data, or 0-based indices given as digits. Every value read must be
    finite. A divergence file is valid input."""


__all__: list[str]
