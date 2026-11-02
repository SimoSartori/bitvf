"""The Back-in-time Void Finder: voids from the optimal transport
displacement field of a tracer catalogue.

Built on otswap, whose Python types bitvf takes and returns as they are:
``otswap.Config``, ``CorrectionConfig``, ``BoxLimits``, ``DistanceTable``,
``MpsProfile``, ``Result`` and ``RealSpaceCatalog``. Coordinates are NumPy
arrays of shape (N, 3), in coordinate units, the unit of the tracers given.
Inputs are converted to float64 (``ids`` to int64); the library works on its
own copy of them. The number of threads follows ``OMP_NUM_THREADS``.

Results are read-only: their arrays are read-only views that keep the
result alive, and none can be made from Python. Every per-cell array of a
box grid has shape (nx, ny, nz), ``DivergenceGrid.n_cells_axis``, C order,
so that ``a[i, j, k]`` is the cell (i, j, k); a cell row stored as a value
(``TracerAssignment.cell_index``, ``Partition.void_seed``) is a flat row,
and ``numpy.unravel_index(row, n_cells_axis)`` gives the cell. Flags are
bool; ``BasinCatalog.edge_flag``, 0 to 2, is uint8. An empty array has shape
(0,), or (0, 3) for triplets.

Every search releases the GIL while it runs and prints its report to
``sys.stdout`` as it runs, each line as its step finishes. Ctrl+C stops a
search at its next report line with ``KeyboardInterrupt``; a silent search,
or a long step, runs on until its next line or its end, and a search called
from a thread other than the main one is not interrupted. An exception
raised by ``sys.stdout`` stops the printing, and the search raises it when
it ends.

Files are written and read by the submodule ``bitvf.io``.
"""

import enum
from typing import Literal, Optional, Union

import numpy as np
import otswap
from numpy.typing import ArrayLike, NDArray
from typing_extensions import TypeAlias

from . import io as io

__version__: str

Verbosity: TypeAlias = Literal["silent", "normal", "detailed"]
"""What a search prints, to ``sys.stdout``, each line starting with
``bitvf:``:

- ``"silent"``: nothing;
- ``"normal"``: what the box limits and a reused grid left out, then one line
  per call with its time::

      bitvf: findVoidsBox: <V> voids and <O> overdense basins from <N> tracers in <T> s; <W> without a valid displacement

- ``"detailed"``: as ``"normal"``, with one line per step of the search: the
  mean particle separation and its source, the grid, the smoothing, the
  basins.

Every count printed can be read from the object returned.
"""

Box: TypeAlias = Union[otswap.BoxLimits, tuple[ArrayLike, ArrayLike]]


class Error(Exception):
    """Raised by every bitvf function on invalid input or failure; its message
    starts with ``otswap:`` when otswap raised it. Unrelated to
    ``otswap.Error``."""


# --------------------------------------------------------------------------
# configuration

class LengthUnit(enum.Enum):
    """Unit of a length the user chooses."""

    MPS = 0
    """Multiples of the mean particle separation."""
    COORDINATE = 1
    """The unit of the tracer coordinates."""


class Smoothing:
    """Gaussian smoothing of the divergence field before the watershed.

    Each cell takes the mean of the divergence of the cells around it,
    weighted by their tracer weight times exp(-r^2 / (2 sigma^2)); in a box
    the cells within floor(3.5 sigma / s) cells along each axis, s the cell
    side. A larger sigma merges small basins; sigma = 0 disables it.
    """

    sigma: float
    """Width of the Gaussian, in ``unit``; finite and non-negative. 1 mps is a
    reasonable starting value."""
    unit: LengthUnit

    def __init__(self, *, sigma: float = 1.0, unit: LengthUnit = LengthUnit.MPS) -> None: ...
    def __eq__(self, other: object) -> bool: ...


class Config:
    """Free parameters of the void finder. Each field is converted to its
    type when set, and checked by the call that reads it."""

    smoothing: Smoothing
    """A reference into this Config: ``config.smoothing.sigma = 2`` changes
    it. Assigning a Smoothing copies it in."""
    divergence_cell_size: Optional[float]
    """Side of the divergence-grid cell, in units of the mean particle
    separation; None for min(1, (20 / n_realizations)^(1/3)). Changes the
    result."""
    keep_divergence: bool
    """Fill ``Result.divergence`` with the grid and the field."""
    keep_unsmoothed_divergence: bool
    """Also fill ``DivergenceGrid.divergence_unsmoothed``; needs
    ``keep_divergence``."""
    keep_assignment: bool
    """Fill ``Result.assignment``, where each tracer ended up."""
    n_realizations_mean: int
    """Realizations averaged into each displacement given: 1 for single
    realizations, the NREC of a mean displacement file
    (``io.DisplacementMean.n_realizations_averaged``). A record, written as
    NRECMEAN."""
    verbosity: Verbosity

    def __init__(self, *, smoothing: Smoothing = Smoothing(), divergence_cell_size: Optional[float] = None,
                 keep_divergence: bool = False, keep_unsmoothed_divergence: bool = False,
                 keep_assignment: bool = False, n_realizations_mean: int = 1,
                 verbosity: Verbosity = "normal") -> None: ...
    def __eq__(self, other: object) -> bool: ...


class BoxRsd:
    """Correction of redshift-space distortions in a box, as
    ``otswap.real_space_box`` applies it."""

    axis: Optional[int]
    """The line of sight: 0, 1 or 2, for x, y or z; None is refused by
    ``run_box``."""
    redshift: float
    """The redshift of the box, at which f is read from the distance table."""
    bias: float
    """The linear bias of the tracers; finite and positive."""
    sigma: float
    """Width of the neighbour average, in coordinate units; 0 for none. For
    coordinates in Mpc/h, 10 is a reasonable starting value, to check in each
    analysis."""
    correction: otswap.CorrectionConfig
    """A reference into this BoxRsd; assigning one copies it in."""
    zspace_voids: bool
    """Also search voids on the redshift-space reconstruction."""

    def __init__(self, *, axis: Optional[int], redshift: float, bias: float, sigma: float,
                 correction: otswap.CorrectionConfig = otswap.CorrectionConfig(),
                 zspace_voids: bool = False) -> None: ...
    def __eq__(self, other: object) -> bool: ...


class BoxConfig:
    """The parameters of ``run_box``. The nested configurations ``ot``,
    ``rsd`` and ``voids`` read as references into this BoxConfig;
    assigning one copies it in."""

    ot: otswap.Config
    """The reconstruction, for both reconstructions of a run with
    correction, with the same seed; 0 draws one."""
    limits: Optional[otswap.BoxLimits]
    """The box of the run; None for the bounding box of the tracers, and of
    the randoms when given. A pair of corners is accepted."""
    mps: Optional[float]
    """The mean particle separation of both reconstructions; None for
    otswap's."""
    rsd: Optional[BoxRsd]
    """The correction of redshift-space distortions; None for none."""
    voids: Config
    """The void search, on every pass; its verbosity also governs the lines
    of the pipeline."""
    keep_reconstructions: bool
    """Keep the otswap results and the corrected catalogue in
    ``BitVFResult``."""

    def __init__(self, *, ot: otswap.Config = otswap.Config(), limits: Optional[Box] = None,
                 mps: Optional[float] = None, rsd: Optional[BoxRsd] = None, voids: Config = Config(),
                 keep_reconstructions: bool = False) -> None: ...
    def __eq__(self, other: object) -> bool: ...


# --------------------------------------------------------------------------
# results

class BasinCatalog:
    """One catalogue of watershed basins, the voids or the overdense basins.
    Row b is basin b, in the order of its seed cell. Read-only."""

    @property
    def n_basins(self) -> int: ...
    @property
    def centre(self) -> NDArray[np.float64]:
        """Cartesian centre, shape (n_basins, 3), in coordinate units."""
    @property
    def centre_sky(self) -> NDArray[np.float64]:
        """Sky coordinates of the centre; shape (0, 3) in a box."""
    @property
    def radius(self) -> NDArray[np.float64]:
        """Radius of the sphere with the volume of the core."""
    @property
    def radius_ext(self) -> NDArray[np.float64]:
        """Radius of the sphere with the volume of the extended basin."""
    @property
    def delta_av(self) -> NDArray[np.float64]:
        """Mean density contrast of the tracers in the core; (0,) after a field
        search."""
    @property
    def delta_2mps(self) -> NDArray[np.float64]:
        """Density contrast within 2 mps of the centre; (0,) after a field
        search."""
    @property
    def edge_flag(self) -> NDArray[np.uint8]:
        """Contact with the boundary: 0 none, 1 only the extended basin, 2 the
        core."""
    @property
    def contamination(self) -> NDArray[np.float64]:
        """Fraction of the basin's surface on the boundary."""
    @property
    def n_tracers(self) -> NDArray[np.uint64]:
        """Tracers starting in the cells of the core; (0,) after a field
        search."""
    @property
    def field_centre(self) -> NDArray[np.float64]:
        """The field at the centre, after a field search; (0,) otherwise."""


class DivergenceGrid:
    """The grid of the void finder and the divergence field on it. Read-only:
    it comes from ``Result.divergence``, ``box_grid`` or
    ``io.read_divergence``, and can be passed back as ``grid``, bit for bit.
    Per-cell arrays have shape ``n_cells_axis`` (``centre``: plus 3), C
    order; (0,) when not filled."""

    @property
    def geometry(self) -> Literal["box", "lightcone"]: ...
    @property
    def angle_unit(self) -> Optional[Literal["deg", "rad"]]:
        """None for a box."""
    @property
    def n_cells(self) -> int: ...
    @property
    def cell_size(self) -> float:
        """Cell side, in coordinate units."""
    @property
    def lower(self) -> tuple[float, float, float]:
        """The lower grid face."""
    @property
    def upper(self) -> tuple[float, float, float]:
        """The upper face of the padded box as placed."""
    @property
    def n_cells_axis(self) -> tuple[int, int, int]: ...
    @property
    def limits(self) -> Optional[otswap.BoxLimits]:
        """The box limits of the search that filled the grid; a record."""
    @property
    def n_realizations(self) -> int: ...
    @property
    def n_realizations_mean(self) -> int: ...
    @property
    def distance_min(self) -> float: ...
    @property
    def distance_max(self) -> float: ...
    @property
    def ra_min(self) -> float: ...
    @property
    def ra_max(self) -> float: ...
    @property
    def dec_min(self) -> float: ...
    @property
    def dec_max(self) -> float: ...
    @property
    def cell_size_mps(self) -> float:
        """The cell size in units of the mps."""
    @property
    def mps(self) -> float:
        """Box: the mean particle separation the grid was built with."""
    @property
    def mps_profile(self) -> Optional[otswap.MpsProfile]:
        """None for a box."""
    @property
    def cell_ra_min(self) -> NDArray[np.float64]: ...
    @property
    def cell_ra_max(self) -> NDArray[np.float64]: ...
    @property
    def cell_dec_min(self) -> NDArray[np.float64]: ...
    @property
    def cell_dec_max(self) -> NDArray[np.float64]: ...
    @property
    def cell_distance_min(self) -> NDArray[np.float64]: ...
    @property
    def cell_distance_max(self) -> NDArray[np.float64]: ...
    @property
    def centre(self) -> NDArray[np.float64]:
        """Centre of each cell, shape (nx, ny, nz, 3)."""
    @property
    def centre_sky(self) -> NDArray[np.float64]: ...
    @property
    def divergence(self) -> NDArray[np.float64]:
        """The divergence the watershed ran on, after smoothing."""
    @property
    def divergence_unsmoothed(self) -> NDArray[np.float64]:
        """The divergence before smoothing, with
        ``keep_unsmoothed_divergence``."""
    @property
    def volume(self) -> NDArray[np.float64]: ...
    @property
    def weight(self) -> NDArray[np.float64]:
        """Tracer weight of each cell, the weight of the smoothing."""
    @property
    def n_displacements(self) -> NDArray[np.uint64]: ...
    @property
    def n_tracers(self) -> NDArray[np.uint64]: ...
    @property
    def void_index(self) -> NDArray[np.int64]:
        """Row of the cell's void; -1 none."""
    @property
    def overdense_index(self) -> NDArray[np.int64]:
        """Row of the cell's overdense basin; -1 none."""
    @property
    def masked(self) -> NDArray[np.bool_]:
        """(0,) in a box."""


class TracerAssignment:
    """Where each input tracer ended up, one row per input row. Read-only."""

    @property
    def cell_index(self) -> NDArray[np.int64]:
        """Flat row of the cell the tracer starts in; -1 when outside the grid
        or left out."""
    @property
    def void_index(self) -> NDArray[np.int64]: ...
    @property
    def overdense_index(self) -> NDArray[np.int64]: ...
    @property
    def distance(self) -> NDArray[np.float64]:
        """Distance from the centre of the tracer's basin; NaN when in none."""
    @property
    def in_mask(self) -> NDArray[np.bool_]: ...
    @property
    def in_box(self) -> NDArray[np.bool_]:
        """False for a tracer outside the box limits, skipped."""
    @property
    def used(self) -> NDArray[np.bool_]:
        """True when a valid displacement of the tracer entered the divergence;
        (0,) after ``find_voids_box_from_divergence``."""


class Counts:
    """What the void finder counted. Read-only."""

    @property
    def tracers(self) -> int: ...
    @property
    def tracers_used(self) -> int: ...
    @property
    def tracers_uncorrected(self) -> int: ...
    @property
    def tracers_outside_box(self) -> int: ...
    @property
    def tracers_outside_grid(self) -> int: ...
    @property
    def tracers_outside_cut(self) -> int: ...
    @property
    def randoms(self) -> int: ...
    @property
    def randoms_outside_grid(self) -> int: ...
    @property
    def displacements(self) -> int: ...
    @property
    def displacements_rejected(self) -> int: ...


class Result:
    """Output of one void search. Read-only. Every field but
    ``elapsed_seconds`` is the same bit for bit whatever the number of
    threads."""

    @property
    def geometry(self) -> Literal["box", "lightcone"]: ...
    @property
    def angle_unit(self) -> Optional[Literal["deg", "rad"]]: ...
    @property
    def config(self) -> Config:
        """The configuration the search ran with, ``divergence_cell_size`` the
        value used: a new copy on each access."""
    @property
    def n_tracers(self) -> int:
        """0 for a field result."""
    @property
    def n_realizations(self) -> int: ...
    @property
    def ids(self) -> NDArray[np.int64]:
        """A new read-only copy on each access; 0..n_tracers-1 when none were
        set."""
    @property
    def random_seed(self) -> Optional[int]:
        """The seed of the randoms ``run_box`` drew; None when it drew none."""
    @property
    def voids(self) -> BasinCatalog: ...
    @property
    def overdense(self) -> BasinCatalog: ...
    @property
    def divergence(self) -> Optional[DivergenceGrid]:
        """With ``keep_divergence``."""
    @property
    def assignment(self) -> Optional[TracerAssignment]:
        """With ``keep_assignment``."""
    @property
    def uncorrected(self) -> NDArray[np.bool_]:
        """Real-space pass of ``run_box`` with a correction: True for a tracer
        the correction left unmoved; (0,) elsewhere, which reads as all
        False."""
    @property
    def outside_grid(self) -> NDArray[np.bool_]:
        """True for a tracer outside the reused grid; (0,) when no grid was
        reused."""
    @property
    def outside_cut(self) -> NDArray[np.bool_]:
        """(0,) in a box."""
    @property
    def limits(self) -> Optional[otswap.BoxLimits]: ...
    @property
    def mps(self) -> float:
        """The mean particle separation the search used; NaN for a field result
        given none."""
    @property
    def mps_profile(self) -> Optional[otswap.MpsProfile]:
        """None for a box."""
    @property
    def distances(self) -> Optional[otswap.DistanceTable]:
        """None for a box."""
    @property
    def counts(self) -> Counts: ...
    @property
    def elapsed_seconds(self) -> float: ...


class Partition:
    """The basins of a field, as ``watershed_box`` finds them. Read-only."""

    @property
    def n_voids(self) -> int: ...
    @property
    def n_overdense(self) -> int: ...
    @property
    def void_index(self) -> NDArray[np.int64]:
        """Row of each cell's void, -1 none; the shape of the field."""
    @property
    def overdense_index(self) -> NDArray[np.int64]: ...
    @property
    def void_seed(self) -> NDArray[np.int64]:
        """Flat row of each void's seed cell, ascending."""
    @property
    def overdense_seed(self) -> NDArray[np.int64]: ...
    @property
    def void_cells(self) -> NDArray[np.uint64]: ...
    @property
    def void_core_cells(self) -> NDArray[np.uint64]: ...
    @property
    def overdense_cells(self) -> NDArray[np.uint64]: ...
    @property
    def overdense_core_cells(self) -> NDArray[np.uint64]: ...


class BitVFResult:
    """Output of ``run_box``. Read-only."""

    @property
    def voids(self) -> Result:
        """The voids of the run, one row per input tracer."""
    @property
    def zspace_voids(self) -> Optional[Result]:
        """With ``BoxRsd.zspace_voids``."""
    @property
    def reconstruction(self) -> Optional[otswap.Result]:
        """With ``keep_reconstructions``: the reconstruction the voids were
        found on."""
    @property
    def zspace_reconstruction(self) -> Optional[otswap.Result]: ...
    @property
    def real_space(self) -> Optional[otswap.RealSpaceCatalog]: ...


# --------------------------------------------------------------------------
# functions

def find_voids_box(tracers: Union[ArrayLike, otswap.Result], displacement: Optional[ArrayLike] = None,
                   valid: Optional[ArrayLike] = None, *, mps: Optional[float] = None, box: Optional[Box] = None,
                   config: Config = Config(), ids: Optional[ArrayLike] = None,
                   grid: Optional[DivergenceGrid] = None) -> Result:
    """Find voids and overdense basins in a box.

    Two forms. With arrays: ``tracers`` (N, 3); ``displacement`` (R, N, 3),
    required; ``valid`` (R, N), bool or 0/1 numbers, or None for all valid;
    ``mps``, the mean particle separation, required; ``box``, the box limits,
    whose outside tracers are skipped and on which the grid is built, or None
    for the extent of the tracers and their displacements' endpoints;
    ``ids``, one 64-bit integer per tracer. With an ``otswap.Result`` of
    ``otswap.reconstruct_box`` as the first argument, every one of those is
    read from it, and none may be given.

    ``grid`` reuses the cells of a grid (``Result.divergence``, ``box_grid``,
    ``io.read_divergence``); tracers outside it are excluded
    (``Result.outside_grid``).
    """


def find_voids_box_from_divergence(grid: DivergenceGrid, tracers: ArrayLike, *, config: Config = Config(),
                                   ids: Optional[ArrayLike] = None) -> Result:
    """Find voids and overdense basins from a divergence computed earlier,
    with the tracers given: the grid's ``divergence_unsmoothed`` smoothed with
    ``config.smoothing``."""


def find_voids_box_from_field(centres: ArrayLike, values: ArrayLike, *, config: Config = Config(),
                              weights: Optional[ArrayLike] = None, mps: Optional[float] = None) -> Result:
    """Find voids and overdense basins in a field given at the centres (n, 3),
    in any order, of the cells of a box lattice; voids are below 0. Without
    ``mps``, a positive sigma must be in ``LengthUnit.COORDINATE``."""


def box_grid(reconstruction: Optional[otswap.Result] = None, *, box: Optional[Box] = None,
             mps: Optional[float] = None, cell_size: Optional[float] = None) -> DivergenceGrid:
    """The grid ``find_voids_box`` builds, bit for bit: on ``box`` with
    ``mps`` and ``cell_size`` (all three required), or from a
    reconstruction, ``cell_size`` then defaulting to its realizations'.
    Only the fields that define the cells are set."""


def watershed_box(field: ArrayLike, *, verbosity: Verbosity = "normal") -> Partition:
    """Split a field of shape (nx, ny, nz), in any memory order, into basins
    by steepest descent; voids below 0. No smoothing, no catalogue."""


def run_box(tracers: ArrayLike, randoms: Optional[ArrayLike] = None, *, config: BoxConfig = BoxConfig(),
            distances: Optional[otswap.DistanceTable] = None, ids: Optional[ArrayLike] = None) -> BitVFResult:
    """Reconstruct a box and find its voids. ``randoms`` (M, 3), or None to
    have bitvf draw them in the box (``Result.random_seed``). With
    ``config.rsd``, the redshift-space distortions are corrected and the
    tracers reconstructed again in real space; ``distances``, carrying the
    growth rate, is then required. Prints otswap's report lines too."""


__all__: list[str]
