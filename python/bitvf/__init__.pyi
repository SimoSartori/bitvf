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
bool; ``BasinCatalog.edge_flag``, 0 to 2, is uint8. A value that was not
computed is NaN; a basin or cell index that does not apply is -1. An empty
array has shape (0,), or (0, 3) for triplets.

Every search releases the GIL while it runs and prints its report to
``sys.stdout`` as it runs, each line as its step finishes. Ctrl+C stops a
search at its next report line with ``KeyboardInterrupt``; a silent search,
or a long step, runs on until its next line or its end, and a search called
from a thread other than the main one is not interrupted. An exception
raised by ``sys.stdout`` stops the printing, and the search raises it when
it ends.

Lightcone searches take sky positions as (right ascension, declination,
redshift) rows in the ``angle_unit`` of the call, ``"deg"`` or ``"rad"``, and
report the sky arrays of what they return in that unit: a ``Result``,
``BitVFResult`` or ``DivergenceGrid`` made by a call keeps its unit, and its
parts (``Result.voids``, ``Result.divergence``, ``BitVFResult.voids``) the
unit of the object that holds them. A grid passed back to a call in another
unit is accepted. The per-cell arrays of a lightcone grid are flat,
``(n_cells,)`` and ``(n_cells, 3)``, in the order of the cells.

Files are written and read by the submodule ``bitvf.io``.
"""

import enum
from typing import Literal, Optional, Union, overload

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

      bitvf: kept <k> of <N> tracers: <b> outside the box [<x0>, <x1>] x [<y0>, <y1>] x [<z0>, <z1>]
      bitvf: kept <k'> of <k> tracers: <g> outside the reused grid
      bitvf: moved <m> tracer(s) onto the box edge
      bitvf: findVoidsBox: <V> voids and <O> overdense basins from <N> tracers in <T> s; <W> without a valid displacement

  The box line is printed when the search has box limits (given, the
  reconstruction's, or the reused grid's); the grid line whenever a grid is
  reused, even with no tracer outside it, and by
  ``find_voids_box_from_divergence`` as "outside the divergence grid". The
  moved line is printed by the search of a reconstruction without box
  limits, when <m> is not 0: a tracer outside the measured box by rounding
  is moved onto its edge. <W>
  counts the tracers read without a valid displacement. The call lines of
  the other searches are
  ``bitvf: findVoidsBoxFromDivergence: <V> voids and <O> overdense basins from <N> tracers in <T> s``
  and
  ``bitvf: findVoidsBoxFromField: <V> voids and <O> overdense basins on <nx> x <ny> x <nz> cells in <T> s``.
  In a lightcone the lines are::

      bitvf: kept <k> of <N> tracers: <m> on unobserved pixels, <c> outside the redshift cut
      bitvf: kept <k> of <M> randoms: <m> on unobserved pixels, <c> outside the redshift cut
      bitvf: grid of <C> cells in <S> shells, cell size <s> mps (<source>), distances [<d0>, <d1>], RA [<a0>, <a1>] deg, Dec [<b0>, <b1>] deg (<from>) in <t> s
      bitvf: kept <k'> of <k> tracers: <g> outside the reused grid
      bitvf: <n> tracers start in a masked cell, in no basin
      bitvf: divergence of the valid displacements of <u> tracers, <r> displacements rejected (masked cells, outside or below the grid) in <t> s
      bitvf: filled <f> cells in <p> passes, <e> left empty in <t> s
      bitvf: findVoidsLightcone: <V> voids and <O> overdense basins from <N> tracers in <T> s; <W> without a valid displacement

  The redshift cut appears in the selection lines only when it has a finite
  bound; the grid line adds ", central sphere included" when the grid has
  one, its source is default, given or reused, and it ends "(from the
  tracers and randoms)" or "(from the reused grid)"; the reused-grid line is
  printed whenever a grid is reused; the masked-cell, divergence and fill
  lines only when <n>, <r> or <e> is not 0
  (``Counts.tracers_in_masked_cells``, ``displacements_rejected``). The
  call line names ``findVoidsLightconeFromCartesian`` for that function; <W>
  counts the tracers kept and inside the grid without a valid displacement.
- ``"detailed"``: as ``"normal"``, with one line per step of the search in
  the order it runs, each with its time: the mean particle separation and
  its source, the grid, its cell side and the box it covers, the smoothing
  sigma, the basins found by the watershed.

Each line is printed as soon as its step finishes. A call that raises
leaves the lines of the steps it finished, and never prints its call line.
Every count printed can be read from the object returned
(``Result.counts``, ``BasinCatalog.n_basins``); step times and descriptive
values (mps, cell side, box, sigma) appear only in the lines. Times have two
decimals.
"""

Box: TypeAlias = Union[otswap.BoxLimits, tuple[ArrayLike, ArrayLike]]
"""Box limits: an ``otswap.BoxLimits``, or a pair of corners
((x0, y0, z0), (x1, y1, z1)), in coordinate units. Each bound must be
finite, with each lower bound below the upper one. Closed: a point on a
face is inside."""


class Error(Exception):
    """Raised by every bitvf function on invalid input or failure.

    An error raised by otswap inside a bitvf call reaches the caller as an
    ``Error`` whose message is ``"otswap: "`` followed by otswap's own.
    ``Error`` is unrelated to ``otswap.Error``: an ``except otswap.Error``
    does not catch it, and catching ``bitvf.Error`` catches every error of
    a bitvf call but a ``TypeError`` for a wrong number of positional
    arguments or an unknown keyword, and ``KeyboardInterrupt``.
    """


# --------------------------------------------------------------------------
# configuration

class LengthUnit(enum.Enum):
    """Unit of a length the user chooses."""

    MPS = 0
    """Multiples of the mean particle separation: the box mps, or in a
    lightcone mps(z) at each cell's redshift."""
    COORDINATE = 1
    """The same length unit as the tracer coordinates."""


class Smoothing:
    """Gaussian smoothing of the divergence field before the watershed.

    Each cell takes the weighted mean of the divergence of the cells around
    it, itself included. The weight of a cell is its tracer weight
    (``DivergenceGrid.weight``) times exp(-r^2 / (2 sigma^2)), r the distance
    between the two cell centres. A cell whose weights sum to at most 1e-12
    times the largest tracer weight of the grid takes 0; this includes a
    filled cell (weight 0) whose surrounding cells hold no weight. A larger
    sigma merges small basins into larger ones; sigma = 0 disables the
    smoothing. With sigma > 0 the smoothed field and the
    smoothed weight are kept beside the divergence and the weight
    (``DivergenceGrid.divergence_smoothed``, ``weight_smoothed``).

    Box: the cells around are those within floor(3.5 sigma / s) cells
    along each axis, s the cell side, cut at the faces of the grid: a cube,
    not a sphere. When 3.5 sigma / s lies within 1e-9 relative of an integer,
    that integer is taken instead, so that a cell side rounded by one ulp
    gives the same cube. Below sigma = s / 3.5 the cube holds the cell alone.
    In ``LengthUnit.MPS``, sigma is multiplied by the box mps.

    Lightcone: the cells around are those whose centre lies closer than
    4 sigma, with the same weights. Masked cells (``DivergenceGrid.masked``)
    enter no mean and take 0. In ``LengthUnit.MPS``, sigma is multiplied by
    mps(z) at the redshift of the centre of the shell of the cell being
    smoothed.

    Each field is a keyword of the constructor and a read/write attribute,
    converted to its type when set; its rules are checked by the search
    that reads it. ``repr`` gives the constructor call; two smoothings
    compare equal when both fields are equal.
    """

    sigma: float
    """Width of the Gaussian, in ``unit``; default 1. Must be finite and
    non-negative; 1 mps is a reasonable starting value. A positive sigma
    must be between about 1e-162 and 1e154 in coordinate units; use 0 to
    disable the smoothing."""
    unit: LengthUnit
    """Unit of ``sigma``; default ``LengthUnit.MPS``."""

    def __init__(self, *, sigma: float = 1.0, unit: LengthUnit = LengthUnit.MPS) -> None: ...
    def __eq__(self, other: object) -> bool: ...


class Config:
    """Free parameters of the void finder: the ``config`` of every search.

    Nothing here depends on the geometry. Voids and overdense basins are
    always both found. Every field is a keyword of the constructor, with the
    default of the C++ struct, and a read/write attribute. A value is
    converted to the field's type when it is stored, and raises ``Error`` if
    it cannot be; its rules are checked by the call that reads it, which
    refuses it with the library's message. Each call works on its own copy,
    so the object passed is never changed. ``repr`` gives the constructor
    call; two configs compare equal when every field is equal.
    """

    smoothing: Smoothing
    """Smoothing of the divergence field; default ``Smoothing()``, sigma
    1 mps. A reference into this Config: ``config.smoothing.sigma = 2``
    changes it. Assigning a Smoothing copies it in."""
    divergence_cell_size: Optional[float]
    """Side of the divergence-grid cell, in units of the mean particle
    separation; None, the default, for min(1, (20 / n_realizations)^(1/3)):
    1 up to 20 realizations, (20 / n_realizations)^(1/3) above. A search on a
    mean displacement file runs with one realization, so its default is 1.
    Must be positive when set, and None when a grid is reused.
    ``Result.config`` holds the value used. Unrelated to
    ``otswap.Config.cell_size``, which affects speed only: this one changes
    the result."""
    keep_divergence: bool
    """Fill ``Result.divergence`` with the grid and the divergence field;
    default False, since it holds one entry per cell in several arrays."""
    keep_assignment: bool
    """Fill ``Result.assignment``, where each tracer ended up; default False,
    since it holds one entry per tracer in several arrays."""
    n_realizations_mean: int
    """Realizations averaged into each displacement given; default 1, for
    the displacements of single realizations. For the ``otswap.Result`` read
    from a mean displacement file (``otswap.io.read_displacement_mean``), its
    ``n_realizations_mean``. A record
    only, written as the keyword NRECMEAN of bitvf's files; no computation
    reads it. Must be at least 1, at most 4294967295."""
    verbosity: Verbosity
    """What the search prints; default ``"normal"``. See ``Verbosity``."""

    def __init__(self, *, smoothing: Smoothing = Smoothing(), divergence_cell_size: Optional[float] = None,
                 keep_divergence: bool = False, keep_assignment: bool = False,
                 n_realizations_mean: int = 1, verbosity: Verbosity = "normal") -> None: ...
    def __eq__(self, other: object) -> bool: ...


class BoxRsd:
    """Correction of redshift-space distortions in a box, as
    ``otswap.real_space_box`` applies it: the ``rsd`` of ``BoxConfig``.

    ``axis``, ``redshift``, ``bias`` and ``sigma`` have no default: each is a
    required keyword of the constructor. Fields are read/write attributes,
    converted when set and checked by ``run_box``. ``repr`` gives the
    constructor call; two compare equal when every field is equal.
    """

    axis: Optional[int]
    """The line of sight: 0, 1 or 2, for x, y or z. None is accepted here
    and refused by ``run_box``."""
    redshift: float
    """The redshift of the box, at which the growth rate f is read from the
    distance table; inside the table's range."""
    bias: float
    """The linear bias of the tracers, constant; finite and positive."""
    sigma: float
    """Width of the neighbour average of the displacements, in coordinate
    units; 0 for none, finite and non-negative. For coordinates in Mpc/h,
    10 is a reasonable starting value; the best value depends on the sample,
    and should be checked in each analysis."""
    correction: otswap.CorrectionConfig
    """The weighting of the neighbour average, and what the correction
    prints, passed to ``otswap.real_space_box`` unchanged; default
    ``otswap.CorrectionConfig()``. A reference into this BoxRsd; assigning
    one copies it in."""
    zspace_voids: bool
    """Also search voids on the first, redshift-space, reconstruction
    (``BitVFResult.zspace_voids``); default False."""

    def __init__(self, *, axis: Optional[int], redshift: float, bias: float, sigma: float,
                 correction: otswap.CorrectionConfig = otswap.CorrectionConfig(),
                 zspace_voids: bool = False) -> None: ...
    def __eq__(self, other: object) -> bool: ...


class BoxConfig:
    """The parameters of ``run_box``.

    Every field is a keyword of the constructor and a read/write attribute.
    The nested configurations ``ot``, ``rsd`` and ``voids`` read as
    references into this BoxConfig, so ``cfg.voids.keep_divergence = True``
    changes ``cfg``; assigning one copies it in. ``repr`` gives the
    constructor call; two compare equal when every field is equal.
    """

    ot: otswap.Config
    """The reconstruction, used for both reconstructions of a run with
    correction, with the same seed; default ``otswap.Config()``. A seed of 0
    is replaced by one drawn at random before the first reconstruction, and
    both use it; the ``config`` of the reconstruction records it."""
    limits: Optional[otswap.BoxLimits]
    """The box of the run, given to ``otswap.reconstruct_box`` as limits on
    both reconstructions; a pair of corners is accepted (see ``Box``). None,
    the default, for the bounding box of the tracers, or of the tracers and
    the randoms when randoms are given."""
    mps: Optional[float]
    """The mean particle separation, in coordinate units, given to both
    reconstructions. None, the default: otswap computes it on each from the
    volume of the box and the tracers inside it."""
    rsd: Optional[BoxRsd]
    """The correction of redshift-space distortions; None, the default, for
    one reconstruction and one search, no correction. A reference into this
    BoxConfig; assigning one copies it in."""
    voids: Config
    """The void search, on every pass; default ``Config()``. Its verbosity
    also governs the lines of the pipeline itself."""
    keep_reconstructions: bool
    """Keep the otswap results and the corrected catalogue in
    ``BitVFResult``; default False."""

    def __init__(self, *, ot: otswap.Config = otswap.Config(), limits: Optional[Box] = None,
                 mps: Optional[float] = None, rsd: Optional[BoxRsd] = None, voids: Config = Config(),
                 keep_reconstructions: bool = False) -> None: ...
    def __eq__(self, other: object) -> bool: ...


class SecondPassRandoms(enum.Enum):
    """The randoms of the second, real-space, reconstruction of
    ``run_lightcone``."""

    NONE = ...
    """The first reconstruction's randoms, unchanged."""
    RESAMPLE = ...
    """The first reconstruction's randoms, with redshifts drawn from the
    corrected tracers."""
    EXTERNAL = ...
    """Randoms given to ``run_lightcone`` for the second reconstruction."""


class LightconeRsd:
    """Correction of redshift-space distortions in a lightcone, as
    ``otswap.real_space_lightcone`` applies it: the ``rsd`` of
    ``LightconeConfig``.

    ``bias`` and ``sigma`` have no default: each is a required keyword of
    the constructor. Fields are read/write attributes, converted when set
    and checked by ``run_lightcone``. ``repr`` gives the constructor call;
    two compare equal when every field is equal.
    """

    bias: otswap.BiasTable
    """The linear bias of the tracers as a function of redshift; positive at
    the redshift of every tracer kept."""
    sigma: float
    """Width of the neighbour average of the displacements, in coordinate
    units; 0 for none, finite and non-negative. For coordinates in Mpc/h,
    10 is a reasonable starting value; the best value depends on the sample,
    and should be checked in each analysis."""
    correction: otswap.CorrectionConfig
    """The weighting of the neighbour average, and what the correction
    prints, passed to ``otswap.real_space_lightcone`` unchanged; default
    ``otswap.CorrectionConfig()``. A reference into this LightconeRsd;
    assigning one copies it in."""
    second_pass_randoms: SecondPassRandoms
    """The randoms of the second reconstruction. ``RESAMPLE``, the default,
    gives them the redshift distribution of the corrected tracers; ``NONE``
    keeps the observed one; ``EXTERNAL`` takes randoms built for the
    real-space sample and given to ``run_lightcone``."""
    zspace_voids: bool
    """Also search voids on the first, redshift-space, reconstruction
    (``BitVFResult.zspace_voids``); default False."""

    def __init__(self, *, bias: otswap.BiasTable, sigma: float,
                 correction: otswap.CorrectionConfig = otswap.CorrectionConfig(),
                 second_pass_randoms: SecondPassRandoms = SecondPassRandoms.RESAMPLE,
                 zspace_voids: bool = False) -> None: ...
    def __eq__(self, other: object) -> bool: ...


class LightconeConfig:
    """The parameters of ``run_lightcone``.

    ``n_bins`` has no default: it is a required keyword of the constructor.
    Every field is a read/write attribute. The nested configurations ``ot``,
    ``cut``, ``rsd`` and ``voids`` read as references into this
    LightconeConfig, so ``cfg.voids.keep_divergence = True`` changes
    ``cfg``; assigning one copies it in. ``repr`` gives the constructor
    call; two compare equal when every field is equal.
    """

    ot: otswap.Config
    """The reconstruction, used for both reconstructions of a run with
    correction, with the same seed; default ``otswap.Config()``. A seed of 0
    is replaced by one drawn at random before the first reconstruction, and
    both use it; the ``config`` of the reconstruction records it."""
    n_bins: int
    """The redshift bins of the mps(z) profile of each reconstruction; at
    least 1. 50 is a reasonable starting value."""
    cut: otswap.RedshiftCut
    """The closed range of observed redshift of the run, given to both
    reconstructions and to every search: the tracers outside keep their row
    and are skipped, the randoms outside are removed. A ``(min, max)`` pair
    is accepted, and None for no cut, the default."""
    rsd: Optional[LightconeRsd]
    """The correction of redshift-space distortions; None, the default, for
    one reconstruction and one search, no correction. A reference into this
    LightconeConfig; assigning one copies it in."""
    voids: Config
    """The void search, on every pass; default ``Config()``. Its verbosity
    also governs the lines of the pipeline itself."""
    keep_reconstructions: bool
    """Keep the otswap results and the corrected catalogue in
    ``BitVFResult``; default False."""

    def __init__(self, *, ot: otswap.Config = otswap.Config(), n_bins: int,
                 cut: Optional[Union[otswap.RedshiftCut, tuple[float, float]]] = None,
                 rsd: Optional[LightconeRsd] = None, voids: Config = Config(),
                 keep_reconstructions: bool = False) -> None: ...
    def __eq__(self, other: object) -> bool: ...


# --------------------------------------------------------------------------
# results

class BasinCatalog:
    """One catalogue of watershed basins: the voids, or the overdense basins.

    The columns are the same for both kinds. Row b is basin b, numbered from
    0 in the order of the index of its seed cell, a local minimum of the
    divergence (smoothed when a smoothing is set): below 0 for a void, at or
    above 0 for an overdense basin. In a box every cell belongs to exactly
    one basin. The order does not depend on the number of threads.

    The core of a void is the part of the void where the divergence is
    below 0; the extended void is the whole basin. In a box the core and the
    extended basin of an overdense basin coincide.

    Each cell flows to the neighbour where the divergence drops the most
    (the largest decrease, not divided by the distance), among its
    neighbours (box: the 26 cells that share a face, an edge or a corner
    with it; lightcone: the unmasked cells whose centres are closer than
    r (s_i + s_j) / 2, r = (sqrt(3) + 2) / 2, s_i and s_j the sides of the
    two cells, about 1.87 times their mean), when that value is strictly
    below its own; among equal smallest values, the neighbour with the
    lowest cell row. A cell with no strictly smaller neighbour is a seed. A
    basin is the seed and every cell whose flow reaches it. In a box the
    grid is not periodic: a cell on a face has fewer neighbours.

    Read-only; every array has one row per basin.
    """

    @property
    def n_basins(self) -> int:
        """Number of basins."""
    @property
    def centre(self) -> NDArray[np.float64]:
        """Cartesian centre, shape (n_basins, 3), in coordinate units: the
        seed cell's centre x0, shifted by
        ``s * sum_n div_n (x_n - x0) / sum_n |div_n|`` over the seed's
        neighbouring cells n, with div the field the watershed split (the smoothed divergence, or the
        divergence when sigma is 0) and s = -1
        for a void, +1 for an overdense basin (towards lower divergence in a
        void, higher in an overdense basin), and kept inside the seed
        cell."""
    @property
    def centre_sky(self) -> NDArray[np.float64]:
        """Sky coordinates of the centre, lightcone only, shape (n_basins, 3):
        right ascension and declination in the ``angle_unit`` of the call
        that made the result, right ascension in [0, 2 pi) or [0, 360), and
        the redshift; shape (0, 3) in a box. In degrees, a read-only copy
        made once, with otswap's conversion."""
    @property
    def radius(self) -> NDArray[np.float64]:
        """Radius of the sphere with the volume of the core, in coordinate
        units, shape (n_basins,)."""
    @property
    def radius_ext(self) -> NDArray[np.float64]:
        """Radius of the sphere with the volume of the extended basin, in
        coordinate units, shape (n_basins,)."""
    @property
    def delta_av(self) -> NDArray[np.float64]:
        """Mean density contrast of the tracers in the core, shape
        (n_basins,): the tracers starting in its cells (``n_tracers``) over
        its volume, divided by the mean density mps^-3 (mps(z)^-3 in a
        lightcone), minus 1. Shape (0,) after a field search."""
    @property
    def delta_2mps(self) -> NDArray[np.float64]:
        """Density contrast within 2 mps of the centre, the sphere closed,
        shape (n_basins,): n_t over its volume, divided by the mean density
        mps^-3, minus 1, with n_t the tracers the search reads (inside the
        box limits and the reused grid) in the sphere. In a lightcone, within
        2 mps(z) at the centre's redshift, n_t / (n_r N_t / N_r) - 1, with
        n_t and n_r the tracers and randoms in the sphere and N_t / N_r the
        ratio of the tracers to the randoms over the whole search
        (``Counts.randoms``); NaN when no random lies in the sphere. At
        least 10 to 20 randoms per tracer is a reasonable start. Shape (0,)
        after a field search."""
    @property
    def edge_flag(self) -> NDArray[np.uint8]:
        """Contact with the boundary, shape (n_basins,): 0 none; 1 only the
        extended basin touches it; 2 the core already touches it. The
        boundary is, in a box, the faces of the grid, with or without box
        limits; in a lightcone, the masked cells and the border of the grid.
        In a box, a basin touches it when one of its cells lies on a face of
        the grid. In a lightcone, when one of its cells shares a face with a
        masked cell, or lies on the border: the first shell unless it is the
        central sphere, the last shell, the first and the last declination
        ring of a shell unless its edge is on a pole, and the first and the
        last sector of a ring that is not a full circle."""
    @property
    def contamination(self) -> NDArray[np.float64]:
        """Fraction of the basin's surface on the boundary, shape
        (n_basins,): in a box, the faces of its cells on the faces of the
        grid, over those faces plus the faces its cells share with cells of
        other basins; 1 for a basin that fills the grid. In a lightcone, the
        area its cells share with masked cells plus their area on the
        border, over that area plus the area they share with cells of other
        basins; 0 when that sum is not positive."""
    @property
    def n_tracers(self) -> NDArray[np.uint64]:
        """Tracers starting in the cells of the core, the count behind
        ``delta_av``, shape (n_basins,). Shape (0,) after a field search."""
    @property
    def field_centre(self) -> NDArray[np.float64]:
        """After ``find_voids_box_from_field`` only, shape (n_basins,): the
        field the watershed split, smoothed when a smoothing was applied, at
        ``centre``, interpolated trilinearly between the eight nearest cell
        centres, the point first clamped per axis to the range of the
        centres; at a cell centre, that cell's value. Shape (0,) otherwise."""


class DivergenceGrid:
    """The grid of the void finder and the divergence field on it.

    Read-only: it comes from ``Result.divergence`` (with
    ``Config.keep_divergence``), ``box_grid`` or ``io.read_divergence``, and
    can be passed back to ``find_voids_box`` as ``grid``, bit for bit. A
    reuse reads only the fields that define the cells (``cell_size``,
    ``lower``, ``upper``, ``n_cells_axis``, ``n_cells``), ``cell_size_mps``
    and ``mps``, checks that they are consistent, and raises ``Error``
    otherwise; the per-cell results, the limits and the realization counts
    of the search that filled the grid are ignored.

    Box: cubic cells of side ``cell_size``; cell (i, j, k) spans lower +
    (i, j, k) * cell_size to lower + (i + 1, j + 1, k + 1) * cell_size. The
    grid is built on the box limits of the search, or without them on the
    bounding box of the tracers read and of the endpoints of their valid
    displacements; from a reconstruction, on its box
    (``otswap.Result.box``). It is padded by 0.1 mps on each side. With E the
    padded extent and s the cell side, an axis needs floor(E / s + 1) cells;
    the largest count is used on all three axes, so the grid is a cube of
    cells even when the box is not. On each axis, with n the count and
    S = n s - E the slack, the padded box's lower face moves down by S / 2,
    which is ``lower``, and its upper face up by S / 4, which is ``upper``.
    With limits the grid overhangs them.

    Per-cell arrays have shape ``n_cells_axis`` (``centre``: plus 3), C
    order, so that ``a[i, j, k]`` is the cell (i, j, k); shape (0,) when not
    filled (every per-cell array of a ``box_grid`` grid). The lightcone
    fields are NaN, (0,) or None in a box.
    """

    @property
    def geometry(self) -> Literal["box", "lightcone"]:
        """The geometry of the grid."""
    @property
    def angle_unit(self) -> Optional[Literal["deg", "rad"]]:
        """The unit of the angle fields (``ra_min`` to ``dec_max``, the
        ``cell_ra_*`` and ``cell_dec_*`` arrays, ``centre_sky``): that of the
        call that made the grid, or of ``io.read_divergence``; None for a
        box. In degrees the fields are copies made with otswap's factor, so
        they match ``numpy.rad2deg``; the library holds radians."""
    @property
    def n_cells(self) -> int:
        """Number of cells."""
    @property
    def cell_size(self) -> float:
        """Defines the cells: the cell side, in coordinate units,
        ``cell_size_mps * mps``. After ``find_voids_box_from_field``, the
        spacing of the centres inferred along the axis with the most nodes.
        NaN in a lightcone."""
    @property
    def lower(self) -> tuple[float, float, float]:
        """Defines the cells: x, y, z of the lower grid face, in coordinate
        units."""
    @property
    def upper(self) -> tuple[float, float, float]:
        """Defines the cells: x, y, z of the upper face of the padded box as
        placed, in coordinate units. It lies below the last cell face,
        lower + n_cells_axis * cell_size; it is stored as built."""
    @property
    def n_cells_axis(self) -> tuple[int, int, int]:
        """Defines the cells: (nx, ny, nz), the cells along x, y, z; equal on
        the three axes but after ``find_voids_box_from_field``, where they
        may differ."""
    @property
    def limits(self) -> Optional[otswap.BoxLimits]:
        """The box limits of the search that filled the grid; None without
        limits. A record: a search reusing the grid does not read it, its
        own limits apply."""
    @property
    def n_realizations(self) -> int:
        """The realizations of the search that filled the grid (keyword
        NREC). A record; 0 after ``find_voids_box_from_field`` and in a
        ``box_grid`` grid."""
    @property
    def n_realizations_mean(self) -> int:
        """``Config.n_realizations_mean`` of the search that filled the grid
        (keyword NRECMEAN). A record; 0 after ``find_voids_box_from_field``
        and in a ``box_grid`` grid."""
    @property
    def distance_min(self) -> float:
        """Lightcone only; NaN in a box."""
    @property
    def distance_max(self) -> float:
        """Lightcone only; NaN in a box."""
    @property
    def ra_min(self) -> float:
        """Lightcone only, in ``angle_unit``; NaN in a box."""
    @property
    def ra_max(self) -> float:
        """Lightcone only, in ``angle_unit``; NaN in a box."""
    @property
    def dec_min(self) -> float:
        """Lightcone only, in ``angle_unit``; NaN in a box."""
    @property
    def dec_max(self) -> float:
        """Lightcone only, in ``angle_unit``; NaN in a box."""
    @property
    def cell_size_mps(self) -> float:
        """``Config.divergence_cell_size`` used, in units of the mps: a
        record of the choice beside ``cell_size``, which a search reusing the
        grid reports in ``Result.config``. NaN after a field search given no
        mps."""
    @property
    def mps(self) -> float:
        """The mean particle separation the grid was built with, in
        coordinate units; the smoothing in mps and the ``delta_2mps`` radius
        of a search that reuses the grid follow it. NaN after a field search
        given no mps."""
    @property
    def mps_profile(self) -> Optional[otswap.MpsProfile]:
        """Lightcone only; None in a box."""
    @property
    def cell_ra_min(self) -> NDArray[np.float64]:
        """Lightcone only, in ``angle_unit``, shape (n_cells,); (0,) in a
        box and in a grid read from a file."""
    @property
    def cell_ra_max(self) -> NDArray[np.float64]:
        """Lightcone only, in ``angle_unit``, shape (n_cells,); (0,) in a
        box and in a grid read from a file."""
    @property
    def cell_dec_min(self) -> NDArray[np.float64]:
        """Lightcone only, in ``angle_unit``, shape (n_cells,); (0,) in a
        box and in a grid read from a file."""
    @property
    def cell_dec_max(self) -> NDArray[np.float64]:
        """Lightcone only, in ``angle_unit``, shape (n_cells,); (0,) in a
        box and in a grid read from a file."""
    @property
    def cell_distance_min(self) -> NDArray[np.float64]:
        """Lightcone only; shape (0,) in a box."""
    @property
    def cell_distance_max(self) -> NDArray[np.float64]:
        """Lightcone only; shape (0,) in a box."""
    @property
    def centre(self) -> NDArray[np.float64]:
        """Cartesian centre of each cell, in coordinate units, shape (nx, ny,
        nz, 3)."""
    @property
    def centre_sky(self) -> NDArray[np.float64]:
        """Lightcone only: right ascension and declination of each cell
        centre in ``angle_unit``, and the redshift, shape (n_cells, 3); shape
        (0, 3) in a box."""
    @property
    def divergence(self) -> NDArray[np.float64]:
        """The divergence of each cell, before smoothing, shape (nx, ny, nz);
        what the watershed splits with sigma = 0."""
    @property
    def divergence_smoothed(self) -> NDArray[np.float64]:
        """The divergence after smoothing (see ``Smoothing``), what the
        watershed splits, shape (nx, ny, nz), with sigma > 0; shape (0,)
        with sigma = 0."""
    @property
    def volume(self) -> NDArray[np.float64]:
        """Volume of each cell, in coordinate units cubed, shape (nx, ny,
        nz); uniform in a box."""
    @property
    def weight(self) -> NDArray[np.float64]:
        """Tracer weight of each cell, the weight of the smoothing, shape
        (nx, ny, nz): n_displacements / n_realizations where a displacement
        contributed. A cell that no displacement contributed to is empty: it
        keeps a weight of 0, and takes the mean divergence of nearby cells
        that hold a value, weighted by 1 / r^2, r the distance between the
        cell centres. Box: the neighbours within one cell along each axis,
        or within two when fewer than three of those hold a value, in passes
        that each read the values left by the previous one, until no cell is
        empty. Lightcone: the neighbouring cells (as ``BasinCatalog``) that
        hold a value, when at least three do; the passes repeat until one
        fills no cell, and a masked cell is never filled. With its weight of
        0, a filled value enters no smoothed mean; with sigma = 0 the
        watershed reads it as it is."""
    @property
    def weight_smoothed(self) -> NDArray[np.float64]:
        """Smoothed tracer weight of each cell, shape (nx, ny, nz), with
        sigma > 0: the mean of ``weight`` over the cells of the smoothing
        (see ``Smoothing``) with the Gaussian alone as weight, sum K w / sum
        K, K = exp(-r^2 / (2 sigma^2)). It shows how much of a cell's
        smoothed divergence comes from displacements; the watershed does not
        read it. The sum of K over the cube is the product of the three sums
        of the 1-D kernel along the axes, each cut at the faces. Shape (0,)
        with sigma = 0."""
    @property
    def n_displacements(self) -> NDArray[np.uint64]:
        """Displacement ends counted in each cell, before the fill, shape
        (nx, ny, nz): a valid displacement that leaves its tracer's cell
        counts 1 there and 1 in the cell it ends in; one that ends in its
        tracer's cell counts 2 there. A point on a face between two cells
        belongs to the upper one, and a displacement that starts or ends on
        a face, an edge or a vertex, or runs along a face, is counted as
        that of the same tracer moved one ulp into its cell. Shape (0,)
        after a field search."""
    @property
    def n_tracers(self) -> NDArray[np.uint64]:
        """Tracers starting in each cell, shape (nx, ny, nz). Shape (0,)
        after a field search."""
    @property
    def void_index(self) -> NDArray[np.int64]:
        """Row of the cell's void in ``Result.voids``, shape (nx, ny, nz);
        -1 none."""
    @property
    def overdense_index(self) -> NDArray[np.int64]:
        """Row of the cell's overdense basin in ``Result.overdense``, shape
        (nx, ny, nz); -1 none."""
    @property
    def masked(self) -> NDArray[np.bool_]:
        """Lightcone only; shape (0,) in a box."""


class TracerAssignment:
    """Where each input tracer ended up, one row per input row, shape (N,)
    each. Filled in ``Result.assignment`` with ``Config.keep_assignment``.
    Read-only."""

    @property
    def cell_index(self) -> NDArray[np.int64]:
        """Flat row of the cell the tracer starts in; -1 when outside the grid
        or left out. ``numpy.unravel_index(row, grid.n_cells_axis)`` gives
        (i, j, k)."""
    @property
    def void_index(self) -> NDArray[np.int64]:
        """Row of the tracer's void in ``Result.voids``; -1 none."""
    @property
    def overdense_index(self) -> NDArray[np.int64]:
        """Row of the tracer's overdense basin in ``Result.overdense``; -1
        none."""
    @property
    def distance(self) -> NDArray[np.float64]:
        """Distance from the centre of the tracer's basin, its void or its
        overdense basin, in coordinate units; NaN for a tracer in no basin
        (not read by the search)."""
    @property
    def in_mask(self) -> NDArray[np.bool_]:
        """False when the tracer falls on an unobserved pixel of the mask.
        Shape (0,) in a box."""
    @property
    def in_cut(self) -> NDArray[np.bool_]:
        """False when the tracer lies outside the redshift cut, and so was
        skipped. Shape (0,) in a box."""
    @property
    def in_box(self) -> NDArray[np.bool_]:
        """False for a tracer outside the box limits, which was skipped; all
        True without limits. Shape (0,) in a lightcone."""
    @property
    def used(self) -> NDArray[np.bool_]:
        """True when at least one valid displacement of the tracer entered
        the divergence. Shape (0,) after ``find_voids_box_from_divergence``,
        which reads no displacement."""


class Counts:
    """What the void finder counted: the counts its report prints. Read-only.
    Every count is 0 in a field result."""

    @property
    def tracers(self) -> int:
        """Tracers given."""
    @property
    def tracers_used(self) -> int:
        """Of them, with at least one displacement used."""
    @property
    def tracers_uncorrected(self) -> int:
        """Real-space pass of ``run_box`` with a correction: left unmoved by
        the correction, and excluded (``Result.uncorrected``)."""
    @property
    def tracers_outside_box(self) -> int:
        """With box limits: outside them, and skipped
        (``TracerAssignment.in_box``)."""
    @property
    def tracers_outside_grid(self) -> int:
        """Reused grid: outside it, and excluded (``Result.outside_grid``).
        In the real-space pass of ``run_box``: corrected, but outside the
        first pass's grid."""
    @property
    def tracers_outside_cut(self) -> int:
        """Real-space pass of ``run_lightcone`` with a correction: corrected,
        but moved outside the redshift cut, and excluded
        (``Result.outside_cut``). 0 in a box."""
    @property
    def tracers_outside_mask(self) -> int:
        """Lightcone: tracers on unobserved pixels of the mask, and so
        skipped (``TracerAssignment.in_mask``), whether or not they lie
        inside the redshift cut. 0 in a box."""
    @property
    def tracers_outside_redshift_cut(self) -> int:
        """Lightcone: of the tracers on observed pixels, those outside the
        redshift cut, and so skipped (``TracerAssignment.in_cut``). A tracer
        both masked and outside the cut counts only in
        ``tracers_outside_mask``. 0 in a box."""
    @property
    def randoms(self) -> int:
        """Lightcone: the randoms the search reads (N_r of ``delta_2mps``),
        kept by the mask and the redshift cut and inside the grid. 0 in a
        box."""
    @property
    def randoms_outside_mask(self) -> int:
        """Lightcone: randoms on unobserved pixels of the mask, removed.
        0 in a box."""
    @property
    def randoms_outside_redshift_cut(self) -> int:
        """Lightcone: of the randoms on observed pixels, those outside the
        redshift cut, removed. 0 in a box."""
    @property
    def tracers_in_masked_cells(self) -> int:
        """Lightcone: tracers on an observed pixel whose cell is masked;
        they keep their cell but are in no basin. 0 in a box."""
    @property
    def randoms_outside_grid(self) -> int:
        """Lightcone, reused grid: randoms outside it, removed. 0 in a
        box."""
    @property
    def displacements(self) -> int:
        """Valid displacements given."""
    @property
    def displacements_rejected(self) -> int:
        """Lightcone: valid displacements of the tracers read that were
        rejected (a start or end in a masked cell, an end outside every cell,
        a path too close to the observer). 0 in a box."""


class Result:
    """Output of one void search. Read-only.

    Given the same inputs, every field but ``elapsed_seconds`` is the same
    bit for bit whatever the number of threads.

    A result of ``find_voids_box_from_field``, a field result, holds no
    tracer: ``n_tracers`` and ``n_realizations`` are 0, which tells it from
    the result of every other search, and so are the counts; the per-tracer
    flags, ``assignment``, and the catalogue columns ``delta_av``,
    ``delta_2mps`` and ``n_tracers`` are empty; ``limits`` is None; ``mps``
    is NaN when none was given.
    """

    @property
    def geometry(self) -> Literal["box", "lightcone"]:
        """The geometry of the search."""
    @property
    def angle_unit(self) -> Optional[Literal["deg", "rad"]]:
        """The unit of the sky arrays, that of the call that made the
        result; None for a box."""
    @property
    def config(self) -> Config:
        """The configuration the search ran with, ``divergence_cell_size``
        set to the value used: a new copy on each access, so changing it
        changes nothing in the result."""
    @property
    def n_tracers(self) -> int:
        """Tracers given, one row each; 0 for a field result."""
    @property
    def n_realizations(self) -> int:
        """Realizations of the displacements searched; 0 for a field
        result."""
    @property
    def ids(self) -> NDArray[np.int64]:
        """The identity of each tracer, shape (n_tracers,), as given: a new
        read-only copy on each access; 0..n_tracers-1 when none were set. No
        computation reads it."""
    @property
    def random_seed(self) -> Optional[int]:
        """The seed ``run_box`` derived and used to draw the randoms, which
        alone reproduces them; None when it drew none."""
    @property
    def voids(self) -> BasinCatalog:
        """The voids."""
    @property
    def overdense(self) -> BasinCatalog:
        """The overdense basins."""
    @property
    def divergence(self) -> Optional[DivergenceGrid]:
        """The grid and the divergence field, with
        ``Config.keep_divergence``; None otherwise. With a reused grid, its
        cells with this search's values."""
    @property
    def assignment(self) -> Optional[TracerAssignment]:
        """Where each tracer ended up, with ``Config.keep_assignment``; None
        otherwise."""
    @property
    def uncorrected(self) -> NDArray[np.bool_]:
        """Real-space pass of ``run_box`` with a correction, shape
        (n_tracers,): True for a tracer the correction left unmoved, excluded
        from the second reconstruction and from this search, row kept. Shape
        (0,) elsewhere, which reads as all False."""
    @property
    def outside_grid(self) -> NDArray[np.bool_]:
        """True for a tracer outside the reused grid, excluded from this
        search, row kept, shape (n_tracers,). In the real-space pass of
        ``run_box``: corrected, but outside the first pass's grid. Shape (0,)
        when no grid was reused, which reads as all False. A tracer carries
        at most one of ``uncorrected`` and ``outside_grid``."""
    @property
    def outside_cut(self) -> NDArray[np.bool_]:
        """Lightcone only; shape (0,) in a box."""
    @property
    def limits(self) -> Optional[otswap.BoxLimits]:
        """The box limits of the search: those given, the reconstruction's,
        or the grid's for ``find_voids_box_from_divergence``; the tracers
        outside them were skipped. None without limits."""
    @property
    def mps(self) -> float:
        """The mean particle separation the search used, in coordinate
        units: the one given, the reconstruction's, or a reused grid's. NaN
        for a field result given none."""
    @property
    def mps_profile(self) -> Optional[otswap.MpsProfile]:
        """Lightcone only; None in a box."""
    @property
    def mask_nside(self) -> int:
        """Lightcone: the NSIDE of the mask the search was given; 0 in a
        box."""
    @property
    def redshift_cut(self) -> Optional[otswap.RedshiftCut]:
        """Lightcone: the closed range of observed redshift the search
        kept; None in a box."""
    @property
    def distances(self) -> Optional[otswap.DistanceTable]:
        """Lightcone only; None in a box."""
    @property
    def counts(self) -> Counts:
        """What the search counted."""
    @property
    def elapsed_seconds(self) -> float:
        """Wall time of the call, in seconds; the only field that differs
        between two runs on the same inputs."""


class Partition:
    """The basins of a field, as ``watershed_box`` finds them, with no
    catalogue. Read-only."""

    @property
    def n_voids(self) -> int:
        """Number of voids."""
    @property
    def n_overdense(self) -> int:
        """Number of overdense basins."""
    @property
    def void_index(self) -> NDArray[np.int64]:
        """Row of each cell's void, -1 none; int64, the shape of the field,
        (nx, ny, nz)."""
    @property
    def overdense_index(self) -> NDArray[np.int64]:
        """Row of each cell's overdense basin, -1 none; int64, the shape of
        the field."""
    @property
    def void_seed(self) -> NDArray[np.int64]:
        """Flat row of each void's seed cell, ascending, shape (n_voids,);
        ``numpy.unravel_index(row, field.shape)`` gives (i, j, k)."""
    @property
    def overdense_seed(self) -> NDArray[np.int64]:
        """Flat row of each overdense basin's seed cell, ascending, shape
        (n_overdense,)."""
    @property
    def void_cells(self) -> NDArray[np.uint64]:
        """Cells of each extended void, shape (n_voids,)."""
    @property
    def void_core_cells(self) -> NDArray[np.uint64]:
        """Cells of each void's core, where the field is below 0, shape
        (n_voids,)."""
    @property
    def overdense_cells(self) -> NDArray[np.uint64]:
        """Cells of each overdense basin, shape (n_overdense,)."""
    @property
    def overdense_core_cells(self) -> NDArray[np.uint64]:
        """Cells of each overdense basin's core, where the field is at or
        above 0, shape (n_overdense,)."""


class BitVFResult:
    """Output of ``run_box`` and ``run_lightcone``. Read-only.

    The otswap objects it holds are made in C++: a lightcone's sky arrays
    (``tracers_sky``, ``lagrangian_sky``, ``sky``) are in radians whatever
    the ``angle_unit`` of the call, and their ``angle_unit`` reads
    ``"rad"``. ``voids`` and ``zspace_voids`` report the call's unit.

    Given the same inputs and a seed other than 0, every field but the
    searches' ``elapsed_seconds`` is the same bit for bit whatever the number
    of threads.
    """

    @property
    def voids(self) -> Result:
        """The voids of the run: of the real-space pass with a correction, of
        the only pass without. Its per-tracer columns have one row per input
        tracer, in input order."""
    @property
    def zspace_voids(self) -> Optional[Result]:
        """The voids of the redshift-space pass, with
        ``BoxRsd.zspace_voids`` or ``LightconeRsd.zspace_voids``; None
        otherwise."""
    @property
    def reconstruction(self) -> Optional[otswap.Result]:
        """With ``keep_reconstructions``: the reconstruction the voids were
        found on. With a correction its rows are the second reconstruction's
        tracers only, and its ids are their identity among the input's."""
    @property
    def zspace_reconstruction(self) -> Optional[otswap.Result]:
        """With ``keep_reconstructions`` and a correction: the first,
        redshift-space, reconstruction, one row per input tracer."""
    @property
    def real_space(self) -> Optional[otswap.RealSpaceCatalog]:
        """With ``keep_reconstructions`` and a correction: the corrected
        catalogue, one row per input tracer."""


# --------------------------------------------------------------------------
# functions

def find_voids_box(tracers: Union[ArrayLike, otswap.Result], displacement: Optional[ArrayLike] = None,
                   valid: Optional[ArrayLike] = None, *, mps: Optional[float] = None, box: Optional[Box] = None,
                   config: Config = Config(), ids: Optional[ArrayLike] = None,
                   grid: Optional[DivergenceGrid] = None) -> Result:
    """Find voids and overdense basins in a box.

    The divergence of the displacement field is computed on a grid of cubic
    cells (``DivergenceGrid``), smoothed (``Config.smoothing``) and divided
    by a watershed into basins around its local minima: the voids, below 0,
    and the overdense basins, at or above 0 (``BasinCatalog``).

    Two forms. With arrays, ``tracers`` and ``displacement`` are required,
    and so is ``mps``. With an ``otswap.Result`` of
    ``otswap.reconstruct_box`` as the first argument, the tracers,
    displacements, validity, mps, box and ids are read from it, and none of
    ``displacement``, ``valid``, ``mps``, ``box`` and ``ids`` may be given:
    the grid is built on ``otswap.Result.box``, and the tracers flagged in
    ``otswap.Result.outside_box`` are skipped. Given the same arrays and
    ``otswap.Result.box`` as ``box``, both forms give the same result.

    The ``otswap.Result`` may be the mean of a reconstruction, read from a
    mean displacement file by ``otswap.io.read_displacement_mean`` (one
    realization, ``n_realizations_mean`` the number averaged; set
    ``Config.n_realizations_mean`` to it, so that bitvf's files record it).
    Its box is the file's box keywords, with or without BOXLIM; a mean file
    without box keywords holds no box limits and is refused. The divergence is not linear in the displacement, so a search
    on the mean of R realizations does not give the result of the search on
    the R realizations.

    Parameters
    ----------
    tracers : array_like, shape (N, 3), or otswap.Result
        Cartesian coordinates of the tracers, the start of their
        displacements, in coordinate units; or a box reconstruction.
    displacement : array_like, shape (R, N, 3)
        Displacements of R realizations, R at least 1; a tracer plus its
        displacement is the random matched to it. Every valid displacement
        of a tracer the search reads must be finite; the others, and those of
        the tracers skipped outside the box, are not read.
    valid : array_like, shape (R, N), optional
        Validity of each displacement, bool or numbers 0 and 1. If omitted,
        every displacement is valid.
    mps : float
        Mean particle separation, in coordinate units; finite and positive.
        Use ``otswap.Result.mps`` of the reconstruction, or pass the
        reconstruction itself. With ``grid``, it must equal ``grid.mps``.
    box : BoxLimits or ((x0, y0, z0), (x1, y1, z1)), optional
        The box limits. Given, the tracers outside them keep their row but
        are skipped (``TracerAssignment.in_box``,
        ``Counts.tracers_outside_box``), and the grid is built on them; pass
        the limits the reconstruction used, or a valid displacement ending
        outside the grid is an error. If omitted, the grid covers the
        tracers and the endpoints of their valid displacements in every
        realization, so that displacements from any source can be searched.
    config : Config, optional
        The parameters of the search.
    ids : array_like of int, shape (N,), optional
        The identity of each tracer, any 64-bit integer, no value twice,
        copied into ``Result.ids``.
    grid : DivergenceGrid, optional
        A box grid to reuse (``Result.divergence``, ``box_grid``,
        ``io.read_divergence``); ``config.divergence_cell_size`` must then be
        None. Its cells are used, and its mps sets the smoothing in mps and
        the ``delta_2mps`` radius; with a reconstruction, the grid's mps wins
        over the reconstruction's. Densities come from this search's own
        counts. Tracers outside every cell are excluded
        (``Result.outside_grid``, ``Counts.tracers_outside_grid``).

    Returns
    -------
    Result
        The voids and overdense basins, ``Result.mps`` the mps used.

    Raises
    ------
    Error
        On the arrays: if ``tracers`` is not of shape (N, 3) or is empty, if
        ``displacement`` is not of shape (R, N, 3) with R at least 1, if
        ``valid`` is given empty, is not of shape (R, N) or holds a value
        other than 0 and 1, if a coordinate of ``tracers`` or a component of
        a displacement read is not finite, if ``ids`` is not one integer per
        tracer or holds a value twice; if ``displacement`` or ``mps`` is left
        out with arrays, or one of ``displacement``, ``valid``, ``mps``,
        ``box`` and ``ids`` is given with a reconstruction.
    Error
        If ``mps`` is not finite and positive; if a bound of ``box`` is not
        finite or a lower bound is not below the upper one.
    Error
        On ``config``: if ``smoothing.sigma`` is not finite and
        non-negative, or is positive but, in coordinate units, below about
        1e-162 or above about 1e154; if ``divergence_cell_size`` is set and
        not finite and positive, or set with ``grid``; if
        ``n_realizations_mean`` is 0.
    Error
        On ``grid``: if it is not a box grid, or a field that defines the
        cells is inconsistent; if ``mps`` differs from ``grid.mps``, bit for
        bit.
    Error
        On the reconstruction: if it is not a box result, or is malformed:
        no objects, ``n_realizations`` or ``n_realizations_mean`` 0, a mean
        (``n_realizations_mean`` above 1) of several realizations,
        ``valid_realizations`` of the wrong length or above
        ``n_realizations`` x ``n_realizations_mean``, or no box.
    Error
        On the data: if the search reads more than 4294967295 tracers; if no
        valid displacement is left to read; with ``box`` or ``grid``, if a
        valid displacement read ends outside every cell of the grid, which
        means the arrays do not come from one reconstruction in that box; if
        the grid would hold more cells than an int64 index can address; if
        the divergence is 0 or positive in every cell, so that no void can
        be found.

    Notes
    -----
    Prints as ``Verbosity`` describes, with the call line
    ``bitvf: findVoidsBox: ...``; with a reconstruction and ``grid``, the
    ``"detailed"`` report states both mps.
    """


def find_voids_box_from_divergence(grid: DivergenceGrid, tracers: ArrayLike, *, config: Config = Config(),
                                   ids: Optional[ArrayLike] = None) -> Result:
    """Find voids and overdense basins in a box from a divergence computed
    earlier, with the tracers given.

    The field before smoothing (``grid.divergence``) is smoothed with
    ``config.smoothing``, in units of the grid's mps, split into basins and
    catalogued as ``find_voids_box`` does; a different sigma is the usual
    reason for the call. With the tracers of the search that produced the
    grid and the same smoothing, the voids, the overdense basins, the
    smoothed divergence, the smoothed weight and the assignment but its
    ``used`` flags equal that search's bit for bit, sigma = 0 included.

    The tracers are read as with a reused grid: those outside the grid's
    box limits (``grid.limits``) are skipped, those inside but outside every
    cell are excluded (``Result.outside_grid``). When the grid holds
    ``n_tracers``, the tracers read in every cell must match it.

    Parameters
    ----------
    grid : DivergenceGrid
        A box grid holding the fields that define its cells, ``divergence``
        and ``weight``: from ``io.read_divergence``, or ``Result.divergence``
        of a search with ``keep_divergence``. Its other per-cell fields are
        not read: ``divergence_smoothed`` and ``weight_smoothed`` of the
        result are computed from these with ``config.smoothing``.
    tracers : array_like, shape (N, 3)
        Cartesian coordinates of the tracers, in coordinate units; at least
        one.
    config : Config, optional
        The parameters of the search; ``divergence_cell_size`` must be None.
    ids : array_like of int, shape (N,), optional
        The identity of each tracer, no value twice.

    Returns
    -------
    Result
        No displacement enters the call: ``TracerAssignment.used`` is empty,
        ``Counts.displacements`` and ``Counts.tracers_used`` are 0. The
        grid's ``n_displacements``, ``n_realizations`` and
        ``n_realizations_mean`` are copied into the result, so the files
        written from it record those of the field; ``Result.limits`` is the
        grid's limits, ``Result.mps`` its mps,
        ``Result.config.divergence_cell_size`` its ``cell_size_mps``.

    Raises
    ------
    Error
        If ``tracers`` is not of shape (N, 3), is empty, or a coordinate is
        not finite; if ``ids`` is not one integer per tracer or holds a value
        twice.
    Error
        On ``config`` as ``find_voids_box``; if ``divergence_cell_size`` is
        set.
    Error
        On ``grid`` as on a grid to reuse in ``find_voids_box``; if
        ``divergence`` or ``weight`` does not hold one entry per cell, or an
        entry of either is not finite; if its limits are not valid.
    Error
        On the tracers: if the search reads more than 4294967295 tracers; if
        none is left to read; if the count of a cell differs from the grid's
        ``n_tracers``. If the smoothed divergence is 0 or positive in every
        cell.

    Notes
    -----
    Prints as ``Verbosity`` describes, the grid line as "outside the
    divergence grid", the call line ``bitvf: findVoidsBoxFromDivergence:
    ...``.
    """


def find_voids_box_from_field(centres: ArrayLike, values: ArrayLike, *, config: Config = Config(),
                              weights: Optional[ArrayLike] = None, mps: Optional[float] = None) -> Result:
    """Find voids and overdense basins in a field given at the centres of
    the cells of a box lattice, without tracers.

    The lattice is inferred from the centres. On each axis, the coordinates
    are sorted and split into nodes at the gaps larger than half the largest
    gap; with n nodes from min to max, the spacing is (max - min) / (n - 1).
    The spacings of the axes with more than one node must agree within 1e-6
    relative: the cell side s is the spacing of the axis with the most
    nodes, x before y before z on a tie. Each centre must lie within 1e-6 s
    of its node on every axis, and each node of the nx * ny * nz must hold
    exactly one centre. The grid spans min - s / 2 to max + s / 2 on each
    axis, without padding; nx, ny, nz need not be equal. The cells of the
    result are the cells of the lattice, not the rows of the input, and the
    order of the input does not change the result.

    The field is smoothed with ``config.smoothing``, as the divergence of
    ``find_voids_box``, with ``weights`` as the tracer weight, and split into
    basins at 0: seeds below 0 start voids, those at or above 0 overdense
    basins. A field with another convention, such as a density or a density
    contrast with voids above 0, must be shifted or flipped before the call.

    Parameters
    ----------
    centres : array_like, shape (n, 3)
        Cartesian coordinates of the cell centres, in any order.
    values : array_like, shape (n,)
        The field at each centre, in the order of ``centres``.
    config : Config, optional
        The parameters of the search. ``divergence_cell_size`` must be None,
        ``keep_assignment`` False and ``n_realizations_mean`` 1: none applies
        without tracers. Without ``mps``, a positive ``smoothing.sigma`` must
        be in ``LengthUnit.COORDINATE``.
    weights : array_like, shape (n,), optional
        The weight of each cell in the smoothing, non-negative, in the order
        of ``centres``, used as given: a cell whose weights sum to at most
        1e-12 times the largest weight takes 0 (see ``Smoothing``). If
        omitted, every cell has weight 1.
    mps : float, optional
        A mean particle separation, in coordinate units. Given, it is the
        unit of a sigma in ``LengthUnit.MPS`` and is recorded in
        ``Result.mps``, with the cell size in units of it in
        ``Result.config.divergence_cell_size``.

    Returns
    -------
    Result
        A field result (see ``Result``). The catalogues hold ``centre``,
        ``radius``, ``radius_ext``, ``edge_flag`` and ``contamination`` as
        ``find_voids_box`` computes them, on the faces of the inferred grid,
        and ``field_centre``. With ``Config.keep_divergence``,
        ``Result.divergence`` holds the grid, the field as given
        (``divergence``) and, with sigma > 0, smoothed
        (``divergence_smoothed``), the weights as given and, with sigma > 0,
        their smoothed mean (``weight_smoothed``), the basin indices, the
        cell centres and volumes.

    Raises
    ------
    Error
        On the arrays: if ``centres`` is not of shape (n, 3), is empty, or a
        coordinate is not finite; if ``values`` does not hold one finite
        entry per centre; if ``weights`` is given empty, does not hold one
        entry per centre, or an entry is negative or not finite.
    Error
        If ``mps`` is given and not finite and positive.
    Error
        On ``config`` as ``find_voids_box``; if ``divergence_cell_size`` is
        set, ``keep_assignment`` is set or ``n_realizations_mean`` is not 1;
        if a positive ``smoothing.sigma`` is in ``LengthUnit.MPS`` without
        ``mps``; if ``smoothing.sigma`` is positive and every weight is 0.
    Error
        On the lattice: if every axis has one node; if the spacings of two
        axes differ by more than 1e-6 relative; if a centre lies farther than
        1e-6 s from its node; if a node holds two centres, or none. If the
        field, smoothed when a smoothing is set, is 0 or positive in every
        cell.

    Notes
    -----
    A bitvf divergence file is valid input, read with ``io.read_field``: its
    cell centres, ``divergence`` and ``weight`` with the smoothing
    of the search that wrote it give back that search's basins. Their
    ``centre``, ``radius`` and ``radius_ext`` can differ from that search's
    by a few ulp: the grid's lower face and s are inferred from the centres
    read and may differ from the original grid's by rounding.

    Prints the call line ``bitvf: findVoidsBoxFromField: ...``.
    """


def box_grid(reconstruction: Optional[otswap.Result] = None, *, box: Optional[Box] = None,
             mps: Optional[float] = None, cell_size: Optional[float] = None) -> DivergenceGrid:
    """The grid ``find_voids_box`` builds, bit for bit, so that several
    searches can share it.

    Two forms: ``box_grid(box=..., mps=..., cell_size=...)``, all three
    required, the grid built on box limits; or ``box_grid(reconstruction)``,
    the grid built on ``otswap.Result.box`` with ``otswap.Result.mps``.

    Only the fields that define the cells are set: ``geometry``, ``lower``,
    ``upper``, ``n_cells_axis``, ``n_cells``, ``cell_size``,
    ``cell_size_mps`` and ``mps``. ``limits``, the realization counts and
    every per-cell field are left None, 0 or empty: records of a search,
    which a reuse ignores. The grid passes the checks of a grid to reuse.

    Parameters
    ----------
    reconstruction : otswap.Result, optional
        A box result of ``otswap.reconstruct_box``.
    box : BoxLimits or ((x0, y0, z0), (x1, y1, z1)), optional
        The box limits the grid is built on, as the ``box`` of
        ``find_voids_box``. Not with ``reconstruction``.
    mps : float, optional
        Mean particle separation, in coordinate units; finite and positive.
        Not with ``reconstruction``.
    cell_size : float, optional
        The cell side in units of mps, as ``Config.divergence_cell_size``;
        finite and positive. With ``reconstruction``, if omitted, the
        default for its realizations, min(1, (20 / n_realizations)^(1/3)).

    Returns
    -------
    DivergenceGrid

    Raises
    ------
    Error
        If neither a reconstruction nor all of ``box``, ``mps`` and
        ``cell_size`` are given, or ``box`` or ``mps`` is given with a
        reconstruction.
    Error
        If a bound of ``box`` is not finite or a lower bound is not below the
        upper one; if ``mps`` or ``cell_size`` is not finite and positive; if
        the grid would hold more cells than an int64 index can address.
    Error
        If the reconstruction is not a box result, has no realization, its
        box is missing or invalid, or its mps is not finite and positive.
    """


def watershed_box(field: ArrayLike, *, verbosity: Verbosity = "normal") -> Partition:
    """Split a field on a box grid into basins by steepest descent.

    The descent is that of the void finder (``BasinCatalog``): each cell
    flows to the neighbour, among the 26 cut at the faces of the grid, where
    the field drops the most, ties to the lowest row; a cell with no
    strictly smaller neighbour is a seed. Seeds below 0 start voids, seeds
    at or above 0 overdense basins; each kind is numbered by ascending seed
    row. A field with another convention must be shifted or flipped before
    the call. No smoothing, no fill and no centres: the field is used as
    given. A field with no cell below 0 gives no void, and no error. The
    result does not depend on the number of threads.

    Parameters
    ----------
    field : array_like, shape (nx, ny, nz)
        The value of each cell, in any memory order; each count at least 1,
        not necessarily equal.
    verbosity : {"silent", "normal", "detailed"}, optional
        At ``"normal"``, one line with its time::

            bitvf: watershedBox: <V> voids and <O> overdense basins on <nx> x <ny> x <nz> cells in <T> s

        and at ``"detailed"``, before it, the cells in voids and in
        overdense basins, and how many of each lie in the cores.

    Returns
    -------
    Partition

    Raises
    ------
    Error
        If ``field`` is not three-dimensional, a count is 0, or an entry is
        not finite.
    """


def run_box(tracers: ArrayLike, randoms: Optional[ArrayLike] = None, *, config: BoxConfig = BoxConfig(),
            distances: Optional[otswap.DistanceTable] = None, ids: Optional[ArrayLike] = None) -> BitVFResult:
    """Reconstruct a box and find its voids, with the correction of
    redshift-space distortions when ``config.rsd`` is set.

    The run:

    1. Checks every argument and every field of ``config`` that can be
       checked before any work, so that an invalid run fails before the
       reconstruction.
    2. Fixes the box: ``config.limits``, or the bounding box of the tracers,
       or of the tracers and the randoms when randoms are given.
    3. Without randoms, draws ``config.ot.n_realizations`` * N of them,
       uniformly in the box, N the tracers inside it.
    4. Reconstructs with ``otswap.reconstruct_box``: the randoms, the box as
       limits, ``config.mps`` when given.
    Without ``config.rsd``, it then searches that reconstruction with
    ``find_voids_box``, which reads its mps. With ``config.rsd``, after
    step 4:

    5. With ``BoxRsd.zspace_voids``, searches the first reconstruction
       (``BitVFResult.zspace_voids``).
    6. Builds the grid of the first reconstruction, ``box_grid(first,
       cell_size=config.voids.divergence_cell_size)``: the grid of every
       search of the run.
    7. Corrects the tracers with ``otswap.real_space_box``.
    8. Selects, in input order, the tracers corrected whose corrected
       position lies in a cell of that grid. The others keep their row and
       are flagged: unmoved by the correction (including the tracers outside
       the box) in ``Result.uncorrected``, moved outside every cell in
       ``Result.outside_grid``.
    9. Reconstructs the selected tracers at their corrected positions, with
       the randoms of the first reconstruction, the same box as limits,
       ``config.mps`` when given, and the same seed.
    10. Searches that reconstruction on the grid of step 6, with the
        smoothing and the ``delta_2mps`` radius from the first
        reconstruction's mps, which is ``Result.mps``. A tracer moved outside
        the box but into a cell is skipped (``TracerAssignment.in_box``).
    11. Returns the search with one row per input tracer.

    The box always reaches otswap as limits, so ``Result.limits`` holds it,
    and the files record it, also when ``config.limits`` is None.

    Parameters
    ----------
    tracers : array_like, shape (N, 3)
        Cartesian coordinates, in coordinate units; in redshift space with a
        correction.
    randoms : array_like, shape (M, 3), optional
        Cartesian coordinates of the randoms, at least
        ``config.ot.n_realizations`` * N inside the box; those outside it are
        dropped by otswap. They go to the reconstructions only. If omitted,
        bitvf draws them: from the seed s = ``config.ot.seed`` + 1299709,
        modulo 2^32, a seed of 0 counting as the one drawn to replace it,
        recorded in ``Result.random_seed`` of every pass and in
        the keyword RANDSEED of bitvf's files; s alone reproduces them.
    config : BoxConfig, optional
        The parameters of the run.
    distances : otswap.DistanceTable, optional
        A table carrying the growth rate at ``BoxRsd.redshift``; required
        with ``config.rsd``, not read without it.
    ids : array_like of int, shape (N,), optional
        The tracers' identity, no value twice; if omitted, the rows. They go
        to the reconstruction and to ``Result.ids``.

    Returns
    -------
    BitVFResult
        The voids, and as ``config`` asks the redshift-space voids, the two
        reconstructions and the corrected catalogue. In ``voids``, a row
        excluded from the second pass carries -1 in every index, NaN in
        ``TracerAssignment.distance``, False in ``used`` and, in ``in_box``,
        the first reconstruction's value (False for a tracer outside the
        box limits, True otherwise), and one flag, ``uncorrected`` before
        ``outside_grid``; ``Result.counts`` counts the input rows.

    Raises
    ------
    Error
        If an argument or a field of ``config`` is invalid, as
        ``find_voids_box`` and ``otswap.reconstruct_box`` require: shapes as
        ``find_voids_box``, ``randoms`` given empty.
    Error
        If ``config.rsd`` is set without ``distances``; if a field of
        ``config.rsd`` is invalid (``axis`` None or not 0, 1 or 2, ``bias``
        not finite and positive, ``sigma`` not finite and non-negative) or
        the distance table cannot give the growth rate at ``redshift``,
        checked before the first reconstruction.
    Error
        If the bounding box spans no volume; if randoms are to be drawn and
        no tracer lies inside the box; with a correction, if no tracer is
        selected for the second reconstruction.
    Error
        From otswap, the message starting with ``"otswap: "``; from each
        search, as ``find_voids_box``.

    Notes
    -----
    Prints otswap's report lines at ``config.ot.verbosity``, the
    correction's at ``BoxRsd.correction.verbosity``, and the searches' and
    the pipeline's at ``config.voids.verbosity``: at ``"detailed"``, after
    drawing the randoms and after step 8::

        bitvf: drew <n> randoms in <box> with seed <s> in <t> s
        bitvf: second pass on <k> of <N> tracers: <u> uncorrected, <g> outside the first pass's grid in <t> s

    and at ``"normal"``, last, with the time of the whole call::

        bitvf: runBox: <V> voids and <O> overdense basins from <N> tracers in <T> s
    """


AngleUnit: TypeAlias = Literal["deg", "rad"]
"""The unit of right ascension and declination: degrees or radians."""

RedshiftCutArg: TypeAlias = Union[otswap.RedshiftCut, tuple[float, float]]
"""A closed range of observed redshift: an ``otswap.RedshiftCut`` or a
``(min, max)`` pair."""


@overload
def find_voids_lightcone(tracers_sky: ArrayLike, displacement: ArrayLike, randoms_sky: ArrayLike,
                         valid: Optional[ArrayLike] = None, *, mask: otswap.Mask,
                         distances: otswap.DistanceTable, mps_profile: otswap.MpsProfile, angle_unit: AngleUnit,
                         config: Config = Config(), redshift_cut: Optional[RedshiftCutArg] = None,
                         ids: Optional[ArrayLike] = None, grid: Optional[DivergenceGrid] = None) -> Result: ...
@overload
def find_voids_lightcone(tracers_sky: otswap.Result, randoms_sky: ArrayLike, *, mask: otswap.Mask,
                         angle_unit: AngleUnit, config: Config = Config(),
                         grid: Optional[DivergenceGrid] = None) -> Result: ...
def find_voids_lightcone(*args, **kwargs):  # type: ignore[no-untyped-def]
    """Find voids and overdense basins in a lightcone, from sky positions or
    from an ``otswap.Result``.

    The divergence of the displacement field is computed on a grid of shells
    in comoving distance, declination rings and right-ascension sectors,
    smoothed (``config.smoothing``) and divided by a watershed into basins
    around its local minima: the voids, below 0, and the overdense basins, at
    or above 0. The observer is at the origin.

    A tracer on an unobserved pixel of the mask, or with a redshift outside
    the cut, keeps its row but is skipped (``TracerAssignment.in_mask`` and
    ``in_cut``, ``Counts.tracers_outside_mask`` and
    ``tracers_outside_redshift_cut``); a random there is removed. A cell
    whose centre falls on an unobserved pixel is masked
    (``DivergenceGrid.masked``): it holds no divergence and belongs to no
    basin. Without ``grid``, the grid covers the tracers and randoms kept,
    with cells ``config.divergence_cell_size`` times mps(z) wide at the
    centre of each shell, by default ``min(1, (20 / R)**(1/3))`` for R
    realizations. Any finite right ascension is accepted and folded into
    [0, 2 pi).

    Parameters
    ----------
    tracers_sky : array_like, shape (N, 3)
        Right ascension and declination in ``angle_unit``, and the observed
        redshift, which must lie inside the distance table. A declination
        must lie in [-90, 90] degrees.
    displacement : array_like, shape (R, N, 3)
        Cartesian displacements of R realizations, in the unit of ``distances``, from the
        tracers' Cartesian positions (``otswap.to_cartesian`` of
        ``tracers_sky``). Every valid displacement of a tracer kept must be
        finite.
    randoms_sky : array_like, shape (M, 3)
        The randoms, as ``tracers_sky``: required. They describe the survey,
        set the extent of a new grid with the tracers, and give the expected
        counts of ``delta_2mps``. Pass those of the reconstruction. With an
        ``otswap.Result``, the second positional argument.
    valid : array_like, shape (R, N), optional
        0 or 1 for each displacement; if omitted, all valid.
    mask : otswap.Mask
        The observed pixels of the survey; that of the reconstruction, or a
        stricter one.
    distances : otswap.DistanceTable
        The redshift-distance table; coordinates are in its distance unit.
        Required with arrays; an ``otswap.Result`` carries its own.
    mps_profile : otswap.MpsProfile
        The mean particle separation as a function of redshift
        (``otswap.Result.mps_profile`` of the reconstruction): the cell size
        of each shell, a smoothing in mps, the densities. Required with
        arrays. With ``grid`` it must equal the grid's profile.
    angle_unit : {"deg", "rad"}
        The unit of ``tracers_sky`` and ``randoms_sky``, and of the sky
        arrays of the result.
    config : Config, optional
        The parameters of the search.
    redshift_cut : otswap.RedshiftCut or (min, max), optional
        The closed range of observed redshift searched; pass the cut of the
        reconstruction. Arrays only; the default cuts nothing.
    ids : array_like of int, shape (N,), optional
        The identity of each tracer; arrays only. As ``find_voids_box``.
    grid : DivergenceGrid, optional
        A lightcone grid to reuse, built with the same distance table, whose
        mps profile sets the smoothing and the densities. Tracers kept but
        outside every cell are excluded (``Result.outside_grid``).
        ``config.divergence_cell_size`` must be None with it. It may have
        been made in the other ``angle_unit``.

    Returns
    -------
    Result
        With ``Result.geometry`` ``"lightcone"`` and ``Result.angle_unit``
        that of the call.

    Raises
    ------
    Error
        If an argument has the wrong type or shape, a declination lies
        outside [-90, 90] degrees, ``angle_unit`` is neither ``"deg"`` nor
        ``"rad"``, an argument that arrays need is missing, or one that the
        form with an ``otswap.Result`` does not take is given.
    Error
        On the data and the configuration, as the library requires: no
        tracer, valid displacement or random left after the mask, the cut
        and a reused grid; a redshift outside the distance table; a grid that
        fails the checks of a grid to reuse; a divergence that is nowhere
        negative.
    """


def find_voids_lightcone_from_cartesian(tracers: ArrayLike, displacement: ArrayLike, randoms: ArrayLike,
                                        valid: Optional[ArrayLike] = None, *, mask: otswap.Mask,
                                        distances: otswap.DistanceTable, mps_profile: otswap.MpsProfile,
                                        angle_unit: AngleUnit, config: Config = Config(),
                                        redshift_cut: Optional[RedshiftCutArg] = None,
                                        ids: Optional[ArrayLike] = None,
                                        grid: Optional[DivergenceGrid] = None) -> Result:
    """Find voids and overdense basins in a lightcone, from Cartesian
    positions with the observer at the origin.

    As ``find_voids_lightcone`` from sky positions. The sky coordinates of
    every row are computed first with ``otswap.to_sky`` and ``distances``;
    they alone decide the mask, the redshift cut and the cell of every
    tracer and random. Every row must be finite, must not lie at the origin
    and must lie within the distance table, rows on unobserved pixels or
    outside the cut included.

    Parameters
    ----------
    tracers : array_like, shape (N, 3)
        x, y, z, in the unit of ``distances``.
    displacement, valid, mask, distances, mps_profile, config, redshift_cut, ids, grid
        As ``find_voids_lightcone``.
    randoms : array_like, shape (M, 3)
        As ``tracers``, for the randoms: required.
    angle_unit : {"deg", "rad"}
        The unit of the sky arrays of the result; no argument is a sky
        array.

    Returns
    -------
    Result
        As ``find_voids_lightcone``.
    """


@overload
def lightcone_grid(tracers_sky: otswap.Result, randoms_sky: ArrayLike, *, mask: otswap.Mask,
                   angle_unit: AngleUnit, cell_size: Optional[float] = None) -> DivergenceGrid: ...
@overload
def lightcone_grid(tracers_sky: ArrayLike, randoms_sky: ArrayLike, *, mask: otswap.Mask,
                   distances: otswap.DistanceTable, mps_profile: otswap.MpsProfile, cell_size: float,
                   angle_unit: AngleUnit, redshift_cut: Optional[RedshiftCutArg] = None) -> DivergenceGrid: ...
def lightcone_grid(*args, **kwargs):  # type: ignore[no-untyped-def]
    """The grid ``find_voids_lightcone`` builds, bit for bit, from an
    ``otswap.Result`` or from arrays, with the same selection of the tracers
    and randoms and the same checks, those of the displacements included
    with a reconstruction.

    One grid can be shared by several searches (``grid``), by the two passes
    of a pipeline or by a set of mocks. Only the fields that define the cells
    are set (``n_cells``, ``distance_min`` to ``dec_max``, ``cell_size_mps``,
    ``mps_profile``); every per-cell array is empty and the realization
    counts are 0.

    Parameters
    ----------
    tracers_sky, randoms_sky, mask, distances, mps_profile, redshift_cut
        As ``find_voids_lightcone``; ``tracers_sky`` may be an
        ``otswap.Result``.
    cell_size : float, optional
        The cell size in units of mps(z), as ``Config.divergence_cell_size``;
        finite and positive. Required with arrays, where a reasonable
        starting value is ``min(1, (20 / R)**(1/3))``, R the number of
        realizations the search will read. With an ``otswap.Result``, empty
        for the default, ``min(1, (20 / R)**(1/3))``, R the realizations of
        the result.
    angle_unit : {"deg", "rad"}
        The unit of ``tracers_sky`` and ``randoms_sky``, and of the angle
        fields of the grid.

    Returns
    -------
    DivergenceGrid
        With ``geometry`` ``"lightcone"``.
    """


def run_lightcone(tracers_sky: ArrayLike, randoms_sky: ArrayLike, randoms_sky_second_pass: Optional[ArrayLike] = None,
                  *, mask: otswap.Mask, distances: otswap.DistanceTable, angle_unit: AngleUnit,
                  config: LightconeConfig, ids: Optional[ArrayLike] = None) -> BitVFResult:
    """Reconstruct a lightcone and find its voids, with the correction of
    redshift-space distortions and a second reconstruction in real space when
    ``config.rsd`` is set.

    The run:

    1. Checks every argument and every field of ``config`` that can be
       checked before any work. With a correction, the factor of the
       correction is computed at the observed redshifts of the tracers kept,
       so that a table without growth rate or an unusable bias fails here.
    2. Reconstructs with ``otswap.reconstruct_lightcone``: the randoms, the
       mask, ``config.n_bins``, the distance table, ``config.ot`` and
       ``config.cut``.
    3. Without ``config.rsd``, searches that reconstruction with
       ``find_voids_lightcone`` and returns.

    With ``config.rsd``, after step 2:

    4. Builds the grid of the first reconstruction with
       ``config.voids.divergence_cell_size``: the grid of every search of
       the run.
    5. With ``LightconeRsd.zspace_voids``, searches the first reconstruction
       on it (``BitVFResult.zspace_voids``).
    6. Corrects the tracers with ``otswap.real_space_lightcone``; the
       correction moves a tracer in redshift only.
    7. Selects, in input order, the tracers corrected whose corrected
       redshift lies inside ``config.cut`` and whose corrected position lies
       in a cell of the grid. The others keep their row and are flagged:
       unmoved by the correction in ``Result.uncorrected``, moved outside
       the cut in ``Result.outside_cut``, moved outside every cell in
       ``Result.outside_grid``.
    8. Reconstructs the selected tracers at their corrected redshifts, with
       the randoms of ``LightconeRsd.second_pass_randoms``, the same mask,
       cut, bins and seed.
    9. Searches that reconstruction on the grid of step 4, with the first
       reconstruction's mps(z) profile for the smoothing and the densities.
    10. Returns the search with one row per input tracer.

    Parameters
    ----------
    tracers_sky : array_like, shape (N, 3)
        Right ascension and declination in ``angle_unit``, and the observed
        redshift, in redshift space.
    randoms_sky : array_like, shape (M, 3)
        The randoms, as ``tracers_sky``: required. Every systematic belongs
        in them. At least ``config.ot.n_realizations`` times the tracers kept
        are needed.
    randoms_sky_second_pass : array_like, shape (K, 3), optional
        The randoms of the second reconstruction, with
        ``SecondPassRandoms.EXTERNAL`` and only then.
    mask : otswap.Mask
        The observed pixels of the survey, for both reconstructions and
        every search.
    distances : otswap.DistanceTable
        The redshift-distance table, carrying the growth rate when
        ``config.rsd`` is set.
    angle_unit : {"deg", "rad"}
        The unit of the sky arrays given and of those of ``BitVFResult.voids``
        and ``zspace_voids``. The otswap objects the result holds are in
        radians whatever it is.
    config : LightconeConfig
        The parameters of the run; required, since ``n_bins`` has no default.
    ids : array_like of int, shape (N,), optional
        The tracers' identity, no value twice; if omitted, the rows.

    Returns
    -------
    BitVFResult
        The voids, and as ``config`` asks the redshift-space voids, the
        reconstructions and the corrected catalogue. In ``voids``, an excluded
        row carries -1 in every index, NaN in ``TracerAssignment.distance``
        and False in ``used``, and one flag, ``uncorrected`` before
        ``outside_cut`` before ``outside_grid``; ``in_mask`` and ``in_cut``
        are the first reconstruction's. ``Result.counts`` counts the input
        rows.

    Raises
    ------
    Error
        If an argument or a field of ``config`` is invalid, as
        ``find_voids_lightcone`` and ``otswap.reconstruct_lightcone``
        require; if ``randoms_sky_second_pass`` is given empty, without a
        correction, or without ``EXTERNAL``; if no tracer is selected at
        step 7; with ``EXTERNAL``, if fewer randoms remain than needed.
    Error
        From otswap, the message starting with ``"otswap: "``.

    Notes
    -----
    Prints otswap's report lines at ``config.ot.verbosity``, the
    correction's at ``LightconeRsd.correction.verbosity``, and the searches'
    and the pipeline's at ``config.voids.verbosity``, the last line, at
    ``"normal"``, being::

        bitvf: runLightcone: <V> voids and <O> overdense basins from <N> tracers in <T> s
    """


__all__: list[str]
