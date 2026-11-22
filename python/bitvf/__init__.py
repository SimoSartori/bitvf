# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""The Back-in-time Void Finder: voids from the optimal transport
displacement field of a tracer catalogue.

Built on otswap, whose Python types bitvf takes and returns as they are.
Coordinates are NumPy arrays of shape (N, 3), in coordinate units, the unit
of the tracers given; sky positions are (right ascension, declination,
redshift) in the ``angle_unit`` of the call. Inputs are converted to float64; the library works on
its own copy of them. The number of threads follows ``OMP_NUM_THREADS``.
Every search releases the GIL while it runs, and prints its report to
``sys.stdout`` line by line as it runs; Ctrl+C stops it at its next line.

Files are written and read by the submodule ``bitvf.io``.
"""

# otswap registers its types on import, so it comes before the compiled
# module that uses them.
import otswap  # noqa: F401

from ._bitvf import (
    BasinCatalog,
    BitVFResult,
    BoxConfig,
    BoxRsd,
    Config,
    Counts,
    DivergenceGrid,
    Error,
    LengthUnit,
    LightconeConfig,
    LightconeRsd,
    Partition,
    Result,
    SecondPassRandoms,
    Smoothing,
    TracerAssignment,
    __version__,
    box_grid,
    find_voids_box,
    find_voids_box_from_divergence,
    find_voids_box_from_field,
    find_voids_lightcone,
    find_voids_lightcone_from_cartesian,
    lightcone_grid,
    run_box,
    run_lightcone,
    watershed_box,
)
from . import io

__all__ = [
    "BasinCatalog",
    "BitVFResult",
    "BoxConfig",
    "BoxRsd",
    "Config",
    "Counts",
    "DivergenceGrid",
    "Error",
    "LengthUnit",
    "LightconeConfig",
    "LightconeRsd",
    "Partition",
    "Result",
    "SecondPassRandoms",
    "Smoothing",
    "TracerAssignment",
    "__version__",
    "box_grid",
    "find_voids_box",
    "find_voids_box_from_divergence",
    "find_voids_box_from_field",
    "find_voids_lightcone",
    "find_voids_lightcone_from_cartesian",
    "io",
    "lightcone_grid",
    "run_box",
    "run_lightcone",
    "watershed_box",
]
