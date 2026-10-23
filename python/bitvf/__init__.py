# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""The Back-in-time Void Finder: voids from the optimal transport
displacement field of a tracer catalogue.

Built on otswap, whose Python types bitvf takes and returns as they are.
"""

# otswap registers its types on import, so it comes before the compiled
# module that uses them.
import otswap  # noqa: F401

from ._bitvf import __version__

__all__ = ["__version__"]
