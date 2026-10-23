# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Tests of the package import."""

from importlib.metadata import version

import otswap
import bitvf


# otswap and bitvf import together, and bitvf reports the version it was
# installed as.
def test_import_and_version():
    assert otswap.Config is not None
    assert bitvf.__version__ == version("bitvf")
