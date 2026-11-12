# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Writing the products of a search, and reading files back.

The format of a file follows its extension: ``.fits``, ``.fit`` and
``.fits.gz`` are FITS, anything else is ASCII. The writing and reading are
the C++ library's, so a file written here is the same, byte for byte, as the
one the C++ library writes from the same result.
"""

from . import _bitvf

Product = _bitvf.io.Product
Field = _bitvf.io.Field
file_name = _bitvf.io.file_name
write_voids = _bitvf.io.write_voids
write_overdense = _bitvf.io.write_overdense
write_divergence = _bitvf.io.write_divergence
write_galaxies = _bitvf.io.write_galaxies
read_divergence = _bitvf.io.read_divergence
read_field = _bitvf.io.read_field

__all__ = [
    "Field",
    "Product",
    "file_name",
    "read_divergence",
    "read_field",
    "write_divergence",
    "write_galaxies",
    "write_overdense",
    "write_voids",
]
