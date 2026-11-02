# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# =================================================================
# Example code: how to find the voids of a field on a box lattice
# =================================================================

"""This example explains how to find voids in a field given at the centres of
the cells of a box lattice, with no tracers and no reconstruction: a density
contrast, a potential, or any field whose voids lie below 0. A field with
voids above 0 must be flipped first.

The field is read from a table with one row per cell, in any order; the
lattice, its cell side and its extent are inferred from the centres. The
field is smoothed, split by the same watershed as the void finder, and
catalogued: centres, radii, contact with the faces of the grid, and the field
at each centre. Densities, which need tracers, are not computed.

The input here is the divergence file written by the full_recons_box example:
a divergence file of bitvf is a valid field, with the cell centres, the
divergence before smoothing and the tracer weight of each cell. Searched with
the smoothing of the search that wrote it, it gives back that search's
basins; this example smooths it with a Gaussian of 15 Mpc/h instead, wider
than the 1 mps of that search, so that small voids merge into larger ones.

Usage: python field.py [field_file]. The field is read from field_file, by
default the divergence file of the full_recons_box example, which must then
be run first; the output is written to the folder output/ next to this file,
created if missing.
"""

import sys
from pathlib import Path

import bitvf


def main():

    # -------------------------------------------------------------
    # ------------ Locate the input and output folders ------------
    # -------------------------------------------------------------

    # The input is the divergence file of the full_recons_box example,
    # unless a first command-line argument replaces it
    here = Path(__file__).resolve().parent
    default = here.parent / "full_recons_box" / "output" / "divergence_box.dat"
    field_file = Path(sys.argv[1]) if len(sys.argv) > 1 else default
    output = here / "output"
    output.mkdir(parents=True, exist_ok=True)

    # -------------------------------------
    # ------------ Read the field ---------
    # -------------------------------------

    # The columns of the cell centres, in Mpc/h, of the field and of the
    # weight of each cell in the smoothing. Without a weight column every
    # cell weighs the same. The field is divergenceUnsmoothed, not the
    # divergence column: the search smooths the field it reads, so the
    # smoothed divergence would be smoothed twice
    field = bitvf.io.read_field(field_file, centres=["cellX", "cellY", "cellZ"], value="divergenceUnsmoothed",
                                weight="weight")

    print(f"Cells: {len(field.values)}")

    # ---------------------------------------------
    # ------------ Find the voids -----------------
    # ---------------------------------------------

    # No mean particle separation is given, so sigma is in coordinate units,
    # here Mpc/h. "normal" prints one line with the lattice found
    config = bitvf.Config(smoothing=bitvf.Smoothing(sigma=15.0, unit=bitvf.LengthUnit.COORDINATE),
                          verbosity="normal")
    result = bitvf.find_voids_box_from_field(field.centres, field.values, config=config, weights=field.weights)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # field_centre is the smoothed field at each void's centre, so the
    # deepest void has the lowest value
    voids = result.voids

    print(f"Voids: {voids.n_basins}")
    print(f"Largest void radius: {voids.radius.max():g} Mpc/h")
    print(f"Lowest field at a void centre: {min(voids.field_centre.min(), 0.0):g}")
    print(f"Overdense basins: {result.overdense.n_basins}")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # The voids and the overdense basins, voids_field.dat and
    # overdense_basins_field.dat, with the column fieldCentre in place of
    # the densities. The writers are the C++ library's, so the files are the
    # C++ example's, byte for byte
    product = bitvf.io.Product
    bitvf.io.write_voids(bitvf.io.file_name(product.VOIDS, output, "field.dat", zspace=False), result)
    bitvf.io.write_overdense(bitvf.io.file_name(product.OVERDENSE, output, "field.dat", zspace=False), result)

    print("Written: output/voids_field.dat, output/overdense_basins_field.dat")


if __name__ == "__main__":
    main()
