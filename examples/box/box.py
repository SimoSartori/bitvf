# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# ===============================================================
# Example code: how to find the voids of a box reconstruction
# ===============================================================

"""This example explains how to find voids in box geometry from an optimal
transport reconstruction. otswap reconstructs the displacement field of halos
in real space; bitvf then computes the divergence of that field on a grid of
cubic cells, smooths it, and splits it by a watershed into basins around its
local minima. Basins whose minimum is below 0, where the reconstructed flow
diverges, are the voids; the others are the overdense basins.

The two steps are separate calls, so the reconstruction can be kept, written,
or searched several times with different parameters. The run_box pipeline of
the full_recons_box example does both in one call.

The box is not periodic: nothing flows through its faces, and the basins that
touch them are flagged (edge_flag), so that they can be left out of an
analysis.

Usage: python box.py [data_dir]. The catalogue is read from data_dir, by
default the examples' data folder; the output is written to the folder
output/ next to this file, created if missing.
"""

import sys
from pathlib import Path

import otswap

import bitvf


def main():

    # -------------------------------------------------------------
    # ------------ Locate the input and output folders ------------
    # -------------------------------------------------------------

    # The data folder is the examples' one, next to this example's folder,
    # unless a first command-line argument replaces it
    here = Path(__file__).resolve().parent
    data = Path(sys.argv[1]) if len(sys.argv) > 1 else here.parent / "data"
    output = here / "output"
    output.mkdir(parents=True, exist_ok=True)

    # ----------------------------------------------------
    # ------------ Read the tracer catalogue -------------
    # ----------------------------------------------------

    # Comoving Cartesian coordinates in Mpc/h, shape (N, 3), in real space
    tracers = otswap.io.read(data / "box_halos_real_space.dat", [0, 1, 2]).values

    print(f"Tracers: {len(tracers)}")

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # The halos fill a cube of side 400 Mpc/h, [269, 669) on each axis.
    # Given as box limits, the cube is where otswap draws the randoms and
    # where bitvf builds its grid; the mean particle separation follows from
    # its volume and the halos inside it
    box = ((269.0, 269.0, 269.0), (669.0, 669.0, 669.0))

    # 8 realizations, a convergence threshold of 1e-2 and a fixed seed, so
    # that the run can be repeated. Without randoms, otswap draws 8 per halo,
    # uniformly in the box
    ot_config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, verbosity="normal")
    reconstruction = otswap.reconstruct_box(tracers, box=box, config=ot_config)

    # ---------------------------------------------
    # ------------ Find the voids -----------------
    # ---------------------------------------------

    # The search reads everything it needs from the reconstruction: the
    # tracers, their displacements in each realization, the mean particle
    # separation (mps) and the box.
    #
    # The two parameters that change the result:
    # - the smoothing: a Gaussian of width sigma, here 1 mps, the default and
    #   a reasonable starting value; a larger sigma merges small basins into
    #   larger ones;
    # - the cell side, left to its default: 1 mps up to 20 realizations,
    #   (20 / n_realizations)^(1/3) mps above.
    #
    # keep_assignment keeps, for each halo, the basin it ended up in, which
    # the galaxies file records. "detailed" prints one line per step
    config = bitvf.Config(smoothing=bitvf.Smoothing(sigma=1.0, unit=bitvf.LengthUnit.MPS),
                          keep_assignment=True, verbosity="detailed")
    result = bitvf.find_voids_box(reconstruction, config=config)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # Each catalogue has one row per basin, in the order of its seed cell.
    # radius is the radius of the sphere with the volume of the void's core,
    # the part where the divergence is below 0; edge_flag is 0 for a void
    # that does not touch the faces of the grid
    voids = result.voids
    inside = int((voids.edge_flag == 0).sum())

    print(f"Mean particle separation: {result.mps:g} Mpc/h")
    print(f"Voids: {voids.n_basins}, {inside} not touching the faces of the grid")
    print(f"Largest void radius: {voids.radius.max():g} Mpc/h")
    print(f"Overdense basins: {result.overdense.n_basins}")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One file per product, ASCII from the extension: the voids and the
    # overdense basins, one row per basin, and where each halo ended up, one
    # row per halo. io.file_name gives the conventional names, voids_box.dat,
    # overdense_basins_box.dat and galaxies_box.dat. The writers are the C++
    # library's, so the files are the C++ example's, byte for byte
    product = bitvf.io.Product
    bitvf.io.write_voids(bitvf.io.file_name(product.VOIDS, output, "box.dat", zspace=False), result)
    bitvf.io.write_overdense(bitvf.io.file_name(product.OVERDENSE, output, "box.dat", zspace=False), result)
    bitvf.io.write_galaxies(bitvf.io.file_name(product.GALAXIES, output, "box.dat", zspace=False), result)

    print("Written: output/voids_box.dat, output/overdense_basins_box.dat, output/galaxies_box.dat")


if __name__ == "__main__":
    main()
