# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# ===================================================================
# Example code: how to find the voids of a lightcone reconstruction
# ===================================================================

"""This example explains how to find voids in lightcone geometry from an
optimal transport reconstruction. otswap reconstructs the displacement field
of tracers observed on a part of the sky, in real space; bitvf then computes
the divergence of that field on a grid of shells in comoving distance,
declination rings and right-ascension sectors, smooths it, and splits it by a
watershed into basins around its local minima. Basins whose minimum is below
0, where the reconstructed flow diverges, are the voids; the others are the
overdense basins.

The two steps are separate calls, so the reconstruction can be kept, written,
or searched several times with different parameters. The run_lightcone
pipeline of the full_recons_lightcone example does both in one call.

The survey is described by the randoms, which must carry its footprint and
redshift distribution, and by the mask, a HEALPix map of the observed sky. A
cell of the grid whose centre falls on an unobserved pixel is masked: it
holds no divergence and belongs to no basin.

Usage: python lightcone.py [data_dir]. The catalogues are read from data_dir,
by default the examples' data folder; the output is written to the folder
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

    # ---------------------------------------------------------------
    # ------------ Read the tracer and random catalogues ------------
    # ---------------------------------------------------------------

    # Right ascension and declination in degrees, and redshift, shape (N, 3).
    # They are passed as they are, with angle_unit="deg": the Python
    # interface converts at its boundary
    tracers_sky = otswap.io.read(data / "lightcone_tracers.dat", [0, 1, 2]).values
    randoms_sky = otswap.io.read(data / "lightcone_randoms.dat", [0, 1, 2]).values

    print(f"Tracers read: {len(tracers_sky)}")
    print(f"Randoms read: {len(randoms_sky)}, {len(randoms_sky) / len(tracers_sky):g} per tracer")

    # ---------------------------------------
    # ------------ Read the mask ------------
    # ---------------------------------------

    # A HEALPix map, RING or NESTED, from a FITS file, here gzip-compressed:
    # a pixel is observed when its value is greater than 0. Tracers on
    # unobserved pixels keep their row in every per-tracer output, flagged,
    # and take no part
    mask = otswap.Mask(str(data / "lightcone_mask.fits.gz"))

    print(f"Mask: NSIDE {mask.nside}, {mask.sky_area_deg2:g} deg^2")

    # ----------------------------------------------------
    # ------------ Set the cosmological model ------------
    # ----------------------------------------------------

    # The comoving distances, from Omega_m, h, w0 and wa of a flat
    # cosmology, tabulated by default from z = 0 to 10. The table must cover
    # the redshifts of tracers and randoms
    distances = otswap.DistanceTable.flat(0.286, 0.7, w0=-1.0, wa=0.0)

    # ------------------------------------------------
    # ------------ Run the reconstruction ------------
    # ------------------------------------------------

    # 8 realizations, a convergence threshold of 1e-2 and a fixed seed, so
    # that the run can be repeated; displacements that cross an unobserved
    # pixel of the mask are rejected. Each realization uses as many randoms
    # as there are tracers, so at least 8 per tracer are needed: 10 per
    # tracer is a good start
    ot_config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, reject_crossings=True,
                              max_unobserved_pixels_crossed=0, verbosity="normal")

    # The mean particle separation (mps) is measured in 30 redshift bins;
    # each needs enough tracers, so use fewer bins for a small sample. The
    # redshift cut spans the whole catalogue: narrow it to search a slice
    reconstruction = otswap.reconstruct_lightcone(tracers_sky, randoms_sky, mask=mask, n_bins=30,
                                                  distances=distances, angle_unit="deg",
                                                  redshift_cut=otswap.RedshiftCut(0.885, 1.10),
                                                  config=ot_config)

    # ---------------------------------------------
    # ------------ Find the voids -----------------
    # ---------------------------------------------

    # The search reads everything it needs from the reconstruction: the
    # tracers, their displacements in each realization, the mps(z) profile,
    # the distance table and the redshift cut. The randoms are given again
    # because the grid and the mean density of each cell are built on them.
    #
    # The two parameters that change the result:
    # - the smoothing: a Gaussian of width sigma, here 1 mps, the default and
    #   a reasonable starting value; a larger sigma merges small basins into
    #   larger ones;
    # - the cell size, left to its default: 1 mps up to 20 realizations,
    #   (20 / n_realizations)^(1/3) mps above, measured at the centre of each
    #   shell.
    #
    # keep_assignment keeps, for each tracer, the basin it ended up in, which
    # the galaxies file records. "detailed" prints one line per step
    config = bitvf.Config(smoothing=bitvf.Smoothing(sigma=1.0, unit=bitvf.LengthUnit.MPS),
                          keep_assignment=True, verbosity="detailed")
    result = bitvf.find_voids_lightcone(reconstruction, randoms_sky, mask=mask, angle_unit="deg", config=config)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # Each catalogue has one row per basin, in the order of its seed cell.
    # radius is the radius of the sphere with the volume of the void's core,
    # the part where the divergence is below 0, in Mpc/h; edge_flag is 0 for
    # a void that does not touch the edges of the grid, where the flow is not
    # known
    voids = result.voids
    inside = int((voids.edge_flag == 0).sum())

    print(f"Tracers on unobserved pixels: {result.counts.tracers_outside_mask}")
    print(f"Voids: {voids.n_basins}, {inside} not touching the edges of the grid")
    print(f"Largest void radius: {voids.radius.max():g} Mpc/h")
    print(f"Overdense basins: {result.overdense.n_basins}")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # One file per product, ASCII from the extension: the voids and the
    # overdense basins, one row per basin, and where each tracer ended up,
    # one row per tracer. Angles are written in degrees. io.file_name gives
    # the conventional names, voids_lightcone.dat,
    # overdense_basins_lightcone.dat and galaxies_lightcone.dat. The writers
    # are the C++ library's, so the files are the C++ example's, byte for
    # byte
    product = bitvf.io.Product
    bitvf.io.write_voids(bitvf.io.file_name(product.VOIDS, output, "lightcone.dat", zspace=False), result)
    bitvf.io.write_overdense(bitvf.io.file_name(product.OVERDENSE, output, "lightcone.dat", zspace=False), result)
    bitvf.io.write_galaxies(bitvf.io.file_name(product.GALAXIES, output, "lightcone.dat", zspace=False), result)

    print("Written: output/voids_lightcone.dat, output/overdense_basins_lightcone.dat, "
          "output/galaxies_lightcone.dat")


if __name__ == "__main__":
    main()
