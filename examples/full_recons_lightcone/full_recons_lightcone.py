# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# ===========================================================================
# Example code: how to find real-space voids from a redshift-space lightcone
# ===========================================================================

"""This example explains how to run the whole BitVF pipeline on a survey
catalogue observed in redshift space, in lightcone geometry, in one call. The
redshift distortions displace each tracer along its line of sight, the
direction from the observer. The pipeline:

1. reconstructs the displacement field of the redshift-space tracers by
   optimal transport (otswap);
2. moves each tracer along its line of sight to the position it would have in
   real space, from the reconstructed displacements and linear theory
   (otswap's correction of redshift-space distortions);
3. reconstructs the corrected tracers again, with randoms that follow the
   redshift distribution of the corrected tracers;
4. finds the voids of that real-space reconstruction, on the grid of the
   first one.

On request it also finds the voids of the first, redshift-space,
reconstruction, which this example compares with the real-space ones.

Tracers on unobserved pixels or outside the redshift cut keep their row in
every per-tracer output, flagged, and take no part. So do those the
correction leaves unmoved, or moves outside the grid or the cut, in the
second reconstruction.

Usage: python full_recons_lightcone.py [data_dir]. The catalogues are read
from data_dir, by default the examples' data folder; the output is written to
the folder output/ next to this file, created if missing.
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

    # Right ascension and declination in degrees, and redshift, shape (N, 3);
    # passed as they are, with angle_unit="deg"
    tracers_sky = otswap.io.read(data / "lightcone_tracers.dat", [0, 1, 2]).values
    randoms_sky = otswap.io.read(data / "lightcone_randoms.dat", [0, 1, 2]).values

    print(f"Tracers read: {len(tracers_sky)}")
    print(f"Randoms read: {len(randoms_sky)}")

    # ---------------------------------------
    # ------------ Read the mask ------------
    # ---------------------------------------

    # A HEALPix map from a FITS file, here gzip-compressed: a pixel is
    # observed when its value is greater than 0
    mask = otswap.Mask(str(data / "lightcone_mask.fits.gz"))

    print(f"Mask: NSIDE {mask.nside}, {mask.sky_area_deg2:g} deg^2")

    # --------------------------------------------------
    # ------------ Configure the pipeline --------------
    # --------------------------------------------------

    # The reconstruction, for both passes: 8 realizations, a convergence
    # threshold of 1e-2 and a fixed seed; displacements that cross an
    # unobserved pixel of the mask are rejected. Each realization uses as
    # many randoms as there are tracers, so at least 8 per tracer are needed:
    # 30 per tracer, the minimum with otswap's default of 30 realizations, is
    # a good start
    ot_config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345, reject_crossings=True,
                              max_unobserved_pixels_crossed=0)

    # The correction: the tracers' linear bias as nodes in redshift,
    # interpolated linearly between them, and sigma, the width of the
    # Gaussian average of the displacements around each tracer, in Mpc/h:
    # 10 Mpc/h is a reasonable starting value, but the best value depends on
    # the sample and should be checked in each analysis. RESAMPLE, the
    # default, gives the randoms of the second reconstruction the redshift
    # distribution of the corrected tracers; NONE would keep the observed
    # one. zspace_voids also searches the redshift-space pass
    rsd = bitvf.LightconeRsd(bias=otswap.io.read_bias_table(data / "lightcone_bias.dat"), sigma=10.0,
                             second_pass_randoms=bitvf.SecondPassRandoms.RESAMPLE, zspace_voids=True)

    # The void search, on both passes: the default smoothing of 1 mps and the
    # default cell size. The divergence grid is kept, with the field before
    # and after smoothing, so that it can be written and searched again later
    # without the reconstruction
    voids = bitvf.Config(keep_divergence=True, verbosity="detailed")

    # The mean particle separation (mps) is measured in 30 redshift bins on
    # each pass; each needs enough tracers, so use fewer bins for a small
    # sample. The redshift cut spans the whole catalogue: narrow it to search
    # a slice
    config = bitvf.LightconeConfig(ot=ot_config, n_bins=30, cut=otswap.RedshiftCut(0.885, 1.10), rsd=rsd,
                                   voids=voids)

    # The comoving distances and the linear growth rate f(z) of the
    # correction, from Omega_m, h, w0 and wa of a flat cosmology. The table
    # must cover the corrected distances too, which can fall slightly outside
    # the catalogue's redshift range: the default range, z = 0 to 10, leaves
    # room for that
    distances = otswap.DistanceTable.flat(0.286, 0.7, w0=-1.0, wa=0.0)

    # ----------------------------------------------
    # ------------ Run the pipeline ----------------
    # ----------------------------------------------

    run = bitvf.run_lightcone(tracers_sky, randoms_sky, mask=mask, distances=distances, angle_unit="deg",
                              config=config)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # run.voids is the real-space search, whose per-tracer outputs keep one
    # row per input tracer, and run.zspace_voids the redshift-space one.
    # Redshift-space distortions stretch voids along the line of sight and
    # change their number
    real, zspace = run.voids, run.zspace_voids

    print(f"Tracers left out of the second pass: {real.counts.tracers_uncorrected} uncorrected, "
          f"{real.counts.tracers_outside_grid} moved outside the grid, "
          f"{real.counts.tracers_outside_cut} moved outside the redshift cut")
    print(f"Voids in redshift space: {zspace.voids.n_basins}")
    print(f"Voids in real space: {real.voids.n_basins}")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # The voids and the divergence grid of each pass, with the mps(z) profile
    # of the grid in a twin file (.mps.dat). The redshift-space files have
    # zspace_ in their names: voids_zspace_lightcone.dat and
    # divergence_zspace_lightcone.dat, then voids_lightcone.dat and
    # divergence_lightcone.dat. The writers are the C++ library's, so the
    # files are the C++ example's, byte for byte
    product = bitvf.io.Product
    bitvf.io.write_voids(bitvf.io.file_name(product.VOIDS, output, "lightcone.dat", zspace=True), zspace)
    bitvf.io.write_divergence(bitvf.io.file_name(product.DIVERGENCE, output, "lightcone.dat", zspace=True), zspace)
    bitvf.io.write_voids(bitvf.io.file_name(product.VOIDS, output, "lightcone.dat", zspace=False), real)
    bitvf.io.write_divergence(bitvf.io.file_name(product.DIVERGENCE, output, "lightcone.dat", zspace=False), real)

    print("Written: output/voids_zspace_lightcone.dat, output/divergence_zspace_lightcone.dat, "
          "output/voids_lightcone.dat, output/divergence_lightcone.dat")


if __name__ == "__main__":
    main()
