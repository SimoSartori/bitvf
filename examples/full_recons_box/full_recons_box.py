# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

# =====================================================================
# Example code: how to find real-space voids from a redshift-space box
# =====================================================================

"""This example explains how to run the whole BitVF pipeline on a box
catalogue observed in redshift space, in one call. The tracers are distorted
along one Cartesian axis, the line of sight. The pipeline:

1. reconstructs the displacement field of the redshift-space tracers by
   optimal transport (otswap);
2. moves each tracer along the line of sight to the position it would have in
   real space, from the reconstructed displacements and linear theory
   (otswap's correction of redshift-space distortions);
3. reconstructs the corrected tracers again, with the same randoms and seed;
4. finds the voids of that real-space reconstruction, on the grid of the
   first one.

On request it also finds the voids of the first, redshift-space,
reconstruction, which this example compares with the real-space ones.

Tracers the correction leaves unmoved, or moves outside the grid, keep their
row in every per-tracer output, flagged, and take no part in the second
reconstruction.

Usage: python full_recons_box.py [data_dir]. The catalogue is read from
data_dir, by default the examples' data folder; the output is written to the
folder output/ next to this file, created if missing. The field example reads
the divergence file it writes.
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

    # Comoving Cartesian coordinates in Mpc/h, shape (N, 3), in redshift
    # space with the line of sight along y
    tracers = otswap.io.read(data / "box_halos_redshift_space.dat", [0, 1, 2]).values

    print(f"Tracers: {len(tracers)}")

    # --------------------------------------------------
    # ------------ Configure the pipeline --------------
    # --------------------------------------------------

    # The reconstruction, for both passes: 8 realizations, a convergence
    # threshold of 1e-2 and a fixed seed. No randoms are given, so bitvf
    # draws 8 per halo, uniformly in the box, from a seed derived from this
    # one (Result.random_seed)
    ot_config = otswap.Config(n_realizations=8, convergence=1e-2, seed=12345)

    # The halos fill a cube of side 400 Mpc/h, [269, 669) on each axis
    box = ((269.0, 269.0, 269.0), (669.0, 669.0, 669.0))

    # The correction: the line of sight (axis 0, 1 or 2 for x, y or z), the
    # box's redshift, at which the growth rate is read from the distance
    # table below, and the tracers' linear bias. sigma is the width of the
    # Gaussian average of the displacements around each tracer, in Mpc/h:
    # 10 Mpc/h is a reasonable starting value, but the best value depends on
    # the sample and should be checked in each analysis. zspace_voids also
    # searches the redshift-space pass
    rsd = bitvf.BoxRsd(axis=1, redshift=1.0, bias=1.91255, sigma=10.0, zspace_voids=True)

    # The void search, on both passes: the default smoothing of 1 mps and the
    # default cell size. The divergence grid is kept, with the field before
    # and after smoothing, so that it can be written and searched again later
    # without the reconstruction (see the field example)
    voids = bitvf.Config(keep_divergence=True, verbosity="detailed")

    config = bitvf.BoxConfig(ot=ot_config, limits=box, rsd=rsd, voids=voids)

    # The correction needs only the linear growth rate f at the box's
    # redshift. The table computes it from Omega_m, h, w0 and wa of a flat
    # cosmology
    distances = otswap.DistanceTable.flat(0.3186, 0.67, w0=-1.0, wa=0.0)

    # ----------------------------------------------
    # ------------ Run the pipeline ----------------
    # ----------------------------------------------

    run = bitvf.run_box(tracers, config=config, distances=distances)

    # ----------------------------------------------
    # ------------ Summarize the result ------------
    # ----------------------------------------------

    # run.voids is the real-space search, whose per-halo outputs keep one
    # row per input halo, and run.zspace_voids the redshift-space one. Redshift-space distortions
    # stretch voids along the line of sight and change their number
    real, zspace = run.voids, run.zspace_voids

    print(f"Randoms drawn with seed {real.random_seed}")
    print(f"Halos left out of the second pass: {real.counts.tracers_uncorrected} uncorrected, "
          f"{real.counts.tracers_outside_grid} moved outside the grid")
    print(f"Voids in redshift space: {zspace.voids.n_basins}")
    print(f"Voids in real space: {real.voids.n_basins}")

    # ------------------------------------------
    # ------------ Write the output ------------
    # ------------------------------------------

    # The voids and the divergence grid of each pass. The redshift-space
    # files have zspace_ in their names: voids_zspace_box.dat and
    # divergence_zspace_box.dat, then voids_box.dat and divergence_box.dat.
    # The writers are the C++ library's, so the files are the C++ example's,
    # byte for byte
    product = bitvf.io.Product
    bitvf.io.write_voids(bitvf.io.file_name(product.VOIDS, output, "box.dat", zspace=True), zspace)
    bitvf.io.write_divergence(bitvf.io.file_name(product.DIVERGENCE, output, "box.dat", zspace=True), zspace)
    bitvf.io.write_voids(bitvf.io.file_name(product.VOIDS, output, "box.dat", zspace=False), real)
    bitvf.io.write_divergence(bitvf.io.file_name(product.DIVERGENCE, output, "box.dat", zspace=False), real)

    print("Written: output/voids_zspace_box.dat, output/divergence_zspace_box.dat, output/voids_box.dat, "
          "output/divergence_box.dat")


if __name__ == "__main__":
    main()
