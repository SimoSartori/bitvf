// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// =====================================================================
// Example code: how to find real-space voids from a redshift-space box
// =====================================================================

/*
  This example explains how to run the whole BitVF pipeline on a box
  catalogue observed in redshift space, in one call. The tracers are
  distorted along one Cartesian axis, the line of sight. The pipeline:

  1. reconstructs the displacement field of the redshift-space tracers
     by optimal transport (otswap);
  2. moves each tracer along the line of sight to the position it would
     have in real space, from the reconstructed displacements and linear
     theory (otswap's correction of redshift-space distortions);
  3. reconstructs the corrected tracers again, with the same randoms and
     seed;
  4. finds the voids of that real-space reconstruction, on the grid of
     the first one.

  On request it also finds the voids of the first, redshift-space,
  reconstruction, which this example compares with the real-space ones.

  Tracers the correction leaves unmoved, or moves outside the grid, keep
  their row in every per-tracer output, flagged, and take no part in the
  second reconstruction.

  Usage: full_recons_box [data_dir]. The catalogue is read from
  data_dir, by default the examples' data folder; the output is written
  to the folder output/ next to this file, created if missing. The field
  example reads the divergence file it writes.
*/

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "bitvf/BitVF.h"
#include "bitvf/io.h"


// ==========================================================================


int main (int argc, char** argv)
{
  try {

    // -------------------------------------------------------------
    // ------------ Locate the input and output folders ------------
    // -------------------------------------------------------------

    // BITVF_EXAMPLES_DATA and BITVF_EXAMPLE_OUTPUT are set by this
    // example's CMakeLists.txt to the data folder and to output/ next to
    // this file. A first command-line argument replaces the data folder
    const std::string data = (argc > 1) ? argv[1] : BITVF_EXAMPLES_DATA;
    const std::string output = BITVF_EXAMPLE_OUTPUT;
    std::filesystem::create_directories(output);


    // ----------------------------------------------------
    // ------------ Read the tracer catalogue -------------
    // ----------------------------------------------------

    // Comoving Cartesian coordinates in Mpc/h, x, y, z for each halo, one
    // row after the other, in redshift space with the line of sight along y
    const std::vector<double> tracers = otswap::io::read(data + "/box_halos_redshift_space.dat",
                                                         {"0", "1", "2"}).values;

    std::cout << "Tracers: " << tracers.size() / 3 << std::endl;


    // --------------------------------------------------
    // ------------ Configure the pipeline --------------
    // --------------------------------------------------

    // The reconstruction, for both passes: 8 realizations, a convergence
    // threshold of 1e-2 and a fixed seed. No randoms are given, so bitvf
    // draws 8 per halo, uniformly in the box, from a seed derived from
    // this one (Result::randomSeed)
    bitvf::BoxConfig config;
    config.ot.nRealizations = 8;
    config.ot.convergence = 1.e-2;
    config.ot.seed = 12345;
    config.ot.log = &std::cout;

    // The halos fill a cube of side 400 Mpc/h, [269, 669) on each axis
    config.limits = otswap::BoxLimits{{269., 269., 269.}, {669., 669., 669.}};

    // The correction: the line of sight (axis 0, 1 or 2 for x, y or z),
    // the box's redshift, at which the growth rate is read from the
    // distance table below, and the tracers' linear bias. sigma is the
    // width of the Gaussian average of the displacements around each
    // tracer, in Mpc/h: 10 Mpc/h is a reasonable starting value, but the
    // best value depends on the sample and should be checked in each
    // analysis. zspaceVoids also searches the redshift-space pass
    bitvf::BoxRsd rsd;
    rsd.axis = 1;
    rsd.redshift = 1.;
    rsd.bias = 1.91255;
    rsd.sigma = 10.;
    rsd.correction.log = &std::cout;
    rsd.zspaceVoids = true;
    config.rsd = rsd;

    // The void search, on both passes: the default smoothing of 1 mps and
    // the default cell size. The divergence grid is kept, with the field
    // before and after smoothing, so that it can be written and searched
    // again later without the reconstruction (see the field example)
    config.voids.keepDivergence = true;
    config.voids.verbosity = otswap::Verbosity::Detailed;
    config.voids.log = &std::cout;

    // The correction needs only the linear growth rate f at the box's
    // redshift. The table computes it from Omega_m, h, w0 and wa of a flat
    // cosmology
    const otswap::DistanceTable distances(0.3186, 0.67, -1., 0.);


    // ----------------------------------------------
    // ------------ Run the pipeline ----------------
    // ----------------------------------------------

    const bitvf::BitVFResult run = bitvf::runBox(tracers, {}, distances, config);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // run.voids is the real-space search, whose per-halo outputs keep
    // one row per input halo, and run.zspaceVoids the redshift-space one. Redshift-space distortions
    // stretch voids along the line of sight and change their number
    const bitvf::Result& real = run.voids;
    const bitvf::Result& zspace = *run.zspaceVoids;

    std::cout << "Randoms drawn with seed " << *real.randomSeed << std::endl;
    std::cout << "Halos left out of the second pass: " << real.counts.tracersUncorrected << " uncorrected, "
              << real.counts.tracersOutsideGrid << " moved outside the grid" << std::endl;
    std::cout << "Voids in redshift space: " << zspace.voids.nBasins << std::endl;
    std::cout << "Voids in real space: " << real.voids.nBasins << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // The voids and the divergence grid of each pass. The redshift-space
    // files have zspace_ in their names: voids_zspace_box.dat and
    // divergence_zspace_box.dat, then voids_box.dat and divergence_box.dat
    using bitvf::io::Product;
    bitvf::io::writeVoids(bitvf::io::fileName(Product::Voids, output, "box.dat", true), zspace);
    bitvf::io::writeDivergence(bitvf::io::fileName(Product::Divergence, output, "box.dat", true), zspace);
    bitvf::io::writeVoids(bitvf::io::fileName(Product::Voids, output, "box.dat", false), real);
    bitvf::io::writeDivergence(bitvf::io::fileName(Product::Divergence, output, "box.dat", false), real);

    std::cout << "Written: output/voids_zspace_box.dat, output/divergence_zspace_box.dat, output/voids_box.dat, "
                 "output/divergence_box.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }
  catch (const bitvf::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
