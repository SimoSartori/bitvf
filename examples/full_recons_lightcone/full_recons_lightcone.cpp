// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// ===========================================================================
// Example code: how to find real-space voids from a redshift-space lightcone
// ===========================================================================

/*
  This example explains how to run the whole BitVF pipeline on a survey
  catalogue observed in redshift space, in lightcone geometry, in one
  call. The redshift distortions displace each tracer along its line of
  sight, the direction from the observer. The pipeline:

  1. reconstructs the displacement field of the redshift-space tracers
     by optimal transport (otswap);
  2. moves each tracer along its line of sight to the position it would
     have in real space, from the reconstructed displacements and linear
     theory (otswap's correction of redshift-space distortions);
  3. reconstructs the corrected tracers again, with randoms that follow
     the redshift distribution of the corrected tracers;
  4. finds the voids of that real-space reconstruction, on the grid of
     the first one.

  On request it also finds the voids of the first, redshift-space,
  reconstruction, which this example compares with the real-space ones.

  Tracers on unobserved pixels or outside the redshift cut keep their row
  in every per-tracer output, flagged, and take no part. So do those the
  correction leaves unmoved, or moves outside the grid or the cut, in the
  second reconstruction.

  Usage: full_recons_lightcone [data_dir]. The catalogues are read from
  data_dir, by default the examples' data folder; the output is written
  to the folder output/ next to this file, created if missing.
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


    // ---------------------------------------------------------------
    // ------------ Read the tracer and random catalogues ------------
    // ---------------------------------------------------------------

    // Right ascension, declination and redshift for each object, one row
    // after the other. The angles are read in degrees and converted to
    // radians by skyToRadians: the C++ library works in radians
    const std::vector<double> tracersSky =
      otswap::skyToRadians(otswap::io::read(data + "/lightcone_tracers.dat", {"0", "1", "2"}).values);
    const std::vector<double> randomsSky =
      otswap::skyToRadians(otswap::io::read(data + "/lightcone_randoms.dat", {"0", "1", "2"}).values);

    std::cout << "Tracers read: " << tracersSky.size() / 3 << std::endl;
    std::cout << "Randoms read: " << randomsSky.size() / 3 << std::endl;


    // ---------------------------------------
    // ------------ Read the mask ------------
    // ---------------------------------------

    // A HEALPix map from a FITS file, here gzip-compressed: a pixel is
    // observed when its value is greater than 0
    const otswap::Mask mask(data + "/lightcone_mask.fits.gz");

    std::cout << "Mask: NSIDE " << mask.nside() << ", " << mask.skyAreaDeg2() << " deg^2" << std::endl;


    // --------------------------------------------------
    // ------------ Configure the pipeline --------------
    // --------------------------------------------------

    // The reconstruction, for both passes: 8 realizations, a convergence
    // threshold of 1e-2 and a fixed seed; displacements that cross an
    // unobserved pixel of the mask are rejected. Each realization uses as
    // many randoms as there are tracers, so at least 8 per tracer are
    // needed: 10 per tracer is a good start
    bitvf::LightconeConfig config;
    config.ot.nRealizations = 8;
    config.ot.convergence = 1.e-2;
    config.ot.seed = 12345;
    config.ot.rejectCrossings = true;
    config.ot.maxUnobservedPixelsCrossed = 0;
    config.ot.log = &std::cout;

    // The mean particle separation (mps) is measured in 30 redshift bins on
    // each pass; each needs enough tracers, so use fewer bins for a small
    // sample. The redshift cut spans the whole catalogue: narrow it to
    // search a slice
    config.nBins = 30;
    config.cut = otswap::RedshiftCut{0.885, 1.10};

    // The correction: the tracers' linear bias as nodes in redshift,
    // interpolated linearly between them, and sigma, the width of the
    // Gaussian average of the displacements around each tracer, in Mpc/h:
    // 10 Mpc/h is a reasonable starting value, but the best value depends
    // on the sample and should be checked in each analysis. Resample, the
    // default, gives the randoms of the second reconstruction the redshift
    // distribution of the corrected tracers; None would keep the observed
    // one. zspaceVoids also searches the redshift-space pass
    bitvf::LightconeRsd rsd;
    rsd.bias = otswap::io::readBiasTable(data + "/lightcone_bias.dat");
    rsd.sigma = 10.;
    rsd.correction.log = &std::cout;
    rsd.secondPassRandoms = bitvf::SecondPassRandoms::Resample;
    rsd.zspaceVoids = true;
    config.rsd = rsd;

    // The void search, on both passes: the default smoothing of 1 mps and
    // the default cell size. The divergence grid is kept, with the field
    // before and after smoothing, so that it can be written and searched
    // again later without the reconstruction
    config.voids.keepDivergence = true;
    config.voids.verbosity = otswap::Verbosity::Detailed;
    config.voids.log = &std::cout;

    // The comoving distances and the linear growth rate f(z) of the
    // correction, from Omega_m, h, w0 and wa of a flat cosmology. The table
    // must cover the corrected distances too, which can fall slightly
    // outside the catalogue's redshift range: the default range, z = 0 to
    // 10, leaves room for that
    const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);


    // ----------------------------------------------
    // ------------ Run the pipeline ----------------
    // ----------------------------------------------

    const bitvf::BitVFResult run = bitvf::runLightcone(tracersSky, randomsSky, mask, distances, config);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // run.voids is the real-space search, whose per-tracer outputs keep
    // one row per input tracer, and run.zspaceVoids the redshift-space one.
    // Redshift-space distortions stretch voids along the line of sight and
    // change their number
    const bitvf::Result& real = run.voids;
    const bitvf::Result& zspace = *run.zspaceVoids;

    std::cout << "Tracers left out of the second pass: " << real.counts.tracersUncorrected << " uncorrected, "
              << real.counts.tracersOutsideGrid << " moved outside the grid, "
              << real.counts.tracersOutsideCut << " moved outside the redshift cut" << std::endl;
    std::cout << "Voids in redshift space: " << zspace.voids.nBasins << std::endl;
    std::cout << "Voids in real space: " << real.voids.nBasins << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // The voids and the divergence grid of each pass, with the mps(z)
    // profile of the grid in a twin file (.mps.dat). The redshift-space
    // files have zspace_ in their names: voids_zspace_lightcone.dat and
    // divergence_zspace_lightcone.dat, then voids_lightcone.dat and
    // divergence_lightcone.dat
    using bitvf::io::Product;
    bitvf::io::writeVoids(bitvf::io::fileName(Product::Voids, output, "lightcone.dat", true), zspace);
    bitvf::io::writeDivergence(bitvf::io::fileName(Product::Divergence, output, "lightcone.dat", true), zspace);
    bitvf::io::writeVoids(bitvf::io::fileName(Product::Voids, output, "lightcone.dat", false), real);
    bitvf::io::writeDivergence(bitvf::io::fileName(Product::Divergence, output, "lightcone.dat", false), real);

    std::cout << "Written: output/voids_zspace_lightcone.dat, output/divergence_zspace_lightcone.dat, "
                 "output/voids_lightcone.dat, output/divergence_lightcone.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }
  catch (const bitvf::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
