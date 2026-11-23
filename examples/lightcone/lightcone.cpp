// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// ===================================================================
// Example code: how to find the voids of a lightcone reconstruction
// ===================================================================

/*
  This example explains how to find voids in lightcone geometry from an
  optimal transport reconstruction. otswap reconstructs the displacement
  field of tracers observed on a part of the sky, in real space; bitvf
  then computes the divergence of that field on a grid of shells in
  comoving distance, declination rings and right-ascension sectors,
  smooths it, and splits it by a watershed into basins around its local
  minima. Basins whose minimum is below 0, where the reconstructed flow
  diverges, are the voids; the others are the overdense basins.

  The two steps are separate calls, so the reconstruction can be kept,
  written, or searched several times with different parameters. The
  runLightcone pipeline of the full_recons_lightcone example does both
  in one call.

  The survey is described by the randoms, which must carry its footprint
  and redshift distribution, and by the mask, a HEALPix map of the
  observed sky. A cell of the grid whose centre falls on an unobserved
  pixel is masked: it holds no divergence and belongs to no basin.

  Usage: lightcone [data_dir]. The catalogues are read from data_dir, by
  default the examples' data folder; the output is written to the folder
  output/ next to this file, created if missing.
*/

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "bitvf/VoidFinder.h"
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

    const std::size_t nTracers = tracersSky.size() / 3;
    const std::size_t nRandoms = randomsSky.size() / 3;

    std::cout << "Tracers read: " << nTracers << std::endl;
    std::cout << "Randoms read: " << nRandoms << ", " << (double)nRandoms / (double)nTracers << " per tracer"
              << std::endl;


    // ---------------------------------------
    // ------------ Read the mask ------------
    // ---------------------------------------

    // A HEALPix map, RING or NESTED, from a FITS file, here gzip-compressed:
    // a pixel is observed when its value is greater than 0. Tracers on
    // unobserved pixels keep their row in every per-tracer output, flagged,
    // and take no part
    const otswap::Mask mask(data + "/lightcone_mask.fits.gz");

    std::cout << "Mask: NSIDE " << mask.nside() << ", " << mask.skyAreaDeg2() << " deg^2" << std::endl;


    // ----------------------------------------------------
    // ------------ Set the cosmological model ------------
    // ----------------------------------------------------

    // The comoving distances, from Omega_m, h, w0 and wa of a flat
    // cosmology, tabulated by default from z = 0 to 10. The table must
    // cover the redshifts of tracers and randoms
    const otswap::DistanceTable distances(0.286, 0.7, -1., 0.);


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // 8 realizations, a convergence threshold of 1e-2 and a fixed seed, so
    // that the run can be repeated; displacements that cross an unobserved
    // pixel of the mask are rejected. Each realization uses as many randoms
    // as there are tracers, so at least 8 per tracer are needed: 30 per
    // tracer, the minimum with otswap's default of 30 realizations, is a
    // good start. The report goes to std::cout, in order with
    // this program's lines
    otswap::Config otConfig;
    otConfig.nRealizations = 8;
    otConfig.convergence = 1.e-2;
    otConfig.seed = 12345;
    otConfig.rejectCrossings = true;
    otConfig.maxUnobservedPixelsCrossed = 0;
    otConfig.verbosity = otswap::Verbosity::Normal;
    otConfig.log = &std::cout;

    // The mean particle separation (mps) is measured in 30 redshift bins;
    // each needs enough tracers, so use fewer bins for a small sample. The
    // redshift cut spans the whole catalogue: narrow it to search a slice
    const unsigned nBins = 30;
    const otswap::RedshiftCut cut{0.885, 1.10};

    const otswap::Result reconstruction =
      otswap::reconstructLightcone(tracersSky, randomsSky, mask, nBins, distances, otConfig, cut);


    // ---------------------------------------------
    // ------------ Find the voids -----------------
    // ---------------------------------------------

    // The search reads everything it needs from the reconstruction: the
    // tracers, their displacements in each realization, the mps(z) profile,
    // the distance table and the redshift cut. The randoms are given again
    // because the grid and the mean density of each cell are built on them.
    //
    // The two parameters that change the result:
    // - the smoothing: a Gaussian of width sigma, here 1 mps, the default
    //   and a reasonable starting value; a larger sigma merges small basins
    //   into larger ones;
    // - the cell size, left to its default: 1 mps up to 20 realizations,
    //   (20 / nRealizations)^(1/3) mps above, measured at the centre of
    //   each shell.
    //
    // keepAssignment keeps, for each tracer, the basin it ended up in, which
    // the galaxies file records. Detailed prints one line per step
    bitvf::Config config;
    config.smoothing.sigma = 1.;
    config.smoothing.unit = bitvf::LengthUnit::Mps;
    config.keepAssignment = true;
    config.verbosity = otswap::Verbosity::Detailed;
    config.log = &std::cout;

    const bitvf::Result result = bitvf::findVoidsLightcone(reconstruction, randomsSky, mask, config);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // Each catalogue has one row per basin, in the order of its seed cell.
    // radius is the radius of the sphere with the volume of the void's
    // core, the part where the divergence is below 0, in Mpc/h; edgeFlag is
    // 0 for a void that does not touch the edges of the grid, where the
    // flow is not known
    const bitvf::BasinCatalog& voids = result.voids;
    std::size_t inside = 0;
    double largest = 0.;
    for (std::size_t b = 0; b < voids.nBasins; ++b) {
      if (voids.edgeFlag[b] == 0) ++inside;
      largest = std::max(largest, voids.radius[b]);
    }

    std::cout << "Tracers on unobserved pixels: " << result.counts.tracersOutsideMask << std::endl;
    std::cout << "Voids: " << voids.nBasins << ", " << inside << " not touching the edges of the grid" << std::endl;
    std::cout << "Largest void radius: " << largest << " Mpc/h" << std::endl;
    std::cout << "Overdense basins: " << result.overdense.nBasins << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // One file per product, ASCII from the extension: the voids and the
    // overdense basins, one row per basin, and where each tracer ended up,
    // one row per tracer. Angles are written in degrees. io::fileName gives
    // the conventional names, voids_lightcone.dat,
    // overdense_basins_lightcone.dat and galaxies_lightcone.dat
    const std::string voidsFile = bitvf::io::fileName(bitvf::io::Product::Voids, output, "lightcone.dat", false);
    const std::string overdenseFile = bitvf::io::fileName(bitvf::io::Product::Overdense, output, "lightcone.dat", false);
    const std::string galaxiesFile = bitvf::io::fileName(bitvf::io::Product::Galaxies, output, "lightcone.dat", false);

    bitvf::io::writeVoids(voidsFile, result);
    bitvf::io::writeOverdense(overdenseFile, result);
    bitvf::io::writeGalaxies(galaxiesFile, result);

    std::cout << "Written: output/voids_lightcone.dat, output/overdense_basins_lightcone.dat, "
                 "output/galaxies_lightcone.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }
  catch (const bitvf::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
