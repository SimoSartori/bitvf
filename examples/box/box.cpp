// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// ===============================================================
// Example code: how to find the voids of a box reconstruction
// ===============================================================

/*
  This example explains how to find voids in box geometry from an
  optimal transport reconstruction. otswap reconstructs the displacement
  field of halos in real space; bitvf then computes the divergence of
  that field on a grid of cubic cells, smooths it, and splits it by a
  watershed into basins around its local minima. Basins whose minimum
  is below 0, where the reconstructed flow diverges, are the voids; the
  others are the overdense basins.

  The two steps are separate calls, so the reconstruction can be kept,
  written, or searched several times with different parameters. The
  runBox pipeline of the full_recons_box example does both in one call.

  The box is not periodic: nothing flows through its faces, and the
  basins that touch them are flagged (edgeFlag), so that they can be
  left out of an analysis.

  Usage: box [data_dir]. The catalogue is read from data_dir, by default
  the examples' data folder; the output is written to the folder output/
  next to this file, created if missing.
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


    // ----------------------------------------------------
    // ------------ Read the tracer catalogue -------------
    // ----------------------------------------------------

    // Comoving Cartesian coordinates in Mpc/h, x, y, z for each halo, one
    // row after the other, in real space
    const std::vector<double> tracers = otswap::io::read(data + "/box_halos_real_space.dat", {"0", "1", "2"}).values;

    std::cout << "Tracers: " << tracers.size() / 3 << std::endl;


    // ------------------------------------------------
    // ------------ Run the reconstruction ------------
    // ------------------------------------------------

    // The halos fill a cube of side 400 Mpc/h, [269, 669) on each axis.
    // Given as box limits, the cube is where otswap draws the randoms and
    // where bitvf builds its grid; the mean particle separation follows
    // from its volume and the halos inside it
    const otswap::BoxLimits box{{269., 269., 269.}, {669., 669., 669.}};

    // 8 realizations, a convergence threshold of 1e-2 and a fixed seed, so
    // that the run can be repeated. Without randoms, otswap draws 8 per
    // halo, uniformly in the box. The report goes to std::cout, in order
    // with this program's lines
    otswap::Config otConfig;
    otConfig.nRealizations = 8;
    otConfig.convergence = 1.e-2;
    otConfig.seed = 12345;
    otConfig.verbosity = otswap::Verbosity::Normal;
    otConfig.log = &std::cout;

    const otswap::Result reconstruction = otswap::reconstructBox(tracers, {}, otConfig, box);


    // ---------------------------------------------
    // ------------ Find the voids -----------------
    // ---------------------------------------------

    // The search reads everything it needs from the reconstruction: the
    // tracers, their displacements in each realization, the mean particle
    // separation (mps) and the box.
    //
    // The two parameters that change the result:
    // - the smoothing: a Gaussian of width sigma, here 1 mps, the default
    //   and a reasonable starting value; a larger sigma merges small basins
    //   into larger ones;
    // - the cell side, left to its default: 1 mps up to 20 realizations,
    //   (20 / nRealizations)^(1/3) mps above.
    //
    // keepAssignment keeps, for each halo, the basin it ended up in, which
    // the galaxies file records. Detailed prints one line per step
    bitvf::Config config;
    config.smoothing.sigma = 1.;
    config.smoothing.unit = bitvf::LengthUnit::Mps;
    config.keepAssignment = true;
    config.verbosity = otswap::Verbosity::Detailed;
    config.log = &std::cout;

    const bitvf::Result result = bitvf::findVoidsBox(reconstruction, config);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // Each catalogue has one row per basin, in the order of its seed cell.
    // radius is the radius of the sphere with the volume of the void's
    // core, the part where the divergence is below 0; edgeFlag is 0 for a
    // void that does not touch the faces of the grid
    const bitvf::BasinCatalog& voids = result.voids;
    std::size_t inside = 0;
    double largest = 0.;
    for (std::size_t b = 0; b < voids.nBasins; ++b) {
      if (voids.edgeFlag[b] == 0) ++inside;
      largest = std::max(largest, voids.radius[b]);
    }

    std::cout << "Mean particle separation: " << result.mps << " Mpc/h" << std::endl;
    std::cout << "Voids: " << voids.nBasins << ", " << inside << " not touching the faces of the grid" << std::endl;
    std::cout << "Largest void radius: " << largest << " Mpc/h" << std::endl;
    std::cout << "Overdense basins: " << result.overdense.nBasins << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // One file per product, ASCII from the extension: the voids and the
    // overdense basins, one row per basin, and where each halo ended up,
    // one row per halo. io::fileName gives the conventional names,
    // voids_box.dat, overdense_basins_box.dat and galaxies_box.dat
    const std::string voidsFile = bitvf::io::fileName(bitvf::io::Product::Voids, output, "box.dat", false);
    const std::string overdenseFile = bitvf::io::fileName(bitvf::io::Product::Overdense, output, "box.dat", false);
    const std::string galaxiesFile = bitvf::io::fileName(bitvf::io::Product::Galaxies, output, "box.dat", false);

    bitvf::io::writeVoids(voidsFile, result);
    bitvf::io::writeOverdense(overdenseFile, result);
    bitvf::io::writeGalaxies(galaxiesFile, result);

    std::cout << "Written: output/voids_box.dat, output/overdense_basins_box.dat, output/galaxies_box.dat"
              << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }
  catch (const bitvf::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
