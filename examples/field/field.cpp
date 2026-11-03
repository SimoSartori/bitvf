// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it

// =================================================================
// Example code: how to find the voids of a field on a box lattice
// =================================================================

/*
  This example explains how to find voids in a field given at the
  centres of the cells of a box lattice, with no tracers and no
  reconstruction: a density contrast, a potential, or any field whose
  voids lie below 0. A field with voids above 0 must be flipped first.

  The field is read from a table with one row per cell, in any order;
  the lattice, its cell side and its extent are inferred from the
  centres. The field is smoothed, split by the same watershed as the
  void finder, and catalogued: centres, radii, contact with the faces
  of the grid, and the field at each centre. Densities, which need
  tracers, are not computed.

  The input here is the divergence file written by the full_recons_box
  example: a divergence file of bitvf is a valid field, with the cell
  centres, the divergence before smoothing and the tracer weight of
  each cell. Searched with the smoothing of the search that wrote it,
  it gives back that search's basins; this example smooths it with a
  Gaussian of 15 Mpc/h instead, wider than the 1 mps of that search, so
  that small voids merge into larger ones.

  Usage: field [field_file]. The field is read from field_file, by
  default the divergence file of the full_recons_box example, which must
  then be run first; the output is written to the folder output/ next
  to this file, created if missing.
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

    // BITVF_EXAMPLE_FIELD and BITVF_EXAMPLE_OUTPUT are set by this
    // example's CMakeLists.txt to the divergence file of the
    // full_recons_box example and to output/ next to this file. A first
    // command-line argument replaces the input file
    const std::string input = (argc > 1) ? argv[1] : BITVF_EXAMPLE_FIELD;
    const std::string output = BITVF_EXAMPLE_OUTPUT;
    std::filesystem::create_directories(output);


    // -------------------------------------
    // ------------ Read the field ---------
    // -------------------------------------

    // The columns of the cell centres, in Mpc/h, of the field and of the
    // weight of each cell in the smoothing. Without a weight column every
    // cell weighs the same. The field is the divergence column, before
    // smoothing, not divergenceSmoothed: the search smooths the field it
    // reads, so the smoothed divergence would be smoothed twice
    const bitvf::io::Field field = bitvf::io::readField(input, {"cellX", "cellY", "cellZ"}, "divergence", "weight");

    std::cout << "Cells: " << field.values.size() << std::endl;


    // ---------------------------------------------
    // ------------ Find the voids -----------------
    // ---------------------------------------------

    // No mean particle separation is given, so sigma is in coordinate
    // units, here Mpc/h. Normal prints one line with the lattice found
    bitvf::Config config;
    config.smoothing.sigma = 15.;
    config.smoothing.unit = bitvf::LengthUnit::Coordinate;
    config.verbosity = otswap::Verbosity::Normal;
    config.log = &std::cout;

    const bitvf::Result result = bitvf::findVoidsBoxFromField(field.centres, field.values, config, field.weights);


    // ----------------------------------------------
    // ------------ Summarize the result ------------
    // ----------------------------------------------

    // fieldCentre is the smoothed field at each void's centre, so the
    // deepest void has the lowest value
    const bitvf::BasinCatalog& voids = result.voids;
    double largest = 0., deepest = 0.;
    for (std::size_t b = 0; b < voids.nBasins; ++b) {
      largest = std::max(largest, voids.radius[b]);
      deepest = std::min(deepest, voids.fieldCentre[b]);
    }

    std::cout << "Voids: " << voids.nBasins << std::endl;
    std::cout << "Largest void radius: " << largest << " Mpc/h" << std::endl;
    std::cout << "Lowest field at a void centre: " << deepest << std::endl;
    std::cout << "Overdense basins: " << result.overdense.nBasins << std::endl;


    // ------------------------------------------
    // ------------ Write the output ------------
    // ------------------------------------------

    // The voids and the overdense basins, voids_field.dat and
    // overdense_basins_field.dat, with the column fieldCentre in place of
    // the densities
    using bitvf::io::Product;
    bitvf::io::writeVoids(bitvf::io::fileName(Product::Voids, output, "field.dat", false), result);
    bitvf::io::writeOverdense(bitvf::io::fileName(Product::Overdense, output, "field.dat", false), result);

    std::cout << "Written: output/voids_field.dat, output/overdense_basins_field.dat" << std::endl;

  }

  catch (const otswap::Error& e) { std::cerr << e.what() << std::endl; return 1; }
  catch (const bitvf::Error& e) { std::cerr << e.what() << std::endl; return 1; }

  return 0;
}
