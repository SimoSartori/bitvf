// SPDX-License-Identifier: GPL-2.0-or-later
/********************************************************************
 * Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it     *
 *                                                                  *
 * This program is free software; you can redistribute it and/or    *
 * modify it under the terms of the GNU General Public License as   *
 * published by the Free Software Foundation; either version 2 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU General Public License for more details.                     *
 *                                                                  *
 * You should have received a copy of the GNU General Public        *
 * License along with this program; if not, write to the Free       *
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,      *
 * Boston, MA 02110-1301 USA.                                       *
 ********************************************************************/

/**
 *  @file tests/executable/reference.cpp
 *
 *  @brief The files each kind of run of the executables must write, made
 *  with the library calls the parameters of its case name.
 *
 *  Usage: executable_reference <data directory> <output directory>. Writes
 *  one directory per case kind, named as the "#! same" lines of the cases
 *  of tests/executable/cases read them; the runs from a file read the
 *  files of pipeline/. Not a test by itself.
 */

#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "otswap/OT.h"
#include "otswap/RSD.h"
#include "otswap/io.h"

#include "bitvf/BitVF.h"
#include "bitvf/VoidFinder.h"
#include "bitvf/io.h"

namespace {

  std::string g_out;

  std::string path (const std::string& kind, const std::string& prefix, const std::string& extension)
  {
    return g_out + "/" + kind + "/" + prefix + "box." + extension;
  }

  void write_search (const std::string& kind, const bitvf::Result& r, const std::string& extension,
                     const bool zspace = false)
  {
    const std::string z = zspace ? "zspace_" : "";
    bitvf::io::writeVoids(path(kind, "voids_" + z, extension), r);
    bitvf::io::writeOverdense(path(kind, "overdense_basins_" + z, extension), r);
    if (r.divergence) bitvf::io::writeDivergence(path(kind, "divergence_" + z, extension), r);
    if (r.assignment) bitvf::io::writeGalaxies(path(kind, "galaxies_" + z, extension), r);
  }

  bitvf::Config search_config ()
  {
    bitvf::Config c;
    c.verbosity = otswap::Verbosity::Silent;
    c.keepDivergence = true;
    c.keepAssignment = true;
    return c;
  }

}

int main (int argc, char** argv)
{
  if (argc != 3) {
    std::cerr << "usage: executable_reference <data directory> <output directory>\n";
    return 1;
  }
  const std::string data = argv[1];
  g_out = argv[2];
  for (const char* kind : {"pipeline", "rsd", "complete", "mean", "divergence", "field"})
    std::filesystem::create_directories(g_out + "/" + kind);

  try {
    const otswap::io::Table table = otswap::io::read(data + "/box_tracers.dat", {"x", "y", "z"}, {"id"});

    // pipeline.ini: Nrec 2, seed 5, IDs, every output, FITS.
    {
      bitvf::BoxConfig c;
      c.ot.nRealizations = 2;
      c.ot.seed = 5;
      c.ot.verbosity = otswap::Verbosity::Silent;
      c.voids = search_config();
      c.keepReconstructions = true;
      const bitvf::BitVFResult r = bitvf::runBox(table.values, {}, c, table.integers);
      otswap::io::writeDisplacementMean(path("pipeline", "displacement_mean_", "fits"), *r.reconstruction);
      otswap::io::writeDisplacementComplete(path("pipeline", "displacement_complete_", "fits"), *r.reconstruction,
                                            true);
      write_search("pipeline", r.voids, "fits");
    }

    // rsd.ini: the correction along z, the redshift-space voids, the
    // corrected catalogue, ASCII, no IDs; verbosity normal, which otswap's
    // files record.
    {
      bitvf::BoxConfig c;
      c.ot.nRealizations = 2;
      c.ot.seed = 9;
      c.ot.log = nullptr;
      c.voids = search_config();
      c.voids.verbosity = otswap::Verbosity::Normal;
      c.voids.log = nullptr;
      c.keepReconstructions = true;
      bitvf::BoxRsd rsd;
      rsd.axis = 2;
      rsd.redshift = 0.5;
      rsd.bias = 1.5;
      rsd.sigma = 5.;
      rsd.zspaceVoids = true;
      rsd.correction.log = nullptr;
      c.rsd = rsd;
      const otswap::DistanceTable distances(0.3, 0.7);
      const bitvf::BitVFResult r = bitvf::runBox(table.values, {}, distances, c);
      otswap::io::writeDisplacementMean(path("rsd", "displacement_mean_zspace_", "dat"), *r.zspaceReconstruction);
      otswap::io::writeDisplacementMean(path("rsd", "displacement_mean_", "dat"), *r.reconstruction);
      otswap::io::writeRealSpaceCatalog(path("rsd", "reconstructed_catalogue_", "dat"), *r.realSpace);
      write_search("rsd", *r.zspaceVoids, "dat", true);
      write_search("rsd", r.voids, "dat");
    }

    // complete.ini: pipeline's complete file on pipeline's grid, sigma 0.5,
    // ASCII.
    {
      const otswap::Result reconstruction =
        otswap::io::readDisplacementComplete(path("pipeline", "displacement_complete_", "fits"));
      const bitvf::DivergenceGrid grid = bitvf::io::readDivergence(path("pipeline", "divergence_", "fits"));
      bitvf::Config c = search_config();
      c.smoothing.sigma = 0.5;
      write_search("complete", bitvf::findVoidsBox(reconstruction, c, &grid), "dat");
    }

    // mean.ini: pipeline's mean file, cell 1.2 mps, sigma 3 in coordinate
    // units, FITS.
    {
      const bitvf::io::DisplacementMean m = bitvf::io::readDisplacementMean(path("pipeline", "displacement_mean_",
                                                                                 "fits"));
      bitvf::Config c = search_config();
      c.divergenceCellSize = 1.2;
      c.smoothing.sigma = 3.;
      c.smoothing.unit = bitvf::LengthUnit::Coordinate;
      c.nRealizationsMean = m.nRealizationsAveraged;
      write_search("mean", bitvf::findVoidsBox(m.tracers, m.displacement, m.valid, m.mps, c, m.limits, m.ids),
                   "fits");
    }

    // divergence.ini: pipeline's divergence with the tracers and their IDs,
    // sigma 2, ASCII.
    {
      const bitvf::DivergenceGrid grid = bitvf::io::readDivergence(path("pipeline", "divergence_", "fits"));
      bitvf::Config c = search_config();
      c.smoothing.sigma = 2.;
      write_search("divergence", bitvf::findVoidsBoxFromDivergence(grid, table.values, c, table.integers), "dat");
    }

    // field.ini: field.dat with its weights, sigma 1 in coordinate units,
    // ASCII.
    {
      const bitvf::io::Field f = bitvf::io::readField(data + "/field.dat", {"x", "y", "z"}, "value", "weight");
      bitvf::Config c;
      c.verbosity = otswap::Verbosity::Silent;
      c.smoothing.unit = bitvf::LengthUnit::Coordinate;
      c.keepDivergence = true;
      write_search("field", bitvf::findVoidsBoxFromField(f.centres, f.values, c, f.weights), "dat");
    }
  }
  catch (const std::exception& e) {
    std::cerr << "executable_reference: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
