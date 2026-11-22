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
#include <fstream>
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
    return g_out + "/" + kind + "/" + prefix + (kind.compare(0, 2, "lc") == 0 ? "lightcone." : "box.") + extension;
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
  for (const char* kind : {"pipeline", "rsd", "complete", "mean", "meanmulti", "meanold", "divergence", "field",
                           "lcpipeline", "lcrsd", "lccomplete", "lcmean", "lcexternal"})
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
      const otswap::Result m = otswap::io::readDisplacementMean(path("pipeline", "displacement_mean_", "fits"));
      bitvf::Config c = search_config();
      c.divergenceCellSize = 1.2;
      c.smoothing.sigma = 3.;
      c.smoothing.unit = bitvf::LengthUnit::Coordinate;
      c.nRealizationsMean = m.nRealizationsMean;
      write_search("mean", bitvf::findVoidsBox(m, c), "fits");
    }

    // meanmulti.ini: pipeline's mean file written again from the Result
    // read (NREC 1, NRECMEAN 2), searched as mean.ini does.
    {
      const otswap::Result m = otswap::io::readDisplacementMean(path("pipeline", "displacement_mean_", "fits"));
      otswap::io::writeDisplacementMean(path("meanmulti", "displacement_mean_", "fits"), m);
      const otswap::Result again = otswap::io::readDisplacementMean(path("meanmulti", "displacement_mean_", "fits"));
      bitvf::Config c = search_config();
      c.divergenceCellSize = 1.2;
      c.smoothing.sigma = 3.;
      c.smoothing.unit = bitvf::LengthUnit::Coordinate;
      c.nRealizationsMean = again.nRealizationsMean;
      write_search("meanmulti", bitvf::findVoidsBox(again, c), "fits");
    }

    // meanold.ini: pipeline's mean file as ASCII without its box keywords,
    // as otswap 0.2 wrote it.
    {
      const otswap::Result m = otswap::io::readDisplacementMean(path("pipeline", "displacement_mean_", "fits"));
      const std::string full = path("meanold", "full_", "dat");
      otswap::io::writeDisplacementMean(full, m);
      std::ifstream in(full);
      std::ofstream out(path("meanold", "displacement_mean_", "dat"));
      for (std::string line; std::getline(in, line); )
        if (line.compare(0, 6, "## BOX") != 0) out << line << '\n';
    }

    // divergence.ini: pipeline's divergence with the tracers and their IDs,
    // sigma 2, ASCII.
    {
      const bitvf::DivergenceGrid grid = bitvf::io::readDivergence(path("pipeline", "divergence_", "fits"));
      bitvf::Config c = search_config();
      c.smoothing.sigma = 2.;
      write_search("divergence", bitvf::findVoidsBoxFromDivergence(grid, table.values, c, table.integers), "dat");
    }

    // lcpipeline.ini: a lightcone with IDs, the redshift cut and every
    // output, FITS.
    const otswap::io::Table skyTracers = otswap::io::read(data + "/lightcone_tracers.dat", {"ra", "dec", "z"}, {"id"});
    const std::vector<double> sky = otswap::skyToRadians(skyTracers.values);
    const std::vector<double> skyRandoms =
      otswap::skyToRadians(otswap::io::read(data + "/lightcone_randoms.dat", {"ra", "dec", "z"}).values);
    const otswap::Mask mask(data + "/lightcone_mask.fits");
    const otswap::DistanceTable lightDistances(0.3, 0.7, -1., 0.);
    auto lightcone_config = [] (const unsigned seed) {
      bitvf::LightconeConfig c;
      c.ot.nRealizations = 2;
      c.ot.seed = seed;
      c.ot.verbosity = otswap::Verbosity::Silent;
      c.nBins = 1;
      c.cut = {0.26, 0.39};
      c.voids = search_config();
      c.voids.divergenceCellSize = 0.5;
      c.keepReconstructions = true;
      return c;
    };
    {
      const bitvf::BitVFResult r =
        bitvf::runLightcone(sky, skyRandoms, mask, lightDistances, lightcone_config(5), {}, skyTracers.integers);
      otswap::io::writeDisplacementMean(path("lcpipeline", "displacement_mean_", "fits"), *r.reconstruction);
      otswap::io::writeDisplacementComplete(path("lcpipeline", "displacement_complete_", "fits"),
                                            *r.reconstruction, true);
      write_search("lcpipeline", r.voids, "fits");
    }

    // lcrsd.ini: the correction with resampled randoms, the redshift-space
    // voids, the corrected catalogue and the profiles, ASCII; verbosity
    // normal, which otswap's files record.
    {
      bitvf::LightconeConfig c = lightcone_config(9);
      c.ot.verbosity = otswap::Verbosity::Normal;
      c.ot.log = nullptr;
      c.voids.verbosity = otswap::Verbosity::Normal;
      c.voids.log = nullptr;
      bitvf::LightconeRsd rsd;
      rsd.bias = otswap::io::readBiasTable(data + "/lightcone_bias.dat");
      rsd.sigma = 20.;
      rsd.zspaceVoids = true;
      rsd.correction.log = nullptr;
      c.rsd = rsd;
      const bitvf::BitVFResult r = bitvf::runLightcone(sky, skyRandoms, mask, lightDistances, c);
      otswap::io::writeMpsProfile(path("lcrsd", "MPS_zspace_", "dat"), r.zspaceReconstruction->mpsProfile);
      otswap::io::writeMpsProfile(path("lcrsd", "MPS_", "dat"), r.reconstruction->mpsProfile);
      otswap::io::writeRealSpaceCatalog(path("lcrsd", "reconstructed_catalogue_", "dat"), *r.realSpace);
      write_search("lcrsd", *r.zspaceVoids, "dat", true);
      write_search("lcrsd", r.voids, "dat");
    }

    // lcexternal.ini: the correction with the randoms of the second
    // reconstruction given, ASCII.
    {
      bitvf::LightconeConfig c = lightcone_config(3);
      c.keepReconstructions = false;
      bitvf::LightconeRsd rsd;
      rsd.correction.verbosity = otswap::Verbosity::Silent;
      rsd.bias = otswap::io::readBiasTable(data + "/lightcone_bias.dat");
      rsd.sigma = 20.;
      rsd.secondPassRandoms = bitvf::SecondPassRandoms::External;
      c.rsd = rsd;
      const bitvf::BitVFResult r = bitvf::runLightcone(sky, skyRandoms, mask, lightDistances, c, skyRandoms);
      write_search("lcexternal", r.voids, "dat");
    }

    // lccomplete.ini: lcpipeline's complete file on lcpipeline's grid,
    // sigma 0.5, ASCII.
    {
      const otswap::Result reconstruction =
        otswap::io::readDisplacementComplete(path("lcpipeline", "displacement_complete_", "fits"));
      const bitvf::DivergenceGrid grid = bitvf::io::readDivergence(path("lcpipeline", "divergence_", "fits"));
      bitvf::Config c = search_config();
      c.smoothing.sigma = 0.5;
      write_search("lccomplete", bitvf::findVoidsLightcone(reconstruction, skyRandoms, mask, c, &grid), "dat");
    }

    // lcmean.ini: lcpipeline's mean file, cell 0.8 mps, sigma 3 in
    // coordinate units, FITS.
    {
      const otswap::Result m = otswap::io::readDisplacementMean(path("lcpipeline", "displacement_mean_", "fits"));
      bitvf::Config c = search_config();
      c.divergenceCellSize = 0.8;
      c.smoothing.sigma = 3.;
      c.smoothing.unit = bitvf::LengthUnit::Coordinate;
      c.nRealizationsMean = m.nRealizationsMean;
      write_search("lcmean", bitvf::findVoidsLightcone(m, skyRandoms, mask, c), "fits");
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
