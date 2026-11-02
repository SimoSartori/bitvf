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
 *  @file app/bitvf_field.cpp
 *
 *  @brief The bitvf-field executable: the voids of a field given at the
 *  centres of the cells of a box lattice, with no tracers
 *  (bitvf::findVoidsBoxFromField), run from a parameter file.
 *
 *  Usage: bitvf-field <parameter-file>, bitvf-field --help, bitvf-field
 *  --version. Exit status 0 when every requested file is written, 1 on
 *  any error, with one line per error on stderr, "bitvf-field: error:
 *  <message>". Every check that needs no data comes before any data is
 *  read. The library's own report lines start with "bitvf: ". The keys
 *  are listed at https://simosartori.github.io/bitvf/executable.html
 *
 *  Uses the public headers only, as an outside program would.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <array>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "otswap/OT.h"
#include "otswap/ParameterFile.h"

#include "bitvf/VoidFinder.h"
#include "bitvf/io.h"

#include "common.h"

namespace {

  const char* const kProgram = "bitvf-field";

  /// Every key, in the order of the typed reads.
  const char* const kKeys[] = {
    "verbosity", "field_file", "column_centre", "column_value", "column_weight", "delimiter", "comment", "mps",
    "sigma", "sigma_unit", "output_dir", "output_file", "output_format", "printOverdense", "printDiv",
    "printDivUnsmoothed",
  };

  /// Why a key does not apply to the run, or empty when it does.
  std::string rule (const std::string& key, const bool printDiv)
  {
    for (const char* k : kKeys)
      if (key == k) return key == "printDivUnsmoothed" && !printDiv ? "does not apply to a run without printDiv" : "";
    return "unknown key";
  }

  /// One run from a parameter file. Every key that applies is read before
  /// the data; one key left after the reads is a bug here, reported as an
  /// unknown key.
  int run_file (const std::string& file)
  {
    const otswap::io::ParameterFile parameters(file);
    const app::Parameters p(parameters, file);
    const otswap::io::ParameterFile& f = parameters;

    const bool printDiv = f.get<bool>("printDiv", false);
    if (app::report_errors(kProgram, app::irrelevant_keys(parameters, file, [printDiv] (const std::string& key) {
      return rule(key, printDiv);
    })))
      return 1;

    bitvf::Config config;
    config.verbosity = p.verbosity();
    const std::string fieldFile = f.get<std::string>("field_file");
    const std::vector<std::string> centre = p.three<std::string>("column_centre");
    const std::string value = f.get<std::string>("column_value");
    const std::optional<std::string> weight = p.optional<std::string>("column_weight");
    const char delimiter = f.get<char>("delimiter", ' ');
    const char comment = f.get<char>("comment", '#');
    const std::optional<double> mps = p.optional<double>("mps");
    if (mps) p.require("mps", *mps > 0., "positive");
    config.smoothing.sigma = f.get<double>("sigma", config.smoothing.sigma);
    p.require("sigma", config.smoothing.sigma >= 0., "non-negative");
    if (f.has("sigma_unit"))
      config.smoothing.unit = p.word("sigma_unit", {"mps", "coordinate"}) == "mps" ? bitvf::LengthUnit::Mps
                                                                                  : bitvf::LengthUnit::Coordinate;
    if (config.smoothing.sigma > 0. && config.smoothing.unit == bitvf::LengthUnit::Mps && !mps) {
      if (f.has("sigma_unit"))
        p.fail("sigma_unit", "mps needs a value of mps; set mps, or sigma_unit = coordinate, or sigma = 0");
      throw app::Error(file + ": sigma is in mps (the default of sigma_unit), but mps has no value; set mps, or "
                       "sigma_unit = coordinate, or sigma = 0");
    }
    const std::string outputDir = f.get<std::string>("output_dir");
    const std::string outputName = p.output_file();
    const std::string extension = p.word("output_format", {"dat", "fits"}, "dat");
    const bool printOverdense = f.get<bool>("printOverdense", false);
    config.keepDivergence = printDiv;
    if (printDiv) config.keepUnsmoothedDivergence = f.get<bool>("printDivUnsmoothed", false);

    if (app::report_unread(kProgram, parameters, file)) return 1;

    app::create_directory(p, outputDir);

    const bitvf::io::Field field = bitvf::io::readField(fieldFile, {centre[0], centre[1], centre[2]}, value,
                                                        weight ? *weight : "", delimiter, comment);
    const bitvf::Result result = bitvf::findVoidsBoxFromField(field.centres, field.values, config, field.weights,
                                                              mps);

    const bool report = config.verbosity != otswap::Verbosity::Silent;
    const std::string name = outputName + "." + extension;
    auto write = [&] (const bitvf::io::Product product, void (*writer) (const std::string&, const bitvf::Result&)) {
      const std::string path = bitvf::io::fileName(product, outputDir, name, false);
      writer(path, result);
      if (report) std::clog << kProgram << ": wrote " << path << std::endl;
    };
    write(bitvf::io::Product::Voids, bitvf::io::writeVoids);
    if (printOverdense) write(bitvf::io::Product::Overdense, bitvf::io::writeOverdense);
    if (printDiv) write(bitvf::io::Product::Divergence, bitvf::io::writeDivergence);

    return 0;
  }

}


// ============================================================================


int main (int argc, char** argv)
{
  return app::main(argc, argv, kProgram, BITVF_VERSION, run_file);
}
