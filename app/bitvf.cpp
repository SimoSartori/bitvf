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
 *  @file app/bitvf.cpp
 *
 *  @brief The bitvf executable: a void search run from a parameter file,
 *  on tracers (the whole pipeline: reconstruction, optional correction of
 *  redshift-space distortions, search) or on a file written earlier (a
 *  complete or mean displacement file of otswap, or a divergence file of
 *  bitvf).
 *
 *  Usage: bitvf <parameter-file>, bitvf --help, bitvf --version. Exit
 *  status 0 when every requested file is written, 1 on any error, with
 *  one line per error on stderr, "bitvf: error: <message>". Every check
 *  that needs no data comes before any data is read. Only box runs are
 *  implemented. The keys are listed at
 *  https://simosartori.github.io/bitvf/executable.html
 *
 *  Uses the public headers only, as an outside program would.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "otswap/OT.h"
#include "otswap/ParameterFile.h"
#include "otswap/RSD.h"
#include "otswap/io.h"

#include "bitvf/BitVF.h"
#include "bitvf/VoidFinder.h"
#include "bitvf/io.h"

#include "common.h"

namespace {

  const char* const kProgram = "bitvf";

  // ==========================================================================
  // The keys, and the runs each applies to
  // ==========================================================================

  /// The kinds of a box run: a pipeline from the tracers (none of the
  /// three input keys), or a search on the file of displacement_file,
  /// displacement_mean_file or divergence_input_file.
  enum Kind : unsigned { Pipeline = 1, Complete = 2, Mean = 4, Divergence = 8, Any = 15 };

  /// The switch a key needs, within the kinds it applies to.
  enum class Needs { Nothing, PerformRecons, RandomFile, NoDivergenceFile };

  struct Key {
    const char* name;
    unsigned kinds;     ///< the kinds of box run it applies to; 0: a key of a lightcone only
    Needs needs;
    bool boxOnly;       ///< does not apply to a lightcone
  };

  /// Every key, in the order of the typed reads, after the switches.
  const Key kKeys[] = {
    {"geometry",                      Any,                   Needs::Nothing,          false},
    {"verbosity",                     Any,                   Needs::Nothing,          false},
    {"tracer_file",                   Pipeline | Divergence, Needs::Nothing,          false},
    {"column_tracer",                 Pipeline | Divergence, Needs::Nothing,          false},
    {"id_column",                     Pipeline | Divergence, Needs::Nothing,          false},
    {"random_file",                   Pipeline,              Needs::Nothing,          false},
    {"column_random",                 Pipeline,              Needs::RandomFile,       false},
    {"units",                         0,                     Needs::Nothing,          false},
    {"random_recons_mode",            0,                     Needs::Nothing,          false},
    {"random_file_recons",            0,                     Needs::Nothing,          false},
    {"delimiter",                     Pipeline | Divergence, Needs::Nothing,          false},
    {"comment",                       Pipeline | Divergence, Needs::Nothing,          false},
    {"veto_mask_filename",            0,                     Needs::Nothing,          false},
    {"z_min",                         0,                     Needs::Nothing,          false},
    {"z_max",                         0,                     Needs::Nothing,          false},
    {"displacement_file",             Any,                   Needs::Nothing,          false},
    {"displacement_mean_file",        Any,                   Needs::Nothing,          true},
    {"divergence_input_file",         Any,                   Needs::Nothing,          true},
    {"divergence_file",               Complete | Mean,       Needs::Nothing,          false},
    {"Omega_matter",                  Pipeline,              Needs::PerformRecons,    false},
    {"little_h",                      Pipeline,              Needs::PerformRecons,    false},
    {"w0",                            Pipeline,              Needs::PerformRecons,    false},
    {"wa",                            Pipeline,              Needs::PerformRecons,    false},
    {"Nrec",                          Pipeline,              Needs::Nothing,          false},
    {"threshold",                     Pipeline,              Needs::Nothing,          false},
    {"seed",                          Pipeline,              Needs::Nothing,          false},
    {"mesh_cellsize",                 Pipeline,              Needs::Nothing,          false},
    {"bins",                          0,                     Needs::Nothing,          false},
    {"reject_mask_crossings",         0,                     Needs::Nothing,          false},
    {"max_unobserved_pixels_crossed", 0,                     Needs::Nothing,          false},
    {"box_min",                       Pipeline,              Needs::Nothing,          true},
    {"box_max",                       Pipeline,              Needs::Nothing,          true},
    {"mps",                           Pipeline,              Needs::Nothing,          true},
    {"sigma",                         Any,                   Needs::Nothing,          false},
    {"sigma_unit",                    Any,                   Needs::Nothing,          false},
    {"cell_size",                     Pipeline | Complete | Mean, Needs::NoDivergenceFile, false},
    {"perform_recons",                Pipeline,              Needs::Nothing,          false},
    {"save_zspace_voids",             Pipeline,              Needs::PerformRecons,    false},
    {"rsd_smoothing_scale",           Pipeline,              Needs::PerformRecons,    false},
    {"rsd_weight_by_realizations",    Pipeline,              Needs::PerformRecons,    false},
    {"redshift",                      Pipeline,              Needs::PerformRecons,    true},
    {"bias",                          Pipeline,              Needs::PerformRecons,    true},
    {"los_axis",                      Pipeline,              Needs::PerformRecons,    true},
    {"bias_file",                     0,                     Needs::Nothing,          false},
    {"output_dir",                    Any,                   Needs::Nothing,          false},
    {"output_file",                   Any,                   Needs::Nothing,          false},
    {"output_format",                 Any,                   Needs::Nothing,          false},
    {"printDisplMean",                Pipeline,              Needs::Nothing,          false},
    {"printDisplComplete",            Pipeline,              Needs::Nothing,          false},
    {"printMps",                      0,                     Needs::Nothing,          false},
    {"printDiv",                      Any,                   Needs::Nothing,          false},
    {"printGal",                      Any,                   Needs::Nothing,          false},
    {"printOverdense",                Any,                   Needs::Nothing,          false},
    {"printRecCat",                   Pipeline,              Needs::PerformRecons,    false},
  };

  /// The keys that choose the kind of a run; at most one is given.
  const char* const kInputKeys[] = {"displacement_file", "displacement_mean_file", "divergence_input_file"};

  /// The switches that decide which keys apply, read before the others.
  struct Switches {
    Kind kind = Pipeline;
    std::string inputKey;
    bool performRecons = false;
    bool randomFile = false;
    bool divergenceFile = false;
  };

  const Key* find_key (const std::string& name)
  {
    for (const Key& k : kKeys)
      if (name == k.name) return &k;
    return nullptr;
  }

  /// Why a key does not apply to a box run, or empty when it does. The
  /// first reason that holds: the geometry, the kind of run, a divergence
  /// file, a switch off.
  std::string box_rule (const std::string& name, const Switches& run)
  {
    const Key* key = find_key(name);
    if (key == nullptr) return "unknown key";
    if (key->kinds == 0) return "does not apply to a box";
    if ((key->kinds & run.kind) == 0)
      return run.kind == Pipeline ? "does not apply to a run without displacement_file or displacement_mean_file"
                                  : "does not apply to a run with " + run.inputKey;
    switch (key->needs) {
      case Needs::Nothing:
        break;
      case Needs::PerformRecons:
        if (!run.performRecons) return "does not apply to a run without perform_recons";
        break;
      case Needs::RandomFile:
        if (!run.randomFile) return "does not apply to a run without random_file";
        break;
      case Needs::NoDivergenceFile:
        if (run.divergenceFile) return "does not apply to a run with divergence_file";
        break;
    }
    return "";
  }

  /// Why a key does not apply to a lightcone run, or empty when it does.
  std::string lightcone_rule (const std::string& name)
  {
    const Key* key = find_key(name);
    if (key == nullptr) return "unknown key";
    if (key->boxOnly) return "does not apply to a lightcone";
    return "";
  }

  // ==========================================================================
  // Typed reads, with the executable's own checks
  // ==========================================================================

  /// Everything the parameter file sets, typed and checked.
  struct Run {
    Switches switches;
    otswap::Verbosity verbosity = otswap::Verbosity::Normal;

    std::string inputFile;
    std::string tracerFile;
    std::vector<std::string> tracerColumns;
    std::optional<std::string> idColumn;
    std::optional<std::string> randomFile;
    std::vector<std::string> randomColumns;
    char delimiter = ' ';
    char comment = '#';
    std::optional<std::string> divergenceFile;

    double omegaM = 0., h = 0., w0 = -1., wa = 0.;

    bitvf::BoxConfig config;

    std::string outputDir, outputName, extension;
    bool printMean = false, printComplete = false, printDiv = false;
    bool printGal = false, printOverdense = false, printCatalog = false;
  };

  /// The typed reads and the value checks, in the order of the keys, the
  /// first error stopping the run. Reads only the keys that apply to the
  /// run.
  Run read_run (const app::Parameters& p, const Switches& switches)
  {
    const otswap::io::ParameterFile& f = p.file();
    const Kind kind = switches.kind;
    const bool pipeline = kind == Pipeline;
    const bool recons = switches.performRecons;

    Run run;
    run.switches = switches;
    run.verbosity = p.verbosity();

    if (pipeline || kind == Divergence) {
      run.tracerFile = f.get<std::string>("tracer_file");
      run.tracerColumns = p.three<std::string>("column_tracer");
      run.idColumn = p.optional<std::string>("id_column");
      if (run.idColumn)
        for (const std::string& c : run.tracerColumns)
          if (app::same_column(*run.idColumn, c, app::is_fits(run.tracerFile)))
            p.fail("id_column", *run.idColumn + " is also a column of column_tracer, at line " +
                   p.line("column_tracer"));
    }
    if (pipeline) {
      run.randomFile = p.optional<std::string>("random_file");
      if (run.randomFile) run.randomColumns = p.three<std::string>("column_random");
    }
    if (pipeline || kind == Divergence) {
      run.delimiter = f.get<char>("delimiter", run.delimiter);
      run.comment = f.get<char>("comment", run.comment);
    }
    if (!pipeline) run.inputFile = f.get<std::string>(switches.inputKey);
    if (kind == Complete || kind == Mean) run.divergenceFile = p.optional<std::string>("divergence_file");

    if (recons) {
      run.omegaM = f.get<double>("Omega_matter");
      p.require("Omega_matter", run.omegaM > 0. && run.omegaM <= 1., "in (0, 1]");
      run.h = f.get<double>("little_h");
      p.require("little_h", run.h > 0., "positive");
      run.w0 = f.get<double>("w0", run.w0);
      run.wa = f.get<double>("wa", run.wa);
    }

    bitvf::BoxConfig& c = run.config;
    if (pipeline) {
      c.ot.nRealizations = f.get<unsigned>("Nrec", c.ot.nRealizations);
      p.require("Nrec", c.ot.nRealizations >= 1, "at least 1");
      c.ot.convergence = f.get<double>("threshold", c.ot.convergence);
      p.require("threshold", c.ot.convergence > 0. && c.ot.convergence < 1., "strictly between 0 and 1");
      c.ot.seed = f.get<unsigned>("seed", c.ot.seed);
      c.ot.cellSize = f.get<double>("mesh_cellsize", c.ot.cellSize);
      p.require("mesh_cellsize", c.ot.cellSize > 0., "positive");

      const std::vector<double> lower = p.optional_three("box_min");
      const std::vector<double> upper = p.optional_three("box_max");
      if (lower.empty() != upper.empty()) {
        const std::string given = lower.empty() ? "box_max" : "box_min";
        const std::string other = lower.empty() ? "box_min" : "box_max";
        if (!f.has(other)) p.fail(given, other + " is missing; give both or neither");
        p.fail(other, "empty while " + given + ", at line " + p.line(given) + ", is not; give both or neither");
      }
      if (!lower.empty()) {
        otswap::BoxLimits limits;
        for (std::size_t k = 0; k < 3; ++k) {
          if (!(upper[k] > lower[k]))
            p.fail("box_max", "element " + std::to_string(k+1) + ", " + p.written("box_max", k) +
                   ", is not above element " + std::to_string(k+1) + " of box_min, " +
                   p.written("box_min", k) + ", at line " + p.line("box_min"));
          limits.lower[k] = lower[k];
          limits.upper[k] = upper[k];
        }
        c.limits = limits;
      }
      c.mps = p.optional<double>("mps");
      if (c.mps) p.require("mps", *c.mps > 0., "positive");
    }

    bitvf::Config& v = c.voids;
    v.smoothing.sigma = f.get<double>("sigma", v.smoothing.sigma);
    p.require("sigma", v.smoothing.sigma >= 0., "non-negative");
    if (f.has("sigma_unit"))
      v.smoothing.unit = p.word("sigma_unit", {"mps", "coordinate"}) == "mps" ? bitvf::LengthUnit::Mps
                                                                             : bitvf::LengthUnit::Coordinate;
    if (kind != Divergence && !switches.divergenceFile) {
      v.divergenceCellSize = p.optional<double>("cell_size");
      if (v.divergenceCellSize) p.require("cell_size", *v.divergenceCellSize > 0., "positive");
    }

    if (recons) {
      bitvf::BoxRsd rsd;
      rsd.zspaceVoids = f.get<bool>("save_zspace_voids", rsd.zspaceVoids);
      rsd.sigma = f.get<double>("rsd_smoothing_scale");
      p.require("rsd_smoothing_scale", rsd.sigma >= 0., "non-negative");
      rsd.correction.weightByRealizations = f.get<bool>("rsd_weight_by_realizations",
                                                        rsd.correction.weightByRealizations);
      rsd.redshift = f.get<double>("redshift");
      p.require("redshift", rsd.redshift >= 0. && rsd.redshift <= 10.,
                "in [0, 10], the range of the distance table");
      rsd.bias = f.get<double>("bias");
      p.require("bias", rsd.bias > 0., "positive");
      const unsigned axis = f.get<unsigned>("los_axis");
      p.require("los_axis", axis <= 2, "0, 1 or 2");
      rsd.axis = axis;
      rsd.correction.verbosity = run.verbosity;
      c.rsd = rsd;
    }

    run.outputDir = f.get<std::string>("output_dir");
    run.outputName = p.output_file();
    run.extension = p.word("output_format", {"dat", "fits"}, "dat");
    if (pipeline) {
      run.printMean = f.get<bool>("printDisplMean", false);
      run.printComplete = f.get<bool>("printDisplComplete", false);
    }
    run.printDiv = f.get<bool>("printDiv", false);
    run.printGal = f.get<bool>("printGal", false);
    run.printOverdense = f.get<bool>("printOverdense", false);
    if (recons) run.printCatalog = f.get<bool>("printRecCat", false);

    c.ot.verbosity = run.verbosity;
    v.verbosity = run.verbosity;
    v.keepDivergence = run.printDiv;
    v.keepAssignment = run.printGal;
    c.keepReconstructions = run.printMean || run.printComplete || run.printCatalog;

    return run;
  }

  // ==========================================================================
  // The run
  // ==========================================================================

  /// The tracers of tracer_file, and their IDs when id_column is set.
  void read_tracers (const Run& run, std::vector<double>& tracers, std::vector<std::int64_t>& ids)
  {
    if (run.idColumn) {
      otswap::io::Table table = app::call_otswap([&] {
        return otswap::io::read(run.tracerFile, run.tracerColumns, {*run.idColumn}, run.delimiter, run.comment);
      });
      tracers = std::move(table.values);
      ids = std::move(table.integers);
      app::check_ids(ids, run.tracerFile, *run.idColumn);
    }
    else
      tracers = app::call_otswap([&] {
        return otswap::io::read(run.tracerFile, run.tracerColumns, run.delimiter, run.comment);
      }).values;
  }

  /// One run from a parameter file. Every key that applies is read before
  /// the data; one key left after the reads is a bug here, reported as an
  /// unknown key. The files are written after the run, the otswap files
  /// first, then the redshift-space pass's, then the final pass's.
  int run_file (const std::string& file)
  {
    const otswap::io::ParameterFile parameters(file);
    const app::Parameters p(parameters, file);

    if (p.word("geometry", {"box", "lightcone"}) == "lightcone") {
      if (app::report_errors(kProgram, app::irrelevant_keys(parameters, file, lightcone_rule))) return 1;
      p.fail("geometry", "lightcone runs are not implemented yet");
    }

    Switches switches;
    std::vector<std::string> given;
    for (const char* key : kInputKeys)
      if (parameters.has(key)) given.push_back(key);
    std::sort(given.begin(), given.end(), [&parameters] (const std::string& a, const std::string& b) {
      return parameters.line(a) < parameters.line(b);
    });
    if (given.size() > 1)
      p.fail(given[1], "cannot be given with " + given[0] + ", at line " + p.line(given[0]) +
             "; a run reads one of displacement_file, displacement_mean_file and divergence_input_file, "
             "or none for a pipeline run");
    if (!given.empty()) {
      switches.inputKey = given[0];
      switches.kind = given[0] == "displacement_file" ? Complete : given[0] == "displacement_mean_file" ? Mean
                                                                                                       : Divergence;
    }
    if (switches.kind == Pipeline) {
      switches.performRecons = parameters.get<bool>("perform_recons", false);
      switches.randomFile = parameters.has("random_file") &&
                            !parameters.get<std::vector<std::string>>("random_file").empty();
    }
    if ((switches.kind == Complete || switches.kind == Mean) && parameters.has("divergence_file"))
      switches.divergenceFile = !parameters.get<std::vector<std::string>>("divergence_file").empty();

    if (app::report_errors(kProgram, app::irrelevant_keys(parameters, file, [&switches] (const std::string& key) {
      return box_rule(key, switches);
    })))
      return 1;

    const Run run = read_run(p, switches);
    if (app::report_unread(kProgram, parameters, file)) return 1;

    app::create_directory(p, run.outputDir);

    const bool report = run.verbosity != otswap::Verbosity::Silent;
    auto wrote = [report] (const std::string& written) {
      if (report) std::clog << kProgram << ": wrote " << written << std::endl;
    };
    auto otswap_path = [&run] (const std::string& prefix) {
      return (std::filesystem::path(run.outputDir) / (prefix + run.outputName + "." + run.extension)).string();
    };
    auto write_search = [&] (const bitvf::Result& result, const bool zspace) {
      const std::string name = run.outputName + "." + run.extension;
      std::string path = bitvf::io::fileName(bitvf::io::Product::Voids, run.outputDir, name, zspace);
      bitvf::io::writeVoids(path, result);
      wrote(path);
      if (run.printOverdense) {
        path = bitvf::io::fileName(bitvf::io::Product::Overdense, run.outputDir, name, zspace);
        bitvf::io::writeOverdense(path, result);
        wrote(path);
      }
      if (run.printDiv) {
        path = bitvf::io::fileName(bitvf::io::Product::Divergence, run.outputDir, name, zspace);
        bitvf::io::writeDivergence(path, result);
        wrote(path);
      }
      if (run.printGal) {
        path = bitvf::io::fileName(bitvf::io::Product::Galaxies, run.outputDir, name, zspace);
        bitvf::io::writeGalaxies(path, result);
        wrote(path);
      }
    };
    auto write_reconstruction = [&] (const otswap::Result& result, const std::string& zspace) {
      if (run.printMean) {
        const std::string path = otswap_path("displacement_mean_" + zspace);
        app::call_otswap([&] { otswap::io::writeDisplacementMean(path, result); });
        wrote(path);
      }
      if (run.printComplete) {
        const std::string path = otswap_path("displacement_complete_" + zspace);
        app::call_otswap([&] { otswap::io::writeDisplacementComplete(path, result, true); });
        wrote(path);
      }
    };

    std::optional<bitvf::DivergenceGrid> grid;
    const bitvf::Config& config = run.config.voids;

    switch (run.switches.kind) {
      case Pipeline: {
        std::vector<double> tracers, randoms;
        std::vector<std::int64_t> ids;
        read_tracers(run, tracers, ids);
        if (run.randomFile)
          randoms = app::call_otswap([&] {
            return otswap::io::read(*run.randomFile, run.randomColumns, run.delimiter, run.comment);
          }).values;

        bitvf::BitVFResult result;
        if (run.switches.performRecons) {
          const otswap::DistanceTable distances = app::call_otswap([&] {
            return otswap::DistanceTable(run.omegaM, run.h, run.w0, run.wa);
          });
          result = bitvf::runBox(tracers, randoms, distances, run.config, ids);
          if (result.zspaceReconstruction) write_reconstruction(*result.zspaceReconstruction, "zspace_");
        }
        else
          result = bitvf::runBox(tracers, randoms, run.config, ids);
        if (result.reconstruction) write_reconstruction(*result.reconstruction, "");
        if (run.printCatalog) {
          const std::string path = otswap_path("reconstructed_catalogue_");
          app::call_otswap([&] { otswap::io::writeRealSpaceCatalog(path, *result.realSpace); });
          wrote(path);
        }
        if (result.zspaceVoids) write_search(*result.zspaceVoids, true);
        write_search(result.voids, false);
        break;
      }
      case Complete: {
        const otswap::Result reconstruction = app::call_otswap([&] {
          return otswap::io::readDisplacementComplete(run.inputFile);
        });
        if (run.divergenceFile) grid = bitvf::io::readDivergence(*run.divergenceFile);
        write_search(bitvf::findVoidsBox(reconstruction, config, grid ? &*grid : nullptr), false);
        break;
      }
      case Mean: {
        const bitvf::io::DisplacementMean mean = bitvf::io::readDisplacementMean(run.inputFile);
        if (run.divergenceFile) grid = bitvf::io::readDivergence(*run.divergenceFile);
        bitvf::Config c = config;
        c.nRealizationsMean = mean.nRealizationsAveraged;
        write_search(bitvf::findVoidsBox(mean.tracers, mean.displacement, mean.valid, mean.mps, c, mean.limits,
                                         mean.ids, grid ? &*grid : nullptr), false);
        break;
      }
      case Divergence: {
        const bitvf::DivergenceGrid divergence = bitvf::io::readDivergence(run.inputFile);
        std::vector<double> tracers;
        std::vector<std::int64_t> ids;
        read_tracers(run, tracers, ids);
        write_search(bitvf::findVoidsBoxFromDivergence(divergence, tracers, config, ids), false);
        break;
      }
      case Any:
        break;
    }

    return 0;
  }

}


// ============================================================================


int main (int argc, char** argv)
{
  return app::main(argc, argv, kProgram, BITVF_VERSION, run_file);
}
