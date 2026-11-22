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
 *  that needs no data comes before any data is read. The geometry is a
 *  box or a lightcone. The keys are listed at
 *  https://simosartori.github.io/bitvf/executable.html
 *
 *  Uses the public headers only, as an outside program would.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
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

  /// The kinds of a run: a pipeline from the tracers (none of the three
  /// input keys), or a search on the file of displacement_file,
  /// displacement_mean_file or divergence_input_file (a box only).
  enum Kind : unsigned { Pipeline = 1, Complete = 2, Mean = 4, Divergence = 8, Any = 15 };

  /// The switch a key needs, within the kinds it applies to.
  enum class Needs { Nothing, PerformRecons, RandomFile, NoDivergenceFile, RejectCrossings, ExternalRandoms };

  /// Where a key applies in one geometry: the kinds of run, 0 for none,
  /// and the switch it needs.
  struct Rule {
    unsigned kinds;
    Needs needs;
  };

  struct Key {
    const char* name;
    Rule box;
    Rule lightcone;
  };

  constexpr unsigned kFiles = Pipeline | Complete | Mean;
  constexpr Rule kNever = {0, Needs::Nothing};

  /// Every key, in the order of the typed reads, after the switches.
  const Key kKeys[] = {
    {"geometry",                      {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"verbosity",                     {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"tracer_file",                   {Pipeline | Divergence, Needs::Nothing}, {Pipeline, Needs::Nothing}},
    {"column_tracer",                 {Pipeline | Divergence, Needs::Nothing}, {Pipeline, Needs::Nothing}},
    {"id_column",                     {Pipeline | Divergence, Needs::Nothing}, {Pipeline, Needs::Nothing}},
    {"random_file",                   {Pipeline, Needs::Nothing},              {kFiles, Needs::Nothing}},
    {"column_random",                 {Pipeline, Needs::RandomFile},           {kFiles, Needs::Nothing}},
    {"units",                         kNever,                                  {kFiles, Needs::Nothing}},
    {"random_recons_mode",            kNever,                                  {Pipeline, Needs::PerformRecons}},
    {"random_file_recons",            kNever,                                  {Pipeline, Needs::ExternalRandoms}},
    {"delimiter",                     {Pipeline | Divergence, Needs::Nothing}, {kFiles, Needs::Nothing}},
    {"comment",                       {Pipeline | Divergence, Needs::Nothing}, {kFiles, Needs::Nothing}},
    {"veto_mask_filename",            kNever,                                  {kFiles, Needs::Nothing}},
    {"z_min",                         kNever,                                  {Pipeline, Needs::Nothing}},
    {"z_max",                         kNever,                                  {Pipeline, Needs::Nothing}},
    {"displacement_file",             {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"displacement_mean_file",        {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"divergence_input_file",         {Any, Needs::Nothing},                   kNever},
    {"divergence_file",               {Complete | Mean, Needs::Nothing},       {Complete | Mean, Needs::Nothing}},
    {"Omega_matter",                  {Pipeline, Needs::PerformRecons},        {Pipeline, Needs::Nothing}},
    {"little_h",                      {Pipeline, Needs::PerformRecons},        {Pipeline, Needs::Nothing}},
    {"w0",                            {Pipeline, Needs::PerformRecons},        {Pipeline, Needs::Nothing}},
    {"wa",                            {Pipeline, Needs::PerformRecons},        {Pipeline, Needs::Nothing}},
    {"Nrec",                          {Pipeline, Needs::Nothing},              {Pipeline, Needs::Nothing}},
    {"threshold",                     {Pipeline, Needs::Nothing},              {Pipeline, Needs::Nothing}},
    {"seed",                          {Pipeline, Needs::Nothing},              {Pipeline, Needs::Nothing}},
    {"mesh_cellsize",                 {Pipeline, Needs::Nothing},              {Pipeline, Needs::Nothing}},
    {"bins",                          kNever,                                  {Pipeline, Needs::Nothing}},
    {"reject_mask_crossings",         kNever,                                  {Pipeline, Needs::Nothing}},
    {"max_unobserved_pixels_crossed", kNever,                                  {Pipeline, Needs::RejectCrossings}},
    {"box_min",                       {Pipeline, Needs::Nothing},              kNever},
    {"box_max",                       {Pipeline, Needs::Nothing},              kNever},
    {"mps",                           {Pipeline, Needs::Nothing},              kNever},
    {"sigma",                         {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"sigma_unit",                    {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"cell_size",                     {kFiles, Needs::NoDivergenceFile},       {kFiles, Needs::NoDivergenceFile}},
    {"perform_recons",                {Pipeline, Needs::Nothing},              {Pipeline, Needs::Nothing}},
    {"save_zspace_voids",             {Pipeline, Needs::PerformRecons},        {Pipeline, Needs::PerformRecons}},
    {"rsd_smoothing_scale",           {Pipeline, Needs::PerformRecons},        {Pipeline, Needs::PerformRecons}},
    {"rsd_weight_by_realizations",    {Pipeline, Needs::PerformRecons},        {Pipeline, Needs::PerformRecons}},
    {"redshift",                      {Pipeline, Needs::PerformRecons},        kNever},
    {"bias",                          {Pipeline, Needs::PerformRecons},        kNever},
    {"los_axis",                      {Pipeline, Needs::PerformRecons},        kNever},
    {"bias_file",                     kNever,                                  {Pipeline, Needs::PerformRecons}},
    {"output_dir",                    {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"output_file",                   {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"output_format",                 {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"printDisplMean",                {Pipeline, Needs::Nothing},              {Pipeline, Needs::Nothing}},
    {"printDisplComplete",            {Pipeline, Needs::Nothing},              {Pipeline, Needs::Nothing}},
    {"printMps",                      kNever,                                  {Pipeline, Needs::Nothing}},
    {"printDiv",                      {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"printGal",                      {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"printOverdense",                {Any, Needs::Nothing},                   {Any, Needs::Nothing}},
    {"printRecCat",                   {Pipeline, Needs::PerformRecons},        {Pipeline, Needs::PerformRecons}},
  };

  /// The keys that choose the kind of a run; at most one is given.
  const char* const kInputKeys[] = {"displacement_file", "displacement_mean_file", "divergence_input_file"};

  /// The switches that decide which keys apply, read before the others.
  struct Switches {
    bool lightcone = false;
    Kind kind = Pipeline;
    std::string inputKey;
    bool performRecons = false;
    bool randomFile = false;
    bool divergenceFile = false;
    bool rejectCrossings = false;
    bitvf::SecondPassRandoms secondPass = bitvf::SecondPassRandoms::Resample;
  };

  const Key* find_key (const std::string& name)
  {
    for (const Key& k : kKeys)
      if (name == k.name) return &k;
    return nullptr;
  }

  /// Why a key does not apply to the run, or empty when it does. The first
  /// reason that holds: the geometry, the kind of run, a divergence file,
  /// a switch off.
  std::string key_rule (const std::string& name, const Switches& run)
  {
    const Key* key = find_key(name);
    if (key == nullptr) return "unknown key";
    const Rule& rule = run.lightcone ? key->lightcone : key->box;
    if (rule.kinds == 0) return run.lightcone ? "does not apply to a lightcone" : "does not apply to a box";
    if ((rule.kinds & run.kind) == 0)
      return run.kind == Pipeline ? "does not apply to a run without displacement_file or displacement_mean_file"
                                  : "does not apply to a run with " + run.inputKey;
    switch (rule.needs) {
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
      case Needs::RejectCrossings:
        if (!run.rejectCrossings) return "does not apply to a run without reject_mask_crossings";
        break;
      case Needs::ExternalRandoms:
        if (!run.performRecons) return "does not apply to a run without perform_recons";
        if (run.secondPass != bitvf::SecondPassRandoms::External) return "does not apply to a run without random_recons_mode = external";
        break;
    }
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

    bool degrees = true;
    std::string maskFile, biasFile, secondRandomFile;

    double omegaM = 0., h = 0., w0 = -1., wa = 0.;

    bitvf::BoxConfig config;
    bitvf::LightconeConfig lightConfig;

    std::string outputDir, outputName, extension;
    bool printMean = false, printComplete = false, printDiv = false, printMps = false;
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
    const bool lightcone = switches.lightcone;

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
    if (lightcone) {
      run.randomFile = f.get<std::string>("random_file");
      run.randomColumns = p.three<std::string>("column_random");
      run.degrees = p.word("units", {"degrees", "radians"}) == "degrees";
      if (switches.secondPass == bitvf::SecondPassRandoms::External)
        run.secondRandomFile = f.get<std::string>("random_file_recons");
    }
    else if (pipeline) {
      run.randomFile = p.optional<std::string>("random_file");
      if (run.randomFile) run.randomColumns = p.three<std::string>("column_random");
    }
    if (pipeline || kind == Divergence || lightcone) {
      run.delimiter = f.get<char>("delimiter", run.delimiter);
      run.comment = f.get<char>("comment", run.comment);
    }
    if (lightcone) run.maskFile = f.get<std::string>("veto_mask_filename");
    if (lightcone && pipeline) {
      otswap::RedshiftCut& cut = run.lightConfig.cut;
      const std::optional<double> lower = p.optional<double>("z_min");
      const std::optional<double> upper = p.optional<double>("z_max");
      if (lower) cut.min = *lower;
      if (upper) cut.max = *upper;
      if (lower && upper && *upper < *lower)
        p.fail("z_max", p.written("z_max") + " is below z_min, " + p.written("z_min") + ", at line " +
               p.line("z_min"));
    }
    if (!pipeline) run.inputFile = f.get<std::string>(switches.inputKey);
    if (kind == Complete || kind == Mean) run.divergenceFile = p.optional<std::string>("divergence_file");

    if (pipeline && (lightcone || recons)) {
      run.omegaM = f.get<double>("Omega_matter");
      p.require("Omega_matter", run.omegaM > 0. && run.omegaM <= 1., "in (0, 1]");
      run.h = f.get<double>("little_h");
      p.require("little_h", run.h > 0., "positive");
      run.w0 = f.get<double>("w0", run.w0);
      run.wa = f.get<double>("wa", run.wa);
    }

    bitvf::BoxConfig& c = run.config;
    bitvf::LightconeConfig& lc = run.lightConfig;
    otswap::Config& ot = lightcone ? lc.ot : c.ot;
    if (pipeline) {
      ot.nRealizations = f.get<unsigned>("Nrec", ot.nRealizations);
      p.require("Nrec", ot.nRealizations >= 1, "at least 1");
      ot.convergence = f.get<double>("threshold", ot.convergence);
      p.require("threshold", ot.convergence > 0. && ot.convergence < 1., "strictly between 0 and 1");
      ot.seed = f.get<unsigned>("seed", ot.seed);
      ot.cellSize = f.get<double>("mesh_cellsize", ot.cellSize);
      p.require("mesh_cellsize", ot.cellSize > 0., "positive");
    }
    if (pipeline && lightcone) {
      lc.nBins = f.get<unsigned>("bins");
      p.require("bins", lc.nBins >= 1, "at least 1");
      ot.rejectCrossings = switches.rejectCrossings;
      if (ot.rejectCrossings)
        ot.maxUnobservedPixelsCrossed = f.get<unsigned>("max_unobserved_pixels_crossed",
                                                        ot.maxUnobservedPixelsCrossed);
    }
    if (pipeline && !lightcone) {
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
          if (!std::isfinite(lower[k]))
            p.fail("box_min", "element " + std::to_string(k+1) + ", " + p.written("box_min", k) + ", is not finite");
          if (!std::isfinite(upper[k]))
            p.fail("box_max", "element " + std::to_string(k+1) + ", " + p.written("box_max", k) + ", is not finite");
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

    bitvf::Config& v = lightcone ? lc.voids : c.voids;
    v.smoothing.sigma = f.get<double>("sigma", v.smoothing.sigma);
    p.require("sigma", v.smoothing.sigma >= 0., "non-negative");
    if (f.has("sigma_unit"))
      v.smoothing.unit = p.word("sigma_unit", {"mps", "coordinate"}) == "mps" ? bitvf::LengthUnit::Mps
                                                                             : bitvf::LengthUnit::Coordinate;
    if (kind != Divergence && !switches.divergenceFile) {
      v.divergenceCellSize = p.optional<double>("cell_size");
      if (v.divergenceCellSize) p.require("cell_size", *v.divergenceCellSize > 0., "positive");
    }

    if (recons && lightcone) {
      bitvf::LightconeRsd rsd;
      rsd.zspaceVoids = f.get<bool>("save_zspace_voids", rsd.zspaceVoids);
      rsd.sigma = f.get<double>("rsd_smoothing_scale");
      p.require("rsd_smoothing_scale", rsd.sigma >= 0., "non-negative");
      rsd.correction.weightByRealizations = f.get<bool>("rsd_weight_by_realizations",
                                                        rsd.correction.weightByRealizations);
      run.biasFile = f.get<std::string>("bias_file");
      rsd.correction.verbosity = run.verbosity;
      rsd.secondPassRandoms = switches.secondPass;
      lc.rsd = rsd;
    }
    else if (recons) {
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
      if (lightcone) run.printMps = f.get<bool>("printMps", false);
    }
    run.printDiv = f.get<bool>("printDiv", false);
    run.printGal = f.get<bool>("printGal", false);
    run.printOverdense = f.get<bool>("printOverdense", false);
    if (recons) run.printCatalog = f.get<bool>("printRecCat", false);

    ot.verbosity = run.verbosity;
    v.verbosity = run.verbosity;
    v.keepDivergence = run.printDiv;
    v.keepAssignment = run.printGal;
    const bool keep = run.printMean || run.printComplete || run.printCatalog || run.printMps;
    c.keepReconstructions = keep;
    lc.keepReconstructions = keep;

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

  /// A sky catalogue of the run, in radians; the columns of the randoms.
  std::vector<double> read_sky (const Run& run, const std::string& file)
  {
    std::vector<double> sky = app::call_otswap([&] {
      return otswap::io::read(file, run.randomColumns, run.delimiter, run.comment);
    }).values;
    if (run.degrees) sky = app::call_otswap([&] { return otswap::skyToRadians(std::move(sky)); });
    return sky;
  }

  /// One run from a parameter file. Every key that applies is read before
  /// the data; one key left after the reads is a bug here, reported as an
  /// unknown key. The files are written after the run, the otswap files
  /// first, then the redshift-space pass's, then the final pass's.
  int run_file (const std::string& file)
  {
    const otswap::io::ParameterFile parameters(file);
    const app::Parameters p(parameters, file);

    Switches switches;
    switches.lightcone = p.word("geometry", {"box", "lightcone"}) == "lightcone";
    std::vector<std::string> given;
    for (const char* key : kInputKeys)
      if (parameters.has(key) && !(switches.lightcone && std::string(key) == "divergence_input_file"))
        given.push_back(key);
    std::sort(given.begin(), given.end(), [&parameters] (const std::string& a, const std::string& b) {
      return parameters.line(a) < parameters.line(b);
    });
    if (given.size() > 1)
      p.fail(given[1], "cannot be given with " + given[0] + ", at line " + p.line(given[0]) +
             "; a run reads one of displacement_file, displacement_mean_file" +
             (switches.lightcone ? "" : " and divergence_input_file") + ", or none for a pipeline run");
    if (!given.empty()) {
      switches.inputKey = given[0];
      switches.kind = given[0] == "displacement_file" ? Complete : given[0] == "displacement_mean_file" ? Mean
                                                                                                       : Divergence;
    }
    if (switches.kind == Pipeline) {
      switches.performRecons = parameters.get<bool>("perform_recons", false);
      if (switches.lightcone) {
        switches.rejectCrossings = parameters.get<bool>("reject_mask_crossings", otswap::Config{}.rejectCrossings);
        if (switches.performRecons) {
          const std::string mode = p.word("random_recons_mode", {"none", "resample", "external"}, "resample");
          switches.secondPass = mode == "none" ? bitvf::SecondPassRandoms::None
                              : mode == "resample" ? bitvf::SecondPassRandoms::Resample
                                                   : bitvf::SecondPassRandoms::External;
        }
      }
      else
        switches.randomFile = parameters.has("random_file") &&
                              !parameters.get<std::vector<std::string>>("random_file").empty();
    }
    if ((switches.kind == Complete || switches.kind == Mean) && parameters.has("divergence_file"))
      switches.divergenceFile = !parameters.get<std::vector<std::string>>("divergence_file").empty();

    if (app::report_errors(kProgram, app::irrelevant_keys(parameters, file, [&switches] (const std::string& key) {
      return key_rule(key, switches);
    })))
      return 1;

    const Run run = read_run(p, switches);
    if (app::report_unread(kProgram, parameters, file)) return 1;

    app::create_directory(p, run.outputDir);

    const bool lightcone = switches.lightcone;
    const bool report = run.verbosity != otswap::Verbosity::Silent;
    auto wrote = [report] (const std::string& written) {
      if (report) std::clog << kProgram << ": wrote " << written << std::endl;
    };
    auto otswap_path = [&run] (const std::string& prefix, const std::string& inserted = "") {
      return (std::filesystem::path(run.outputDir) / (prefix + run.outputName + inserted + "." + run.extension))
        .string();
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
        if (lightcone && run.extension == "dat") wrote(path.substr(0, path.size() - 4) + ".mps.dat");
      }
      if (run.printGal) {
        path = bitvf::io::fileName(bitvf::io::Product::Galaxies, run.outputDir, name, zspace);
        bitvf::io::writeGalaxies(path, result);
        wrote(path);
      }
    };
    auto wrote_twins = [&] (const std::string& prefix) {
      if (lightcone && run.extension == "dat") {
        wrote(otswap_path(prefix, ".mps"));
        wrote(otswap_path(prefix, ".distances"));
      }
    };
    auto write_reconstruction = [&] (const otswap::Result& result, const std::string& zspace) {
      if (run.printMean) {
        const std::string path = otswap_path("displacement_mean_" + zspace);
        app::call_otswap([&] { otswap::io::writeDisplacementMean(path, result); });
        wrote(path);
        wrote_twins("displacement_mean_" + zspace);
      }
      if (run.printComplete) {
        const std::string path = otswap_path("displacement_complete_" + zspace);
        app::call_otswap([&] { otswap::io::writeDisplacementComplete(path, result, true); });
        wrote(path);
        wrote_twins("displacement_complete_" + zspace);
      }
      if (run.printMps && !run.printMean && !run.printComplete) {
        const std::string path = otswap_path("MPS_" + zspace);
        app::call_otswap([&] { otswap::io::writeMpsProfile(path, result.mpsProfile); });
        wrote(path);
      }
    };

    std::optional<bitvf::DivergenceGrid> grid;
    const bitvf::Config& config = lightcone ? run.lightConfig.voids : run.config.voids;

    switch (run.switches.kind) {
      case Pipeline: {
        std::vector<double> tracers, randoms;
        std::vector<std::int64_t> ids;
        read_tracers(run, tracers, ids);
        if (lightcone && run.degrees)
          tracers = app::call_otswap([&] { return otswap::skyToRadians(std::move(tracers)); });
        if (run.randomFile)
          randoms = lightcone ? read_sky(run, *run.randomFile)
                              : app::call_otswap([&] {
                                  return otswap::io::read(*run.randomFile, run.randomColumns, run.delimiter,
                                                          run.comment);
                                }).values;

        bitvf::BitVFResult result;
        if (lightcone) {
          const otswap::Mask mask = app::call_otswap([&] { return otswap::Mask(run.maskFile); });
          const otswap::DistanceTable distances = app::call_otswap([&] {
            return otswap::DistanceTable(run.omegaM, run.h, run.w0, run.wa);
          });
          bitvf::LightconeConfig lc = run.lightConfig;
          std::vector<double> second;
          if (lc.rsd) {
            lc.rsd->bias = app::call_otswap([&] {
              return otswap::io::readBiasTable(run.biasFile, run.delimiter, run.comment);
            });
            if (switches.secondPass == bitvf::SecondPassRandoms::External)
              second = read_sky(run, run.secondRandomFile);
          }
          result = bitvf::runLightcone(tracers, randoms, mask, distances, lc, second, ids);
        }
        else if (run.switches.performRecons) {
          const otswap::DistanceTable distances = app::call_otswap([&] {
            return otswap::DistanceTable(run.omegaM, run.h, run.w0, run.wa);
          });
          result = bitvf::runBox(tracers, randoms, distances, run.config, ids);
        }
        else
          result = bitvf::runBox(tracers, randoms, run.config, ids);
        if (result.zspaceReconstruction) write_reconstruction(*result.zspaceReconstruction, "zspace_");
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
      case Complete:
      case Mean: {
        const bool mean = run.switches.kind == Mean;
        const otswap::Result reconstruction = app::call_otswap([&] {
          return mean ? otswap::io::readDisplacementMean(run.inputFile)
                      : otswap::io::readDisplacementComplete(run.inputFile);
        });
        if (mean && reconstruction.geometry == otswap::Geometry::Box && !reconstruction.box)
          throw bitvf::Error(run.inputFile + ": the file holds no box limits; write it with otswap 0.3 or later");
        if (run.divergenceFile) grid = bitvf::io::readDivergence(*run.divergenceFile);
        bitvf::Config c = config;
        if (mean) c.nRealizationsMean = reconstruction.nRealizationsMean;
        if (lightcone) {
          const otswap::Mask mask = app::call_otswap([&] { return otswap::Mask(run.maskFile); });
          write_search(bitvf::findVoidsLightcone(reconstruction, read_sky(run, *run.randomFile), mask, c,
                                                 grid ? &*grid : nullptr), false);
        }
        else
          write_search(bitvf::findVoidsBox(reconstruction, c, grid ? &*grid : nullptr), false);
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
