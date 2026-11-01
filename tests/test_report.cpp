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
 *  @file tests/test_report.cpp
 *
 *  @brief The report of the box entry points: the lines of each level,
 *  their order, what a throwing call leaves, and that the text does not
 *  depend on the number of threads.
 */

#include <cstdint>
#include <functional>
#include <iostream>
#include <locale>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "check.h"
#include "fixtures.h"

using namespace bitvf;
using namespace fixtures;

namespace {

  /// Every time "in <digits>.<2 digits> s" replaced with "in T s".
  std::string masked (const std::string& text)
  {
    return std::regex_replace(text, std::regex("in [0-9]+\\.[0-9]{2} s"), "in T s");
  }

  /// The text written by f to log, masked; a bitvf::Error from f is
  /// caught and sets threw.
  std::string written (const std::function<void (std::ostream*)>& f, bool* threw = nullptr)
  {
    std::ostringstream log;
    if (threw != nullptr) *threw = false;
    try {
      f(&log);
    }
    catch (const Error&) {
      if (threw == nullptr) throw;
      *threw = true;
    }
    return masked(log.str());
  }

  /// The text written to std::clog by f.
  std::string to_clog (const std::function<void ()>& f)
  {
    std::ostringstream sink;
    std::streambuf* old = std::clog.rdbuf(sink.rdbuf());
    f();
    std::clog.rdbuf(old);
    return sink.str();
  }

  Config with (Config c, const otswap::Verbosity verbosity, std::ostream* log)
  {
    c.verbosity = verbosity;
    c.log = log;
    return c;
  }

  /// A stream for expected text, as the report composes it.
  std::ostringstream text ()
  {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    return out;
  }

  std::string box (const std::array<double, 3>& lower, const std::array<double, 3>& upper)
  {
    std::ostringstream out = text();
    for (std::size_t k = 0; k < 3; ++k)
      out << (k ? " x [" : "[") << lower[k] << ", " << upper[k] << "]";
    return out.str();
  }

  std::string basins (const Result& r)
  {
    return std::to_string(r.voids.nBasins) + " voids and " + std::to_string(r.overdense.nBasins) +
      " overdense basins";
  }

  /// The grid line of a search, without its time.
  std::string grid_line (const DivergenceGrid& g, const std::string& source, const std::string& where)
  {
    std::ostringstream out = text();
    out << "bitvf: grid of " << g.nCellsAxis[0] << " x " << g.nCellsAxis[1] << " x " << g.nCellsAxis[2]
        << " cells of side " << g.cellSize << " in coordinate units (" << g.cellSizeMps << " mps, " << source
        << ") on " << box(g.lower, g.upper) << where;
    return out.str();
  }

  /// The reconstruction of the Translated search: its arrays, the limits
  /// as the box, the tracers outside them flagged.
  otswap::Result reconstruction (const Translated& in, const double mps)
  {
    const std::size_t n = Translated::kN + Translated::kOutside;
    otswap::Result rec;
    rec.geometry = otswap::Geometry::Box;
    rec.nObjects = n;
    rec.nRealizations = 2;
    rec.tracers = in.tracers;
    rec.displacement = in.displacement;
    rec.valid.assign(2 * n, 1);
    rec.outsideBox.assign(n, 0);
    for (std::size_t i = Translated::kN; i < n; ++i) rec.outsideBox[i] = 1;
    rec.ids = in.ids;
    rec.mps = mps;
    rec.selection.boxLimits = in.limits;
    rec.box = in.limits;
    return rec;
  }

}

int main ()
{
  set_threads(4);
  const Translated in;
  const Result r = in.run();
  const DivergenceGrid& d = *r.divergence;
  const otswap::Result rec = reconstruction(in, Translated::kMps);
  const LatticeField field;
  Config fieldConfig;
  fieldConfig.smoothing.sigma = 1.5;
  fieldConfig.smoothing.unit = LengthUnit::Coordinate;
  fieldConfig.keepDivergence = true;
  fieldConfig.log = nullptr;
  std::vector<double> ramp(8);
  for (std::size_t c = 0; c < 8; ++c) ramp[c] = (double)c - 0.5;

  const std::string limitsBox = box(in.limits->lower, in.limits->upper);
  const std::string boxLine = "bitvf: kept 3000 of 3004 tracers: 4 outside the box " + limitsBox + "\n";

  const std::vector<std::pair<std::string, std::function<void (otswap::Verbosity, std::ostream*)>>> calls = {
    {"findVoidsBox, arrays", [&] (otswap::Verbosity v, std::ostream* log) {
      findVoidsBox(in.tracers, in.displacement, in.valid, Translated::kMps, with(in.config, v, log), in.limits, in.ids);
    }},
    {"findVoidsBox, reconstruction", [&] (otswap::Verbosity v, std::ostream* log) {
      findVoidsBox(rec, with(in.config, v, log));
    }},
    {"findVoidsBoxFromDivergence", [&] (otswap::Verbosity v, std::ostream* log) {
      findVoidsBoxFromDivergence(d, in.tracers, with(in.config, v, log), in.ids);
    }},
    {"findVoidsBoxFromField", [&] (otswap::Verbosity v, std::ostream* log) {
      findVoidsBoxFromField(field.centres, field.values, with(fieldConfig, v, log), field.weights);
    }},
    {"watershedBox", [&] (otswap::Verbosity v, std::ostream* log) {
      watershedBox(field.values, field.n, v, log);
    }},
  };

  group("Silent writes nothing, and a null log writes nothing at Detailed, not even to std::clog; every entry point");
  for (const auto& call : calls) {
    check(written([&] (std::ostream* log) { call.second(otswap::Verbosity::Silent, log); }).empty(),
          call.first + ", Silent");
    check(to_clog([&] { call.second(otswap::Verbosity::Detailed, nullptr); }).empty(), call.first + ", null log");
  }

  group("Normal: the selection lines and the call line, exactly");
  {
    const std::string call = "bitvf: findVoidsBox: " + basins(r) + " from 3004 tracers in T s; 0 without a valid "
      "displacement\n";
    check(written([&] (std::ostream* log) { calls[0].second(otswap::Verbosity::Normal, log); }) == boxLine + call,
          "arrays, with limits");

    const Result open = findVoidsBox(in.tracers, in.displacement, in.valid, Translated::kMps, in.config,
                                     std::nullopt, in.ids);
    check(written([&] (std::ostream* log) {
            findVoidsBox(in.tracers, in.displacement, in.valid, Translated::kMps,
                         with(in.config, otswap::Verbosity::Normal, log), std::nullopt, in.ids);
          }) == "bitvf: findVoidsBox: " + basins(open) + " from 3004 tracers in T s; 0 without a valid displacement\n",
          "arrays, without limits: the call line alone");

    check(written([&] (std::ostream* log) {
            findVoidsBox(in.tracers, in.displacement, in.valid, Translated::kMps,
                         with(in.config, otswap::Verbosity::Normal, log), in.limits, in.ids, &d);
          }) == boxLine + "bitvf: kept 3000 of 3000 tracers: 0 outside the reused grid\n" + call,
          "arrays, with a reused grid: the grid line also with 0 outside");

    check(written([&] (std::ostream* log) { calls[1].second(otswap::Verbosity::Normal, log); }) == boxLine + call,
          "reconstruction with limits: the box line repeated");

    const std::string fromDivergence = "bitvf: findVoidsBoxFromDivergence: " + basins(r) + " from 3004 tracers in T s\n";
    check(written([&] (std::ostream* log) { calls[2].second(otswap::Verbosity::Normal, log); }) ==
          boxLine + "bitvf: kept 3000 of 3000 tracers: 0 outside the divergence grid\n" + fromDivergence,
          "from a divergence with limits");
    DivergenceGrid unlimited = d;
    unlimited.limits.reset();
    check(written([&] (std::ostream* log) {
            findVoidsBoxFromDivergence(unlimited, in.tracers, with(in.config, otswap::Verbosity::Normal, log), in.ids);
          }) == "bitvf: kept 3000 of 3004 tracers: 4 outside the divergence grid\n" + fromDivergence,
          "from a divergence without limits: the tracers beyond the grid on the grid line");

    const Result f = findVoidsBoxFromField(field.centres, field.values, fieldConfig, field.weights);
    check(written([&] (std::ostream* log) { calls[3].second(otswap::Verbosity::Normal, log); }) ==
          "bitvf: findVoidsBoxFromField: " + basins(f) + " on 12 x 17 x 9 cells in T s\n", "field search");

    const Partition p = watershedBox(field.values, field.n, otswap::Verbosity::Silent);
    check(written([&] (std::ostream* log) { calls[4].second(otswap::Verbosity::Normal, log); }) ==
          "bitvf: watershedBox: " + std::to_string(p.nVoids) + " voids and " + std::to_string(p.nOverdense) +
          " overdense basins on 12 x 17 x 9 cells in T s\n", "watershedBox");

    std::vector<double> positive(8);
    for (std::size_t c = 0; c < 8; ++c) positive[c] = 1. + (double)c;
    check(written([&] (std::ostream* log) { watershedBox(positive, {2, 2, 2}, otswap::Verbosity::Normal, log); }) ==
          "bitvf: watershedBox: 0 voids and 1 overdense basin on 2 x 2 x 2 cells in T s\n",
          "watershedBox, no cell below 0: no error, the call line written");
  }

  group("singular forms: 1 void, 1 overdense basin, 1 tracer, 1 cell, 1 realization");
  {
    Config c = in.config;
    c.keepAssignment = false;
    DivergenceGrid one = d;
    one.nTracers.clear();
    one.nRealizations = 1;
    one.nRealizationsMean = 1;
    const std::vector<double> tracer = {50., 50., 50.};
    const Result f = findVoidsBoxFromDivergence(one, tracer, c);
    const std::string t = written([&] (std::ostream* log) {
      findVoidsBoxFromDivergence(one, tracer, with(c, otswap::Verbosity::Detailed, log));
    });
    check(t.find("bitvf: kept 1 of 1 tracer: 0 outside the box ") != std::string::npos &&
          t.find("bitvf: kept 1 of 1 tracer: 0 outside the divergence grid\n") != std::string::npos &&
          t.find("bitvf: checked 1 tracer in T s\n") != std::string::npos &&
          t.find(", from 1 realization\n") != std::string::npos &&
          t.find("bitvf: placed 1 tracer in the cells in T s\n") != std::string::npos &&
          t.find(" from 1 tracer in T s\n") != std::string::npos, "one tracer, one realization");

    Config sharp;
    sharp.smoothing.sigma = 0.;
    const std::vector<double> centres = {0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 1., 1.,
                                         1., 0., 0., 1., 0., 1., 1., 1., 0., 1., 1., 1.};
    check(written([&] (std::ostream* log) {
            findVoidsBoxFromField(centres, ramp, with(sharp, otswap::Verbosity::Normal, log));
          }) == "bitvf: findVoidsBoxFromField: 1 void and 0 overdense basins on 2 x 2 x 2 cells in T s\n",
          "one void from the field search");
    check(written([&] (std::ostream* log) {
            watershedBox({-1.}, {1, 1, 1}, otswap::Verbosity::Detailed, log);
          }) == "bitvf: 1 cell in voids, 1 of them below 0; 0 in overdense basins, 0 of them at or above 0\n"
                "bitvf: watershedBox: 1 void and 0 overdense basins on 1 x 1 x 1 cells in T s\n", "one cell");
  }

  group("Detailed: the Normal lines, the step lines in the order of the steps, the call line last");
  {
    std::ostringstream expected = text();
    expected << boxLine
             << "bitvf: checked 3004 tracers and " << r.counts.displacements << " valid displacements in T s\n"
             << "bitvf: mean particle separation 7 in coordinate units (given)\n"
             << grid_line(d, "default", " (the box limits) in T s\n")
             << "bitvf: divergence of the valid displacements of 3000 tracers in T s\n"
             << "bitvf: smoothed with sigma 3.5 in coordinate units (0.5 mps) in T s\n"
             << "bitvf: watershed into " << basins(r) << " in T s\n"
             << "bitvf: catalogues of the basins in T s\n"
             << "bitvf: kept the tracer assignment and the divergence grid in T s\n"
             << "bitvf: findVoidsBox: " << basins(r) << " from 3004 tracers in T s; 0 without a valid displacement\n";
    check(written([&] (std::ostream* log) { calls[0].second(otswap::Verbosity::Detailed, log); }) == expected.str(),
          "arrays, with limits, every output kept");

    Config plain = in.config;
    plain.keepAssignment = false;
    plain.keepDivergence = false;
    plain.keepUnsmoothedDivergence = false;
    plain.divergenceCellSize = 1.;
    const std::string t = written([&] (std::ostream* log) {
      findVoidsBox(in.tracers, in.displacement, in.valid, Translated::kMps,
                   with(plain, otswap::Verbosity::Detailed, log), std::nullopt, in.ids);
    });
    check(t.find("mps, given) on ") != std::string::npos &&
          t.find(" (the tracers and their displacement endpoints) in T s\n") != std::string::npos &&
          t.find("bitvf: kept ") == std::string::npos,
          "arrays without limits, cell size given, nothing kept: no selection line, no kept line");

    const otswap::Result other = reconstruction(in, 6.5);
    std::ostringstream reused = text();
    reused << boxLine << "bitvf: kept 3000 of 3000 tracers: 0 outside the reused grid\n"
           << "bitvf: checked 3004 tracers and " << r.counts.displacements << " valid displacements in T s\n"
           << "bitvf: mean particle separation 7 in coordinate units (from the reused grid; the reconstruction's "
              "is 6.5)\n"
           << grid_line(d, "reused", " (the reused grid) in T s\n");
    const std::string u = written([&] (std::ostream* log) {
      findVoidsBox(other, with(in.config, otswap::Verbosity::Detailed, log), &d);
    });
    check(u.compare(0, reused.str().size(), reused.str()) == 0, "reconstruction with a reused grid: both mps");
    check(written([&] (std::ostream* log) { calls[1].second(otswap::Verbosity::Detailed, log); })
            .find("bitvf: mean particle separation 7 in coordinate units (from the reconstruction)\n"
                  + grid_line(d, "default", " (the box limits) in T s\n")) != std::string::npos,
          "reconstruction: its mps, the grid on its box limits");

    std::ostringstream fromDivergence = text();
    fromDivergence << boxLine << "bitvf: kept 3000 of 3000 tracers: 0 outside the divergence grid\n"
                   << "bitvf: checked 3004 tracers in T s\n"
                   << "bitvf: mean particle separation 7 in coordinate units (from the divergence)\n"
                   << grid_line(d, "from the divergence", ", from 2 realizations, each the mean of 3\n")
                   << "bitvf: placed 3000 tracers in the cells in T s\n"
                   << "bitvf: smoothed with sigma 3.5 in coordinate units (0.5 mps) in T s\n"
                   << "bitvf: watershed into " << basins(r) << " in T s\n"
                   << "bitvf: catalogues of the basins in T s\n"
                   << "bitvf: kept the tracer assignment and the divergence grid in T s\n"
                   << "bitvf: findVoidsBoxFromDivergence: " << basins(r) << " from 3004 tracers in T s\n";
    check(written([&] (std::ostream* log) { calls[2].second(otswap::Verbosity::Detailed, log); }) ==
          fromDivergence.str(), "from a divergence");

    const Result f = findVoidsBoxFromField(field.centres, field.values, fieldConfig, field.weights);
    const DivergenceGrid& g = *f.divergence;
    std::ostringstream fromField = text();
    fromField << "bitvf: lattice inferred: 12 x 17 x 9 cells of side " << g.cellSize << " in coordinate units on "
              << box(g.lower, g.upper) << " in T s\n"
              << "bitvf: no mean particle separation given\n"
              << "bitvf: smoothed with sigma 1.5 in coordinate units in T s\n"
              << "bitvf: watershed into " << basins(f) << " in T s\n"
              << "bitvf: catalogues of the basins in T s\n"
              << "bitvf: kept the divergence grid in T s\n"
              << "bitvf: findVoidsBoxFromField: " << basins(f) << " on 12 x 17 x 9 cells in T s\n";
    check(written([&] (std::ostream* log) { calls[3].second(otswap::Verbosity::Detailed, log); }) ==
          fromField.str(), "field search, no mps");

    Config bare = fieldConfig;
    bare.smoothing.sigma = 0.;
    bare.keepDivergence = false;
    const std::string b = written([&] (std::ostream* log) {
      findVoidsBoxFromField(field.centres, field.values, with(bare, otswap::Verbosity::Detailed, log), field.weights,
                            2.);
    });
    check(b.find("bitvf: mean particle separation 2 in coordinate units (given)\n"
                 "bitvf: not smoothed (sigma 0)\nbitvf: watershed into ") != std::string::npos &&
          b.find("bitvf: kept") == std::string::npos, "field search with mps, sigma 0: no time, no kept line");

    const Partition p = watershedBox(field.values, field.n, otswap::Verbosity::Silent);
    std::uint64_t sums[4] = {0, 0, 0, 0};
    for (std::size_t k = 0; k < p.nVoids; ++k) { sums[0] += p.voidCells[k]; sums[1] += p.voidCoreCells[k]; }
    for (std::size_t k = 0; k < p.nOverdense; ++k) {
      sums[2] += p.overdenseCells[k];
      sums[3] += p.overdenseCoreCells[k];
    }
    check(written([&] (std::ostream* log) { calls[4].second(otswap::Verbosity::Detailed, log); }) ==
          "bitvf: " + std::to_string(sums[0]) + " cells in voids, " + std::to_string(sums[1]) + " of them below 0; " +
          std::to_string(sums[2]) + " in overdense basins, " + std::to_string(sums[3]) + " of them at or above 0\n" +
          "bitvf: watershedBox: " + std::to_string(p.nVoids) + " voids and " + std::to_string(p.nOverdense) +
          " overdense basins on 12 x 17 x 9 cells in T s\n", "watershedBox: the cell counts, then the call line");
  }

  group("<W>: the tracers read without a valid displacement, the tracers read minus tracersUsed");
  {
    std::vector<std::uint8_t> valid(2 * (Translated::kN + Translated::kOutside), 1);
    for (std::size_t i = 0; i < 5; ++i)
      for (std::size_t s = 0; s < 2; ++s) valid[s * (Translated::kN + Translated::kOutside) + 100 * i] = 0;
    valid[7] = 0;
    const Result w = findVoidsBox(in.tracers, in.displacement, valid, Translated::kMps, in.config, in.limits, in.ids);
    const std::size_t read = w.counts.tracers - w.counts.tracersOutsideBox - w.counts.tracersOutsideGrid;
    check(read - w.counts.tracersUsed == 5, "5 tracers with every realization invalid; one with one left is used");
    const std::string t = written([&] (std::ostream* log) {
      findVoidsBox(in.tracers, in.displacement, valid, Translated::kMps,
                   with(in.config, otswap::Verbosity::Normal, log), in.limits, in.ids);
    });
    check(t.find("from 3004 tracers in T s; 5 without a valid displacement\n") != std::string::npos,
          "the call line counts 5");
  }

  group("a throwing call: a check error writes nothing; no void found leaves the step lines up to the watershed, "
        "no call line; at Normal only the selection lines");
  {
    bool threw = false;
    const std::string t = written([&] (std::ostream* log) {
      findVoidsBox(in.tracers, in.displacement, in.valid, -1., with(in.config, otswap::Verbosity::Detailed, log),
                   in.limits, in.ids);
    }, &threw);
    check(threw && t.empty(), "check error, Detailed: nothing");

    DivergenceGrid flat = d;
    flat.divergenceUnsmoothed.assign(flat.nCells, 1.);
    flat.nTracers.clear();
    const std::string fd = written([&] (std::ostream* log) {
      findVoidsBoxFromDivergence(flat, in.tracers, with(in.config, otswap::Verbosity::Detailed, log), in.ids);
    }, &threw);
    const std::size_t watershed = fd.find("bitvf: watershed into 0 voids and ");
    check(threw && watershed != std::string::npos && fd.find('\n', watershed) + 1 == fd.size() &&
          fd.find("bitvf: placed 3000 tracers in the cells in T s\nbitvf: smoothed with sigma ") != std::string::npos,
          "from a divergence, Detailed: the lines end with the watershed's");
    check(written([&] (std::ostream* log) {
            findVoidsBoxFromDivergence(flat, in.tracers, with(in.config, otswap::Verbosity::Normal, log), in.ids);
          }, &threw) == boxLine + "bitvf: kept 3000 of 3000 tracers: 0 outside the divergence grid\n" && threw,
          "from a divergence, Normal: the selection lines only");

    std::vector<double> up(field.values.size(), 1.);
    const std::string ff = written([&] (std::ostream* log) {
      findVoidsBoxFromField(field.centres, up, with(fieldConfig, otswap::Verbosity::Detailed, log), field.weights);
    }, &threw);
    const std::size_t fw = ff.find("bitvf: watershed into 0 voids and ");
    check(threw && ff.rfind("bitvf: lattice inferred: ", 0) == 0 && fw != std::string::npos &&
          ff.find('\n', fw) + 1 == ff.size(), "field search, Detailed: the lines end with the watershed's");
    check(written([&] (std::ostream* log) {
            findVoidsBoxFromField(field.centres, up, with(fieldConfig, otswap::Verbosity::Normal, log), field.weights);
          }, &threw).empty() && threw, "field search, Normal: nothing");
  }

  group("the masked text is equal with 1 and 4 threads, at Detailed, on every entry point");
  for (const auto& call : calls) {
    set_threads(1);
    const std::string one = written([&] (std::ostream* log) { call.second(otswap::Verbosity::Detailed, log); });
    set_threads(4);
    const std::string four = written([&] (std::ostream* log) { call.second(otswap::Verbosity::Detailed, log); });
    check(!one.empty() && one == four, call.first);
  }

  return report("test_report");
}
