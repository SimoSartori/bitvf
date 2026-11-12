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
 *  @file tests/test_readers.cpp
 *
 *  @brief The readers of bitvf::io: the divergence file and otswap's mean
 *  displacement file read back, reuse, the default cell size, errors.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
#include <iomanip>
#include <limits>
#include <locale>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "otswap/Maths.h"
#include "otswap/OT.h"
#include "otswap/io.h"

#include "bitvf/io.h"
#include "check.h"
#include "fixtures.h"

using namespace bitvf;

namespace {

  constexpr std::size_t kN = 200;
  constexpr std::size_t kOut = 4;
  constexpr std::size_t kR = 2;
  constexpr double kMps = 10.;

  const char* const kFormats[2] = {".fits", ".txt"};

  otswap::BoxLimits box (const double lo, const double hi)
  {
    otswap::BoxLimits limits;
    limits.lower = {lo, lo, lo};
    limits.upper = {hi, hi, hi};
    return limits;
  }

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "bitvf_test_readers_" + name;
  }

  bool same_bits (const double a, const double b)
  {
    if (std::isnan(a) && std::isnan(b)) return true;
    return std::memcmp(&a, &b, sizeof(double)) == 0;
  }

  bool same_bits (const std::vector<double>& a, const std::vector<double>& b)
  {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
      if (!same_bits(a[i], b[i])) return false;
    return true;
  }

  bool same_limits (const std::optional<otswap::BoxLimits>& a, const std::optional<otswap::BoxLimits>& b)
  {
    if (!a || !b) return !a && !b;
    for (int c = 0; c < 3; ++c)
      if (!same_bits(a->lower[c], b->lower[c]) || !same_bits(a->upper[c], b->upper[c])) return false;
    return true;
  }

  /// The double an ASCII file holds for @p value written with @p digits
  /// significant digits.
  double as_written (const double value, const int digits)
  {
    if (std::isnan(value)) return value;
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::defaultfloat << std::setprecision(digits) << value;
    return std::strtod(out.str().c_str(), nullptr);
  }

  std::vector<double> as_written (const std::vector<double>& values, const int digits)
  {
    std::vector<double> w;
    for (const double v : values) w.push_back(as_written(v, digits));
    return w;
  }

  /// A search on tracers uniform in [0, 100)^3, kOut more beyond the
  /// limits [-2, 102]^3, displaced at random and towards the centre, so
  /// that both kinds of basin appear on cells of cbrt(10) mps; divergence
  /// kept, three realizations averaged.
  struct Input {
    std::vector<double> tracers, displacement;
    std::vector<std::uint8_t> valid;
    Config config;
    std::optional<otswap::BoxLimits> limits = box(-2., 102.);
    std::vector<std::int64_t> ids;
    const DivergenceGrid* grid = nullptr;

    explicit Input (const std::size_t nRealizations = kR)
    {
      std::mt19937 rng(2026);
      auto uniform = [&rng] (const double lo, const double hi) {
        return lo + (hi - lo) * (rng() / 4294967296.);
      };
      const std::size_t n = kN + kOut;
      for (std::size_t i = 0; i < 3 * kN; ++i) tracers.push_back(uniform(0., 100.));
      for (std::size_t i = 0; i < 3 * kOut; ++i) tracers.push_back(150. + i);
      for (std::size_t r = 0; r < nRealizations; ++r)
        for (std::size_t i = 0; i < 3 * n; ++i)
          displacement.push_back(uniform(-2., 2.) - (tracers[i] - 50.) / 50.);
      config.log = nullptr;
      config.keepDivergence = true;
      config.nRealizationsMean = 3;
      config.divergenceCellSize = otswap::maths::cbrt(10.);
    }

    Result run () const
    {
      return findVoidsBox(tracers, displacement, valid, kMps, config, limits, ids, grid);
    }
  };

  /// Every field of two grids equal, bit for bit.
  bool same_grid (const DivergenceGrid& a, const DivergenceGrid& b, const std::string& what)
  {
    bool equal = true;
    auto field = [&] (const bool same, const std::string& name) {
      check(same, what + ": " + name);
      equal = equal && same;
    };
    field(a.geometry == b.geometry && a.nCells == b.nCells && a.nCellsAxis == b.nCellsAxis, "geometry and counts");
    field(same_bits(a.cellSize, b.cellSize) && same_bits(a.cellSizeMps, b.cellSizeMps) && same_bits(a.mps, b.mps),
          "cellSize, cellSizeMps, mps");
    bool bounds = true;
    for (int c = 0; c < 3; ++c) bounds = bounds && same_bits(a.lower[c], b.lower[c]) && same_bits(a.upper[c], b.upper[c]);
    field(bounds, "lower and upper");
    field(same_limits(a.limits, b.limits), "limits");
    field(a.nRealizations == b.nRealizations && a.nRealizationsMean == b.nRealizationsMean,
          "nRealizations and nRealizationsMean");
    field(same_bits(a.distanceMin, b.distanceMin) && same_bits(a.distanceMax, b.distanceMax) &&
          same_bits(a.raMin, b.raMin) && same_bits(a.raMax, b.raMax) && same_bits(a.decMin, b.decMin) &&
          same_bits(a.decMax, b.decMax) && a.mpsProfile.redshift.empty() && b.mpsProfile.redshift.empty() &&
          a.mpsProfile.mps.empty() && b.mpsProfile.mps.empty() && a.mpsProfile.count.empty() &&
          b.mpsProfile.count.empty() && a.cellRaMin.empty() && b.cellRaMin.empty() && a.cellRaMax.empty() &&
          b.cellRaMax.empty() && a.cellDecMin.empty() && b.cellDecMin.empty() && a.cellDecMax.empty() &&
          b.cellDecMax.empty() && a.cellDistanceMin.empty() && b.cellDistanceMin.empty() &&
          a.cellDistanceMax.empty() && b.cellDistanceMax.empty() && a.centreSky.empty() && b.centreSky.empty() &&
          a.masked.empty() && b.masked.empty(), "the lightcone fields NaN or empty");
    field(same_bits(a.centre, b.centre), "centre");
    field(same_bits(a.divergence, b.divergence), "divergence");
    field(same_bits(a.divergenceSmoothed, b.divergenceSmoothed), "divergenceSmoothed");
    field(same_bits(a.volume, b.volume), "volume");
    field(same_bits(a.weight, b.weight), "weight");
    field(same_bits(a.weightSmoothed, b.weightSmoothed), "weightSmoothed");
    field(a.nDisplacements == b.nDisplacements && a.nTracers == b.nTracers, "nDisplacements and nTracers");
    field(a.voidIndex == b.voidIndex && a.overdenseIndex == b.overdenseIndex, "voidIndex and overdenseIndex");
    return equal;
  }

  bool same_catalogue (const BasinCatalog& a, const BasinCatalog& b)
  {
    return a.nBasins == b.nBasins && same_bits(a.centre, b.centre) && same_bits(a.centreSky, b.centreSky) &&
           same_bits(a.radius, b.radius) && same_bits(a.radiusExt, b.radiusExt) &&
           same_bits(a.deltaAv, b.deltaAv) && same_bits(a.delta2Mps, b.delta2Mps) && a.edgeFlag == b.edgeFlag &&
           same_bits(a.contamination, b.contamination) && a.nTracers == b.nTracers;
  }

  /// writeDivergence, then readDivergence, in @p format.
  DivergenceGrid divergence_read_back (const Result& r, const std::string& format)
  {
    const std::string file = temporary(std::string("divergence") + format);
    io::writeDivergence(file, r);
    DivergenceGrid g = io::readDivergence(file);
    std::remove(file.c_str());
    return g;
  }

  otswap::io::Column real_column (const std::string& name, std::vector<double> data)
  {
    otswap::io::Column c;
    c.name = name;
    c.type = 'D';
    c.data = std::move(data);
    return c;
  }

  otswap::io::Column flag_column (const std::string& name, std::vector<double> data)
  {
    otswap::io::Column c;
    c.name = name;
    c.type = 'J';
    c.data = std::move(data);
    return c;
  }

  otswap::io::Column integer_column (const std::string& name, std::vector<std::int64_t> integers)
  {
    otswap::io::Column c;
    c.name = name;
    c.type = 'K';
    c.integers = std::move(integers);
    return c;
  }

  /// A table and its keywords, written by otswap::io::write.
  struct Table {
    std::vector<otswap::io::Column> columns;
    std::vector<otswap::io::Keyword> keywords;

    otswap::io::Column& column (const std::string& name)
    {
      for (otswap::io::Column& c : columns)
        if (c.name == name) return c;
      throw std::runtime_error("no column " + name);
    }

    void erase_column (const std::string& name)
    {
      for (std::size_t i = 0; i < columns.size(); ++i)
        if (columns[i].name == name) columns.erase(columns.begin() + (std::ptrdiff_t)i);
    }

    void set (const std::string& name, const std::string& value)
    {
      for (otswap::io::Keyword& k : keywords)
        if (k.name == name) {
          k.value = value;
          return;
        }
      keywords.push_back({name, value, ""});
    }

    void erase (const std::string& name)
    {
      for (std::size_t i = 0; i < keywords.size(); ++i)
        if (keywords[i].name == name) keywords.erase(keywords.begin() + (std::ptrdiff_t)i);
    }

    void write (const std::string& file) const
    {
      otswap::io::WriteOptions options;
      options.precision = 17;
      options.keywords = keywords;
      otswap::io::write(file, columns, options);
    }
  };

  /// A small mean displacement file: three tracers, NREC 2, MPS 10.
  Table mean_table ()
  {
    Table t;
    t.columns = {real_column("tracX", {1., 2., 3.}), real_column("tracY", {4., 5., 6.}),
                 real_column("tracZ", {7., 8., 9.}), real_column("displX", {0.1, 0.2, 0.3}),
                 real_column("displY", {0.4, 0.5, 0.6}), real_column("displZ", {0.7, 0.8, 0.9}),
                 flag_column("nValidRec", {2., 1., 0.})};
    t.keywords = {{"PRODUCT", "displacements", ""}, {"GEOMETRY", "box", ""}, {"NREC", "2", ""},
                  {"MPS", "10", ""}};
    return t;
  }

  /// The divergence file of @p r, as writeDivergence writes it.
  Table divergence_table (const Result& r)
  {
    const DivergenceGrid& g = *r.divergence;
    const std::string ascii = temporary("keywords.txt");
    io::writeDivergence(ascii, r);
    Table t;
    t.keywords = otswap::io::readKeywords(ascii);
    std::remove(ascii.c_str());
    std::vector<std::int64_t> index;
    std::vector<double> x, y, z;
    for (std::size_t c = 0; c < g.nCells; ++c) {
      index.push_back((std::int64_t)c);
      x.push_back(g.centre[3*c]);
      y.push_back(g.centre[3*c+1]);
      z.push_back(g.centre[3*c+2]);
    }
    std::vector<std::int64_t> nDispl(g.nDisplacements.begin(), g.nDisplacements.end());
    std::vector<std::int64_t> nTracers(g.nTracers.begin(), g.nTracers.end());
    t.columns = {integer_column("index", index), real_column("cellX", x), real_column("cellY", y),
                 real_column("cellZ", z), real_column("divergence", g.divergence)};
    if (!g.divergenceSmoothed.empty()) t.columns.push_back(real_column("divergenceSmoothed", g.divergenceSmoothed));
    t.columns.push_back(real_column("volume", g.volume));
    t.columns.push_back(real_column("weight", g.weight));
    if (!g.weightSmoothed.empty()) t.columns.push_back(real_column("weightSmoothed", g.weightSmoothed));
    t.columns.push_back(integer_column("nDispl", nDispl));
    t.columns.push_back(integer_column("nTracers", nTracers));
    t.columns.push_back(integer_column("voidIndex", g.voidIndex));
    t.columns.push_back(integer_column("overdenseIndex", g.overdenseIndex));
    return t;
  }

  /// @p read must raise a bitvf::Error naming @p message on @p table,
  /// written in both formats.
  template <typename R>
  void refused (const Table& table, R&& read, const std::string& message, const std::string& what)
  {
    for (const char* format : kFormats) {
      const std::string file = temporary(std::string("refused") + format);
      table.write(file);
      throws_naming([&] { read(file); }, message, what + " " + format);
      std::remove(file.c_str());
    }
  }

  void read_mean (const std::string& f) { io::readDisplacementMean(f); }
  void read_divergence (const std::string& f) { io::readDivergence(f); }

  /// @p f must raise a bitvf::Error starting with "otswap: ", an
  /// otswap::Error nested.
  template <typename F>
  void otswap_nested (F&& f, const std::string& what)
  {
    bool nested = false;
    try {
      f();
    }
    catch (const Error& e) {
      const bool prefixed = std::string(e.what()).compare(0, 8, "otswap: ") == 0;
      try {
        std::rethrow_if_nested(e);
      }
      catch (const otswap::Error&) {
        nested = prefixed;
      }
      catch (...) {
      }
    }
    catch (...) {
    }
    check(nested, what + ": \"otswap: \", otswap's error nested");
  }

}

int main ()
{
  const Input in;
  const Result r = in.run();
  check(r.voids.nBasins > 0 && r.overdense.nBasins > 0 && r.limits && r.divergence &&
        !r.divergence->divergenceSmoothed.empty() && !r.divergence->weightSmoothed.empty() &&
        r.divergence->nRealizations == kR && r.divergence->nRealizationsMean == 3,
        "the fixture: basins of both kinds, limits, smoothed divergence and weight");

  Input plain = in;
  plain.config.smoothing.sigma = 0.;
  plain.limits.reset();
  plain.config.nRealizationsMean = 1;
  const Result p = plain.run();
  check(!p.limits && p.divergence->divergenceSmoothed.empty() && p.divergence->weightSmoothed.empty(),
        "the fixture: no limits, sigma 0, no smoothed divergence or weight");

  group("divergence round trip, FITS and ASCII: every field bit for bit, with and without the optional parts");
  for (const char* format : kFormats) {
    same_grid(divergence_read_back(r, format), *r.divergence, std::string("limits and smoothed ") + format);
    same_grid(divergence_read_back(p, format), *p.divergence, std::string("neither ") + format);
  }

  group("reuse: the grid read back, passed as grid, gives the catalogues of the search that wrote it");
  for (const Result* written : {&r, &p})
    for (const char* format : kFormats) {
      const DivergenceGrid g = divergence_read_back(*written, format);
      Input reuse = written == &r ? in : plain;
      reuse.grid = &g;
      reuse.config.divergenceCellSize.reset();
      const Result s = reuse.run();
      check(same_catalogue(s.voids, written->voids) && same_catalogue(s.overdense, written->overdense) &&
            same_bits(s.divergence->divergence, written->divergence->divergence) &&
            same_bits(s.divergence->divergenceSmoothed, written->divergence->divergenceSmoothed) &&
            s.divergence->voidIndex == written->divergence->voidIndex,
            std::string(written == &r ? "limits " : "no limits ") + format + ": catalogues and grid equal");
    }

  group("mean file, FITS: tracers, displacements, validity, mps, limits, ids, NREC; the search on it uses R = 1");
  {
    std::vector<double> tracers(in.tracers.begin(), in.tracers.begin() + 3 * kN);
    otswap::Config ot;
    ot.nRealizations = 3;
    ot.seed = 1;
    ot.log = nullptr;
    for (const bool limited : {true, false})
      for (const bool named : {true, false}) {
        const std::string what = std::string(limited ? "limits" : "no limits") + (named ? ", ids" : ", no ids");
        const std::vector<double> all = limited ? in.tracers : tracers;
        otswap::Result rec = otswap::reconstructBox(all, {}, ot, limited ? std::optional<otswap::BoxLimits>(box(-2., 102.))
                                                                          : std::nullopt);
        if (named)
          for (std::size_t i = 0; i < rec.nObjects; ++i) rec.ids.push_back(i % 2 == 0 ? (std::int64_t)i : -(std::int64_t)i);
        for (const char* format : kFormats) {
          const bool fits = std::string(format) == ".fits";
          const std::string file = temporary(std::string("mean") + format);
          otswap::io::writeDisplacementMean(file, rec);
          const io::DisplacementMean m = io::readDisplacementMean(file);
          std::remove(file.c_str());

          std::vector<std::uint8_t> valid;
          for (const unsigned v : rec.validRealizations) valid.push_back(v > 0 ? 1 : 0);
          const std::string tag = what + " " + format;
          if (fits) {
            check(same_bits(m.tracers, rec.tracers) && same_bits(m.displacement, rec.meanDisplacement),
                  tag + ": tracers and displacement bit for bit");
          }
          else {
            check(same_bits(m.tracers, as_written(rec.tracers, 9)) &&
                  same_bits(m.displacement, as_written(rec.meanDisplacement, 9)),
                  tag + ": tracers and displacement as the 9-digit text");
          }
          check(m.valid == valid, tag + ": valid where nValidRec > 0");
          check(same_bits(m.mps, rec.mps), tag + ": mps");
          check(limited ? same_limits(m.limits, rec.box) && m.limits.has_value() : !m.limits, tag + ": limits");
          check(m.ids == rec.ids, tag + ": ids");
          check(m.nRealizationsAveraged == 3, tag + ": nRealizationsAveraged");
          if (limited) {
            bool invalidOutside = true;
            for (std::size_t i = kN; i < kN + kOut; ++i) invalidOutside = invalidOutside && m.valid[i] == 0;
            check(invalidOutside, tag + ": the tracers outside the limits not valid");
          }

          if (!fits) continue;
          Config config;
          config.log = nullptr;
          config.nRealizationsMean = m.nRealizationsAveraged;
          const Result s = findVoidsBox(m.tracers, m.displacement, m.valid, m.mps, config, m.limits, m.ids);
          check(s.nRealizations == 1 && s.config.divergenceCellSize == 1. && s.voids.nBasins > 0 &&
                same_limits(s.limits, m.limits) && s.ids == m.ids,
                tag + ": the search runs with one realization, default cell 1");
        }
      }
  }

  group("default cell size: 1 up to 20 realizations, cbrt(20 / R) above, recorded in Result::config");
  for (const std::size_t nRec : {1, 2, 19, 20, 21, 30}) {
    Input many(nRec);
    many.config.divergenceCellSize.reset();
    many.config.keepDivergence = false;
    const Result s = many.run();
    const double expected = nRec <= 20 ? 1. : otswap::maths::cbrt(20. / (double)nRec);
    check(s.nRealizations == nRec && s.config.divergenceCellSize == expected &&
          (nRec <= 20 || expected < 1.), "R = " + std::to_string(nRec));
  }

  group("errors RD1-RD10, on files written with otswap::io::write; complete and lightcone files named");
  {
    for (const char* format : kFormats) {
      const std::string missing = std::string("/nonexistent_bitvf_directory/x") + format;
      otswap_nested([&] { io::readDisplacementMean(missing); }, std::string("RD1, mean, no file ") + format);
      otswap_nested([&] { io::readDivergence(missing); }, std::string("RD1, divergence, no file ") + format);
    }
    {
      Table t = mean_table();
      t.erase_column("displZ");
      for (const char* format : kFormats) {
        const std::string file = temporary(std::string("rd1") + format);
        t.write(file);
        otswap_nested([&] { io::readDisplacementMean(file); }, std::string("RD1, a missing column ") + format);
        std::remove(file.c_str());
      }
      t = mean_table();
      t.column("nValidRec").type = 'D';
      t.column("nValidRec").data[0] = 1.5;
      for (const char* format : kFormats) {
        const std::string file = temporary(std::string("rd1") + format);
        t.write(file);
        otswap_nested([&] { io::readDisplacementMean(file); }, std::string("RD1, not an integer ") + format);
        std::remove(file.c_str());
      }
      Table d = divergence_table(p);
      d.erase_column("weight");
      for (const char* format : kFormats) {
        const std::string file = temporary(std::string("rd1") + format);
        d.write(file);
        otswap_nested([&] { io::readDivergence(file); }, std::string("RD1, divergence, a missing column ") + format);
        std::remove(file.c_str());
      }
    }

    const Table mean = mean_table();
    const Table div = divergence_table(r);
    {
      Table t = mean;
      t.set("PRODUCT", "voids");
      refused(t, read_mean, " is not a mean displacement file: its PRODUCT is 'voids', where 'displacements' is "
              "expected", "RD2, mean");
      t.set("PRODUCT", "complete displacements");
      refused(t, read_mean, "; read it with otswap::io::readDisplacementComplete and pass the otswap::Result to "
              "findVoidsBox", "RD2, mean, a complete file");
      t.erase("PRODUCT");
      refused(t, read_mean, " has no PRODUCT keyword: not an otswap mean displacement file", "RD2, mean, none");

      Table d = div;
      d.set("PRODUCT", "displacements");
      refused(d, read_divergence, " is not a divergence file: its PRODUCT is 'displacements', where 'divergence' "
              "is expected", "RD2, divergence");
      d.erase("PRODUCT");
      refused(d, read_divergence, " has no PRODUCT keyword: not a bitvf divergence file", "RD2, divergence, none");
    }
    {
      otswap::Config ot;
      ot.nRealizations = 2;
      ot.seed = 1;
      ot.log = nullptr;
      const otswap::Result rec = otswap::reconstructBox(std::vector<double>(in.tracers.begin(),
                                                                            in.tracers.begin() + 3 * kN), {}, ot);
      for (const char* format : kFormats) {
        const std::string file = temporary(std::string("complete") + format);
        otswap::io::writeDisplacementComplete(file, rec);
        throws_naming([&] { io::readDisplacementMean(file); },
                      "its PRODUCT is 'complete displacements', where 'displacements' is expected; read it with "
                      "otswap::io::readDisplacementComplete", std::string("RD2, a real complete file ") + format);
        throws_naming([&] { io::readDivergence(file); },
                      "is not a divergence file: its PRODUCT is 'complete displacements'",
                      std::string("RD2, a real complete file as a divergence ") + format);
        std::remove(file.c_str());
      }
    }
    {
      Table t = mean;
      t.set("GEOMETRY", "lightcone");
      refused(t, read_mean, "reading a lightcone mean displacement file is not implemented", "RD3, mean");
      t.set("GEOMETRY", "sphere");
      refused(t, read_mean, ": GEOMETRY is 'sphere'; it must be box or lightcone", "RD4, mean");
      Table d = div;
      d.set("GEOMETRY", "lightcone");
      refused(d, read_divergence, "reading a lightcone divergence file is not implemented", "RD3, divergence");
      d.set("GEOMETRY", "Box");
      refused(d, read_divergence, ": GEOMETRY is 'Box'; it must be box or lightcone", "RD4, divergence");
    }
    {
      Table t = mean;
      t.erase("MPS");
      refused(t, read_mean, ": keyword MPS is missing", "RD5, MPS missing");
      t = mean;
      t.set("MPS", "ten");
      refused(t, read_mean, ": keyword MPS is 'ten', not a number", "RD5, MPS not a number");
      t.set("MPS", "10 m");
      refused(t, read_mean, ": keyword MPS is '10 m', not a number", "RD5, MPS with trailing text");
      t = mean;
      t.set("NREC", "1.5");
      refused(t, read_mean, ": keyword NREC is '1.5', not a whole number in [0, 4294967295]", "RD5, NREC 1.5");
      t.set("NREC", "-1");
      refused(t, read_mean, ": keyword NREC is '-1', not a whole number", "RD5, NREC negative");
      t.set("NREC", "4294967296");
      refused(t, read_mean, ": keyword NREC is '4294967296', not a whole number", "RD5, NREC too large");

      for (const char* key : {"NCELLS", "NRECMEAN", "CELLSIZE", "CELLSIDE", "LOWERY", "UPPERZ", "NCELLX"}) {
        Table d = div;
        d.erase(key);
        refused(d, read_divergence, std::string(": keyword ") + key + " is missing",
                std::string("RD5, divergence ") + key + " missing");
      }
      Table d = div;
      d.set("NCELLY", "x");
      refused(d, read_divergence, ": keyword NCELLY is 'x', not a whole number", "RD5, NCELLY");
    }
    {
      Table t = mean;
      t.set("BOXLIM", "2");
      refused(t, read_mean, ": keyword BOXLIM is '2', where only 1 is written", "RD6, mean BOXLIM 2");
      t.set("BOXLIM", "1");
      const char* const boxKeys[6] = {"BOXXLO", "BOXXHI", "BOXYLO", "BOXYHI", "BOXZLO", "BOXZHI"};
      for (int k = 0; k < 5; ++k) t.set(boxKeys[k], k % 2 == 0 ? "-1" : "11");
      refused(t, read_mean, ": BOXLIM = 1 records box limits, but keyword BOXZHI is missing", "RD6, mean BOXZHI");
      t.set("BOXZHI", "11");
      t.set("BOXYLO", "abc");
      refused(t, read_mean, ": keyword BOXYLO is 'abc', not a number", "RD6, mean BOXYLO not a number");

      Table loose = mean;
      for (int k = 0; k < 6; ++k) loose.set(boxKeys[k], k % 2 == 0 ? "-1" : "11");
      for (const char* format : kFormats) {
        const std::string file = temporary(std::string("loose") + format);
        loose.write(file);
        const io::DisplacementMean m = io::readDisplacementMean(file);
        check(!m.limits, std::string("RD6: mean BOX* without BOXLIM ignored ") + format);
        std::remove(file.c_str());
      }

      Table d = div;
      d.set("BOXLIM", "0");
      refused(d, read_divergence, ": keyword BOXLIM is '0', where only 1 is written", "RD6, divergence BOXLIM 0");
      d = div;
      d.erase("BOXXHI");
      refused(d, read_divergence, ": BOXLIM = 1 records box limits, but keyword BOXXHI is missing",
              "RD6, divergence BOXXHI");
      d = div;
      d.erase("BOXLIM");
      refused(d, read_divergence, ": keyword BOXXLO is present without BOXLIM = 1", "RD6, divergence no BOXLIM");
    }
    {
      Table t = mean;
      t.set("NREC", "0");
      refused(t, read_mean, ": NREC is 0; it must be at least 1", "RD7, mean");
      Table d = div;
      d.set("NREC", "0");
      refused(d, read_divergence, ": NREC is 0; it must be at least 1", "RD7, divergence NREC");
      d = div;
      d.set("NRECMEAN", "0");
      refused(d, read_divergence, ": NRECMEAN is 0; it must be at least 1", "RD7, divergence NRECMEAN");
    }
    {
      Table d = div;
      for (otswap::io::Column& c : d.columns) {
        c.data.pop_back();
        if (!c.integers.empty()) c.integers.pop_back();
      }
      refused(d, read_divergence, ": it holds " + std::to_string(r.divergence->nCells - 1) + " rows, but NCELLS is " +
              std::to_string(r.divergence->nCells), "RD8, rows");
      d = div;
      std::swap(d.column("index").integers[3], d.column("index").integers[4]);
      refused(d, read_divergence, ": data row 4 has index 4; the rows must hold the cells 0 to NCELLS - 1 in order",
              "RD8, index");
      d = div;
      d.column("nTracers").integers[5] = -1;
      refused(d, read_divergence, ": data row 6 has a negative count (nDispl ", "RD8, nTracers negative");
      d = div;
      d.column("nDispl").integers[0] = -3;
      refused(d, read_divergence, ": data row 1 has a negative count (nDispl -3", "RD8, nDispl negative");
      d = div;
      d.column("overdenseIndex").integers[7] = -2;
      refused(d, read_divergence, ": data row 8 has a basin index below -1 (voidIndex ", "RD8, overdenseIndex");
      d = div;
      d.column("voidIndex").integers[2] = -5;
      refused(d, read_divergence, ": data row 3 has a basin index below -1 (voidIndex -5", "RD8, voidIndex");
    }
    {
      std::ostringstream side;
      side.imbue(std::locale::classic());
      side << std::setprecision(17) << std::nextafter(r.divergence->cellSize, 100.);
      Table d = div;
      d.set("CELLSIDE", side.str());
      refused(d, read_divergence, ": the grid to reuse has a cell side of", "RD9, G9");
      d = div;
      d.set("NCELLS", std::to_string(r.divergence->nCells + 1));
      refused(d, read_divergence, ": the grid to reuse has", "RD9, G5");
    }
    {
      Table t = mean;
      t.column("nValidRec").data[1] = 3.;
      refused(t, read_mean, ": data row 2 has nValidRec 3, outside 0..2", "RD10, above NREC");
      t.column("nValidRec").data[1] = -1.;
      refused(t, read_mean, ": data row 2 has nValidRec -1, outside 0..2", "RD10, negative");
    }
  }

  group("readField: centres, values and weights in row order, FITS and ASCII, by name and by ASCII index; no "
        "weight column gives no weights");
  for (const char* format : kFormats) {
    const std::string file = temporary(std::string("field") + format);
    io::writeDivergence(file, r);
    const DivergenceGrid& g = *r.divergence;
    const io::Field f = io::readField(file, {"cellX", "cellY", "cellZ"}, "divergence", "weight");
    check(fixtures::same_bits(f.centres, g.centre) && fixtures::same_bits(f.values, g.divergence) &&
          fixtures::same_bits(f.weights, g.weight), std::string("by name, bit for bit ") + format);
    const io::Field n = io::readField(file, {"cellX", "cellY", "cellZ"}, "divergenceSmoothed");
    check(fixtures::same_bits(n.values, g.divergenceSmoothed) && n.weights.empty(),
          std::string("no weights ") + format);
    if (std::string(format) == ".txt") {
      const io::Field i = io::readField(file, {"1", "2", "3"}, "4", "7");
      check(fixtures::same_bits(i.centres, f.centres) && fixtures::same_bits(i.values, f.values) &&
            fixtures::same_bits(i.weights, f.weights), "ASCII by 0-based index");
    }
    otswap_nested([&] { io::readField(file, {"cellX", "cellY", "nope"}, "divergence"); },
                  std::string("RF1, a missing column ") + format);
    otswap_nested([&] { io::readField(file, {"cellX", "cellX", "cellZ"}, "divergence"); },
                  std::string("RF1, a column named twice ") + format);
    std::remove(file.c_str());
    otswap_nested([&] { io::readField(file, {"cellX", "cellY", "cellZ"}, "divergence"); },
                  std::string("RF1, no file ") + format);
  }

  group("readField with a delimiter and a comment character: a comma-separated file with % comments, its columns in "
        "another order, reads as its space-separated twin; without the delimiter its names are not found, without "
        "the comment character a comment line that is not data is refused");
  {
    const fixtures::LatticeField lf;
    const std::size_t commentRow = 7;
    const std::string spaced = temporary("field_spaced.txt");
    const std::string comma = temporary("field_comma.txt");
    {
      std::FILE* s = std::fopen(spaced.c_str(), "w");
      std::FILE* c = std::fopen(comma.c_str(), "w");
      std::fprintf(s, "# 12 x 17 x 9 cells\n### x y z value weight\n");
      std::fprintf(c, "### weight,z,value,x,y\n");
      for (std::size_t i = 0; i < lf.values.size(); ++i) {
        std::fprintf(s, "%.17g %.17g %.17g %.17g %.17g\n", lf.centres[3*i], lf.centres[3*i+1], lf.centres[3*i+2],
                     lf.values[i], lf.weights[i]);
        if (i == commentRow) std::fprintf(c, "%% 1,2,3,4,5\n");
        std::fprintf(c, "%.17g,%.17g,%.17g,%.17g,%.17g\n", lf.weights[i], lf.centres[3*i+2], lf.values[i],
                     lf.centres[3*i], lf.centres[3*i+1]);
      }
      std::fclose(s);
      std::fclose(c);
    }
    const io::Field a = io::readField(spaced, {"x", "y", "z"}, "value", "weight");
    const io::Field b = io::readField(comma, {"x", "y", "z"}, "value", "weight", ',', '%');
    check(fixtures::same_bits(a.centres, lf.centres) && fixtures::same_bits(a.values, lf.values) &&
          fixtures::same_bits(a.weights, lf.weights), "the space-separated file, with the defaults");
    check(fixtures::same_bits(b.centres, a.centres) && fixtures::same_bits(b.values, a.values) &&
          fixtures::same_bits(b.weights, a.weights), "the comma-separated twin, bit for bit");
    otswap_nested([&] { io::readField(comma, {"x", "y", "z"}, "value", "weight"); },
                  "the default delimiter: the names are not found");
    otswap_nested([&] { io::readField(comma, {"x", "y", "z"}, "value", "weight", ','); },
                  "the default comment character: the comment line is not data");
    std::remove(spaced.c_str());
    std::remove(comma.c_str());
  }

  group("readField refuses a value that is not finite, naming the file, the data row from 1 and the column as "
        "selected");
  {
    const std::string file = temporary("field_not_finite.txt");
    const struct { const char* row; const char* column; const char* word; const char* selector; } cases[] = {
      {"0.5 1.5 0.5 nan 1", "value", "nan", "value"},
      {"0.5 1.5 0.5 3 inf", "weight", "inf", "weight"},
      {"0.5 -inf 0.5 3 1", "y", "-inf", "y"},
      {"0.5 1.5 0.5 -inf 1", "3", "-inf", "value"},
    };
    for (const auto& c : cases) {
      std::FILE* out = std::fopen(file.c_str(), "w");
      std::fprintf(out, "### x y z value weight\n0.5 0.5 0.5 1 1\n1.5 0.5 0.5 2 1\n%s\n1.5 1.5 0.5 4 1\n", c.row);
      std::fclose(out);
      const bool byIndex = std::string(c.column) == "3";
      throws_naming([&] {
        if (byIndex) io::readField(file, {"0", "1", "2"}, "3", "4");
        else io::readField(file, {"x", "y", "z"}, "value", "weight");
      }, file + ": data row 3 has " + c.word + " in column " + c.column + "; every value must be finite",
      std::string(c.word) + " in " + c.selector + (byIndex ? ", by index" : ""));
    }
    std::remove(file.c_str());
  }

  group("the divergence file of a field result: readField reads it back, readDivergence refuses it");
  {
    const fixtures::LatticeField lf;
    Config c;
    c.log = nullptr;
    c.keepDivergence = true;
    c.smoothing.unit = LengthUnit::Coordinate;
    c.smoothing.sigma = 1.5;
    const Result fr = findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights);
    for (const char* format : kFormats) {
      const std::string file = temporary(std::string("field_result") + format);
      io::writeDivergence(file, fr);
      const io::Field f = io::readField(file, {"cellX", "cellY", "cellZ"}, "divergence", "weight");
      const Result back = findVoidsBoxFromField(f.centres, f.values, c, f.weights);
      check(fixtures::same_bits(f.centres, fr.divergence->centre) &&
            fixtures::same_bits(f.values, fr.divergence->divergence) &&
            back.divergence->voidIndex == fr.divergence->voidIndex &&
            back.divergence->overdenseIndex == fr.divergence->overdenseIndex,
            std::string("read back, the same partition ") + format);
      throws_naming([&] { io::readDivergence(file); }, "keyword NREC is missing",
                    std::string("readDivergence refuses it ") + format);
      std::remove(file.c_str());
    }
  }

  return report("test_readers");
}
