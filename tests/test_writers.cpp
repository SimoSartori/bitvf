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
 *  @file tests/test_writers.cpp
 *
 *  @brief The writers of bitvf::io: columns, keywords, precision, errors
 *  and file names.
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iterator>
#include <limits>
#include <locale>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "otswap/Maths.h"
#include "otswap/io.h"

#include "bitvf/io.h"
#include "internal.h"
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

  void threads (const int n)
  {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
  }

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "bitvf_test_writers_" + name;
  }

  bool exists (const std::string& file)
  {
    std::ifstream in(file);
    return in.good();
  }

  std::string bytes_of (const std::string& file)
  {
    std::ifstream in(file, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  }

  /// A search on tracers uniform in [0, 100)^3, kOut more beyond the
  /// limits [-2, 102]^3, displaced at random and towards the centre, so
  /// that both kinds of basin appear on cells of cbrt(10) mps; every
  /// output kept, IDs of both signs and three realizations averaged.
  struct Input {
    std::vector<double> tracers, displacement;
    std::vector<std::uint8_t> valid;
    Config config;
    std::optional<otswap::BoxLimits> limits = box(-2., 102.);
    std::vector<std::int64_t> ids;
    const DivergenceGrid* grid = nullptr;

    Input ()
    {
      std::mt19937 rng(2026);
      auto uniform = [&rng] (const double lo, const double hi) {
        return lo + (hi - lo) * (rng() / 4294967296.);
      };
      const std::size_t n = kN + kOut;
      for (std::size_t i = 0; i < 3 * kN; ++i) tracers.push_back(uniform(0., 100.));
      for (std::size_t i = 0; i < 3 * kOut; ++i) tracers.push_back(150. + i);
      for (std::size_t r = 0; r < kR; ++r)
        for (std::size_t i = 0; i < 3 * n; ++i)
          displacement.push_back(uniform(-2., 2.) - (tracers[i] - 50.) / 50.);
      for (std::size_t i = 0; i < n; ++i)
        ids.push_back(i % 2 == 0 ? (std::int64_t)(1000 + i) : -(std::int64_t)(1000 + i));
      config.log = nullptr;
      config.keepDivergence = true;
      config.keepUnsmoothedDivergence = true;
      config.keepAssignment = true;
      config.nRealizationsMean = 3;
      config.divergenceCellSize = otswap::maths::cbrt(10.);
    }

    Result run () const
    {
      return findVoidsBox(tracers, displacement, valid, kMps, config, limits, ids, grid);
    }
  };

  /// A column as the header of io.h lists it, with the value of each row.
  struct Expected {
    std::string name;
    char type;
    std::function<double(std::size_t)> real;
    std::function<std::int64_t(std::size_t)> integer;
  };

  Expected real (const std::string& name, std::function<double(std::size_t)> f)
  {
    return {name, 'D', std::move(f), {}};
  }

  Expected integer (const std::string& name, const char type, std::function<std::int64_t(std::size_t)> f)
  {
    return {name, type, {}, std::move(f)};
  }

  /// The columns of a catalogue; @p field for a field result.
  std::vector<Expected> catalogue_columns (const BasinCatalog& c, const bool field = false)
  {
    std::vector<Expected> e = {
      integer("index", 'K', [] (std::size_t b) { return (std::int64_t)b; }),
      real("centreX", [&c] (std::size_t b) { return c.centre[3*b]; }),
      real("centreY", [&c] (std::size_t b) { return c.centre[3*b+1]; }),
      real("centreZ", [&c] (std::size_t b) { return c.centre[3*b+2]; }),
      real("radius", [&c] (std::size_t b) { return c.radius[b]; }),
      real("radiusExt", [&c] (std::size_t b) { return c.radiusExt[b]; })};
    if (!field) {
      e.push_back(real("deltaAv", [&c] (std::size_t b) { return c.deltaAv[b]; }));
      e.push_back(real("delta2Mps", [&c] (std::size_t b) { return c.delta2Mps[b]; }));
    }
    e.push_back(integer("edgeFlag", 'J', [&c] (std::size_t b) { return (std::int64_t)c.edgeFlag[b]; }));
    e.push_back(real("contamination", [&c] (std::size_t b) { return c.contamination[b]; }));
    if (field) e.push_back(real("fieldCentre", [&c] (std::size_t b) { return c.fieldCentre[b]; }));
    else e.push_back(integer("nTracers", 'K', [&c] (std::size_t b) { return (std::int64_t)c.nTracers[b]; }));
    return e;
  }

  /// The columns of a divergence grid; @p field for a field result.
  std::vector<Expected> divergence_columns (const DivergenceGrid& g, const bool field = false)
  {
    std::vector<Expected> e = {
      integer("index", 'K', [] (std::size_t c) { return (std::int64_t)c; }),
      real("cellX", [&g] (std::size_t c) { return g.centre[3*c]; }),
      real("cellY", [&g] (std::size_t c) { return g.centre[3*c+1]; }),
      real("cellZ", [&g] (std::size_t c) { return g.centre[3*c+2]; }),
      real("divergence", [&g] (std::size_t c) { return g.divergence[c]; })};
    if (!g.divergenceUnsmoothed.empty())
      e.push_back(real("divergenceUnsmoothed", [&g] (std::size_t c) { return g.divergenceUnsmoothed[c]; }));
    e.push_back(real("volume", [&g] (std::size_t c) { return g.volume[c]; }));
    e.push_back(real("weight", [&g] (std::size_t c) { return g.weight[c]; }));
    if (!field) {
      e.push_back(integer("nDispl", 'K', [&g] (std::size_t c) { return (std::int64_t)g.nDisplacements[c]; }));
      e.push_back(integer("nTracers", 'K', [&g] (std::size_t c) { return (std::int64_t)g.nTracers[c]; }));
    }
    e.push_back(integer("voidIndex", 'K', [&g] (std::size_t c) { return g.voidIndex[c]; }));
    e.push_back(integer("overdenseIndex", 'K', [&g] (std::size_t c) { return g.overdenseIndex[c]; }));
    return e;
  }

  std::vector<Expected> galaxies_columns (const Result& r)
  {
    const TracerAssignment& a = *r.assignment;
    std::vector<Expected> e = {
      integer("id", 'K', [&r] (std::size_t i) { return r.ids.empty() ? (std::int64_t)i : r.ids[i]; }),
      integer("cellIndex", 'K', [&a] (std::size_t i) { return a.cellIndex[i]; }),
      integer("voidIndex", 'K', [&a] (std::size_t i) { return a.voidIndex[i]; }),
      integer("overdenseIndex", 'K', [&a] (std::size_t i) { return a.overdenseIndex[i]; }),
      real("distance", [&a] (std::size_t i) { return a.distance[i]; }),
      integer("inMask", 'J', [&a] (std::size_t i) { return (std::int64_t)a.inMask[i]; })};
    if (r.limits) e.push_back(integer("inBox", 'J', [&a] (std::size_t i) { return (std::int64_t)a.inBox[i]; }));
    e.push_back(integer("used", 'J', [&a] (std::size_t i) { return (std::int64_t)a.used[i]; }));
    if (!r.uncorrected.empty())
      e.push_back(integer("uncorrected", 'J', [&r] (std::size_t i) { return (std::int64_t)r.uncorrected[i]; }));
    if (!r.outsideGrid.empty())
      e.push_back(integer("outsideGrid", 'J', [&r] (std::size_t i) { return (std::int64_t)r.outsideGrid[i]; }));
    return e;
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

  bool same_bits (const double a, const double b)
  {
    if (std::isnan(a) && std::isnan(b)) return true;
    return std::memcmp(&a, &b, sizeof(double)) == 0;
  }

  /// The column names of a file, in order: TTYPEn in FITS, the ### line in
  /// ASCII.
  std::vector<std::string> names_of (const std::string& file)
  {
    std::vector<std::string> names;
    if (file.size() > 5 && file.compare(file.size() - 5, 5, ".fits") == 0) {
      const std::vector<otswap::io::Keyword> k = otswap::io::readKeywords(file);
      for (std::size_t c = 1; ; ++c) {
        bool found = false;
        for (const otswap::io::Keyword& w : k)
          if (w.name == "TTYPE" + std::to_string(c)) {
            names.push_back(w.value);
            found = true;
          }
        if (!found) break;
      }
      return names;
    }
    std::ifstream in(file);
    std::string line;
    while (std::getline(in, line))
      if (line.compare(0, 3, "###") == 0) {
        std::istringstream words(line.substr(3));
        std::string w;
        while (words >> w) names.push_back(w);
        break;
      }
    return names;
  }

  /// The column descriptions of a file, in order: TCOMMn in FITS, the
  /// "# name description" lines in ASCII, which hold them whole.
  std::vector<std::string> descriptions_of (const std::string& file)
  {
    std::vector<std::string> descriptions;
    if (file.size() > 5 && file.compare(file.size() - 5, 5, ".fits") == 0) {
      const std::vector<otswap::io::Keyword> k = otswap::io::readKeywords(file);
      const std::size_t n = names_of(file).size();
      for (std::size_t c = 1; c <= n; ++c) {
        std::string text;
        for (const otswap::io::Keyword& w : k)
          if (w.name == "TCOMM" + std::to_string(c)) text = w.value;
        descriptions.push_back(text);
      }
      return descriptions;
    }
    std::ifstream in(file);
    std::string line;
    while (std::getline(in, line) && line.compare(0, 3, "###") != 0) {
      if (line.compare(0, 2, "# ") != 0) continue;
      const std::size_t end = line.find(' ', 2);
      const std::size_t start = end == std::string::npos ? end : line.find_first_not_of(' ', end);
      descriptions.push_back(start == std::string::npos ? "" : line.substr(start));
    }
    return descriptions;
  }

  /// The file holds exactly @p columns, in order, with @p nRows rows equal
  /// to theirs: bit for bit at @p digits significant digits in ASCII (17
  /// is exact), and always in FITS.
  bool reads_back (const std::string& file, const std::vector<Expected>& columns, const std::size_t nRows,
                   const int digits, const std::string& what)
  {
    std::vector<std::string> expectedNames, reals, integers;
    for (const Expected& e : columns) {
      expectedNames.push_back(e.name);
      (e.type == 'D' ? reals : integers).push_back(e.name);
    }
    const std::vector<std::string> names = names_of(file);
    if (names != expectedNames) {
      check(false, what + ": the column names and their order");
      return false;
    }
    const otswap::io::Table t = otswap::io::read(file, reals, integers);
    if (t.nRows != nRows) {
      check(false, what + ": " + std::to_string(t.nRows) + " rows, " + std::to_string(nRows) + " expected");
      return false;
    }
    const bool fits = file.size() > 5 && file.compare(file.size() - 5, 5, ".fits") == 0;
    std::size_t d = 0, k = 0;
    for (const Expected& e : columns) {
      bool equal = true;
      for (std::size_t row = 0; equal && row < nRows; ++row) {
        if (e.type == 'D') {
          const double expected = fits ? e.real(row) : as_written(e.real(row), digits);
          equal = same_bits(t.values[row * t.nColumns + d], expected);
        }
        else {
          equal = t.integers[row * t.nIntegerColumns + k] == e.integer(row);
        }
      }
      check(equal, what + ": column " + e.name);
      if (e.type == 'D') ++d; else ++k;
    }
    return true;
  }

  /// A keyword as written: a double compared bit for bit after strtod, or
  /// a text compared as is.
  struct Key {
    std::string name;
    bool number;
    double value;
    std::string text;
  };

  Key number (const std::string& name, const double value) { return {name, true, value, ""}; }
  Key text (const std::string& name, const std::string& value) { return {name, false, 0., value}; }

  std::vector<Key> expected_keywords (const std::string& product, const Result& r, const std::string& count,
                                      const std::size_t n)
  {
    const bool field = r.nTracers == 0;
    std::vector<Key> k = {
      text("PRODUCT", product), text("BITVFV", internal::version), text("GEOMETRY", "box"),
      text(count, std::to_string(n))};
    if (!field) {
      k.push_back(text("NREC", std::to_string(r.nRealizations)));
      k.push_back(text("NRECMEAN", std::to_string(r.config.nRealizationsMean)));
    }
    k.push_back(number("SIGMA", r.config.smoothing.sigma));
    k.push_back(text("SIGUNIT", r.config.smoothing.unit == LengthUnit::Mps ? "mps" : "coordinate"));
    if (!field || !std::isnan(r.mps)) {
      k.push_back(number("CELLSIZE", *r.config.divergenceCellSize));
      k.push_back(number("MPS", r.mps));
    }
    if (r.randomSeed) k.push_back(text("RANDSEED", std::to_string(*r.randomSeed)));
    if (r.limits) {
      k.push_back(text("BOXLIM", "1"));
      const char* axis[3] = {"X", "Y", "Z"};
      for (int c = 0; c < 3; ++c) {
        k.push_back(number(std::string("BOX") + axis[c] + "LO", r.limits->lower[c]));
        k.push_back(number(std::string("BOX") + axis[c] + "HI", r.limits->upper[c]));
      }
    }
    if (product == "divergence") {
      const DivergenceGrid& g = *r.divergence;
      const char* axis[3] = {"X", "Y", "Z"};
      for (int c = 0; c < 3; ++c) k.push_back(number(std::string("LOWER") + axis[c], g.lower[c]));
      for (int c = 0; c < 3; ++c) k.push_back(number(std::string("UPPER") + axis[c], g.upper[c]));
      k.push_back(number("CELLSIDE", g.cellSize));
      for (int c = 0; c < 3; ++c) k.push_back(text(std::string("NCELL") + axis[c], std::to_string(g.nCellsAxis[c])));
    }
    return k;
  }

  /// The keywords of the file from PRODUCT on are exactly @p expected, in
  /// order.
  void keywords_are (const std::string& file, const std::vector<Key>& expected, const std::string& what)
  {
    const std::vector<otswap::io::Keyword> all = otswap::io::readKeywords(file);
    std::size_t first = 0;
    while (first < all.size() && all[first].name != "PRODUCT") ++first;
    if (all.size() - first != expected.size()) {
      check(false, what + ": " + std::to_string(all.size() - first) + " keywords from PRODUCT on, " +
                   std::to_string(expected.size()) + " expected");
      return;
    }
    for (std::size_t i = 0; i < expected.size(); ++i) {
      const otswap::io::Keyword& w = all[first + i];
      const Key& e = expected[i];
      const bool equal = w.name == e.name &&
                         (e.number ? same_bits(std::strtod(w.value.c_str(), nullptr), e.value) : w.value == e.text);
      check(equal, what + ": keyword " + e.name + " (read " + w.name + " = " + w.value + ")");
    }
  }

  /// Each of the four products, written to <stem><format>, read back and
  /// removed; a field result without the galaxies.
  void round_trip (const Result& r, const std::string& stem, const std::string& what)
  {
    const bool field = r.nTracers == 0;
    for (const char* format : kFormats) {
      const bool fits = std::string(format) == ".fits";
      const std::string tag = what + " " + format;
      const int digits = fits ? 17 : 9;

      std::string file = temporary(stem + "_voids" + format);
      io::writeVoids(file, r);
      reads_back(file, catalogue_columns(r.voids, field), r.voids.nBasins, digits, tag + " voids");
      keywords_are(file, expected_keywords("voids", r, "NBASINS", r.voids.nBasins), tag + " voids");
      std::remove(file.c_str());

      file = temporary(stem + "_overdense" + format);
      io::writeOverdense(file, r);
      reads_back(file, catalogue_columns(r.overdense, field), r.overdense.nBasins, digits, tag + " overdense");
      keywords_are(file, expected_keywords("overdense_basins", r, "NBASINS", r.overdense.nBasins),
                   tag + " overdense");
      std::remove(file.c_str());

      file = temporary(stem + "_divergence" + format);
      io::writeDivergence(file, r);
      reads_back(file, divergence_columns(*r.divergence, field), r.divergence->nCells, 17, tag + " divergence");
      keywords_are(file, expected_keywords("divergence", r, "NCELLS", r.divergence->nCells), tag + " divergence");
      std::remove(file.c_str());

      if (field) continue;

      file = temporary(stem + "_galaxies" + format);
      io::writeGalaxies(file, r);
      reads_back(file, galaxies_columns(r), r.nTracers, digits, tag + " galaxies");
      keywords_are(file, expected_keywords("galaxies", r, "NOBJECTS", r.nTracers), tag + " galaxies");
      std::remove(file.c_str());
    }
  }

  /// @p write must raise a bitvf::Error naming @p message, and leave no
  /// file.
  template <typename W>
  void refused (W&& write, const Result& r, const std::string& message, const std::string& what)
  {
    for (const char* format : kFormats) {
      const std::string file = temporary(std::string("refused") + format);
      std::remove(file.c_str());
      throws_naming([&] { write(file, r); }, message, what + " " + format);
      check(!exists(file), what + " " + format + ": no file left");
      std::remove(file.c_str());
    }
  }

  void write_voids (const std::string& f, const Result& r) { io::writeVoids(f, r); }
  void write_overdense (const std::string& f, const Result& r) { io::writeOverdense(f, r); }
  void write_divergence (const std::string& f, const Result& r) { io::writeDivergence(f, r); }
  void write_galaxies (const std::string& f, const Result& r) { io::writeGalaxies(f, r); }

  /// The bytes of the four files of @p r in @p format.
  std::string written_bytes (const Result& r, const std::string& format)
  {
    std::string all;
    const std::string file = temporary("threads" + format);
    for (auto write : {write_voids, write_overdense, write_divergence, write_galaxies}) {
      write(file, r);
      all += bytes_of(file);
      std::remove(file.c_str());
    }
    return all;
  }

}

int main ()
{
  const Input in;
  const Result r = in.run();
  check(r.voids.nBasins > 0 && r.overdense.nBasins > 0 && r.counts.tracersOutsideBox == kOut &&
        r.limits && r.divergence && !r.divergence->divergenceUnsmoothed.empty() && r.assignment,
        "the fixture: basins of both kinds, tracers outside the limits, every output kept");

  group("FITS and ASCII round trips of every product: names, order, values; keywords in order");
  round_trip(r, "full", "full");

  group("keywords: RANDSEED exactly when set, BOXLIM and BOX* exactly with limits, NRECMEAN always");
  {
    Result seeded = r;
    seeded.randomSeed = 4294967295u;
    round_trip(seeded, "seeded", "RANDSEED");

    Input free = in;
    free.limits.reset();
    free.config.nRealizationsMean = 1;
    free.config.smoothing.unit = LengthUnit::Coordinate;
    free.config.smoothing.sigma = 0.1 + 0.2;
    const Result f = free.run();
    check(!f.limits && !f.divergence->limits, "no limits recorded without limits");
    round_trip(f, "free", "no limits");
  }

  group("optional columns: divergenceUnsmoothed, inBox, uncorrected, outsideGrid exactly when they apply");
  {
    Input plain = in;
    plain.config.keepUnsmoothedDivergence = false;
    const Result p = plain.run();
    check(p.divergence->divergenceUnsmoothed.empty(), "no unsmoothed divergence kept");
    round_trip(p, "plain", "no divergenceUnsmoothed");

    Result flagged = r;
    flagged.uncorrected.assign(r.nTracers, 0);
    flagged.uncorrected[3] = 1;
    flagged.outsideGrid.assign(r.nTracers, 0);
    flagged.outsideGrid[5] = 1;
    round_trip(flagged, "flagged", "uncorrected and outsideGrid");

    Input reuse = in;
    reuse.grid = &*r.divergence;
    reuse.config.keepDivergence = false;
    reuse.config.keepUnsmoothedDivergence = false;
    reuse.config.divergenceCellSize.reset();
    reuse.limits.reset();
    const Result g = reuse.run();
    check(!g.outsideGrid.empty() && !g.limits, "a reused grid: outsideGrid set, no limits");
    for (const char* format : kFormats) {
      const std::string file = temporary(std::string("reuse") + format);
      io::writeGalaxies(file, g);
      reads_back(file, galaxies_columns(g), g.nTracers, std::string(format) == ".fits" ? 17 : 9,
                 std::string("a reused grid ") + format);
      std::remove(file.c_str());
    }
  }

  group("descriptions: every column of every product described, FITS equal to ASCII, none truncated");
  {
    Result flagged = r;
    flagged.uncorrected.assign(r.nTracers, 0);
    flagged.outsideGrid.assign(r.nTracers, 0);
    const char* const product[4] = {"voids", "overdense", "divergence", "galaxies"};
    int p = 0;
    for (auto write : {write_voids, write_overdense, write_divergence, write_galaxies}) {
      const std::string fits = temporary("descriptions.fits"), ascii = temporary("descriptions.txt");
      write(fits, flagged);
      write(ascii, flagged);
      const std::vector<std::string> f = descriptions_of(fits), a = descriptions_of(ascii);
      const std::vector<std::string> names = names_of(fits);
      check(a.size() == names.size() && f == a, std::string("FITS descriptions whole: ") + product[p]);
      for (std::size_t c = 0; c < a.size(); ++c)
        check(!a[c].empty(), std::string(product[p]) + " column " + names[c] + " described");
      std::remove(fits.c_str());
      std::remove(ascii.c_str());
      ++p;
    }
  }

  group("id: the IDs, negative ones included, or the row from 0 without");
  {
    Input unnamed = in;
    unnamed.ids.clear();
    const Result u = unnamed.run();
    for (const char* format : kFormats) {
      const std::string file = temporary(std::string("ids") + format);
      io::writeGalaxies(file, r);
      otswap::io::Table t = otswap::io::read(file, {}, {"id"});
      check(t.integers == r.ids && r.ids[1] < 0, std::string("the IDs ") + format);
      io::writeGalaxies(file, u);
      t = otswap::io::read(file, {}, {"id"});
      bool rows = t.nRows == u.nTracers;
      for (std::size_t i = 0; rows && i < t.nRows; ++i) rows = t.integers[i] == (std::int64_t)i;
      check(rows, std::string("the rows without IDs ") + format);
      std::remove(file.c_str());
    }
  }

  group("K columns: values at and above 2^31, and the largest 64-bit integer, written exactly");
  {
    Result big = r;
    const std::int64_t top = std::numeric_limits<std::int64_t>::max();
    big.voids.nTracers[0] = (std::size_t)1 << 31;
    big.voids.nTracers.back() = (std::size_t)top;
    big.overdense.nTracers[0] = ((std::size_t)1 << 40) + 7;
    big.ids[0] = top;
    big.ids[1] = std::numeric_limits<std::int64_t>::min();
    big.assignment->cellIndex[0] = (std::int64_t)1 << 33;
    big.assignment->voidIndex[0] = ((std::int64_t)1 << 31) + 1;
    big.assignment->overdenseIndex[1] = (std::int64_t)1 << 50;
    big.divergence->nDisplacements[0] = (std::uint64_t)top;
    big.divergence->nTracers[1] = (std::uint64_t)1 << 32;
    big.divergence->voidIndex[2] = (std::int64_t)1 << 45;
    big.divergence->overdenseIndex[3] = ((std::int64_t)1 << 31) * 3;
    round_trip(big, "big", "K values");
  }

  group("errors WR1-WR10, each before the file is opened");
  {
    Result lightcone = r;
    lightcone.geometry = otswap::Geometry::Lightcone;
    for (auto write : {write_voids, write_overdense, write_divergence, write_galaxies})
      refused(write, lightcone, "writing a lightcone result is not implemented", "WR1");

    auto malformed = [&r] (auto&& change, auto write, const std::string& message, const std::string& what) {
      Result bad = r;
      change(bad);
      refused(write, bad, message, what);
    };

    malformed([] (Result& b) { b.voids.centre.pop_back(); }, write_voids,
              "cannot write voids.centre: it holds", "WR2, voids centre");
    malformed([] (Result& b) { b.voids.radius.push_back(1.); }, write_voids,
              "cannot write voids.radius: it holds", "WR2, voids radius");
    malformed([] (Result& b) { b.overdense.centreSky = {1., 2.}; }, write_overdense,
              "cannot write overdense.centreSky: it holds 2 entries where 0 are expected", "WR2, centreSky");
    malformed([] (Result& b) { b.overdense.nTracers.pop_back(); }, write_overdense,
              "cannot write overdense.nTracers", "WR2, overdense nTracers");
    malformed([] (Result& b) { b.voids.nTracers[1] = (std::size_t)std::numeric_limits<std::int64_t>::max() + 1; },
              write_voids, "cannot write voids.nTracers: entry 1 is 9223372036854775808, above the largest",
              "WR2, a count above 2^63 - 1");
    malformed([] (Result& b) { b.divergence->nDisplacements[2] = std::numeric_limits<std::uint64_t>::max(); },
              write_divergence, "cannot write divergence.nDisplacements: entry 2 is", "WR2, nDispl too large");
    malformed([] (Result& b) { b.ids.pop_back(); }, write_galaxies, "cannot write ids: it holds", "WR2, ids");
    malformed([] (Result& b) { b.uncorrected = {0}; }, write_galaxies, "cannot write uncorrected", "WR2, uncorrected");
    malformed([] (Result& b) { b.outsideGrid = {0, 1}; }, write_galaxies, "cannot write outsideGrid",
              "WR2, outsideGrid");
    malformed([] (Result& b) { b.assignment->distance.pop_back(); }, write_galaxies,
              "cannot write assignment.distance", "WR2, assignment");
    malformed([] (Result& b) { b.assignment->inBox.clear(); }, write_galaxies,
              "cannot write assignment.inBox", "WR2, inBox");
    malformed([] (Result& b) { b.divergence->centre.push_back(0.); }, write_divergence,
              "cannot write divergence.centre", "WR2, grid centre");
    malformed([] (Result& b) { b.divergence->weight.pop_back(); }, write_divergence,
              "cannot write divergence.weight", "WR2, grid weight");
    malformed([] (Result& b) { b.divergence->divergenceUnsmoothed.pop_back(); }, write_divergence,
              "cannot write divergence.divergenceUnsmoothed", "WR2, divergenceUnsmoothed");

    malformed([] (Result& b) { b.divergence.reset(); }, write_divergence,
              "the result holds no divergence grid: set Config::keepDivergence", "WR3");
    malformed([] (Result& b) { b.assignment.reset(); }, write_galaxies,
              "the result holds no tracer assignment: set Config::keepAssignment", "WR4");

    malformed([] (Result& b) { b.divergence->cellSizeMps = std::nextafter(b.divergence->cellSizeMps, 10.); },
              write_divergence, "cannot write the divergence: the grid to reuse has a cell side of",
              "WR5, G9");
    malformed([] (Result& b) { b.divergence->nCells += 1; }, write_divergence,
              "cannot write the divergence: the grid to reuse has", "WR5, G5");
    malformed([] (Result& b) { b.divergence->raMin = 0.; }, write_divergence,
              "cannot write the divergence: the grid to reuse is a box grid with a field of a lightcone grid",
              "WR5, G7");

    malformed([] (Result& b) { b.mps = std::nextafter(b.mps, 20.); }, write_divergence,
              "the result is malformed: its divergence grid has an mps of", "WR6, mps");
    malformed([] (Result& b) { b.config.divergenceCellSize = 2. * *b.config.divergenceCellSize; },
              write_divergence, "the result is malformed: its divergence grid has a cell size of", "WR6, cell size");
    malformed([] (Result& b) { b.limits.reset(); }, write_divergence,
              "the result is malformed: the box limits of its divergence grid differ", "WR6, limits absent");
    malformed([] (Result& b) { b.limits->upper[1] = std::nextafter(b.limits->upper[1], 200.); }, write_divergence,
              "the result is malformed: the box limits of its divergence grid differ", "WR6, one bound");
    malformed([] (Result& b) { b.divergence->nRealizations += 1; }, write_divergence,
              "the result is malformed: its divergence grid records 3 realizations, but the search's are 2",
              "WR6, nRealizations");
    malformed([] (Result& b) { b.config.nRealizationsMean = 4; }, write_divergence,
              "the result is malformed: its divergence grid records 3 realizations averaged, but its "
              "configuration records 4", "WR6, nRealizationsMean");

    for (auto write : {write_voids, write_overdense, write_divergence, write_galaxies}) {
      malformed([] (Result& b) { b.config.divergenceCellSize.reset(); }, write,
                "the result is malformed: it records no divergence cell size", "WR7");
      malformed([] (Result& b) { b.config.nRealizationsMean = 0; }, write,
                "the result is malformed: it records 0 realizations averaged", "WR10");
    }

    throws_naming([] { io::fileName(io::Product::Voids, "dir", "", false); }, "the file name is empty", "WR8");

    for (const char* format : kFormats) {
      const std::string file = std::string("/nonexistent_bitvf_directory/x") + format;
      bool nested = false;
      try {
        io::writeVoids(file, r);
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
      check(nested, std::string("WR9: an unwritable path, otswap's error nested ") + format);
    }
  }

  group("fileName: every product, zspace, directory with and without '/', empty directory");
  {
    check(io::fileName(io::Product::Voids, "out", "run.fits", false) == "out/voids_run.fits", "voids");
    check(io::fileName(io::Product::Overdense, "out/", "run.fits", false) == "out/overdense_basins_run.fits",
          "overdense, directory ending in '/'");
    check(io::fileName(io::Product::Divergence, "", "run.txt", false) == "divergence_run.txt",
          "divergence, empty directory");
    check(io::fileName(io::Product::Galaxies, "/a/b", "run", false) == "/a/b/galaxies_run", "galaxies");
    check(io::fileName(io::Product::Voids, "out", "run.fits", true) == "out/voids_zspace_run.fits", "voids zspace");
    check(io::fileName(io::Product::Overdense, "", "r", true) == "overdense_basins_zspace_r", "overdense zspace");
    check(io::fileName(io::Product::Divergence, "d/", "r", true) == "d/divergence_zspace_r", "divergence zspace");
    check(io::fileName(io::Product::Galaxies, "/", "r", true) == "/galaxies_zspace_r", "galaxies zspace, root");
  }

  group("a field result, without and with mps: catalogues with fieldCentre and without deltaAv, delta2Mps, "
        "nTracers; divergence without nDispl, nTracers; no NREC, NRECMEAN; CELLSIZE and MPS only with mps");
  {
    const fixtures::LatticeField lf;
    Config c;
    c.log = nullptr;
    c.keepDivergence = true;
    c.keepUnsmoothedDivergence = true;
    c.smoothing.unit = LengthUnit::Coordinate;
    c.smoothing.sigma = 1.5;
    const Result f = findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights);
    check(f.nTracers == 0 && f.voids.nBasins > 0 && f.overdense.nBasins > 0 && std::isnan(f.mps),
          "the field result: both kinds of basin, no mps");
    round_trip(f, "field", "field, no mps");
    const Result m = findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights, 3.);
    round_trip(m, "field_mps", "field, mps 3");
    c.keepUnsmoothedDivergence = false;
    round_trip(findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights), "field_plain",
               "field, no divergenceUnsmoothed");

    const char* const product[3] = {"voids", "overdense", "divergence"};
    int p = 0;
    for (auto write : {write_voids, write_overdense, write_divergence}) {
      const std::string fits = temporary("field_descriptions.fits"), ascii = temporary("field_descriptions.txt");
      write(fits, f);
      write(ascii, f);
      const std::vector<std::string> d = descriptions_of(fits), a = descriptions_of(ascii);
      const std::vector<std::string> names = names_of(fits);
      bool described = a.size() == names.size() && d == a;
      for (const std::string& text : a) described = described && !text.empty();
      check(described, std::string("field ") + product[p] + ": every column described, FITS equal to ASCII");
      std::remove(fits.c_str());
      std::remove(ascii.c_str());
      ++p;
    }
  }

  group("errors on a field result: writeGalaxies refused (WR11); its grid checks (WR12); its columns and "
        "records (WR13); fieldCentre in a finder's result refused");
  {
    const fixtures::LatticeField lf;
    Config c;
    c.log = nullptr;
    c.keepDivergence = true;
    c.smoothing.unit = LengthUnit::Coordinate;
    c.smoothing.sigma = 1.5;
    const Result f = findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights);
    const Result m = findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights, 3.);
    refused(write_galaxies, f, "the result is a search on a field (findVoidsBoxFromField), which has no tracers",
            "WR11");

    auto changed = [] (const Result& base, auto&& change) {
      Result bad = base;
      change(bad);
      return bad;
    };
    refused(write_divergence, changed(f, [] (Result& b) { b.divergence->cellSize = -1.; }),
            "cannot write the divergence: the grid has a cell size of", "WR12 cell side");
    refused(write_divergence, changed(f, [] (Result& b) {
              b.divergence->upper[1] = std::numeric_limits<double>::infinity(); }),
            "cannot write the divergence: the grid has a bound that is not finite on axis y", "WR12 bound");
    refused(write_divergence, changed(f, [] (Result& b) { b.divergence->nCellsAxis[2] += 1; }),
            "cells for (", "WR12 counts");
    refused(write_voids, changed(f, [] (Result& b) { b.voids.deltaAv.assign(b.voids.nBasins, 0.); }),
            "cannot write voids.deltaAv", "WR13 deltaAv");
    refused(write_overdense, changed(f, [] (Result& b) { b.overdense.nTracers.assign(b.overdense.nBasins, 1); }),
            "cannot write overdense.nTracers", "WR13 nTracers");
    refused(write_voids, changed(f, [] (Result& b) { b.voids.fieldCentre.pop_back(); }),
            "cannot write voids.fieldCentre", "WR13 fieldCentre short");
    refused(write_voids, changed(m, [] (Result& b) { b.config.divergenceCellSize.reset(); }),
            "it records no divergence cell size", "WR13 with mps, no cell size");
    refused(write_voids, changed(r, [] (Result& b) { b.voids.fieldCentre.assign(b.voids.nBasins, 0.); }),
            "where 0 are expected", "fieldCentre in a finder's result");
  }

  group("the same bytes with 1 and 4 threads, in both formats");
  {
    threads(1);
    const Result one = in.run();
    threads(4);
    const Result four = in.run();
    for (const char* format : kFormats) {
      const std::string a = written_bytes(one, format);
      const std::string b = written_bytes(four, format);
      check(!a.empty() && a == b, std::string("identical bytes ") + format);
    }
  }

  return report("test_writers");
}
