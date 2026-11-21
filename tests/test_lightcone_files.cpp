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
 *  @file tests/test_lightcone_files.cpp
 *
 *  @brief The divergence file of a lightcone: its columns, record and
 *  profile part, reading it back, the checks of the writer and of the
 *  reader.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "otswap/Maths.h"
#include "otswap/OT.h"
#include "otswap/io.h"

#include "bitvf/io.h"
#include "bitvf/VoidFinder.h"

#include "check.h"
#include "fixtures.h"
#include "internal.h"

using namespace bitvf;
using namespace fixtures;

namespace {

  const char* const kFormats[2] = {".fits", ".txt"};
  constexpr double kPi = 3.14159265358979323846;

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  const otswap::Mask& full_mask ()
  {
    static const otswap::Mask m(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    return m;
  }

  otswap::MpsProfile profile ()
  {
    otswap::MpsProfile p;
    p.redshift = {0.1, 0.5};
    p.mps = {5., 6.};
    p.count = {40, 60};
    p.redshiftMin = 0.1;
    p.redshiftMax = 0.5;
    p.representative = 5.5;
    return p;
  }

  std::vector<double> objects (const std::size_t n, const unsigned seed)
  {
    std::mt19937 rng(seed);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < n; ++i) {
      sky.push_back(uniform(0.5, 1.2));
      sky.push_back(uniform(-0.3, 0.3));
      sky.push_back(uniform(0.25, 0.33));
    }
    return sky;
  }

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "bitvf_test_lightcone_files_" + name;
  }

  std::string bytes_of (const std::string& file)
  {
    std::ifstream in(file, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  }

  bool same (const double a, const double b)
  {
    return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0;
  }

  bool close (const double a, const double b)
  {
    return std::fabs(a - b) <= 4e-16 * std::max(1., std::fabs(a));
  }

  /// A search on 400 tracers displaced at random, every output kept.
  Result search (const unsigned nThreads = 1, const DivergenceGrid* grid = nullptr)
  {
    set_threads((int)nThreads);
    const std::size_t n = 400;
    const std::vector<double> sky = objects(n, 7);
    const std::vector<double> cart = otswap::toCartesian(sky, table());
    std::mt19937 rng(11);
    std::vector<double> displacement;
    for (std::size_t i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k)
        displacement.push_back(-0.02 * cart[3*i+k] + 30. * ((rng() / 4294967296.) - 0.5));
    Config config;
    config.verbosity = otswap::Verbosity::Silent;
    if (grid == nullptr) config.divergenceCellSize = 8.;
    config.keepDivergence = true;
    config.keepAssignment = true;
    config.nRealizationsMean = 1;
    otswap::RedshiftCut cut;
    cut.min = 0.26;
    return findVoidsLightcone(sky, displacement, {}, objects(1500, 8), full_mask(), table(), profile(), config, cut, {},
                              grid);
  }

  std::vector<std::string> names_of (const std::string& file)
  {
    return otswap::io::columnNames(file);
  }

  /// @p write must raise a bitvf::Error naming @p message, and leave no file.
  template <typename W>
  void refused (W&& write, const Result& r, const std::string& message, const std::string& what)
  {
    for (const char* format : kFormats) {
      const std::string file = temporary(std::string("refused") + format);
      std::remove(file.c_str());
      throws_naming([&] { write(file, r); }, message, what + " " + format);
      check(!std::ifstream(file).good(), what + " " + format + ": no file left");
      std::remove(file.c_str());
    }
  }

  void write_divergence (const std::string& f, const Result& r) { io::writeDivergence(f, r); }

  /// The value of keyword @p key of an ASCII @p file replaced by @p value; false when absent.
  bool set_keyword (const std::string& file, const std::string& key, const std::string& value)
  {
    std::string text = bytes_of(file);
    const std::string head = "## " + key + " = ";
    const std::size_t at = text.find(head);
    if (at == std::string::npos) return false;
    const std::size_t end = text.find(" /", at + head.size());
    text.replace(at + head.size(), end - at - head.size(), value);
    std::ofstream(file, std::ios::binary) << text;
    return true;
  }

  /// Field @p field (from 0) of the first data row of an ASCII @p file replaced by @p value.
  bool set_first_row (const std::string& file, const std::size_t field, const std::string& value)
  {
    const std::string text = bytes_of(file);
    std::size_t at = 0;
    while (at < text.size() && text[at] == '#') at = text.find('\n', at) + 1;
    const std::size_t end = text.find('\n', at);
    if (at >= text.size() || end == std::string::npos) return false;
    std::istringstream in(text.substr(at, end - at));
    std::vector<std::string> tokens;
    for (std::string t; in >> t;) tokens.push_back(t);
    if (field >= tokens.size()) return false;
    tokens[field] = value;
    std::string row;
    for (const std::string& t : tokens) row += (row.empty() ? "" : " ") + t;
    std::ofstream(file, std::ios::binary) << text.substr(0, at) + row + text.substr(end);
    return true;
  }

  /// The main table of an ASCII divergence @p from copied to a FITS @p to, keywords included, without the profile.
  void ascii_to_fits (const std::string& from, const std::string& to)
  {
    const std::vector<std::string> names = otswap::io::columnNames(from);
    const std::vector<std::string> integers = {"index", "nDispl", "nTracers", "voidIndex", "overdenseIndex"};
    std::vector<std::string> reals;
    for (const std::string& n : names)
      if (std::find(integers.begin(), integers.end(), n) == integers.end()) reals.push_back(n);
    const otswap::io::Table t = otswap::io::read(from, reals, integers);
    std::vector<otswap::io::Column> columns;
    for (const std::string& n : names) {
      otswap::io::Column c;
      c.name = n;
      const auto i = std::find(integers.begin(), integers.end(), n);
      if (i != integers.end()) {
        c.type = 'K';
        const std::size_t k = (std::size_t)(i - integers.begin());
        for (std::size_t r = 0; r < t.nRows; ++r) c.integers.push_back(t.integers[t.nIntegerColumns*r + k]);
      }
      else {
        const std::size_t k = (std::size_t)(std::find(reals.begin(), reals.end(), n) - reals.begin());
        c.type = n == "masked" ? 'J' : 'D';
        for (std::size_t r = 0; r < t.nRows; ++r) c.data.push_back(t.values[t.nColumns*r + k]);
      }
      columns.push_back(c);
    }
    otswap::io::WriteOptions options;
    bool started = false;
    for (const otswap::io::Keyword& k : otswap::io::readKeywords(from)) {
      if (k.name == "PRODUCT") started = true;
      if (started) options.keywords.push_back(k);
    }
    std::remove(to.c_str());
    otswap::io::write(to, columns, options);
  }

  /// The first occurrence of @p from in @p file replaced by @p to; false when absent.
  bool edit (const std::string& file, const std::string& from, const std::string& to)
  {
    std::string text = bytes_of(file);
    const std::size_t at = text.find(from);
    if (at == std::string::npos) return false;
    text.replace(at, from.size(), to);
    std::ofstream out(file, std::ios::binary);
    out << text;
    return true;
  }

}

int main ()
{
  const Result r = search();
  check(r.geometry == otswap::Geometry::Lightcone && r.divergence && r.voids.nBasins > 0 && r.overdense.nBasins > 0 &&
        r.divergence->nCells > 100 && r.divergence->nCells < 3000 && r.nRealizations == 1 &&
        r.divergence->nRealizations == 1 && r.divergence->nRealizations == r.nRealizations &&
        r.divergence->nRealizationsMean == 1 && r.redshiftCut.min == 0.26 && std::isinf(r.redshiftCut.max),
        "the fixture: a lightcone with basins of both kinds, its grid recording the realizations");
  std::cout << "  cells " << r.divergence->nCells << ", voids " << r.voids.nBasins << ", overdense "
            << r.overdense.nBasins << std::endl;

  group("columns: cellDist right after cellRed, the per-cell limits in degrees, masked last; the record keywords");
  for (const char* format : kFormats) {
    const std::string file = temporary(std::string("grid") + format);
    io::writeDivergence(file, r);
    const std::vector<std::string> expected = {
      "index", "cellX", "cellY", "cellZ", "cellRA", "cellDec", "cellRed", "cellDist", "cellRAMin", "cellRAMax",
      "cellDecMin", "cellDecMax", "cellDistMin", "cellDistMax", "divergence", "divergenceSmoothed", "volume",
      "weight", "weightSmoothed", "nDispl", "nTracers", "voidIndex", "overdenseIndex", "masked"};
    check(names_of(file) == expected, std::string("the columns and their order ") + format);

    const std::vector<otswap::io::Keyword> all = otswap::io::readKeywords(file);
    std::vector<std::string> names;
    bool started = false;
    for (const otswap::io::Keyword& k : all) {
      if (k.name == "PRODUCT") started = true;
      if (started) names.push_back(k.name);
    }
    const std::vector<std::string> keys = {"PRODUCT", "BITVFV", "GEOMETRY", "NCELLS", "NREC", "NRECMEAN", "SIGMA",
                                           "SIGUNIT", "CELLSIZE", "ZCUTMIN", "NRANDUSE", "MASKNSID", "DISTMIN",
                                           "DISTMAX", "RAMIN", "RAMAX", "DECMIN", "DECMAX"};
    check(names == keys, std::string("the keywords and their order ") + format);
    std::remove(file.c_str());
  }

  group("the profile part: PRODUCT divergence mps profile, no GEOMETRY, every column and keyword, 17 digits");
  for (const char* format : kFormats) {
    const std::string file = temporary(std::string("profile") + format);
    io::writeDivergence(file, r);
    const bool fits = std::string(format) == ".fits";
    const std::string part = fits ? file : otswap::io::twinFile(file, "mps");
    const std::vector<otswap::io::Keyword> k = otswap::io::readKeywords(part, fits ? "MPSPROF" : "");
    bool product = false, geometry = false, zmin = false, repr = false;
    for (const otswap::io::Keyword& w : k) {
      if (w.name == "PRODUCT") product = w.value == "divergence mps profile";
      if (w.name == "GEOMETRY") geometry = true;
      if (w.name == "ZMIN") zmin = std::strtod(w.value.c_str(), nullptr) == 0.1;
      if (w.name == "MPSREPR") repr = std::strtod(w.value.c_str(), nullptr) == 5.5;
    }
    check(product && !geometry && zmin && repr, std::string("the profile's keywords ") + format);
    const otswap::io::Table t = otswap::io::readPart(file, "MPSPROF", "mps", {"redshift", "MPS"}, {"nTracers"});
    check(t.nRows == 2 && t.values[0] == 0.1 && t.values[1] == 5. && t.values[2] == 0.5 && t.values[3] == 6. &&
          t.integers[0] == 40 && t.integers[1] == 60, std::string("the profile's rows ") + format);
    std::remove(file.c_str());
    if (!fits) std::remove(part.c_str());
  }

  group("readDivergence of a lightcone file: the record, the profile and the cells back; the angles to rounding");
  for (const char* format : kFormats) {
    const std::string file = temporary(std::string("back") + format);
    io::writeDivergence(file, r);
    const DivergenceGrid g = io::readDivergence(file);
    const DivergenceGrid& w = *r.divergence;
    check(g.geometry == otswap::Geometry::Lightcone && g.nCells == w.nCells && g.nRealizations == 1 &&
          g.nRealizationsMean == 1, std::string("geometry and counts ") + format);
    check(same(g.distanceMin, w.distanceMin) && same(g.distanceMax, w.distanceMax) && same(g.cellSizeMps, w.cellSizeMps),
          std::string("distances and cell size, bit for bit ") + format);
    check(close(g.raMin, w.raMin) && close(g.raMax, w.raMax) && close(g.decMin, w.decMin) && close(g.decMax, w.decMax),
          std::string("angles to rounding ") + format);
    check(same_bits(g.mpsProfile.redshift, w.mpsProfile.redshift) && same_bits(g.mpsProfile.mps, w.mpsProfile.mps) &&
          g.mpsProfile.count == w.mpsProfile.count && same(g.mpsProfile.redshiftMin, w.mpsProfile.redshiftMin) &&
          same(g.mpsProfile.representative, w.mpsProfile.representative), std::string("the profile ") + format);
    check(same_bits(g.centre, w.centre) && same_bits(g.divergence, w.divergence) &&
          same_bits(g.divergenceSmoothed, w.divergenceSmoothed) && same_bits(g.volume, w.volume) &&
          same_bits(g.weight, w.weight) && same_bits(g.weightSmoothed, w.weightSmoothed) &&
          g.nDisplacements == w.nDisplacements && g.nTracers == w.nTracers && g.voidIndex == w.voidIndex &&
          g.overdenseIndex == w.overdenseIndex && g.masked == w.masked, std::string("the cells ") + format);
    check(g.cellRaMin.empty() && g.cellDistanceMax.empty(), std::string("the per-cell limits are not read ") + format);
    std::remove(file.c_str());
    if (std::string(format) == ".txt") std::remove(otswap::io::twinFile(file, "mps").c_str());
  }

  group("a grid read from a lightcone file is a grid to reuse: the search on it keeps the cells and the basins");
  for (const char* format : kFormats) {
    const std::string file = temporary(std::string("reuse") + format);
    io::writeDivergence(file, r);
    const DivergenceGrid g = io::readDivergence(file);
    const Result again = search(1, &g);
    check(again.divergence && again.divergence->nCells == r.divergence->nCells, std::string("the same cells ") + format);
    bool near = again.divergence && again.divergence->voidIndex == r.divergence->voidIndex &&
                again.divergence->overdenseIndex == r.divergence->overdenseIndex;
    for (std::size_t c = 0; near && c < r.divergence->nCells; ++c)
      near = std::fabs(again.divergence->divergence[c] - r.divergence->divergence[c]) <=
             1e-12 * std::fabs(r.divergence->divergence[c]);
    check(near, std::string("the same basins, the divergence within the rounding of the angles ") + format);
    std::remove(file.c_str());
    if (std::string(format) == ".txt") std::remove(otswap::io::twinFile(file, "mps").c_str());
  }

  group("the grid kept in memory, reused, gives the same divergence and basins bit for bit");
  {
    const Result again = search(1, &*r.divergence);
    check(again.divergence && again.divergence->divergence == r.divergence->divergence &&
          again.voids.centre == r.voids.centre && again.overdense.centre == r.overdense.centre,
          "the same divergence and centres");
  }

  group("the same bytes with 1 and 4 threads, in both formats, the profile part included");
  {
    const Result four = search(4);
    for (const char* format : kFormats) {
      std::string bytes[2];
      for (int which = 0; which < 2; ++which) {
        const std::string file = temporary(std::string("threads") + format);
        io::writeDivergence(file, which == 0 ? r : four);
        bytes[which] = bytes_of(file);
        if (std::string(format) == ".txt") {
          bytes[which] += bytes_of(otswap::io::twinFile(file, "mps"));
          std::remove(otswap::io::twinFile(file, "mps").c_str());
        }
        std::remove(file.c_str());
      }
      check(!bytes[0].empty() && bytes[0] == bytes[1], std::string("identical files ") + format);
    }
    set_threads(1);
  }

  group("errors of the writer, each before the file is opened: the result, the grid, the grid against the result");
  {
    auto bad = [&r] (auto&& change, const std::string& message, const std::string& what) {
      Result b = r;
      change(b);
      refused(write_divergence, b, message, what);
    };
    bad([] (Result& b) { b.distances.reset(); }, "writeDivergence: the result is malformed: it holds no distance table",
        "no distance table");
    bad([] (Result& b) { b.mpsProfile = otswap::MpsProfile(); },
        "writeDivergence: the result is malformed: it holds no mps profile",
        "no mps profile");
    bad([] (Result& b) { b.maskNside = 0; },
        "writeDivergence: the result is malformed: maskNside is 0; it must lie in [1, 536870912]",
        "maskNside 0");
    bad([] (Result& b) { b.maskNside = 536870913; }, "maskNside is 536870913; it must lie in [1, 536870912]",
        "maskNside too large");
    bad([] (Result& b) { b.divergence->masked.pop_back(); }, "writeDivergence: cannot write divergence.masked: it holds",
        "a short masked column");
    bad([] (Result& b) { b.divergence->cellDistanceMax.push_back(1.); }, "cannot write divergence.cellDistanceMax",
        "a long limit column");
    bad([] (Result& b) { b.divergence->centreSky.pop_back(); }, "cannot write divergence.centreSky", "centreSky");
    bad([] (Result& b) { b.divergence->distanceMin = -1.; },
        "writeDivergence: the divergence grid has a distance range [-1, ", "a negative distance");
    bad([] (Result& b) { b.divergence->cellSizeMps = 0.; },
        "writeDivergence: the divergence grid has a cell size of 0 mps",
        "a zero cell size");
    bad([] (Result& b) { b.divergence->mps = 5.; },
        "writeDivergence: the divergence grid is a lightcone grid but holds the box field mps",
        "a box field");
    bad([] (Result& b) { b.divergence->mpsProfile.mps[0] = -1.; },
        "writeDivergence: the mps profile of the divergence grid holds a mean particle separation of -1",
        "a negative profile");
    bad([] (Result& b) { b.distances = std::make_shared<const otswap::DistanceTable>(0.9, 0.1, -1., 0., 0., 2., 1000); },
        " cells but the rebuild from its record gives",
        "another distance table");
    bad([] (Result& b) { b.distances = std::make_shared<const otswap::DistanceTable>(0.3, 0.7, -1., 0., 0., 0.2, 1000); },
        "writeDivergence: the divergence grid reaches distance ", "a distance table too short for the rebuild");
    bad([] (Result& b) { b.mpsProfile.mps[1] += 1.; },
        "the mps profile of its divergence grid differs from the search's at entry 1 (from 0): mps", "another profile");
    bad([] (Result& b) { b.config.divergenceCellSize = 9.; },
        "its divergence grid has a cell size of 8 mps, but its configuration records 9", "another cell size");
    bad([] (Result& b) { b.nRealizations = 2; }, "its divergence grid records 1 realizations, but the search's are 2",
        "another number of realizations");
    bad([] (Result& b) { b.config.nRealizationsMean = 2; },
        "its divergence grid records 1 realizations averaged, but its configuration records 2", "another mean");
    for (auto writer : {io::writeVoids, io::writeOverdense, io::writeGalaxies}) {
      Result b = r;
      b.maskNside = 0;
      refused([&writer] (const std::string& f, const Result& x) { writer(f, x); }, b,
              "the result is malformed: maskNside is 0; it must lie in [1, 536870912]", "maskNside in every writer");
    }
  }

  group("errors of the reader: the profile part missing or another product, a box keyword in a lightcone file and "
        "the reverse, a record that fails the checks of a grid to reuse, a missing or zero keyword, a bad row");
  {
    const std::string file = temporary("reader.txt");
    const std::string twin = otswap::io::twinFile(file, "mps");
    auto fresh = [&] {
      io::writeDivergence(file, r);
    };
    auto cleanup = [&] {
      std::remove(file.c_str());
      std::remove(twin.c_str());
    };
    fresh();
    std::remove(twin.c_str());
    throws_naming([&] { io::readDivergence(file); },
                  ": the mps profile part is missing; a lightcone divergence file holds it beside the grid",
                  "no profile part");
    cleanup();

    fresh();
    check(edit(twin, "divergence mps profile", "complete displacements mps profile"), "the profile's PRODUCT edited");
    throws_naming([&] { io::readDivergence(file); }, "is not a divergence mps profile file", "another profile product");
    cleanup();

    fresh();
    check(edit(file, "## NCELLS", "## MPS = 5\n## NCELLS"), "a box keyword added");
    throws_naming([&] { io::readDivergence(file); }, "keyword MPS is present, which a lightcone file does not hold",
                  "a box keyword in a lightcone file");
    cleanup();

    fresh();
    check(edit(file, "GEOMETRY = lightcone", "GEOMETRY = box"), "GEOMETRY edited");
    throws_naming([&] { io::readDivergence(file); }, "keyword DISTMIN is present, which a box file does not hold",
                  "a lightcone keyword in a box file");
    cleanup();

    fresh();
    check(edit(file, "DISTMIN = ", "DISTMIN = -5 / "), "DISTMIN edited");
    throws_naming([&] { io::readDivergence(file); }, ": the grid to reuse has a distance range [-5, ",
                  "a record failing the checks");
    cleanup();

    fresh();
    check(edit(file, "NRECMEAN = 1", "NRECMEAN = 0"), "NRECMEAN edited");
    throws_naming([&] { io::readDivergence(file); }, ": NRECMEAN is 0; it must be at least 1", "NRECMEAN 0");
    cleanup();

    fresh();
    check(edit(file, "## NREC = 1", "## NREC = 0"), "NREC edited");
    throws_naming([&] { io::readDivergence(file); }, ": NREC is 0; it must be at least 1", "NREC 0");
    cleanup();

    fresh();
    check(edit(file, "## DISTMAX = ", "## DISTMAY = "), "DISTMAX renamed");
    throws_naming([&] { io::readDivergence(file); }, ": keyword DISTMAX is missing", "DISTMAX missing");
    cleanup();

    fresh();
    check(set_first_row(file, names_of(file).size() - 1, "2"), "masked edited");
    throws_naming([&] { io::readDivergence(file); }, ": data row 1 (from 1) has masked 2; it must be 0 or 1",
                  "masked 2");
    cleanup();

    fresh();
    check(set_first_row(file, 0, "1"), "index edited");
    throws_naming([&] { io::readDivergence(file); },
                  ": data row 1 (from 1) has index 1; the rows must hold the cells 0 to NCELLS - 1 in order",
                  "a row out of order");
    cleanup();

    fresh();
    const std::string fits = temporary("reader.fits");
    ascii_to_fits(file, fits);
    throws_naming([&] { io::readDivergence(fits); },
                  ": the mps profile part is missing; a lightcone divergence file holds it beside the grid",
                  "no profile HDU in FITS");
    std::remove(fits.c_str());
    cleanup();
  }

  group("a right-ascension range of 360 degrees reads back within 2 pi");
  {
    const std::string file = temporary("ra.txt");
    io::writeDivergence(file, r);
    check(set_keyword(file, "RAMIN", "0.1") && set_keyword(file, "RAMAX", "360.1"), "RAMIN and RAMAX edited");
    const DivergenceGrid g = io::readDivergence(file);
    check(g.raMax - g.raMin <= 2. * kPi && close(g.raMin, 0.1 * kPi / 180.) && close(g.raMax, 360.1 * kPi / 180.),
          "the width kept within 2 pi, each bound within rounding");
    std::remove(file.c_str());
    std::remove(otswap::io::twinFile(file, "mps").c_str());
  }

  group("a grid read from a file is written again, the limits of its cells rebuilt from its record");
  for (const char* format : kFormats) {
    const std::string first = temporary(std::string("first") + format);
    const std::string second = temporary(std::string("second") + format);
    io::writeDivergence(first, r);
    Result b = r;
    b.divergence = io::readDivergence(first);
    check(b.divergence->cellRaMin.empty(), std::string("the limits are not read ") + format);
    io::writeDivergence(second, b);
    check(names_of(second) == names_of(first), std::string("the same columns ") + format);
    const std::vector<std::string> limits = {"cellDist", "cellRAMin", "cellRAMax", "cellDecMin", "cellDecMax",
                                             "cellDistMin", "cellDistMax"};
    const otswap::io::Table t1 = otswap::io::read(first, limits);
    const otswap::io::Table t2 = otswap::io::read(second, limits);
    bool near = t1.nRows == t2.nRows && t1.values.size() == t2.values.size();
    for (std::size_t i = 0; near && i < t1.values.size(); ++i)
      near = std::fabs(t1.values[i] - t2.values[i]) <= 1e-12 * std::max(1., std::fabs(t1.values[i]));
    check(near, std::string("the limits within the rounding of the angles ") + format);
    for (const std::string& f : {first, second}) {
      std::remove(f.c_str());
      if (std::string(format) == ".txt") std::remove(otswap::io::twinFile(f, "mps").c_str());
    }
  }

  group("lightconeGrid records no realizations; a search reusing its grid records its own");
  {
    const DivergenceGrid g = lightconeGrid(objects(400, 7), objects(1500, 8), full_mask(), table(), profile(), 8.,
                                           otswap::RedshiftCut{0.26, std::numeric_limits<double>::infinity()});
    check(g.nRealizations == 0 && g.nRealizationsMean == 0, "lightconeGrid gives 0 and 0");
    const Result again = search(1, &g);
    check(again.nRealizations == 1 && again.divergence && again.divergence->nRealizations == 1 &&
          again.divergence->nRealizationsMean == 1, "the kept grid records 1 realization, 1 averaged");
  }

  return report("test_lightcone_files");
}
