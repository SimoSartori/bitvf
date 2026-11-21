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
 *  @file tests/test_lightcone_pipeline.cpp
 *
 *  @brief The real-space pass of the lightcone pipeline: the tracers it
 *  moves, their flags and ids, the grid every search shares, the checks
 *  made before the first reconstruction and the report lines.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "otswap/Maths.h"
#include "otswap/RSD.h"

#include "check.h"
#include "internal.h"

using namespace bitvf;

namespace {

  constexpr std::uint64_t kRecordedHash = 0x89e3a8a3e04a4e10ull;
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

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

  std::vector<double> patch (const std::size_t n, const unsigned seed)
  {
    std::mt19937 rng(seed);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < n; ++i) {
      sky.push_back(uniform(0.5, 1.));
      sky.push_back(otswap::maths::asin(uniform(otswap::maths::sin(-0.25), otswap::maths::sin(0.25))));
      sky.push_back(uniform(0.25, 0.4));
    }
    return sky;
  }

  const std::vector<double>& tracers ()
  {
    static const std::vector<double> t = patch(500, 11);
    return t;
  }

  const std::vector<double>& randoms ()
  {
    static const std::vector<double> r = patch(1500, 12);
    return r;
  }

  std::vector<std::int64_t> ids ()
  {
    std::vector<std::int64_t> v(500);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = 1000 + 3 * (std::int64_t)i;
    return v;
  }

  LightconeConfig base_config (std::ostream* log)
  {
    LightconeConfig c;
    c.ot.nRealizations = 2;
    c.ot.seed = 7;
    c.ot.verbosity = otswap::Verbosity::Silent;
    c.nBins = 1;
    c.cut = {0.26, 0.39};
    LightconeRsd rsd;
    rsd.bias = {{0., 1.}, {1.5, 1.5}};
    rsd.sigma = 20.;
    rsd.correction.verbosity = otswap::Verbosity::Silent;
    rsd.secondPassRandoms = SecondPassRandoms::None;
    c.rsd = rsd;
    c.voids.divergenceCellSize = 0.5;
    c.voids.verbosity = otswap::Verbosity::Detailed;
    c.voids.log = log;
    return c;
  }

  struct Chain {
    otswap::Result first;
    DivergenceGrid grid;
    otswap::RealSpaceCatalog real;
    internal::LightconeSecondPass pass;
  };

  Chain chain (const LightconeConfig& c)
  {
    Chain out;
    out.first = otswap::reconstructLightcone(tracers(), randoms(), full_mask(), c.nBins, table(), c.ot, c.cut);
    out.first.ids = ids();
    out.grid = lightconeGrid(out.first, randoms(), full_mask(), c.voids.divergenceCellSize);
    out.real = otswap::realSpaceLightcone(out.first, table(), c.rsd->bias, c.rsd->sigma, c.rsd->correction);
    const internal::LightconeGrid lc = internal::lightcone_grid(out.grid, table(), full_mask(), true);
    out.pass = internal::lightcone_second_pass(out.real, lc, table(), c.cut, ids());
    return out;
  }

  std::string run_log (const LightconeConfig& c, std::string& message)
  {
    std::ostringstream log;
    LightconeConfig cc = c;
    cc.voids.log = &log;
    try {
      runLightcone(tracers(), randoms(), full_mask(), table(), cc, {}, ids());
      message = "returned";
    }
    catch (const Error& e) {
      message = e.what();
    }
    return log.str();
  }

  std::string grid_line (const std::string& log, const std::size_t from = 0)
  {
    const std::size_t a = log.find("grid of ", from);
    if (a == std::string::npos) return "";
    const std::size_t b = log.find("(from the reused grid)", a);
    if (b == std::string::npos) return "";
    return log.substr(a, b - a);
  }

  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  std::uint64_t pass_hash ()
  {
    const Chain c = chain(base_config(nullptr));
    const internal::LightconeSecondPass& p = c.pass;
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const std::size_t r : p.rows) {
      const std::uint64_t v = r;
      hash_bytes(h, &v, sizeof v);
    }
    hash_bytes(h, p.tracersSky.data(), p.tracersSky.size() * sizeof(double));
    hash_bytes(h, p.ids.data(), p.ids.size() * sizeof(std::int64_t));
    hash_bytes(h, p.uncorrected.data(), p.uncorrected.size());
    hash_bytes(h, p.outsideCut.data(), p.outsideCut.size());
    hash_bytes(h, p.outsideGrid.data(), p.outsideGrid.size());
    return h;
  }

}

int main ()
{
  // The randoms kept are those on observed pixels inside the cut; a random outside both is counted once, under the mask.
  group("the selection of the pipeline's inputs: kept redshifts and randoms, with a partial mask and a cut");
  {
    std::vector<double> pixels(12, 1.);
    pixels[4] = 0.;
    const otswap::Mask partial(pixels, otswap::PixelOrdering::Ring);
    check(!partial.allows(0.2, 0.1) && partial.allows(4., 0.1), "the test mask hides RA 0.2 and shows RA 4");
    const std::vector<double> tSky = {0.2, 0.1, 0.3,   4., 0.1, 0.3,   4., 0.1, 0.25};
    const std::vector<double> rSky = {0.2, 0.1, 0.3,   0.2, 0.1, 0.5,   4., 0.1, 0.3,
                                      4., 0.1, 0.3,   4., 0.1, 0.5,   4., 0.1, 0.1};
    const internal::PipelineSelection sel =
      internal::check_pipeline_sky(tSky, rSky, partial, table(), otswap::RedshiftCut{0.2, 0.4}, {}, "test");
    check(sel.redshift == std::vector<double>({0.3, 0.25}), "the kept tracers' redshifts are those on observed pixels");
    check(sel.nRandomsKept == 2, "2 of 6 randoms are kept: 2 masked (one also outside the cut), 2 outside the cut");
    const internal::PipelineSelection open =
      internal::check_pipeline_sky(tSky, rSky, full_mask(), table(), otswap::RedshiftCut{}, {}, "test");
    check(open.redshift.size() == 3 && open.nRandomsKept == 6, "without a cut and with a full mask, all are kept");
    check(sel.randomRows == std::vector<std::size_t>({2, 3}) &&
          open.randomRows == std::vector<std::size_t>({0, 1, 2, 3, 4, 5}),
          "the rows of the randoms kept are listed in order");
  }

  // The resampled randoms keep each kept random's sky position and take a redshift of the tracers given, from a serial draw of the seed alone.
  group("the resampled randoms");
  {
    const std::vector<double> rSky = {0.1, 0.01, 0.30,   0.2, 0.02, 0.31,   0.3, 0.03, 0.32,
                                      0.4, 0.04, 0.33,   0.5, 0.05, 0.34};
    const std::vector<double> tSky = {9., 9., 0.25,   9., 9., 0.26,   9., 9., 0.27};
    const std::vector<std::size_t> rows = {0, 2, 3, 4};
    const std::vector<double> a = internal::resample_randoms(rSky, rows, tSky, 104736u);
    check(a.size() == 12, "3 values per row chosen");
    bool sky = true, redshift = true;
    for (std::size_t j = 0; j < rows.size(); ++j) {
      sky = sky && std::memcmp(&a[3*j], &rSky[3*rows[j]], 2 * sizeof(double)) == 0;
      redshift = redshift && (a[3*j+2] == 0.25 || a[3*j+2] == 0.26 || a[3*j+2] == 0.27);
    }
    check(sky, "right ascension and declination are those of the row, bit for bit");
    check(redshift, "every redshift is one of the tracers'");

    std::mt19937 engine(104736u);
    bool drawn = true;
    for (std::size_t j = 0; j < rows.size(); ++j) {
      std::uint64_t product = (std::uint64_t)engine() * 3u;
      while ((std::uint32_t)product < 3u && (std::uint32_t)product < (0u - 3u) % 3u)
        product = (std::uint64_t)engine() * 3u;
      drawn = drawn && a[3*j+2] == tSky[3 * (product >> 32) + 2];
    }
    check(drawn, "one multiply-and-reject draw per row, in order, from an engine seeded with the seed");
    check(a == internal::resample_randoms(rSky, rows, tSky, 104736u), "the same seed gives the same randoms");
    bool differs = false;
    for (unsigned seed = 1; seed < 20 && !differs; ++seed) differs = internal::resample_randoms(rSky, rows, tSky, seed) != a;
    check(differs, "another seed gives other redshifts");
  }

  // A corrected tracer keeps its sky position bit for bit; excluded rows carry one flag, uncorrected before outsideCut before outsideGrid.
  group("the tracers chosen, their flags in order, and their ids");
  {
    otswap::MpsProfile flat;
    flat.redshift = {0., 0.5};
    flat.mps = {40., 40.};
    const DivergenceGrid record = lightconeGrid(tracers(), randoms(), full_mask(), table(), flat, 1.);
    const internal::LightconeGrid lc = internal::lightcone_grid(record, table(), full_mask(), true);
    check(lc.cells.nCells <= 3000, "the grid is small");

    otswap::RealSpaceCatalog real;
    real.nObjects = 6;
    real.geometry = otswap::Geometry::Lightcone;
    using S = otswap::CorrectionStatus;
    real.status = {S::Corrected, S::NoValidNeighbour, S::LeftOut, S::Corrected, S::MovedByNeighbours,
                   S::MovedByNeighbours};
    real.sky = {0.7, 0.01, 0.3,   kNaN, kNaN, kNaN,   kNaN, kNaN, kNaN,
                0.7, 0.0, 0.45,   0.8, 0.1, 0.41,   0.9, -0.1, 0.33};
    const otswap::RedshiftCut cut {0.2, 0.42};

    const std::vector<std::int64_t> given {10, 11, 12, 13, 14, 15};
    const internal::LightconeSecondPass p = internal::lightcone_second_pass(real, lc, table(), cut, given);
    check(p.rows == std::vector<std::size_t>({0, 5}), "rows 0 and 5 are chosen");
    check(p.ids == std::vector<std::int64_t>({10, 15}), "the ids of the chosen rows are those given");
    check(p.tracersSky.size() == 6 && std::memcmp(p.tracersSky.data(), &real.sky[0], 3 * sizeof(double)) == 0 &&
          std::memcmp(p.tracersSky.data() + 3, &real.sky[15], 3 * sizeof(double)) == 0,
          "the chosen tracers keep the corrected sky position bit for bit");
    check(p.uncorrected == std::vector<std::uint8_t>({0, 1, 1, 0, 0, 0}) &&
          p.outsideCut == std::vector<std::uint8_t>({0, 0, 0, 1, 0, 0}) &&
          p.outsideGrid == std::vector<std::uint8_t>({0, 0, 0, 0, 1, 0}),
          "one flag per excluded row: a tracer beyond both the cut and the grid is outside the cut");
    check(p.nUncorrected == 2 && p.nOutsideCut == 1 && p.nOutsideGrid == 1, "the counts match the flags");

    const internal::LightconeSecondPass q = internal::lightcone_second_pass(real, lc, table(), cut, {});
    check(q.ids == std::vector<std::int64_t>({0, 5}), "without ids, the ids are the rows");

    otswap::RealSpaceCatalog none = real;
    none.status = {S::LeftOut, S::LeftOut, S::NoValidNeighbour, S::Corrected, S::Corrected, S::LeftOut};
    throws_naming([&] { internal::lightcone_second_pass(none, lc, table(), cut, {}); },
                  "no tracer is left for the real-space reconstruction: of 6 tracers, 4 were left uncorrected, "
                  "1 were moved outside the redshift cut and 1 were moved outside the first pass's grid",
                  "no tracer chosen is an error giving the three counts");
  }

  // Invalid inputs and an unusable correction fail before the first reconstruction does any work.
  group("the checks before the first reconstruction");
  {
    std::ostringstream otLog;
    LightconeConfig c = base_config(nullptr);
    c.ot.verbosity = otswap::Verbosity::Normal;
    c.ot.log = &otLog;
    auto run = [&] (const LightconeConfig& cc, const std::vector<double>& second = {}) {
      runLightcone(tracers(), randoms(), full_mask(), table(), cc, second, ids());
    };

    const otswap::DistanceTable noGrowth(table().redshift(), table().distance());
    throws_naming([&] { runLightcone(tracers(), randoms(), full_mask(), noGrowth, c, {}, ids()); }, "otswap: ",
                  "a table without growth rate fails");
    LightconeConfig negative = c;
    negative.rsd->bias = {{0., 1.}, {1.5, -1.}};
    throws_naming([&] { run(negative); }, "otswap: ", "a bias that is not positive fails");
    LightconeConfig noBias = c;
    noBias.rsd->bias = {};
    throws_naming([&] { run(noBias); }, "otswap: ", "a missing bias table fails");
    LightconeConfig sigma = c;
    sigma.rsd->sigma = kNaN;
    throws_naming([&] { run(sigma); }, "the smoothing width is nan in coordinate units; it must be finite and non-negative",
                  "an unset smoothing width fails");
    LightconeConfig bins = c;
    bins.nBins = 0;
    throws_naming([&] { run(bins); }, "the number of bins is 0; it must be at least 1", "0 bins fail");
    throws_naming([&] { run(c, randoms()); }, "but their mode is None or Resample",
                  "second-pass randoms given with None fail");
    LightconeConfig external = c;
    external.rsd->secondPassRandoms = SecondPassRandoms::External;
    throws_naming([&] { run(external); }, "the second-pass randoms are External but none were given",
                  "External without randoms fails");
    LightconeConfig plain = c;
    plain.rsd.reset();
    throws_naming([&] { run(plain, randoms()); }, "no correction of redshift-space distortions is requested",
                  "second-pass randoms without a correction fail");
    std::vector<double> farRandoms = randoms();
    for (std::size_t i = 2; i < farRandoms.size(); i += 3) farRandoms[i] = 0.5;
    throws_naming([&] { runLightcone(tracers(), farRandoms, full_mask(), table(), c, {}, ids()); },
                  "no random is left: of 1500 randoms, 0 lie on unobserved pixels and 1500 outside the redshift cut",
                  "no random inside the cut fails");
    check(otLog.str().empty(), "the first reconstruction never ran");
  }

  // The redshift-space search and the real-space search both reuse the grid of the first reconstruction.
  group("one grid shared by every search");
  {
    LightconeConfig c = base_config(nullptr);
    std::string zMessage, rMessage, refMessage;
    LightconeConfig z = c;
    z.rsd->zspaceVoids = true;
    const std::string zLog = run_log(z, zMessage);
    const std::string rLog = run_log(c, rMessage);
    check(zMessage.find("findVoidsLightcone") != std::string::npos &&
          rMessage.find("findVoidsLightcone") != std::string::npos,
          "both runs reach a search, which stops at the steps not written yet");
    const Chain ref = chain(c);
    std::ostringstream refLog;
    Config v = c.voids;
    v.log = &refLog;
    v.divergenceCellSize.reset();
    try {
      findVoidsLightcone(ref.first, randoms(), full_mask(), v, &ref.grid);
    }
    catch (const Error& e) {
      refMessage = e.what();
    }
    const std::string expected = grid_line(refLog.str());
    check(!expected.empty() && grid_line(zLog) == expected,
          "the redshift-space search is the search of the first reconstruction on its own grid");
    check(grid_line(rLog) == expected, "the real-space search reuses that grid");
  }

  // The tracers chosen, their ids and the report lines are those of the corrected catalogue on the first pass's grid.
  group("the second pass's tracers, ids and report lines");
  {
    const LightconeConfig c = base_config(nullptr);
    const Chain ref = chain(c);
    const internal::LightconeSecondPass& p = ref.pass;
    bool idsMatch = p.ids.size() == p.rows.size();
    const std::vector<std::int64_t> given = ids();
    for (std::size_t j = 0; idsMatch && j < p.rows.size(); ++j) idsMatch = p.ids[j] == given[p.rows[j]];
    check(idsMatch, "the second reconstruction's ids are the ids given for the rows chosen");
    check(p.nUncorrected > 0, "the tracers outside the cut are uncorrected");

    std::string message;
    const std::string log = run_log(c, message);
    std::ostringstream line;
    line << "second pass on " << p.rows.size() << " of 500 tracers: " << p.nUncorrected << " uncorrected, "
         << p.nOutsideCut << " outside the redshift cut, " << p.nOutsideGrid << " outside the first pass's grid in ";
    check(log.find(line.str()) != std::string::npos, "the second-pass line gives the three counts");
    std::size_t kept = 0;
    for (std::size_t i = 0; i < randoms().size() / 3; ++i)
      if (randoms()[3*i+2] >= c.cut.min && randoms()[3*i+2] <= c.cut.max) ++kept;
    check(kept > 0 && kept < 1500, "the cut removes some randoms");
    const otswap::SelectionCounts& sel = ref.first.selection;
    check(kept == sel.randoms - sel.randomsOutsideRedshiftCut - sel.randomsOutsideMask + sel.randomsOutsideBoth,
          "the randoms kept are those the reconstruction kept");
    check(log.find("bitvf: second-pass randoms: the first pass's " + std::to_string(kept) + "\n") != std::string::npos,
          "the randoms line gives the randoms kept by the mask and the cut");
    check(log.find(line.str()) < log.find("second-pass randoms"), "the second-pass line comes first");

    LightconeConfig resample = c;
    resample.rsd->secondPassRandoms = SecondPassRandoms::Resample;
    const std::string rLog = run_log(resample, message);
    check(rLog.find(line.str()) != std::string::npos && rLog.find("second-pass randoms") == std::string::npos,
          "the tracers are chosen the same way whatever the randoms");
    check(rLog.find("bitvf: resampled the redshifts of " + std::to_string(kept) + " randoms from " +
                    std::to_string(p.rows.size()) + " tracers with seed " + std::to_string(c.ot.seed + 104729u) + "\n") !=
          std::string::npos, "the Resample line gives the randoms kept, the tracers and the seed");
    check(message.find("findVoidsLightcone") != std::string::npos,
          "a Resample run reaches the second search");

    auto seed_printed = [&] (const unsigned seed, unsigned& value) {
      LightconeConfig cs = resample;
      cs.ot.seed = seed;
      const std::string text = run_log(cs, message);
      const std::size_t at = text.find(" with seed ");
      if (at == std::string::npos) return false;
      value = (unsigned)std::stoul(text.substr(at + 11));
      return true;
    };
    unsigned printed = 0;
    check(seed_printed(4294967295u, printed) && printed == 104728u,
          "the seed is the reconstruction's plus 104729 modulo 2^32, as printed");
    check(seed_printed(0u, printed) && (unsigned)(printed - 104729u) != 0,
          "a seed of 0 is replaced by a non-zero one before the offset is added");
  }

  // The tracers chosen are identical bit for bit with any number of threads and match the recorded hash.
  group("the second pass is the same with 1 and 4 threads, and matches the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = pass_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = pass_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "second pass hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_pipeline");
}
