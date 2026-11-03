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
 *  @file tests/test_runbox.cpp
 *
 *  @brief The box pipeline, runBox: equal to its steps run by hand, the
 *  randoms it draws, the rows of the real-space pass, the identities, the
 *  kept outputs, its errors, its report, and hashes with 1 and 4 threads.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "otswap/RSD.h"

#include "bitvf/BitVF.h"

#include "check.h"
#include "fixtures.h"
#include "internal.h"

using namespace bitvf;
using namespace fixtures;

namespace {

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

  otswap::BoxLimits box (const double lo, const double hi)
  {
    otswap::BoxLimits b;
    b.lower = {lo, lo, lo};
    b.upper = {hi, hi, hi};
    return b;
  }

  /// 2000 tracers uniform in [0, 100)^3 and 4 beyond x = 100; two
  /// realizations, seed 1; IDs given; sigma 0.5 mps; assignment and
  /// divergence kept; the correction along x at z = 0.5 with a small bias,
  /// so that it moves some tracers out of the box and of the grid.
  struct Pipeline {

    static constexpr std::size_t kN = 2000;
    static constexpr std::size_t kOutside = 4;

    std::vector<double> tracers;
    std::vector<std::int64_t> ids;
    BoxConfig config;

    Pipeline ()
    {
      std::mt19937 rng(2026);
      for (std::size_t i = 0; i < 3 * kN; ++i) tracers.push_back(100. * (rng() / 4294967296.));
      for (std::size_t i = 0; i < kOutside; ++i) {
        const double x[3] = {110. + (double)i, 50., 50.};
        tracers.insert(tracers.end(), x, x + 3);
      }
      for (std::size_t i = 0; i < kN + kOutside; ++i) ids.push_back(i % 2 == 0 ? (std::int64_t)(700 + i) : -(std::int64_t)(700 + i));
      config.ot.nRealizations = 2;
      config.ot.seed = 1;
      config.ot.log = nullptr;
      config.voids.smoothing.sigma = 0.5;
      config.voids.keepAssignment = true;
      config.voids.keepDivergence = true;
      config.voids.log = nullptr;
    }

    BoxRsd rsd () const
    {
      BoxRsd r;
      r.axis = 0;
      r.redshift = 0.5;
      r.bias = 0.02;
      r.sigma = 0.5;
      r.correction.log = nullptr;
      return r;
    }

    /// The tracers inside [0, 100]^3 only, for runs without limits.
    std::vector<double> inner () const
    {
      return std::vector<double>(tracers.begin(), tracers.begin() + 3 * kN);
    }
  };

  /// The randoms runBox draws, written from the formula of its header.
  std::vector<double> header_randoms (const std::size_t n, const otswap::BoxLimits& b, const unsigned seed)
  {
    std::mt19937 rng(seed);
    std::vector<double> r(3 * n);
    for (std::size_t j = 0; j < r.size(); ++j) {
      const std::uint64_t a = rng() >> 5;
      const std::uint64_t c = rng() >> 6;
      const double u = (double)(a * 67108864ull + c) / 9007199254740992.;
      r[j] = b.lower[j%3] + u * (b.upper[j%3] - b.lower[j%3]);
      if (r[j] >= b.upper[j%3]) r[j] = std::nextafter(b.upper[j%3], b.lower[j%3]);
    }
    return r;
  }

  otswap::BoxLimits bounding_box (const std::vector<double>& a)
  {
    otswap::BoxLimits b = box(std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
    for (std::size_t j = 0; j < a.size(); ++j) {
      b.lower[j%3] = std::min(b.lower[j%3], a[j]);
      b.upper[j%3] = std::max(b.upper[j%3], a[j]);
    }
    return b;
  }

  /// Every field of a search but elapsedSeconds and config's log.
  std::uint64_t result_hash (const Result& r)
  {
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const BasinCatalog* c : {&r.voids, &r.overdense}) {
      fnv(h, &c->nBasins, sizeof(c->nBasins));
      fnv_vector(h, c->centre);
      fnv_vector(h, c->radius);
      fnv_vector(h, c->radiusExt);
      fnv_vector(h, c->deltaAv);
      fnv_vector(h, c->delta2Mps);
      fnv_vector(h, c->edgeFlag);
      fnv_vector(h, c->contamination);
      fnv_vector(h, c->nTracers);
    }
    fnv(h, &r.nTracers, sizeof(r.nTracers));
    fnv(h, &r.nRealizations, sizeof(r.nRealizations));
    fnv(h, &r.mps, sizeof(r.mps));
    const double cell = r.config.divergenceCellSize.value_or(-1.);
    fnv(h, &cell, sizeof(cell));
    const unsigned seed = r.randomSeed.value_or(0);
    fnv(h, &seed, sizeof(seed));
    fnv_vector(h, r.ids);
    fnv_vector(h, r.uncorrected);
    fnv_vector(h, r.outsideGrid);
    if (r.limits) {
      fnv(h, r.limits->lower.data(), sizeof(r.limits->lower));
      fnv(h, r.limits->upper.data(), sizeof(r.limits->upper));
    }
    const std::size_t counts[] = {r.counts.tracers, r.counts.tracersUsed, r.counts.tracersUncorrected,
                                  r.counts.tracersOutsideBox, r.counts.tracersOutsideGrid, r.counts.displacements};
    fnv(h, counts, sizeof(counts));
    if (r.divergence) {
      const DivergenceGrid& g = *r.divergence;
      fnv(h, g.lower.data(), sizeof(g.lower));
      fnv(h, &g.cellSize, sizeof(g.cellSize));
      fnv(h, &g.nCells, sizeof(g.nCells));
      fnv_vector(h, g.divergence);
      fnv_vector(h, g.weight);
      fnv_vector(h, g.nTracers);
      fnv_vector(h, g.voidIndex);
      fnv_vector(h, g.overdenseIndex);
    }
    if (r.assignment) {
      const TracerAssignment& a = *r.assignment;
      fnv_vector(h, a.cellIndex);
      fnv_vector(h, a.voidIndex);
      fnv_vector(h, a.overdenseIndex);
      fnv_vector(h, a.distance);
      fnv_vector(h, a.inMask);
      fnv_vector(h, a.inBox);
      fnv_vector(h, a.used);
    }
    return h;
  }

  /// A search with randomSeed cleared, for a comparison with a search run
  /// by hand.
  std::uint64_t hash_without_seed (Result r)
  {
    r.randomSeed.reset();
    return result_hash(r);
  }

  bool same_reconstruction (const otswap::Result& a, const otswap::Result& b)
  {
    return a.nObjects == b.nObjects && a.config.seed == b.config.seed && same_bits(a.tracers, b.tracers) &&
           same_bits(a.displacement, b.displacement) && a.valid == b.valid && a.ids == b.ids &&
           a.mps == b.mps && a.outsideBox == b.outsideBox;
  }

  /// The error's nested exception is an otswap::Error.
  template <typename F>
  bool nests_otswap (F&& f)
  {
    try {
      f();
    }
    catch (const Error& e) {
      try {
        std::rethrow_if_nested(e);
      }
      catch (const otswap::Error&) {
        return true;
      }
      catch (...) {
        return false;
      }
    }
    return false;
  }

  std::string masked (const std::string& text)
  {
    return std::regex_replace(text, std::regex("in [0-9]+\\.[0-9]{2} s"), "in T s");
  }

}

int main ()
{
  set_threads(4);
  const otswap::DistanceTable distances(0.3, 0.7);
  const Pipeline in;
  const unsigned derived = 1u + 1299709u;

  group("without correction: equal to otswap::reconstructBox then findVoidsBox, randoms drawn and given");
  {
    const std::vector<double> t = in.inner();
    const otswap::BoxLimits b = bounding_box(t);
    BoxConfig c = in.config;
    c.keepReconstructions = true;
    const BitVFResult run = runBox(t, {}, c, std::vector<std::int64_t>(in.ids.begin(), in.ids.begin() + Pipeline::kN));

    const std::vector<double> randoms = header_randoms(2 * Pipeline::kN, b, derived);
    otswap::Result rec = otswap::reconstructBox(t, randoms, c.ot, b);
    rec.ids.assign(in.ids.begin(), in.ids.begin() + Pipeline::kN);
    const Result byHand = findVoidsBox(rec, c.voids);
    check(run.voids.randomSeed == derived, "randomSeed is the seed plus 1299709");
    check(hash_without_seed(run.voids) == result_hash(byHand), "drawn: the voids equal the steps by hand");
    check(run.reconstruction && same_reconstruction(*run.reconstruction, rec),
          "drawn: the reconstruction equals reconstructBox on the header's randoms in the bounding box");
    check(run.reconstruction->box->lower == b.lower && run.reconstruction->box->upper == b.upper &&
          run.reconstruction->selection.boxLimits && run.voids.limits, "the bounding box passed as limits");
    check(!run.zspaceVoids && !run.zspaceReconstruction && !run.realSpace, "no pass of the correction");

    const BitVFResult given = runBox(t, randoms, c, rec.ids);
    check(!given.voids.randomSeed && result_hash(given.voids) == result_hash(byHand),
          "given the header's randoms: identical, without randomSeed");

    const otswap::Result plain = otswap::reconstructBox(t, randoms, c.ot);
    check(plain.mps == rec.mps, "the mps equals otswap's without limits, bit for bit");

    BoxConfig limited = in.config;
    limited.limits = box(0., 100.);
    limited.mps = 7.5;
    const BitVFResult l = runBox(in.tracers, randoms, limited, in.ids);
    otswap::Result lrec = otswap::reconstructBox(in.tracers, randoms, 7.5, limited.ot, *limited.limits);
    lrec.ids = in.ids;
    check(result_hash(l.voids) == result_hash(findVoidsBox(lrec, limited.voids)) && l.voids.mps == 7.5 &&
          l.voids.counts.tracersOutsideBox == Pipeline::kOutside, "given limits, randoms and mps: equal by hand");
  }

  group("ot.seed 0: one seed drawn for both reconstructions and the randoms");
  {
    BoxConfig c = in.config;
    c.ot.seed = 0;
    c.limits = box(0., 100.);
    c.rsd = in.rsd();
    c.keepReconstructions = true;
    const BitVFResult r = runBox(in.tracers, {}, distances, c, in.ids);
    const unsigned s = r.zspaceReconstruction->config.seed;
    check(s != 0 && r.reconstruction->config.seed == s && r.voids.randomSeed == s + 1299709u,
          "the drawn seed recorded in both reconstructions and in randomSeed");
    BoxConfig again = c;
    again.ot.seed = s;
    check(result_hash(runBox(in.tracers, {}, distances, again, in.ids).voids) == result_hash(r.voids),
          "the run reproduced from the recorded seed");
  }

  BoxConfig rsdConfig = in.config;
  rsdConfig.limits = box(0., 100.);
  rsdConfig.rsd = in.rsd();
  rsdConfig.rsd->zspaceVoids = true;
  rsdConfig.keepReconstructions = true;
  const BitVFResult full = runBox(in.tracers, {}, distances, rsdConfig, in.ids);
  const otswap::Result& first = *full.zspaceReconstruction;
  const otswap::Result& second = *full.reconstruction;
  const otswap::RealSpaceCatalog& real = *full.realSpace;
  const std::size_t N = Pipeline::kN + Pipeline::kOutside;

  group("with correction: the redshift-space pass and the grid of every search");
  {
    check(hash_without_seed(*full.zspaceVoids) == result_hash(findVoidsBox(first, rsdConfig.voids)) &&
          full.zspaceVoids->randomSeed == derived, "zspaceVoids is findVoidsBox on zspaceReconstruction");
    const DivergenceGrid& a = *full.zspaceVoids->divergence;
    const DivergenceGrid& b = *full.voids.divergence;
    const DivergenceGrid g = boxGrid(first, rsdConfig.voids.divergenceCellSize);
    check(a.lower == b.lower && a.upper == b.upper && a.nCellsAxis == b.nCellsAxis && a.cellSize == b.cellSize &&
          a.mps == b.mps && b.lower == g.lower && b.cellSize == g.cellSize && b.nCells == g.nCells,
          "the second search's grid is the first pass's");
    check(full.voids.mps == first.mps && full.voids.config.divergenceCellSize == g.cellSizeMps,
          "the second search's mps and cell size are the first pass's");
    check(second.mps != first.mps, "the second reconstruction's own mps differs (fewer tracers)");
    check(second.config.seed == first.config.seed && same_bits(first.tracers, in.tracers),
          "the same seed on both reconstructions; the first on the input");
  }

  group("with correction: rows aligned, flags in order, counts over the input rows");
  {
    const Result& v = full.voids;
    const TracerAssignment& a = *v.assignment;
    const DivergenceGrid grid = boxGrid(first, rsdConfig.voids.divergenceCellSize);
    const Result byHand = findVoidsBox(second, rsdConfig.voids, &grid);
    const TracerAssignment& h = *byHand.assignment;

    std::size_t u = 0, g = 0, j = 0;
    bool excludedRight = true, selectedRight = true, flagsRight = v.uncorrected.size() == N && v.outsideGrid.size() == N;
    for (std::size_t i = 0; i < N; ++i) {
      const bool moved = real.status[i] == otswap::CorrectionStatus::Corrected ||
                         real.status[i] == otswap::CorrectionStatus::MovedByNeighbours;
      const bool inGrid = moved && bitvf::internal::box_cell(grid, &real.cartesian[3*i]) >= 0;
      if (!moved) ++u;
      else if (!inGrid) ++g;
      flagsRight = flagsRight && v.uncorrected[i] == !moved && v.outsideGrid[i] == (moved && !inGrid);
      if (!inGrid) {
        excludedRight = excludedRight && a.cellIndex[i] == -1 && a.voidIndex[i] == -1 && a.overdenseIndex[i] == -1 &&
                        std::isnan(a.distance[i]) && a.inMask[i] == 1 && a.inBox[i] == 1 && a.used[i] == 0;
        continue;
      }
      selectedRight = selectedRight && std::memcmp(&second.tracers[3*j], &real.cartesian[3*i], 3 * sizeof(double)) == 0 &&
                      second.ids[j] == in.ids[i] && a.cellIndex[i] == h.cellIndex[j] && a.voidIndex[i] == h.voidIndex[j] &&
                      a.overdenseIndex[i] == h.overdenseIndex[j] &&
                      std::memcmp(&a.distance[i], &h.distance[j], sizeof(double)) == 0 && a.inMask[i] == h.inMask[j] &&
                      a.inBox[i] == h.inBox[j] && a.used[i] == h.used[j];
      ++j;
    }
    check(u >= Pipeline::kOutside && g > 0 && byHand.counts.tracersOutsideBox > 0,
          "the fixture has tracers uncorrected, moved out of the grid, and moved out of the box into the grid");
    check(second.nObjects == j && j == N - u - g, "the second reconstruction holds the selected rows, in input order");
    check(flagsRight, "uncorrected for the unmoved (LeftOut included), else outsideGrid for those outside the grid");
    check(excludedRight, "excluded rows: -1, NaN, inMask and inBox 1, used 0");
    check(selectedRight, "selected rows: the corrected position, the input ID, the second search's assignment");
    check(same_catalog(v.voids, byHand.voids) && same_catalog(v.overdense, byHand.overdense),
          "the catalogues are the second search's on the first pass's grid");

    std::size_t inBoxZeros = 0;
    for (const std::uint8_t b : a.inBox) inBoxZeros += b == 0;
    check(v.nTracers == N && v.ids == in.ids && v.counts.tracers == N && v.counts.tracersUncorrected == u &&
          v.counts.tracersOutsideGrid == g && v.counts.tracersUsed == byHand.counts.tracersUsed &&
          v.counts.tracersOutsideBox == byHand.counts.tracersOutsideBox &&
          v.counts.displacements == byHand.counts.displacements && inBoxZeros == v.counts.tracersOutsideBox,
          "nTracers, ids, counts; the inBox zeros equal tracersOutsideBox");
    check(v.outsideCut.empty() && v.randomSeed == derived, "no outsideCut in a box; randomSeed set");
  }

  group("identities: the input's in the first reconstruction, the corrected catalogue and the searches");
  {
    check(first.ids == in.ids && real.ids == in.ids && full.zspaceVoids->ids == in.ids, "IDs given");
    BoxConfig c = rsdConfig;
    c.rsd->zspaceVoids = false;
    const BitVFResult r = runBox(in.tracers, {}, distances, c);
    bool rows = true;
    for (std::size_t j = 0; j < r.reconstruction->nObjects; ++j)
      rows = rows && same_bits(std::vector<double>(r.reconstruction->tracers.begin() + 3 * j,
                                                   r.reconstruction->tracers.begin() + 3 * j + 3),
                               std::vector<double>(r.realSpace->cartesian.begin() + 3 * r.reconstruction->ids[j],
                                                   r.realSpace->cartesian.begin() + 3 * r.reconstruction->ids[j] + 3));
    check(r.zspaceReconstruction->ids.empty() && r.realSpace->ids.empty() && r.voids.ids.empty() &&
          r.reconstruction->ids.size() == r.reconstruction->nObjects && rows,
          "no IDs: the second reconstruction's ids are the input rows, every other empty");
    check(!r.zspaceVoids, "no zspace search unless asked");
  }

  group("keepReconstructions off: nothing kept, the same voids");
  {
    BoxConfig c = rsdConfig;
    c.keepReconstructions = false;
    const BitVFResult r = runBox(in.tracers, {}, distances, c, in.ids);
    check(!r.reconstruction && !r.zspaceReconstruction && !r.realSpace, "none kept");
    check(result_hash(r.voids) == result_hash(full.voids) && result_hash(*r.zspaceVoids) == result_hash(*full.zspaceVoids),
          "the searches unchanged");
    BoxConfig plain = in.config;
    const BitVFResult p = runBox(in.inner(), {}, plain);
    check(!p.reconstruction && !p.zspaceReconstruction && !p.realSpace && !p.zspaceVoids, "none kept, no correction");
  }

  group("the errors, in order, before the first reconstruction up to the box");
  {
    std::ostringstream otLog;
    BoxConfig c = rsdConfig;
    c.ot.log = &otLog;
    const std::vector<double>& t = in.tracers;
    auto fails = [&] (const BoxConfig& conf, const std::vector<double>& tr, const std::vector<double>& rd,
                      const std::vector<std::int64_t>& id, const std::string& text, const std::string& what) {
      throws_naming([&] { runBox(tr, rd, distances, conf, id); }, text, what);
    };
    fails(c, {}, {}, {}, "tracers is empty", "PB1 empty");
    fails(c, {1., 2.}, {}, {}, "tracers holds 2 entries, which is not a multiple of three", "PB1 size");
    std::vector<double> bad = t;
    bad[5] = kNaN;
    fails(c, bad, {}, {}, "tracers holds a non-finite value at entry 5", "PB1 finite");
    fails(c, t, {1., 2., 3., 4.}, {}, "the random array holds 4 entries, which is not a multiple of three", "PB2 size");
    fails(c, t, {1., 2., kNaN}, {}, "the random array holds a non-finite value at entry 2", "PB2 finite");
    fails(c, t, {}, {1, 2}, "ids holds 2 entries for 2004 tracers", "PB3 size");
    std::vector<std::int64_t> twice = in.ids;
    twice[3] = twice[1];
    fails(c, t, {}, twice, "twice, at entries 1 and 3", "PB3 twice");
    BoxConfig l = c;
    l.limits->upper[1] = -1.;
    fails(l, t, {}, {}, "the box limits are [0.000000, -1.000000] on axis y", "PB4");
    BoxConfig m = c;
    m.mps = -1.;
    fails(m, t, {}, {}, "mps is -1.000000 in coordinate units; it must be finite and positive", "PB5");
    BoxConfig r0 = c;
    r0.ot.nRealizations = 0;
    fails(r0, t, {}, {}, "the number of realizations is 0; it must be at least 1", "PB6");
    BoxConfig s = c;
    s.voids.smoothing.sigma = -1.;
    fails(s, t, {}, {}, "the smoothing sigma is", "PB7 C1");
    throws_naming([&] { runBox(t, {}, c, in.ids); }, "needs a distance table carrying the growth rate; call runBox with distances",
                  "PB8");
    BoxConfig a = c;
    a.rsd->axis.reset();
    fails(a, t, {}, {}, "the line-of-sight axis (BoxRsd::axis) is not set; it must be 0, 1 or 2", "PB9 empty");
    a.rsd->axis = 3;
    fails(a, t, {}, {}, "the line-of-sight axis is 3; it must be 0, 1 or 2", "PB9 above 2");
    BoxConfig w = c;
    w.rsd->sigma = kNaN;
    fails(w, t, {}, {}, "the smoothing width is nan in coordinate units; it must be finite and non-negative", "PB10");
    BoxConfig b = c;
    b.rsd->bias = 0.;
    fails(b, t, {}, {}, "otswap: ", "PB11 bias 0");
    check(nests_otswap([&] { runBox(t, {}, distances, b); }), "PB11: the otswap::Error nested");
    b.rsd->bias = 1.;
    b.rsd->redshift = kNaN;
    fails(b, t, {}, {}, "otswap: ", "PB11 redshift NaN");
    const otswap::DistanceTable noGrowth({0., 1.}, {0., 3000.});
    throws_naming([&] { runBox(t, {}, noGrowth, c); }, "otswap: ", "PB11 no growth rate");
    std::vector<double> flat = in.inner();
    for (std::size_t i = 0; i < Pipeline::kN; ++i) flat[3*i+2] = 4.;
    BoxConfig f = in.config;
    fails(f, flat, {}, {}, "the tracers span no volume (all of them in one plane, on one line or at one point), so they "
          "fix no box; give the box limits (BoxConfig::limits)", "PB12");
    fails(f, flat, {1., 1., 4., 2., 2., 4.}, {}, "the tracers and the given randoms span no volume", "PB12 randoms");
    check(otLog.str().empty(), "no reconstruction ran before these errors");
  }

  group("the errors of the run");
  {
    BoxConfig c = in.config;
    c.limits = box(200., 300.);
    throws_naming([&] { runBox(in.tracers, {}, c); }, "keeps 0 of 2004 tracers; at least one is needed to draw the randoms",
                  "PB13");
    c = in.config;
    throws_naming([&] { runBox(in.inner(), {1., 2., 3.}, c); }, "otswap: ", "PB15 first reconstruction");
    check(nests_otswap([&] { runBox(in.inner(), {1., 2., 3.}, c); }), "PB15: the otswap::Error nested");
    c.voids.smoothing.sigma = 1e200;
    throws_naming([&] { runBox(in.inner(), {}, c); }, "too large to build the kernel", "PB16 C4");

    std::vector<double> corner, far;
    std::mt19937 rng(5);
    for (std::size_t i = 0; i < 3 * 150; ++i) corner.push_back(5. * (rng() / 4294967296.));
    for (std::size_t i = 0; i < 3 * 300; ++i) far.push_back(95. + 5. * (rng() / 4294967296.));
    BoxConfig e = in.config;
    e.limits = box(0., 100.);
    e.mps = 10.;
    e.rsd = in.rsd();
    e.rsd->bias = 0.01;
    e.rsd->sigma = 0.;
    throws_naming([&] { runBox(corner, far, distances, e); },
                  "no tracer is left for the real-space reconstruction: of 150 tracers, 0 were left uncorrected and 150 "
                  "were moved outside the first pass's grid", "PB17");
  }

  group("the report: the pipeline's lines at config.voids.verbosity, the call line last");
  {
    std::ostringstream log;
    BoxConfig c = rsdConfig;
    c.voids.verbosity = otswap::Verbosity::Detailed;
    c.voids.log = &log;
    const BitVFResult r = runBox(in.tracers, {}, distances, c, in.ids);
    const std::string text = masked(log.str());
    std::ostringstream drew, pass, call;
    drew << "bitvf: drew 4000 randoms in [0, 100] x [0, 100] x [0, 100] with seed " << derived << " in T s\n";
    pass << "bitvf: second pass on " << r.reconstruction->nObjects << " of 2004 tracers: " << r.voids.counts.tracersUncorrected
         << " uncorrected, " << r.voids.counts.tracersOutsideGrid << " outside the first pass's grid in T s\n";
    call << "bitvf: runBox: " << r.voids.voids.nBasins << " voids and " << r.voids.overdense.nBasins
         << " overdense basins from 2004 tracers in T s\n";
    const std::size_t p = text.find(pass.str());
    check(text.rfind(drew.str(), 0) == 0, "the draw line first");
    check(p != std::string::npos && text.find("findVoidsBox") < p && text.find("findVoidsBox", p) != std::string::npos,
          "the second-pass line between the two searches");
    check(text.size() >= call.str().size() && text.compare(text.size() - call.str().size(), std::string::npos, call.str()) == 0,
          "the call line last");

    std::ostringstream normal;
    c.voids.verbosity = otswap::Verbosity::Normal;
    c.voids.log = &normal;
    runBox(in.tracers, {}, distances, c, in.ids);
    check(masked(normal.str()).find("drew") == std::string::npos && masked(normal.str()).find("second pass") == std::string::npos &&
          masked(normal.str()).find(call.str()) != std::string::npos, "Normal: the call line, not the Detailed ones");

    std::ostringstream silent;
    c.voids.verbosity = otswap::Verbosity::Silent;
    c.voids.log = &silent;
    runBox(in.tracers, {}, distances, c, in.ids);
    check(silent.str().empty(), "Silent writes nothing");

    std::ostringstream thrown;
    BoxConfig e = c;
    e.voids.verbosity = otswap::Verbosity::Normal;
    e.voids.log = &thrown;
    e.voids.smoothing.sigma = 1e200;
    throws_naming([&] { runBox(in.tracers, {}, distances, e, in.ids); }, "too large", "a throwing run");
    check(thrown.str().find("runBox") == std::string::npos, "a throwing run writes no call line");
  }

  group("hashes, equal with 1 and 4 threads");
  {
    std::uint64_t hashes[2][2];
    for (const int threads : {1, 4}) {
      set_threads(threads);
      const BitVFResult r = runBox(in.tracers, {}, distances, rsdConfig, in.ids);
      const BitVFResult p = runBox(in.inner(), {}, in.config);
      hashes[threads == 4][0] = result_hash(r.voids) ^ (result_hash(*r.zspaceVoids) * 3);
      hashes[threads == 4][1] = result_hash(p.voids);
    }
    check(hashes[0][0] == hashes[1][0] && hashes[0][1] == hashes[1][1], "1 and 4 threads");
    check(hashes[0][0] == 0x08d7a96a4e990a1full, "recorded hash, with correction");
    check(hashes[0][1] == 0x01ab24d1c86b5ef4ull, "recorded hash, without");
  }

  return report("test_runbox");
}
