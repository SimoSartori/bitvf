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
 *  @file tests/test_checks.cpp
 *
 *  @brief The argument checks of findVoidsBox, both overloads.
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "internal.h"
#include "check.h"

using namespace bitvf;

namespace {

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
  constexpr std::size_t kN = 200;
  constexpr std::size_t kR = 2;
  constexpr double kMps = 10.;

  /// A valid input: tracers uniform in [0, 100)^3, displacements within 2
  /// of zero, all valid, no randoms.
  struct Input {
    std::vector<double> tracers, displacement, randoms;
    std::vector<std::uint8_t> valid;
    double mps = kMps;
    Config config;
    std::optional<otswap::BoxLimits> limits;
    std::vector<std::int64_t> ids;
    const DivergenceGrid* grid = nullptr;

    Input ()
    {
      std::mt19937 rng(2026);
      auto uniform = [&rng] (const double lo, const double hi) {
        return lo + (hi - lo) * (rng() / 4294967296.);
      };
      for (std::size_t i = 0; i < 3 * kN; ++i) tracers.push_back(uniform(0., 100.));
      for (std::size_t i = 0; i < 3 * kR * kN; ++i) displacement.push_back(uniform(-2., 2.));
      valid.assign(kR * kN, 1);
      config.log = nullptr;
    }

    Result run () const
    {
      return findVoidsBox(tracers, displacement, valid, randoms, mps, config, limits, ids, grid);
    }

    internal::BoxSearch search () const
    {
      return internal::check_box(tracers, displacement, valid, randoms, mps, config, limits, ids,
                                 grid);
    }

    otswap::Result reconstruction () const
    {
      otswap::Result rec;
      rec.nObjects = kN;
      rec.nRealizations = kR;
      rec.geometry = otswap::Geometry::Box;
      rec.tracers = tracers;
      rec.displacement = displacement;
      rec.valid = valid;
      rec.mps = mps;
      return rec;
    }
  };

  /// A call that passes every check ends in the error of the unbuilt search.
  template <typename F>
  void accepted (F&& f, const std::string& what)
  {
    throws_naming(f, "not implemented", what);
  }

  otswap::BoxLimits box (const double lo, const double hi)
  {
    otswap::BoxLimits limits;
    limits.lower = {lo, lo, lo};
    limits.upper = {hi, hi, hi};
    return limits;
  }

}

int main ()
{
  group("a valid input passes every check");
  {
    accepted([] { Input().run(); }, "the base input");

    Input in;
    in.valid.clear();
    accepted([&] { in.run(); }, "empty valid means all valid");

    Input logged;
    logged.config.log = nullptr;
    logged.config.verbosity = otswap::Verbosity::Detailed;
    accepted([&] { logged.run(); }, "a null log");

    Input withRandoms;
    for (std::size_t i = 0; i < 3 * kN; ++i) withRandoms.randoms.push_back(50.);
    accepted([&] { withRandoms.run(); }, "randoms given");
  }

  group("the tracers: empty, not a multiple of three, a non-finite coordinate on any row");
  {
    Input in;
    in.tracers.clear();
    throws_naming([&] { in.run(); }, "tracers is empty", "T1");

    Input bad;
    bad.tracers.push_back(1.);
    throws_naming([&] { bad.run(); }, "not a multiple of three", "T2");

    Input nan;
    nan.tracers[5] = kNaN;
    throws_naming([&] { nan.run(); }, "tracers holds a non-finite value at entry 5", "T3");

    Input outside;
    outside.limits = box(0., 50.);
    outside.tracers[0] = 200.;
    outside.tracers[1] = kNaN;
    throws_naming([&] { outside.run(); }, "tracers holds a non-finite",
                  "T3 reads the rows outside the limits");
  }

  group("the displacements: their count, and finite where read");
  {
    Input empty;
    empty.displacement.clear();
    throws_naming([&] { empty.run(); }, "displacement holds 0 entries", "D1, empty");

    Input partial;
    partial.displacement.push_back(0.);
    throws_naming([&] { partial.run(); }, "whole nRealizations", "D1, not whole");

    Input nan;
    nan.displacement[3 * kN + 4] = kNaN;
    throws_naming([&] { nan.run(); }, "realization 1, tracer 1, component 1", "D3");

    Input invalid;
    invalid.displacement[3 * kN + 4] = kNaN;
    invalid.valid[kN + 1] = 0;
    accepted([&] { invalid.run(); }, "D3 skips an invalid displacement");

    Input skipped;
    skipped.limits = box(-1., 101.);
    skipped.tracers[3] = 200.;
    skipped.displacement[3] = kNaN;
    accepted([&] { skipped.run(); }, "D3 skips a tracer outside the limits");
  }

  group("valid: its size and its values");
  {
    Input in;
    in.valid.pop_back();
    throws_naming([&] { in.run(); }, "valid holds 399 entries", "V1");

    Input two;
    two.valid[7] = 2;
    throws_naming([&] { two.run(); }, "valid holds 2 at entry 7", "V2");
  }

  group("the randoms: a multiple of three, finite");
  {
    Input in;
    in.randoms = {1., 2., 3., 4.};
    throws_naming([&] { in.run(); }, "randoms holds 4 entries", "R1");

    Input nan;
    nan.randoms = {1., kNaN, 3.};
    throws_naming([&] { nan.run(); }, "randoms holds a non-finite", "R2");
  }

  group("mps: finite and positive");
  for (const double mps : {0., -1., kNaN, std::numeric_limits<double>::infinity()}) {
    Input in;
    in.mps = mps;
    throws_naming([&] { in.run(); }, "mps is", "M1");
  }

  group("the limits: finite, lower below upper; a point on a face is inside");
  {
    Input nan;
    nan.limits = box(0., 100.);
    nan.limits->upper[1] = kNaN;
    throws_naming([&] { nan.run(); }, "every bound must be finite", "L1, not finite");

    Input flat;
    flat.limits = box(0., 100.);
    flat.limits->upper[2] = 0.;
    throws_naming([&] { flat.run(); }, "the lower bound must be below the upper one", "L1, order");

    Input face;
    face.limits = box(0., 100.);
    face.tracers[0] = 100.;
    face.tracers[1] = 0.;
    face.randoms = {100., 0., 50., 100.5, 50., 50.};
    const internal::BoxSearch search = face.search();
    check(search.inBox[0] == 1 && search.counts.tracersOutsideBox == 0,
          "a tracer on a face is inside");
    check(search.counts.randomsOutsideBox == 1 && search.randoms.size() == 3,
          "a random on a face is kept, one outside removed");
  }

  group("ids: one per tracer, no value twice");
  {
    Input in;
    in.ids.assign(kN - 1, 0);
    throws_naming([&] { in.run(); }, "ids holds 199 entries", "I1");

    Input twice;
    for (std::size_t i = 0; i < kN; ++i) twice.ids.push_back((std::int64_t)i - 100);
    twice.ids[150] = twice.ids[3];
    throws_naming([&] { twice.run(); }, "twice, at entries 3 and 150", "I2");
  }

  group("the configuration");
  {
    for (const double sigma : {-1., kNaN}) {
      Input in;
      in.config.smoothing.sigma = sigma;
      throws_naming([&] { in.run(); }, "smoothing sigma", "C1");
    }
    Input zero;
    zero.config.smoothing.sigma = 0.;
    accepted([&] { zero.run(); }, "C1 accepts sigma 0");

    for (const double cell : {0., -1., std::numeric_limits<double>::infinity()}) {
      Input in;
      in.config.divergenceCellSize = cell;
      throws_naming([&] { in.run(); }, "divergence cell size is", "C2");
    }

    Input unsmoothed;
    unsmoothed.config.keepUnsmoothedDivergence = true;
    throws_naming([&] { unsmoothed.run(); }, "requires keeping the divergence", "C3");
    unsmoothed.config.keepDivergence = true;
    accepted([&] { unsmoothed.run(); }, "C3 accepts both");
  }

  group("a reused grid: its defining fields, and its combinations with the arguments");
  {
    const DivergenceGrid good = internal::box_grid(Input().search(), 1.);

    auto with = [&good] (auto&& change, const std::string& text, const std::string& what) {
      DivergenceGrid grid = good;
      change(grid);
      Input in;
      in.grid = &grid;
      throws_naming([&] { in.run(); }, text, what);
    };

    with([] (DivergenceGrid& g) { g.geometry = otswap::Geometry::Lightcone; },
         "is not a box grid", "G1");
    with([] (DivergenceGrid& g) { g.cellSize = 0.; }, "has a cell size of", "G2");
    with([] (DivergenceGrid& g) { g.lower[0] = kNaN; }, "every bound must be finite", "G3");
    with([] (DivergenceGrid& g) { g.upper[2] = g.lower[2]; }, "lower bound must be below", "G3");
    with([] (DivergenceGrid& g) { g.nCellsAxis[1] += 1; }, "three counts must be equal", "G4");
    with([] (DivergenceGrid& g) { g.nCellsAxis = {0, 0, 0}; g.nCells = 0; },
         "three counts must be equal and positive", "G4, zero");
    with([] (DivergenceGrid& g) { g.nCells += 1; }, "it must have their cube", "G5");
    with([] (DivergenceGrid& g) { g.mps = kNaN; }, "has an mps of", "G6");
    with([] (DivergenceGrid& g) { g.raMin = 0.; }, "field of a lightcone grid", "G7");
    with([] (DivergenceGrid& g) { g.mpsProfile.redshift = {0.5}; }, "field of a lightcone grid",
         "G7, the profile");
    with([] (DivergenceGrid& g) { g.mps = std::nextafter(g.mps, 20.); }, "they must be equal",
         "X2, one ulp");

    DivergenceGrid reused = good;
    reused.divergence = {1., 2.};
    reused.voidIndex = {7};
    Input in;
    in.grid = &reused;
    accepted([&] { in.run(); }, "per-cell fields are ignored");

    in.config.divergenceCellSize = 1.;
    throws_naming([&] { in.run(); }, "set together with a grid", "X1");

    Input limited;
    limited.grid = &good;
    limited.limits = box(0., 100.);
    accepted([&] { limited.run(); }, "X3: limits with a grid");

    Input small;
    DivergenceGrid half = good;
    half.nCellsAxis = {good.nCellsAxis[0] / 2, good.nCellsAxis[0] / 2, good.nCellsAxis[0] / 2};
    half.nCells = half.nCellsAxis[0] * half.nCellsAxis[0] * half.nCellsAxis[0];
    small.grid = &half;
    small.randoms = {1., 1., 1., 99., 99., 99.};
    const internal::BoxSearch search = small.search();
    check(search.counts.tracersOutsideGrid > 0 && search.counts.tracersOutsideGrid < kN,
          "tracers outside a reused grid are flagged");
    check(search.counts.randomsOutsideGrid == 1 && search.randoms.size() == 3,
          "randoms outside a reused grid are removed");
  }

  group("the data: something left to read, the size of the grid");
  {
    Input none;
    none.valid.assign(kR * kN, 0);
    throws_naming([&] { none.run(); }, "no valid displacement is left", "E1, all invalid");

    Input away;
    away.limits = box(200., 300.);
    throws_naming([&] { away.run(); }, "no valid displacement is left", "E1, all outside");

    Input randoms;
    randoms.limits = box(-10., 110.);
    randoms.randoms = {200., 200., 200.};
    throws_naming([&] { randoms.run(); }, "no random is left", "E2");

    Input huge;
    huge.config.divergenceCellSize = 1e-9;
    throws_naming([&] { huge.run(); }, "an int64 cell index can address", "E3");
  }

  group("the otswap::Result overload: its own checks, and the same selection");
  {
    const Input in;
    const Config config = in.config;

    auto with = [&] (auto&& change, const std::string& text, const std::string& what) {
      otswap::Result rec = in.reconstruction();
      change(rec);
      throws_naming([&] { findVoidsBox(rec, config); }, text, what);
    };

    accepted([&] { findVoidsBox(in.reconstruction(), config); }, "a valid reconstruction");

    with([] (otswap::Result& r) { r.geometry = otswap::Geometry::Lightcone; },
         "not a box result", "O1");
    with([] (otswap::Result& r) { r.nObjects = 0; }, "holds no objects", "O2, no objects");
    with([] (otswap::Result& r) { r.nRealizations = 0; }, "holds no realizations",
         "O2, no realizations");
    with([] (otswap::Result& r) { r.tracers.pop_back(); }, "malformed: tracers", "O2, tracers");
    with([] (otswap::Result& r) { r.displacement.pop_back(); }, "malformed: displacement",
         "O2, displacement");
    with([] (otswap::Result& r) { r.valid.clear(); }, "malformed: valid", "O2, valid");
    with([] (otswap::Result& r) { r.outsideBox.assign(3, 0); }, "malformed: outsideBox",
         "O2, outsideBox");
    with([] (otswap::Result& r) { r.ids.assign(3, 0); }, "malformed: ids", "O2, ids");
    with([] (otswap::Result& r) { r.valid[0] = 3; }, "the reconstruction's valid holds 3", "O3");
    with([] (otswap::Result& r) { r.ids.assign(kN, 5); }, "twice", "O3, ids");
    with([] (otswap::Result& r) { r.tracers[4] = kNaN; }, "the reconstruction's tracers",
         "O4, tracers");
    with([] (otswap::Result& r) { r.displacement[4] = kNaN; },
         "the reconstruction's displacement", "O4, displacement");
    with([] (otswap::Result& r) { r.mps = kNaN; }, "the reconstruction's mps", "O5");
    with([] (otswap::Result& r) { r.selection.boxLimits = box(5., 1.); },
         "the reconstruction's box limits", "O6");

    otswap::Result flagged = in.reconstruction();
    flagged.selection.boxLimits = box(-10., 110.);
    flagged.outsideBox.assign(kN, 0);
    flagged.outsideBox[1] = 1;
    flagged.tracers[4] = kNaN;
    flagged.displacement[4] = kNaN;
    accepted([&] { findVoidsBox(flagged, config); }, "O4 skips a flagged row unread");

    const DivergenceGrid grid = internal::box_grid(in.search(), 1.);
    otswap::Result other = in.reconstruction();
    other.mps = kNaN;
    accepted([&] { findVoidsBox(other, config, {}, &grid); },
             "O5 and X2: with a grid, the grid's mps wins");

    const otswap::Result rec = in.reconstruction();
    const internal::BoxSearch a = in.search();
    const internal::BoxSearch b = internal::check_box(rec, config, {}, nullptr);
    check(a.nTracers == b.nTracers && a.nRealizations == b.nRealizations && a.mps == b.mps &&
          a.inBox == b.inBox && a.randoms == b.randoms &&
          a.counts.displacements == b.counts.displacements &&
          a.counts.randoms == b.counts.randoms,
          "both overloads select the same from the same arrays");
    check(internal::box_grid(a, 1.).lower == internal::box_grid(b, 1.).lower,
          "both overloads build the same grid");
  }

  return report("test_checks");
}
