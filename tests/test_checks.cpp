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
 *  @brief The argument checks of findVoidsBox, both overloads, and the
 *  Result they return.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "otswap/Maths.h"

#include "internal.h"
#include "check.h"

using namespace bitvf;

namespace {

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
  constexpr double kInf = std::numeric_limits<double>::infinity();
  constexpr std::size_t kN = 200;
  constexpr std::size_t kR = 2;
  constexpr double kMps = 10.;

  otswap::BoxLimits box (const double lo, const double hi)
  {
    otswap::BoxLimits limits;
    limits.lower = {lo, lo, lo};
    limits.upper = {hi, hi, hi};
    return limits;
  }

  /// A valid input: tracers uniform in [0, 100)^3, displacements within 2
  /// of zero, all valid. Its reconstruction is recorded as otswap records
  /// one given the limits [-2, 102]^3, which hold every endpoint.
  struct Input {
    std::vector<double> tracers, displacement;
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
      return findVoidsBox(tracers, displacement, valid, mps, config, limits, ids, grid);
    }

    internal::BoxSearch search () const
    {
      return internal::check_box(tracers, displacement, valid, mps, config, limits, ids, grid);
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
      rec.selection.boxLimits = box(-2., 102.);
      rec.box = rec.selection.boxLimits;
      return rec;
    }
  };

  void threads (const int n)
  {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
  }

  /// FNV-1a over raw bytes.
  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  template <typename T>
  void hash_vector (std::uint64_t& h, const std::vector<T>& v)
  {
    hash_bytes(h, v.data(), v.size() * sizeof(T));
  }

  /// The hash of every column of both catalogues.
  std::uint64_t catalogue_hash (const Result& r)
  {
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const BasinCatalog* c : {&r.voids, &r.overdense}) {
      hash_bytes(h, &c->nBasins, sizeof(c->nBasins));
      hash_vector(h, c->centre);
      hash_vector(h, c->centreSky);
      hash_vector(h, c->radius);
      hash_vector(h, c->radiusExt);
      hash_vector(h, c->deltaAv);
      hash_vector(h, c->delta2Mps);
      hash_vector(h, c->edgeFlag);
      hash_vector(h, c->contamination);
      hash_vector(h, c->nTracers);
    }
    return h;
  }

  /// The hash of everything a Result holds except elapsedSeconds.
  std::uint64_t result_hash (const Result& r)
  {
    std::uint64_t h = catalogue_hash(r);
    hash_bytes(h, &r.nTracers, sizeof(r.nTracers));
    hash_bytes(h, &r.nRealizations, sizeof(r.nRealizations));
    hash_bytes(h, &r.mps, sizeof(r.mps));
    const double cellSize = r.config.divergenceCellSize.value_or(-1.);
    hash_bytes(h, &cellSize, sizeof(cellSize));
    hash_vector(h, r.ids);
    hash_vector(h, r.outsideGrid);
    const std::size_t counts[] = {r.counts.tracers, r.counts.tracersUsed, r.counts.tracersOutsideBox,
                                  r.counts.tracersOutsideGrid, r.counts.displacements,
                                  r.counts.displacementsRejected};
    hash_bytes(h, counts, sizeof(counts));
    if (r.divergence) {
      const DivergenceGrid& g = *r.divergence;
      hash_bytes(h, g.lower.data(), sizeof(g.lower));
      hash_bytes(h, &g.cellSize, sizeof(g.cellSize));
      hash_bytes(h, &g.nCells, sizeof(g.nCells));
      hash_vector(h, g.centre);
      hash_vector(h, g.volume);
      hash_vector(h, g.divergence);
      hash_vector(h, g.weight);
      hash_vector(h, g.nDisplacements);
      hash_vector(h, g.nTracers);
      hash_vector(h, g.voidIndex);
      hash_vector(h, g.overdenseIndex);
    }
    if (r.assignment) {
      const TracerAssignment& a = *r.assignment;
      hash_vector(h, a.cellIndex);
      hash_vector(h, a.voidIndex);
      hash_vector(h, a.overdenseIndex);
      hash_vector(h, a.distance);
      hash_vector(h, a.inMask);
      hash_vector(h, a.inBox);
      hash_vector(h, a.used);
    }
    return h;
  }

  /// Every column of a catalogue of the right size.
  bool sized (const BasinCatalog& c)
  {
    const std::size_t B = c.nBasins;
    return c.centre.size() == 3 * B && c.centreSky.empty() && c.radius.size() == B &&
           c.radiusExt.size() == B && c.deltaAv.size() == B && c.delta2Mps.size() == B &&
           c.edgeFlag.size() == B && c.contamination.size() == B && c.nTracers.size() == B;
  }

  /// A call that passes every check returns a Result.
  template <typename F>
  void accepted (F&& f, const std::string& what)
  {
    try {
      f();
    }
    catch (const std::exception& e) {
      check(false, what + " (threw \"" + e.what() + "\")");
    }
  }

  /// The call passes the smoothing checks: it returns, or fails only after
  /// the kernel is built, finding no void in the smoothed field.
  template <typename F>
  void kernel_built (F&& f, const std::string& what)
  {
    try {
      f();
    }
    catch (const std::exception& e) {
      check(std::string(e.what()).find("no void found") == 0, what + " (threw \"" + e.what() + "\")");
    }
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
  }

  group("the tracers: empty, not a multiple of three, a non-finite coordinate on any row");
  {
    Input in;
    in.tracers.clear();
    throws_naming([&] { in.run(); }, "findVoidsBox: tracers is empty", "T1");

    Input bad;
    bad.tracers.push_back(1.);
    throws_naming([&] { bad.run(); }, "findVoidsBox: tracers holds", "T2");

    Input nan;
    nan.tracers[5] = kNaN;
    throws_naming([&] { nan.run(); }, "the tracers hold a non-finite value at entry 5", "T3");

    Input outside;
    outside.limits = box(0., 50.);
    outside.tracers[0] = 200.;
    outside.tracers[1] = kNaN;
    throws_naming([&] { outside.run(); }, "the tracers hold a non-finite",
                  "T3 reads the rows outside the limits");
  }

  group("the displacements: their count, and finite where read");
  {
    Input empty;
    empty.displacement.clear();
    throws_naming([&] { empty.run(); }, "findVoidsBox: displacement holds 0 entries", "D1, empty");

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
    throws_naming([&] { in.run(); }, "findVoidsBox: valid holds 399 entries", "V1");

    Input two;
    two.valid[7] = 2;
    throws_naming([&] { two.run(); }, "the validity flags hold 2 at entry 7", "V2");
  }

  group("mps: finite and positive");
  for (const double mps : {0., -1., kNaN, std::numeric_limits<double>::infinity()}) {
    Input in;
    in.mps = mps;
    throws_naming([&] { in.run(); }, "the mean particle separation is", "M1");
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
    const internal::BoxSearch search = face.search();
    check(search.inBox[0] == 1 && search.counts.tracersOutsideBox == 0,
          "a tracer on a face is inside");
  }

  group("ids: one per tracer, no value twice");
  {
    Input in;
    in.ids.assign(kN - 1, 0);
    throws_naming([&] { in.run(); }, "findVoidsBox: ids holds 199 entries", "I1");

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

    Input mps;
    mps.config.smoothing.sigma = 0.7;
    check(mps.search().sigmaLength == 0.7 * kMps, "sigmaLength is sigma times mps in Mps");
    Input coordinate;
    coordinate.config.smoothing.sigma = 0.7;
    coordinate.config.smoothing.unit = LengthUnit::Coordinate;
    check(coordinate.search().sigmaLength == 0.7, "sigmaLength is sigma in Coordinate");

    Input averaged;
    averaged.config.nRealizationsMean = 0;
    throws_naming([&] { averaged.run(); }, "the number of realizations averaged is 0", "C5");
    throws_naming([&] { findVoidsBox(averaged.reconstruction(), averaged.config); },
                  "the number of realizations averaged is 0", "C5, reconstruction");
    averaged.config.nRealizationsMean = 7;
    check(averaged.run().config.nRealizationsMean == 7, "nRealizationsMean copied into the Result");
  }

  group("the smoothing sigma in coordinate units, (2 sigma) sigma positive and finite: both "
        "overloads, both units, both sides of both limits (C4)");
  {
    const double tiny = std::ldexp(1., -538);
    double big = std::sqrt(std::numeric_limits<double>::max() / 2);
    while (std::isfinite(2 * std::nextafter(big, kInf) * std::nextafter(big, kInf)))
      big = std::nextafter(big, kInf);
    while (!std::isfinite(2 * big * big)) big = std::nextafter(big, 0.);
    check(1e-162 < tiny && tiny < 1.2e-162 && 9e153 < big && big < 1e154,
          "the limits are about 1e-162 and 1e154");

    struct Sigma { double sigma; LengthUnit unit; double mps; const char* fails; };
    const Sigma cases[] = {
      {tiny, LengthUnit::Coordinate, kMps, "too small"},
      {std::nextafter(tiny, 1.), LengthUnit::Coordinate, kMps, nullptr},
      {1e-200, LengthUnit::Coordinate, kMps, "too small"},
      {1e-150, LengthUnit::Coordinate, kMps, nullptr},
      {big, LengthUnit::Coordinate, kMps, nullptr},
      {std::nextafter(big, kInf), LengthUnit::Coordinate, kMps, "too large"},
      {std::numeric_limits<double>::max(), LengthUnit::Coordinate, kMps, "too large"},
      {1e-163, LengthUnit::Mps, 1., "too small"},
      {tiny, LengthUnit::Mps, 1., "too small"},
      {std::nextafter(tiny, 1.), LengthUnit::Mps, 1., nullptr},
      {1e-160, LengthUnit::Mps, 1e-10, "too small"},
      {big, LengthUnit::Mps, 1., nullptr},
      {std::nextafter(big, kInf), LengthUnit::Mps, 1., "too large"},
      {9e152, LengthUnit::Mps, kMps, nullptr},
      {1e153, LengthUnit::Mps, kMps, "too large"},
      {1e300, LengthUnit::Mps, kMps, "too large"},
    };
    for (const Sigma& c : cases) {
      Input in;
      in.mps = c.mps;
      in.config.smoothing.sigma = c.sigma;
      in.config.smoothing.unit = c.unit;
      const otswap::Result rec = in.reconstruction();
      std::ostringstream label;
      label << std::setprecision(17) << "C4, sigma " << c.sigma << " ("
            << (c.unit == LengthUnit::Mps ? "Mps" : "Coordinate") << "), mps " << c.mps;
      const std::string what = label.str();
      if (c.fails == nullptr) {
        kernel_built([&] { in.run(); }, what + ", arrays: accepted");
        kernel_built([&] { findVoidsBox(rec, in.config); }, what + ", reconstruction: accepted");
        continue;
      }
      const std::string text = std::string(c.fails) + " to build the kernel";
      throws_naming([&] { in.run(); }, "smoothing sigma is", what + ", arrays");
      throws_naming([&] { in.run(); }, text, what + ", arrays");
      throws_naming([&] { findVoidsBox(rec, in.config); }, "smoothing sigma is",
                    what + ", reconstruction");
      throws_naming([&] { findVoidsBox(rec, in.config); }, text, what + ", reconstruction");
    }
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
    for (const double c : {0., -1., kNaN, kInf})
      with([c] (DivergenceGrid& g) { g.cellSizeMps = c; }, "mps (cellSizeMps); it must be finite and positive",
           "G8");
    with([] (DivergenceGrid& g) { g.cellSizeMps = std::nextafter(g.cellSizeMps, 2.); },
         "but cellSizeMps * mps is", "G9, one ulp of cellSizeMps");
    with([] (DivergenceGrid& g) { g.cellSize = std::nextafter(g.cellSize, 0.); },
         "but cellSizeMps * mps is", "G9, one ulp of cellSize");
    with([] (DivergenceGrid& g) { g.raMin = 0.; }, "field of a lightcone grid", "G7");
    with([] (DivergenceGrid& g) { g.mpsProfile.redshift = {0.5}; }, "field of a lightcone grid",
         "G7, the profile");
    with([] (DivergenceGrid& g) { g.mps = std::nextafter(g.mps, 20.); }, "they must be equal",
         "X2, one ulp");

    DivergenceGrid reused = good;
    reused.divergence = {1., 2.};
    reused.voidIndex = {7};
    check(good.cellSizeMps == 1. && good.cellSize == kMps, "box_grid records cellSizeMps beside cellSize");

    Input in;
    in.grid = &good;
    accepted([&] { in.run(); }, "a box grid with cellSizeMps set is no G7");

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
    const internal::BoxSearch search = small.search();
    check(search.counts.tracersOutsideGrid > 0 && search.counts.tracersOutsideGrid < kN,
          "tracers outside a reused grid are flagged");
  }

  group("the data: something left to read, the size of the grid");
  {
    Input none;
    none.valid.assign(kR * kN, 0);
    throws_naming([&] { none.run(); }, "no valid displacement is left", "E1, all invalid");

    Input away;
    away.limits = box(200., 300.);
    throws_naming([&] { away.run(); }, "no valid displacement is left", "E1, all outside");

    Input huge;
    huge.config.divergenceCellSize = 1e-9;
    throws_naming([&] { huge.run(); }, "an int64 cell index can address", "E3");
  }

  group("the divergence: an endpoint outside every cell, a fill that stalls");
  {
    DivergenceGrid cells;
    cells.geometry = otswap::Geometry::Box;
    cells.mps = 1.;
    cells.cellSize = 1.;
    cells.cellSizeMps = 1.;
    cells.lower = {0., 0., 0.};
    cells.upper = {4., 4., 4.};
    cells.nCellsAxis = {4, 4, 4};
    cells.nCells = 64;
    Config config;
    config.log = nullptr;
    const std::vector<double> tracers = {0.5, 0.5, 0.5, 1.5, 0.5, 0.5};

    auto run = [&] (const std::vector<double>& displacement) {
      findVoidsBox(tracers, displacement, {}, 1., config, std::nullopt, {}, &cells);
    };
    throws_naming([&] { run({1., 1., 0., 2.5, 0., 0.}); },
                  "displacement of tracer 1 in realization 0 ends outside every cell", "E6");
    throws_naming([&] { run({1., 0., 0., 0.1, 0., 0.}); },
                  "cannot reach 62 empty cells: the last pass filled none", "E4");
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
    with([] (otswap::Result& r) { r.valid[0] = 3; }, "findVoidsBox: the reconstruction's validity flags hold 3", "O3");
    with([] (otswap::Result& r) { r.ids.assign(kN, 5); }, "twice", "O3, ids");
    with([] (otswap::Result& r) { r.tracers[4] = kNaN; }, "the reconstruction's tracers",
         "O4, tracers");
    with([] (otswap::Result& r) { r.displacement[4] = kNaN; },
         "findVoidsBox: the reconstruction's displacements hold a non-finite", "O4, displacement");
    with([] (otswap::Result& r) { r.mps = kNaN; }, "the reconstruction's mps", "O5");
    with([] (otswap::Result& r) { r.selection.boxLimits = box(5., 1.); },
         "the reconstruction's box limits", "O6");
    with([] (otswap::Result& r) { r.box.reset(); }, "it holds no box", "O6, no box");
    with([] (otswap::Result& r) { r.box->upper[1] = kNaN; }, "every bound must be finite",
         "O6, box not finite");
    with([] (otswap::Result& r) { r.box = box(5., 1.); }, "must not be above the upper one",
         "O6, box reversed");
    with([] (otswap::Result& r) { r.selection.boxLimits = box(-10., 110.); },
         "differs from its box limits", "O6, box not the limits");
    with([] (otswap::Result& r) { r.tracers[3] = r.box->upper[0] + 1.; },
         "tracer 1 lies outside its box", "O7");

    otswap::Result plane = in.reconstruction();
    for (std::size_t i = 0; i < kN; ++i) plane.tracers[3*i+2] = 50.;
    plane.selection.boxLimits.reset();
    plane.box->lower[2] = plane.box->upper[2] = 50.;
    check(internal::check_box(plane, config, nullptr).box->upper[2] == 50.,
          "O6 accepts a flat box");

    otswap::Result flagged = in.reconstruction();
    flagged.selection.boxLimits = box(-10., 110.);
    flagged.box = flagged.selection.boxLimits;
    flagged.outsideBox.assign(kN, 0);
    flagged.outsideBox[1] = 1;
    flagged.tracers[4] = kNaN;
    flagged.displacement[4] = kNaN;
    accepted([&] { findVoidsBox(flagged, config); }, "O4 skips a flagged row unread");

    const DivergenceGrid grid = internal::box_grid(in.search(), 1.);
    otswap::Result other = in.reconstruction();
    other.mps = kNaN;
    accepted([&] { findVoidsBox(other, config, &grid); },
             "O5 and X2: with a grid, the grid's mps wins");

    Input seven;
    seven.mps = 7.;
    const DivergenceGrid grid7 = internal::box_grid(seven.search(), 1.);
    Config sigma = config;
    sigma.smoothing.sigma = 0.7;
    check(internal::check_box(in.reconstruction(), sigma, &grid7).sigmaLength == 0.7 * 7.,
          "sigmaLength in Mps uses the reused grid's mps, not the reconstruction's");

    const otswap::Result rec = in.reconstruction();
    Input limited = in;
    limited.limits = rec.box;
    const internal::BoxSearch a = limited.search();
    const internal::BoxSearch b = internal::check_box(rec, config, nullptr);
    check(a.nTracers == b.nTracers && a.nRealizations == b.nRealizations && a.mps == b.mps &&
          a.sigmaLength == b.sigmaLength &&
          a.maxDisplacement == b.maxDisplacement && a.inBox == b.inBox &&
          a.box->lower == b.box->lower && a.box->upper == b.box->upper &&
          a.counts.displacements == b.counts.displacements,
          "both overloads select the same from the same arrays, the box as the limits");
    check(internal::box_grid(a, 1.).lower == internal::box_grid(b, 1.).lower &&
          internal::box_grid(a, 1.).upper == internal::box_grid(b, 1.).upper,
          "both overloads build the same grid");

    otswap::Result wide = in.reconstruction();
    wide.box = wide.selection.boxLimits = box(-30., 130.);
    Input wideLimits = in;
    wideLimits.limits = wide.box;
    const DivergenceGrid g = internal::box_grid(internal::check_box(wide, config, nullptr), 1.);
    check(g.lower == internal::box_grid(wideLimits.search(), 1.).lower && g.lower[0] < -30. &&
          g.upper[0] > 130.,
          "the grid is built on the reconstruction's box, wider than its tracers");
  }

  group("an endpoint beyond the box: E6 from a reconstruction, none from arrays without limits");
  {
    Input far;
    far.displacement[0] = 150. - far.tracers[0];
    otswap::Result rec = far.reconstruction();
    throws_naming([&] { findVoidsBox(rec, far.config); },
                  "ends outside every cell of the grid, which spans the box limits", "E6, box");
    accepted([&] { far.run(); }, "arrays without limits cover the endpoints");
    far.limits = rec.box;
    throws_naming([&] { far.run(); }, "pass the box limits the reconstruction used", "E6, limits");
  }

  group("an otswap reconstruction on given randoms beyond the tracers, without limits: no E6");
  {
    const Input in;
    std::mt19937 rng(7);
    std::vector<double> randoms;
    for (std::size_t i = 0; i < 3 * kR * kN; ++i) randoms.push_back(-30. + 160. * (rng() / 4294967296.));
    otswap::BoxLimits tracersBox = box(kInf, -kInf);
    for (std::size_t i = 0; i < kN; ++i)
      for (std::size_t a = 0; a < 3; ++a) {
        tracersBox.lower[a] = std::min(tracersBox.lower[a], in.tracers[3*i+a]);
        tracersBox.upper[a] = std::max(tracersBox.upper[a], in.tracers[3*i+a]);
      }
    otswap::Config ot;
    ot.nRealizations = kR;
    ot.seed = 1;
    ot.log = nullptr;

    for (const bool givenMps : {true, false}) {
      const otswap::Result rec = givenMps ? otswap::reconstructBox(in.tracers, randoms, kMps, ot)
                                          : otswap::reconstructBox(in.tracers, randoms, ot);
      const std::string what = givenMps ? "mps given" : "mps computed";
      bool inside = true, beyond = false;
      for (std::size_t j = 0; j < kR * kN; ++j)
        for (std::size_t a = 0; a < 3; ++a) {
          const double end = rec.tracers[3*(j%kN)+a] + rec.displacement[3*j+a];
          inside = inside && end >= rec.box->lower[a] && end <= rec.box->upper[a];
          beyond = beyond || end < tracersBox.lower[a] || end > tracersBox.upper[a];
        }
      check(!rec.selection.boxLimits && inside && beyond,
            what + ": endpoints beyond the tracers' bounding box, all inside otswap::Result::box");

      Input arrays = in;
      arrays.displacement = rec.displacement;
      arrays.valid = rec.valid;
      arrays.mps = rec.mps;
      arrays.limits = tracersBox;
      throws_naming([&] { arrays.run(); }, "ends outside every cell of the grid",
                    what + ": E6 on the tracers' bounding box");
      accepted([&] { findVoidsBox(rec, in.config); }, what + ": the reconstruction accepted");
      arrays.limits = rec.box;
      check(result_hash(arrays.run()) == result_hash(findVoidsBox(rec, in.config)),
            what + ": equal to the arrays with otswap::Result::box as the limits");
    }
  }

  group("the Result: its fields, with and without the kept outputs");
  {
    const Result r = Input().run();
    check(r.geometry == otswap::Geometry::Box && r.nTracers == kN && r.nRealizations == kR && r.mps == kMps,
          "geometry, nTracers, nRealizations, mps");
    check(r.config.divergenceCellSize == 1., "the default cell size, recorded");
    check(r.counts.tracers == kN && r.counts.tracersUsed == kN && r.counts.displacements == kR * kN &&
          r.counts.tracersOutsideBox == 0 && r.counts.tracersOutsideGrid == 0, "counts");
    check(r.ids.empty() && !r.randomSeed && r.outsideGrid.empty() && r.uncorrected.empty() &&
          r.outsideCut.empty() && !r.distances, "the fields a box search without ids leaves empty");
    check(!r.divergence && !r.assignment, "no divergence or assignment unless kept");
    check(r.voids.nBasins > 0 && sized(r.voids) && sized(r.overdense), "both catalogues, every column");
    check(std::isfinite(r.elapsedSeconds) && r.elapsedSeconds >= 0., "elapsedSeconds");

    Input kept;
    kept.config.keepDivergence = true;
    kept.config.keepAssignment = true;
    const Result k = kept.run();
    check(catalogue_hash(k) == catalogue_hash(r), "the kept outputs change no catalogue");
    check(k.divergence && k.assignment, "divergence and assignment kept");
    const DivergenceGrid& g = *k.divergence;
    check(g.nRealizations == kR && g.nRealizationsMean == kept.config.nRealizationsMean,
          "the grid's nRealizations and nRealizationsMean, those of the search");
    DivergenceGrid centres = g;
    internal::box_centres(centres);
    bool volume = g.volume.size() == g.nCells;
    for (const double v : g.volume) volume = volume && v == g.cellSize * g.cellSize * g.cellSize;
    check(g.centre.size() == 3 * g.nCells && g.centre == centres.centre && volume,
          "the grid's centre and volume");
    check(g.divergence.size() == g.nCells && g.divergenceSmoothed.size() == g.nCells &&
          g.weight.size() == g.nCells && g.weightSmoothed.size() == g.nCells && g.voidIndex.size() == g.nCells &&
          g.overdenseIndex.size() == g.nCells, "the grid's per-cell fields");
    bool seeds = true;
    for (std::size_t b = 0; b < k.voids.nBasins; ++b)
      seeds = seeds && g.voidIndex[(std::size_t)internal::box_cell(g, &k.voids.centre[3*b])] == (std::int64_t)b;
    for (std::size_t b = 0; b < k.overdense.nBasins; ++b)
      seeds = seeds &&
              g.overdenseIndex[(std::size_t)internal::box_cell(g, &k.overdense.centre[3*b])] == (std::int64_t)b;
    check(seeds, "every centre in a cell of its basin");

    const TracerAssignment& a = *k.assignment;
    bool assigned = a.cellIndex.size() == kN && a.distance.size() == kN;
    for (std::size_t i = 0; assigned && i < kN; ++i) {
      const std::int64_t c = internal::box_cell(g, &kept.tracers[3*i]);
      const std::int64_t v = a.voidIndex[i], o = a.overdenseIndex[i];
      const double* x = v >= 0 ? &k.voids.centre[3*(std::size_t)v] : &k.overdense.centre[3*(std::size_t)o];
      const double dx = kept.tracers[3*i] - x[0];
      const double dy = kept.tracers[3*i+1] - x[1];
      const double dz = kept.tracers[3*i+2] - x[2];
      assigned = a.cellIndex[i] == c && v == g.voidIndex[(std::size_t)c] && o == g.overdenseIndex[(std::size_t)c] &&
                 a.distance[i] == std::sqrt(dx*dx + dy*dy + dz*dz) && a.inMask[i] == 1 && a.inBox[i] == 1 &&
                 a.used[i] == 1;
    }
    check(assigned, "every tracer's cell, basin and distance");

    Input named = kept;
    for (std::size_t i = 0; i < kN; ++i) named.ids.push_back((std::int64_t)(1000 - 3 * i));
    const Result n = named.run();
    check(n.ids == named.ids && catalogue_hash(n) == catalogue_hash(r), "ids copied, no catalogue changed");
  }

  group("both overloads give the same Result from the same arrays");
  {
    Input in;
    in.config.keepDivergence = true;
    in.config.keepAssignment = true;
    const otswap::Result rec = in.reconstruction();
    in.limits = rec.box;
    check(result_hash(in.run()) == result_hash(findVoidsBox(rec, in.config)), "equal hashes");
  }

  group("a reused grid: the same catalogues, the tracers outside it not read, the cell size recorded");
  {
    Input in;
    in.config.keepDivergence = true;
    in.config.keepAssignment = true;
    const Result first = in.run();

    const std::size_t extra = 5;
    Input more = in;
    more.tracers.clear();
    more.displacement.clear();
    for (std::size_t i = 0; i < kN + extra; ++i)
      for (std::size_t a = 0; a < 3; ++a) more.tracers.push_back(i < kN ? in.tracers[3*i+a] : 150. + i);
    for (std::size_t r = 0; r < kR; ++r)
      for (std::size_t i = 0; i < kN + extra; ++i)
        for (std::size_t a = 0; a < 3; ++a)
          more.displacement.push_back(i < kN ? in.displacement[3*(r*kN+i)+a] : 1.);
    more.valid.clear();
    more.config.keepDivergence = false;
    more.grid = &*first.divergence;
    const Result r = more.run();

    std::vector<std::uint8_t> outside(kN + extra, 0);
    for (std::size_t i = kN; i < kN + extra; ++i) outside[i] = 1;
    check(catalogue_hash(r) == catalogue_hash(first), "the catalogues of the grid's own search");
    check(r.outsideGrid == outside && r.counts.tracersOutsideGrid == extra && r.counts.tracersUsed == kN,
          "outsideGrid and the counts");
    check(r.config.divergenceCellSize == first.divergence->cellSizeMps, "the cell size of the grid, in mps");

    Input fewer = in;
    fewer.displacement.resize(3 * kN);
    fewer.valid.clear();
    fewer.config.keepDivergence = true;
    fewer.config.nRealizationsMean = 7;
    fewer.grid = &*first.divergence;
    const Result f = fewer.run();
    check(first.divergence->nRealizations == kR && first.divergence->nRealizationsMean == 1 &&
          f.divergence->nRealizations == 1 && f.divergence->nRealizationsMean == 7,
          "a reused grid records the realizations of the search reusing it, not those of the grid");

    DivergenceGrid skew = *first.divergence;
    const double c = 1.001799;
    skew.mps = 31.94434690031884;
    check(skew.mps > first.divergence->cellSize && (c * skew.mps) / skew.mps != c,
          "a cell size whose cellSize / mps differs from it, on a grid larger than the first");
    skew.cellSizeMps = c;
    skew.cellSize = c * skew.mps;
    Input onSkew = in;
    onSkew.mps = skew.mps;
    onSkew.config.keepDivergence = false;
    onSkew.grid = &skew;
    const Result s = onSkew.run();
    check(s.config.divergenceCellSize == c && *s.config.divergenceCellSize != skew.cellSize / skew.mps,
          "cellSizeMps recorded, not cellSize / mps, where the two differ");
    const TracerAssignment& a = *r.assignment;
    bool unread = true;
    for (std::size_t i = kN; i < kN + extra; ++i)
      unread = unread && a.cellIndex[i] == -1 && a.voidIndex[i] == -1 && a.overdenseIndex[i] == -1 &&
               std::isnan(a.distance[i]) && a.used[i] == 0;
    check(unread, "the tracers outside it: no cell, no basin, distance NaN, not used");
    bool same = true;
    for (std::size_t i = 0; i < kN; ++i)
      same = same && a.cellIndex[i] == first.assignment->cellIndex[i] &&
             a.distance[i] == first.assignment->distance[i];
    check(same, "the tracers inside it as in the grid's own search");
  }

  group("the limits recorded: on the grid and the Result, as given, otswap's, none, the call's on reuse");
  {
    Input in;
    in.config.keepDivergence = true;
    in.limits = box(-2.5, 102.5);
    const Result given = in.run();
    check(given.limits && given.limits->lower == in.limits->lower && given.limits->upper == in.limits->upper &&
          given.divergence->limits && given.divergence->limits->lower == in.limits->lower &&
          given.divergence->limits->upper == in.limits->upper, "the array limits");

    const otswap::Result rec = in.reconstruction();
    const Result fromRec = findVoidsBox(rec, in.config);
    const otswap::BoxLimits& sel = *rec.selection.boxLimits;
    check(fromRec.limits && fromRec.limits->lower == sel.lower && fromRec.limits->upper == sel.upper &&
          fromRec.divergence->limits && fromRec.divergence->limits->lower == sel.lower &&
          fromRec.divergence->limits->upper == sel.upper, "otswap::SelectionCounts::boxLimits");

    Input none = in;
    none.limits.reset();
    const Result without = none.run();
    check(!without.limits && !without.divergence->limits, "empty without limits");

    Input reuse = none;
    reuse.grid = &*without.divergence;
    reuse.limits = box(1., 99.);
    const Result reused = reuse.run();
    check(reused.limits && reused.limits->lower == reuse.limits->lower &&
          reused.limits->upper == reuse.limits->upper && reused.divergence->limits &&
          reused.divergence->limits->lower == reuse.limits->lower &&
          reused.divergence->limits->upper == reuse.limits->upper, "the call's limits on a reused grid");
    Input reuseNone = reuse;
    reuseNone.grid = &*given.divergence;
    reuseNone.limits.reset();
    const Result reusedNone = reuseNone.run();
    check(!reusedNone.limits && !reusedNone.divergence->limits, "none on a grid built with limits, reused without");
  }

  group("the same Result with 1 and 4 threads");
  {
    Input in;
    in.config.keepDivergence = true;
    in.config.keepAssignment = true;
    threads(1);
    const std::uint64_t one = result_hash(in.run());
    threads(4);
    const std::uint64_t four = result_hash(in.run());
    const otswap::Result rec = in.reconstruction();
    threads(1);
    const std::uint64_t recOne = result_hash(findVoidsBox(rec, in.config));
    threads(4);
    const std::uint64_t recFour = result_hash(findVoidsBox(rec, in.config));
    check(one == four, "arrays");
    check(recOne == recFour, "reconstruction");
  }

  return report("test_checks");
}
