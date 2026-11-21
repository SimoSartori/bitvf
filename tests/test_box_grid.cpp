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
 *  @file tests/test_box_grid.cpp
 *
 *  @brief The box grid as DivergenceGrid describes it, boxGrid, and
 *  their independence of the number of threads.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "otswap/Maths.h"
#include "otswap/OT.h"

#include "internal.h"
#include "check.h"
#include "fixtures.h"

using namespace bitvf;

namespace {

  constexpr std::size_t kN = 2000;
  constexpr double kMps = 3.;

  // The hash of the grid and the tracer cells of the base input, recorded
  // once. A different value on some platform means results are not
  // portable there.
  constexpr std::uint64_t kRecordedHash = 0xe03897a49b9e7605ull;

  /// Tracers uniform in [0, 100) x [0, 60) x [0, 30), so that the axes need
  /// different counts; one realization of displacements within 4 of zero.
  struct Input {
    std::vector<double> tracers, displacement;
    std::vector<std::uint8_t> valid;
    std::optional<otswap::BoxLimits> limits;
    Config config;

    Input ()
    {
      std::mt19937 rng(7);
      auto uniform = [&rng] (const double lo, const double hi) {
        return lo + (hi - lo) * (rng() / 4294967296.);
      };
      const double side[3] = {100., 60., 30.};
      for (std::size_t i = 0; i < kN; ++i)
        for (std::size_t k = 0; k < 3; ++k) tracers.push_back(uniform(0., side[k]));
      for (std::size_t i = 0; i < 3 * kN; ++i) displacement.push_back(uniform(-4., 4.));
      config.log = nullptr;
    }

    internal::BoxSearch search () const
    {
      return internal::check_box(tracers, displacement, valid, kMps, config, limits, {},
                                 nullptr);
    }
  };

  /// Extent of the points read: the tracers and the displacement
  /// endpoints.
  void extent (const Input& in, double lo[3], double hi[3])
  {
    for (std::size_t k = 0; k < 3; ++k) { lo[k] = 1e300; hi[k] = -1e300; }
    auto add = [&] (const double* x) {
      for (std::size_t k = 0; k < 3; ++k) {
        lo[k] = std::fmin(lo[k], x[k]);
        hi[k] = std::fmax(hi[k], x[k]);
      }
    };
    for (std::size_t i = 0; i < kN; ++i) add(&in.tracers[3*i]);
    for (std::size_t i = 0; i < kN; ++i) {
      const double end[3] = {in.tracers[3*i] + in.displacement[3*i],
                             in.tracers[3*i+1] + in.displacement[3*i+1],
                             in.tracers[3*i+2] + in.displacement[3*i+2]};
      add(end);
    }
  }

  /// The placement of DivergenceGrid: one count, floor(E / s + 1) for the
  /// largest padded extent E, lower down by S / 2, upper up by S / 4.
  void check_placement (const DivergenceGrid& g, const double lo[3], const double hi[3],
                        const std::string& what)
  {
    const double s = g.cellSize;
    double largest = 0.;
    for (std::size_t k = 0; k < 3; ++k)
      largest = std::fmax(largest, (hi[k] - lo[k]) + 0.2 * kMps);
    const std::size_t n = (std::size_t)std::floor(largest / s + 1);

    check(g.nCellsAxis[0] == n && g.nCellsAxis[1] == n && g.nCellsAxis[2] == n,
          what + ": one count, from the largest extent");
    check(g.nCells == n * n * n, what + ": nCells is the cube of the count");

    for (std::size_t k = 0; k < 3; ++k) {
      const double E = (hi[k] - lo[k]) + 0.2 * kMps;
      const double S = n * s - E;
      check(std::fabs((lo[k] - 0.1 * kMps) - g.lower[k] - S / 2) < 1e-9,
            what + ": lower moves down by half the slack");
      check(std::fabs(g.upper[k] - (hi[k] + 0.1 * kMps) - S / 4) < 1e-9,
            what + ": upper moves up by a quarter of the slack");
      check(g.upper[k] < g.lower[k] + n * s, what + ": upper lies below the last cell face");
    }
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

  /// The hash of the grid, its centres and the cell of each tracer.
  std::uint64_t grid_hash (const Input& in)
  {
    DivergenceGrid g = internal::box_grid(in.search(), 1.);
    internal::box_centres(g);
    std::vector<std::int64_t> cells(kN);
    for (std::size_t i = 0; i < kN; ++i) cells[i] = internal::box_cell(g, &in.tracers[3*i]);

    std::uint64_t h = 0xcbf29ce484222325ull;
    hash_bytes(h, &g.cellSize, sizeof g.cellSize);
    hash_bytes(h, g.lower.data(), sizeof g.lower);
    hash_bytes(h, g.upper.data(), sizeof g.upper);
    hash_bytes(h, g.nCellsAxis.data(), sizeof g.nCellsAxis);
    hash_bytes(h, g.centre.data(), g.centre.size() * sizeof(double));
    hash_bytes(h, cells.data(), cells.size() * sizeof(std::int64_t));
    return h;
  }

}

int main ()
{
  group("the grid on the tracers and the displacement endpoints");
  {
    const Input in;
    const DivergenceGrid g = internal::box_grid(in.search(), 1.);
    double lo[3], hi[3];
    extent(in, lo, hi);
    check(g.geometry == otswap::Geometry::Box && g.mps == kMps && g.cellSize == kMps,
          "geometry, mps, cell side = divergenceCellSize * mps");
    check_placement(g, lo, hi, "endpoints");
  }

  group("without limits an endpoint beyond the tracers widens the grid");
  {
    Input in;
    in.displacement[0] = -20. - in.tracers[0];
    const DivergenceGrid g = internal::box_grid(in.search(), 0.5);
    double lo[3], hi[3];
    extent(in, lo, hi);
    check(lo[0] == -20., "the endpoint sets the lower extent");
    check_placement(g, lo, hi, "endpoint beyond the tracers");
  }

  group("with limits the grid is built on them, and overhangs them");
  {
    Input in;
    otswap::BoxLimits box;
    box.lower = {10., 10., 10.};
    box.upper = {50., 40., 20.};
    in.limits = box;
    const DivergenceGrid g = internal::box_grid(in.search(), 1.);
    const double lo[3] = {10., 10., 10.}, hi[3] = {50., 40., 20.};
    check_placement(g, lo, hi, "limits");
    for (std::size_t k = 0; k < 3; ++k)
      check(g.lower[k] < box.lower[k] && g.upper[k] > box.upper[k], "the grid overhangs");
  }

  group("cell rows, centres and the cell of every tracer");
  {
    const Input in;
    DivergenceGrid g = internal::box_grid(in.search(), 1.);
    internal::box_centres(g);
    const std::size_t n = g.nCellsAxis[0];
    const double s = g.cellSize;

    check(g.centre.size() == 3 * g.nCells, "one centre per cell");
    bool rows = true, centres = true;
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        for (std::size_t k = 0; k < n; ++k) {
          const std::size_t c = (i * n + j) * n + k;
          centres = centres && g.centre[3*c] == g.lower[0] + s * (i + 0.5) &&
                    g.centre[3*c+1] == g.lower[1] + s * (j + 0.5) &&
                    g.centre[3*c+2] == g.lower[2] + s * (k + 0.5);
          rows = rows && internal::box_cell(g, &g.centre[3*c]) == (std::int64_t)c;
        }
    check(centres, "the centre of cell (i, j, k) is lower + (i, j, k, + 0.5) * cellSize");
    check(rows, "cell (i, j, k) is row (i * nY + j) * nZ + k");

    bool inside = true;
    for (std::size_t t = 0; t < kN; ++t) {
      const std::int64_t c = internal::box_cell(g, &in.tracers[3*t]);
      if (c < 0) { inside = false; continue; }
      const std::size_t idx[3] = {(std::size_t)c / (n * n), ((std::size_t)c / n) % n,
                                  (std::size_t)c % n};
      for (std::size_t k = 0; k < 3; ++k) {
        const double x = in.tracers[3*t+k];
        inside = inside && x >= g.lower[k] + s * idx[k] - 1e-9 &&
                 x <= g.lower[k] + s * (idx[k] + 1) + 1e-9;
      }
    }
    check(inside, "every tracer lies in the cell box_cell gives");

    const double below[3] = {g.lower[0] - 1e-6, g.centre[1], g.centre[2]};
    const double top[3] = {g.lower[0] + n * s, g.centre[1], g.centre[2]};
    check(internal::box_cell(g, below) == -1 && internal::box_cell(g, top) == -1,
          "points below the lower face or on the last face are outside");
  }

  group("the same grid with 1 and 4 threads, and the recorded hash");
  {
    const Input in;
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = grid_hash(in);
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = grid_hash(in);
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash)
      std::cerr << "grid hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  group("boxGrid on limits: the grid of a search with the same limits, mps and cell size, bit for bit; only "
        "the fields that define the cells set");
  const fixtures::Translated tr;
  for (const double cell : {1., 0.7}) {
    Config c = tr.config;
    c.divergenceCellSize = cell;
    const Result r = tr.run(c);
    const DivergenceGrid& d = *r.divergence;
    const DivergenceGrid g = boxGrid(*tr.limits, fixtures::Translated::kMps, cell);
    const std::string what = "cell " + std::to_string(cell) + " mps: ";
    check(g.geometry == otswap::Geometry::Box && g.lower == d.lower && g.upper == d.upper &&
          g.nCellsAxis == d.nCellsAxis && g.nCells == d.nCells && g.cellSize == d.cellSize &&
          g.cellSizeMps == d.cellSizeMps && g.mps == d.mps, what + "the defining fields");
    check(!g.limits && g.nRealizations == 0 && g.nRealizationsMean == 0 && g.divergence.empty() &&
          g.weight.empty() && g.nTracers.empty() && g.nDisplacements.empty() && g.voidIndex.empty() &&
          g.centre.empty() && g.volume.empty(), what + "no record, no per-cell field");
  }

  group("boxGrid from a reconstruction: the grid findVoidsBox builds on it, with the default and a given cell "
        "size; the default min(1, cbrt(20 / nRealizations))");
  otswap::Config ot;
  ot.nRealizations = 3;
  ot.seed = 1;
  ot.log = nullptr;
  std::vector<double> few(tr.tracers.begin(), tr.tracers.begin() + 3 * 400);
  otswap::BoxLimits limits;
  limits.lower = {-2., -2., -2.};
  limits.upper = {102., 102., 102.};
  const otswap::Result rec = otswap::reconstructBox(few, {}, ot, limits);
  {
    Config c = tr.config;
    c.keepAssignment = false;
    c.nRealizationsMean = 1;
    for (const std::optional<double> cell : {std::optional<double>(), std::optional<double>(0.8)}) {
      c.divergenceCellSize = cell;
      const Result searched = findVoidsBox(rec, c);
      const DivergenceGrid& d = *searched.divergence;
      const DivergenceGrid g = boxGrid(rec, cell);
      check(g.lower == d.lower && g.upper == d.upper && g.nCellsAxis == d.nCellsAxis && g.nCells == d.nCells &&
            g.cellSize == d.cellSize && g.cellSizeMps == d.cellSizeMps && g.mps == d.mps && g.mps == rec.mps &&
            g.divergence.empty() && !g.limits,
            cell ? "cell 0.8 mps" : "default cell, 3 realizations: 1 mps");
    }
    otswap::Result many = rec;
    many.nRealizations = 40;
    const DivergenceGrid g = boxGrid(many);
    const DivergenceGrid h = boxGrid(*rec.box, rec.mps, otswap::maths::cbrt(20. / 40.));
    check(g.cellSizeMps == otswap::maths::cbrt(0.5) && g.cellSize == h.cellSize && g.lower == h.lower &&
          g.nCellsAxis == h.nCellsAxis, "40 realizations: cbrt(20 / 40) mps, the limits form's grid");
  }

  group("reuse: a grid from boxGrid passes the reuse checks, and the searches reusing it give the catalogues of "
        "the searches that built their own, both kinds of basin");
  {
    const Result own = tr.run();
    const DivergenceGrid g = boxGrid(*tr.limits, fixtures::Translated::kMps, *own.config.divergenceCellSize);
    const Result reused = tr.run(&g);
    check(own.overdense.nBasins > 0 && fixtures::same_catalog(reused.voids, own.voids) &&
          fixtures::same_catalog(reused.overdense, own.overdense), "arrays with limits");
    Config c = tr.config;
    c.keepAssignment = false;
    c.nRealizationsMean = 1;
    const DivergenceGrid gr = boxGrid(rec);
    const Result a = findVoidsBox(rec, c);
    const Result b = findVoidsBox(rec, c, &gr);
    check(fixtures::same_catalog(a.voids, b.voids) && fixtures::same_catalog(a.overdense, b.overdense),
          "a reconstruction");
  }

  group("errors of the box grid builders");
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    otswap::BoxLimits bad = limits;
    bad.lower[1] = nan;
    throws_naming([&] { boxGrid(bad, 1., 1.); }, "the box limits hold a bound of", "not finite");
    bad = limits;
    bad.upper[2] = bad.lower[2];
    throws_naming([&] { boxGrid(bad, 1., 1.); }, "the box limits are [", "lower not below upper");
    throws_naming([&] { boxGrid(limits, -1., 1.); }, "the mean particle separation is ", "the mean particle separation is ");
    throws_naming([&] { boxGrid(limits, 1., 0.); }, "the divergence cell size is ", "the divergence cell size is ");
    auto from = [&rec] (auto&& change) {
      otswap::Result r = rec;
      change(r);
      return r;
    };
    throws_naming([&] { boxGrid(from([] (otswap::Result& r) { r.geometry = otswap::Geometry::Lightcone; })); },
                  "boxGrid: the reconstruction is not a box result; use lightconeGrid",
                  "boxGrid: the reconstruction is not a box result");
    throws_naming([&] { boxGrid(from([] (otswap::Result& r) { r.nRealizations = 0; })); },
                  "boxGrid: the reconstruction is malformed: it holds no realizations", "boxGrid: the reconstruction is malformed: it holds no realizations");
    throws_naming([&] { boxGrid(from([] (otswap::Result& r) { r.box.reset(); })); },
                  "boxGrid: the reconstruction is malformed: it holds no box", "no box");
    throws_naming([&] { boxGrid(from([&nan] (otswap::Result& r) { r.box->upper[0] = nan; })); },
                  "has a bound of", "not finite");
    throws_naming([&] { boxGrid(from([] (otswap::Result& r) { r.box->lower[2] = r.box->upper[2] + 1.;
                                                               r.selection.boxLimits.reset(); })); },
                  "the lower bound must not be above the upper one", "lower above upper");
    throws_naming([&] { boxGrid(from([] (otswap::Result& r) { r.box->upper[1] += 1.; })); },
                  "differs from its box limits", "differs from the limits");
    throws_naming([&] { boxGrid(from([&nan] (otswap::Result& r) { r.mps = nan; })); },
                  "boxGrid: the reconstruction's mps is ", "boxGrid: the reconstruction's mps is ");
    throws_naming([&] { boxGrid(rec, -2.); }, "the divergence cell size is ", "the divergence cell size is ");
    otswap::BoxLimits huge;
    huge.lower = {0., 0., 0.};
    huge.upper = {1e7, 1., 1.};
    throws_naming([&] { boxGrid(huge, 1., 1.); }, "the grid would need ", "limits form");
    throws_naming([&] { boxGrid(from([] (otswap::Result& r) { r.mps = 1e-6; })); },
                  "the grid would need ", "reconstruction form");
  }

  return report("test_box_grid");
}
