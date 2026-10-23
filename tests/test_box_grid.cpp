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
 *  @brief The box grid as DivergenceGrid describes it, and its
 *  independence of the number of threads.
 */

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "internal.h"
#include "check.h"

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
    std::vector<double> tracers, displacement, randoms;
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
      return internal::check_box(tracers, displacement, valid, randoms, kMps, config, limits, {},
                                 nullptr);
    }
  };

  /// Extent of the points read: the tracers and, given none, the
  /// randoms, else the displacement endpoints.
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
    if (!in.randoms.empty())
      for (std::size_t i = 0; i < in.randoms.size() / 3; ++i) add(&in.randoms[3*i]);
    else
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

  group("the grid on the tracers and the randoms");
  {
    Input in;
    in.randoms = {-20., 10., 10., 50., 80., 40.};
    const DivergenceGrid g = internal::box_grid(in.search(), 0.5);
    double lo[3], hi[3];
    extent(in, lo, hi);
    check(lo[0] == -20. && hi[1] == 80., "the randoms extend the extent");
    check_placement(g, lo, hi, "randoms");
  }

  group("with limits the grid is built on them, and overhangs them");
  {
    Input in;
    otswap::BoxLimits box;
    box.lower = {10., 10., 10.};
    box.upper = {50., 40., 20.};
    in.limits = box;
    in.randoms = {15., 15., 15.};
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

  return report("test_box_grid");
}
