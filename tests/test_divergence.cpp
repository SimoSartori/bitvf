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
 *  @file tests/test_divergence.cpp
 *
 *  @brief The divergence of a box search: counts, terms, weights, the
 *  fill, its errors, the scale of the accumulation, and its independence
 *  of the number of threads.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
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

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

  // The hashes of the divergence of the generated input, cells of 1 and 0.5
  // mps, recorded once. A different value on some platform means results
  // are not portable there.
  constexpr std::uint64_t kRecordedHash1 = 0x562871ab2052e046ull;
  constexpr std::uint64_t kRecordedHash05 = 0xe535a4e780425c9bull;

  void threads (const int n)
  {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
  }

  /// A box grid with lower corner 0, n cells of side s along each axis, mps
  /// 1; the per-cell arrays empty.
  DivergenceGrid cube (const std::size_t n, const double s)
  {
    DivergenceGrid g;
    g.geometry = otswap::Geometry::Box;
    g.mps = 1.;
    g.cellSize = s;
    g.cellSizeMps = s;
    g.lower = {0., 0., 0.};
    g.upper = {n * s, n * s, n * s};
    g.nCellsAxis = {n, n, n};
    g.nCells = n * n * n;
    return g;
  }

  std::size_t row (const DivergenceGrid& g, const std::size_t i, const std::size_t j,
                   const std::size_t k)
  {
    return (i * g.nCellsAxis[1] + j) * g.nCellsAxis[2] + k;
  }

  using Vec = std::array<double, 3>;

  /// Tracers and their displacements, one per realization, searched on a
  /// reused grid through check_box; a NaN displacement is invalid.
  struct Case {
    DivergenceGrid grid;
    unsigned R;
    std::vector<double> tracers, displacement;
    std::vector<std::uint8_t> valid;
    std::vector<std::int64_t> ids;
    std::vector<Vec> x;
    std::vector<std::vector<Vec>> d;
    Config config;

    Case (const std::size_t n, const double s, const unsigned R = 1) : grid(cube(n, s)), R(R)
    {
      config.log = nullptr;
    }

    void add (const Vec& position, const std::vector<Vec>& disp)
    {
      x.push_back(position);
      d.push_back(disp);
    }

    internal::BoxSearch search ()
    {
      const std::size_t N = x.size();
      tracers.clear();
      displacement.assign(3 * R * N, 0.);
      valid.assign(R * N, 1);
      for (std::size_t i = 0; i < N; ++i) {
        tracers.insert(tracers.end(), x[i].begin(), x[i].end());
        for (std::size_t r = 0; r < R; ++r) {
          for (std::size_t k = 0; k < 3; ++k) displacement[3 * (r * N + i) + k] = d[i][r][k];
          if (std::isnan(d[i][r][0])) valid[r * N + i] = 0;
        }
      }
      return internal::check_box(tracers, displacement, valid, 1., config, std::nullopt, ids,
                                 &grid);
    }

    /// The divergence on a copy of the grid, left in out.
    internal::BoxDivergence run (DivergenceGrid& out)
    {
      const internal::BoxSearch s = search();
      out = grid;
      return internal::box_divergence(s, out);
    }
  };

  const Vec kInvalid = {kNaN, kNaN, kNaN};

  /// The 1 / r^2 means of value and weight over the cells flagged in
  /// source within h cells of cell c along each axis, in ascending row
  /// order; false when fewer than three are flagged.
  bool window (const DivergenceGrid& g, const std::vector<double>& value,
               const std::vector<double>& weight, const std::vector<std::uint8_t>& source,
               const std::size_t c, const std::int64_t h, double& v, double& w)
  {
    const std::int64_t n = (std::int64_t)g.nCellsAxis[0];
    const std::int64_t i = (std::int64_t)c / (n * n), j = ((std::int64_t)c / n) % n,
                       k = (std::int64_t)c % n;
    const double s = g.cellSize;
    double num = 0., den = 0., wn = 0.;
    unsigned count = 0;
    for (std::int64_t ii = i - h; ii <= i + h; ++ii)
      for (std::int64_t jj = j - h; jj <= j + h; ++jj)
        for (std::int64_t kk = k - h; kk <= k + h; ++kk) {
          if (ii < 0 || jj < 0 || kk < 0 || ii >= n || jj >= n || kk >= n) continue;
          const std::size_t cc = (std::size_t)((ii * n + jj) * n + kk);
          if (!source[cc]) continue;
          ++count;
          const double inv =
            1. / ((double)((i - ii) * (i - ii) + (j - jj) * (j - jj) + (k - kk) * (k - kk)) * s * s);
          num += value[cc] * inv;
          den += inv;
          wn += weight[cc] * inv;
        }
    if (count < 3) return false;
    v = num / den;
    w = wn / den;
    return true;
  }

  /// A grid for box_fill: every per-cell array zero.
  DivergenceGrid empty_grid (const std::size_t n)
  {
    DivergenceGrid g = cube(n, 1.);
    g.nDisplacements.assign(g.nCells, 0);
    g.divergence.assign(g.nCells, 0.);
    g.weight.assign(g.nCells, 0.);
    return g;
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

  /// Generated input: tracers uniform in [0, 100)^3, R = 4, displacements
  /// uniform within 6 of zero, about 5% invalid; limits the bounding box
  /// of the tracers and of randoms uniform in [-8, 108)^3, which covers
  /// every endpoint.
  struct Generated {
    static constexpr std::size_t N = 20000;
    static constexpr unsigned R = 4;
    std::vector<double> tracers, displacement, randoms;
    std::vector<std::uint8_t> valid;
    std::vector<std::int64_t> ids;
    otswap::BoxLimits limits;
    Config config;

    Generated ()
    {
      std::mt19937 rng(11);
      auto uniform = [&rng] (const double lo, const double hi) {
        return lo + (hi - lo) * (rng() / 4294967296.);
      };
      for (std::size_t i = 0; i < 3 * N; ++i) tracers.push_back(uniform(0., 100.));
      for (std::size_t i = 0; i < 3 * R * N; ++i) displacement.push_back(uniform(-6., 6.));
      for (std::size_t i = 0; i < R * N; ++i) valid.push_back(rng() % 20 != 0);
      for (std::size_t i = 0; i < 3 * N; ++i) randoms.push_back(uniform(-8., 108.));
      limits.lower = {tracers[0], tracers[1], tracers[2]};
      limits.upper = limits.lower;
      for (const std::vector<double>* points : {&tracers, &randoms})
        for (std::size_t i = 0; i < points->size(); ++i) {
          limits.lower[i%3] = std::min(limits.lower[i%3], (*points)[i]);
          limits.upper[i%3] = std::max(limits.upper[i%3], (*points)[i]);
        }
      config.log = nullptr;
    }

    /// The hash of the divergence on cells of cellSizeMps; the passes of
    /// the fill in passes.
    std::uint64_t hash (const double cellSizeMps, unsigned& passes) const
    {
      const internal::BoxSearch s = internal::check_box(tracers, displacement, valid, 4., config,
                                                        limits, ids, nullptr);
      DivergenceGrid g = internal::box_grid(s, cellSizeMps);
      const internal::BoxDivergence b = internal::box_divergence(s, g);
      passes = b.fillPasses;

      std::uint64_t h = 0xcbf29ce484222325ull;
      hash_bytes(h, g.nTracers.data(), g.nTracers.size() * sizeof(std::uint64_t));
      hash_bytes(h, g.nDisplacements.data(), g.nDisplacements.size() * sizeof(std::uint64_t));
      hash_bytes(h, g.divergence.data(), g.divergence.size() * sizeof(double));
      hash_bytes(h, g.weight.data(), g.weight.size() * sizeof(double));
      hash_bytes(h, b.cell.data(), b.cell.size() * sizeof(std::int64_t));
      hash_bytes(h, b.used.data(), b.used.size());
      hash_bytes(h, &b.scaleExponent, sizeof b.scaleExponent);
      hash_bytes(h, &b.fillPasses, sizeof b.fillPasses);
      return h;
    }
  };

}

int main ()
{
  threads(4);

  group("counts: a displacement leaving its cell counts 1 + 1, one ending in it 2; nTracers, used, cell");
  {
    Case c(4, 1., 2);
    c.add({0.5, 0.5, 0.5}, {{1., 0., 0.}, {0., 0.2, 0.}});
    c.add({0.25, 0.75, 0.5}, {kInvalid, kInvalid});
    c.add({2.5, 2.5, 2.5}, {{0., 0., -1.}, kInvalid});
    c.add({5., 5., 5.}, {{0.1, 0., 0.}, {0.1, 0., 0.}});
    DivergenceGrid g;
    const internal::BoxDivergence b = c.run(g);

    std::vector<std::uint64_t> displacements(g.nCells, 0), tracers(g.nCells, 0);
    displacements[row(g, 0, 0, 0)] = 3;
    displacements[row(g, 1, 0, 0)] = 1;
    displacements[row(g, 2, 2, 2)] = 1;
    displacements[row(g, 2, 2, 1)] = 1;
    tracers[row(g, 0, 0, 0)] = 2;
    tracers[row(g, 2, 2, 2)] = 1;
    check(g.nDisplacements == displacements, "nDisplacements: 1 + 1 leaving, 2 staying");
    check(g.nTracers == tracers, "nTracers counts the tracers read in their start cell");
    check(b.used == std::vector<std::uint8_t>({1, 0, 1, 0}) && b.nUsed == 2,
          "used: a valid displacement read; nUsed");
    check(b.cell == std::vector<std::int64_t>({(std::int64_t)row(g, 0, 0, 0),
                                               (std::int64_t)row(g, 0, 0, 0),
                                               (std::int64_t)row(g, 2, 2, 2), -1}),
          "cell: the start cell, -1 for a tracer outside the grid");
  }

  group("divergence: the outward component over the cell side, divided by the count; "
        "intermediate cells receive nothing");
  {
    Case c(10, 2.);
    auto at = [] (const double i, const double j, const double k) {
      return Vec{2. * i + 1., 2. * j + 1., 2. * k + 1.};
    };
    c.add(at(1, 1, 1), {{1.5, 0., 0.}});
    c.add(at(1, 5, 1), {{-1.5, 0., 0.}});
    c.add(at(5, 1, 1), {{0., 1.5, 0.}});
    c.add(at(5, 5, 1), {{0., -1.5, 0.}});
    c.add(at(1, 1, 5), {{0., 0., 1.5}});
    c.add(at(5, 5, 5), {{0., 0., -1.5}});
    c.add(at(8, 8, 8), {{1.5, 0., 0.}});
    c.add(at(8, 8, 8), {{0., 0.3, 0.}});
    c.add(at(1, 8, 1), {{4.5, 0., 0.}});
    DivergenceGrid g;
    const internal::BoxDivergence b = c.run(g);

    const std::size_t from[6][3] = {{1, 1, 1}, {1, 5, 1}, {5, 1, 1}, {5, 5, 1}, {1, 1, 5}, {5, 5, 5}};
    const std::size_t to[6][3] = {{2, 1, 1}, {0, 5, 1}, {5, 2, 1}, {5, 4, 1}, {1, 1, 6}, {5, 5, 4}};
    bool faces = true;
    for (std::size_t t = 0; t < 6; ++t) {
      const std::size_t a = row(g, from[t][0], from[t][1], from[t][2]);
      const std::size_t e = row(g, to[t][0], to[t][1], to[t][2]);
      faces = faces && g.divergence[a] == 0.75 && g.divergence[e] == -0.75 &&
              g.nDisplacements[a] == 1 && g.nDisplacements[e] == 1;
    }
    check(faces, "each axis and sign: +|d|/s in the start cell, -|d|/s in the end cell");
    check(g.divergence[row(g, 8, 8, 8)] == 0.25 && g.nDisplacements[row(g, 8, 8, 8)] == 3,
          "a displacement ending in its cell adds nothing and counts 2: (0.75 + 0) / 3");
    check(g.divergence[row(g, 1, 8, 1)] == 2.25 && g.divergence[row(g, 3, 8, 1)] == -2.25 &&
          g.nDisplacements[row(g, 2, 8, 1)] == 0,
          "a crossed intermediate cell counts nothing");
    int e = 0;
    std::frexp(4.5 / 2., &e);
    check(b.scaleExponent == 61 - e - 20, "q = 61 - e - 20, B = maxDisplacement / s < 2^e");
  }

  group("faces: a tracer or an end on a face, edge or vertex, or a displacement along a face, gives "
        "the cells and terms of the same tracer moved one ulp into its cell");
  {
    struct Face {
      const char* what;
      Vec x, d;
      std::array<bool, 3> move;
    };
    const Face faces[] = {
      {"start on a face, leaving through it", {2.5, 2.5, 2.}, {0.3, -0.2, -0.7}, {false, false, true}},
      {"start on a face, leaving through another", {2.5, 2.5, 2.}, {0.9, 0.1, 0.4}, {false, false, true}},
      {"start on an edge, y faster", {2.5, 2., 2.}, {0.1, -0.9, -0.3}, {false, true, true}},
      {"start on an edge, z faster", {2.5, 2., 2.}, {0.1, -0.3, -0.9}, {false, true, true}},
      {"start on a vertex, x fastest", {2., 2., 2.}, {-0.9, -0.3, -0.1}, {true, true, true}},
      {"start on a vertex, z fastest", {2., 2., 2.}, {-0.1, -0.3, -0.9}, {true, true, true}},
      {"start on a vertex, y into the cell", {2., 2., 2.}, {-0.9, 0.3, -0.1}, {true, true, true}},
      {"start on a face, along it", {2.5, 2.5, 2.}, {0.9, 0., 0.}, {false, false, true}},
      {"start on an edge, along it", {2.5, 2., 2.}, {0.9, 0., 0.}, {false, true, true}},
      {"start on an edge, along one face", {2.5, 2., 2.}, {0.9, -0.6, 0.}, {false, true, true}},
      {"end on a face", {4.5, 4.5, 4.5}, {1.5, 0.2, 0.1}, {true, false, false}},
      {"end on an edge", {5.5, 4.5, 4.5}, {0.5, 1.5, 0.1}, {true, true, false}},
    };
    const Vec stay = {0.01, 0., 0.};
    auto run = [&] (const Vec& x, const Vec& d, DivergenceGrid& g) {
      Case c(8, 1.);
      for (std::size_t i = 0; i < 8; ++i)
        for (std::size_t j = 0; j < 8; ++j)
          for (std::size_t k = 0; k < 8; ++k) c.add({i + 0.5, j + 0.5, k + 0.5}, {stay});
      c.add(x, {d});
      return c.run(g);
    };
    for (const Face& f : faces) {
      Vec moved = f.x;
      for (std::size_t a = 0; a < 3; ++a)
        if (f.move[a]) moved[a] = std::nextafter(moved[a], std::numeric_limits<double>::infinity());
      DivergenceGrid g, h;
      const internal::BoxDivergence b = run(f.x, f.d, g);
      const internal::BoxDivergence m = run(moved, f.d, h);
      check(g.nDisplacements == h.nDisplacements && g.nTracers == h.nTracers &&
            g.divergence == h.divergence && g.weight == h.weight && b.cell == m.cell &&
            b.used == m.used && b.scaleExponent == m.scaleExponent,
            std::string(f.what) + ": as the tracer moved one ulp into its cell");
    }

    DivergenceGrid g;
    run({2.5, 2.5, 2.5}, {1., 1., 0.1}, g);
    std::uint64_t total = 0;
    for (const std::uint64_t n : g.nDisplacements) total += n;
    check(total == 2 * 513 && g.nDisplacements[row(g, 2, 2, 2)] == 3 &&
          g.nDisplacements[row(g, 3, 3, 2)] == 3 && g.divergence[row(g, 2, 2, 2)] == 1. / 3. &&
          g.divergence[row(g, 3, 3, 2)] == -1. / 3.,
          "a segment through an edge: one face of each cell, the lower axis on equal t and d");

    struct Tiny {
      const char* what;
      Vec x, d;
      std::array<std::size_t, 3> from, to;
      int exit, entry;
    };
    for (const double t : {1e-10, 1e-12}) {
      const double in = 3. - 0.5 * t;
      const Tiny tiny[] = {
        {"on a face", {3., 2.5, 2.5}, {-t, 0.2 * t, 0.3 * t}, {3, 2, 2}, {2, 2, 2}, 0, 0},
        {"inside, across a face", {in, 2.5, 2.5}, {t, 0.2 * t, -0.3 * t}, {2, 2, 2}, {3, 2, 2}, 0, 0},
        {"on an edge", {3., 3., 2.5}, {-t, -2. * t, 0.3 * t}, {3, 3, 2}, {2, 2, 2}, 1, 0},
        {"inside, across an edge", {in, in, 2.5}, {t, 1.5 * t, 0.3 * t}, {2, 2, 2}, {3, 3, 2}, 1, 0},
        {"on a vertex", {3., 3., 3.}, {-t, -2. * t, -3. * t}, {3, 3, 3}, {2, 2, 2}, 2, 0},
        {"inside, across a vertex", {in, in, in}, {t, 2. * t, 3. * t}, {2, 2, 2}, {3, 3, 3}, 2, 0},
      };
      for (const Tiny& c : tiny) {
        DivergenceGrid h;
        const internal::BoxDivergence b = run(c.x, c.d, h);
        const std::size_t from = row(h, c.from[0], c.from[1], c.from[2]);
        const std::size_t to = row(h, c.to[0], c.to[1], c.to[2]);
        auto term = [&] (const int a) {
          return std::ldexp((double)std::llrint(std::ldexp(std::abs(c.d[a]), b.scaleExponent)),
                            -b.scaleExponent) / 3.;
        };
        std::uint64_t sum = 0;
        for (const std::uint64_t n : h.nDisplacements) sum += n;
        check(b.cell[512] == (std::int64_t)from && sum == 2 * 513 && h.nDisplacements[from] == 3 &&
              h.nDisplacements[to] == 3 && h.divergence[from] == term(c.exit) &&
              h.divergence[to] == -term(c.entry),
              "|d| = " + std::to_string(t) + ", " + c.what +
              ": one exit from the start cell, one entry into the end cell, each through one face");
      }
    }
  }

  group("weight: nDisplacements / R, with R = 3");
  {
    Case c(4, 1., 3);
    c.add({1.5, 1.5, 1.5}, {{1., 0., 0.}, {0.1, 0., 0.}, {0., 0.1, 0.}});
    c.add({0.5, 2.5, 0.5}, {{0.1, 0., 0.}, {0.1, 0., 0.}, {0.1, 0., 0.}});
    DivergenceGrid g;
    c.run(g);
    check(g.weight[row(g, 1, 1, 1)] == 5. / 3. && g.weight[row(g, 2, 1, 1)] == 1. / 3. &&
          g.weight[row(g, 0, 2, 0)] == 6. / 3.,
          "weight = n / R in every counted cell");
  }

  group("fill: the 1 / r^2 mean of the 3x3x3 window");
  {
    DivergenceGrid g = empty_grid(3);
    for (std::size_t c = 0; c < g.nCells; ++c) {
      g.nDisplacements[c] = 1;
      g.divergence[c] = 0.37 * c - 3.;
      g.weight[c] = 1. + 0.125 * c;
    }
    const std::size_t centre = row(g, 1, 1, 1);
    g.nDisplacements[centre] = 0;
    g.divergence[centre] = 0.;
    g.weight[centre] = 0.;
    std::vector<std::uint8_t> source(g.nCells, 1);
    source[centre] = 0;
    double v = 0., w = 0.;
    window(g, g.divergence, g.weight, source, centre, 1, v, w);
    const unsigned passes = internal::box_fill(g);
    check(passes == 1 && g.divergence[centre] == v && g.weight[centre] == w,
          "the empty centre takes the 1 / r^2 means of its 26 neighbours");
  }

  group("fill: the 5x5x5 window when the 3x3x3 one holds fewer than three values");
  {
    DivergenceGrid g = empty_grid(5);
    const std::size_t sources[3] = {row(g, 0, 0, 0), row(g, 4, 4, 4), row(g, 0, 4, 0)};
    std::vector<std::uint8_t> source(g.nCells, 0);
    for (std::size_t t = 0; t < 3; ++t) {
      g.nDisplacements[sources[t]] = 1;
      g.divergence[sources[t]] = 1. + t;
      g.weight[sources[t]] = 0.5 + t;
      source[sources[t]] = 1;
    }
    const std::size_t centre = row(g, 2, 2, 2);
    double v = 0., w = 0.;
    window(g, g.divergence, g.weight, source, centre, 2, v, w);
    internal::box_fill(g);
    check(g.divergence[centre] == v && g.weight[centre] == w,
          "the centre, with no value within one cell, takes the means over two cells");
  }

  group("fill: a pass reads only the values of the previous passes");
  {
    const std::size_t n = 7;
    DivergenceGrid g = empty_grid(n);
    std::vector<std::uint8_t> plane0(g.nCells, 0), plane2(g.nCells, 0);
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t k = 0; k < n; ++k) {
        const std::size_t c = row(g, 0, j, k);
        g.nDisplacements[c] = 1;
        g.divergence[c] = (double)(j * n + k) - 20.;
        g.weight[c] = 0.25 * (j + k + 1);
        plane0[c] = 1;
        plane2[row(g, 2, j, k)] = 1;
      }

    std::vector<double> value = g.divergence, weight = g.weight;
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t k = 0; k < n; ++k) {
        const std::size_t c = row(g, 2, j, k);
        window(g, g.divergence, g.weight, plane0, c, 2, value[c], weight[c]);
      }
    const std::size_t target = row(g, 3, 3, 3);
    double v = 0., w = 0.;
    window(g, value, weight, plane2, target, 1, v, w);
    double v1 = 0., w1 = 0.;
    window(g, g.divergence, g.weight, plane0, row(g, 1, 3, 3), 1, v1, w1);

    const unsigned passes = internal::box_fill(g);
    check(passes == 3, "the planes fill in three passes");
    check(g.divergence[row(g, 1, 3, 3)] == v1 && g.weight[row(g, 1, 3, 3)] == w1,
          "pass 1, plane 1: the 3x3x3 window of plane 0");
    check(g.divergence[target] == v && g.weight[target] == w,
          "pass 2, plane 3: plane 2's pass-1 values, not plane 4's pass-2 ones");
  }

  group("fill: a counted cell whose terms cancel to 0 keeps it and is a source");
  {
    Case c(4, 1.);
    c.add({1.5, 1.5, 1.5}, {{0.75, 0., 0.}});
    c.add({1.5, 0.5, 1.5}, {{0., 0.75, 0.}});
    DivergenceGrid g;
    c.run(g);
    const std::size_t zero = row(g, 1, 1, 1);
    check(g.nDisplacements[zero] == 2 && g.divergence[zero] == 0. && g.weight[zero] == 2.,
          "the cell keeps divergence 0 and weight 2");

    std::vector<std::uint8_t> source(g.nCells, 0);
    std::vector<double> value(g.nCells, 0.), weight(g.nCells, 0.);
    source[zero] = 1;
    weight[zero] = 2.;
    source[row(g, 2, 1, 1)] = 1;
    value[row(g, 2, 1, 1)] = -0.75;
    weight[row(g, 2, 1, 1)] = 1.;
    source[row(g, 1, 0, 1)] = 1;
    value[row(g, 1, 0, 1)] = 0.75;
    weight[row(g, 1, 0, 1)] = 1.;
    double v = 0., w = 0.;
    window(g, value, weight, source, row(g, 2, 1, 2), 1, v, w);
    check(g.divergence[row(g, 2, 1, 2)] == v && g.weight[row(g, 2, 1, 2)] == w,
          "a neighbour fills from the three counted cells in its 3x3x3 window");
  }

  group("errors: E6, naming the first case in tracer-major order; E4");
  {
    const Vec stay = {0.1, 0., 0.};
    Case outside(4, 1., 2);
    outside.add({0.5, 0.5, 0.5}, {{1., 1., 0.}, stay});
    outside.add({1.5, 0.5, 0.5}, {stay, {2.5, 0., 0.}});
    outside.add({3.5, 0.5, 0.5}, {{0.5, 0., 0.}, stay});
    throws_naming([&] { DivergenceGrid g; outside.run(g); },
                  "the valid displacement of tracer 1 in realization 1 ends outside every cell",
                  "E6, on the last face of the grid, the first in tracer-major order");

    Case two(4, 1.);
    two.add({0.5, 0.5, 0.5}, {{1., 0., 0.}});
    throws_naming([&] { DivergenceGrid g; two.run(g); },
                  "cannot reach 62 empty cells: the last pass filled none",
                  "E4, two counted cells only");
  }

  group("scale: a cell counting 2^20 + 2 reruns at q = 61 - e - 21, one counting 2^20 does not");
  {
    int e = 0;
    std::frexp(0.7, &e);
    for (const std::size_t many : {(std::size_t{1} << 19) + 1, std::size_t{1} << 19}) {
      const bool rerun = many > (std::size_t{1} << 19);
      Case c(4, 1.);
      c.x.assign(many, Vec{0.5, 0.5, 0.5});
      c.d.assign(many, std::vector<Vec>{Vec{0.1, 0., 0.}});
      c.add({2.5, 2.5, 2.5}, {{0.7, 0., 0.}});
      c.add({1.5, 2.5, 0.5}, {{0., -0.7, 0.}});
      c.add({2.5, 0.5, 2.5}, {{0., 0., 0.7}});

      DivergenceGrid g1, g4;
      threads(1);
      const internal::BoxDivergence b = c.run(g1);
      threads(4);
      c.run(g4);

      const std::string what = rerun ? "2^20 + 2: " : "2^20: ";
      const int q = 61 - e - (rerun ? 21 : 20);
      check(b.scaleExponent == q, what + "scaleExponent");
      check(g1.nDisplacements[0] == 2 * many, what + "the count of the crowded cell, counted once");

      const std::size_t start[3] = {row(g1, 2, 2, 2), row(g1, 1, 2, 0), row(g1, 2, 0, 2)};
      const std::size_t end[3] = {row(g1, 3, 2, 2), row(g1, 1, 1, 0), row(g1, 2, 0, 3)};
      const double bound = std::ldexp(1., -(q + 1));
      bool within = true;
      for (std::size_t t = 0; t < 3; ++t)
        within = within && std::fabs(g1.divergence[start[t]] - 0.7) <= bound &&
                 std::fabs(g1.divergence[end[t]] + 0.7) <= bound;
      check(within, what + "the crossed cells are within 2^-(q+1) of +-0.7");
      check(g1.divergence == g4.divergence && g1.weight == g4.weight &&
            g1.nDisplacements == g4.nDisplacements,
            what + "the same with 1 and 4 threads");
    }
  }

  group("the same divergence with 1, 3 and 4 threads, and the recorded hashes");
  {
    const Generated in;
    for (const double cellSizeMps : {1., 0.5}) {
      const std::string what = cellSizeMps == 1. ? "cell 1 mps: " : "cell 0.5 mps: ";
      unsigned passes = 0;
      threads(1);
      const std::uint64_t one = in.hash(cellSizeMps, passes);
      threads(3);
      const std::uint64_t three = in.hash(cellSizeMps, passes);
      threads(4);
      const std::uint64_t four = in.hash(cellSizeMps, passes);
      check(one == three && one == four, what + "1, 3 and 4 threads give the same hash");

      const std::uint64_t recorded = cellSizeMps == 1. ? kRecordedHash1 : kRecordedHash05;
      if (one != recorded)
        std::cerr << what << "hash 0x" << std::hex << one << std::dec << ", " << passes
                  << " fill passes" << std::endl;
      check(one == recorded, what + "the hash equals the recorded one");
    }
  }

  return report("test_divergence");
}
