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
 *  @file tests/test_watershed.cpp
 *
 *  @brief The watershed of a box search: the flow and its tie rule, seeds
 *  and their sign, the faces cut, the numbering, the counts, the error
 *  when there is no void, walks longer than the path buffer, agreement
 *  with a serial reference walk, independence of the number of
 *  threads, and watershedBox.
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
#include "fixtures.h"

using namespace bitvf;

namespace {

  // The hashes of the watershed of the generated input, cells of 1 and 0.5
  // mps, sigma 1 mps, recorded once. A different value on some platform
  // means results are not portable there.
  constexpr std::uint64_t kRecordedHash1 = 0x4a0e42d5a02d1affull;
  constexpr std::uint64_t kRecordedHash05 = 0x055f50cab4240dbaull;

  void threads (const int n)
  {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
  }

  /// A box grid of nCellsAxis cells, divergence v everywhere, no tracers.
  DivergenceGrid box (const std::array<std::size_t, 3>& n, const double v)
  {
    DivergenceGrid g;
    g.geometry = otswap::Geometry::Box;
    g.mps = 1.;
    g.cellSize = 1.;
    g.lower = {0., 0., 0.};
    g.upper = {(double)n[0], (double)n[1], (double)n[2]};
    g.nCellsAxis = n;
    g.nCells = n[0] * n[1] * n[2];
    g.divergence.assign(g.nCells, v);
    g.nTracers.assign(g.nCells, 0);
    return g;
  }

  DivergenceGrid cube (const std::size_t n, const double v)
  {
    return box({n, n, n}, v);
  }

  std::int64_t row (const std::size_t n, const std::size_t i, const std::size_t j,
                    const std::size_t k)
  {
    return (std::int64_t)((i * n + j) * n + k);
  }

  /// The watershed of a copy of g.
  struct Run {
    DivergenceGrid grid;
    internal::BoxBasins basins;
  };

  Run run (const DivergenceGrid& g)
  {
    Run r;
    r.grid = g;
    internal::box_watershed(r.grid, r.basins);
    return r;
  }

  /// The serial reference: CBL's flow loop and CBL's walk, which labels
  /// only the cell it starts from, with the seeds numbered by ascending row.
  struct Reference {
    std::vector<std::int64_t> flow, voidIndex, overdenseIndex;
    internal::BoxBasins basins;
  };

  Reference reference (const DivergenceGrid& g)
  {
    const std::int64_t nx = (std::int64_t)g.nCellsAxis[0];
    const std::int64_t ny = (std::int64_t)g.nCellsAxis[1];
    const std::int64_t nz = (std::int64_t)g.nCellsAxis[2];
    const std::int64_t C = (std::int64_t)g.nCells;
    const std::vector<double>& D = g.divergence;
    Reference r;
    r.flow.resize(C);
    for (std::int64_t i = 0; i < nx; ++i)
      for (std::int64_t j = 0; j < ny; ++j)
        for (std::int64_t k = 0; k < nz; ++k) {
          const std::int64_t c = (i * ny + j) * nz + k;
          double minVal = D[c];
          std::int64_t target = c;
          for (std::int64_t ii = std::max<std::int64_t>(0, i - 1); ii <= std::min(nx - 1, i + 1); ++ii)
            for (std::int64_t jj = std::max<std::int64_t>(0, j - 1); jj <= std::min(ny - 1, j + 1); ++jj)
              for (std::int64_t kk = std::max<std::int64_t>(0, k - 1); kk <= std::min(nz - 1, k + 1); ++kk) {
                if (ii == i && jj == j && kk == k) continue;
                const std::int64_t nb = (ii * ny + jj) * nz + kk;
                if (D[nb] < minVal) {
                  minVal = D[nb];
                  target = nb;
                }
              }
          r.flow[c] = target;
        }

    std::vector<std::int64_t> label(C, -1);
    for (std::int64_t c = 0; c < C; ++c)
      if (r.flow[c] == c) {
        internal::BasinSet& set = D[c] < 0. ? r.basins.voids : r.basins.overdense;
        label[c] = D[c] < 0. ? (std::int64_t)set.seed.size() : -2 - (std::int64_t)set.seed.size();
        set.seed.push_back(c);
      }
    for (std::int64_t c = 0; c < C; ++c) {
      if (label[c] != -1) continue;
      std::int64_t cur = c;
      while (label[cur] == -1) cur = r.flow[cur];
      label[c] = label[cur];
    }

    for (internal::BasinSet* set : {&r.basins.voids, &r.basins.overdense}) {
      set->nCells.assign(set->seed.size(), 0);
      set->nCoreCells.assign(set->seed.size(), 0);
      set->nCoreTracers.assign(set->seed.size(), 0);
    }
    r.voidIndex.assign(C, -1);
    r.overdenseIndex.assign(C, -1);
    for (std::int64_t c = 0; c < C; ++c) {
      const bool isVoid = label[c] >= 0;
      const std::int64_t b = isVoid ? label[c] : -2 - label[c];
      internal::BasinSet& set = isVoid ? r.basins.voids : r.basins.overdense;
      (isVoid ? r.voidIndex : r.overdenseIndex)[c] = b;
      ++set.nCells[b];
      if ((D[c] < 0.) == isVoid) {
        ++set.nCoreCells[b];
        set.nCoreTracers[b] += g.nTracers[c];
      }
    }
    return r;
  }

  bool same_set (const internal::BasinSet& a, const internal::BasinSet& b)
  {
    return a.seed == b.seed && a.nCells == b.nCells && a.nCoreCells == b.nCoreCells &&
           a.nCoreTracers == b.nCoreTracers;
  }

  /// Every output of r equal to the reference.
  bool matches (const Run& r, const Reference& ref)
  {
    return r.grid.voidIndex == ref.voidIndex && r.grid.overdenseIndex == ref.overdenseIndex &&
           same_set(r.basins.voids, ref.basins.voids) &&
           same_set(r.basins.overdense, ref.basins.overdense);
  }

  /// The header's guarantees on the output alone: seeds ascending, of the
  /// right sign and labelled with their position; every cell in exactly one
  /// basin; the counts equal to a direct count over the cells; the cells
  /// of an overdense basin all in its core.
  bool consistent (const Run& r)
  {
    const DivergenceGrid& g = r.grid;
    const internal::BasinSet& v = r.basins.voids;
    const internal::BasinSet& o = r.basins.overdense;
    bool ok = g.voidIndex.size() == g.nCells && g.overdenseIndex.size() == g.nCells;
    for (std::size_t b = 0; ok && b < v.seed.size(); ++b)
      ok = (b == 0 || v.seed[b - 1] < v.seed[b]) && g.divergence[v.seed[b]] < 0. &&
           g.voidIndex[v.seed[b]] == (std::int64_t)b;
    for (std::size_t b = 0; ok && b < o.seed.size(); ++b)
      ok = (b == 0 || o.seed[b - 1] < o.seed[b]) && g.divergence[o.seed[b]] >= 0. &&
           g.overdenseIndex[o.seed[b]] == (std::int64_t)b;
    if (!ok) return false;

    internal::BasinSet count;
    std::vector<std::uint64_t> oCells(o.seed.size(), 0), oCore(o.seed.size(), 0), oTracers(o.seed.size(), 0);
    count.nCells.assign(v.seed.size(), 0);
    count.nCoreCells.assign(v.seed.size(), 0);
    count.nCoreTracers.assign(v.seed.size(), 0);
    std::uint64_t total = 0;
    for (std::size_t c = 0; c < g.nCells; ++c) {
      const std::int64_t a = g.voidIndex[c], b = g.overdenseIndex[c];
      const bool inVoid = a >= 0 && a < (std::int64_t)v.seed.size() && b == -1;
      const bool inOverdense = b >= 0 && b < (std::int64_t)o.seed.size() && a == -1;
      if (inVoid == inOverdense) return false;
      if (inVoid) {
        ++count.nCells[a];
        if (g.divergence[c] < 0.) {
          ++count.nCoreCells[a];
          count.nCoreTracers[a] += g.nTracers[c];
        }
      }
      else {
        ++oCells[b];
        if (g.divergence[c] >= 0.) {
          ++oCore[b];
          oTracers[b] += g.nTracers[c];
        }
      }
    }
    for (const std::uint64_t x : v.nCells) total += x;
    for (const std::uint64_t x : o.nCells) total += x;
    return total == g.nCells && count.nCells == v.nCells && count.nCoreCells == v.nCoreCells &&
           count.nCoreTracers == v.nCoreTracers && oCells == o.nCells && oCore == o.nCoreCells &&
           oTracers == o.nCoreTracers && o.nCoreCells == o.nCells;
  }

  /// The output equal to the reference with 1, 3 and 4 threads.
  bool matches_threads (const DivergenceGrid& g)
  {
    const Reference ref = reference(g);
    bool ok = true;
    for (const int t : {1, 3, 4}) {
      threads(t);
      const Run r = run(g);
      ok = ok && matches(r, ref) && consistent(r);
    }
    threads(4);
    return ok;
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

    /// The hash of the watershed on cells of cellSizeMps, sigma 1 mps:
    /// voidIndex, overdenseIndex, then each set's seeds and counts; the
    /// output checked against the reference in agrees.
    std::uint64_t hash (const double cellSizeMps, bool& agrees) const
    {
      const internal::BoxSearch s = internal::check_box(tracers, displacement, valid, 4., config,
                                                        limits, ids, nullptr);
      DivergenceGrid g = internal::box_grid(s, cellSizeMps);
      internal::box_divergence(s, g);
      internal::box_smooth(g, s.sigmaLength, nullptr);
      const Run r = run(g);
      agrees = matches(r, reference(g)) && consistent(r);

      std::uint64_t h = 0xcbf29ce484222325ull;
      hash_vector(h, r.grid.voidIndex);
      hash_vector(h, r.grid.overdenseIndex);
      for (const internal::BasinSet* set : {&r.basins.voids, &r.basins.overdense}) {
        hash_vector(h, set->seed);
        hash_vector(h, set->nCells);
        hash_vector(h, set->nCoreCells);
        hash_vector(h, set->nCoreTracers);
      }
      return h;
    }
  };

  /// A field of n = 40 whose cells, in the even planes i, form one corridor
  /// of width 1 running along k, turning at the ends of each even row j and
  /// stepping to the next even plane at the end of each plane, between walls
  /// of width 1 at divergence 1; the divergence decreases along the
  /// corridor. The corridor's first cell is row 0.
  DivergenceGrid serpentine ()
  {
    const std::size_t n = 40;
    DivergenceGrid g = cube(n, 1.);
    std::vector<std::int64_t> path;
    std::size_t k = 0;
    bool up = true;
    for (std::size_t p = 0; p < n / 2; ++p) {
      const std::size_t i = 2 * p;
      for (std::size_t q = 0; q < n / 2; ++q) {
        const std::size_t j = p % 2 == 0 ? 2 * q : n - 2 - 2 * q;
        if (q > 0) path.push_back(row(n, i, p % 2 == 0 ? j - 1 : j + 1, k));
        else if (p > 0) path.push_back(row(n, i - 1, j, k));
        for (std::size_t s = 0; s < n; ++s) path.push_back(row(n, i, j, up ? s : n - 1 - s));
        k = up ? n - 1 : 0;
        up = !up;
      }
    }
    for (std::size_t q = 0; q < path.size(); ++q) g.divergence[path[q]] = -(double)(q + 1);
    return g;
  }

}

int main ()
{
  threads(4);
  std::mt19937 rng(2026);
  auto uniform = [&rng] (const double lo, const double hi) {
    return lo + (hi - lo) * (rng() / 4294967296.);
  };

  group("strict descent: a cell flows to its strictly smallest neighbour; of two equal smallest, to "
        "the lower row");
  {
    const std::size_t n = 5;
    const std::int64_t c = row(n, 2, 2, 2), a = row(n, 1, 1, 1), b = row(n, 3, 3, 3);
    const double below = std::nextafter(-3., -4.);
    const struct { double a, b; std::int64_t to; const char* what; } cases[] = {
      {-2., -3., 1, "a later neighbour strictly smaller"},
      {-3., -3., 0, "two equal smallest: the lower row"},
      {-3., below, 1, "smaller by one ulp: the smaller"},
    };
    for (const auto& t : cases) {
      DivergenceGrid g = cube(n, 5.);
      g.divergence[c] = 1.;
      g.divergence[a] = t.a;
      g.divergence[b] = t.b;
      const Run r = run(g);
      check(r.basins.voids.seed == std::vector<std::int64_t>{a, b} && r.grid.voidIndex[c] == t.to &&
            consistent(r), t.what);
    }
  }

  group("plateau: two adjacent cells of equal value with no smaller neighbour are both seeds");
  {
    const std::size_t n = 4;
    DivergenceGrid g = cube(n, 1.);
    const std::int64_t a = row(n, 1, 1, 1), b = row(n, 1, 1, 2);
    g.divergence[a] = g.divergence[b] = -1.;
    const Run r = run(g);
    check(r.basins.voids.seed == std::vector<std::int64_t>{a, b} && r.grid.voidIndex[a] == 0 &&
          r.grid.voidIndex[b] == 1 && consistent(r), "two voids");
  }

  group("a seed at 0 or -0 starts an overdense basin, at the smallest negative subnormal a void");
  for (const double v : {0., -0., -std::numeric_limits<double>::denorm_min()}) {
    const std::size_t n = 5;
    DivergenceGrid g = cube(n, 1.);
    const std::int64_t c = row(n, 2, 2, 2);
    g.divergence[0] = -1.;
    g.divergence[c] = v;
    const Run r = run(g);
    const bool isVoid = v < 0.;
    const std::string what = v != 0. ? "-denorm_min" : std::signbit(v) ? "-0" : "0";
    if (isVoid)
      check(r.basins.voids.seed == std::vector<std::int64_t>{0, c} && r.grid.voidIndex[c] == 1 &&
            r.grid.overdenseIndex[c] == -1 && consistent(r), what);
    else {
      const std::vector<std::int64_t>& seed = r.basins.overdense.seed;
      bool found = false;
      for (std::size_t b = 0; b < seed.size(); ++b)
        found = found || (seed[b] == c && r.grid.overdenseIndex[c] == (std::int64_t)b);
      check(r.basins.voids.seed == std::vector<std::int64_t>{0} && found &&
            r.grid.voidIndex[c] == -1 && consistent(r), what);
    }
  }

  group("faces, edges and corners: neighbours cut, no periodicity");
  {
    const std::size_t n = 5;
    const struct { std::int64_t a, image; const char* what; } cases[] = {
      {row(n, 0, 2, 2), row(n, 4, 2, 2), "face"},
      {row(n, 0, 0, 2), row(n, 4, 4, 2), "edge"},
      {row(n, 0, 0, 0), row(n, 4, 4, 4), "corner"},
    };
    for (const auto& t : cases) {
      DivergenceGrid g = cube(n, 1.);
      g.divergence[t.a] = -1.;
      g.divergence[t.image] = -2.;
      const Run r = run(g);
      check(r.basins.voids.seed == std::vector<std::int64_t>{t.a, t.image} && consistent(r),
            std::string(t.what) + ": a seed whose periodic image is lower stays a seed");
    }
    DivergenceGrid g = cube(n, 1.);
    g.divergence[row(n, 0, 2, 2)] = -0.5;
    g.divergence[row(n, 1, 2, 2)] = -1.;
    g.divergence[row(n, 4, 2, 2)] = -2.;
    const Run r = run(g);
    check(r.grid.voidIndex[row(n, 0, 2, 2)] == r.grid.voidIndex[row(n, 1, 2, 2)] &&
          r.basins.voids.seed.size() == 2 && consistent(r),
          "a face cell flows inside the grid, not to its periodic image");
  }

  group("no void: a divergence of 0 or above in every cell gives an empty void set, voidIndex -1 in every "
        "cell, no error");
  {
    const auto noVoid = [](const DivergenceGrid& g, const std::string& what) {
      const Run r = run(g);
      bool none = true;
      for (const std::int64_t b : r.grid.voidIndex) none = none && b == -1;
      check(r.basins.voids.seed.empty() && r.basins.voids.nCells.empty() && !r.basins.overdense.seed.empty() &&
            none && consistent(r), what);
    };
    for (const double v : {1., 0., -0.}) noVoid(cube(6, v), "constant " + std::to_string(v));
    DivergenceGrid g = cube(7, 0.);
    const double levels[] = {0., -0., 1., 2.};
    for (double& d : g.divergence) d = levels[rng() % 4];
    noVoid(g, "levels 0, -0, 1, 2");
    noVoid(cube(1, 0.), "one cell at 0");
  }

  group("no overdense seed: an empty overdense set, overdenseIndex -1 in every cell");
  {
    DivergenceGrid g = cube(7, 0.);
    for (double& d : g.divergence) d = uniform(-2., -1.);
    const Run r = run(g);
    bool none = true;
    for (const std::int64_t b : r.grid.overdenseIndex) none = none && b == -1;
    check(r.basins.overdense.seed.empty() && r.basins.overdense.nCells.empty() && none &&
          consistent(r), "negative field");
  }

  group("random fields, continuous and with many exact ties: equal to the serial reference walk "
        "with 1, 3 and 4 threads; numbering, partition and counts");
  for (const std::size_t n : {1, 2, 3, 7, 16, 50})
    for (const bool ties : {false, true}) {
      DivergenceGrid g = cube(n, 0.);
      const double levels[] = {-2., -1., -0., 0., 1., 2.};
      for (std::size_t c = 0; c < g.nCells; ++c) {
        g.divergence[c] = ties ? levels[rng() % 6] : uniform(-1., 1.);
        g.nTracers[c] = rng() % 5;
      }
      if (n == 1) g.divergence[0] = -1.;
      check(matches_threads(g), "n " + std::to_string(n) + (ties ? ", ties" : ", continuous"));
    }

  group("any cell counts: non-cubic grids, axes of 1 cell included, equal to the serial reference walk "
        "with 1, 3 and 4 threads");
  for (const std::array<std::size_t, 3> n : {std::array<std::size_t, 3>{5, 9, 13}, {13, 5, 9}, {9, 13, 5},
                                             {1, 7, 11}, {7, 1, 11}, {7, 11, 1}, {1, 1, 40}, {2, 3, 4}})
    for (const bool ties : {false, true}) {
      DivergenceGrid g = box(n, 0.);
      const double levels[] = {-2., -1., -0., 0., 1., 2.};
      for (std::size_t c = 0; c < g.nCells; ++c) {
        g.divergence[c] = ties ? levels[rng() % 6] : uniform(-1., 1.);
        g.nTracers[c] = rng() % 5;
      }
      check(matches_threads(g), "(" + std::to_string(n[0]) + ", " + std::to_string(n[1]) + ", " +
            std::to_string(n[2]) + ")" + (ties ? ", ties" : ", continuous"));
    }

  group("any cell counts: a line of cells, along each axis, is the 1-D descent (left neighbour first)");
  for (const std::size_t axis : {0, 1, 2}) {
    const std::size_t L = 60;
    std::array<std::size_t, 3> n = {1, 1, 1};
    n[axis] = L;
    DivergenceGrid g = box(n, 0.);
    const double levels[] = {-2., -1., 0., 1., 2.};
    for (double& d : g.divergence) d = rng() % 2 ? levels[rng() % 5] : uniform(-2., 2.);
    const std::vector<double>& D = g.divergence;
    std::vector<std::size_t> flow(L);
    for (std::size_t c = 0; c < L; ++c) {
      flow[c] = c;
      if (c > 0 && D[c - 1] < D[c]) flow[c] = c - 1;
      if (c + 1 < L && D[c + 1] < D[flow[c]]) flow[c] = c + 1;
    }
    std::vector<std::int64_t> voidIndex(L, -1), overdenseIndex(L, -1);
    std::int64_t nVoid = 0, nOverdense = 0;
    for (std::size_t c = 0; c < L; ++c)
      if (flow[c] == c) (D[c] < 0. ? voidIndex[c] = nVoid++ : overdenseIndex[c] = nOverdense++);
    for (std::size_t c = 0; c < L; ++c) {
      std::size_t cur = c;
      while (flow[cur] != cur) cur = flow[cur];
      voidIndex[c] = voidIndex[cur];
      overdenseIndex[c] = overdenseIndex[cur];
    }
    const Run r = run(g);
    check(r.grid.voidIndex == voidIndex && r.grid.overdenseIndex == overdenseIndex && consistent(r),
          "axis " + std::to_string(axis));
  }

  group("no nTracers: the same partition and cell counts, nCoreTracers empty");
  {
    DivergenceGrid g = box({6, 4, 9}, 0.);
    for (std::size_t c = 0; c < g.nCells; ++c) {
      g.divergence[c] = uniform(-1., 1.);
      g.nTracers[c] = rng() % 5;
    }
    const Run with = run(g);
    g.nTracers.clear();
    const Run without = run(g);
    bool ok = without.grid.voidIndex == with.grid.voidIndex &&
              without.grid.overdenseIndex == with.grid.overdenseIndex;
    const internal::BasinSet* x[2] = {&with.basins.voids, &with.basins.overdense};
    const internal::BasinSet* y[2] = {&without.basins.voids, &without.basins.overdense};
    for (int kind = 0; kind < 2; ++kind)
      ok = ok && x[kind]->seed == y[kind]->seed && x[kind]->nCells == y[kind]->nCells &&
           x[kind]->nCoreCells == y[kind]->nCoreCells && y[kind]->nCoreTracers.empty();
    check(ok, "(6, 4, 9)");
  }

  group("a walk longer than kPathCells: equal to the serial reference walk with 1, 3 and 4 threads");
  {
    const DivergenceGrid g = serpentine();
    const Reference ref = reference(g);
    std::size_t steps = 0;
    for (std::int64_t c = 0; ref.flow[c] != c; c = ref.flow[c]) ++steps;
    std::cout << "    descent from row 0: " << steps << " steps" << std::endl;
    check(steps > internal::kPathCells, "the descent from row 0 is longer than kPathCells");
    check(ref.basins.voids.seed.size() == 1 && ref.basins.overdense.seed.empty(),
          "one void, no overdense basin");
    check(matches_threads(g), "labels and counts");
  }

  group("the same watershed with 1, 3 and 4 threads, equal to the reference, and the recorded hashes");
  {
    const Generated in;
    for (const double cellSizeMps : {1., 0.5}) {
      const std::string what = cellSizeMps == 1. ? "cell 1 mps: " : "cell 0.5 mps: ";
      bool agrees1 = false, agrees3 = false, agrees4 = false;
      threads(1);
      const std::uint64_t one = in.hash(cellSizeMps, agrees1);
      threads(3);
      const std::uint64_t three = in.hash(cellSizeMps, agrees3);
      threads(4);
      const std::uint64_t four = in.hash(cellSizeMps, agrees4);
      check(agrees1 && agrees3 && agrees4, what + "equal to the reference, consistent");
      check(one == three && one == four, what + "1, 3 and 4 threads give the same hash");

      const std::uint64_t recorded = cellSizeMps == 1. ? kRecordedHash1 : kRecordedHash05;
      if (one != recorded)
        std::cerr << what << "hash 0x" << std::hex << one << std::dec << std::endl;
      check(one == recorded, what + "the hash equals the recorded one");
    }
  }

  group("watershedBox: equal to the serial reference walk with 1 and 4 threads, on cubic and non-cubic grids, "
        "axes of 1 cell included; seeds, labels and cell counts");
  for (const std::array<std::size_t, 3> n : {std::array<std::size_t, 3>{16, 16, 16}, {5, 9, 13}, {2, 3, 4},
                                             {1, 7, 11}, {7, 11, 1}, {1, 1, 40}, {1, 1, 1}})
    for (const bool ties : {false, true}) {
      DivergenceGrid g = box(n, 0.);
      const double levels[] = {-2., -1., -0., 0., 1., 2.};
      for (std::size_t c = 0; c < g.nCells; ++c) g.divergence[c] = ties ? levels[rng() % 6] : uniform(-1., 1.);
      const Reference ref = reference(g);
      bool ok = true;
      for (const int t : {1, 4}) {
        threads(t);
        const Partition p = watershedBox(g.divergence, n);
        ok = ok && p.nVoids == ref.basins.voids.seed.size() && p.nOverdense == ref.basins.overdense.seed.size() &&
             p.voidIndex == ref.voidIndex && p.overdenseIndex == ref.overdenseIndex &&
             p.voidSeed == ref.basins.voids.seed && p.overdenseSeed == ref.basins.overdense.seed &&
             p.voidCells == ref.basins.voids.nCells && p.voidCoreCells == ref.basins.voids.nCoreCells &&
             p.overdenseCells == ref.basins.overdense.nCells &&
             p.overdenseCoreCells == ref.basins.overdense.nCoreCells;
      }
      threads(4);
      check(ok, "(" + std::to_string(n[0]) + ", " + std::to_string(n[1]) + ", " + std::to_string(n[2]) + ")" +
            (ties ? ", ties" : ", continuous"));
    }

  group("watershedBox on a search's smoothed divergence: the search's labels, seeds and basin cell counts, "
        "both kinds of basin");
  {
    const fixtures::Translated in;
    const Result r = in.run();
    const DivergenceGrid& d = *r.divergence;
    const Partition p = watershedBox(d.divergence, d.nCellsAxis);
    check(r.overdense.nBasins > 0 && p.nVoids == r.voids.nBasins && p.nOverdense == r.overdense.nBasins,
          "basin counts");
    check(p.voidIndex == d.voidIndex && p.overdenseIndex == d.overdenseIndex, "labels");
    std::vector<std::uint64_t> voidCells(p.nVoids, 0), voidCore(p.nVoids, 0);
    std::vector<std::uint64_t> overdenseCells(p.nOverdense, 0);
    for (std::size_t c = 0; c < d.nCells; ++c) {
      if (d.voidIndex[c] >= 0) {
        ++voidCells[(std::size_t)d.voidIndex[c]];
        if (d.divergence[c] < 0.) ++voidCore[(std::size_t)d.voidIndex[c]];
      }
      else ++overdenseCells[(std::size_t)d.overdenseIndex[c]];
    }
    bool seeds = true;
    for (std::size_t b = 0; b < p.nVoids; ++b)
      seeds = seeds && d.voidIndex[(std::size_t)p.voidSeed[b]] == (std::int64_t)b &&
              (b == 0 || p.voidSeed[b - 1] < p.voidSeed[b]);
    for (std::size_t b = 0; b < p.nOverdense; ++b)
      seeds = seeds && d.overdenseIndex[(std::size_t)p.overdenseSeed[b]] == (std::int64_t)b &&
              d.divergence[(std::size_t)p.overdenseSeed[b]] >= 0.;
    check(seeds, "seeds ascending, in their basins, overdense seeds at or above 0");
    check(p.voidCells == voidCells && p.voidCoreCells == voidCore && p.overdenseCells == overdenseCells &&
          p.overdenseCoreCells == overdenseCells, "cell counts");
  }

  group("watershedBox: no cell below 0 gives no void and no error");
  {
    std::vector<double> field(2 * 3 * 4);
    for (double& v : field) v = uniform(0., 1.);
    field[5] = -0.;
    const Partition p = watershedBox(field, {2, 3, 4});
    check(p.nVoids == 0 && p.voidSeed.empty() && p.voidCells.empty() && p.nOverdense > 0 &&
          p.voidIndex == std::vector<std::int64_t>(field.size(), -1), "no void, every cell overdense");
  }

  group("watershedBox errors WS1-WS3");
  {
    const std::vector<double> field(2 * 3 * 4, -1.);
    throws_naming([&] { watershedBox(field, {2, 0, 4}); }, "nCellsAxis is (2, 0, 4); every count must be at least 1",
                  "WS1");
    throws_naming([&] { watershedBox(field, {2, 3, 5}); }, "field holds 24 entries for (2, 3, 5) cells", "WS2");
    throws_naming([&] { watershedBox(field, {(std::size_t)1 << 40, (std::size_t)1 << 40, 1}); }, "field holds 24",
                  "WS2, a product that overflows");
    std::vector<double> bad = field;
    bad[9] = std::numeric_limits<double>::infinity();
    throws_naming([&] { watershedBox(bad, {2, 3, 4}); }, "field holds a non-finite value at entry 9", "WS3");
  }

  return report("test_watershed");
}
