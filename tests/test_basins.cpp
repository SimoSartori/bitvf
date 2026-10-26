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
 *  @file tests/test_basins.cpp
 *
 *  @brief The catalogues of a box search: centres and their clamp to the
 *  seed cell, radii, deltaAv, nTracers, delta2Mps and its sphere, edge
 *  flags and contamination, the assignment of the tracers, and
 *  independence of the number of threads.
 */

#include <algorithm>
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

#include "otswap/Maths.h"

#include "internal.h"
#include "check.h"

using namespace bitvf;

namespace {

  // The hashes of the catalogues and the assignment of the generated
  // input, cells of 1 and 0.5 mps, sigma 1 mps, recorded once. A different
  // value on some platform means results are not portable there.
  constexpr std::uint64_t kRecordedHash1 = 0x635ad8bd9ce6e272ull;
  constexpr std::uint64_t kRecordedHash05 = 0x6896aeaadc5576faull;

  constexpr double kPi = 3.141592653589793;
  constexpr double kInf = std::numeric_limits<double>::infinity();

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

  /// A box grid of n cells per axis of side s from lower, divergence v
  /// everywhere, no tracers, mps s.
  DivergenceGrid cube (const std::size_t n, const double v, const double lower = 0., const double s = 1.)
  {
    DivergenceGrid g;
    g.geometry = otswap::Geometry::Box;
    g.mps = s;
    g.cellSize = s;
    g.lower = {lower, lower, lower};
    g.upper = {lower + s * n, lower + s * n, lower + s * n};
    g.nCellsAxis = {n, n, n};
    g.nCells = n * n * n;
    g.divergence.assign(g.nCells, v);
    g.nTracers.assign(g.nCells, 0);
    return g;
  }

  std::int64_t row (const std::size_t n, const std::size_t i, const std::size_t j, const std::size_t k)
  {
    return (std::int64_t)((i * n + j) * n + k);
  }

  /// A grid split into basins by hand: label >= 0 is that void, label < 0
  /// the overdense basin -1 - label. The seed of a basin is its first core
  /// cell in row order, else its first cell, unless given in seeds.
  struct Hand {
    DivergenceGrid grid;
    internal::BoxBasins basins;
  };

  Hand partition (const DivergenceGrid& g, const std::vector<std::int64_t>& label,
                  const std::vector<std::int64_t>& seeds = {})
  {
    Hand h;
    h.grid = g;
    std::int64_t nVoid = 0, nOver = 0;
    for (const std::int64_t l : label) {
      if (l >= 0) nVoid = std::max(nVoid, l + 1);
      else nOver = std::max(nOver, -l);
    }
    internal::BasinSet* sets[2] = {&h.basins.voids, &h.basins.overdense};
    for (int kind = 0; kind < 2; ++kind) {
      const std::size_t B = (std::size_t)(kind == 0 ? nVoid : nOver);
      sets[kind]->seed.assign(B, -1);
      sets[kind]->nCells.assign(B, 0);
      sets[kind]->nCoreCells.assign(B, 0);
      sets[kind]->nCoreTracers.assign(B, 0);
    }
    h.grid.voidIndex.assign(g.nCells, -1);
    h.grid.overdenseIndex.assign(g.nCells, -1);
    std::vector<std::int64_t> first[2] = {std::vector<std::int64_t>(nVoid, -1),
                                          std::vector<std::int64_t>(nOver, -1)};
    for (std::size_t c = 0; c < g.nCells; ++c) {
      const bool isVoid = label[c] >= 0;
      const std::size_t b = (std::size_t)(isVoid ? label[c] : -1 - label[c]);
      internal::BasinSet& set = isVoid ? h.basins.voids : h.basins.overdense;
      (isVoid ? h.grid.voidIndex : h.grid.overdenseIndex)[c] = (std::int64_t)b;
      ++set.nCells[b];
      if (first[isVoid ? 0 : 1][b] < 0) first[isVoid ? 0 : 1][b] = (std::int64_t)c;
      if (isVoid ? g.divergence[c] < 0. : g.divergence[c] >= 0.) {
        ++set.nCoreCells[b];
        set.nCoreTracers[b] += g.nTracers[c];
        if (set.seed[b] < 0) set.seed[b] = (std::int64_t)c;
      }
    }
    for (int kind = 0; kind < 2; ++kind)
      for (std::size_t b = 0; b < sets[kind]->seed.size(); ++b)
        if (sets[kind]->seed[b] < 0) sets[kind]->seed[b] = first[kind][b];
    for (const std::int64_t c : seeds) {
      const bool isVoid = h.grid.voidIndex[c] >= 0;
      (isVoid ? h.basins.voids : h.basins.overdense)
        .seed[isVoid ? h.grid.voidIndex[c] : h.grid.overdenseIndex[c]] = c;
    }
    return h;
  }

  /// The whole grid as one basin with the given seed: a void when the seed
  /// is below 0, else an overdense basin.
  Hand one_basin (const DivergenceGrid& g, const std::int64_t seed)
  {
    return partition(g, std::vector<std::int64_t>(g.nCells, g.divergence[seed] < 0. ? 0 : -1), {seed});
  }

  struct Catalogues {
    BasinCatalog voids, overdense;
  };

  /// box_basins on a grid split by hand, with the given tracers; cell
  /// gives which are read (cell -1 not read), box_cell of each by default.
  Catalogues catalogues (const Hand& h, const std::vector<double>& tracers = {},
                         std::vector<std::int64_t> cell = {})
  {
    internal::BoxSearch search;
    search.nTracers = tracers.size() / 3;
    search.nRealizations = 1;
    search.tracers = &tracers;
    search.mps = h.grid.mps;
    internal::BoxDivergence div;
    if (cell.empty())
      for (std::size_t i = 0; i < search.nTracers; ++i) cell.push_back(internal::box_cell(h.grid, &tracers[3*i]));
    div.cell = cell;
    Catalogues out;
    internal::box_basins(search, div, h.grid, h.basins, out.voids, out.overdense);
    return out;
  }

  /// The centre CBL computes before any clamp, evaluated directly.
  std::array<double, 3> unclamped (const DivergenceGrid& g, const std::int64_t c)
  {
    const std::int64_t n = (std::int64_t)g.nCellsAxis[0];
    const std::int64_t i[3] = {c / (n * n), (c / n) % n, c % n};
    const double s = g.cellSize;
    double num[3] = {0., 0., 0.}, den = 0.;
    for (std::int64_t a = -1; a <= 1; ++a)
      for (std::int64_t b = -1; b <= 1; ++b)
        for (std::int64_t d = -1; d <= 1; ++d) {
          const std::int64_t ii = i[0] + a, jj = i[1] + b, kk = i[2] + d;
          if (ii < 0 || jj < 0 || kk < 0 || ii >= n || jj >= n || kk >= n) continue;
          if (a == 0 && b == 0 && d == 0) continue;
          const double D = g.divergence[(ii * n + jj) * n + kk];
          num[0] += D * s * (i[0] - ii);
          num[1] += D * s * (i[1] - jj);
          num[2] += D * s * (i[2] - kk);
          den += std::abs(D);
        }
    if (den == 0.) den = 1.;
    std::array<double, 3> x;
    for (int k = 0; k < 3; ++k) x[k] = g.lower[k] + s * (i[k] + 0.5) + num[k] / den;
    return x;
  }

  /// The edge flag and contamination of every basin recounted from the
  /// definitions: a cell's faces on the faces of the grid are the
  /// coordinates at 0 or n - 1; a face shared with another basin is one
  /// whose neighbour has a different label.
  bool brute_force (const Hand& h, const Catalogues& cat)
  {
    const DivergenceGrid& g = h.grid;
    const std::int64_t n = (std::int64_t)g.nCellsAxis[0];
    const std::size_t nVoid = h.basins.voids.seed.size();
    const std::size_t nAll = nVoid + h.basins.overdense.seed.size();
    std::vector<std::uint64_t> out(nAll, 0), outCore(nAll, 0), phy(nAll, 0);
    auto basin = [&] (const std::int64_t c) {
      return g.voidIndex[c] >= 0 ? g.voidIndex[c] : (std::int64_t)nVoid + g.overdenseIndex[c];
    };
    for (std::int64_t i = 0; i < n; ++i)
      for (std::int64_t j = 0; j < n; ++j)
        for (std::int64_t k = 0; k < n; ++k) {
          const std::int64_t c = (i * n + j) * n + k;
          const std::int64_t b = basin(c);
          const bool core = g.voidIndex[c] >= 0 ? g.divergence[c] < 0. : g.divergence[c] >= 0.;
          std::uint64_t onFace = 0;
          for (const std::int64_t x : {i, j, k}) onFace += (x == 0) + (x == n - 1);
          out[b] += onFace;
          if (core) outCore[b] += onFace;
          if (i > 0 && basin(c - n * n) != b) ++phy[b];
          if (i < n - 1 && basin(c + n * n) != b) ++phy[b];
          if (j > 0 && basin(c - n) != b) ++phy[b];
          if (j < n - 1 && basin(c + n) != b) ++phy[b];
          if (k > 0 && basin(c - 1) != b) ++phy[b];
          if (k < n - 1 && basin(c + 1) != b) ++phy[b];
        }
    bool ok = cat.voids.edgeFlag.size() == nVoid && cat.overdense.edgeFlag.size() == nAll - nVoid;
    for (std::size_t b = 0; ok && b < nAll; ++b) {
      const BasinCatalog& c = b < nVoid ? cat.voids : cat.overdense;
      const std::size_t r = b < nVoid ? b : b - nVoid;
      const std::uint8_t flag = outCore[b] > 0 ? 2 : out[b] > 0 ? 1 : 0;
      const double cont = out[b] + phy[b] > 0 ? (double)out[b] / (double)(out[b] + phy[b]) : 0.;
      ok = c.edgeFlag[r] == flag && c.contamination[r] == cont;
    }
    return ok;
  }

  /// CBL's delta2Mps for n_t tracers in the sphere.
  double delta2 (const std::size_t nt, const double m)
  {
    return ((double)nt / (4 * kPi * 8 * m * m * m / 3.)) / (1. / (m * m * m)) - 1.;
  }

  /// Every centre of the catalogues in its seed cell.
  bool in_seed (const DivergenceGrid& g, const internal::BoxBasins& basins, const Catalogues& cat)
  {
    bool ok = true;
    for (std::size_t b = 0; b < basins.voids.seed.size(); ++b)
      ok = ok && internal::box_cell(g, &cat.voids.centre[3*b]) == basins.voids.seed[b];
    for (std::size_t b = 0; b < basins.overdense.seed.size(); ++b)
      ok = ok && internal::box_cell(g, &cat.overdense.centre[3*b]) == basins.overdense.seed[b];
    return ok;
  }

  /// box_cell of the centre of every cell of g is that cell.
  bool centres_in_cells (DivergenceGrid g)
  {
    internal::box_centres(g);
    bool ok = true;
    for (std::size_t c = 0; ok && c < g.nCells; ++c)
      ok = internal::box_cell(g, &g.centre[3*c]) == (std::int64_t)c;
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
  /// of the tracers and of points uniform in [-8, 108)^3, which covers
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

    /// The hash of the catalogues and the assignment on cells of
    /// cellSizeMps, sigma 1 mps: every column of the voids, then of the
    /// overdense basins, then of the assignment. The header's guarantees
    /// checked on the output in agrees; with countSpheres, delta2Mps also
    /// against a direct count over every tracer read.
    std::uint64_t hash (const double cellSizeMps, const bool countSpheres, bool& agrees) const
    {
      const internal::BoxSearch s = internal::check_box(tracers, displacement, valid, 4., config,
                                                        limits, ids, nullptr);
      Hand h;
      h.grid = internal::box_grid(s, cellSizeMps);
      const internal::BoxDivergence div = internal::box_divergence(s, h.grid);
      internal::box_smooth(h.grid, s.sigmaLength, nullptr);
      internal::box_watershed(h.grid, h.basins);
      Catalogues cat;
      internal::box_basins(s, div, h.grid, h.basins, cat.voids, cat.overdense);
      const TracerAssignment as = internal::box_assignment(s, div, h.grid, cat.voids, cat.overdense);

      const internal::BasinSet& o = h.basins.overdense;
      bool overdense = o.nCoreCells == o.nCells;
      for (std::size_t b = 0; b < o.seed.size(); ++b)
        overdense = overdense && cat.overdense.edgeFlag[b] != 1 &&
                    cat.overdense.radius[b] == cat.overdense.radiusExt[b];
      agrees = in_seed(h.grid, h.basins, cat) && brute_force(h, cat) && overdense;

      if (countSpheres) {
        const double r2 = (2 * s.mps) * (2 * s.mps);
        for (const BasinCatalog* c : {&cat.voids, &cat.overdense})
          for (std::size_t b = 0; agrees && b < c->nBasins; ++b) {
            std::size_t nt = 0;
            for (std::size_t i = 0; i < N; ++i) {
              if (div.cell[i] < 0) continue;
              const double dx = tracers[3*i] - c->centre[3*b];
              const double dy = tracers[3*i+1] - c->centre[3*b+1];
              const double dz = tracers[3*i+2] - c->centre[3*b+2];
              nt += dx*dx + dy*dy + dz*dz <= r2;
            }
            agrees = c->delta2Mps[b] == delta2(nt, s.mps);
          }
      }

      std::uint64_t hash = 0xcbf29ce484222325ull;
      for (const BasinCatalog* c : {&cat.voids, &cat.overdense}) {
        hash_bytes(hash, &c->nBasins, sizeof(c->nBasins));
        hash_vector(hash, c->centre);
        hash_vector(hash, c->radius);
        hash_vector(hash, c->radiusExt);
        hash_vector(hash, c->deltaAv);
        hash_vector(hash, c->delta2Mps);
        hash_vector(hash, c->edgeFlag);
        hash_vector(hash, c->contamination);
        hash_vector(hash, c->nTracers);
      }
      hash_vector(hash, as.cellIndex);
      hash_vector(hash, as.voidIndex);
      hash_vector(hash, as.overdenseIndex);
      hash_vector(hash, as.distance);
      hash_vector(hash, as.inMask);
      hash_vector(hash, as.inBox);
      hash_vector(hash, as.used);
      return hash;
    }
  };

}

int main ()
{
  threads(4);
  std::mt19937 rng(2026);
  auto uniform = [&rng] (const double lo, const double hi) {
    return lo + (hi - lo) * (rng() / 4294967296.);
  };

  group("centre: CBL's expression over the neighbours, bit for bit where no clamp binds; faces, edges "
        "and corners cut; a zero neighbourhood gives the cell centre; voids and overdense seeds alike");
  {
    const std::size_t n = 5;
    const std::size_t places[][3] = {{2, 2, 2}, {0, 2, 2}, {4, 2, 2}, {2, 0, 4}, {0, 0, 2}, {4, 4, 0},
                                     {0, 0, 0}, {4, 4, 4}, {0, 4, 0}};
    std::size_t free = 0, total = 0;
    for (int trial = 0; trial < 200; ++trial)
      for (const auto& p : places)
        for (const bool isVoid : {true, false}) {
          DivergenceGrid g = cube(n, 0.);
          for (double& d : g.divergence) d = uniform(-1., 1.);
          const std::int64_t c = row(n, p[0], p[1], p[2]);
          g.divergence[c] = isVoid ? -2. : 2.;
          const Hand h = one_basin(g, c);
          const Catalogues cat = catalogues(h);
          const double* x = &(isVoid ? cat.voids : cat.overdense).centre[0];
          const std::array<double, 3> u = unclamped(g, c);
          bool inside = true;
          for (int a = 0; a < 3; ++a) {
            const double t = (u[a] - g.lower[a]) * (1. / g.cellSize);
            inside = inside && t >= (double)p[a] && t < (double)p[a] + 1.;
          }
          ++total;
          if (inside) {
            ++free;
            if (!(x[0] == u[0] && x[1] == u[1] && x[2] == u[2]))
              check(false, "trial " + std::to_string(trial) + ": the unclamped centre bit for bit");
          }
          if (internal::box_cell(g, x) != c)
            check(false, "trial " + std::to_string(trial) + ": the centre in its seed cell");
        }
    std::cout << "    " << free << " of " << total << " centres unclamped" << std::endl;
    check(free > 0 && free < total, "both clamped and unclamped centres occur");

    for (const auto& p : places) {
      DivergenceGrid g = cube(n, 0.);
      const std::int64_t c = row(n, p[0], p[1], p[2]);
      g.divergence[c] = -1.;
      const Catalogues cat = catalogues(one_basin(g, c));
      check(cat.voids.centre[0] == p[0] + 0.5 && cat.voids.centre[1] == p[1] + 0.5 &&
            cat.voids.centre[2] == p[2] + 0.5, "all neighbours 0: the cell centre");
    }
  }

  group("clamp: a centre pulled past each of the six faces of its seed cell stays in it, box_cell of it "
        "the seed, at the closest point; faces at which the cell lookup rounds across are stepped");
  {
    const std::size_t n = 5;
    for (int a = 0; a < 3; ++a)
      for (const int side : {-1, 1}) {
        DivergenceGrid g = cube(n, 0.);
        const std::int64_t c = row(n, 2, 2, 2);
        g.divergence[c] = -1.;
        std::size_t q[3] = {2, 2, 2};
        q[a] = (std::size_t)(2 - side);
        g.divergence[row(n, q[0], q[1], q[2])] = 1.;
        const Catalogues cat = catalogues(one_basin(g, c));
        const double* x = &cat.voids.centre[0];
        double expected[3] = {2.5, 2.5, 2.5};
        expected[a] = side < 0 ? 2. : std::nextafter(3., 0.);
        check(x[0] == expected[0] && x[1] == expected[1] && x[2] == expected[2] &&
              internal::box_cell(g, x) == c,
              "axis " + std::to_string(a) + (side < 0 ? ", lower face" : ", upper face"));
      }

    struct Face { double lower, s; std::size_t i; bool upper; };
    std::vector<Face> found;
    bool haveLower = false, haveUpper = false;
    for (int trial = 0; trial < 2000000 && !(haveLower && haveUpper); ++trial) {
      const double lower = uniform(-1000., 1000.);
      const double s = uniform(0.01, 10.);
      const std::size_t i = rng() % 29 + 1;
      const double inverse = 1. / s;
      const double lo = lower + s * i;
      const double hi = std::nextafter(lower + s * (i + 1.), -kInf);
      if (!haveLower && (lo - lower) * inverse < (double)i) {
        found.push_back({lower, s, i, false});
        haveLower = true;
      }
      if (!haveUpper && (hi - lower) * inverse >= (double)i + 1.) {
        found.push_back({lower, s, i, true});
        haveUpper = true;
      }
    }
    check(haveLower, "a lower face whose lookup falls in the cell below was found");
    check(haveUpper, "an upper face one ulp below which the lookup still falls in the cell above was found");

    for (const Face& f : found) {
      const std::size_t m = f.i + 2;
      for (int a = 0; a < 3; ++a) {
        DivergenceGrid g = cube(m, 0., f.lower, f.s);
        const std::int64_t c = row(m, f.i, f.i, f.i);
        g.divergence[c] = -1.;
        std::size_t q[3] = {f.i, f.i, f.i};
        q[a] = f.upper ? f.i - 1 : f.i + 1;
        g.divergence[row(m, q[0], q[1], q[2])] = 1.;
        const Catalogues cat = catalogues(one_basin(g, c));
        const double* x = &cat.voids.centre[0];
        double beyond[3] = {x[0], x[1], x[2]};
        beyond[a] = std::nextafter(x[a], f.upper ? kInf : -kInf);
        const std::string what = std::string(f.upper ? "upper" : "lower") + " face, axis " + std::to_string(a);
        check(internal::box_cell(g, x) == c, what + ": in the seed cell");
        check(internal::box_cell(g, beyond) != c, what + ": one more ulp leaves it");
        check(centres_in_cells(g), what + ": every cell centre in its cell");
      }
    }
    check(centres_in_cells(cube(5, 0.)) && centres_in_cells(cube(7, 0., -0.3, 0.7)),
          "every cell centre in its cell, other grids");
  }

  group("radii, deltaAv and nTracers from the counts of the basins");
  {
    const std::size_t n = 6;
    DivergenceGrid g = cube(n, 0.);
    for (double& d : g.divergence) d = uniform(-1., 1.);
    for (std::uint64_t& t : g.nTracers) t = rng() % 7;
    std::vector<std::int64_t> label(g.nCells);
    for (std::size_t c = 0; c < g.nCells; ++c) label[c] = (std::int64_t)(c % 5) - 2;
    for (std::size_t c = 0; c < 5; ++c) g.divergence[c] = label[c] >= 0 ? -1. : 1.;
    g.mps = 1.7;
    g.cellSize = 0.9;
    const Hand h = partition(g, label);
    const Catalogues cat = catalogues(h);
    const double s = g.cellSize, m = g.mps;
    bool ok = true, ulp = true;
    for (int kind = 0; kind < 2; ++kind) {
      const internal::BasinSet& set = kind == 0 ? h.basins.voids : h.basins.overdense;
      const BasinCatalog& c = kind == 0 ? cat.voids : cat.overdense;
      ok = ok && c.nBasins == set.seed.size();
      for (std::size_t b = 0; ok && b < set.seed.size(); ++b) {
        const double r = s * otswap::maths::cbrt(3.0 * (double)set.nCoreCells[b] / (4.0 * kPi));
        const double re = s * otswap::maths::cbrt(3.0 * (double)set.nCells[b] / (4.0 * kPi));
        const double dAv = ((double)set.nCoreTracers[b] / ((double)set.nCoreCells[b] * s * s * s)) /
                           (1. / (m * m * m)) - 1.;
        ok = c.radius[b] == r && c.radiusExt[b] == re && c.deltaAv[b] == dAv &&
             c.nTracers[b] == set.nCoreTracers[b];
        const double exact = s * std::cbrt(3.0 * (double)set.nCells[b] / (4.0 * kPi));
        ulp = ulp && (c.radiusExt[b] == exact || std::nextafter(c.radiusExt[b], exact) == exact);
      }
    }
    check(ok, "the formulas");
    check(ulp, "radiusExt within one ulp of std::cbrt's");
  }

  group("delta2Mps: the sphere closed at 2 mps, the tracers not read not counted (outside the limits, "
        "outside a reused grid), -1 with no tracer in it");
  {
    const double m = 10.;
    Config config;
    config.log = nullptr;
    config.smoothing.sigma = 0.;
    const std::vector<std::uint8_t> allValid;
    const std::vector<std::int64_t> noIds;
    auto sphere = [&] (const std::vector<double>& tracers, const std::optional<otswap::BoxLimits>& limits,
                       const DivergenceGrid* reused, const std::array<std::size_t, 3>& seed,
                       const std::size_t expected, const std::string& what) {
      std::vector<double> displacement(tracers.size(), 0.);
      for (std::size_t i = 0; i < tracers.size(); i += 3) displacement[i] = 0.01;
      const internal::BoxSearch s = internal::check_box(tracers, displacement, allValid, m, config, limits,
                                                        noIds, reused);
      Hand h;
      h.grid = reused ? *reused : internal::box_grid(s, 1.);
      const internal::BoxDivergence div = internal::box_divergence(s, h.grid);
      const std::size_t n = h.grid.nCellsAxis[0];
      const std::int64_t c = row(n, seed[0], seed[1], seed[2]);
      h.grid.divergence.assign(h.grid.nCells, 1.);
      for (std::size_t i = seed[0] - 1; i <= seed[0] + 1; ++i)
        for (std::size_t j = seed[1] - 1; j <= seed[1] + 1; ++j)
          for (std::size_t k = seed[2] - 1; k <= seed[2] + 1; ++k)
            h.grid.divergence[row(n, i, j, k)] = 0.;
      h.grid.divergence[c] = -1.;
      internal::box_watershed(h.grid, h.basins);
      Catalogues cat;
      internal::box_basins(s, div, h.grid, h.basins, cat.voids, cat.overdense);
      const double nt = (cat.voids.delta2Mps[0] + 1.) * (4 * kPi * 8 / 3.);
      check(cat.voids.nBasins == 1 && cat.voids.delta2Mps[0] == delta2(expected, m),
            what + " (n_t " + std::to_string(nt) + ")");
      return cat;
    };

    std::vector<double> tracers;
    auto add = [&tracers] (const double x, const double y, const double z) {
      tracers.insert(tracers.end(), {x, y, z});
    };
    auto beyond = [] (const double centre, const double at) {
      double x = at;
      while ((x - centre) * (x - centre) <= 400.) x = std::nextafter(x, at < centre ? -kInf : kInf);
      return x;
    };
    auto lattice = [&add] (const std::size_t cells, const double x, const double y, const double z) {
      for (std::size_t i = 0; i < cells; ++i)
        for (std::size_t j = 0; j < cells; ++j)
          for (std::size_t k = 0; k < cells; ++k) {
            const double p[3] = {1.5 + 10. * i, 1.5 + 10. * j, 1.5 + 10. * k};
            if ((p[0] - x) * (p[0] - x) + (p[1] - y) * (p[1] - y) + (p[2] - z) * (p[2] - z) > 25. * 25.)
              add(p[0], p[1], p[2]);
          }
    };
    const otswap::BoxLimits limits = box(0., 93.);
    add(0., 0., 0.);
    add(93., 93., 93.);
    add(61.5, 41.5, 41.5);
    add(beyond(81.5, 61.5), 41.5, 41.5);
    add(81.5, 61.5, 41.5);
    add(81.5, 41.5, beyond(41.5, 21.5));
    add(95., 41.5, 41.5);
    add(81.5, 41.5, 41.5);
    lattice(10, 81.5, 41.5, 41.5);
    const DivergenceGrid built = internal::box_grid(
      internal::check_box(tracers, std::vector<double>(tracers.size(), 0.01), allValid, m, config, limits,
                          noIds, nullptr), 1.);
    check(built.lower[0] == -3.5 && built.cellSize == 10. && built.nCellsAxis[0] == 10,
          "the grid of the sphere tests: lower -3.5, 10 cells of 10");
    const Catalogues cat = sphere(tracers, limits, nullptr, {8, 4, 4}, 3,
                                  "centre (81.5, 41.5, 41.5): at 2 mps counted, one ulp beyond not, outside "
                                  "the limits not");
    check(cat.voids.centre[0] == 81.5 && cat.voids.centre[1] == 41.5 && cat.voids.centre[2] == 41.5,
          "the centre is the cell centre");

    tracers.clear();
    add(0., 0., 0.);
    add(93., 93., 93.);
    add(11.5, 31.5, 31.5);
    add(beyond(31.5, 11.5), 31.5, 31.5);
    add(31.5, 31.5, 11.5);
    add(47., 31.5, 31.5);
    add(31.5, 31.5, 31.5);
    lattice(5, 31.5, 31.5, 31.5);
    DivergenceGrid half = built;
    half.nCellsAxis = {5, 5, 5};
    half.nCells = 125;
    half.upper = {46.5, 46.5, 46.5};
    sphere(tracers, std::nullopt, &half, {3, 3, 3}, 3,
           "reused grid, centre (31.5, 31.5, 31.5): at 2 mps counted, one ulp beyond not, outside the grid "
           "not");

    tracers.clear();
    add(0., 0., 0.);
    add(93., 93., 93.);
    lattice(10, 41.5, 41.5, 41.5);
    const Catalogues empty = sphere(tracers, limits, nullptr, {4, 4, 4}, 0, "no tracer in the sphere");
    check(empty.voids.delta2Mps[0] == -1., "n_t = 0 gives -1 exactly");
  }

  group("edge flags 0, 1, 2 and contamination on hand-built partitions");
  {
    const std::size_t n = 4;
    DivergenceGrid g = cube(n, 1.);
    std::vector<std::int64_t> label(g.nCells, -1);
    for (std::size_t i = 1; i <= 2; ++i)
      for (std::size_t j = 1; j <= 2; ++j)
        for (std::size_t k = 1; k <= 2; ++k) {
          label[row(n, i, j, k)] = 0;
          g.divergence[row(n, i, j, k)] = -1.;
        }
    Hand h = partition(g, label);
    Catalogues cat = catalogues(h);
    check(cat.voids.edgeFlag[0] == 0 && cat.voids.contamination[0] == 0.,
          "an interior void: flag 0, contamination 0");
    check(cat.overdense.edgeFlag[0] == 2 && cat.overdense.contamination[0] == 96. / 120.,
          "the overdense shell around it: flag 2, 96 faces on the grid over 96 + 24");

    g = cube(n, 1.);
    label.assign(g.nCells, -1);
    label[row(n, 1, 1, 1)] = label[row(n, 0, 1, 1)] = 0;
    g.divergence[row(n, 1, 1, 1)] = -1.;
    g.divergence[row(n, 0, 1, 1)] = 0.5;
    cat = catalogues(partition(g, label));
    check(cat.voids.edgeFlag[0] == 1 && cat.voids.contamination[0] == 1. / 10.,
          "a void touching a face with its extended part only: flag 1, 1 face over 1 + 9");

    g.divergence[row(n, 0, 1, 1)] = -0.5;
    cat = catalogues(partition(g, label));
    check(cat.voids.edgeFlag[0] == 2 && cat.voids.contamination[0] == 1. / 10.,
          "a void touching a face with its core: flag 2");

    g = cube(n, 1.);
    g.divergence[row(n, 2, 1, 2)] = -1.;
    cat = catalogues(one_basin(g, row(n, 2, 1, 2)));
    check(cat.voids.edgeFlag[0] == 1 && cat.voids.contamination[0] == 1.,
          "a void filling the grid, core inside: flag 1, contamination 1");
    g.divergence[row(n, 0, 0, 0)] = -1.;
    cat = catalogues(one_basin(g, row(n, 2, 1, 2)));
    check(cat.voids.edgeFlag[0] == 2 && cat.voids.contamination[0] == 1.,
          "a void filling the grid, core on a face: flag 2, contamination 1");
    cat = catalogues(one_basin(cube(1, -1.), 0));
    check(cat.voids.edgeFlag[0] == 2 && cat.voids.contamination[0] == 1., "one cell: flag 2, contamination 1");
  }

  group("random partitions and watershed partitions: edge flags and contamination equal to a brute-force "
        "count with 1, 3 and 4 threads; overdense core = basin, its flag 0 or 2");
  for (const std::size_t n : {1, 2, 3, 6, 13}) {
    DivergenceGrid g = cube(n, 0.);
    for (double& d : g.divergence) d = uniform(-1., 1.);
    g.divergence[0] = -1.;
    std::vector<std::int64_t> label(g.nCells);
    for (std::size_t c = 0; c < g.nCells; ++c) label[c] = (std::int64_t)(rng() % 4) - 2;
    const std::int64_t firsts[4] = {0, 1, -1, -2};
    for (std::size_t c = 0; c < std::min<std::size_t>(4, g.nCells); ++c) label[c] = firsts[c];
    bool ok = true;
    for (const int t : {1, 3, 4}) {
      threads(t);
      const Hand h = partition(g, label);
      ok = ok && brute_force(h, catalogues(h));
    }
    Hand w;
    w.grid = g;
    internal::box_watershed(w.grid, w.basins);
    for (const int t : {1, 3, 4}) {
      threads(t);
      const Catalogues cat = catalogues(w);
      ok = ok && brute_force(w, cat) && in_seed(w.grid, w.basins, cat) &&
           w.basins.overdense.nCoreCells == w.basins.overdense.nCells;
      for (std::size_t b = 0; b < cat.overdense.nBasins; ++b)
        ok = ok && cat.overdense.edgeFlag[b] != 1 && cat.overdense.radius[b] == cat.overdense.radiusExt[b];
    }
    threads(4);
    check(ok, "n " + std::to_string(n));
  }

  group("assignment: cell, void or overdense basin, distance from its centre (NaN when not read), inMask, "
        "inBox, used");
  {
    const std::size_t n = 4;
    DivergenceGrid g = cube(n, 1.);
    std::vector<std::int64_t> label(g.nCells, -1);
    for (std::size_t c = 0; c < 32; ++c) {
      label[c] = 0;
      g.divergence[c] = -1.;
    }
    g.divergence[5] = -2.;
    const Hand h = partition(g, label, {5});
    const std::vector<double> tracers = {0.25, 1.5, 3.75, 3.5, 2.25, 0.125, 1.75, 1.25, 1.5};
    internal::BoxSearch search;
    search.nTracers = 3;
    search.nRealizations = 1;
    search.tracers = &tracers;
    search.inBox = {1, 1, 0};
    internal::BoxDivergence div;
    div.cell = {internal::box_cell(g, &tracers[0]), internal::box_cell(g, &tracers[3]), -1};
    div.used = {1, 0, 0};
    Catalogues cat;
    internal::box_basins(search, div, h.grid, h.basins, cat.voids, cat.overdense);
    const TracerAssignment as = internal::box_assignment(search, div, h.grid, cat.voids, cat.overdense);
    auto distance = [&] (const std::size_t i, const double* centre) {
      const double dx = tracers[3*i] - centre[0];
      const double dy = tracers[3*i+1] - centre[1];
      const double dz = tracers[3*i+2] - centre[2];
      return std::sqrt(dx*dx + dy*dy + dz*dz);
    };
    check(as.cellIndex == std::vector<std::int64_t>{div.cell[0], div.cell[1], -1} && div.cell[0] < 32 &&
          div.cell[1] >= 32, "cellIndex");
    check(as.voidIndex == std::vector<std::int64_t>{0, -1, -1} &&
          as.overdenseIndex == std::vector<std::int64_t>{-1, 0, -1}, "voidIndex, overdenseIndex");
    check(as.distance[0] == distance(0, &cat.voids.centre[0]) &&
          as.distance[1] == distance(1, &cat.overdense.centre[0]) && std::isnan(as.distance[2]),
          "distance from the void, from the overdense basin, NaN when not read");
    check(as.inMask == std::vector<std::uint8_t>{1, 1, 1} && as.inBox == search.inBox && as.used == div.used,
          "inMask all 1, inBox and used passed through");
  }

  group("the same catalogues and assignment with 1, 3 and 4 threads, and the recorded hashes");
  {
    const Generated in;
    for (const double cellSizeMps : {1., 0.5}) {
      const std::string what = cellSizeMps == 1. ? "cell 1 mps: " : "cell 0.5 mps: ";
      const bool count = cellSizeMps == 1.;
      bool agrees1 = false, agrees3 = false, agrees4 = false;
      threads(1);
      const std::uint64_t one = in.hash(cellSizeMps, count, agrees1);
      threads(3);
      const std::uint64_t three = in.hash(cellSizeMps, false, agrees3);
      threads(4);
      const std::uint64_t four = in.hash(cellSizeMps, false, agrees4);
      check(agrees1 && agrees3 && agrees4, what + "centres in their seeds, flags and contamination equal to "
            "the brute force, overdense core = basin" + (count ? ", delta2Mps equal to a direct count" : ""));
      check(one == three && one == four, what + "1, 3 and 4 threads give the same hash");

      const std::uint64_t recorded = cellSizeMps == 1. ? kRecordedHash1 : kRecordedHash05;
      if (one != recorded)
        std::cerr << what << "hash 0x" << std::hex << one << std::dec << std::endl;
      check(one == recorded, what + "the hash equals the recorded one");
    }
  }

  return report("test_basins");
}
