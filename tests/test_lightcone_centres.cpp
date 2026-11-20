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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bitvf/VoidFinder.h"
#include "otswap/Maths.h"

#include "check.h"
#include "internal.h"

using namespace bitvf;
using internal::LightconeAdjacency;
using internal::LightconeBasins;
using internal::LightconeGrid;

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;
  constexpr double kHalfPi = kPi / 2.;

  constexpr std::uint64_t kRecordedHash = 0xa8d414c75c8e3d49ull;

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  enum class Sky { Patch, Full, AcrossZero };

  LightconeGrid built_grid (const Sky kind)
  {
    std::mt19937 rng(1);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < 4000; ++i) {
      if (kind == Sky::Full) {
        sky.push_back(uniform(0., kTwoPi));
        sky.push_back(otswap::maths::asin(uniform(-1., 1.)));
        sky.push_back(uniform(0.001, 0.05));
      }
      else {
        const double ra = kind == Sky::Patch ? uniform(0.5, 2.) : uniform(-0.6, 0.6);
        sky.push_back(ra < 0. ? ra + kTwoPi : ra);
        sky.push_back(otswap::maths::asin(uniform(otswap::maths::sin(-0.3), otswap::maths::sin(0.3))));
        sky.push_back(uniform(0.25, 0.4));
      }
    }
    const internal::LightconeExtent e = internal::lightcone_extent(sky, {}, {}, {}, table());
    DivergenceGrid record;
    record.geometry = otswap::Geometry::Lightcone;
    record.distanceMin = e.distanceMin;
    record.distanceMax = e.distanceMax;
    record.raMin = e.raMin;
    record.raMax = e.raMax;
    record.decMin = e.decMin;
    record.decMax = e.decMax;
    record.cellSizeMps = kind == Sky::Full ? 8. : 9.;
    record.mpsProfile.redshift = {0., 0.5};
    record.mpsProfile.mps = {5., 6.};
    const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    return internal::lightcone_grid(record, table(), mask);
  }

  std::int64_t ring_of (const LightconeGrid& g, const std::int64_t c)
  {
    return std::upper_bound(g.ringFirstCell.begin(), g.ringFirstCell.end(), c) - g.ringFirstCell.begin() - 1;
  }

  std::int64_t shell_of (const LightconeGrid& g, const std::int64_t c)
  {
    const std::int64_t ring = ring_of(g, c);
    return std::upper_bound(g.shellFirstRing.begin(), g.shellFirstRing.end(), ring) - g.shellFirstRing.begin() - 1;
  }

  struct Centre {
    double x[3];
    double distance, ra, dec, z;
  };

  /// The centre of a single basin seeded at @p seed when only @p neighbour has a non-zero divergence, @p value.
  Centre forced (LightconeGrid& g, const LightconeAdjacency& a, const std::int64_t seed, const std::int64_t neighbour,
                 const double value, const bool isVoid)
  {
    g.cells.divergence.assign(g.cells.nCells, 0.);
    g.cells.divergenceSmoothed.clear();
    g.cells.divergence[(std::size_t)neighbour] = value;
    LightconeBasins b;
    (isVoid ? b.voids : b.overdense).seed.push_back(seed);
    BasinCatalog voids, overdense;
    internal::lightcone_centres(g, a, b, table(), voids, overdense);
    const BasinCatalog& cat = isVoid ? voids : overdense;
    Centre out;
    for (int i = 0; i < 3; ++i) out.x[i] = cat.centre[(std::size_t)i];
    internal::lightcone_sky(out.x, out.distance, out.ra, out.dec);
    out.z = cat.centreSky[2];
    check(cat.centreSky[0] == out.ra && cat.centreSky[1] == out.dec, "the sky centre is the sky of the centre");
    return out;
  }

  bool is_neighbour (const LightconeAdjacency& a, const std::int64_t c, const std::int64_t n)
  {
    return std::binary_search(a.radiusNeighbour.begin() + a.radiusFirst[(std::size_t)c],
                              a.radiusNeighbour.begin() + a.radiusFirst[(std::size_t)c + 1], n);
  }

  double dot_shift (const DivergenceGrid& cells, const Centre& c, const std::int64_t seed, const std::int64_t n)
  {
    const double* x0 = &cells.centre[3 * (std::size_t)seed];
    const double* xn = &cells.centre[3 * (std::size_t)n];
    double s = 0.;
    for (int i = 0; i < 3; ++i) s += (c.x[i] - x0[i]) * (xn[i] - x0[i]);
    return s;
  }

  void smooth_field (DivergenceGrid& cells)
  {
    cells.divergence.assign(cells.nCells, 0.);
    for (std::size_t c = 0; c < cells.nCells; ++c) {
      const double* x = &cells.centre[3*c];
      cells.divergence[c] = otswap::maths::sin(0.021 * x[0]) * otswap::maths::cos(0.017 * x[1]) +
                            otswap::maths::sin(0.013 * x[2] + 0.5) +
                            0.8 * otswap::maths::sin(0.11 * x[0]) * otswap::maths::sin(0.13 * x[1]) + 0.3;
    }
    cells.divergenceSmoothed.clear();
  }

  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  std::uint64_t centres_hash ()
  {
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const Sky kind : {Sky::Patch, Sky::Full, Sky::AcrossZero}) {
      LightconeGrid g = built_grid(kind);
      for (std::size_t c = 5; c < g.cells.nCells; c += 29) g.cells.masked[c] = 1;
      smooth_field(g.cells);
      const LightconeAdjacency a = internal::lightcone_adjacency(g);
      LightconeBasins b = internal::lightcone_watershed(a, g.cells);
      BasinCatalog voids, overdense;
      internal::lightcone_centres(g, a, b, table(), voids, overdense);
      for (const BasinCatalog* cat : {&voids, &overdense}) {
        hash_bytes(h, cat->centre.data(), cat->centre.size() * sizeof(double));
        hash_bytes(h, cat->centreSky.data(), cat->centreSky.size() * sizeof(double));
      }
    }
    return h;
  }

}


int main ()
{
  // A void's centre moves towards its lower neighbours and away from its higher ones; an overdense basin's centre moves towards its higher neighbours and away from its lower ones.
  group("the direction of the centre of a void and of an overdense basin");
  {
    LightconeGrid g = built_grid(Sky::Patch);
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    const std::int64_t seed = (std::int64_t)g.cells.nCells / 2;
    const std::int64_t n = a.radiusNeighbour[(std::size_t)a.radiusFirst[(std::size_t)seed]];
    check(dot_shift(g.cells, forced(g, a, seed, n, -1., true), seed, n) > 0., "a void moves towards a lower cell");
    check(dot_shift(g.cells, forced(g, a, seed, n, 1., true), seed, n) < 0., "a void moves away from a higher cell");
    check(dot_shift(g.cells, forced(g, a, seed, n, 1., false), seed, n) > 0.,
          "an overdense basin moves towards a higher cell");
    check(dot_shift(g.cells, forced(g, a, seed, n, -1., false), seed, n) < 0.,
          "an overdense basin moves away from a lower cell");
  }

  // Every centre lies in its seed cell; a centre the neighbour sum already places there is that sum, bit for bit, and its redshift is the distance table's.
  group("the centres of a watershed: the neighbour sum, kept inside the seed cell");
  {
    for (const Sky kind : {Sky::Patch, Sky::Full, Sky::AcrossZero}) {
      LightconeGrid g = built_grid(kind);
      smooth_field(g.cells);
      const LightconeAdjacency a = internal::lightcone_adjacency(g);
      LightconeBasins b = internal::lightcone_watershed(a, g.cells);
      BasinCatalog voids, overdense;
      internal::lightcone_centres(g, a, b, table(), voids, overdense);
      std::size_t kept = 0, inSeed = 0, total = 0, sameSky = 0, sameZ = 0;
      for (int kind2 = 0; kind2 < 2; ++kind2) {
        const internal::BasinSet& set = kind2 == 0 ? b.voids : b.overdense;
        const BasinCatalog& cat = kind2 == 0 ? voids : overdense;
        const double sign = kind2 == 0 ? -1. : 1.;
        for (std::size_t k = 0; k < set.seed.size(); ++k) {
          const std::int64_t c = set.seed[k];
          const double* x0 = &g.cells.centre[3 * (std::size_t)c];
          double num[3] = {0., 0., 0.}, den = 0.;
          for (std::int64_t p = a.radiusFirst[(std::size_t)c]; p < a.radiusFirst[(std::size_t)c + 1]; ++p) {
            const std::size_t n = (std::size_t)a.radiusNeighbour[(std::size_t)p];
            for (int i = 0; i < 3; ++i) num[i] += g.cells.divergence[n] * (g.cells.centre[3 * n + (std::size_t)i] - x0[i]);
            den += std::abs(g.cells.divergence[n]);
          }
          if (den == 0.) den = 1.;
          double expected[3];
          for (int i = 0; i < 3; ++i) expected[i] = x0[i] + sign * (num[i] / den);
          const double* x = &cat.centre[3 * k];
          ++total;
          inSeed += internal::lightcone_cell_at(g, x) == c;
          if (internal::lightcone_cell_at(g, expected) == c)
            kept += x[0] == expected[0] && x[1] == expected[1] && x[2] == expected[2];
          else
            kept += !(x[0] == expected[0] && x[1] == expected[1] && x[2] == expected[2]);
          double r, ra, dec;
          internal::lightcone_sky(x, r, ra, dec);
          sameSky += cat.centreSky[3 * k] == ra && cat.centreSky[3 * k + 1] == dec && ra >= 0. && ra < kTwoPi;
          sameZ += cat.centreSky[3 * k + 2] == table().redshiftAt(r);
        }
      }
      check(total > 0 && inSeed == total, "every centre is in its seed cell");
      check(kept == total, "a sum inside the seed cell is kept bit for bit, one outside is moved");
      check(sameSky == total, "the sky centre is the sky of the centre, right ascension in [0, 2 pi)");
      check(sameZ == total, "the redshift is the distance table's at the centre's distance");
    }
  }

  // A centre pulled beyond the seed cell's outer distance edge stops just inside it, its direction from the observer unchanged.
  group("a centre beyond the seed cell's distance range is clamped in distance only");
  {
    LightconeGrid g = built_grid(Sky::Patch);
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    std::int64_t seed = -1, n = -1;
    for (std::int64_t c = 0; c < (std::int64_t)g.cells.nCells && seed < 0; ++c)
      for (std::int64_t p = a.radiusFirst[(std::size_t)c]; p < a.radiusFirst[(std::size_t)c + 1]; ++p) {
        const std::int64_t q = a.radiusNeighbour[(std::size_t)p];
        const double* xs = &g.cells.centreSky[3 * (std::size_t)q];
        const std::size_t cc = (std::size_t)c;
        if (shell_of(g, q) == shell_of(g, c) + 1 && xs[1] >= g.cells.cellDecMin[cc] && xs[1] < g.cells.cellDecMax[cc] &&
            xs[0] >= g.cells.cellRaMin[cc] && xs[0] < g.cells.cellRaMax[cc]) {
          seed = c;
          n = q;
          break;
        }
      }
    check(seed >= 0, "the grid has a cell with a neighbour straight behind it");
    if (seed >= 0) {
      const Centre c = forced(g, a, seed, n, -1., true);
      const double rHi = g.cells.cellDistanceMax[(std::size_t)seed];
      check(internal::lightcone_cell_at(g, c.x) == seed, "the centre is in the seed cell");
      check(c.distance < rHi && c.distance > rHi * (1. - 1e-12), "the distance is just inside the outer edge");
      check(std::abs(c.ra - g.cells.centreSky[3 * (std::size_t)n]) < 1e-12 &&
            std::abs(c.dec - g.cells.centreSky[3 * (std::size_t)n + 1]) < 1e-12, "the direction is the neighbour's");
    }
  }

  // In the central sphere only the distance is clamped: the centre keeps the direction it was pulled towards.
  group("the central sphere clamps the distance only");
  {
    LightconeGrid g = built_grid(Sky::Full);
    check(g.centralSphere, "the full-sky grid has a central sphere");
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    const std::int64_t n = a.radiusNeighbour[(std::size_t)a.radiusFirst[0]];
    const Centre c = forced(g, a, 0, n, -1., true);
    double rn, ran, decn;
    internal::lightcone_sky(&g.cells.centre[3 * (std::size_t)n], rn, ran, decn);
    check(internal::lightcone_cell_at(g, c.x) == 0, "the centre is in the central sphere");
    check(c.distance < g.shellEdge[1] && c.distance > g.shellEdge[1] * (1. - 1e-12),
          "the distance is just inside the sphere");
    check(std::abs(c.ra - ran) < 1e-12 && std::abs(c.dec - decn) < 1e-12, "the direction is the neighbour's");
  }

  // In a polar cap only the declination is clamped, to the cap's edge; the right ascension is kept, and the centre stays in the cap.
  group("a polar cap clamps the declination only, north and south");
  {
    LightconeGrid g = built_grid(Sky::Full);
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    std::size_t tried = 0;
    for (const bool north : {true, false}) {
      std::int64_t seed = -1, n = -1;
      for (std::size_t ring = 0; ring + 1 < g.ringFirstCell.size() && seed < 0; ++ring) {
        if (!g.ringPole[ring] || (g.centralSphere && (std::int64_t)ring < g.shellFirstRing[1])) continue;
        if (north ? g.ringDecMax[ring] != kHalfPi : g.ringDecMin[ring] != -kHalfPi) continue;
        const std::int64_t c = g.ringFirstCell[ring];
        for (std::int64_t p = a.radiusFirst[(std::size_t)c]; p < a.radiusFirst[(std::size_t)c + 1]; ++p) {
          const std::int64_t q = a.radiusNeighbour[(std::size_t)p];
          if (shell_of(g, q) == shell_of(g, c) && ring_of(g, q) != (std::int64_t)ring) {
            seed = c;
            n = q;
            break;
          }
        }
      }
      if (seed < 0) continue;
      ++tried;
      const Centre c = forced(g, a, seed, n, -1., true);
      const std::size_t ring = (std::size_t)ring_of(g, seed);
      double rn, ran, decn;
      internal::lightcone_sky(&g.cells.centre[3 * (std::size_t)n], rn, ran, decn);
      check(internal::lightcone_cell_at(g, c.x) == seed, "the centre is in the cap");
      const double edge = north ? g.ringDecMin[ring] : g.ringDecMax[ring];
      check(north ? c.dec >= edge && c.dec < edge + 1e-9 : c.dec < edge && c.dec > edge - 1e-9,
            "the declination is at the cap's edge");
      check(std::abs(c.ra - ran) < 1e-9, "the right ascension is the neighbour's");
    }
    check(tried == 2, "the full-sky grid has a north and a south cap");
  }

  // A centre pulled across right ascension 0 is clamped to the seed sector's edge on its own side, and its right ascension is reported in [0, 2 pi).
  group("a centre pulled across right ascension 0, on a full circle and on a partial ring");
  {
    LightconeGrid g = built_grid(Sky::Full);
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    std::size_t tried = 0;
    for (std::size_t ring = 0; ring + 1 < g.ringFirstCell.size() && tried == 0; ++ring) {
      const std::int64_t firstCell = g.ringFirstCell[ring], lastCell = g.ringFirstCell[ring + 1] - 1;
      if (!g.ringFull[ring] || lastCell - firstCell < 3 || !is_neighbour(a, firstCell, lastCell)) continue;
      ++tried;
      const Centre low = forced(g, a, firstCell, lastCell, -1., true);
      check(internal::lightcone_cell_at(g, low.x) == firstCell, "the first sector keeps its centre");
      check(low.ra >= 0. && low.ra < 1e-9, "pulled below 0, the right ascension stops at 0");
      const Centre high = forced(g, a, lastCell, firstCell, -1., true);
      check(internal::lightcone_cell_at(g, high.x) == lastCell, "the last sector keeps its centre");
      check(high.ra < kTwoPi && high.ra > kTwoPi - 1e-9, "pulled above 2 pi, the right ascension stops below 2 pi");
    }
    check(tried == 1, "the full-sky grid has a full ring with its first and last sectors as neighbours");

    LightconeGrid p = built_grid(Sky::AcrossZero);
    const LightconeAdjacency pa = internal::lightcone_adjacency(p);
    std::size_t beyond = 0, beyondOk = 0;
    for (std::size_t ring = 0; ring + 1 < p.ringFirstCell.size(); ++ring) {
      if (p.ringFull[ring]) continue;
      const std::int64_t c = p.ringFirstCell[ring + 1] - 1;
      const std::size_t cc = (std::size_t)c;
      if (!(p.cells.cellRaMin[cc] >= kTwoPi) || !is_neighbour(pa, c, c - 1)) continue;
      ++beyond;
      const Centre centre = forced(p, pa, c, c - 1, 1., true);
      beyondOk += internal::lightcone_cell_at(p, centre.x) == c && centre.ra >= p.cells.cellRaMin[cc] - kTwoPi &&
                  centre.ra < p.cells.cellRaMax[cc] - kTwoPi;
    }
    check(beyond > 0 && beyondOk == beyond,
          "on a ring across 0, a centre in a sector past 2 pi is reported in [0, 2 pi), inside its sector");
  }

  // Centres and sky centres of fixed grids are identical bit for bit with any number of threads and match the recorded hash.
  group("the centres are the same with 1 and 4 threads, and match the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = centres_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = centres_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "centres hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_centres");
}
