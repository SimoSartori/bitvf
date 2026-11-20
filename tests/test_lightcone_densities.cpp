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
using internal::LightconeBasins;
using internal::LightconeGrid;

namespace {

  constexpr std::uint64_t kRecordedHash = 0xc636a64d8eb7918aull;

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  LightconeGrid built_grid ()
  {
    std::mt19937 rng(1);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < 4000; ++i) {
      sky.push_back(uniform(0.5, 2.));
      sky.push_back(otswap::maths::asin(uniform(otswap::maths::sin(-0.3), otswap::maths::sin(0.3))));
      sky.push_back(uniform(0.25, 0.4));
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
    record.cellSizeMps = 9.;
    record.mpsProfile.redshift = {0., 0.5};
    record.mpsProfile.mps = {5., 6.};
    const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    return internal::lightcone_grid(record, table(), mask);
  }

  /// Adds a basin centred on cell @p c, with @p nTracers core tracers and core volume @p volume.
  void add_basin (const LightconeGrid& g, const std::int64_t c, const std::size_t nTracers, const double volume,
                  BasinCatalog& cat, std::vector<double>& coreVolume)
  {
    const double* x = &g.cells.centre[3 * (std::size_t)c];
    double distance = 0., ra = 0., dec = 0.;
    internal::lightcone_sky(x, distance, ra, dec);
    cat.centre.insert(cat.centre.end(), x, x + 3);
    cat.centreSky.insert(cat.centreSky.end(), {ra, dec, table().redshiftAt(distance)});
    cat.nTracers.push_back(nTracers);
    coreVolume.push_back(volume);
    ++cat.nBasins;
  }

  /// The point at @p fraction of @p length from @p x along the unit vector @p u.
  void push_at (std::vector<double>& points, const double* x, const double* u, const double length,
                const double fraction)
  {
    for (int a = 0; a < 3; ++a) points.push_back(x[a] + fraction * length * u[a]);
  }

  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  std::uint64_t densities_hash ()
  {
    const LightconeGrid g = built_grid();
    std::mt19937 rng(7);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> tracers, randoms;
    for (std::size_t c = 0; c < g.cells.nCells; ++c)
      for (int k = 0; k < 40; ++k) {
        std::vector<double>& points = k < 4 ? tracers : randoms;
        const double f = uniform(0.998, 1.002);
        for (int a = 0; a < 3; ++a) points.push_back(f * g.cells.centre[3*c + (std::size_t)a] + uniform(-2., 2.));
      }
    LightconeBasins basins;
    BasinCatalog voids, overdense;
    for (std::size_t c = 0; c < g.cells.nCells; c += 7) {
      if (c % 2 == 0) add_basin(g, (std::int64_t)c, c % 11, g.cells.volume[c], voids, basins.voidCoreVolume);
      else add_basin(g, (std::int64_t)c, c % 13, g.cells.volume[c], overdense, basins.overdenseCoreVolume);
    }
    internal::lightcone_densities(g, basins, tracers, randoms, voids, overdense);
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const BasinCatalog* cat : {&voids, &overdense}) {
      hash_bytes(h, cat->deltaAv.data(), cat->deltaAv.size() * sizeof(double));
      hash_bytes(h, cat->delta2Mps.data(), cat->delta2Mps.size() * sizeof(double));
    }
    return h;
  }

}


int main ()
{
  const LightconeGrid g = built_grid();
  std::int64_t outer = 0, inner = 0;
  for (std::size_t c = 0; c < g.cells.nCells; ++c) {
    if (g.cells.cellDistanceMax[c] > g.cells.cellDistanceMax[(std::size_t)outer]) outer = (std::int64_t)c;
    if (g.cells.cellDistanceMin[c] < g.cells.cellDistanceMin[(std::size_t)inner]) inner = (std::int64_t)c;
  }

  // deltaAv is the core tracers over the core volume, over mps(z)^-3 at the centre's redshift, minus 1; delta2Mps is n_t / (alpha n_r) - 1 within 2 mps(z) of the centre, alpha the ratio of all tracers to all randoms, and NaN with no random in the sphere.
  group("deltaAv and delta2Mps follow mps at the centre's redshift and the counts within 2 mps");
  {
    LightconeBasins basins;
    BasinCatalog voids, overdense;
    add_basin(g, inner, 7, 1234.5, voids, basins.voidCoreVolume);
    add_basin(g, outer, 0, 99., overdense, basins.overdenseCoreVolume);
    const double m = g.cells.mpsProfile.at(voids.centreSky[2]);
    const double R = 2 * m;
    const double* x = voids.centre.data();
    const double ux[3] = {1., 0., 0.}, uy[3] = {0., 1., 0.}, uz[3] = {0., 0., 1.};
    std::vector<double> tracers, randoms;
    push_at(tracers, x, ux, R, 0.5);
    push_at(tracers, x, uy, R, -0.9);
    push_at(tracers, x, uz, R, 0.99);
    push_at(tracers, x, ux, R, 1.01);
    push_at(tracers, x, uy, R, -3.);
    push_at(randoms, x, ux, R, 0.);
    push_at(randoms, x, ux, R, -0.2);
    push_at(randoms, x, uy, R, 0.7);
    push_at(randoms, x, uz, R, -0.98);
    for (int k = 0; k < 5; ++k) push_at(randoms, x, uz, R, 1.02 + k);
    internal::lightcone_densities(g, basins, tracers, randoms, voids, overdense);

    check(voids.deltaAv[0] == (7. / 1234.5) / (1. / (m * m * m)) - 1., "deltaAv of the void, exactly");
    const double mo = g.cells.mpsProfile.at(overdense.centreSky[2]);
    check(overdense.deltaAv[0] == (0. / 99.) / (1. / (mo * mo * mo)) - 1., "deltaAv of the overdense basin");
    const double alpha = 5. / 9.;
    check(voids.delta2Mps[0] == 3. / (alpha * 4.) - 1., "delta2Mps from 3 tracers and 4 randoms in the sphere");
    check(std::isnan(overdense.delta2Mps[0]), "no random in the sphere gives NaN");
  }

  // A centre in the outermost cell of the grid, far from every tracer and random, is counted without error.
  group("a centre beyond the data, inside the grid, is counted");
  {
    LightconeBasins basins;
    BasinCatalog voids, overdense;
    add_basin(g, outer, 3, 50., voids, basins.voidCoreVolume);
    const double* x = &g.cells.centre[3 * (std::size_t)inner];
    std::vector<double> tracers(x, x + 3), randoms(x, x + 3);
    internal::lightcone_densities(g, basins, tracers, randoms, voids, overdense);
    check(voids.delta2Mps.size() == 1 && std::isnan(voids.delta2Mps[0]), "the outer centre has an empty sphere");
    check(overdense.deltaAv.empty() && overdense.delta2Mps.empty(), "no overdense basin, empty columns");
  }

  // deltaAv and delta2Mps are identical bit for bit with any number of threads and match the recorded hash.
  group("the densities are the same with 1 and 4 threads, and match the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = densities_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = densities_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "densities hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_densities");
}
