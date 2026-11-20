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
 *  @file src/LightconeCentres.cpp
 *
 *  @brief The centres of the basins of a lightcone search.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <string>
#include <vector>

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;
  constexpr unsigned kMaxSteps = 1u << 20;
  constexpr double kInfinity = std::numeric_limits<double>::infinity();

  void cartesian (const double distance, const double ra, const double dec, double* x)
  {
    const double cosDec = otswap::maths::cos(dec);
    x[0] = distance * cosDec * otswap::maths::cos(ra);
    x[1] = distance * cosDec * otswap::maths::sin(ra);
    x[2] = distance * otswap::maths::sin(dec);
  }

}


// ============================================================================


void bitvf::internal::lightcone_centres (const LightconeGrid& grid, const LightconeAdjacency& adjacency,
                                         const LightconeBasins& basins, const otswap::DistanceTable& distances,
                                         BasinCatalog& voids, BasinCatalog& overdense)
{
  const DivergenceGrid& cells = grid.cells;
  const double* D = (cells.divergenceSmoothed.empty() ? cells.divergence : cells.divergenceSmoothed).data();
  const double* X = cells.centre.data();
  const std::int64_t* first = adjacency.radiusFirst.data();
  const std::int64_t* neighbour = adjacency.radiusNeighbour.data();
  const std::vector<std::int64_t>& ringFirstCell = grid.ringFirstCell;
  const std::vector<std::int64_t>& shellFirstRing = grid.shellFirstRing;

  const BasinSet* sets[2] = {&basins.voids, &basins.overdense};
  BasinCatalog* catalogues[2] = {&voids, &overdense};
  const char* names[2] = {"void", "overdense basin"};

  for (int kind = 0; kind < 2; ++kind) {
    const std::vector<std::int64_t>& seed = sets[kind]->seed;
    BasinCatalog& cat = *catalogues[kind];
    const std::int64_t B = (std::int64_t)seed.size();
    const double sign = kind == 0 ? -1. : 1.;
    cat.centre.assign(3 * (std::size_t)B, 0.);
    cat.centreSky.assign(3 * (std::size_t)B, 0.);

    std::int64_t notPlaced = B, noRedshift = B;
    std::exception_ptr redshiftFailure;

#pragma omp parallel for schedule(static)
    for (std::int64_t b = 0; b < B; ++b) {
      const std::int64_t c = seed[(std::size_t)b];
      const std::size_t cc = (std::size_t)c;
      const double* x0 = X + 3 * cc;

      double num[3] = {0., 0., 0.};
      double den = 0.;
      for (std::int64_t p = first[c]; p < first[c + 1]; ++p) {
        const std::size_t n = (std::size_t)neighbour[p];
        const double d = D[n];
        num[0] += d * (X[3 * n] - x0[0]);
        num[1] += d * (X[3 * n + 1] - x0[1]);
        num[2] += d * (X[3 * n + 2] - x0[2]);
        den += std::abs(d);
      }
      if (den == 0.) den = 1.;

      double x[3];
      for (int a = 0; a < 3; ++a) x[a] = x0[a] + sign * (num[a] / den);

      if (lightcone_cell_at(grid, x) != c) {
        const std::size_t ring = (std::size_t)(std::upper_bound(ringFirstCell.begin(), ringFirstCell.end(), c) -
                                               ringFirstCell.begin()) - 1;
        const std::size_t shell = (std::size_t)(std::upper_bound(shellFirstRing.begin(), shellFirstRing.end(),
                                                                 (std::int64_t)ring) - shellFirstRing.begin()) - 1;
        const bool central = grid.centralSphere && shell == 0;
        const bool cap = !central && grid.ringPole[ring];

        double r = 0., ra = 0., dec = 0.;
        lightcone_sky(x, r, ra, dec);
        const double raLo = cells.cellRaMin[cc], raHi = cells.cellRaMax[cc];
        ra += kTwoPi * std::round((0.5 * (raLo + raHi) - ra) / kTwoPi);
        r = std::max(cells.cellDistanceMin[cc], std::min(cells.cellDistanceMax[cc], r));
        if (!central) dec = std::max(cells.cellDecMin[cc], std::min(cells.cellDecMax[cc], dec));
        if (!central && !cap) ra = std::max(raLo, std::min(raHi, ra));
        cartesian(r, ra, dec, x);

        unsigned steps = 0;
        while (lightcone_cell_at(grid, x) != c) {
          if (++steps > kMaxSteps) {
#pragma omp critical(bitvf_lightcone_centres)
            notPlaced = std::min(notPlaced, b);
            break;
          }
          double r2 = 0., ra2 = 0., dec2 = 0.;
          lightcone_sky(x, r2, ra2, dec2);
          const LightconeOutside out = lightcone_outside(grid, c, r2, ra2, dec2);
          if (out.side[0]) r = std::nextafter(r, out.side[0] > 0 ? -kInfinity : kInfinity);
          if (out.side[1]) dec = std::nextafter(dec, out.side[1] > 0 ? -kInfinity : kInfinity);
          if (out.side[2]) ra = std::nextafter(ra, out.side[2] > 0 ? -kInfinity : kInfinity);
          cartesian(r, ra, dec, x);
        }
      }

      double distance = 0., ra = 0., dec = 0.;
      lightcone_sky(x, distance, ra, dec);
      double z = 0.;
      try {
        z = distances.redshiftAt(distance);
      }
      catch (...) {
#pragma omp critical(bitvf_lightcone_centres)
        if (b < noRedshift) {
          noRedshift = b;
          redshiftFailure = std::current_exception();
        }
      }

      double* centre = &cat.centre[3 * (std::size_t)b];
      double* sky = &cat.centreSky[3 * (std::size_t)b];
      centre[0] = x[0];
      centre[1] = x[1];
      centre[2] = x[2];
      sky[0] = ra;
      sky[1] = dec;
      sky[2] = z;
    }

    if (notPlaced != B)
      throw Error("internal error: the centre of " + std::string(names[kind]) + " " + std::to_string(notPlaced) +
                  " (from 0) cannot be placed in its seed cell " + std::to_string(seed[(std::size_t)notPlaced]) +
                  " (from 0); please report");
    if (redshiftFailure) {
      try {
        std::rethrow_exception(redshiftFailure);
      }
      catch (...) {
        std::throw_with_nested(Error("internal error: the redshift lookup failed at the centre of " +
                                     std::string(names[kind]) + " " + std::to_string(noRedshift) +
                                     " (from 0); please report"));
      }
    }
  }
}
