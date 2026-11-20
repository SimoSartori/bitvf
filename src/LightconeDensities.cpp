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
 *  @file src/LightconeDensities.cpp
 *
 *  @brief The density contrasts of the basins of a lightcone search.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "meshsearch/MeshGrid.h"
#include "otswap/Maths.h"

#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kHalfPi = kPi / 2.;

  std::unique_ptr<meshsearch::MeshGrid> mesh_of (const std::vector<double>& points, const double cellSize,
                                                 const std::vector<std::vector<double>>& limits)
  {
    const std::size_t N = points.size() / 3;
    std::vector<double> X(N), Y(N), Z(N);
    for (std::size_t i = 0; i < N; ++i) {
      X[i] = points[3*i];
      Y[i] = points[3*i+1];
      Z[i] = points[3*i+2];
    }
    try {
      return std::unique_ptr<meshsearch::MeshGrid>(new meshsearch::MeshGrid(X, Y, Z, cellSize, limits));
    }
    catch (const meshsearch::Error& e) {
      std::throw_with_nested(bitvf::Error(std::string("internal error: meshsearch: ") + e.what()));
    }
  }

}


// ============================================================================


void bitvf::internal::lightcone_densities (const LightconeGrid& grid, const LightconeBasins& basins,
                                           const std::vector<double>& tracers, const std::vector<double>& randoms,
                                           BasinCatalog& voids, BasinCatalog& overdense)
{
  const DivergenceGrid& cells = grid.cells;
  const std::size_t nVoid = voids.nBasins;
  const std::size_t nAll = nVoid + overdense.nBasins;
  const std::vector<double>* coreVolume[2] = {&basins.voidCoreVolume, &basins.overdenseCoreVolume};
  BasinCatalog* catalogues[2] = {&voids, &overdense};
  for (BasinCatalog* cat : catalogues) {
    cat->deltaAv.assign(cat->nBasins, 0.);
    cat->delta2Mps.assign(cat->nBasins, 0.);
  }

  constexpr std::int64_t none = std::numeric_limits<std::int64_t>::max();
  std::int64_t noMps = none;
  std::vector<double> radius(nAll, 0.);

#pragma omp parallel for schedule(static) reduction(min:noMps)
  for (std::int64_t g = 0; g < (std::int64_t)nAll; ++g) {
    const int kind = (std::size_t)g < nVoid ? 0 : 1;
    const std::size_t b = kind == 0 ? (std::size_t)g : (std::size_t)g - nVoid;
    BasinCatalog& cat = *catalogues[kind];
    double m = 0.;
    try {
      m = cells.mpsProfile.at(cat.centreSky[3*b+2]);
    }
    catch (...) {
      noMps = std::min(noMps, g);
      continue;
    }
    const double avDens = 1. / (m * m * m);
    cat.deltaAv[b] = ((double)cat.nTracers[b] / (*coreVolume[kind])[b]) / avDens - 1.;
    radius[(std::size_t)g] = 2 * m;
  }

  if (noMps != none) {
    const bool isVoid = (std::size_t)noMps < nVoid;
    const std::size_t b = isVoid ? (std::size_t)noMps : (std::size_t)noMps - nVoid;
    const double z = (isVoid ? voids : overdense).centreSky[3*b+2];
    try {
      cells.mpsProfile.at(z);
    }
    catch (const std::exception&) {
      std::throw_with_nested(Error(std::string("internal error: the mps lookup failed at the centre of ") +
                                   (isVoid ? "void " : "overdense basin ") + std::to_string(b) +
                                   " (from 0); please report"));
    }
    throw Error(std::string("internal error: the mps lookup failed at the centre of ") +
                (isVoid ? "void " : "overdense basin ") + std::to_string(b) + " (from 0); please report");
  }

  const double pad = nAll == 0 ? 0. : *std::max_element(radius.begin(), radius.end());
  double lo[3] = {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                  std::numeric_limits<double>::infinity()};
  double hi[3] = {-lo[0], -lo[1], -lo[2]};
  for (std::size_t c = 0; c < cells.nCells; ++c) {
    const double r[2] = {cells.cellDistanceMin[c], cells.cellDistanceMax[c]};
    double dec[3] = {cells.cellDecMin[c], cells.cellDecMax[c], 0.};
    const int nDec = dec[0] < 0. && dec[1] > 0. ? 3 : 2;
    double ra[7] = {cells.cellRaMin[c], cells.cellRaMax[c]};
    int nRa = 2;
    for (double k = std::ceil(ra[0] / kHalfPi); k * kHalfPi < ra[1] && nRa < 7; ++k) ra[nRa++] = k * kHalfPi;
    for (int i = 0; i < nDec; ++i) {
      const double cosDec = otswap::maths::cos(dec[i]);
      const double sinDec = otswap::maths::sin(dec[i]);
      for (int j = 0; j < nRa; ++j) {
        const double cosRa = otswap::maths::cos(ra[j]);
        const double sinRa = otswap::maths::sin(ra[j]);
        for (const double d : r) {
          const double x[3] = {d * cosDec * cosRa, d * cosDec * sinRa, d * sinDec};
          for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], x[a]);
            hi[a] = std::max(hi[a], x[a]);
          }
        }
      }
    }
  }
  std::vector<std::vector<double>> limits(3);
  for (int a = 0; a < 3; ++a) limits[(std::size_t)a] = {lo[a] - pad, hi[a] + pad};

  const double cellSize = 4 * grid.mpsMin;
  const char* objects[2] = {"tracers", "randoms"};
  std::vector<std::size_t> count[2];
  for (int set = 0; set < 2; ++set) {
    const std::unique_ptr<meshsearch::MeshGrid> mesh = mesh_of(set == 0 ? tracers : randoms, cellSize, limits);
    count[set].assign(nAll, 0);
    std::int64_t failed = none;

#pragma omp parallel for schedule(dynamic) reduction(min:failed)
    for (std::int64_t g = 0; g < (std::int64_t)nAll; ++g) {
      const BasinCatalog& cat = (std::size_t)g < nVoid ? voids : overdense;
      const std::size_t b = (std::size_t)g < nVoid ? (std::size_t)g : (std::size_t)g - nVoid;
      const double* x = &cat.centre[3*b];
      try {
        count[set][(std::size_t)g] = mesh->closeObjects(x[0], x[1], x[2], radius[(std::size_t)g]).size();
      }
      catch (...) {
        failed = std::min(failed, g);
      }
    }

    if (failed != none) {
      const bool isVoid = (std::size_t)failed < nVoid;
      const BasinCatalog& cat = isVoid ? voids : overdense;
      const std::size_t b = isVoid ? (std::size_t)failed : (std::size_t)failed - nVoid;
      const double* x = &cat.centre[3*b];
      const std::string what = std::string("internal error: counting the ") + objects[set] +
        " within 2 mps of the centre of " + (isVoid ? "void " : "overdense basin ") + std::to_string(b) +
        " (from 0)";
      try {
        mesh->closeObjects(x[0], x[1], x[2], radius[(std::size_t)failed]);
      }
      catch (const meshsearch::Error& e) {
        std::throw_with_nested(Error(what + ": meshsearch: " + e.what()));
      }
      throw Error(what + " failed inside the parallel region");
    }
  }

  const double alpha = (double)(tracers.size() / 3) / (double)(randoms.size() / 3);
  for (std::size_t g = 0; g < nAll; ++g) {
    BasinCatalog& cat = g < nVoid ? voids : overdense;
    const std::size_t b = g < nVoid ? g : g - nVoid;
    cat.delta2Mps[b] = count[1][g] == 0 ? std::numeric_limits<double>::quiet_NaN()
                                        : (double)count[0][g] / (alpha * (double)count[1][g]) - 1.;
  }
}
