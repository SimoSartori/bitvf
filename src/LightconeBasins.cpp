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
 *  @file src/LightconeBasins.cpp
 *
 *  @brief The volumes, boundary contact and contamination of the basins of
 *  a lightcone search.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;
  constexpr double kHalfPi = kPi / 2.;

}


// ============================================================================


double bitvf::internal::lightcone_shared_area (const DivergenceGrid& grid, const std::size_t c, const std::size_t n)
{
  const double rLo = grid.cellDistanceMin[c], rHi = grid.cellDistanceMax[c];
  const double nrLo = grid.cellDistanceMin[n], nrHi = grid.cellDistanceMax[n];
  const double dLo = grid.cellDecMin[c], dHi = grid.cellDecMax[c];
  const double ndLo = grid.cellDecMin[n], ndHi = grid.cellDecMax[n];
  const double aLo = grid.cellRaMin[c], aHi = grid.cellRaMax[c];
  double naLo = grid.cellRaMin[n], naHi = grid.cellRaMax[n];
  const double shift = kTwoPi * std::round((0.5 * (aLo + aHi) - 0.5 * (naLo + naHi)) / kTwoPi);
  naLo += shift;
  naHi += shift;

  const double tolR = 1e-9 * (rHi - rLo), tolD = 1e-9 * (dHi - dLo), tolA = 1e-9 * (aHi - aLo);
  const double roLo = std::max(rLo, nrLo), roHi = std::min(rHi, nrHi);
  const double doLo = std::max(dLo, ndLo), doHi = std::min(dHi, ndHi);
  const double aoLo = std::max(aLo, naLo), aoHi = std::min(aHi, naHi);
  const double Lr = roHi - roLo, Ld = doHi - doLo, La = aoHi - aoLo;

  if (Lr < -tolR || Ld < -tolD || La < -tolA) return 0.;
  if (Lr <= tolR && Ld > tolD && La > tolA) {
    const double R = std::abs(rHi - nrLo) < tolR ? rHi : rLo;
    return R * R * (otswap::maths::sin(doHi) - otswap::maths::sin(doLo)) * La;
  }
  if (Ld <= tolD && Lr > tolR && La > tolA) {
    const double D = std::abs(dHi - ndLo) < tolD ? dHi : dLo;
    return 0.5 * (roHi * roHi - roLo * roLo) * otswap::maths::cos(D) * La;
  }
  if (La <= tolA && Lr > tolR && Ld > tolD) return 0.5 * (roHi * roHi - roLo * roLo) * Ld;
  return 0.;
}


// ============================================================================


void bitvf::internal::lightcone_basins (const LightconeGrid& grid, const LightconeAdjacency& adjacency,
                                        LightconeBasins& basins, BasinCatalog& voids, BasinCatalog& overdense)
{
  const DivergenceGrid& cells = grid.cells;
  const std::int64_t C = (std::int64_t)cells.nCells;
  const double* D = (cells.divergenceSmoothed.empty() ? cells.divergence : cells.divergenceSmoothed).data();
  const std::int64_t* voidIndex = cells.voidIndex.data();
  const std::int64_t* overdenseIndex = cells.overdenseIndex.data();
  const bool withTracers = !cells.nTracers.empty();
  const std::size_t nVoid = basins.voids.seed.size();
  const std::size_t nAll = nVoid + basins.overdense.seed.size();
  const std::size_t nShells = grid.shellCentre.size();
  const std::size_t nRings = grid.ringDecMin.size();

  std::vector<std::size_t> ringShell(nRings);
  for (std::size_t i = 0; i < nShells; ++i)
    for (std::int64_t r = grid.shellFirstRing[i]; r < grid.shellFirstRing[i + 1]; ++r) ringShell[(std::size_t)r] = i;

  std::vector<std::int64_t> first(nAll + 1, 0);
  for (std::int64_t c = 0; c < C; ++c) {
    if (cells.masked[(std::size_t)c]) continue;
    const std::int64_t g = voidIndex[c] >= 0 ? voidIndex[c] : (std::int64_t)nVoid + overdenseIndex[c];
    ++first[(std::size_t)g + 1];
  }
  for (std::size_t g = 0; g < nAll; ++g) first[g + 1] += first[g];
  std::vector<std::int64_t> member((std::size_t)first[nAll]);
  {
    std::vector<std::int64_t> next(first.begin(), first.end() - 1);
    for (std::int64_t c = 0; c < C; ++c) {
      if (cells.masked[(std::size_t)c]) continue;
      const std::int64_t g = voidIndex[c] >= 0 ? voidIndex[c] : (std::int64_t)nVoid + overdenseIndex[c];
      member[(std::size_t)next[(std::size_t)g]++] = c;
    }
  }

  BasinSet* sets[2] = {&basins.voids, &basins.overdense};
  std::vector<double>* coreVolume[2] = {&basins.voidCoreVolume, &basins.overdenseCoreVolume};
  std::vector<double>* volume[2] = {&basins.voidVolume, &basins.overdenseVolume};
  BasinCatalog* catalogues[2] = {&voids, &overdense};
  for (int kind = 0; kind < 2; ++kind) {
    const std::size_t B = sets[kind]->seed.size();
    sets[kind]->nCells.assign(B, 0);
    sets[kind]->nCoreCells.assign(B, 0);
    sets[kind]->nCoreTracers.assign(withTracers ? B : 0, 0);
    coreVolume[kind]->assign(B, 0.);
    volume[kind]->assign(B, 0.);
    BasinCatalog& cat = *catalogues[kind];
    cat.nBasins = B;
    cat.centre.clear();
    cat.centreSky.clear();
    cat.deltaAv.clear();
    cat.delta2Mps.clear();
    cat.radius.resize(B);
    cat.radiusExt.resize(B);
    cat.edgeFlag.resize(B);
    cat.contamination.resize(B);
    cat.nTracers.assign(withTracers ? B : 0, 0);
  }

#pragma omp parallel for schedule(dynamic)
  for (std::int64_t g = 0; g < (std::int64_t)nAll; ++g) {
    const bool isVoid = (std::size_t)g < nVoid;
    const int kind = isVoid ? 0 : 1;
    const std::size_t b = isVoid ? (std::size_t)g : (std::size_t)g - nVoid;
    std::uint64_t nCells = 0, nCore = 0, nCoreTracers = 0;
    double vCore = 0., vExt = 0., lost = 0., phys = 0., total = 0.;
    bool touch = false, touchCore = false;

    for (std::int64_t m = first[(std::size_t)g]; m < first[(std::size_t)g + 1]; ++m) {
      const std::int64_t c = member[(std::size_t)m];
      const std::size_t cu = (std::size_t)c;
      const bool core = isVoid ? D[c] < 0. : D[c] >= 0.;
      const std::size_t ring = (std::size_t)(std::upper_bound(grid.ringFirstCell.begin(), grid.ringFirstCell.end(), c) -
                                             grid.ringFirstCell.begin()) - 1;
      const std::size_t shell = ringShell[ring];
      const std::int64_t ringFirst = grid.ringFirstCell[ring], ringEnd = grid.ringFirstCell[ring + 1];
      const double decLo = cells.cellDecMin[cu], decHi = cells.cellDecMax[cu];

      bool edge = (shell == 0 && !grid.centralSphere) || shell + 1 == nShells ||
                  ((std::int64_t)ring == grid.shellFirstRing[shell] && decLo != -kHalfPi) ||
                  ((std::int64_t)ring + 1 == grid.shellFirstRing[shell + 1] && decHi != kHalfPi) ||
                  (!grid.ringFull[ring] && (c == ringFirst || c + 1 == ringEnd));

      const double rIn = cells.cellDistanceMin[cu], rOut = cells.cellDistanceMax[cu];
      const double width = cells.cellRaMax[cu] - cells.cellRaMin[cu];
      const double sinBand = otswap::maths::sin(decHi) - otswap::maths::sin(decLo);
      const double half = 0.5 * (rOut * rOut - rIn * rIn);
      double area = rIn * rIn * sinBand * width;
      area += rOut * rOut * sinBand * width;
      if (decLo != -kHalfPi) area += half * otswap::maths::cos(decLo) * width;
      if (decHi != kHalfPi) area += half * otswap::maths::cos(decHi) * width;
      if (!(grid.ringFull[ring] && ringEnd - ringFirst == 1)) area += 2. * half * (decHi - decLo);

      double aInt = 0., aMsk = 0., aPhy = 0.;
      for (std::int64_t k = adjacency.contactFirst[cu]; k < adjacency.contactFirst[cu + 1]; ++k) {
        const std::int64_t n = adjacency.contactNeighbour[(std::size_t)k];
        const double a = lightcone_shared_area(cells, cu, (std::size_t)n);
        if (a <= 0.) continue;
        if (cells.masked[(std::size_t)n]) {
          aMsk += a;
          edge = true;
        }
        else if (isVoid ? voidIndex[n] == (std::int64_t)b : overdenseIndex[n] == (std::int64_t)b) aInt += a;
        else aPhy += a;
      }
      lost += aMsk + (area - aInt - aMsk - aPhy);
      phys += aPhy;
      total += area;

      if (edge) {
        touch = true;
        if (core) touchCore = true;
      }
      ++nCells;
      vExt += cells.volume[cu];
      if (core) {
        ++nCore;
        vCore += cells.volume[cu];
        if (withTracers) nCoreTracers += cells.nTracers[cu];
      }
    }

    BasinSet& set = *sets[kind];
    set.nCells[b] = nCells;
    set.nCoreCells[b] = nCore;
    if (withTracers) set.nCoreTracers[b] = nCoreTracers;
    (*coreVolume[kind])[b] = vCore;
    (*volume[kind])[b] = vExt;

    BasinCatalog& cat = *catalogues[kind];
    cat.radius[b] = otswap::maths::cbrt(3. * vCore / (4. * kPi));
    cat.radiusExt[b] = otswap::maths::cbrt(3. * vExt / (4. * kPi));
    cat.edgeFlag[b] = touchCore ? 2 : touch ? 1 : 0;
    const double bound = lost + phys;
    const double snap = 1024. * std::numeric_limits<double>::epsilon() * total;
    const double shown = lost < 0. && lost > -snap ? 0. : lost;
    cat.contamination[b] = bound > 0. ? shown / bound : 0.;
    if (withTracers) cat.nTracers[b] = (std::size_t)nCoreTracers;
  }
}


// ============================================================================
