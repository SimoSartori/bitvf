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
 *  @file src/LightconeWatershed.cpp
 *
 *  @brief The watershed of the divergence of a lightcone search into voids
 *  and overdense basins.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "internal.h"

namespace {

  constexpr std::int64_t kMaskedLabel = std::numeric_limits<std::int64_t>::min();

}


// ============================================================================


bitvf::internal::LightconeBasins bitvf::internal::lightcone_watershed (const LightconeAdjacency& adjacency,
                                                                       DivergenceGrid& grid)
{
  LightconeBasins basins;
  const std::int64_t C = (std::int64_t)grid.nCells;
  const double* D = (grid.divergenceSmoothed.empty() ? grid.divergence : grid.divergenceSmoothed).data();

  grid.voidIndex.assign((std::size_t)C, -1);
  grid.overdenseIndex.assign((std::size_t)C, -1);
  std::int64_t* label = grid.voidIndex.data();
  std::unique_ptr<std::int64_t[]> flow(new std::int64_t[(std::size_t)C]);

#pragma omp parallel for schedule(dynamic, 1024)
  for (std::int64_t c = 0; c < C; ++c) {
    flow[c] = c;
    if (grid.masked[(std::size_t)c]) {
      label[c] = kMaskedLabel;
      continue;
    }
    label[c] = -1;
    double minVal = D[c];
    for (std::int64_t k = adjacency.radiusFirst[(std::size_t)c]; k < adjacency.radiusFirst[(std::size_t)c + 1]; ++k) {
      const std::int64_t nb = adjacency.radiusNeighbour[(std::size_t)k];
      if (D[nb] < minVal) {
        minVal = D[nb];
        flow[c] = nb;
      }
    }
  }

  const std::int64_t nBlocks = (C + kSeedBlockRows - 1) / kSeedBlockRows;
  std::vector<std::int64_t> voidOffset((std::size_t)nBlocks + 1, 0), overdenseOffset((std::size_t)nBlocks + 1, 0);

#pragma omp parallel for schedule(static)
  for (std::int64_t b = 0; b < nBlocks; ++b) {
    const std::int64_t end = std::min(C, (b + 1) * kSeedBlockRows);
    std::int64_t nVoid = 0, nOverdense = 0;
    for (std::int64_t c = b * kSeedBlockRows; c < end; ++c)
      if (label[c] == -1 && flow[c] == c) {
        if (D[c] < 0.) ++nVoid;
        else ++nOverdense;
      }
    voidOffset[(std::size_t)b + 1] = nVoid;
    overdenseOffset[(std::size_t)b + 1] = nOverdense;
  }

  for (std::size_t b = 0; b < (std::size_t)nBlocks; ++b) {
    voidOffset[b + 1] += voidOffset[b];
    overdenseOffset[b + 1] += overdenseOffset[b];
  }

  basins.voids.seed.resize((std::size_t)voidOffset[(std::size_t)nBlocks]);
  basins.overdense.seed.resize((std::size_t)overdenseOffset[(std::size_t)nBlocks]);
  std::int64_t* voidSeed = basins.voids.seed.data();
  std::int64_t* overdenseSeed = basins.overdense.seed.data();

#pragma omp parallel for schedule(static)
  for (std::int64_t b = 0; b < nBlocks; ++b) {
    const std::int64_t end = std::min(C, (b + 1) * kSeedBlockRows);
    std::int64_t v = voidOffset[(std::size_t)b], o = overdenseOffset[(std::size_t)b];
    for (std::int64_t c = b * kSeedBlockRows; c < end; ++c)
      if (label[c] == -1 && flow[c] == c) {
        if (D[c] < 0.) {
          voidSeed[v] = c;
          label[c] = v++;
        }
        else {
          overdenseSeed[o] = c;
          label[c] = -2 - o++;
        }
      }
  }

#ifdef _OPENMP
  const std::size_t nThreads = (std::size_t)omp_get_max_threads();
#else
  const std::size_t nThreads = 1;
#endif
  std::vector<std::int64_t> paths(nThreads * kPathCells);

#pragma omp parallel
  {
#ifdef _OPENMP
    std::int64_t* path = &paths[(std::size_t)omp_get_thread_num() * kPathCells];
#else
    std::int64_t* path = paths.data();
#endif

#pragma omp for schedule(dynamic, 4096)
    for (std::int64_t c = 0; c < C; ++c) {
      std::int64_t L;
#pragma omp atomic read
      L = label[c];
      if (L != -1) continue;

      std::size_t length = 0;
      std::int64_t cur = c;
      while (true) {
#pragma omp atomic read
        L = label[cur];
        if (L != -1) break;
        if (length < kPathCells) path[length++] = cur;
        cur = flow[cur];
      }
      for (std::size_t p = 0; p < length; ++p) {
#pragma omp atomic write
        label[path[p]] = L;
      }
    }
  }

  std::int64_t* overdense = grid.overdenseIndex.data();
#pragma omp parallel for schedule(static)
  for (std::int64_t c = 0; c < C; ++c) {
    const std::int64_t L = label[c];
    if (L == kMaskedLabel) label[c] = -1;
    else if (L < 0) {
      overdense[c] = -2 - L;
      label[c] = -1;
    }
  }
  return basins;
}


// ============================================================================
