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
 *  @file src/Watershed.cpp
 *
 *  @brief The watershed of the divergence of a box search into voids and
 *  overdense basins.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <initializer_list>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "internal.h"

namespace {

  /// The rows of a block of the seed search: blocks of fixed size, so that
  /// the seeds are gathered in the same order with any number of threads.
  constexpr std::int64_t kSeedBlockRows = std::int64_t(1) << 16;

}


// ============================================================================


void bitvf::internal::box_watershed (DivergenceGrid& grid, BoxBasins& basins)
{
  const std::int64_t n = (std::int64_t)grid.nCellsAxis[0];
  const std::int64_t C = (std::int64_t)grid.nCells;
  const double* D = grid.divergence.data();

  grid.voidIndex.resize((std::size_t)C);
  grid.overdenseIndex.resize((std::size_t)C);
  std::int64_t* label = grid.voidIndex.data();

  std::unique_ptr<std::int64_t[]> flow(new std::int64_t[(std::size_t)C]);

#pragma omp parallel for schedule(static)
  for (std::int64_t line = 0; line < n * n; ++line) {
    const std::int64_t i = line / n;
    const std::int64_t j = line % n;
    const std::int64_t i0 = std::max<std::int64_t>(0, i - 1), i1 = std::min(n - 1, i + 1);
    const std::int64_t j0 = std::max<std::int64_t>(0, j - 1), j1 = std::min(n - 1, j + 1);
    for (std::int64_t k = 0; k < n; ++k) {
      const std::int64_t c = line * n + k;
      const std::int64_t k0 = std::max<std::int64_t>(0, k - 1), k1 = std::min(n - 1, k + 1);
      double minVal = D[c];
      std::int64_t target = c;
      for (std::int64_t ii = i0; ii <= i1; ++ii)
        for (std::int64_t jj = j0; jj <= j1; ++jj)
          for (std::int64_t kk = k0; kk <= k1; ++kk) {
            const std::int64_t nb = (ii * n + jj) * n + kk;
            if (nb == c) continue;
            if (D[nb] < minVal) {
              minVal = D[nb];
              target = nb;
            }
          }
      flow[c] = target;
      label[c] = -1;
    }
  }

  const std::int64_t nBlocks = (C + kSeedBlockRows - 1) / kSeedBlockRows;
  std::vector<std::int64_t> voidOffset((std::size_t)nBlocks + 1, 0), overdenseOffset((std::size_t)nBlocks + 1, 0);

#pragma omp parallel for schedule(static)
  for (std::int64_t b = 0; b < nBlocks; ++b) {
    const std::int64_t end = std::min(C, (b + 1) * kSeedBlockRows);
    std::int64_t nVoid = 0, nOverdense = 0;
    for (std::int64_t c = b * kSeedBlockRows; c < end; ++c)
      if (flow[c] == c) {
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

  if (voidOffset[(std::size_t)nBlocks] == 0)
    throw Error("no void found: the divergence is 0 or positive in every cell");

  basins.voids.seed.resize((std::size_t)voidOffset[(std::size_t)nBlocks]);
  basins.overdense.seed.resize((std::size_t)overdenseOffset[(std::size_t)nBlocks]);
  std::int64_t* voidSeed = basins.voids.seed.data();
  std::int64_t* overdenseSeed = basins.overdense.seed.data();

#pragma omp parallel for schedule(static)
  for (std::int64_t b = 0; b < nBlocks; ++b) {
    const std::int64_t end = std::min(C, (b + 1) * kSeedBlockRows);
    std::int64_t v = voidOffset[(std::size_t)b], o = overdenseOffset[(std::size_t)b];
    for (std::int64_t c = b * kSeedBlockRows; c < end; ++c)
      if (flow[c] == c) {
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
  std::vector<std::int64_t>().swap(paths);
  flow.reset();

  for (BasinSet* set : {&basins.voids, &basins.overdense}) {
    set->nCells.assign(set->seed.size(), 0);
    set->nCoreCells.assign(set->seed.size(), 0);
    set->nCoreTracers.assign(set->seed.size(), 0);
  }
  std::uint64_t* voidCells = basins.voids.nCells.data();
  std::uint64_t* voidCoreCells = basins.voids.nCoreCells.data();
  std::uint64_t* voidCoreTracers = basins.voids.nCoreTracers.data();
  std::uint64_t* overdenseCells = basins.overdense.nCells.data();
  std::uint64_t* overdenseCoreCells = basins.overdense.nCoreCells.data();
  std::uint64_t* overdenseCoreTracers = basins.overdense.nCoreTracers.data();
  const std::uint64_t* nTracers = grid.nTracers.data();
  std::int64_t* overdense = grid.overdenseIndex.data();

#pragma omp parallel for schedule(static)
  for (std::int64_t c = 0; c < C; ++c) {
    const std::int64_t L = label[c];
    if (L >= 0) {
      overdense[c] = -1;
#pragma omp atomic
      voidCells[L] += 1;
      if (D[c] < 0.) {
#pragma omp atomic
        voidCoreCells[L] += 1;
#pragma omp atomic
        voidCoreTracers[L] += nTracers[c];
      }
    }
    else {
      const std::int64_t b = -2 - L;
      overdense[c] = b;
      label[c] = -1;
#pragma omp atomic
      overdenseCells[b] += 1;
      if (D[c] >= 0.) {
#pragma omp atomic
        overdenseCoreCells[b] += 1;
#pragma omp atomic
        overdenseCoreTracers[b] += nTracers[c];
      }
    }
  }
}
