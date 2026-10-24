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

  const std::size_t nVoid = basins.voids.seed.size();
  const std::size_t nOverdense = basins.overdense.seed.size();
  const std::size_t stride = 3 * (nVoid + nOverdense);
  std::unique_ptr<std::uint64_t[]> counts(new std::uint64_t[nThreads * stride]);
  std::size_t nTeam = 1;
  const std::uint64_t* nTracers = grid.nTracers.data();
  std::int64_t* overdense = grid.overdenseIndex.data();

#pragma omp parallel
  {
#ifdef _OPENMP
    std::uint64_t* voidCells = &counts[(std::size_t)omp_get_thread_num() * stride];
#pragma omp single nowait
    nTeam = (std::size_t)omp_get_num_threads();
#else
    std::uint64_t* voidCells = counts.get();
#endif
    std::fill(voidCells, voidCells + stride, std::uint64_t(0));
    std::uint64_t* voidCoreCells = voidCells + nVoid;
    std::uint64_t* voidCoreTracers = voidCoreCells + nVoid;
    std::uint64_t* overdenseCells = voidCoreTracers + nVoid;
    std::uint64_t* overdenseCoreCells = overdenseCells + nOverdense;
    std::uint64_t* overdenseCoreTracers = overdenseCoreCells + nOverdense;

#pragma omp for schedule(static)
    for (std::int64_t c = 0; c < C; ++c) {
      const std::int64_t L = label[c];
      if (L >= 0) {
        overdense[c] = -1;
        voidCells[L] += 1;
        if (D[c] < 0.) {
          voidCoreCells[L] += 1;
          voidCoreTracers[L] += nTracers[c];
        }
      }
      else {
        const std::int64_t b = -2 - L;
        overdense[c] = b;
        label[c] = -1;
        overdenseCells[b] += 1;
        if (D[c] >= 0.) {
          overdenseCoreCells[b] += 1;
          overdenseCoreTracers[b] += nTracers[c];
        }
      }
    }
  }

  std::size_t offset = 0;
  for (BasinSet* set : {&basins.voids, &basins.overdense})
    for (std::vector<std::uint64_t>* count : {&set->nCells, &set->nCoreCells, &set->nCoreTracers}) {
      const std::int64_t m = (std::int64_t)set->seed.size();
      count->resize((std::size_t)m);
      std::uint64_t* sum = count->data();
#pragma omp parallel for schedule(static)
      for (std::int64_t b = 0; b < m; ++b) {
        std::uint64_t s = 0;
        for (std::size_t t = 0; t < nTeam; ++t)
          s += counts[t * stride + offset + (std::size_t)b];
        sum[b] = s;
      }
      offset += (std::size_t)m;
    }
}
