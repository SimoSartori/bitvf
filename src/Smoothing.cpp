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
 *  @file src/Smoothing.cpp
 *
 *  @brief The Gaussian smoothing of the divergence of a box search.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <initializer_list>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  /// The number of lines of a tile: the lines a thread convolves together.
  constexpr std::size_t kTileLines = 64;

  /**
   *  @brief Replace num and den by their 1-D convolution with K along one
   *  axis of a cubic grid of n cells per axis, cut at the faces.
   *
   *  Each entry becomes the sum, from 0. and over the cells x' of its line
   *  within h of it in ascending order, of K[x' - x + h] times the entry at
   *  x'. The lines are convolved in tiles of up to kTileLines adjacent
   *  lines, each tile by one thread; the result does not depend on the
   *  tiles or the threads.
   *
   *  @param axis 0, 1 or 2: the axis of the row index (i, j, k), k the
   *  contiguous one
   *  @param buffer 4 n kTileLines entries per thread of the parallel region
   */
  void convolve_axis (std::vector<double>& num, std::vector<double>& den, const std::size_t n,
                      const int axis, const std::vector<double>& K, const std::size_t h,
                      std::vector<double>& buffer)
  {
    const std::size_t nOuter = axis == 0 ? 1 : n;
    const std::size_t nInner = axis == 0 ? n * n : n;
    const std::size_t outerStride = n * n;
    const std::size_t lineStride = axis == 2 ? n : 1;
    const std::size_t stride = axis == 0 ? n * n : (axis == 1 ? n : 1);
    const std::size_t nBlocks = (nInner + kTileLines - 1) / kTileLines;
    const std::size_t nTiles = nOuter * nBlocks;
    const std::size_t tileSize = n * kTileLines;

#pragma omp parallel
    {
#ifdef _OPENMP
      const std::size_t thread = (std::size_t)omp_get_thread_num();
#else
      const std::size_t thread = 0;
#endif
      double* inNum = &buffer[thread * 4 * tileSize];
      double* inDen = inNum + tileSize;
      double* outNum = inDen + tileSize;
      double* outDen = outNum + tileSize;

#pragma omp for schedule(static)
      for (std::size_t t = 0; t < nTiles; ++t) {
        const std::size_t p0 = (t % nBlocks) * kTileLines;
        const std::size_t B = std::min(kTileLines, nInner - p0);
        const std::size_t base = (t / nBlocks) * outerStride + p0 * lineStride;

        for (std::size_t x = 0; x < n; ++x)
          for (std::size_t b = 0; b < B; ++b) {
            inNum[x*B+b] = num[base + x * stride + b * lineStride];
            inDen[x*B+b] = den[base + x * stride + b * lineStride];
          }

        for (std::size_t x = 0; x < n; ++x) {
          double* oNum = outNum + x * B;
          double* oDen = outDen + x * B;
          for (std::size_t b = 0; b < B; ++b) {
            oNum[b] = 0.;
            oDen[b] = 0.;
          }
          const std::size_t hi = std::min(n - 1, x + h);
          for (std::size_t xx = (x >= h ? x - h : 0); xx <= hi; ++xx) {
            const double k = K[xx + h - x];
            const double* iNum = inNum + xx * B;
            const double* iDen = inDen + xx * B;
            for (std::size_t b = 0; b < B; ++b) {
              oNum[b] += k * iNum[b];
              oDen[b] += k * iDen[b];
            }
          }
        }

        for (std::size_t x = 0; x < n; ++x)
          for (std::size_t b = 0; b < B; ++b) {
            num[base + x * stride + b * lineStride] = outNum[x*B+b];
            den[base + x * stride + b * lineStride] = outDen[x*B+b];
          }
      }
    }
  }

}


// ============================================================================


void bitvf::internal::box_smooth (DivergenceGrid& grid, const double sigmaLength,
                                  std::vector<double>* unsmoothed)
{
  if (sigmaLength == 0.) {
    if (unsmoothed != nullptr) *unsmoothed = grid.divergence;
    return;
  }

  const std::size_t n = grid.nCellsAxis[0];
  const std::size_t C = grid.nCells;
  const double s = grid.cellSize;

  const double hd = 3.5 * sigmaLength / s;
  const std::size_t h = (std::size_t)std::min(hd, (double)(n - 1));
  const double two = 2 * sigmaLength * sigmaLength;
  std::vector<double> K(2 * h + 1);
  for (std::size_t i = 0; i < K.size(); ++i) {
    const double dist = (double)((long long)i - (long long)h) * s;
    K[i] = otswap::maths::exp(-(dist * dist) / two);
  }

  if (unsmoothed != nullptr) *unsmoothed = std::move(grid.divergence);
  const std::vector<double>& field = unsmoothed != nullptr ? *unsmoothed : grid.divergence;

  std::vector<double> num(C), den(C);
#pragma omp parallel for schedule(static)
  for (std::size_t c = 0; c < C; ++c) {
    num[c] = grid.weight[c] * field[c];
    den[c] = grid.weight[c];
  }

#ifdef _OPENMP
  const std::size_t nThreads = (std::size_t)omp_get_max_threads();
#else
  const std::size_t nThreads = 1;
#endif
  std::vector<double> buffer(nThreads * 4 * n * kTileLines);
  for (const int axis : {2, 1, 0})
    convolve_axis(num, den, n, axis, K, h, buffer);
  std::vector<double>().swap(buffer);

#pragma omp parallel for schedule(static)
  for (std::size_t c = 0; c < C; ++c)
    num[c] = den[c] > 1e-12 ? num[c] / den[c] : 0.;

  grid.divergence = std::move(num);
}
