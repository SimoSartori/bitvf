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
#include <array>
#include <cmath>
#include <initializer_list>
#include <memory>
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
   *  axis of a grid of nCellsAxis cells, cut at the faces.
   *
   *  Each entry becomes the sum, from 0. and over the cells x' of its line
   *  within h of it in ascending order, of K[x' - x + hK] times the entry
   *  at x'. The lines are convolved in tiles of up to kTileLines adjacent
   *  lines, each tile by one thread; the result does not depend on the
   *  tiles or the threads.
   *
   *  @param axis 0, 1 or 2: the axis of the row index (i, j, k), k the
   *  contiguous one
   *  @param h the half-width along axis, at most hK
   *  @param hK the half-width of K, which holds 2 hK + 1 entries
   *  @param buffer 4 max(nCellsAxis) kTileLines entries per thread of the
   *  parallel region
   */
  void convolve_axis (double* num, double* den, const std::array<std::size_t, 3>& nCellsAxis,
                      const int axis, const std::vector<double>& K, const std::size_t h, const std::size_t hK,
                      std::vector<double>& buffer)
  {
    const std::size_t n = nCellsAxis[(std::size_t)axis];
    const std::size_t ny = nCellsAxis[1], nz = nCellsAxis[2];
    const std::size_t nOuter = axis == 0 ? 1 : nCellsAxis[0];
    const std::size_t nInner = axis == 0 ? ny * nz : (axis == 1 ? nz : ny);
    const std::size_t outerStride = ny * nz;
    const std::size_t lineStride = axis == 2 ? nz : 1;
    const std::size_t stride = axis == 0 ? ny * nz : (axis == 1 ? nz : 1);
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

#pragma omp for schedule(dynamic)
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
            const double k = K[xx + hK - x];
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


void bitvf::internal::box_smooth (DivergenceGrid& grid, const double sigmaLength, const double threshold)
{
  if (sigmaLength == 0.) {
    grid.divergenceSmoothed.clear();
    grid.weightSmoothed.clear();
    return;
  }

  const std::array<std::size_t, 3> n = grid.nCellsAxis;
  const std::size_t C = grid.nCells;
  const double s = grid.cellSize;

  double hd = 3.5 * sigmaLength / s;
  const double hn = std::floor(hd + 0.5);
  if (std::abs(hd - hn) <= 1e-9 * hn) hd = hn;
  std::array<std::size_t, 3> h;
  for (std::size_t a = 0; a < 3; ++a)
    h[a] = (std::size_t)std::min(hd, (double)(n[a] - 1));
  const std::size_t hK = *std::max_element(h.begin(), h.end());
  const double two = 2 * sigmaLength * sigmaLength;
  std::vector<double> K(2 * hK + 1);
  for (std::size_t i = 0; i < K.size(); ++i) {
    const double dist = (double)((long long)i - (long long)hK) * s;
    K[i] = otswap::maths::exp(-(dist * dist) / two);
  }

  std::array<std::vector<double>, 3> kernelSum;
  for (std::size_t a = 0; a < 3; ++a) {
    kernelSum[a].assign(n[a], 0.);
    for (std::size_t x = 0; x < n[a]; ++x) {
      const std::size_t hi = std::min(n[a] - 1, x + h[a]);
      for (std::size_t xx = (x >= h[a] ? x - h[a] : 0); xx <= hi; ++xx)
        kernelSum[a][x] += K[xx + hK - x];
    }
  }

  const std::vector<double>& field = grid.divergence;

  std::unique_ptr<double[]> num(new double[C]), den(new double[C]);
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
  std::vector<double> buffer(nThreads * 4 * *std::max_element(n.begin(), n.end()) * kTileLines);
  for (const int axis : {2, 1, 0})
    convolve_axis(num.get(), den.get(), n, axis, K, h[(std::size_t)axis], hK, buffer);
  std::vector<double>().swap(buffer);

  grid.divergenceSmoothed.resize(C);
  grid.weightSmoothed.resize(C);
  const std::size_t nyz = n[1] * n[2];
#pragma omp parallel for schedule(static)
  for (std::size_t c = 0; c < C; ++c) {
    grid.divergenceSmoothed[c] = den[c] > threshold ? num[c] / den[c] : 0.;
    grid.weightSmoothed[c] =
      den[c] / (kernelSum[0][c / nyz] * kernelSum[1][(c / n[2]) % n[1]] * kernelSum[2][c % n[2]]);
  }
}
