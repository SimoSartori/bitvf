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
 *  @file src/Divergence.cpp
 *
 *  @brief The divergence of a box search: the face crossed by a
 *  displacement, the fixed-point accumulation, the fill of the empty
 *  cells.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "internal.h"

namespace {

  /// The base-2 logarithm of the largest count of a cell the scale of the
  /// accumulation assumes; a larger count reruns the accumulation at the
  /// scale of the real count.
  constexpr int kAssumedCountBits = 20;

  /// The largest exponent of the scale, so that 2^-q times a sum stays a
  /// normal double.
  constexpr int kMaxScaleExponent = 1000;

  /// No offending displacement.
  constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

  /// The first offending displacements of a sweep, as tracer * R +
  /// realization; kNone for none.
  struct SweepErrors {
    std::size_t outside = kNone;
    std::size_t noFace = kNone;
  };

  /**
   *  @brief The scale exponent q = 61 - e - L, with B < 2^e, capped at
   *  kMaxScaleExponent; 0 when B is 0.
   *  @param B the largest component of a displacement, in cell sides
   *  @param L the base-2 logarithm of the largest count of a cell
   */
  int scale_exponent (const double B, const int L)
  {
    if (!(B > 0.)) return 0;
    int e = 0;
    std::frexp(B, &e);
    return std::min(61 - e - L, kMaxScaleExponent);
  }

  /**
   *  @brief One sweep over the tracers read: start and end cells, faces,
   *  and the terms of the divergence added to numerator as multiples of
   *  2^-q, modulo 2^64.
   *  @param count whether to also count DivergenceGrid::nTracers and
   *  nDisplacements, and fill result.cell and result.used
   */
  SweepErrors sweep (const bitvf::internal::BoxSearch& search, bitvf::DivergenceGrid& grid,
                     bitvf::internal::BoxDivergence& result, std::vector<std::uint64_t>& numerator,
                     const int q, const bool count)
  {
    const std::size_t N = search.nTracers;
    const std::size_t R = search.nRealizations;
    const std::vector<double>& tracers = *search.tracers;
    const std::vector<double>& displacement = *search.displacement;
    const std::vector<std::uint8_t>& valid = *search.valid;
    const std::size_t n1 = grid.nCellsAxis[1];
    const std::size_t n2 = grid.nCellsAxis[2];
    const double inverse = 1. / grid.cellSize;

    std::size_t outside = kNone;
    std::size_t noFace = kNone;

#pragma omp parallel for schedule(static) reduction(min:outside,noFace)
    for (std::size_t i = 0; i < N; ++i) {
      if (!search.inBox[i] || (!search.outsideGrid.empty() && search.outsideGrid[i])) continue;

      const double* x = &tracers[3*i];
      const std::int64_t c0 = bitvf::internal::box_cell(grid, x);
      if (count) {
        result.cell[i] = c0;
#pragma omp atomic
        grid.nTracers[c0] += 1;
      }

      for (std::size_t r = 0; r < R; ++r) {
        const std::size_t row = r * N + i;
        if (!valid.empty() && !valid[row]) continue;

        const double* d = &displacement[3*row];
        const double end[3] = {x[0] + d[0], x[1] + d[1], x[2] + d[2]};
        const std::int64_t c1 = bitvf::internal::box_cell(grid, end);
        if (c1 < 0) {
          outside = std::min(outside, i * R + r);
          continue;
        }
        if (count) result.used[i] = 1;

        if (c1 == c0) {
          if (count) {
#pragma omp atomic
            grid.nDisplacements[c0] += 2;
          }
          continue;
        }

        for (const std::int64_t c : {c0, c1}) {
          const std::size_t u = (std::size_t)c;
          const bitvf::internal::BoxFace face =
            bitvf::internal::box_face(grid, {u / (n1 * n2), (u / n2) % n1, u % n2}, x, d);
          if (face.axis < 0) {
            noFace = std::min(noFace, i * R + r);
            continue;
          }
          if (count) {
#pragma omp atomic
            grid.nDisplacements[c] += 1;
          }
          const double t = (face.sign * d[face.axis]) * inverse;
          const std::uint64_t term = (std::uint64_t)(std::int64_t)std::llrint(std::ldexp(t, q));
#pragma omp atomic
          numerator[c] += term;
        }
      }
    }

    return {outside, noFace};
  }

}


// ============================================================================


bitvf::internal::BoxFace bitvf::internal::box_face (const DivergenceGrid& grid,
                                                    const std::array<std::size_t, 3>& cell,
                                                    const double* x, const double* d)
{
  const std::array<double, 3>& lo = grid.lower;
  const double s = grid.cellSize;

  for (std::size_t fac = 0; fac < 2; ++fac)
    for (int a = 0; a < 3; ++a) {
      if (!(std::abs(d[a]) > 1e-9)) continue;
      const int b = (a == 0) ? 1 : 0;
      const int c = (a == 2) ? 1 : 2;
      const double t = (lo[a] + s * (double)(cell[a] + fac) - x[a]) / d[a];
      const double yb = x[b] + t * d[b];
      const double yc = x[c] + t * d[c];
      if (t > 0 && t < 1 &&
          yb > (lo[b] + s * (double)cell[b]) && yb < (lo[b] + s * (double)(cell[b] + 1)) &&
          yc > (lo[c] + s * (double)cell[c]) && yc < (lo[c] + s * (double)(cell[c] + 1)))
        return {a, -1 + 2 * (int)fac};
    }

  return {-1, 0};
}


// ============================================================================


bitvf::internal::BoxDivergence bitvf::internal::box_divergence (const BoxSearch& search,
                                                                DivergenceGrid& grid)
{
  const std::size_t N = search.nTracers;
  const std::size_t R = search.nRealizations;
  const std::size_t C = grid.nCells;

  grid.nTracers.assign(C, 0);
  grid.nDisplacements.assign(C, 0);

  BoxDivergence result;
  result.cell.assign(N, -1);
  result.used.assign(N, 0);

  const double B = search.maxDisplacement * (1. / grid.cellSize);
  int q = scale_exponent(B, kAssumedCountBits);

  std::vector<std::uint64_t> numerator(C, 0);
  const SweepErrors errors = sweep(search, grid, result, numerator, q, true);

  if (errors.outside != kNone)
    throw Error("the valid displacement of tracer " + std::to_string(errors.outside / R) +
                " in realization " + std::to_string(errors.outside % R) + " ends outside every "
                "cell of the grid, which spans the box limits (with an otswap result, its box) or "
                "the reused grid; the tracers and the displacement do not come from one "
                "reconstruction in that box: pass the box limits the reconstruction used, or none");
  if (errors.noFace != kNone)
    throw Error("internal error: no face of its cell is crossed by the displacement of tracer " +
                std::to_string(errors.noFace / R) + " in realization " +
                std::to_string(errors.noFace % R) + ", which leaves that cell");

  std::uint64_t K = 0;
#pragma omp parallel for schedule(static) reduction(max:K)
  for (std::size_t c = 0; c < C; ++c)
    K = std::max(K, grid.nDisplacements[c]);

  if (K > (std::uint64_t{1} << kAssumedCountBits)) {
    int L = kAssumedCountBits;
    while ((std::uint64_t{1} << L) < K) ++L;
    q = scale_exponent(B, L);
    std::fill(numerator.begin(), numerator.end(), 0);
    sweep(search, grid, result, numerator, q, false);
  }
  result.scaleExponent = q;

  grid.divergence.assign(C, 0.);
  grid.weight.assign(C, 0.);
#pragma omp parallel for schedule(static)
  for (std::size_t c = 0; c < C; ++c) {
    const std::uint64_t n = grid.nDisplacements[c];
    grid.divergence[c] = n > 0 ? std::ldexp((double)(std::int64_t)numerator[c], -q) / (double)n : 0.;
    grid.weight[c] = (double)n / (double)R;
  }
  std::vector<std::uint64_t>().swap(numerator);

  std::size_t nUsed = 0;
#pragma omp parallel for schedule(static) reduction(+:nUsed)
  for (std::size_t i = 0; i < N; ++i)
    nUsed += result.used[i];
  result.nUsed = nUsed;

  result.fillPasses = box_fill(grid);

  return result;
}


// ============================================================================


unsigned bitvf::internal::box_fill (DivergenceGrid& grid)
{
  const std::int64_t n0 = (std::int64_t)grid.nCellsAxis[0];
  const std::int64_t n1 = (std::int64_t)grid.nCellsAxis[1];
  const std::int64_t n2 = (std::int64_t)grid.nCellsAxis[2];
  const double s = grid.cellSize;

  std::vector<std::size_t> empty;
  for (std::size_t c = 0; c < grid.nCells; ++c)
    if (grid.nDisplacements[c] == 0) {
      empty.push_back(c);
      grid.divergence[c] = std::numeric_limits<double>::quiet_NaN();
    }

  std::vector<double> value(empty.size());
  std::vector<double> weight(empty.size());

  unsigned passes = 0;
  while (!empty.empty()) {
    ++passes;
    const std::size_t E = empty.size();
    std::size_t nFilled = 0;

#pragma omp parallel for schedule(dynamic, 64) reduction(+:nFilled)
    for (std::size_t e = 0; e < E; ++e) {
      const std::int64_t c = (std::int64_t)empty[e];
      const std::int64_t i = c / (n1 * n2);
      const std::int64_t j = (c / n2) % n1;
      const std::int64_t k = c % n2;

      value[e] = std::numeric_limits<double>::quiet_NaN();
      for (std::int64_t h = 1; h <= 2 && std::isnan(value[e]); ++h) {
        double num = 0., den = 0., w = 0.;
        unsigned counterValues = 0;
        for (std::int64_t ii = std::max<std::int64_t>(0, i - h); ii < std::min(n0, i + h + 1); ++ii)
          for (std::int64_t jj = std::max<std::int64_t>(0, j - h); jj < std::min(n1, j + h + 1); ++jj)
            for (std::int64_t kk = std::max<std::int64_t>(0, k - h); kk < std::min(n2, k + h + 1); ++kk) {
              const std::size_t cc = (std::size_t)((ii * n1 + jj) * n2 + kk);
              if (std::isnan(grid.divergence[cc])) continue;
              ++counterValues;
              const double inv_d2 =
                1. / ((double)((i - ii) * (i - ii) + (j - jj) * (j - jj) + (k - kk) * (k - kk)) * s * s);
              num += grid.divergence[cc] * inv_d2;
              den += inv_d2;
              w += grid.weight[cc] * inv_d2;
            }
        if (!(den == 0 || counterValues < 3)) {
          value[e] = num / den;
          weight[e] = w / den;
          ++nFilled;
        }
      }
    }

    if (nFilled == 0) {
      for (const std::size_t c : empty) grid.divergence[c] = 0.;
      throw Error("the fill of the divergence cannot reach " + std::to_string(E) +
                  " empty cells: the last pass filled none");
    }

    std::size_t kept = 0;
    for (std::size_t e = 0; e < E; ++e)
      if (!std::isnan(value[e])) {
        grid.divergence[empty[e]] = value[e];
        grid.weight[empty[e]] = weight[e];
      }
      else
        empty[kept++] = empty[e];
    empty.resize(kept);
  }

  return passes;
}


// ============================================================================


bitvf::internal::BoxDivergence bitvf::internal::box_tracer_cells (const BoxSearch& search, DivergenceGrid& grid)
{
  const std::size_t N = search.nTracers;
  const std::vector<double>& tracers = *search.tracers;

  grid.nTracers.assign(grid.nCells, 0);

  BoxDivergence result;
  result.cell.assign(N, -1);

#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < N; ++i) {
    if (!search.inBox[i] || (!search.outsideGrid.empty() && search.outsideGrid[i])) continue;
    const std::int64_t c = bitvf::internal::box_cell(grid, &tracers[3*i]);
    result.cell[i] = c;
#pragma omp atomic
    grid.nTracers[c] += 1;
  }

  return result;
}
