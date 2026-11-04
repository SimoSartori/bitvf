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
 *  @file src/BoxGrid.cpp
 *
 *  @brief The box grid: its construction, boxGrid, the cell of a point,
 *  the cell centres.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>

#include "otswap/Maths.h"

#include "internal.h"


// ============================================================================


bitvf::DivergenceGrid bitvf::internal::box_grid (const BoxSearch& search, const double cellSizeMps)
{
  std::array<double, 3> lo, hi;

  if (search.box) {
    lo = search.box->lower;
    hi = search.box->upper;
  }
  else {
    lo.fill(std::numeric_limits<double>::infinity());
    hi.fill(-std::numeric_limits<double>::infinity());

    const std::vector<double>& tracers = *search.tracers;
    const std::size_t N = search.nTracers;
    for (std::size_t i = 0; i < N; ++i)
      if (search.inBox[i])
        for (std::size_t k = 0; k < 3; ++k) {
          lo[k] = std::min(lo[k], tracers[3*i+k]);
          hi[k] = std::max(hi[k], tracers[3*i+k]);
        }

    const std::vector<double>& displacement = *search.displacement;
    const std::vector<std::uint8_t>& valid = *search.valid;
    for (std::size_t r = 0; r < search.nRealizations; ++r)
      for (std::size_t i = 0; i < N; ++i) {
        const std::size_t row = r * N + i;
        if (!search.inBox[i] || (!valid.empty() && !valid[row])) continue;
        for (std::size_t k = 0; k < 3; ++k) {
          const double end = tracers[3*i+k] + displacement[3*row+k];
          lo[k] = std::min(lo[k], end);
          hi[k] = std::max(hi[k], end);
        }
      }
  }

  const double mps = search.mps;

  DivergenceGrid grid;
  grid.geometry = otswap::Geometry::Box;
  grid.mps = mps;
  grid.cellSize = cellSizeMps * mps;
  grid.cellSizeMps = cellSizeMps;

  double count = 0.;
  for (std::size_t k = 0; k < 3; ++k) {
    grid.lower[k] = lo[k] - 0.1*mps;
    grid.upper[k] = hi[k] + 0.1*mps;
    count = std::max(count, (grid.upper[k] - grid.lower[k]) / grid.cellSize + 1);
  }

  if (!(count < (double)(kMaxCellsAxis + 1)))
    throw Error("the grid would need " + std::to_string(count) + " cells along each axis, more "
                "than the " + std::to_string(kMaxCellsAxis) + " an int64 cell index can address;"
                " increase the divergence cell size");

  const std::size_t n = (std::size_t)count;
  grid.nCellsAxis = {n, n, n};
  grid.nCells = n * n * n;

  for (std::size_t k = 0; k < 3; ++k) {
    grid.lower[k] = grid.lower[k] - (grid.cellSize * n - (grid.upper[k] - grid.lower[k])) / 2;
    grid.upper[k] = grid.upper[k] + (grid.cellSize * n - (grid.upper[k] - grid.lower[k])) / 2;
  }

  return grid;
}


// ============================================================================


std::int64_t bitvf::internal::box_cell (const DivergenceGrid& grid, const double* x)
{
  const double inverse = 1. / grid.cellSize;

  std::int64_t row = 0;
  for (std::size_t k = 0; k < 3; ++k) {
    const double t = (x[k] - grid.lower[k]) * inverse;
    if (!(t >= 0. && t < (double)grid.nCellsAxis[k])) return -1;
    row = row * (std::int64_t)grid.nCellsAxis[k] + (std::int64_t)t;
  }

  return row;
}


// ============================================================================


void bitvf::internal::box_centres (DivergenceGrid& grid)
{
  const std::array<std::size_t, 3>& n = grid.nCellsAxis;
  const double s = grid.cellSize;

  grid.centre.resize(3 * grid.nCells);
  for (std::size_t i = 0; i < n[0]; ++i)
    for (std::size_t j = 0; j < n[1]; ++j)
      for (std::size_t k = 0; k < n[2]; ++k) {
        const std::size_t c = (i * n[1] + j) * n[2] + k;
        grid.centre[3*c]   = grid.lower[0] + s * (i + 0.5);
        grid.centre[3*c+1] = grid.lower[1] + s * (j + 0.5);
        grid.centre[3*c+2] = grid.lower[2] + s * (k + 0.5);
      }
}


// ============================================================================


bitvf::DivergenceGrid bitvf::boxGrid (const otswap::BoxLimits& box, const double mps,
                                      const double divergenceCellSize)
{
  internal::check_limits(box, "the box limits");

  if (!std::isfinite(mps) || mps <= 0.)
    throw Error("the mean particle separation is " + std::to_string(mps) +
                " in coordinate units; it must be finite and positive");

  if (!std::isfinite(divergenceCellSize) || divergenceCellSize <= 0.)
    throw Error("the divergence cell size is " + std::to_string(divergenceCellSize) +
                " mean particle separations; it must be finite and positive");

  internal::BoxSearch search;
  search.box = box;
  search.mps = mps;
  return internal::box_grid(search, divergenceCellSize);
}


// ============================================================================


bitvf::DivergenceGrid bitvf::boxGrid (const otswap::Result& reconstruction,
                                      const std::optional<double> divergenceCellSize)
{
  const otswap::Result& rec = reconstruction;

  if (rec.geometry != otswap::Geometry::Box)
    throw Error("boxGrid: the reconstruction is not a box result; use the lightcone grid builder");

  if (rec.nRealizations == 0)
    throw Error("boxGrid: the reconstruction is malformed: it holds no realizations");

  internal::check_reconstruction_box(rec, "boxGrid");

  if (!std::isfinite(rec.mps) || rec.mps <= 0.)
    throw Error("boxGrid: the reconstruction's mps is " + std::to_string(rec.mps) +
                " in coordinate units; it must be finite and positive");

  if (divergenceCellSize && (!std::isfinite(*divergenceCellSize) || *divergenceCellSize <= 0.))
    throw Error("the divergence cell size is " + std::to_string(*divergenceCellSize) +
                " mean particle separations; it must be finite and positive");

  const double cellSizeMps = divergenceCellSize
    ? *divergenceCellSize
    : std::min(1., otswap::maths::cbrt(20. / rec.nRealizations));

  internal::BoxSearch search;
  search.box = *rec.box;
  search.mps = rec.mps;
  return internal::box_grid(search, cellSizeMps);
}
