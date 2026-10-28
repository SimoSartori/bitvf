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
 *  @file src/VoidFinderBox.cpp
 *
 *  @brief The box void finder, findVoidsBox.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <chrono>
#include <utility>

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  /**
   *  @brief The search shared by both overloads, once the arguments are
   *  checked: the grid, then the steps of the search, then the result.
   *  @param grid the grid to reuse, or null to build one
   *  @param start the time the public call began, for
   *  Result::elapsedSeconds
   */
  bitvf::Result search_box (const bitvf::internal::BoxSearch& search, const bitvf::Config& config,
                            const bitvf::DivergenceGrid* grid,
                            const std::chrono::steady_clock::time_point start)
  {
    bitvf::Result result;
    result.config = config;

    bitvf::DivergenceGrid cells;
    if (grid == nullptr) {
      const double cellSizeMps = config.divergenceCellSize
        ? *config.divergenceCellSize
        : std::min(1., otswap::maths::cbrt(20. / search.nRealizations));
      result.config.divergenceCellSize = cellSizeMps;
      cells = bitvf::internal::box_grid(search, cellSizeMps);
    }
    else {
      result.config.divergenceCellSize = grid->cellSizeMps;
      cells.geometry = grid->geometry;
      cells.mps = grid->mps;
      cells.cellSize = grid->cellSize;
      cells.cellSizeMps = grid->cellSizeMps;
      cells.lower = grid->lower;
      cells.upper = grid->upper;
      cells.nCellsAxis = grid->nCellsAxis;
      cells.nCells = grid->nCells;
    }
    cells.limits = search.limits;
    cells.nRealizations = search.nRealizations;
    cells.nRealizationsMean = config.nRealizationsMean;

    const bitvf::internal::BoxDivergence div = bitvf::internal::box_divergence(search, cells);
    bitvf::internal::box_smooth(cells, search.sigmaLength,
                                config.keepUnsmoothedDivergence ? &cells.divergenceUnsmoothed : nullptr);
    bitvf::internal::BoxBasins basins;
    bitvf::internal::box_watershed(cells, basins);
    bitvf::internal::box_basins(search, div, cells, basins, result.voids, result.overdense);

    result.nTracers = search.nTracers;
    result.nRealizations = search.nRealizations;
    result.geometry = otswap::Geometry::Box;
    result.ids = *search.ids;
    result.outsideGrid = search.outsideGrid;
    result.limits = search.limits;
    result.mps = search.mps;
    result.counts = search.counts;
    result.counts.tracersUsed = div.nUsed;

    if (config.keepAssignment)
      result.assignment = bitvf::internal::box_assignment(search, div, cells, result.voids, result.overdense);

    if (config.keepDivergence) {
      bitvf::internal::box_centres(cells);
      cells.volume.assign(cells.nCells, cells.cellSize * cells.cellSize * cells.cellSize);
      result.divergence = std::move(cells);
    }

    result.elapsedSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    return result;
  }

}


// ============================================================================


bitvf::Result bitvf::findVoidsBox (const std::vector<double>& tracers,
                                   const std::vector<double>& displacement,
                                   const std::vector<std::uint8_t>& valid,
                                   const double mps,
                                   const Config& config,
                                   const std::optional<otswap::BoxLimits>& limits,
                                   const std::vector<std::int64_t>& ids,
                                   const DivergenceGrid* grid)
{
  const auto start = std::chrono::steady_clock::now();
  const internal::BoxSearch search =
    internal::check_box(tracers, displacement, valid, mps, config, limits, ids, grid);
  return search_box(search, config, grid, start);
}


// ============================================================================


bitvf::Result bitvf::findVoidsBox (const otswap::Result& reconstruction,
                                   const Config& config,
                                   const DivergenceGrid* grid)
{
  const auto start = std::chrono::steady_clock::now();
  const internal::BoxSearch search = internal::check_box(reconstruction, config, grid);
  return search_box(search, config, grid, start);
}
