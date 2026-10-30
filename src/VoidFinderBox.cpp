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
   *  @brief The steps of a box search from the smoothing on: smoothing,
   *  watershed, catalogues, then the fields of the result.
   *  @param cells the grid with the divergence before smoothing, weight
   *  and the records of the search set; moved into Result::divergence
   *  when kept
   *  @param result the result, with config set
   *  @param start the time the public call began, for
   *  Result::elapsedSeconds
   */
  void finish_box (const bitvf::internal::BoxSearch& search, const bitvf::internal::BoxDivergence& div,
                   bitvf::DivergenceGrid& cells, bitvf::Result& result,
                   const std::chrono::steady_clock::time_point start)
  {
    const bitvf::Config& config = result.config;

    bitvf::internal::box_smooth(cells, search.sigmaLength,
                                config.keepUnsmoothedDivergence ? &cells.divergenceUnsmoothed : nullptr);
    bitvf::internal::BoxBasins basins;
    bitvf::internal::box_watershed(cells, basins);
    if (basins.voids.seed.empty())
      throw bitvf::Error("no void found: the divergence is 0 or positive in every cell");
    bitvf::internal::box_basins(&search, &div, cells, basins, result.voids, result.overdense);

    result.nTracers = search.nTracers;
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
  }

  /**
   *  @brief The search shared by both overloads of findVoidsBox, once the
   *  arguments are checked: the grid, the divergence, then finish_box.
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
    result.nRealizations = search.nRealizations;
    finish_box(search, div, cells, result, start);
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


// ============================================================================


bitvf::Result bitvf::findVoidsBoxFromDivergence (const DivergenceGrid& divergence,
                                                 const std::vector<double>& tracers,
                                                 const Config& config,
                                                 const std::vector<std::int64_t>& ids)
{
  const auto start = std::chrono::steady_clock::now();
  const internal::BoxSearch search = internal::check_box_from_divergence(divergence, tracers, config, ids);

  Result result;
  result.config = config;
  result.config.divergenceCellSize = divergence.cellSizeMps;
  result.config.nRealizationsMean = divergence.nRealizationsMean;
  result.nRealizations = divergence.nRealizations;

  DivergenceGrid cells;
  cells.geometry = divergence.geometry;
  cells.mps = divergence.mps;
  cells.cellSize = divergence.cellSize;
  cells.cellSizeMps = divergence.cellSizeMps;
  cells.lower = divergence.lower;
  cells.upper = divergence.upper;
  cells.nCellsAxis = divergence.nCellsAxis;
  cells.nCells = divergence.nCells;
  cells.limits = divergence.limits;
  cells.nRealizations = divergence.nRealizations;
  cells.nRealizationsMean = divergence.nRealizationsMean;
  cells.divergence = divergence.divergenceUnsmoothed;
  cells.weight = divergence.weight;
  cells.nDisplacements = divergence.nDisplacements;

  const internal::BoxDivergence div = internal::box_tracer_cells(search, cells);
  if (!divergence.nTracers.empty())
    for (std::size_t c = 0; c < cells.nCells; ++c)
      if (cells.nTracers[c] != divergence.nTracers[c])
        throw Error("the tracers do not match the divergence: cell " + std::to_string(c) + " holds " +
                    std::to_string(cells.nTracers[c]) + " of the tracers read, but the grid's nTracers is " +
                    std::to_string(divergence.nTracers[c]) + "; pass the tracers of the search that produced the field");

  finish_box(search, div, cells, result, start);
  return result;
}
