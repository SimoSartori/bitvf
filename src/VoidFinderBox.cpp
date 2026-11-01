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
#include <sstream>
#include <string>
#include <utility>

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  /**
   *  @brief Write the selection lines of a box search, at Normal: the
   *  tracers outside the box limits, when the search has limits, then those
   *  outside the grid, when grid names one.
   *  @param grid the name of the grid in the line, or null for no line
   */
  void report_selection (const bitvf::internal::Report& report, const bitvf::internal::BoxSearch& search,
                         const char* grid)
  {
    if (!report.normal()) return;
    const bitvf::Counts& counts = search.counts;
    const std::size_t kept = counts.tracers - counts.tracersOutsideBox;

    if (search.limits) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "kept " << kept << " of " << counts.tracers << (counts.tracers == 1 ? " tracer: " : " tracers: ")
          << counts.tracersOutsideBox << " outside the box ";
      bitvf::internal::Report::box(out, search.limits->lower, search.limits->upper);
      report.write(out.str());
    }

    if (grid != nullptr) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "kept " << kept - counts.tracersOutsideGrid << " of " << kept << (kept == 1 ? " tracer: " : " tracers: ")
          << counts.tracersOutsideGrid << " outside " << grid;
      report.write(out.str());
    }
  }

  /**
   *  @brief Write the Detailed lines of the checks of findVoidsBox, the
   *  first step: the tracers and displacements checked, then the mps.
   *  @param mpsSource the source of the mps, in the parentheses of its line
   */
  void report_checked (bitvf::internal::Report& report, const bitvf::internal::BoxSearch& search,
                       const std::string& mpsSource)
  {
    if (!report.detailed()) return;
    const bitvf::Counts& counts = search.counts;

    std::ostringstream out = bitvf::internal::Report::line();
    out << "checked " << counts.tracers << (counts.tracers == 1 ? " tracer and " : " tracers and ")
        << counts.displacements << (counts.displacements == 1 ? " valid displacement" : " valid displacements");
    report.step(out);

    std::ostringstream mps = bitvf::internal::Report::line();
    mps << "mean particle separation " << search.mps << " in coordinate units (" << mpsSource << ")";
    report.write(mps.str());
  }

  /**
   *  @brief Write the call line of a box search, at Normal.
   *  @param name the public function
   *  @param without whether the line ends with the tracers read without a
   *  valid displacement
   */
  void report_call (const bitvf::internal::Report& report, const char* name, const bitvf::Result& result,
                    const bool without)
  {
    if (!report.normal()) return;
    const bitvf::Counts& counts = result.counts;

    std::ostringstream out = bitvf::internal::Report::line();
    out << name << ": ";
    bitvf::internal::Report::basins(out, result.voids.nBasins, result.overdense.nBasins);
    out << " from " << counts.tracers << (counts.tracers == 1 ? " tracer " : " tracers ");
    bitvf::internal::Report::time(out, result.elapsedSeconds);
    if (without)
      out << "; " << counts.tracers - counts.tracersOutsideBox - counts.tracersOutsideGrid - counts.tracersUsed
          << " without a valid displacement";
    report.write(out.str());
  }

  /**
   *  @brief The steps of a box search from the smoothing on: smoothing,
   *  watershed, catalogues, then the fields of the result.
   *  @param cells the grid with the divergence before smoothing, weight
   *  and the records of the search set; moved into Result::divergence
   *  when kept
   *  @param result the result, with config set
   *  @param report the report of the call, whose start is that of the
   *  public call, for Result::elapsedSeconds
   */
  void finish_box (const bitvf::internal::BoxSearch& search, const bitvf::internal::BoxDivergence& div,
                   bitvf::DivergenceGrid& cells, bitvf::Result& result, bitvf::internal::Report& report)
  {
    const bitvf::Config& config = result.config;

    bitvf::internal::box_smooth(cells, search.sigmaLength,
                                config.keepUnsmoothedDivergence ? &cells.divergenceUnsmoothed : nullptr);
    if (report.detailed()) {
      std::ostringstream out = bitvf::internal::Report::line();
      if (search.sigmaLength > 0.) {
        out << "smoothed with sigma " << search.sigmaLength << " in coordinate units ("
            << search.sigmaLength / search.mps << " mps)";
        report.step(out);
      }
      else {
        report.lap();
        report.write("not smoothed (sigma 0)");
      }
    }

    bitvf::internal::BoxBasins basins;
    bitvf::internal::box_watershed(cells, basins);
    if (report.detailed()) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "watershed into ";
      bitvf::internal::Report::basins(out, basins.voids.seed.size(), basins.overdense.seed.size());
      report.step(out);
    }
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
    if (report.detailed()) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "catalogues of the basins";
      report.step(out);
    }

    if (config.keepAssignment)
      result.assignment = bitvf::internal::box_assignment(search, div, cells, result.voids, result.overdense);

    if (config.keepDivergence) {
      bitvf::internal::box_centres(cells);
      cells.volume.assign(cells.nCells, cells.cellSize * cells.cellSize * cells.cellSize);
      result.divergence = std::move(cells);
    }
    if (report.detailed() && (config.keepAssignment || config.keepDivergence)) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "kept the " << (!config.keepDivergence ? "tracer assignment"
                             : !config.keepAssignment ? "divergence grid"
                             : "tracer assignment and the divergence grid");
      report.step(out);
    }

    result.elapsedSeconds = report.elapsed();
  }

  /**
   *  @brief The search shared by both overloads of findVoidsBox, once the
   *  arguments are checked: the grid, the divergence, then finish_box.
   *  @param grid the grid to reuse, or null to build one
   *  @param boxSource what the grid is built on, without grid, in the
   *  parentheses of the grid line
   *  @param report the report of the call, with the checks reported
   */
  bitvf::Result search_box (const bitvf::internal::BoxSearch& search, const bitvf::Config& config,
                            const bitvf::DivergenceGrid* grid, const char* boxSource,
                            bitvf::internal::Report& report)
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
    if (report.detailed()) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "grid of ";
      bitvf::internal::Report::shape(out, cells.nCellsAxis);
      out << " of side " << cells.cellSize << " in coordinate units (" << cells.cellSizeMps << " mps, "
          << (grid != nullptr ? "reused" : config.divergenceCellSize ? "given" : "default") << ") on ";
      bitvf::internal::Report::box(out, cells.lower, cells.upper);
      out << " (" << (grid != nullptr ? "the reused grid" : boxSource) << ")";
      report.step(out);
    }

    const bitvf::internal::BoxDivergence div = bitvf::internal::box_divergence(search, cells);
    if (report.detailed()) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "divergence of the valid displacements of " << div.nUsed << (div.nUsed == 1 ? " tracer" : " tracers");
      report.step(out);
    }

    result.nRealizations = search.nRealizations;
    finish_box(search, div, cells, result, report);
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
  internal::Report report(config.verbosity, config.log);
  const internal::BoxSearch search =
    internal::check_box(tracers, displacement, valid, mps, config, limits, ids, grid);
  report_selection(report, search, grid != nullptr ? "the reused grid" : nullptr);
  report_checked(report, search, "given");

  Result result = search_box(search, config, grid,
                             limits ? "the box limits" : "the tracers and their displacement endpoints", report);
  report_call(report, "findVoidsBox", result, true);
  return result;
}


// ============================================================================


bitvf::Result bitvf::findVoidsBox (const otswap::Result& reconstruction,
                                   const Config& config,
                                   const DivergenceGrid* grid)
{
  internal::Report report(config.verbosity, config.log);
  const internal::BoxSearch search = internal::check_box(reconstruction, config, grid);
  report_selection(report, search, grid != nullptr ? "the reused grid" : nullptr);
  if (grid == nullptr) report_checked(report, search, "from the reconstruction");
  else {
    std::ostringstream source = internal::Report::line();
    source << "from the reused grid; the reconstruction's is " << reconstruction.mps;
    report_checked(report, search, source.str());
  }

  Result result = search_box(search, config, grid,
                             reconstruction.selection.boxLimits ? "the box limits" : "the reconstruction's box", report);
  report_call(report, "findVoidsBox", result, true);
  return result;
}


// ============================================================================


bitvf::Result bitvf::findVoidsBoxFromDivergence (const DivergenceGrid& divergence,
                                                 const std::vector<double>& tracers,
                                                 const Config& config,
                                                 const std::vector<std::int64_t>& ids)
{
  internal::Report report(config.verbosity, config.log);
  const internal::BoxSearch search = internal::check_box_from_divergence(divergence, tracers, config, ids);
  report_selection(report, search, "the divergence grid");
  if (report.detailed()) {
    std::ostringstream out = internal::Report::line();
    out << "checked " << search.counts.tracers << (search.counts.tracers == 1 ? " tracer" : " tracers");
    report.step(out);

    std::ostringstream mps = internal::Report::line();
    mps << "mean particle separation " << divergence.mps << " in coordinate units (from the divergence)";
    report.write(mps.str());

    std::ostringstream grid = internal::Report::line();
    grid << "grid of ";
    internal::Report::shape(grid, divergence.nCellsAxis);
    grid << " of side " << divergence.cellSize << " in coordinate units (" << divergence.cellSizeMps
         << " mps, from the divergence) on ";
    internal::Report::box(grid, divergence.lower, divergence.upper);
    if (divergence.nRealizations > 0) {
      grid << ", from " << divergence.nRealizations
           << (divergence.nRealizations == 1 ? " realization" : " realizations");
      if (divergence.nRealizationsMean > 1) grid << ", each the mean of " << divergence.nRealizationsMean;
    }
    report.write(grid.str());
  }

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
  if (report.detailed()) {
    const std::size_t placed =
      search.counts.tracers - search.counts.tracersOutsideBox - search.counts.tracersOutsideGrid;
    std::ostringstream out = internal::Report::line();
    out << "placed " << placed << (placed == 1 ? " tracer" : " tracers") << " in the cells";
    report.step(out);
  }

  finish_box(search, div, cells, result, report);
  report_call(report, "findVoidsBoxFromDivergence", result, false);
  return result;
}
