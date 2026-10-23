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
 *  @file src/Checks.cpp
 *
 *  @brief The checks of the arguments of the void finder, and the
 *  selection of what a search reads.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

#include "internal.h"

namespace {

  const char* const kAxes[3] = {"x", "y", "z"};

  /// A double written with 17 significant digits, so that two values that
  /// differ are written differently.
  std::string exact (const double x)
  {
    std::ostringstream out;
    out << std::setprecision(17) << x;
    return out.str();
  }

  /**
   *  @brief Check a flat array of Cartesian coordinates, every entry.
   *  @param name the name of the array in the messages
   *  @param allowEmpty whether an empty array is accepted
   *  @return the number of rows
   */
  std::size_t check_coordinates (const std::vector<double>& a, const std::string& name,
                                 const bool allowEmpty)
  {
    if (a.empty()) {
      if (allowEmpty) return 0;
      throw bitvf::Error(name + " is empty");
    }

    if (a.size() % 3 != 0)
      throw bitvf::Error(name + " holds " + std::to_string(a.size()) +
                         " entries, which is not a multiple of three");

    for (std::size_t i = 0; i < a.size(); ++i)
      if (!std::isfinite(a[i]))
        throw bitvf::Error(name + " holds a non-finite value at entry " + std::to_string(i) +
                           " (object " + std::to_string(i/3) + ", component " +
                           std::to_string(i%3) + ")");

    return a.size() / 3;
  }

  /**
   *  @brief Check that every entry of a validity array is 0 or 1.
   *  @param name the name of the array in the messages
   */
  void check_valid (const std::vector<std::uint8_t>& valid, const std::string& name)
  {
    for (std::size_t i = 0; i < valid.size(); ++i)
      if (valid[i] > 1)
        throw bitvf::Error(name + " holds " + std::to_string(valid[i]) + " at entry " +
                           std::to_string(i) + "; each entry must be 0 or 1");
  }

  /**
   *  @brief Check box limits.
   *  @param name the name of the limits in the messages
   */
  void check_limits (const otswap::BoxLimits& limits, const std::string& name)
  {
    for (std::size_t k = 0; k < 3; ++k)
      for (const double bound : {limits.lower[k], limits.upper[k]})
        if (!std::isfinite(bound))
          throw bitvf::Error(name + " hold a bound of " + std::to_string(bound) + " on axis " +
                             kAxes[k] + "; every bound must be finite");

    for (std::size_t k = 0; k < 3; ++k)
      if (!(limits.lower[k] < limits.upper[k]))
        throw bitvf::Error(name + " are [" + std::to_string(limits.lower[k]) + ", " +
                           std::to_string(limits.upper[k]) + "] on axis " + kAxes[k] +
                           "; the lower bound must be below the upper one");
  }

  /**
   *  @brief Check the IDs of the tracers: one per tracer, or none, and no
   *  value twice.
   *  @param name the name of the IDs in the messages
   */
  void check_ids (const std::vector<std::int64_t>& ids, const std::size_t nTracers,
                  const std::string& name)
  {
    if (ids.empty()) return;
    if (ids.size() != nTracers)
      throw bitvf::Error(name + " holds " + std::to_string(ids.size()) + " entries for " +
                         std::to_string(nTracers) + " tracers; give one per tracer, or none");

    std::vector<std::pair<std::int64_t, std::size_t>> sorted(ids.size());
    for (std::size_t i = 0; i < ids.size(); ++i) sorted[i] = {ids[i], i};
    std::sort(sorted.begin(), sorted.end());
    std::size_t start = 0, first = 0, second = ids.size();
    for (std::size_t k = 1; k < sorted.size(); ++k) {
      if (sorted[k].first != sorted[k-1].first) start = k;
      else if (k == start + 1 && sorted[k].second < second) {
        first = sorted[start].second;
        second = sorted[k].second;
      }
    }
    if (second < ids.size())
      throw bitvf::Error(name + " holds " + std::to_string(ids[first]) + " twice, at entries " +
                         std::to_string(first) + " and " + std::to_string(second) +
                         "; each tracer needs its own ID");
  }

  /**
   *  @brief The checks and the selection shared by both overloads, once
   *  the arrays, mps, the limits and the IDs are checked and the tracers
   *  outside the limits are known.
   *
   *  Removes the randoms outside the limits; with a grid, flags the
   *  tracers and removes the randoms outside it; checks the valid
   *  displacements of the tracers read; counts.
   *
   *  @param search the search, with every field but outsideGrid, randoms
   *  and counts other than tracers and tracersOutsideBox set
   *  @param randoms the randoms given
   *  @param displacementName the name of the displacements in the messages
   */
  void select_box (bitvf::internal::BoxSearch& search, const std::vector<double>& randoms,
                   const bitvf::Config& config, const bitvf::DivergenceGrid* grid,
                   const bool mpsMustMatch, const std::string& displacementName)
  {
    bitvf::internal::check_config(config, grid != nullptr);
    if (grid != nullptr) bitvf::internal::check_box_grid(*grid, search.mps, mpsMustMatch);

    const std::size_t N = search.nTracers;
    const std::size_t nRandoms = randoms.size() / 3;
    search.randomsGiven = nRandoms > 0;

    std::vector<std::uint8_t> keep(nRandoms, 1);
    if (search.limits) {
      const otswap::BoxLimits& box = *search.limits;
      for (std::size_t i = 0; i < nRandoms; ++i)
        for (std::size_t k = 0; k < 3; ++k)
          if (!(randoms[3*i+k] >= box.lower[k] && randoms[3*i+k] <= box.upper[k])) {
            keep[i] = 0;
            ++search.counts.randomsOutsideBox;
            break;
          }
    }

    if (grid != nullptr) {
      search.outsideGrid.assign(N, 0);
      for (std::size_t i = 0; i < N; ++i)
        if (search.inBox[i] && bitvf::internal::box_cell(*grid, &(*search.tracers)[3*i]) < 0) {
          search.outsideGrid[i] = 1;
          ++search.counts.tracersOutsideGrid;
        }
      for (std::size_t i = 0; i < nRandoms; ++i)
        if (keep[i] && bitvf::internal::box_cell(*grid, &randoms[3*i]) < 0) {
          keep[i] = 0;
          ++search.counts.randomsOutsideGrid;
        }
    }

    const std::size_t nKept = nRandoms - search.counts.randomsOutsideBox -
                              search.counts.randomsOutsideGrid;
    search.randoms.reserve(3 * nKept);
    for (std::size_t i = 0; i < nRandoms; ++i)
      if (keep[i]) search.randoms.insert(search.randoms.end(), &randoms[3*i], &randoms[3*i+3]);

    const std::vector<double>& displacement = *search.displacement;
    const std::vector<std::uint8_t>& valid = *search.valid;
    std::size_t nRead = 0;
    for (std::size_t r = 0; r < search.nRealizations; ++r)
      for (std::size_t i = 0; i < N; ++i) {
        const std::size_t row = r * N + i;
        if (!valid.empty() && !valid[row]) continue;
        ++search.counts.displacements;
        if (!search.inBox[i] || (!search.outsideGrid.empty() && search.outsideGrid[i])) continue;
        for (std::size_t k = 0; k < 3; ++k) {
          if (!std::isfinite(displacement[3*row+k]))
            throw bitvf::Error(displacementName + " holds a non-finite value at entry " +
                               std::to_string(3*row+k) + " (realization " + std::to_string(r) +
                               ", tracer " + std::to_string(i) + ", component " +
                               std::to_string(k) + "), a valid displacement of a tracer the "
                               "search reads");
          search.maxDisplacement = std::max(search.maxDisplacement, std::abs(displacement[3*row+k]));
        }
        ++nRead;
      }

    if (nRead == 0)
      throw bitvf::Error("no valid displacement is left to read: of " +
                         std::to_string(search.counts.displacements) +
                         " valid displacements, none belongs to a tracer inside the box limits"
                         " and the grid (" + std::to_string(search.counts.tracersOutsideBox) +
                         " tracers outside the limits, " +
                         std::to_string(search.counts.tracersOutsideGrid) +
                         " outside the grid)");

    if (search.randomsGiven && nKept == 0)
      throw bitvf::Error("no random is left: of " + std::to_string(nRandoms) + " randoms, " +
                         std::to_string(search.counts.randomsOutsideBox) +
                         " lie outside the box limits and " +
                         std::to_string(search.counts.randomsOutsideGrid) +
                         " outside the grid");

    search.counts.randoms = search.randomsGiven ? nKept : nRead;
  }

}


// ============================================================================


bitvf::internal::BoxSearch bitvf::internal::check_box (const std::vector<double>& tracers,
                                                       const std::vector<double>& displacement,
                                                       const std::vector<std::uint8_t>& valid,
                                                       const std::vector<double>& randoms,
                                                       const double mps,
                                                       const Config& config,
                                                       const std::optional<otswap::BoxLimits>& limits,
                                                       const std::vector<std::int64_t>& ids,
                                                       const DivergenceGrid* grid)
{
  const std::size_t N = check_coordinates(tracers, "tracers", false);

  if (displacement.empty() || displacement.size() % (3 * N) != 0)
    throw Error("displacement holds " + std::to_string(displacement.size()) + " entries for " +
                std::to_string(N) + " tracers; it must hold 3 * nRealizations * " +
                std::to_string(N) + " entries, with a whole nRealizations of at least 1");

  const std::size_t R = displacement.size() / (3 * N);
  if (R > UINT_MAX)
    throw Error("displacement holds " + std::to_string(R) + " realizations; at most " +
                std::to_string(UINT_MAX) + " are supported");

  if (!valid.empty() && valid.size() != R * N)
    throw Error("valid holds " + std::to_string(valid.size()) + " entries for " +
                std::to_string(R) + " realizations of " + std::to_string(N) +
                " tracers; give one per displacement, or none");
  check_valid(valid, "valid");

  check_coordinates(randoms, "randoms", true);

  if (!std::isfinite(mps) || mps <= 0.)
    throw Error("mps is " + std::to_string(mps) + " Mpc/h; it must be finite and positive");

  if (limits) check_limits(*limits, "the box limits");

  check_ids(ids, N, "ids");

  BoxSearch search;
  search.nTracers = N;
  search.nRealizations = (unsigned)R;
  search.tracers = &tracers;
  search.displacement = &displacement;
  search.valid = &valid;
  search.ids = &ids;
  search.mps = mps;
  search.limits = limits;
  search.counts.tracers = N;

  search.inBox.assign(N, 1);
  if (limits)
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t k = 0; k < 3; ++k)
        if (!(tracers[3*i+k] >= limits->lower[k] && tracers[3*i+k] <= limits->upper[k])) {
          search.inBox[i] = 0;
          ++search.counts.tracersOutsideBox;
          break;
        }

  select_box(search, randoms, config, grid, true, "displacement");
  return search;
}


// ============================================================================


bitvf::internal::BoxSearch bitvf::internal::check_box (const otswap::Result& reconstruction,
                                                       const Config& config,
                                                       const std::vector<double>& randoms,
                                                       const DivergenceGrid* grid)
{
  const otswap::Result& rec = reconstruction;

  if (rec.geometry != otswap::Geometry::Box)
    throw Error("the reconstruction is not a box result; use the lightcone void finder");

  const std::size_t N = rec.nObjects;
  const std::size_t R = rec.nRealizations;
  const std::string malformed = "the reconstruction is malformed: ";
  if (N == 0) throw Error(malformed + "it holds no objects");
  if (R == 0) throw Error(malformed + "it holds no realizations");

  auto check_size = [&] (const std::size_t size, const std::size_t expected,
                         const std::string& name, const bool allowEmpty) {
    if (size != expected && !(allowEmpty && size == 0))
      throw Error(malformed + name + " holds " + std::to_string(size) + " entries for " +
                  std::to_string(R) + " realizations of " + std::to_string(N) + " objects");
  };
  check_size(rec.tracers.size(), 3 * N, "tracers", false);
  check_size(rec.displacement.size(), 3 * R * N, "displacement", false);
  check_size(rec.valid.size(), R * N, "valid", false);
  check_size(rec.outsideBox.size(), N, "outsideBox", true);
  check_size(rec.ids.size(), N, "ids", true);

  check_valid(rec.valid, "the reconstruction's valid");
  check_ids(rec.ids, N, "the reconstruction's ids");

  BoxSearch search;
  search.nTracers = N;
  search.nRealizations = (unsigned)R;
  search.tracers = &rec.tracers;
  search.displacement = &rec.displacement;
  search.valid = &rec.valid;
  search.ids = &rec.ids;
  search.limits = rec.selection.boxLimits;
  search.counts.tracers = N;

  search.inBox.assign(N, 1);
  for (std::size_t i = 0; i < rec.outsideBox.size(); ++i)
    if (rec.outsideBox[i]) {
      search.inBox[i] = 0;
      ++search.counts.tracersOutsideBox;
    }

  for (std::size_t i = 0; i < N; ++i)
    if (search.inBox[i])
      for (std::size_t k = 0; k < 3; ++k)
        if (!std::isfinite(rec.tracers[3*i+k]))
          throw Error("the reconstruction's tracers hold a non-finite value at entry " +
                      std::to_string(3*i+k) + " (object " + std::to_string(i) + ", component " +
                      std::to_string(k) + ")");

  check_coordinates(randoms, "randoms", true);

  if (grid == nullptr) {
    if (!std::isfinite(rec.mps) || rec.mps <= 0.)
      throw Error("the reconstruction's mps is " + std::to_string(rec.mps) +
                  " Mpc/h; it must be finite and positive");
    search.mps = rec.mps;
  }
  else search.mps = grid->mps;

  if (search.limits) check_limits(*search.limits, "the reconstruction's box limits");

  select_box(search, randoms, config, grid, false, "the reconstruction's displacement");
  return search;
}


// ============================================================================


void bitvf::internal::check_config (const Config& config, const bool gridGiven)
{
  const double sigma = config.smoothing.sigma;
  if (!std::isfinite(sigma) || sigma < 0.)
    throw Error("the smoothing sigma is " + std::to_string(sigma) +
                "; it must be finite and non-negative");

  if (config.divergenceCellSize) {
    const double cell = *config.divergenceCellSize;
    if (!std::isfinite(cell) || cell <= 0.)
      throw Error("the divergence cell size is " + std::to_string(cell) +
                  " mean particle separations; it must be finite and positive");
  }

  if (config.keepUnsmoothedDivergence && !config.keepDivergence)
    throw Error("keeping the unsmoothed divergence requires keeping the divergence");

  if (gridGiven && config.divergenceCellSize)
    throw Error("the divergence cell size is set together with a grid to reuse; leave it empty, "
                "the cells are the grid's");
}


// ============================================================================


void bitvf::internal::check_box_grid (const DivergenceGrid& grid, const double mps,
                                      const bool mpsMustMatch)
{
  const std::string bad = "the grid to reuse ";

  if (grid.geometry != otswap::Geometry::Box)
    throw Error(bad + "is not a box grid");

  if (!std::isfinite(grid.cellSize) || grid.cellSize <= 0.)
    throw Error(bad + "has a cell size of " + std::to_string(grid.cellSize) +
                " Mpc/h; it must be finite and positive");

  for (std::size_t k = 0; k < 3; ++k) {
    if (!std::isfinite(grid.lower[k]) || !std::isfinite(grid.upper[k]))
      throw Error(bad + "has a bound of " +
                  std::to_string(std::isfinite(grid.lower[k]) ? grid.upper[k] : grid.lower[k]) +
                  " on axis " + kAxes[k] + "; every bound must be finite");
    if (!(grid.lower[k] < grid.upper[k]))
      throw Error(bad + "spans [" + std::to_string(grid.lower[k]) + ", " +
                  std::to_string(grid.upper[k]) + "] on axis " + kAxes[k] +
                  "; the lower bound must be below the upper one");
  }

  const std::size_t n = grid.nCellsAxis[0];
  if (n == 0 || grid.nCellsAxis[1] != n || grid.nCellsAxis[2] != n)
    throw Error(bad + "has " + std::to_string(grid.nCellsAxis[0]) + ", " +
                std::to_string(grid.nCellsAxis[1]) + " and " + std::to_string(grid.nCellsAxis[2]) +
                " cells along x, y and z; the three counts must be equal and positive");

  if (n > kMaxCellsAxis || grid.nCells != n * n * n)
    throw Error(bad + "has " + std::to_string(grid.nCells) + " cells for " + std::to_string(n) +
                " along each axis; it must have their cube, and at most " +
                std::to_string(kMaxCellsAxis) + " along each axis");

  if (!std::isfinite(grid.mps) || grid.mps <= 0.)
    throw Error(bad + "has an mps of " + std::to_string(grid.mps) +
                " Mpc/h; it must be finite and positive");

  const bool lightconeSet =
    !std::isnan(grid.distanceMin) || !std::isnan(grid.distanceMax) ||
    !std::isnan(grid.raMin) || !std::isnan(grid.raMax) ||
    !std::isnan(grid.decMin) || !std::isnan(grid.decMax) ||
    !std::isnan(grid.cellSizeMps) || !grid.mpsProfile.redshift.empty() || !grid.mpsProfile.mps.empty() ||
    !grid.mpsProfile.count.empty();
  if (lightconeSet)
    throw Error(bad + "is a box grid with a field of a lightcone grid set");

  if (mpsMustMatch && mps != grid.mps)
    throw Error("mps is " + exact(mps) + " Mpc/h, but the grid to reuse was built with " +
                exact(grid.mps) + "; they must be equal");
}
