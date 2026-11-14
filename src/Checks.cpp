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
#include <limits>
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
   *  @brief Check that every entry of a validity array is 0 or 1.
   *  @param fn the public function that starts the message of a reconstruction's flags
   *  @param fromResult whether the flags are a reconstruction's
   */
  void check_valid (const std::vector<std::uint8_t>& valid, const std::string& fn, const bool fromResult)
  {
    for (std::size_t i = 0; i < valid.size(); ++i)
      if (valid[i] > 1)
        throw bitvf::Error((fromResult ? fn + ": the reconstruction's validity flags hold "
                                       : std::string("the validity flags hold ")) + std::to_string(valid[i]) + " at entry " +
                           std::to_string(i) + "; each must be 0 or 1");
  }

  /**
   *  @brief The smoothing sigma in coordinate units: sigma, times mps in
   *  Mps.
   *  @exception Error if sigma is positive and (2 sigmaLength) sigmaLength
   *  is 0 or not finite, too small or too large to build the kernel
   */
  double sigma_length (const bitvf::Smoothing& smoothing, const double mps)
  {
    const bool inMps = smoothing.unit == bitvf::LengthUnit::Mps;
    const double sigmaLength = inMps ? smoothing.sigma * mps : smoothing.sigma;
    const double twoSigmaSq = 2 * sigmaLength * sigmaLength;
    if (smoothing.sigma > 0. && (twoSigmaSq == 0. || !std::isfinite(twoSigmaSq)))
      throw bitvf::Error("the smoothing sigma is " + exact(smoothing.sigma) +
                         (inMps ? " mps (" : " in coordinate units (") + exact(sigmaLength) +
                         " in coordinate units), too " + (twoSigmaSq == 0. ? "small" : "large") +
                         " to build the kernel; use 0 to disable the smoothing");
    return sigmaLength;
  }

  /**
   *  @brief With a grid, flag and count the tracers inside the limits but
   *  outside every cell; then check the number of tracers read.
   *  @param search the search, with inBox and counts.tracersOutsideBox set
   *  @param grid the grid the tracers are read in, or null for none
   *  @exception Error if the search reads more tracers than the count of
   *  delta2Mps supports
   */
  void select_tracers (bitvf::internal::BoxSearch& search, const bitvf::DivergenceGrid* grid)
  {
    const std::size_t N = search.nTracers;
    if (grid != nullptr) {
      search.outsideGrid.assign(N, 0);
      for (std::size_t i = 0; i < N; ++i)
        if (search.inBox[i] && bitvf::internal::box_cell(*grid, &(*search.tracers)[3*i]) < 0) {
          search.outsideGrid[i] = 1;
          ++search.counts.tracersOutsideGrid;
        }
    }

    const std::size_t nTracersRead = N - search.counts.tracersOutsideBox - search.counts.tracersOutsideGrid;
    if (nTracersRead > std::numeric_limits<unsigned>::max())
      throw bitvf::Error("the search reads " + std::to_string(nTracersRead) + " tracers, more than " +
                         std::to_string(std::numeric_limits<unsigned>::max()) +
                         ", the most that the density contrast within 2 mean particle separations of a centre supports");
  }

  /**
   *  @brief The checks and the selection shared by both overloads, once
   *  the arrays, mps, the limits and the IDs are checked and the tracers
   *  outside the limits are known.
   *
   *  Sets and checks sigmaLength; with a grid, flags the tracers outside
   *  it; checks the valid displacements of the tracers read; counts.
   *
   *  @param search the search, with every field but outsideGrid,
   *  sigmaLength and counts other than tracers and tracersOutsideBox set
   *  @param fn the public function that starts the messages about a grid and a reconstruction
   *  @param fromResult whether the displacements are a reconstruction's
   */
  void select_box (bitvf::internal::BoxSearch& search, const bitvf::Config& config, const bitvf::DivergenceGrid* grid,
                   const bool mpsMustMatch, const std::string& fn, const bool fromResult)
  {
    bitvf::internal::check_config(config, grid != nullptr);
    if (grid != nullptr) bitvf::internal::check_box_grid(*grid, search.mps, mpsMustMatch, fn);

    search.sigmaLength = sigma_length(config.smoothing, search.mps);

    select_tracers(search, grid);

    const std::size_t N = search.nTracers;
    const std::vector<double>& displacement = *search.displacement;
    const std::vector<std::uint8_t>& valid = *search.valid;
    const std::size_t rows = (std::size_t)search.nRealizations * N;
    constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
    std::size_t nValid = 0, nRead = 0, firstNonFinite = none;
    double maxDisplacement = 0.;
#pragma omp parallel for schedule(static) reduction(+:nValid, nRead) \
  reduction(min:firstNonFinite) reduction(max:maxDisplacement)
    for (std::size_t row = 0; row < rows; ++row) {
      const std::size_t i = row % N;
      if (!valid.empty() && !valid[row]) continue;
      ++nValid;
      if (!search.inBox[i] || (!search.outsideGrid.empty() && search.outsideGrid[i])) continue;
      for (std::size_t k = 0; k < 3; ++k) {
        const double d = displacement[3*row+k];
        if (!std::isfinite(d)) firstNonFinite = std::min(firstNonFinite, 3*row+k);
        else maxDisplacement = std::max(maxDisplacement, std::abs(d));
      }
      ++nRead;
    }

    if (firstNonFinite != none) {
      const std::size_t row = firstNonFinite / 3;
      throw bitvf::Error((fromResult ? fn + ": the reconstruction's displacements hold " : std::string("the displacements hold ")) +
                         "a non-finite value at entry " +
                         std::to_string(firstNonFinite) + " (realization " +
                         std::to_string(row / N) + ", tracer " + std::to_string(row % N) +
                         ", component " + std::to_string(firstNonFinite % 3) +
                         "), a valid displacement of a tracer the search reads");
    }
    search.counts.displacements = nValid;
    search.maxDisplacement = maxDisplacement;

    if (nRead == 0)
      throw bitvf::Error("no valid displacement is left to read: of " +
                         std::to_string(search.counts.displacements) +
                         " valid displacements, none belongs to a tracer inside the box limits"
                         " and the grid (" + std::to_string(search.counts.tracersOutsideBox) +
                         " tracers outside the limits, " +
                         std::to_string(search.counts.tracersOutsideGrid) +
                         " outside the grid)");
  }

}


// ============================================================================


void bitvf::internal::check_limits (const otswap::BoxLimits& limits, const std::string& name, const std::string& fn)
{
  const std::string at = fn.empty() ? "" : fn + ": ";
  for (std::size_t k = 0; k < 3; ++k)
    for (const double bound : {limits.lower[k], limits.upper[k]})
      if (!std::isfinite(bound))
        throw Error(at + name + " hold a bound of " + std::to_string(bound) + " on axis " +
                    kAxes[k] + "; every bound must be finite");

  for (std::size_t k = 0; k < 3; ++k)
    if (!std::isfinite(limits.upper[k] - limits.lower[k]))
      throw Error(at + name + " are [" + std::to_string(limits.lower[k]) + ", " +
                  std::to_string(limits.upper[k]) + "] on axis " + kAxes[k] +
                  "; their width must be finite");

  for (std::size_t k = 0; k < 3; ++k)
    if (!(limits.lower[k] < limits.upper[k]))
      throw Error(at + name + " are [" + std::to_string(limits.lower[k]) + ", " +
                  std::to_string(limits.upper[k]) + "] on axis " + kAxes[k] +
                  "; the lower bound must be below the upper one");
}


// ============================================================================


std::size_t bitvf::internal::check_coordinates (const std::vector<double>& a, const std::string& name,
                                                const bool allowEmpty, const std::string& fn)
{
  if (a.empty()) {
    if (allowEmpty) return 0;
    throw Error(fn + ": " + name + " is empty");
  }

  if (a.size() % 3 != 0)
    throw Error(fn + ": " + name + " holds " + std::to_string(a.size()) +
                " entries, which is not a multiple of three");

  for (std::size_t i = 0; i < a.size(); ++i)
    if (!std::isfinite(a[i]))
      throw Error((name.compare(0, 4, "the ") == 0 ? name : "the " + name) + (name.back() == 's' ? " hold" : " holds") +
                  " a non-finite value at entry " + std::to_string(i) + " (object " + std::to_string(i/3) + ", component " +
                  std::to_string(i%3) + ")");

  return a.size() / 3;
}


// ============================================================================


void bitvf::internal::check_ids (const std::vector<std::int64_t>& ids, const std::size_t nTracers,
                                 const std::string& fn, const bool fromResult)
{
  if (ids.empty()) return;
  if (ids.size() != nTracers)
    throw Error(fn + ": " + (fromResult ? "the reconstruction's ids" : "ids") + " holds " + std::to_string(ids.size()) + " entries for " +
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
    throw Error((fromResult ? fn + ": the reconstruction's ids hold " : std::string("the ids hold ")) + std::to_string(ids[first]) + " twice, at entries " +
                std::to_string(first) + " and " + std::to_string(second) +
                "; each tracer needs its own ID");
}


// ============================================================================


void bitvf::internal::check_reconstruction_box (const otswap::Result& rec, const std::string& fn)
{
  const std::string malformed = fn + ": the reconstruction is malformed: ";

  if (!rec.box) throw Error(malformed + "it holds no box (otswap::Result::box)");
  const otswap::BoxLimits& box = *rec.box;
  for (std::size_t k = 0; k < 3; ++k) {
    if (!std::isfinite(box.lower[k]) || !std::isfinite(box.upper[k]))
      throw Error(malformed + "its box (otswap::Result::box) has a bound of " +
                  std::to_string(std::isfinite(box.lower[k]) ? box.upper[k] : box.lower[k]) +
                  " on axis " + kAxes[k] + "; every bound must be finite");
    if (box.lower[k] > box.upper[k])
      throw Error(malformed + "its box (otswap::Result::box) spans [" + exact(box.lower[k]) + ", " +
                  exact(box.upper[k]) + "] on axis " + kAxes[k] +
                  "; the lower bound must not be above the upper one");
  }
  if (rec.selection.boxLimits && (rec.selection.boxLimits->lower != box.lower ||
                                  rec.selection.boxLimits->upper != box.upper))
    throw Error(malformed + "its box (otswap::Result::box) differs from its box limits "
                "(otswap::SelectionCounts::boxLimits); a reconstruction with limits uses them as its box");
}


// ============================================================================


void bitvf::internal::check_reconstruction_shape (const otswap::Result& rec, const std::string& fn)
{
  const std::size_t N = rec.nObjects;
  const std::size_t R = rec.nRealizations;
  const std::string malformed = fn + ": the reconstruction is malformed: ";
  if (N == 0) throw Error(malformed + "it holds no objects");
  if (R == 0) throw Error(malformed + "it holds no realizations");
  if (rec.nRealizationsMean == 0) throw Error(malformed + "nRealizationsMean is 0; it must be at least 1");
  if (rec.nRealizationsMean > 1 && R != 1)
    throw Error(malformed + "it holds " + std::to_string(R) + " realizations, each the mean of " +
                std::to_string(rec.nRealizationsMean) + "; a mean holds one realization");
  if (!rec.validRealizations.empty()) {
    if (rec.validRealizations.size() != N)
      throw Error(malformed + "validRealizations holds " + std::to_string(rec.validRealizations.size()) +
                  " entries for " + std::to_string(N) + " objects");
    for (std::size_t i = 0; i < N; ++i)
      if ((std::size_t)rec.validRealizations[i] > R * rec.nRealizationsMean)
        throw Error(malformed + "validRealizations holds " + std::to_string(rec.validRealizations[i]) +
                    " for object " + std::to_string(i) + ", above nRealizations * nRealizationsMean = " +
                    std::to_string(R * rec.nRealizationsMean));
  }

  auto check_size = [&] (const std::size_t size, const std::size_t expected,
                         const std::string& name, const bool allowEmpty) {
    if (size != expected && !(allowEmpty && size == 0))
      throw Error(malformed + name + " holds " + std::to_string(size) + " entries for " +
                  std::to_string(R) + " realizations of " + std::to_string(N) + " objects");
  };
  check_size(rec.tracers.size(), 3 * N, "tracers", false);
  check_size(rec.displacement.size(), 3 * R * N, "displacement", false);
  check_size(rec.valid.size(), R * N, "valid", false);
  check_size(rec.ids.size(), N, "ids", true);

  check_valid(rec.valid, fn, true);
  check_ids(rec.ids, N, fn, true);
}


// ============================================================================


std::size_t bitvf::internal::check_displacement_shape (const std::vector<double>& displacement,
                                                       const std::vector<std::uint8_t>& valid,
                                                       const std::size_t N, const std::string& fn)
{
  if (displacement.empty() || displacement.size() % (3 * N) != 0)
    throw Error(fn + ": displacement holds " + std::to_string(displacement.size()) + " entries for " +
                std::to_string(N) + " tracers; it must hold 3 * nRealizations * " +
                std::to_string(N) + " entries, with a whole nRealizations of at least 1");

  const std::size_t R = displacement.size() / (3 * N);
  if (R > UINT_MAX)
    throw Error("the displacements hold " + std::to_string(R) + " realizations; at most " +
                std::to_string(UINT_MAX) + " are supported");

  if (!valid.empty() && valid.size() != R * N)
    throw Error(fn + ": valid holds " + std::to_string(valid.size()) + " entries for " +
                std::to_string(R) + " realizations of " + std::to_string(N) +
                " tracers; give one per displacement, or none");
  check_valid(valid, fn, false);

  return R;
}


// ============================================================================


bitvf::internal::BoxSearch bitvf::internal::check_box (const std::vector<double>& tracers,
                                                       const std::vector<double>& displacement,
                                                       const std::vector<std::uint8_t>& valid,
                                                       const double mps,
                                                       const Config& config,
                                                       const std::optional<otswap::BoxLimits>& limits,
                                                       const std::vector<std::int64_t>& ids,
                                                       const DivergenceGrid* grid)
{
  const std::string fn = "findVoidsBox";
  const std::size_t N = check_coordinates(tracers, "tracers", false, fn);

  const std::size_t R = check_displacement_shape(displacement, valid, N, fn);

  if (!std::isfinite(mps) || mps <= 0.)
    throw Error("the mean particle separation is " + std::to_string(mps) +
                " in coordinate units; it must be finite and positive");

  if (limits) check_limits(*limits, "the box limits", "");

  check_ids(ids, N, fn, false);

  BoxSearch search;
  search.nTracers = N;
  search.nRealizations = (unsigned)R;
  search.tracers = &tracers;
  search.displacement = &displacement;
  search.valid = &valid;
  search.ids = &ids;
  search.mps = mps;
  search.box = limits;
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

  select_box(search, config, grid, true, fn, false);
  return search;
}


// ============================================================================


bitvf::internal::BoxSearch bitvf::internal::check_box (const otswap::Result& reconstruction,
                                                       const Config& config,
                                                       const DivergenceGrid* grid,
                                                       std::vector<double>& moved)
{
  const otswap::Result& rec = reconstruction;
  const std::string fn = "findVoidsBox";

  if (rec.geometry != otswap::Geometry::Box)
    throw Error(fn + ": the reconstruction is not a box result; use the lightcone void finder");

  const std::size_t N = rec.nObjects;
  const std::size_t R = rec.nRealizations;
  const std::string malformed = fn + ": the reconstruction is malformed: ";
  check_reconstruction_shape(rec, fn);
  if (!rec.outsideBox.empty() && rec.outsideBox.size() != N)
    throw Error(malformed + "outsideBox holds " + std::to_string(rec.outsideBox.size()) + " entries for " +
                std::to_string(R) + " realizations of " + std::to_string(N) + " objects");

  BoxSearch search;
  search.nTracers = N;
  search.nRealizations = (unsigned)R;
  search.tracers = &rec.tracers;
  search.displacement = &rec.displacement;
  search.valid = &rec.valid;
  search.ids = &rec.ids;
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
          throw Error(fn + ": the reconstruction's tracers hold a non-finite value at entry " +
                      std::to_string(3*i+k) + " (object " + std::to_string(i) + ", component " +
                      std::to_string(k) + ")");

  if (grid == nullptr) {
    if (!std::isfinite(rec.mps) || rec.mps <= 0.)
      throw Error(fn + ": the reconstruction's mps is " + std::to_string(rec.mps) +
                  " in coordinate units; it must be finite and positive");
    search.mps = rec.mps;
  }
  else search.mps = grid->mps;

  if (rec.selection.boxLimits) check_limits(*rec.selection.boxLimits, "the reconstruction's box limits", fn);

  check_reconstruction_box(rec, fn);
  const otswap::BoxLimits& box = *rec.box;

  for (std::size_t i = 0; i < N; ++i)
    if (search.inBox[i]) {
      bool movedTracer = false;
      for (std::size_t k = 0; k < 3; ++k) {
        const double x = rec.tracers[3*i+k];
        if (x >= box.lower[k] && x <= box.upper[k]) continue;
        const double bound = x < box.lower[k] ? box.lower[k] : box.upper[k];
        const double tolerance = 1e-8 * std::max(std::abs(x), std::abs(bound));
        if (rec.selection.boxLimits || !(std::abs(x - bound) <= tolerance))
          throw Error(malformed + "tracer " + std::to_string(i) + " lies outside its box "
                      "(otswap::Result::box) and is not flagged in outsideBox");
        if (moved.empty()) moved = rec.tracers;
        moved[3*i+k] = x < box.lower[k] ? box.lower[k] : std::nextafter(box.upper[k], box.lower[k]);
        movedTracer = true;
      }
      if (movedTracer) ++search.tracersMoved;
    }
  if (!moved.empty()) search.tracers = &moved;
  search.box = box;
  search.limits = rec.selection.boxLimits;

  select_box(search, config, grid, false, fn, true);
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

  if (config.nRealizationsMean == 0)
    throw Error("the number of realizations averaged is 0; it must be at least 1 "
                "(1 for displacements of single realizations)");

  if (gridGiven && config.divergenceCellSize)
    throw Error("the divergence cell size is set together with a grid to reuse; leave it empty, "
                "the cells are the grid's");
}


// ============================================================================


void bitvf::internal::check_box_grid (const DivergenceGrid& grid, const double mps,
                                      const bool mpsMustMatch, const std::string& fn, const std::string& name)
{
  const std::string bad = (fn.empty() ? "" : fn + ": ") + name + " ";

  if (grid.geometry != otswap::Geometry::Box)
    throw Error(bad + "is not a box grid");

  if (!std::isfinite(grid.cellSize) || grid.cellSize <= 0.)
    throw Error(bad + "has a cell size of " + std::to_string(grid.cellSize) +
                " in coordinate units; it must be finite and positive");

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
                " in coordinate units; it must be finite and positive");

  if (!std::isfinite(grid.cellSizeMps) || grid.cellSizeMps <= 0.)
    throw Error(bad + "has a cell size of " + std::to_string(grid.cellSizeMps) +
                " mps (cellSizeMps); it must be finite and positive");

  if (grid.cellSize != grid.cellSizeMps * grid.mps)
    throw Error(bad + "has a cell side of " + exact(grid.cellSize) +
                " in coordinate units, but cellSizeMps * mps is " + exact(grid.cellSizeMps * grid.mps) +
                "; they must be equal");

  const bool lightconeSet =
    !std::isnan(grid.distanceMin) || !std::isnan(grid.distanceMax) ||
    !std::isnan(grid.raMin) || !std::isnan(grid.raMax) ||
    !std::isnan(grid.decMin) || !std::isnan(grid.decMax) ||
    !grid.mpsProfile.redshift.empty() || !grid.mpsProfile.mps.empty() ||
    !grid.mpsProfile.count.empty();
  if (lightconeSet)
    throw Error(bad + "is a box grid with a field of a lightcone grid set");

  if (mpsMustMatch && mps != grid.mps)
    throw Error((fn.empty() ? "" : fn + ": ") + "the mean particle separation is " + exact(mps) + " in coordinate units, but " + name + " was built with " +
                exact(grid.mps) + "; they must be equal");
}


// ============================================================================


bitvf::internal::BoxSearch bitvf::internal::check_box_from_divergence (const DivergenceGrid& grid,
                                                                       const std::vector<double>& tracers,
                                                                       const Config& config,
                                                                       const std::vector<std::int64_t>& ids)
{
  const std::string fn = "findVoidsBoxFromDivergence";
  const std::size_t N = check_coordinates(tracers, "tracers", false, fn);
  check_ids(ids, N, fn, false);

  check_config(config, false);
  if (config.divergenceCellSize)
    throw Error("the divergence cell size is set for a search from a divergence; leave it empty, "
                "the cells are the grid's");

  check_box_grid(grid, grid.mps, false, fn, "the divergence grid");
  if (grid.limits) check_limits(*grid.limits, "the grid's box limits", fn);

  const std::size_t C = grid.nCells;
  const std::pair<const std::vector<double>*, const char*> fields[2] = {
    {&grid.divergence, "divergence"}, {&grid.weight, "weight"}};
  for (const auto& field : fields) {
    if (field.first->size() != C)
      throw Error(fn + ": the grid's " + std::string(field.second) + " holds " + std::to_string(field.first->size()) +
                  " entries for " + std::to_string(C) + " cells; a search from a divergence needs one per cell"
                  " (keep it with Config::keepDivergence, or read it from a divergence file)");
    for (std::size_t c = 0; c < C; ++c)
      if (!std::isfinite((*field.first)[c]))
        throw Error(fn + ": the grid's " + std::string(field.second) + " holds a non-finite value at cell " +
                    std::to_string(c));
  }
  if (!grid.nTracers.empty() && grid.nTracers.size() != C)
    throw Error(fn + ": the grid's nTracers holds " + std::to_string(grid.nTracers.size()) + " entries for " +
                std::to_string(C) + " cells; give one per cell, or none");

  BoxSearch search;
  search.nTracers = N;
  search.nRealizations = grid.nRealizations;
  search.tracers = &tracers;
  search.ids = &ids;
  search.mps = grid.mps;
  search.limits = grid.limits;
  search.counts.tracers = N;
  search.sigmaLength = sigma_length(config.smoothing, grid.mps);

  search.inBox.assign(N, 1);
  if (grid.limits)
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t k = 0; k < 3; ++k)
        if (!(tracers[3*i+k] >= grid.limits->lower[k] && tracers[3*i+k] <= grid.limits->upper[k])) {
          search.inBox[i] = 0;
          ++search.counts.tracersOutsideBox;
          break;
        }

  select_tracers(search, &grid);

  if (search.counts.tracersOutsideBox + search.counts.tracersOutsideGrid == N)
    throw Error("no tracer is left to read: of " + std::to_string(N) + " tracers, " +
                std::to_string(search.counts.tracersOutsideBox) + " lie outside the grid's box limits and " +
                std::to_string(search.counts.tracersOutsideGrid) + " outside the grid");

  return search;
}


// ============================================================================


double bitvf::internal::check_box_from_field (const std::vector<double>& centres,
                                              const std::vector<double>& values,
                                              const Config& config,
                                              const std::vector<double>& weights,
                                              const std::optional<double>& mps)
{
  const std::string fn = "findVoidsBoxFromField";
  const std::size_t N = check_coordinates(centres, "centres", false, fn);

  if (values.size() != N)
    throw Error(fn + ": values holds " + std::to_string(values.size()) + " entries for " + std::to_string(N) +
                " centres; it must hold one per centre");
  for (std::size_t c = 0; c < N; ++c)
    if (!std::isfinite(values[c]))
      throw Error("the field values hold a non-finite value at entry " + std::to_string(c));

  if (!weights.empty() && weights.size() != N)
    throw Error(fn + ": weights holds " + std::to_string(weights.size()) + " entries for " + std::to_string(N) +
                " centres; it must hold one per centre, or none");
  for (std::size_t c = 0; c < weights.size(); ++c)
    if (!std::isfinite(weights[c]) || weights[c] < 0.)
      throw Error("the weights hold " + exact(weights[c]) + " at entry " + std::to_string(c) +
                  "; every weight must be finite and non-negative");

  if (mps && (!std::isfinite(*mps) || *mps <= 0.))
    throw Error("the mean particle separation is " + std::to_string(*mps) +
                " in coordinate units; it must be finite and positive");

  check_config(config, false);
  if (config.divergenceCellSize)
    throw Error("the divergence cell size is set for a search on a field; leave it empty, "
                "the cells are the lattice's");
  if (config.keepAssignment)
    throw Error("a search on a field has no tracers to assign to cells; do not ask for the assignment");
  if (config.nRealizationsMean != 1)
    throw Error("the number of realizations averaged is " + std::to_string(config.nRealizationsMean) +
                " for a search on a field, which has no realizations; leave it 1");

  if (!mps && config.smoothing.unit == LengthUnit::Mps && config.smoothing.sigma > 0.)
    throw Error("the smoothing sigma is in mps, but no mps is given; give the mps, or the sigma in "
                "coordinate units (LengthUnit::Coordinate)");

  const double sigmaLength = sigma_length(config.smoothing, mps ? *mps : 1.);

  if (config.smoothing.sigma > 0. && !weights.empty() &&
      std::all_of(weights.begin(), weights.end(), [] (const double w) { return w == 0.; }))
    throw Error("every weight is 0, so no cell can be smoothed; give a positive weight, or set "
                "the smoothing sigma to 0");

  return sigmaLength;
}
