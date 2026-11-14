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
 *  @file src/VoidFinderField.cpp
 *
 *  @brief The box search on a user field, findVoidsBoxFromField: the
 *  lattice of the centres, the field at the basin centres.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

#include "internal.h"

namespace {

  const char* const kAxes[3] = {"x", "y", "z"};

  /// The box lattice of a set of cell centres.
  struct Lattice {

    /// Nodes along x, y, z.
    std::array<std::size_t, 3> n = {0, 0, 0};

    /// The smallest and the largest coordinate along each axis.
    std::array<double, 3> min = {0., 0., 0.};
    std::array<double, 3> max = {0., 0., 0.};

    /// The cell side.
    double s = 0.;

    /// [node]: the row of the centre at node (i, j, k), row
    /// (i * ny + j) * nz + k of the lattice.
    std::vector<std::size_t> row;
  };

  /**
   *  @brief Infer the lattice of the centres, as findVoidsBoxFromField
   *  describes, serially.
   *
   *  On each axis with more than one node, a centre's node is
   *  floor((x - min) / spacing + 0.5) with the spacing of that axis, so
   *  that min and max are nodes.
   *
   *  @param centres finite coordinates, [centre][xyz], at least one centre
   *  @exception Error on every lattice case listed by findVoidsBoxFromField
   */
  Lattice infer_lattice (const std::vector<double>& centres)
  {
    const std::size_t N = centres.size() / 3;
    Lattice L;
    std::array<double, 3> spacing = {0., 0., 0.};

    std::vector<double> x(N);
    for (std::size_t a = 0; a < 3; ++a) {
      for (std::size_t i = 0; i < N; ++i) x[i] = centres[3*i+a];
      std::sort(x.begin(), x.end());
      double largest = 0.;
      for (std::size_t i = 1; i < N; ++i) largest = std::max(largest, x[i] - x[i-1]);
      std::size_t n = 1;
      for (std::size_t i = 1; i < N; ++i)
        if (x[i] - x[i-1] > largest / 2) ++n;
      L.n[a] = n;
      L.min[a] = x[0];
      L.max[a] = x[N-1];
      if (n > 1) {
        spacing[a] = (x[N-1] - x[0]) / (double)(n - 1);
        if (!std::isfinite(spacing[a]))
          throw bitvf::Error(std::string("the centres span [") + std::to_string(x[0]) + ", " +
                             std::to_string(x[N-1]) + "] along " + kAxes[a] +
                             ", too wide a range for the cell side to be finite");
      }
    }
    std::vector<double>().swap(x);

    std::size_t axis = 3;
    for (std::size_t a = 0; a < 3; ++a)
      if (L.n[a] > 1 && (axis == 3 || L.n[a] > L.n[axis])) axis = a;
    if (axis == 3)
      throw bitvf::Error("the centres span one node on every axis: a single cell, whose side cannot be "
                         "inferred; give at least two cells");
    L.s = spacing[axis];

    for (std::size_t a = 0; a < 3; ++a)
      if (L.n[a] > 1 && !(std::abs(spacing[a] - L.s) <= 1e-6 * L.s))
        throw bitvf::Error(std::string("the cells are not cubic: the centres are spaced ") +
                           std::to_string(spacing[a]) + " along " + kAxes[a] + " and " +
                           std::to_string(L.s) + " along " + kAxes[axis] +
                           "; the spacings must agree within 1e-6 relative");

    const std::string shape = "(" + std::to_string(L.n[0]) + ", " + std::to_string(L.n[1]) + ", " +
      std::to_string(L.n[2]) + ")";
    std::size_t nodes = 1;
    bool fits = true;
    for (const std::size_t n : L.n) {
      if (n > N / nodes) fits = false;
      else nodes *= n;
    }
    if (!fits)
      throw bitvf::Error("the centres form a lattice of " + shape + " nodes but are only " + std::to_string(N) +
                         "; every node must hold one centre, and some hold none");

    const std::size_t none = std::numeric_limits<std::size_t>::max();
    L.row.assign(nodes, none);
    for (std::size_t r = 0; r < N; ++r) {
      std::size_t node[3] = {0, 0, 0};
      for (std::size_t a = 0; a < 3; ++a) {
        if (L.n[a] == 1) continue;
        const double t = std::floor((centres[3*r+a] - L.min[a]) / spacing[a] + 0.5);
        const std::size_t i = t <= 0. ? 0 : t >= (double)(L.n[a] - 1) ? L.n[a] - 1 : (std::size_t)t;
        const double offset = std::abs(centres[3*r+a] - (L.min[a] + spacing[a] * (double)i));
        if (!(offset <= 1e-6 * L.s))
          throw bitvf::Error("the centre of cell " + std::to_string(r) + " (from 0) lies " + std::to_string(offset) +
                             " from its node along " + kAxes[a] + ", more than 1e-6 of the cell side " +
                             std::to_string(L.s) + "; the centres must lie on a lattice");
        node[a] = i;
      }
      const std::size_t c = (node[0] * L.n[1] + node[1]) * L.n[2] + node[2];
      if (L.row[c] != none)
        throw bitvf::Error("the centres of cells " + std::to_string(L.row[c]) + " and " + std::to_string(r) +
                           " (from 0) lie at the same node (" + std::to_string(node[0]) + ", " + std::to_string(node[1]) +
                           ", " + std::to_string(node[2]) + ") of the " + shape +
                           " lattice, counted from 0 along each axis; every node must hold one centre");
      L.row[c] = r;
    }

    return L;
  }

}


// ============================================================================


double bitvf::internal::box_field_at (const DivergenceGrid& grid, const std::vector<double>& field, const double* x)
{
  const double s = grid.cellSize;
  std::size_t index[3][2] = {{0, 0}, {0, 0}, {0, 0}};
  double weight[3][2] = {{1., 0.}, {1., 0.}, {1., 0.}};
  std::size_t count[3] = {1, 1, 1};

  for (std::size_t a = 0; a < 3; ++a) {
    const std::size_t n = grid.nCellsAxis[a];
    if (n == 1) continue;
    const double lower = grid.lower[a];
    const auto centre = [lower, s] (const std::size_t i) { return lower + s * (i + 0.5); };
    const double y = std::max(centre(0), std::min(centre(n - 1), x[a]));
    const double t = std::floor((y - lower) / s - 0.5);
    std::size_t i = t <= 0. ? 0 : t >= (double)(n - 2) ? n - 2 : (std::size_t)t;
    if (i > 0 && y < centre(i)) --i;
    else if (i + 2 < n && y >= centre(i + 1)) ++i;
    const double f = y >= centre(i + 1) ? 1. : std::min(1., (y - centre(i)) / s);
    index[a][0] = i;
    index[a][1] = i + 1;
    weight[a][0] = 1. - f;
    weight[a][1] = f;
    count[a] = 2;
  }

  const std::size_t ny = grid.nCellsAxis[1], nz = grid.nCellsAxis[2];
  double value = 0.;
  for (std::size_t p = 0; p < count[0]; ++p)
    for (std::size_t q = 0; q < count[1]; ++q)
      for (std::size_t r = 0; r < count[2]; ++r)
        value += weight[0][p] * weight[1][q] * weight[2][r] *
          field[(index[0][p] * ny + index[1][q]) * nz + index[2][r]];
  return value;
}


// ============================================================================


bitvf::Result bitvf::findVoidsBoxFromField (const std::vector<double>& centres,
                                            const std::vector<double>& values,
                                            const Config& config,
                                            const std::vector<double>& weights,
                                            const std::optional<double> mps)
{
  internal::Report report(config.verbosity, config.log);
  const double sigmaLength = internal::check_box_from_field(centres, values, config, weights, mps);
  const Lattice lattice = infer_lattice(centres);
  const std::int64_t C = (std::int64_t)values.size();
  const double s = lattice.s;
  const double nan = std::numeric_limits<double>::quiet_NaN();

  Result result;
  result.config = config;
  result.geometry = otswap::Geometry::Box;
  result.mps = mps ? *mps : nan;
  if (mps) result.config.divergenceCellSize = s / *mps;

  DivergenceGrid cells;
  cells.geometry = otswap::Geometry::Box;
  cells.nCells = (std::size_t)C;
  cells.nCellsAxis = lattice.n;
  cells.cellSize = s;
  for (std::size_t a = 0; a < 3; ++a) {
    cells.lower[a] = lattice.min[a] - s / 2;
    cells.upper[a] = lattice.max[a] + s / 2;
  }
  cells.mps = result.mps;
  cells.cellSizeMps = mps ? s / *mps : nan;

  double wMax = weights.empty() ? 1. : 0.;
  for (const double w : weights) wMax = std::max(wMax, w);

  cells.divergence.resize((std::size_t)C);
  cells.weight.resize((std::size_t)C);
#pragma omp parallel for schedule(static)
  for (std::int64_t c = 0; c < C; ++c) {
    const std::size_t r = lattice.row[(std::size_t)c];
    cells.divergence[(std::size_t)c] = values[r];
    cells.weight[(std::size_t)c] = weights.empty() ? 1. : weights[r];
  }
  if (report.detailed()) {
    std::ostringstream out = internal::Report::line();
    out << "lattice inferred: ";
    internal::Report::shape(out, cells.nCellsAxis);
    out << " of side " << s << " in coordinate units on ";
    internal::Report::box(out, cells.lower, cells.upper);
    report.step(out);

    std::ostringstream given = internal::Report::line();
    if (mps) given << "mean particle separation " << *mps << " in coordinate units (given)";
    else given << "no mean particle separation given";
    report.write(given.str());
  }

  internal::box_smooth(cells, sigmaLength, 1e-12 * wMax);
  if (report.detailed()) {
    std::ostringstream out = internal::Report::line();
    if (sigmaLength > 0.) {
      out << "smoothed with sigma " << sigmaLength << " in coordinate units";
      if (mps) out << " (" << sigmaLength / *mps << " mps)";
      report.step(out);
    }
    else {
      report.lap();
      report.write("not smoothed (sigma 0)");
    }
  }

  internal::BoxBasins basins;
  internal::box_watershed(cells, basins);
  if (report.detailed()) {
    std::ostringstream out = internal::Report::line();
    out << "watershed into ";
    internal::Report::basins(out, basins.voids.seed.size(), basins.overdense.seed.size());
    report.step(out);
  }
  if (basins.voids.seed.empty())
    throw Error("no void found: the field is 0 or positive in every cell");
  internal::box_basins(nullptr, nullptr, cells, basins, result.voids, result.overdense);

  for (BasinCatalog* cat : {&result.voids, &result.overdense}) {
    const std::int64_t B = (std::int64_t)cat->nBasins;
    cat->fieldCentre.resize((std::size_t)B);
#pragma omp parallel for schedule(static)
    for (std::int64_t b = 0; b < B; ++b)
      cat->fieldCentre[(std::size_t)b] =
        internal::box_field_at(cells, cells.divergenceSmoothed.empty() ? cells.divergence : cells.divergenceSmoothed,
                               &cat->centre[3 * (std::size_t)b]);
  }
  if (report.detailed()) {
    std::ostringstream out = internal::Report::line();
    out << "catalogues of the basins";
    report.step(out);
  }

  const std::array<std::size_t, 3> shape = cells.nCellsAxis;
  if (config.keepDivergence) {
    internal::box_centres(cells);
    cells.volume.assign(cells.nCells, s * s * s);
    result.divergence = std::move(cells);
    if (report.detailed()) {
      std::ostringstream out = internal::Report::line();
      out << "kept the divergence grid";
      report.step(out);
    }
  }

  result.elapsedSeconds = report.elapsed();
  if (report.normal()) {
    std::ostringstream out = internal::Report::line();
    out << "findVoidsBoxFromField: ";
    internal::Report::basins(out, result.voids.nBasins, result.overdense.nBasins);
    out << " on ";
    internal::Report::shape(out, shape);
    out << ' ';
    internal::Report::time(out, result.elapsedSeconds);
    report.write(out.str());
  }
  return result;
}
