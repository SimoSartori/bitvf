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
 *  @file src/internal.h
 *
 *  @brief Declarations shared by the bitvf translation units. Nothing
 *  here is part of the public API.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef BITVF_INTERNAL_H
#define BITVF_INTERNAL_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "bitvf/VoidFinder.h"

namespace bitvf {
  namespace internal {

    /**
     *  @brief The bitvf version, as "major.minor.patch": the version the
     *  writers record in the files they write.
     */
    extern const char* const version;

    //==========================================================================

    /**
     *  @brief The inputs of one box search after the checks, whichever
     *  overload of findVoidsBox received them.
     *
     *  The arrays are not copied: the pointers refer to the caller's
     *  arguments, which outlive the search. A tracer is read when inBox is 1
     *  and, with a reused grid, outsideGrid is 0; only the valid
     *  displacements of the tracers read enter the search.
     */
    struct BoxSearch {

      std::size_t nTracers = 0;
      unsigned nRealizations = 0;

      /// [tracer][xyz], [realization][tracer][xyz], [realization][tracer];
      /// valid empty means all valid.
      const std::vector<double>* tracers = nullptr;
      const std::vector<double>* displacement = nullptr;
      const std::vector<std::uint8_t>* valid = nullptr;

      /// The identity of each tracer, as given; empty for none.
      const std::vector<std::int64_t>* ids = nullptr;

      /// [tracer]: 0 for a tracer outside the box limits; all 1 without
      /// them.
      std::vector<std::uint8_t> inBox;

      /// [tracer]: 1 for a tracer inside the limits but outside every cell
      /// of the reused grid; empty without a reused grid.
      std::vector<std::uint8_t> outsideGrid;

      /// The randoms kept, inside the limits and the reused grid, [random][xyz];
      /// empty when none were given.
      std::vector<double> randoms;

      /// Whether randoms were given; when not, the displacement endpoints
      /// take their place.
      bool randomsGiven = false;

      /// The mean particle separation of the search, Mpc/h.
      double mps = 0.;

      /// The largest absolute value of a component of the valid
      /// displacements read, Mpc/h; it sets the scale of the fixed-point
      /// accumulation of the divergence.
      double maxDisplacement = 0.;

      /// The box limits, when given.
      std::optional<otswap::BoxLimits> limits;

      /// The counts known after the checks: tracers, tracersOutsideBox,
      /// tracersOutsideGrid, randoms, randomsOutsideBox, randomsOutsideGrid,
      /// displacements.
      Counts counts;
    };

    /**
     *  @brief Check the arguments of the array overload of findVoidsBox
     *  and select what the search reads.
     *
     *  The checks run in a fixed order, so that the same defects always
     *  raise the same error: the arrays, mps, the limits, ids, the
     *  configuration, the reused grid, then the data.
     *
     *  @exception Error on every case listed by findVoidsBox, except those
     *  found by building the grid and by the divergence.
     */
    BoxSearch check_box (const std::vector<double>& tracers,
                         const std::vector<double>& displacement,
                         const std::vector<std::uint8_t>& valid,
                         const std::vector<double>& randoms,
                         double mps,
                         const Config& config,
                         const std::optional<otswap::BoxLimits>& limits,
                         const std::vector<std::int64_t>& ids,
                         const DivergenceGrid* grid);

    /**
     *  @brief Check the arguments of the otswap::Result overload of
     *  findVoidsBox and select what the search reads.
     *
     *  The tracers flagged in otswap::Result::outsideBox are skipped
     *  unread. With a grid, its mps is the search's, and the
     *  reconstruction's is not read.
     *
     *  @exception Error as the overload above, and on a reconstruction that
     *  is not a box result or is malformed.
     */
    BoxSearch check_box (const otswap::Result& reconstruction,
                         const Config& config,
                         const std::vector<double>& randoms,
                         const DivergenceGrid* grid);

    /**
     *  @brief Check the configuration of a search.
     *  @param gridGiven whether a grid is reused, with which
     *  Config::divergenceCellSize must be empty
     *  @exception Error on an invalid field or combination
     */
    void check_config (const Config& config, bool gridGiven);

    /**
     *  @brief Check the fields of a box grid to reuse that define its cells.
     *  @param mps the mps of the search
     *  @param mpsMustMatch whether mps must equal grid.mps, bit for bit
     *  @exception Error on an inconsistent field, or on a different mps
     *  when it must match
     */
    void check_box_grid (const DivergenceGrid& grid, double mps, bool mpsMustMatch);

    //==========================================================================

    /// The most cells a box grid may have along an axis: the largest n with
    /// n^3 below 2^63, so that every cell row fits an int64 index.
    constexpr std::size_t kMaxCellsAxis = 2097151;

    /**
     *  @brief Build the box grid of a search.
     *
     *  The extent is the box limits when given, otherwise the bounding box
     *  of the tracers read and of the randoms kept, or of the endpoints of
     *  the valid displacements read when no randoms were given. It is
     *  padded by 0.1 mps on each side and placed in a cube of cells as
     *  DivergenceGrid describes. Only the fields that define the cells, and
     *  mps, are set; the per-cell arrays are left empty.
     *
     *  @param cellSizeMps the cell side in units of search.mps
     *  @exception Error if the grid would hold more cells than an int64
     *  index can address
     */
    DivergenceGrid box_grid (const BoxSearch& search, double cellSizeMps);

    /**
     *  @brief The row of the cell of a box grid that holds a point.
     *
     *  Every lookup of a point in a box grid goes through this function,
     *  so that a point falls in the same cell everywhere. A point on the
     *  lower face of a cell is in that cell; a point on the upper face of
     *  the grid is outside.
     *
     *  @param x the point, x, y, z
     *  @return the row of the cell, or -1 when the point lies outside every
     *  cell
     */
    std::int64_t box_cell (const DivergenceGrid& grid, const double* x);

    /// Fill DivergenceGrid::centre of a box grid from its defining fields.
    void box_centres (DivergenceGrid& grid);

    //==========================================================================

    /// What the divergence of a box search finds besides the grid fields.
    struct BoxDivergence {

      /// [tracer]: the row of the cell the tracer starts in; -1 for a
      /// tracer not read.
      std::vector<std::int64_t> cell;

      /// [tracer]: 1 when a valid displacement of the tracer entered the
      /// divergence.
      std::vector<std::uint8_t> used;

      /// The tracers with used 1 (Counts::tracersUsed).
      std::size_t nUsed = 0;

      /// The exponent q of the fixed-point accumulation: each term of the
      /// divergence is rounded to a multiple of 2^-q.
      int scaleExponent = 0;

      /// The passes the fill of the empty cells ran.
      unsigned fillPasses = 0;
    };

    /// The face of a cell crossed by a displacement: axis 0, 1 or 2 and
    /// sign -1 (lower face) or +1 (upper face); axis -1 when none is.
    struct BoxFace {
      int axis;
      int sign;
    };

    /**
     *  @brief The face of a box cell through which the segment from x to
     *  x + d crosses.
     *
     *  The lower faces are tried before the upper ones, the axes x, y, z in
     *  that order within each, an axis only when the component of d along
     *  it exceeds 1e-9 in absolute value. A face is crossed when the
     *  segment meets its plane at 0 < t < 1 (x + t d) strictly inside the
     *  face; the first face crossed is returned.
     *
     *  @param cell the cell, (i, j, k)
     *  @param x the start of the segment, x, y, z
     *  @param d the displacement, x, y, z
     */
    BoxFace box_face (const DivergenceGrid& grid, const std::array<std::size_t, 3>& cell,
                      const double* x, const double* d);

    /**
     *  @brief The divergence of the displacement field on a box grid,
     *  before smoothing.
     *
     *  Sets DivergenceGrid::nTracers, nDisplacements, divergence and
     *  weight from the valid displacements read, then fills the empty
     *  cells (box_fill). A displacement that leaves its tracer's cell adds
     *  the outward component, in units of the cell side, through the face
     *  it crosses in the start cell and in the end cell, and counts 1 in
     *  each; one that ends in its tracer's cell adds nothing and counts 2
     *  there. The divergence of a cell is the sum over its count, its
     *  weight the count over the number of realizations.
     *
     *  The sums are accumulated in fixed point, so that the result does not
     *  depend on the number of threads: each term is rounded to a multiple
     *  of 2^-q (BoxDivergence::scaleExponent), with q set by
     *  BoxSearch::maxDisplacement over the cell side and by the largest
     *  count of a cell, so that no sum overflows.
     *
     *  @param search the checked search, from check_box
     *  @param grid a grid with its defining fields set; its per-cell
     *  arrays other than these four are left unchanged
     *  @exception Error if a valid displacement read ends outside every
     *  cell; then if no face of a cell is found for a displacement leaving
     *  it; then if the fill cannot reach every empty cell
     */
    BoxDivergence box_divergence (const BoxSearch& search, DivergenceGrid& grid);

    /**
     *  @brief Fill the empty cells (nDisplacements 0) of a box grid's
     *  divergence and weight, as DivergenceGrid::weight describes.
     *
     *  Reads only nDisplacements, divergence, weight and the fields that
     *  define the cells.
     *
     *  @return the number of passes run, 0 with no empty cell
     *  @exception Error if a pass fills no cell while some are empty
     */
    unsigned box_fill (DivergenceGrid& grid);

  }
}

#endif
