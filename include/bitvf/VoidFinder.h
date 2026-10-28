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
 *  @file include/bitvf/VoidFinder.h
 *
 *  @brief The BitVF void finder: watershed voids and overdense basins
 *  in the divergence of an optimal transport displacement field.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef BITVF_VOIDFINDER_H
#define BITVF_VOIDFINDER_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "otswap/OT.h"

namespace bitvf {

  // ==========================================================================
  // Errors
  // ==========================================================================

  /// Base class of every exception bitvf throws on invalid input or
  /// failure. std::bad_alloc is not converted, and propagates unchanged.
  /// An error raised by otswap inside a bitvf call, an otswap::Error or an
  /// exception derived from it, reaches the caller as an Error whose
  /// message is "otswap: " followed by its own, with the original nested
  /// in it (std::throw_with_nested), so std::rethrow_if_nested recovers
  /// it. Catching Error therefore catches every error of a bitvf call.
  ///
  /// Error is not derived from otswap::Error: the two are unrelated, and a
  /// handler for otswap::Error does not catch it.
  ///
  /// An Error records where it was constructed: the source file, the
  /// function and the line. The location is not part of what(), which
  /// holds the message alone. An error rethrown with context, after a
  /// parallel region or from an otswap call, carries the location of the
  /// rethrow.
  class Error : public std::runtime_error {
  public:

    /// @param what the message, returned by what()
    /// @param file,function,line the location; left to their defaults,
    /// they are those of the expression that constructs the Error
    explicit Error (const std::string& what,
                    const char* file = __builtin_FILE(),
                    const char* function = __builtin_FUNCTION(),
                    int line = __builtin_LINE())
      : std::runtime_error(what), m_file(file), m_function(function), m_line(line) {}

    /// Source file in which the Error was constructed, as the compiler
    /// named it.
    const char* file () const noexcept { return m_file; }

    /// Name of the function in which the Error was constructed.
    const char* function () const noexcept { return m_function; }

    /// Line at which the Error was constructed.
    int line () const noexcept { return m_line; }

  private:

    const char* m_file;
    const char* m_function;
    int m_line;
  };

  // ==========================================================================
  // Array layout
  // ==========================================================================

  /**
   *  @defgroup array_layout Array layout
   *
   *  @brief How the arrays of the void finder are laid out.
   *
   *  The conventions are otswap's (see otswap's coordinate layout). Every
   *  array is flat and row-major. Cartesian arrays hold [object][xyz] in
   *  coordinate units; in a lightcone these are Mpc/h when the distance
   *  table comes from otswap's cosmology. Sky arrays hold [object][ra,
   *  dec, z], angles in radians; displacements are
   *  [realization][object][xyz] and their validity
   *  [realization][object], uint8. The arrays of a result:
   *
   *      [basin][3]   BasinCatalog::centre, centreSky
   *      [basin]      BasinCatalog::radius, radiusExt, deltaAv, delta2Mps,
   *                   contamination (double); edgeFlag (uint8);
   *                   nTracers (size_t)
   *      [cell][3]    DivergenceGrid::centre, centreSky
   *      [cell]       DivergenceGrid::divergence, divergenceUnsmoothed,
   *                   volume, weight, cellRaMin .. cellDistanceMax
   *                   (double); nDisplacements, nTracers (uint64);
   *                   voidIndex, overdenseIndex (int64); masked (uint8)
   *      [tracer]     Result::ids (int64); Result::uncorrected,
   *                   outsideGrid, outsideCut (uint8);
   *                   TracerAssignment::cellIndex, voidIndex,
   *                   overdenseIndex (int64), distance (double), inMask,
   *                   inBox, used (uint8)
   *
   *  Flags are uint8, 0 or 1. Indices of tracers are their input rows;
   *  indices of cells and basins are their rows in the result, and the
   *  fields named cellIndex, voidIndex and overdenseIndex hold them. A
   *  tracer's identity is its ID (Result::ids), never an index. A value
   *  that was not computed is NaN; a basin or cell index that does not
   *  apply is -1.
   */

  // ==========================================================================
  // Configuration
  // ==========================================================================

  /// Unit of a length the user chooses.
  enum class LengthUnit {
    /// Multiples of the mean particle separation: the box mps, or in a
    /// lightcone mps(z) at each cell's redshift.
    Mps,

    /// The same length unit as the tracer coordinates.
    Coordinate
  };

  /**
   *  @brief Gaussian smoothing of the divergence field before the
   *  watershed.
   *
   *  Each cell takes the weighted mean of the divergence of the cells
   *  around it, itself included. The weight of a cell is its tracer
   *  weight (DivergenceGrid::weight) times exp(-r^2 / (2 sigma^2)), r the
   *  distance between the two cell centres. A cell whose weights sum to
   *  1e-12 or less takes 0. A larger sigma merges small basins into
   *  larger ones; sigma = 0 disables the smoothing.
   *
   *  Box: the cells around are those within floor(3.5 sigma / s) cells
   *  along each axis, s the cell side, cut at the faces of the grid: a
   *  cube, not a sphere. Below sigma = s / 3.5 the cube holds the cell
   *  alone. In Mps, sigma is multiplied by the box mps.
   *
   *  Lightcone: the cells around are those whose centre lies closer than
   *  4 sigma, and each weight is also multiplied by the cell's volume.
   *  Masked cells (DivergenceGrid::masked) enter no mean and take 0. In
   *  Mps, sigma is multiplied by mps(z) at the redshift of the centre of
   *  the shell of the cell being smoothed.
   */
  struct Smoothing {

    /// Width of the Gaussian, in unit; must be finite and non-negative.
    /// 1 mps is a reasonable starting value. A positive sigma must be
    /// between about 1e-162 and 1e154 in coordinate units; use 0 to
    /// disable the smoothing.
    double sigma = 1.;

    /// Unit of sigma.
    LengthUnit unit = LengthUnit::Mps;
  };

  /**
   *  @brief Free parameters of the void finder.
   *
   *  Nothing here depends on the geometry: the inputs proper to a
   *  lightcone (mask, distance table, mps(z) profile) are arguments of
   *  the lightcone functions. Voids and overdense basins are always both
   *  found.
   */
  struct Config {

    /// Smoothing of the divergence field.
    Smoothing smoothing;

    /// Side of the divergence-grid cell, in units of the mean particle
    /// separation: the box mps, or in a lightcone mps(z) at each shell.
    /// Empty: (20 / nRealizations)^(1/3), which is 1 at 20 realizations,
    /// the same rule in both geometries. Must be positive when set, and
    /// left empty when a grid is reused. Result::config holds the value
    /// used. Unrelated to otswap::Config::cellSize, the cell of the
    /// neighbour search of the reconstruction, which affects speed only:
    /// this one changes the result.
    std::optional<double> divergenceCellSize;

    /// Fill Result::divergence with the grid and the divergence field.
    /// Off by default: it holds one entry per cell in several arrays.
    bool keepDivergence = false;

    /// Also fill DivergenceGrid::divergenceUnsmoothed, the field before
    /// smoothing. Requires keepDivergence.
    bool keepUnsmoothedDivergence = false;

    /// Fill Result::assignment, where each tracer ended up. Off by
    /// default: it holds one entry per tracer in several arrays.
    bool keepAssignment = false;

    /// Realizations averaged into each displacement given: 1 for the
    /// displacements of single realizations; for the displacements of a
    /// mean displacement file, the NREC of that file. A record only,
    /// written as the keyword NRECMEAN of bitvf's files; no computation
    /// reads it. Must be at least 1.
    unsigned nRealizationsMean = 1;

    /// What the call reports, at otswap's levels; each line starts with
    /// "bitvf: " and ends in '\n'. Every number reported is also in
    /// Result::counts.
    otswap::Verbosity verbosity = otswap::Verbosity::Normal;

    /// Where the report is written; null writes nothing, as Silent.
    std::ostream* log = &std::clog;
  };

  // ==========================================================================
  // Results
  // ==========================================================================

  /**
   *  @brief One catalogue of watershed basins: the voids, or the
   *  overdense basins.
   *
   *  The columns are the same for both kinds and in both geometries.
   *  Row b is basin b, numbered from 0 in the order of the index of its
   *  seed cell, a local minimum of the divergence (smoothed when a
   *  smoothing is set): below 0 for a void, at or above 0 for an
   *  overdense basin. In a box every cell belongs to exactly one basin, a
   *  void or an overdense basin. The order does not depend on the number
   *  of threads.
   *
   *  The core of a void is the part of the void where the divergence is
   *  below 0; the extended void is the whole basin. In a box the core and
   *  the extended basin of an overdense basin coincide, since every cell
   *  of a basin has a divergence at or above that of its seed.
   *
   *  Each cell flows to the neighbour where the divergence drops the most
   *  (the largest decrease, not divided by the distance), among the cells
   *  that share a face, an edge or a corner with it, when that value is
   *  strictly below its own; among equal smallest values, the neighbour
   *  with the lowest cell index. A cell with no strictly smaller
   *  neighbour is a seed, so two adjacent cells of equal value with no
   *  smaller neighbour are both seeds. A basin is the seed and every cell
   *  whose flow reaches it. In a box the grid is not periodic: a cell on
   *  a face has fewer neighbours.
   */
  struct BasinCatalog {

    std::size_t nBasins = 0;

    /// Cartesian centre, in coordinate units: the seed cell's centre x0,
    /// shifted by -sum_n div_n (x_n - x0) / sum_n |div_n| over the seed's
    /// neighbouring cells n, with div the smoothed divergence, and kept
    /// inside the seed cell. Flat, [basin][xyz].
    std::vector<double> centre;

    /// Sky coordinates of centre: right ascension, declination, redshift;
    /// radians. Flat, [basin][3] in a lightcone; empty in a box.
    std::vector<double> centreSky;

    /// Radius of the sphere with the volume of the core, in coordinate units.
    std::vector<double> radius;

    /// Radius of the sphere with the volume of the extended basin, in
    /// coordinate units.
    std::vector<double> radiusExt;

    /// Mean density contrast of the tracers in the core: the tracers
    /// starting in its cells (nTracers) over its volume, divided by the
    /// mean density mps^-3 (mps(z)^-3 in a lightcone), minus 1.
    std::vector<double> deltaAv;

    /// Density contrast within 2 mps of the centre, the sphere closed. In
    /// a box, n_t over its volume, divided by the mean density mps^-3,
    /// minus 1, with n_t the tracers the search reads (inside the limits
    /// and the reused grid) in the sphere. In a lightcone, within 2 mps(z)
    /// at the centre's redshift, n_t / (n_r N_t / N_r) - 1, with n_t and
    /// n_r the tracers and randoms in the sphere and N_t / N_r one ratio of
    /// the tracers to the randoms over the whole search (Counts::randoms).
    std::vector<double> delta2Mps;

    /// Contact with the boundary: 0 none; 1 only the extended basin
    /// touches it; 2 the core already touches it. The boundary is, in a
    /// box, the faces of the grid, with or without box limits; in a
    /// lightcone, the masked cells, the cells outside the survey and the
    /// borders of the grid. In a box, a basin touches the boundary when one
    /// of its cells lies on a face of the grid.
    std::vector<std::uint8_t> edgeFlag;

    /// Fraction of the basin's surface that lies on the boundary of
    /// edgeFlag. In a box, the faces of its cells on the faces of the grid,
    /// over those faces plus the faces its cells share with cells of other
    /// basins; 1 for a basin that fills the grid.
    std::vector<double> contamination;

    /// Tracers starting in the cells of the core: the count behind
    /// deltaAv.
    std::vector<std::size_t> nTracers;
  };

  /**
   *  @brief The grid of the void finder and the divergence field on it,
   *  one row per cell.
   *
   *  An output: filled in Result::divergence when Config::keepDivergence
   *  is set, or read from a divergence file. The void finder takes one
   *  back as a grid to reuse. A reuse reads only the fields that define
   *  the cells, marked so below, cellSizeMps, and mps or mpsProfile; it
   *  checks that they are consistent, and throws Error otherwise. The
   *  per-cell results of the search that filled the grid (divergence,
   *  weights, counts, basin indices, cell limits) and its box limits are
   *  ignored, whatever their size.
   *  Cell c is row c; basins are ordered by the index of their seed cell.
   *
   *  Box. Cubic cells of side cellSize: cell (i, j, k) spans
   *  lower + (i, j, k) * cellSize to lower + (i + 1, j + 1, k + 1) *
   *  cellSize on the three axes, and is row (i * nY + j) * nZ + k, with
   *  (nX, nY, nZ) = nCellsAxis. The grid is built from the box limits of
   *  the call, or without them from the bounding box of the tracers read
   *  and of the endpoints of their valid displacements; from a
   *  reconstruction, from its box (otswap::Result::box). It is padded by
   *  0.1 mps on each side. With E the padded extent
   *  and s the cell side, an axis needs floor(E / s + 1) cells; the
   *  largest count is used on all three axes, so the grid is a cube of
   *  cells even when the box is not. The padded box is then placed in the
   *  cells: on each axis, with n the count and S = n s - E the slack, its
   *  lower face moves down by S / 2, which is lower, and its upper face up
   *  by S / 4, which is upper. With limits the grid overhangs them.
   *
   *  Lightcone. Cells are numbered by shell, then declination ring, then
   *  right-ascension sector. The grid is rebuilt from its construction
   *  parameters (distanceMin, distanceMax, raMin, raMax, decMin, decMax,
   *  cellSizeMps), its mpsProfile and the distance table of the call,
   *  through which each shell's width is evaluated: it is the same grid
   *  only with the same distance table, which nothing checks. The shells
   *  start at a central sphere or at a first edge just below distanceMin,
   *  and are cellSizeMps * mps(z) wide at their centres until past
   *  distanceMax; each is divided into declination rings over the
   *  footprint, with polar caps when the footprint comes near a pole, and
   *  each ring into right-ascension sectors over the footprint, or over
   *  the full circle when the footprint nearly spans it.
   *
   *  The values that define the cells are reproduced bit for bit when the
   *  grid is passed back in memory. A divergence file stores them at 17
   *  significant digits, angles in degrees, so a grid read from one
   *  equals the built one up to the rounding of its right-ascension and
   *  declination edges in a lightcone: a tracer within rounding of such an
   *  edge can change cell.
   */
  struct DivergenceGrid {

    /// The geometry of the grid.
    otswap::Geometry geometry = otswap::Geometry::Box;

    /// Number of cells.
    std::size_t nCells = 0;

    /// Box, defines the cells: cell side, in coordinate units,
    /// cellSizeMps times mps. NaN in a lightcone.
    double cellSize = std::numeric_limits<double>::quiet_NaN();

    /// Box, defines the cells: x, y, z of the lower grid face, in
    /// coordinate units. NaN in a lightcone.
    std::array<double, 3> lower = {std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN()};

    /// Box, defines the cells: x, y, z of the upper face of the padded box
    /// as placed, in coordinate units. It lies below the last cell face,
    /// lower + nCellsAxis * cellSize; it is stored as built and never
    /// recomputed from the other fields. NaN in a lightcone.
    std::array<double, 3> upper = {std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN(),
                                   std::numeric_limits<double>::quiet_NaN()};

    /// Box, defines the cells: cells along x, y, z, the same count on the
    /// three axes. 0 in a lightcone.
    std::array<std::size_t, 3> nCellsAxis = {0, 0, 0};

    /// Box: the box limits of the search that filled the grid; the tracers
    /// outside them were skipped. Empty without limits and in a lightcone.
    /// Not part of the cells: a search reusing the grid does not read it,
    /// its own limits apply.
    std::optional<otswap::BoxLimits> limits;

    /// Lightcone, defines the cells: comoving distance range of the
    /// tracers and randoms the grid covers, in coordinate units. NaN in a
    /// box.
    double distanceMin = std::numeric_limits<double>::quiet_NaN();
    double distanceMax = std::numeric_limits<double>::quiet_NaN();

    /// Lightcone, defines the cells: right-ascension range of the
    /// footprint, radians, unwrapped at its largest gap, so raMin may be
    /// negative. NaN in a box.
    double raMin = std::numeric_limits<double>::quiet_NaN();
    double raMax = std::numeric_limits<double>::quiet_NaN();

    /// Lightcone, defines the cells: declination range of the footprint,
    /// radians. NaN in a box.
    double decMin = std::numeric_limits<double>::quiet_NaN();
    double decMax = std::numeric_limits<double>::quiet_NaN();

    /// Config::divergenceCellSize used. In a box, in units of mps, a record
    /// of the choice beside cellSize, which is cellSizeMps * mps as built,
    /// and which a search reusing the grid reports in Result::config; in a
    /// lightcone, in units of mps(z), and it defines the cells.
    double cellSizeMps = std::numeric_limits<double>::quiet_NaN();

    /// Box: the mean particle separation the grid was built with, in
    /// coordinate units; the smoothing in mps and the delta2Mps radius of
    /// a search that reuses the grid follow it. NaN in a lightcone.
    double mps = std::numeric_limits<double>::quiet_NaN();

    /// Lightcone: the mps(z) profile the grid was built with; the shells,
    /// the smoothing in mps and the delta2Mps radius of a search that
    /// reuses the grid follow it. Empty in a box.
    otswap::MpsProfile mpsProfile;

    /// Lightcone, descriptive: the limits of each cell as built, radians
    /// and coordinate units; never read by a reuse. Right-ascension edges
    /// are continuous across RA = 0 within a ring; a ring of a single
    /// sector spans 2 pi. The central sphere, when present, has distance
    /// limits [0, r] and the full sky; a polar cap reaches declination
    /// pi/2 (or -pi/2) over the full circle. Flat, [cell]; empty in a box.
    std::vector<double> cellRaMin, cellRaMax;
    std::vector<double> cellDecMin, cellDecMax;
    std::vector<double> cellDistanceMin, cellDistanceMax;

    /// Cartesian centre of each cell, in coordinate units. Flat, [cell][xyz].
    std::vector<double> centre;

    /// Sky coordinates of each cell centre: right ascension, declination,
    /// redshift; radians. Flat, [cell][3] in a lightcone; empty in a box.
    std::vector<double> centreSky;

    /// The divergence the watershed ran on, after smoothing.
    std::vector<double> divergence;

    /// The divergence before smoothing, when
    /// Config::keepUnsmoothedDivergence is set; empty otherwise.
    std::vector<double> divergenceUnsmoothed;

    /// Volume of each cell, in coordinate units cubed: uniform in a box.
    std::vector<double> volume;

    /// Tracer weight of each cell, the weight of the smoothing:
    /// nDisplacements / nRealizations where a displacement contributed.
    /// A cell that no displacement contributed to (nDisplacements 0) is
    /// empty, and the fill gives it the divergence and the weight of nearby
    /// cells that hold a value, a contribution or a value filled in an
    /// earlier pass, each mean weighted by 1 / r^2, r the distance between
    /// the cell centres (times the neighbour's volume in a lightcone). The
    /// fill runs in passes, each reading the values left by the previous
    /// one.
    ///
    /// Box: the neighbours within one cell along each axis, or within two
    /// when fewer than three of those hold a value; the passes repeat
    /// until no cell is empty (an Error when a pass fills none).
    ///
    /// Lightcone: the adjacent cells that are not masked, when there are
    /// at least three; the passes repeat until one fills no cell. A cell
    /// left empty keeps a divergence and a weight of 0, and a masked cell
    /// is never filled.
    std::vector<double> weight;

    /// Displacement ends counted in each cell, before the fill: a valid
    /// displacement that leaves its tracer's cell counts 1 in that cell and
    /// 1 in the cell it ends in; one that ends in its tracer's cell counts 2
    /// there. In a lightcone, displacements rejected by the mask count
    /// nowhere (Counts::displacementsRejected).
    std::vector<std::uint64_t> nDisplacements;

    /// Tracers starting in each cell.
    std::vector<std::uint64_t> nTracers;

    /// Row of the void each cell belongs to in Result::voids; -1 none.
    std::vector<std::int64_t> voidIndex;

    /// Row of the overdense basin each cell belongs to in
    /// Result::overdense; -1 none.
    std::vector<std::int64_t> overdenseIndex;

    /// Lightcone: 1 when the cell centre falls on an unobserved pixel of
    /// the mask of the call. Empty in a box.
    std::vector<std::uint8_t> masked;
  };

  /**
   *  @brief Where each input tracer ended up, one row per input row.
   *  Filled in Result::assignment when Config::keepAssignment is set.
   */
  struct TracerAssignment {

    /// Row of the cell the tracer starts in; -1 when outside the grid or
    /// left out.
    std::vector<std::int64_t> cellIndex;

    /// Row of the tracer's void in Result::voids; -1 none.
    std::vector<std::int64_t> voidIndex;

    /// Row of the tracer's overdense basin in Result::overdense; -1 none.
    std::vector<std::int64_t> overdenseIndex;

    /// Distance from the centre of the tracer's basin, its void or its
    /// overdense basin, in coordinate units; NaN for a tracer not in a
    /// basin (not read by the search).
    std::vector<double> distance;

    /// 0 when the tracer falls on an unobserved pixel of the mask; all 1 in
    /// a box.
    std::vector<std::uint8_t> inMask;

    /// 0 when the tracer lies outside the box limits, and so was skipped;
    /// all 1 without limits and in a lightcone.
    std::vector<std::uint8_t> inBox;

    /// 1 when at least one valid displacement of the tracer entered the
    /// divergence.
    std::vector<std::uint8_t> used;
  };

  /// What the void finder counted. Every number it reports is here.
  struct Counts {

    std::size_t tracers = 0;             ///< tracers given
    std::size_t tracersUsed = 0;         ///< of them, with at least one displacement used

    /// Real-space pass of a pipeline with redshift-space correction: left
    /// unmoved by the correction, and excluded (Result::uncorrected).
    std::size_t tracersUncorrected = 0;

    /// Box with limits: outside them, and skipped (TracerAssignment::inBox).
    std::size_t tracersOutsideBox = 0;

    /// Reused grid: outside it, and excluded (Result::outsideGrid). In the
    /// real-space pass of a pipeline: corrected and, in a lightcone, inside
    /// the redshift cut, but outside the first pass's grid.
    std::size_t tracersOutsideGrid = 0;

    /// Real-space pass of a lightcone pipeline with redshift-space
    /// correction: corrected, but moved outside the redshift cut, and
    /// excluded (Result::outsideCut).
    std::size_t tracersOutsideCut = 0;

    /// Lightcone: N_r of delta2Mps, the randoms given and kept. 0 in a box.
    std::size_t randoms = 0;

    std::size_t randomsOutsideGrid = 0;  ///< lightcone, reused grid: randoms outside it, removed
    std::size_t displacements = 0;       ///< valid displacements given

    /// Lightcone: valid displacements rejected because their endpoint, or a
    /// point of their path, falls on a masked cell or outside the survey.
    /// 0 in a box.
    std::size_t displacementsRejected = 0;
  };

  /**
   *  @brief Output of one void search.
   *
   *  Given the same inputs, every field but elapsedSeconds is the same bit
   *  for bit whatever the number of threads.
   */
  struct Result {

    std::size_t nTracers = 0;
    unsigned    nRealizations = 0;

    /// The geometry of the search.
    otswap::Geometry geometry = otswap::Geometry::Box;

    /// The configuration the search ran with, divergenceCellSize set to
    /// the value used.
    Config config;

    /// The identity of each input row, as given: any 64-bit integer,
    /// negative ones included, each row its own. Empty means no IDs, and
    /// the writers then write each row's index, from 0. Flat, [tracer];
    /// size nTracers or 0. No computation reads it.
    std::vector<std::int64_t> ids;

    /// The seed bitvf derived and used to draw the randoms of a box
    /// pipeline given none, which alone reproduces them; empty when bitvf
    /// drew no randoms.
    std::optional<unsigned> randomSeed;

    /// The voids.
    BasinCatalog voids;

    /// The overdense basins.
    BasinCatalog overdense;

    /// The grid and the divergence field, when Config::keepDivergence is
    /// set.
    std::optional<DivergenceGrid> divergence;

    /// Where each tracer ended up, when Config::keepAssignment is set.
    std::optional<TracerAssignment> assignment;

    /// Real-space pass of a pipeline with redshift-space correction: 1 for
    /// a tracer the correction left unmoved, excluded from the second
    /// reconstruction and from this search, row kept. Flat, [tracer]; empty
    /// elsewhere, which reads as all 0.
    std::vector<std::uint8_t> uncorrected;

    /// 1 for a tracer outside the reused grid, excluded from this search,
    /// row kept. In the real-space pass of a pipeline: corrected and, in a
    /// lightcone, inside the redshift cut, but outside the first pass's
    /// grid, so excluded from the second reconstruction too. Flat,
    /// [tracer]; filled only when a grid was reused, empty elsewhere, which
    /// reads as all 0.
    std::vector<std::uint8_t> outsideGrid;

    /// Real-space pass of a lightcone pipeline with redshift-space
    /// correction: 1 for a tracer the correction moved outside the
    /// redshift cut, excluded from the second reconstruction and from this
    /// search, row kept. Flat, [tracer]; empty elsewhere, which reads as
    /// all 0. A tracer carries at most one of uncorrected, outsideCut and
    /// outsideGrid, checked in that order.
    std::vector<std::uint8_t> outsideCut;

    /// Box: the box limits of the search, the limits given or
    /// otswap::SelectionCounts::boxLimits; the tracers outside them were
    /// skipped. Empty without limits and in a lightcone.
    std::optional<otswap::BoxLimits> limits;

    /// Box: the mean particle separation the search used, in coordinate
    /// units: the one given, otswap::Result::mps, or a reused grid's. NaN
    /// in a lightcone.
    double mps = std::numeric_limits<double>::quiet_NaN();

    /// Lightcone: the mps(z) profile the search used, a reused grid's when
    /// one was given. Empty in a box.
    otswap::MpsProfile mpsProfile;

    /// Lightcone: the distance table of the search; coordinates are in
    /// Mpc/h when it comes from otswap's cosmology. Null in a box.
    std::shared_ptr<const otswap::DistanceTable> distances;

    /// What the search counted.
    Counts counts;

    /// Wall time of the call, in seconds. Two runs on the same inputs give
    /// the same result but for this field.
    double elapsedSeconds = 0.;
  };

  // ==========================================================================
  // Void finder
  // ==========================================================================

  /**
   *  @brief Find voids and overdense basins in a box.
   *
   *  The divergence of the displacement field is computed on a grid of
   *  cubic cells (DivergenceGrid), smoothed (Config::smoothing) and
   *  divided by a watershed into basins around its local minima: the
   *  voids, below 0, and the overdense basins, at or above 0
   *  (BasinCatalog).
   *
   *  @param tracers Cartesian coordinates of the tracers, the start of
   *  their displacements, 3 * nTracers entries.
   *
   *  @param displacement displacements, 3 * nRealizations * nTracers
   *  entries, [realization][tracer][xyz]; nRealizations is their size over
   *  that of tracers, at least 1. A tracer plus its displacement is the
   *  random matched to it. Every valid displacement of a tracer the search
   *  reads must be finite; the others, and those of the tracers skipped
   *  outside the limits, are not read.
   *
   *  @param valid validity of each displacement, nRealizations * nTracers
   *  entries, [realization][tracer], each 0 or 1; or empty, which means
   *  all valid.
   *
   *  @param mps mean particle separation, in coordinate units; must be
   *  positive. Use otswap::Result::mps of the reconstruction, or the
   *  overload below. With grid, it must equal grid->mps.
   *
   *  @param config the parameters of the search.
   *
   *  @param limits the box limits, or none. Given, the tracers outside
   *  them keep their row but are skipped (TracerAssignment::inBox,
   *  Counts::tracersOutsideBox), and the grid is built on them. A point on
   *  a face is inside. Each bound must be finite, with lower[k] <
   *  upper[k]. Without them, the grid covers the tracers and the endpoints
   *  of their valid displacements, tracer plus displacement, in every
   *  realization, so that displacements from any source can be searched.
   *  With them, pass the limits the reconstruction used, or a valid
   *  displacement ending outside the grid is an error.
   *
   *  @param ids the identity of each tracer, copied into Result::ids: one
   *  per tracer, no value twice; or empty.
   *
   *  @param grid a grid to reuse, or null to build one. Given, its cells
   *  are used, and its mps sets the smoothing in mps and the delta2Mps
   *  radius; densities come from this search's own counts. Tracers outside
   *  every cell are excluded (Result::outsideGrid,
   *  Counts::tracersOutsideGrid). Result::divergence, when kept, has the
   *  reused cells with this search's values.
   *
   *  @exception Error
   *  - on the arrays: if tracers is empty or its size is not a multiple of
   *    three; if the size of displacement is not 3 * nRealizations *
   *    nTracers for a whole nRealizations of at least 1 that fits in an
   *    unsigned; if valid is neither empty nor of size nRealizations *
   *    nTracers, or holds a value other than 0 and 1; if a coordinate of
   *    tracers, or a component of a displacement read, is not finite; if
   *    ids is neither empty nor of size nTracers, or holds a value twice;
   *  - if mps is not finite and positive; if a bound of the limits is not
   *    finite, or a lower bound is not below the upper one;
   *  - on config: if smoothing.sigma is not finite and non-negative; if
   *    divergenceCellSize is set and not finite and positive; if
   *    nRealizationsMean is 0; if keepUnsmoothedDivergence is set without keepDivergence; if
   *    smoothing.sigma is positive but, in coordinate units (times mps in
   *    Mps), below about 1e-162 or above about 1e154, too small or too
   *    large to build the kernel;
   *  - on grid: if it is not a box grid; if a field that defines the cells
   *    is inconsistent: cellSize, mps or cellSizeMps not finite and
   *    positive, cellSize not equal to cellSizeMps * mps bit for bit, a
   *    bound of lower or upper not finite or lower[k] not below upper[k], a
   *    count of nCellsAxis 0 or the three counts different, nCells not
   *    their cube, a field of the lightcone set; if
   *    config.divergenceCellSize is set;
   *    if mps differs from grid->mps, compared bit for bit;
   *  - on the data: if the search reads more than 4294967295 tracers, the
   *    most the count of delta2Mps supports; if no valid displacement is
   *    left to read; with limits
   *    or a reused grid, if a valid displacement read ends outside every
   *    cell of the grid, which means the arrays do not come from one
   *    reconstruction in that box; if the
   *    grid would hold more cells than an int64 index can address;
   *    if empty cells remain that the fill of the divergence (see
   *    DivergenceGrid::weight) cannot reach; if the divergence is 0 or
   *    positive in every cell, so that no void can be found.
   *
   *  @return the voids and overdense basins, with Result::mps set to mps.
   */
  Result findVoidsBox (const std::vector<double>& tracers,
                       const std::vector<double>& displacement,
                       const std::vector<std::uint8_t>& valid,
                       double mps,
                       const Config& config,
                       const std::optional<otswap::BoxLimits>& limits = std::nullopt,
                       const std::vector<std::int64_t>& ids = {},
                       const DivergenceGrid* grid = nullptr);

  /**
   *  @brief Find voids and overdense basins in a box, from a
   *  reconstruction.
   *
   *  As the overload above, with tracers, displacement, valid, mps and ids
   *  taken from the reconstruction. The grid is built on
   *  otswap::Result::box, the box the reconstruction used: its limits
   *  (otswap::SelectionCounts::boxLimits), or without them the bounding box
   *  of the tracers and the given randoms (of the tracers alone when the
   *  randoms were drawn). Every endpoint otswap computes lies inside it.
   *  The tracers flagged in otswap::Result::outsideBox are skipped without
   *  being read. Given the same arrays and otswap::Result::box as the
   *  limits, the result equals that of the overload above.
   *
   *  @param reconstruction a box result of otswap::reconstructBox.
   *
   *  @param config the parameters of the search.
   *
   *  @param grid a grid to reuse, or null, as above. Its mps wins over the
   *  reconstruction's without an error; at Detailed verbosity the report
   *  states both.
   *
   *  @exception Error if the reconstruction is not a box result; if it is
   *  malformed: nObjects or nRealizations 0, tracers, displacement, valid,
   *  outsideBox or ids of a size that does not match them, valid holding a
   *  value other than 0 and 1, ids holding a value twice, box missing, with
   *  a bound not finite or a lower bound above the upper one, or different
   *  from selection.boxLimits when those are set, a tracer not flagged in
   *  outsideBox lying outside box; if its mps is not finite and positive
   *  and no grid is given; and as above, but for
   *  the mps of the grid, which wins over the reconstruction's. A positive
   *  smoothing.sigma in Mps is checked against the mps of the search: the
   *  grid's when given, otherwise the reconstruction's.
   */
  Result findVoidsBox (const otswap::Result& reconstruction,
                       const Config& config,
                       const DivergenceGrid* grid = nullptr);

}

#endif
