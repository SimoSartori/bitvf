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

    /// @return the source file in which the Error was constructed, as the
    /// compiler named it.
    const char* file () const noexcept { return m_file; }

    /// @return the name of the function in which the Error was
    /// constructed.
    const char* function () const noexcept { return m_function; }

    /// @return the line at which the Error was constructed.
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
   *      [cell]       DivergenceGrid::divergence, divergenceSmoothed,
   *                   volume, weight, weightSmoothed, cellRaMin ..
   *                   cellDistanceMax
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
   *  at most 1e-12 times the largest tracer weight of the grid takes 0:
   *  so does a cell filled from its neighbours (weight 0) when the cells
   *  around hold it alone. A larger sigma merges small basins into larger
   *  ones; sigma = 0 disables the smoothing. With sigma > 0 the smoothed field and the smoothed weight
   *  are kept beside the divergence and the weight
   *  (DivergenceGrid::divergenceSmoothed, weightSmoothed).
   *
   *  Box: the cells around are those within floor(3.5 sigma / s) cells
   *  along each axis, s the cell side, cut at the faces of the grid: a
   *  cube, not a sphere. When 3.5 sigma / s lies within 1e-9 relative of
   *  an integer, that integer is taken instead, so that a cell side
   *  rounded by one ulp gives the same cube. Below sigma = s / 3.5 the
   *  cube holds the cell alone. In Mps, sigma is multiplied by the box
   *  mps.
   *
   *  Lightcone: the cells around are those whose centre lies closer than
   *  4 sigma, with the same weights. Masked cells (DivergenceGrid::masked) enter no mean and take 0. In
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
    /// Empty: min(1, (20 / nRealizations)^(1/3)), that is 1 up to 20
    /// realizations and (20 / nRealizations)^(1/3) above, the same rule in
    /// both geometries. A search on a mean displacement file runs with one
    /// realization, so its default is 1. Must be positive when set, and
    /// left empty when a grid is reused. Result::config holds the value
    /// used. Unrelated to otswap::Config::cellSize, the cell of the
    /// neighbour search of the reconstruction, which affects speed only:
    /// this one changes the result.
    std::optional<double> divergenceCellSize;

    /// Fill Result::divergence with the grid and the divergence field.
    /// Off by default: it holds one entry per cell in several arrays.
    bool keepDivergence = false;

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
    /// "bitvf: " and ends in '\n'.
    ///
    /// - Silent: nothing.
    /// - Normal: the selection, then one line per call with its time:
    ///
    ///       bitvf: kept <k> of <N> tracers: <b> outside the box [<x0>, <x1>] x [<y0>, <y1>] x [<z0>, <z1>]
    ///       bitvf: kept <k'> of <k> tracers: <g> outside the reused grid
    ///       bitvf: findVoidsBox: <V> voids and <O> overdense basins from <N> tracers in <T> s; <W> without a valid displacement
    ///       bitvf: findVoidsBoxFromDivergence: <V> voids and <O> overdense basins from <N> tracers in <T> s
    ///       bitvf: findVoidsBoxFromField: <V> voids and <O> overdense basins on <nx> x <ny> x <nz> cells in <T> s
    ///
    ///   The box line is written when the search has box limits: those
    ///   given, the reconstruction's (otswap::SelectionCounts::boxLimits)
    ///   or the grid's (DivergenceGrid::limits); the grid line whenever a
    ///   grid is reused, even with no tracer outside it, and by
    ///   findVoidsBoxFromDivergence, as "outside the divergence grid". \<W\>
    ///   counts the tracers read without a valid displacement. A count of 1
    ///   takes the singular.
    /// - Detailed: as Normal, with one line per step of the search in the
    ///   order it runs, each with its time, and the values that describe
    ///   it: the mean particle separation and its source, the grid, its
    ///   cell side and the box it covers, the smoothing sigma, the basins
    ///   found by the watershed.
    ///
    /// Each line is written as soon as its step finishes, by the calling
    /// thread, never inside a parallel region. A call that throws leaves
    /// the lines of the steps it finished, and never writes its call line.
    /// Every count printed can be read from the object returned
    /// (Result::counts, the catalogues' nBasins); step times and
    /// descriptive values (mps, cell side, box, sigma) appear only in the
    /// lines. Times have two decimals.
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
    /// shifted by s sum_n div_n (x_n - x0) / sum_n |div_n| over the seed's
    /// neighbouring cells n, with div the field the watershed split and
    /// s = -1 for a void, +1 for an overdense basin (towards lower divergence
    /// in a void, higher in an overdense basin), and kept inside the seed
    /// cell. Flat, [basin][xyz].
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
    /// the tracers to the randoms over the whole search (Counts::randoms);
    /// NaN when no random lies in the sphere. The expected count
    /// n_r N_t / N_r has a relative noise of about 1 / sqrt(33.5 q f), with
    /// q the randoms per tracer and f the fraction of the sphere inside the
    /// survey: at least 10 to 20 randoms per tracer is a reasonable start.
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

    /// findVoidsBoxFromField only, empty otherwise: the field the
    /// watershed split, smoothed when a smoothing was applied, at centre,
    /// interpolated trilinearly between the eight nearest cell centres.
    /// The point is first clamped, per axis, to the range of the cell
    /// centres; an axis with one cell takes that cell's value. At a cell
    /// centre it is that cell's value.
    std::vector<double> fieldCentre;
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
   *  weights, counts, basin indices, cell limits), its box limits and its
   *  realization counts are ignored, whatever their size.
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
    /// cellSizeMps times mps. In a result of findVoidsBoxFromField, the
    /// spacing (max - min) / (n - 1) of the centres along the axis with the
    /// most nodes, x before y before z on a tie, with cellSizeMps that
    /// spacing divided by mps. NaN in a lightcone.
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
    /// three axes but in a result of findVoidsBoxFromField, where they may
    /// differ. 0 in a lightcone.
    std::array<std::size_t, 3> nCellsAxis = {0, 0, 0};

    /// Box: the box limits of the search that filled the grid; the tracers
    /// outside them were skipped. Empty without limits and in a lightcone.
    /// Not part of the cells: a search reusing the grid does not read it,
    /// its own limits apply.
    std::optional<otswap::BoxLimits> limits;

    /// Box: the realizations of the search that filled the grid (keyword
    /// NREC). A record: a search reusing the grid does not read it. 0 in a
    /// lightcone and in a result of findVoidsBoxFromField.
    unsigned nRealizations = 0;

    /// Box: Config::nRealizationsMean of the search that filled the grid
    /// (keyword NRECMEAN). A record, as nRealizations. 0 in a lightcone and
    /// in a result of findVoidsBoxFromField.
    unsigned nRealizationsMean = 0;

    /// Lightcone, defines the cells: comoving distance range of the
    /// tracers and randoms the grid covers, in coordinate units. NaN in a
    /// box.
    double distanceMin = std::numeric_limits<double>::quiet_NaN();
    double distanceMax = std::numeric_limits<double>::quiet_NaN();

    /// Lightcone, defines the cells: right-ascension range of the
    /// footprint, radians, cut at its largest gap on the circle: raMin, the
    /// first value after the gap, lies in [0, 2 pi), and raMax may exceed
    /// 2 pi when the footprint crosses RA = 0. NaN in a box.
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
    /// a search that reuses the grid follow it. NaN in a lightcone, and in
    /// a result of findVoidsBoxFromField given no mps (cellSizeMps too).
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

    /// The divergence of each cell, before smoothing; what the watershed
    /// splits with sigma = 0.
    std::vector<double> divergence;

    /// The divergence after smoothing (see Smoothing), what the watershed
    /// splits, with sigma > 0; empty with sigma = 0.
    std::vector<double> divergenceSmoothed;

    /// Volume of each cell, in coordinate units cubed: uniform in a box.
    std::vector<double> volume;

    /// Tracer weight of each cell, the weight of the smoothing:
    /// nDisplacements / nRealizations where a displacement contributed.
    /// A cell that no displacement contributed to (nDisplacements 0) is
    /// empty: it keeps a weight of 0, and the fill gives it the mean
    /// divergence of nearby cells that hold a value, a contribution or a
    /// value filled in an earlier pass, weighted by 1 / r^2, r the
    /// distance between the cell centres. The fill runs in passes, each
    /// reading the values left by the previous one. With its weight of 0,
    /// a filled value enters no smoothed mean; with sigma = 0 the
    /// watershed reads it as it is.
    ///
    /// Box: the neighbours within one cell along each axis, or within two
    /// when fewer than three of those hold a value; the passes repeat
    /// until no cell is empty (an Error when a pass fills none).
    ///
    /// Lightcone: the adjacent cells that are not masked, when there are
    /// at least three; the passes repeat until one fills no cell. A cell
    /// left empty keeps a divergence of 0, and a masked cell is never
    /// filled.
    std::vector<double> weight;

    /// Smoothed tracer weight of each cell, with sigma > 0: the mean of
    /// weight over the cells of the smoothing (see Smoothing) with the
    /// Gaussian alone as weight, sum K w / sum K, K = exp(-r^2 /
    /// (2 sigma^2)). It shows how much of a cell's smoothed divergence
    /// comes from displacements; the watershed does not read it. Box: the
    /// sum of K over the cube is the product of the three sums of the 1-D
    /// kernel along the axes, each cut at the faces. Empty with sigma = 0.
    std::vector<double> weightSmoothed;

    /// Displacement ends counted in each cell, before the fill: a valid
    /// displacement that leaves its tracer's cell counts 1 in that cell and
    /// 1 in the cell it ends in; one that ends in its tracer's cell counts 2
    /// there. In a box, a point on a face between two cells belongs to the
    /// upper one, and a displacement that starts or ends on a face, an edge
    /// or a vertex, or runs along a face, is counted as that of the same
    /// tracer moved one ulp into its cell. In a lightcone, displacements
    /// rejected count nowhere (Counts::displacementsRejected).
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

    /// Row of the tracer's void in Result::voids; -1 none. In a
    /// lightcone, a tracer that starts in a masked cell
    /// (DivergenceGrid::masked) keeps that cell in cellIndex but is in no
    /// basin (Counts::tracersInMaskedCells).
    std::vector<std::int64_t> voidIndex;

    /// Row of the tracer's overdense basin in Result::overdense; -1 none.
    std::vector<std::int64_t> overdenseIndex;

    /// Distance from the centre of the tracer's basin, its void or its
    /// overdense basin, in coordinate units; NaN for a tracer not in a
    /// basin (not read by the search).
    std::vector<double> distance;

    /// Lightcone: 0 when the tracer falls on an unobserved pixel of the
    /// mask. Empty in a box.
    std::vector<std::uint8_t> inMask;

    /// Lightcone: 0 when the tracer lies outside the redshift cut of the
    /// search, and so was skipped; all 1 without a cut. Empty in a box.
    std::vector<std::uint8_t> inCut;

    /// Box: 0 when the tracer lies outside the box limits, and so was
    /// skipped; all 1 without limits.
    std::vector<std::uint8_t> inBox;

    /// 1 when at least one valid displacement of the tracer entered the
    /// divergence. Empty after findVoidsBoxFromDivergence, which reads no
    /// displacement.
    std::vector<std::uint8_t> used;
  };

  /// What the void finder counted: the counts its report prints, with the
  /// catalogues' nBasins (Config::verbosity).
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

    /// Lightcone: N_r of delta2Mps, the randoms the search reads: kept by
    /// the mask and the redshift cut, and inside the grid. 0 in a box.
    std::size_t randoms = 0;

    /// Lightcone: tracers on unobserved pixels of the mask, and so skipped
    /// (TracerAssignment::inMask), whether or not they lie inside the
    /// redshift cut. 0 in a box.
    std::size_t tracersOutsideMask = 0;

    /// Lightcone: of the tracers on observed pixels, those outside the
    /// redshift cut, and so skipped (TracerAssignment::inCut is 0). A
    /// tracer both masked and outside the cut counts only in
    /// tracersOutsideMask, where its flag TracerAssignment::inMask is 0:
    /// the counts give the mask priority, while its inCut is 0 as well.
    /// 0 in a box.
    std::size_t tracersOutsideRedshiftCut = 0;

    /// Lightcone: randoms on unobserved pixels of the mask, removed, whether
    /// or not they lie inside the redshift cut. 0 in a box.
    std::size_t randomsOutsideMask = 0;

    /// Lightcone: of the randoms on observed pixels, those outside the
    /// redshift cut, removed. A random both masked and outside the cut counts
    /// only in randomsOutsideMask. 0 in a box.
    std::size_t randomsOutsideRedshiftCut = 0;

    /// Lightcone: of the tracers read, those on an observed pixel whose
    /// cell is masked (its centre falls on an unobserved pixel,
    /// DivergenceGrid::masked). They keep their cell in
    /// TracerAssignment::cellIndex but are in no basin, and their
    /// displacements are rejected. A smaller cell size makes them fewer.
    /// At Normal verbosity, when not 0, the search reports them on a line
    /// of their own:
    ///
    ///     bitvf: <n> tracers start in a masked cell, in no basin
    ///
    /// 0 in a box.
    std::size_t tracersInMaskedCells = 0;

    std::size_t randomsOutsideGrid = 0;  ///< lightcone, reused grid: randoms outside it, removed
    std::size_t displacements = 0;       ///< valid displacements given

    /// Lightcone: valid displacements of the tracers read that were
    /// rejected (masked cells, outside or below the grid): the start or the
    /// end lies in a masked cell, the end lies outside every cell, or,
    /// without a central sphere, the path passes closer to the observer
    /// than the grid's first edge. A path through unobserved pixels is the
    /// reconstruction's to reject, with its own mask and threshold. 0 in a
    /// box.
    std::size_t displacementsRejected = 0;
  };

  /**
   *  @brief Output of one void search.
   *
   *  Given the same inputs, every field but elapsedSeconds is the same bit
   *  for bit whatever the number of threads.
   *
   *  A result of findVoidsBoxFromField, a field result, holds no tracer:
   *  nTracers and nRealizations are 0, which tells it from the result of
   *  every other search, and so are the counts; ids, the per-tracer flags,
   *  assignment, and the catalogue columns deltaAv, delta2Mps and nTracers
   *  are empty; limits is empty; mps is NaN when none was given.
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
   *  (BasinCatalog). Reports as Config::verbosity describes.
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
   *    nRealizationsMean is 0; if smoothing.sigma is positive but, in
   *    coordinate units (times mps in Mps), below about 1e-162 or above
   *    about 1e154, too small or too large to build the kernel;
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
   *  limits, the result equals that of the overload above. Reports as
   *  Config::verbosity describes.
   *
   *  @param reconstruction a box result of otswap::reconstructBox, or the
   *  mean of one: the Result read from a mean displacement file by
   *  otswap::io::readDisplacementMean (nRealizations 1,
   *  nRealizationsMean the number averaged; set Config::nRealizationsMean
   *  to it, so that bitvf's files record it). The box is the file's box
   *  keywords, with or without BOXLIM; a file without them (written by
   *  otswap 0.2) holds no box limits and is refused. The divergence is not
   *  linear in the displacement, so a search on the mean of R realizations
   *  does not give the result of the search on the R realizations.
   *
   *  Without box limits, a tracer rounded just outside the box, as in an
   *  ASCII mean file, is moved onto it: on an axis where it lies outside by
   *  at most 1e-8 times the larger of its coordinate and the bound, it is
   *  put on the lower bound, or on the largest value below the upper one.
   *  The displacement is kept. The search and the files written from the
   *  result use the moved positions; the report states how many tracers
   *  were moved. With box limits, the tracers outside are excluded.
   *
   *  @param config the parameters of the search.
   *
   *  @param grid a grid to reuse, or null, as above. Its mps wins over the
   *  reconstruction's without an error; at Detailed verbosity the report
   *  states both.
   *
   *  @exception Error if the reconstruction is not a box result; if it is
   *  malformed: nObjects, nRealizations or nRealizationsMean 0,
   *  nRealizationsMean above 1 with nRealizations not 1, validRealizations
   *  neither empty nor of nObjects entries or holding more than
   *  nRealizations * nRealizationsMean, tracers, displacement, valid,
   *  outsideBox or ids of a size that does not match them, valid holding a
   *  value other than 0 and 1, ids holding a value twice, box missing, with
   *  a bound not finite or a lower bound above the upper one, or different
   *  from selection.boxLimits when those are set; selection.boxLimits with
   *  a bound not finite, a width not finite, or a lower bound not below
   *  the upper one; a tracer not flagged in
   *  outsideBox lying outside box by more than the tolerance above, or at all when box
   *  limits are set; if its mps is not finite and positive
   *  and no grid is given; and as above, but for
   *  the mps of the grid, which wins over the reconstruction's. A positive
   *  smoothing.sigma in Mps is checked against the mps of the search: the
   *  grid's when given, otherwise the reconstruction's.
   *
   *  @return the voids and overdense basins, with Result::mps the mps of
   *  the search: the grid's when given, otherwise the reconstruction's.
   */
  Result findVoidsBox (const otswap::Result& reconstruction,
                       const Config& config,
                       const DivergenceGrid* grid = nullptr);

  /**
   *  @brief Find voids and overdense basins in a box from a divergence
   *  computed earlier, with the tracers given.
   *
   *  The field before smoothing (DivergenceGrid::divergence) is
   *  smoothed with config.smoothing, in units of the grid's mps, split into
   *  basins and catalogued as findVoidsBox does; a different sigma is the
   *  usual reason for the call. With the tracers of the search that
   *  produced the grid and the same smoothing, the voids, the overdense
   *  basins, the smoothed divergence, the smoothed weight and the
   *  assignment but its used flags equal that search's bit for bit,
   *  sigma = 0 included. A grid kept in memory
   *  (Config::keepDivergence) serves as well as one read from a divergence
   *  file (io::readDivergence). Reports as Config::verbosity describes.
   *
   *  The tracers are read as with a reused grid: those outside the grid's
   *  box limits (DivergenceGrid::limits) are skipped
   *  (TracerAssignment::inBox, Counts::tracersOutsideBox), those inside
   *  but outside every cell are excluded (Result::outsideGrid,
   *  Counts::tracersOutsideGrid). The tracers read in each cell are
   *  counted; when the grid holds nTracers, every cell's count must equal
   *  it, or the tracers are not those of the field.
   *
   *  No displacement enters this call: TracerAssignment::used is empty,
   *  Counts::displacements and Counts::tracersUsed are 0. The grid's
   *  nDisplacements, nRealizations and nRealizationsMean are copied into
   *  Result::divergence, Result::nRealizations and
   *  Result::config.nRealizationsMean, so the files written from the result
   *  record those of the field. Result::limits is the grid's limits,
   *  Result::mps its mps, Result::config.divergenceCellSize its
   *  cellSizeMps.
   *
   *  @param divergence a box grid holding the fields that define its cells,
   *  divergence and weight, and optionally nTracers, nDisplacements,
   *  limits, nRealizations and nRealizationsMean; its other per-cell
   *  fields are not read: divergenceSmoothed and weightSmoothed of the
   *  result are computed from these with config.smoothing.
   *
   *  @param tracers Cartesian coordinates of the tracers, 3 * nTracers
   *  entries; at least one.
   *
   *  @param config the parameters of the search; divergenceCellSize must
   *  be empty.
   *
   *  @param ids the identity of each tracer, copied into Result::ids: one
   *  per tracer, no value twice; or empty.
   *
   *  @exception Error
   *  - if tracers is empty, its size is not a multiple of three, or a
   *    coordinate is not finite; if ids is neither empty nor of size
   *    nTracers, or holds a value twice;
   *  - on config as findVoidsBox; if divergenceCellSize is set;
   *  - on the grid as on a grid to reuse in findVoidsBox (a lightcone grid,
   *    or nCellsAxis not equal on the three axes, among others); if
   *    divergence or weight does not hold nCells entries, or an
   *    entry of either is not finite; if nTracers is neither empty nor of
   *    nCells entries; if a bound of its limits is not finite, or a lower
   *    bound is not below the upper one;
   *  - if smoothing.sigma is positive but too small or too large in
   *    coordinate units to build the kernel, as findVoidsBox;
   *  - on the tracers: if the search reads more than 4294967295 tracers;
   *    if none is left to read; if the count of a cell differs from the
   *    grid's nTracers;
   *  - if the smoothed divergence is 0 or positive in every cell, so that
   *    no void can be found.
   *
   *  @return the voids and overdense basins of the divergence smoothed
   *  with config.smoothing, with Result::mps the grid's mps.
   */
  Result findVoidsBoxFromDivergence (const DivergenceGrid& divergence,
                                     const std::vector<double>& tracers,
                                     const Config& config,
                                     const std::vector<std::int64_t>& ids = {});

  /**
   *  @brief Find voids and overdense basins in a field given at the
   *  centres of the cells of a box lattice, without tracers.
   *
   *  The lattice is inferred from the centres. On each axis, the
   *  coordinates are sorted and split into nodes at the gaps larger than
   *  half the largest gap; with n nodes from min to max, the spacing is
   *  (max - min) / (n - 1). The spacings of the axes with more than one
   *  node must agree within 1e-6 relative: the cell side s is the spacing
   *  of the axis with the most nodes, x before y before z on a tie, and is
   *  DivergenceGrid::cellSize of the result; the other axes keep their
   *  own spacing only to place the centres on nodes. Each centre
   *  must lie within 1e-6 s of its node on every axis, and each node of
   *  the nx * ny * nz must hold exactly one centre. The grid spans min -
   *  s / 2 to max + s / 2 on each axis, without padding, with cubic cells
   *  of side s; nx, ny, nz need not be equal. The cells of the result are
   *  the rows of the lattice, cell (i, j, k) row (i * ny + j) * nz + k,
   *  not the rows of the input, and no map back is returned: the order of
   *  the input does not change the result. Reports as Config::verbosity
   *  describes.
   *
   *  The field is smoothed with config.smoothing, as the divergence of
   *  findVoidsBox, with weights as the tracer weight; split into basins at
   *  0 as BasinCatalog describes: seeds below 0 start voids, those at or
   *  above 0 overdense basins. A field with another convention, such as a
   *  density or a density contrast with voids above 0, must be shifted or
   *  flipped before the call. The catalogues hold centre, radius,
   *  radiusExt, edgeFlag and contamination as findVoidsBox computes them,
   *  on the faces of the inferred grid, and fieldCentre; deltaAv,
   *  delta2Mps and nTracers, which need tracers, are empty. Result is a
   *  field result (see Result). With Config::keepDivergence,
   *  Result::divergence holds the grid, the field as given (divergence)
   *  and, with sigma > 0, smoothed (divergenceSmoothed), the weights as
   *  given and, with sigma > 0, their smoothed mean (weightSmoothed), the
   *  basin indices, the cell centres and volumes; nDisplacements and
   *  nTracers are empty.
   *
   *  A bitvf divergence file is valid input, read with io::readField: its
   *  cell centres, divergence and weight with the smoothing of
   *  the search that wrote it give back that search's basins. Their
   *  centre, radius and radiusExt can differ from that search's by a few
   *  ulp: they are computed from the grid's lower and s, which are
   *  inferred here from the centres read and may differ from the original
   *  grid's by rounding.
   *
   *  @param centres Cartesian coordinates of the cell centres, 3 * nCells
   *  entries, [cell][xyz], in any order.
   *
   *  @param values the field at each centre, nCells entries, in the order
   *  of centres.
   *
   *  @param config the parameters of the search. divergenceCellSize must
   *  be empty, keepAssignment unset and nRealizationsMean 1: none applies
   *  without tracers. Without mps, a positive smoothing.sigma must be in
   *  LengthUnit::Coordinate.
   *
   *  @param weights the weight of each cell in the smoothing, as
   *  DivergenceGrid::weight in findVoidsBox, nCells entries in the order
   *  of centres; or empty, which gives every cell the weight 1. Used as
   *  given: a cell whose weights sum to at most 1e-12 times the largest
   *  weight takes 0 (see Smoothing).
   *
   *  @param mps a mean particle separation, in coordinate units, or none.
   *  Given, it is the unit of a sigma in LengthUnit::Mps and is recorded
   *  in Result::mps, with the cell size in units of it in
   *  Result::config.divergenceCellSize.
   *
   *  @exception Error
   *  - on the arrays: if centres is empty, its size is not a multiple of
   *    three, or a coordinate is not finite; if values does not hold one
   *    entry per centre, or an entry is not finite; if weights is neither
   *    empty nor of one entry per centre, or an entry is negative or not
   *    finite;
   *  - if mps is given and not finite and positive;
   *  - on config as findVoidsBox; if divergenceCellSize is set,
   *    keepAssignment is set or nRealizationsMean is not 1; if
   *    smoothing.sigma is positive in LengthUnit::Mps without mps; if
   *    smoothing.sigma is positive but too small or too large in
   *    coordinate units to build the kernel, as findVoidsBox; if
   *    smoothing.sigma is positive and every weight is 0;
   *  - on the lattice: if every axis has one node (a single cell has no
   *    side); if the spacings of two axes differ by more than 1e-6
   *    relative; if a centre lies farther than 1e-6 s from its node; if a
   *    node holds two centres, or none;
   *  - if the field, smoothed when a smoothing is set, is 0 or positive in
   *    every cell, so that no void can be found.
   *
   *  @return the voids and overdense basins of the field, with
   *  Result::nTracers 0, and Result::mps the mps given, or NaN without
   *  one.
   */
  Result findVoidsBoxFromField (const std::vector<double>& centres,
                                const std::vector<double>& values,
                                const Config& config,
                                const std::vector<double>& weights = {},
                                std::optional<double> mps = std::nullopt);

  // ==========================================================================
  // Grid and watershed
  // ==========================================================================

  /**
   *  @brief The grid findVoidsBox builds on box limits, with this mps and
   *  cell size, bit for bit.
   *
   *  Only the fields that define the cells are set: geometry, lower, upper,
   *  nCellsAxis, nCells, cellSize, cellSizeMps and mps. limits,
   *  nRealizations, nRealizationsMean and every per-cell field are left
   *  empty or 0: records of a search, which a reuse ignores. The grid
   *  passes the checks of a grid to reuse, so that several searches can
   *  share it.
   *
   *  @param box the box limits the grid is built on, as the limits of
   *  findVoidsBox; each bound finite, with lower[k] < upper[k].
   *
   *  @param mps mean particle separation, in coordinate units; finite and
   *  positive.
   *
   *  @param divergenceCellSize the cell side in units of mps, as
   *  Config::divergenceCellSize; finite and positive.
   *
   *  @exception Error if a bound of box is not finite or a lower bound is
   *  not below the upper one; if mps or divergenceCellSize is not finite
   *  and positive; if the grid would hold more cells than an int64 index
   *  can address.
   *
   *  @return the grid, with only the fields that define the cells set.
   */
  DivergenceGrid boxGrid (const otswap::BoxLimits& box, double mps, double divergenceCellSize);

  /**
   *  @brief The grid findVoidsBox builds from this reconstruction, bit for
   *  bit: on otswap::Result::box, with otswap::Result::mps.
   *
   *  As the overload above. Without divergenceCellSize, the cell side is
   *  the default for the reconstruction's realizations,
   *  min(1, (20 / nRealizations)^(1/3)) mps.
   *
   *  @param reconstruction a box result of otswap::reconstructBox, or the
   *  mean of one: the Result read from a mean displacement file by
   *  otswap::io::readDisplacementMean (nRealizations 1,
   *  nRealizationsMean the number averaged; set Config::nRealizationsMean
   *  to it, so that bitvf's files record it). The box is the file's box
   *  keywords, with or without BOXLIM; a file without them (written by
   *  otswap 0.2) holds no box limits and is refused. The divergence is not
   *  linear in the displacement, so a search on the mean of R realizations
   *  does not give the result of the search on the R realizations.
   *
   *  @param divergenceCellSize the cell side in units of mps, as
   *  Config::divergenceCellSize; finite and positive. Empty for the
   *  default.
   *
   *  @exception Error if the reconstruction is not a box result; if its
   *  nRealizations is 0; if its box is missing, has a bound not finite or
   *  a lower bound above the upper one, or differs from
   *  selection.boxLimits when those are set; if its mps is not finite and
   *  positive; if divergenceCellSize is set and not finite and positive;
   *  if the grid would hold more cells than an int64 index can address.
   *
   *  @return the grid, with only the fields that define the cells set.
   */
  DivergenceGrid boxGrid (const otswap::Result& reconstruction,
                          std::optional<double> divergenceCellSize = std::nullopt);

  /// The basins of a field on a box grid, as watershedBox finds them,
  /// with no catalogue.
  struct Partition {

    std::size_t nVoids = 0;
    std::size_t nOverdense = 0;

    /// [cell]: row of the cell's void, -1 none.
    std::vector<std::int64_t> voidIndex;

    /// [cell]: row of the cell's overdense basin, -1 none.
    std::vector<std::int64_t> overdenseIndex;

    /// [void]: row of its seed cell, ascending.
    std::vector<std::int64_t> voidSeed;

    /// [overdense basin]: row of its seed cell, ascending.
    std::vector<std::int64_t> overdenseSeed;

    /// [void]: cells of the extended void, and of its core (field below 0).
    std::vector<std::uint64_t> voidCells, voidCoreCells;

    /// [overdense basin]: cells of the basin, and of its core (field at or
    /// above 0).
    std::vector<std::uint64_t> overdenseCells, overdenseCoreCells;
  };

  /**
   *  @brief Split a field on a box grid into basins by steepest descent.
   *
   *  The descent is that of the void finder (BasinCatalog): each cell
   *  flows to the neighbour, among the 26 cut at the faces of the grid,
   *  where the field drops the most, ties to the lowest row; a cell with
   *  no strictly smaller neighbour is a seed. Seeds below 0 start voids,
   *  seeds at or above 0 overdense basins; each kind is numbered by
   *  ascending seed row. A field with another convention, such as a
   *  density or a density contrast with voids above 0, must be shifted or
   *  flipped before the call. No smoothing, no fill and no centres: the
   *  field is used as given. A field with no cell below 0 gives no void,
   *  and no error. The result does not depend on the number of threads.
   *
   *  @param field the value of each cell, nx * ny * nz entries; cell
   *  (i, j, k) is row (i * ny + j) * nz + k.
   *
   *  @param nCellsAxis the cells along x, y, z, (nx, ny, nz), each at
   *  least 1; they need not be equal.
   *
   *  @param verbosity,log what the call reports and where, as
   *  Config::verbosity and Config::log: at Normal, one line with its time,
   *
   *      bitvf: watershedBox: <V> voids and <O> overdense basins on <nx> x <ny> x <nz> cells in <T> s
   *
   *  and at Detailed, before it, the cells in voids and in overdense
   *  basins, and how many of each lie in the cores.
   *
   *  @exception Error if a count of nCellsAxis is 0; if field does not
   *  hold nx * ny * nz entries; if an entry of field is not finite.
   *
   *  @return the voids and overdense basins: the basin of each cell, the
   *  seed and the cells of each basin.
   */
  Partition watershedBox (const std::vector<double>& field, const std::array<std::size_t, 3>& nCellsAxis,
                          otswap::Verbosity verbosity = otswap::Verbosity::Normal,
                          std::ostream* log = &std::clog);

}

#endif
