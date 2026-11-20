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
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "bitvf/VoidFinder.h"

namespace bitvf {
  namespace internal {

    /**
     *  @brief The bitvf version, as "major.minor.patch": the version the
     *  writers record in the files they write.
     */
    extern const char* const version;

    /**
     *  @brief Call otswap, rethrowing an otswap::Error as an Error
     *  prefixed "otswap: ", the original nested.
     */
    template <typename F>
    auto call_otswap (F&& f)
    {
      try {
        return f();
      }
      catch (const otswap::Error& e) {
        std::throw_with_nested(Error(std::string("otswap: ") + e.what()));
      }
    }

    //==========================================================================

    /**
     *  @brief The report of one public call, as Config::verbosity
     *  describes, with the clock of its steps.
     *
     *  Built at the start of the call. Lines are written only from the
     *  calling thread, outside parallel regions, each as one insertion
     *  followed by a flush. The static members compose the parts of the
     *  text shared by several lines, in the stream returned by line().
     */
    class Report {
    public:

      /// @param verbosity,log the level and the stream of the call; a null
      /// log writes nothing, as Silent
      Report (otswap::Verbosity verbosity, std::ostream* log);

      /// Whether the Normal lines are written: Normal or Detailed, with a
      /// log.
      bool normal () const { return m_log != nullptr && m_verbosity != otswap::Verbosity::Silent; }

      /// Whether the Detailed lines are written.
      bool detailed () const { return m_log != nullptr && m_verbosity == otswap::Verbosity::Detailed; }

      /// Write "bitvf: ", text and '\n', then flush; nothing when not
      /// normal().
      void write (const std::string& text) const;

      /// Seconds since the previous lap, or since the start for the first.
      double lap ();

      /// Seconds since the start.
      double elapsed () const;

      /// Append " in <t> s" to out, with t the lap() just taken, and write
      /// the line: the line of a step that has finished.
      void step (std::ostringstream& out);

      /// An empty stream for one line, in the classic locale, at the
      /// default precision of 6 significant digits.
      static std::ostringstream line ();

      /// "in <t> s", t fixed with 2 decimals; out's format is restored.
      static void time (std::ostream& out, double t);

      /// "[x0, x1] x [y0, y1] x [z0, z1]", as otswap writes a box.
      static void box (std::ostream& out, const std::array<double, 3>& lower, const std::array<double, 3>& upper);

      /// "<nx> x <ny> x <nz> cells".
      static void shape (std::ostream& out, const std::array<std::size_t, 3>& n);

      /// "<V> voids and <O> overdense basins", singular for 1.
      static void basins (std::ostream& out, std::size_t nVoids, std::size_t nOverdense);

    private:

      otswap::Verbosity m_verbosity;
      std::ostream* m_log;
      std::chrono::steady_clock::time_point m_start;
      std::chrono::steady_clock::time_point m_last;
    };

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

      /// [tracer]: 0 for a tracer outside the box limits, or flagged in
      /// otswap::Result::outsideBox; all 1 otherwise.
      std::vector<std::uint8_t> inBox;

      /// [tracer]: 1 for a tracer inside the limits but outside every cell
      /// of the reused grid; empty without a reused grid.
      std::vector<std::uint8_t> outsideGrid;

      /// The mean particle separation of the search, in coordinate units.
      double mps = 0.;

      /// The smoothing sigma in coordinate units: Smoothing::sigma, times
      /// mps in Mps; 0 for no smoothing.
      double sigmaLength = 0.;

      /// The largest absolute value of a component of the valid
      /// displacements read, in coordinate units; it sets the scale of the
      /// fixed-point accumulation of the divergence.
      double maxDisplacement = 0.;

      /// The box the grid is built on: the limits given, or the
      /// reconstruction's box (otswap::Result::box). Absent, the grid
      /// covers the tracers read and the endpoints of their valid
      /// displacements.
      std::optional<otswap::BoxLimits> box;

      /// The box limits of the search: the limits given, or
      /// otswap::SelectionCounts::boxLimits; empty without them.
      std::optional<otswap::BoxLimits> limits;

      /// The number of tracers moved onto the box edge; 0 outside the
      /// otswap::Result overload.
      std::size_t tracersMoved = 0;

      /// The counts known after the checks: tracers, tracersOutsideBox,
      /// tracersOutsideGrid, displacements.
      Counts counts;
    };

    /**
     *  @brief Check the arguments of the array overload of findVoidsBox
     *  and select what the search reads.
     *
     *  The checks run in a fixed order, so that the same defects always
     *  raise the same error: the arrays, mps, the limits, ids, the
     *  configuration, the reused grid, the smoothing sigma in coordinate
     *  units, then the data.
     *
     *  @exception Error on every case listed by findVoidsBox, except those
     *  found by building the grid and by the divergence.
     */
    BoxSearch check_box (const std::vector<double>& tracers,
                         const std::vector<double>& displacement,
                         const std::vector<std::uint8_t>& valid,
                         double mps,
                         const Config& config,
                         const std::optional<otswap::BoxLimits>& limits,
                         const std::vector<std::int64_t>& ids,
                         const DivergenceGrid* grid);

    /**
     *  @brief Check the arguments of the otswap::Result overload of
     *  findVoidsBox and select what the search reads.
     *
     *  The box is otswap::Result::box. The tracers flagged in
     *  otswap::Result::outsideBox are skipped unread. With a grid, its mps is the search's, and the
     *  reconstruction's is not read.
     *
     *  Without box limits, a tracer outside the box by at most 1e-8 times
     *  the larger of its coordinate and the bound, on an axis, is moved onto
     *  the box on that axis; the positions are then copied into moved, which
     *  BoxSearch::tracers refers to, and BoxSearch::tracersMoved counts the
     *  tracers moved. moved stays empty when none is moved.
     *
     *  @exception Error as the overload above, and on a reconstruction that
     *  is not a box result or is malformed.
     */
    BoxSearch check_box (const otswap::Result& reconstruction,
                         const Config& config,
                         const DivergenceGrid* grid,
                         std::vector<double>& moved);

    /**
     *  @brief Check the arguments of findVoidsBoxFromDivergence and select
     *  the tracers it reads.
     *
     *  The checks run in a fixed order: the tracers, ids, the
     *  configuration, the grid's defining fields and limits, its per-cell
     *  fields, the smoothing sigma in coordinate units, then the tracers
     *  read. The tracers outside grid.limits are skipped (inBox 0), those
     *  inside but outside every cell flagged in outsideGrid. displacement
     *  and valid are left null; nRealizations is the grid's.
     *
     *  @exception Error on every case listed by findVoidsBoxFromDivergence
     *  but those of the tracer counts and of the watershed
     */
    BoxSearch check_box_from_divergence (const DivergenceGrid& grid,
                                         const std::vector<double>& tracers,
                                         const Config& config,
                                         const std::vector<std::int64_t>& ids);

    /**
     *  @brief Check the arguments of findVoidsBoxFromField but the
     *  lattice.
     *
     *  The checks run in a fixed order: centres, values, weights, mps, the
     *  configuration, the unit of the smoothing, the smoothing sigma in
     *  coordinate units, then the weights against the smoothing.
     *
     *  @return the smoothing sigma in coordinate units, 0 for no smoothing
     *  @exception Error on every case listed by findVoidsBoxFromField but
     *  those of the lattice and of the watershed
     */
    double check_box_from_field (const std::vector<double>& centres,
                                 const std::vector<double>& values,
                                 const Config& config,
                                 const std::vector<double>& weights,
                                 const std::optional<double>& mps);

    /**
     *  @brief Check box limits: every bound finite, each lower bound below
     *  the upper one.
     *  @param name the name of the limits in the messages, plural
     *  @param fn the public function that starts the messages; empty for limits the user gave
     *  @exception Error on the first bound that fails
     */
    void check_limits (const otswap::BoxLimits& limits, const std::string& name, const std::string& fn = "");

    /**
     *  @brief Check a flat array of Cartesian coordinates, every entry.
     *  @param name the name of the array in the messages
     *  @param fn the public function that starts the messages about the size
     *  @param allowEmpty whether an empty array is accepted
     *  @return the number of rows
     *  @exception Error if the array is empty and may not be, if its size
     *  is not a multiple of three, or on the first entry not finite
     */
    std::size_t check_coordinates (const std::vector<double>& a, const std::string& name, bool allowEmpty,
                                   const std::string& fn);

    /**
     *  @brief Check the shape of the displacements of N tracers and of their
     *  validity flags.
     *  @param fn the public function that starts the messages
     *  @return the number of realizations
     *  @exception Error if displacement is empty or not a whole number of
     *  realizations of N tracers, if that number exceeds UINT_MAX, if valid
     *  is neither empty nor one per displacement, or holds a value other than 0 and 1
     */
    std::size_t check_displacement_shape (const std::vector<double>& displacement,
                                          const std::vector<std::uint8_t>& valid, std::size_t N,
                                          const std::string& fn);

    /**
     *  @brief Check the IDs of the tracers: one per tracer, or none, and no
     *  value twice.
     *  @param fn the public function that starts the messages
     *  @param fromResult whether the IDs are a reconstruction's (the message then names it)
     *  @exception Error on a wrong size, or on the first value held twice
     */
    void check_ids (const std::vector<std::int64_t>& ids, std::size_t nTracers, const std::string& fn, bool fromResult);

    /**
     *  @brief Check the box of a box reconstruction (otswap::Result::box), for the public function fn:
     *  present, every bound finite, no lower bound above the upper one,
     *  and equal to otswap::SelectionCounts::boxLimits when those are set.
     *  @exception Error, as a malformed reconstruction
     */
    void check_reconstruction_box (const otswap::Result& reconstruction, const std::string& fn);

    /**
     *  @brief Check what a box and a lightcone reconstruction share, for
     *  the public function fn: at least one object and one realization;
     *  nRealizationsMean at least 1, and above 1 only with one
     *  realization; validRealizations empty or one per object, each at
     *  most nRealizations * nRealizationsMean; tracers, displacement and
     *  valid of their layout's size, ids empty or one per object; every
     *  valid entry 0 or 1, no ID twice. The fields of one geometry are
     *  left to its caller.
     *  @exception Error, as a malformed reconstruction
     */
    void check_reconstruction_shape (const otswap::Result& reconstruction, const std::string& fn);

    /**
     *  @brief Check the configuration of a search.
     *  @param gridGiven whether a grid is reused, with which
     *  Config::divergenceCellSize must be empty
     *  @exception Error on an invalid field or combination
     */
    void check_config (const Config& config, bool gridGiven);

    /**
     *  @brief Check the fields of a box grid to reuse that define its cells,
     *  and cellSizeMps against them; also the record of a divergence file
     *  before it is written.
     *  @param mps the mps of the search
     *  @param mpsMustMatch whether mps must equal grid.mps, bit for bit
     *  @param fn the public function that starts the messages; empty when the caller adds its own
     *  @param name the name of the grid in the messages
     *  @exception Error on an inconsistent field, or on a different mps
     *  when it must match
     */
    void check_box_grid (const DivergenceGrid& grid, double mps, bool mpsMustMatch, const std::string& fn,
                         const std::string& name = "the grid to reuse");

    //==========================================================================

    /// The most cells a box grid may have along an axis, and a lightcone
    /// grid at each level of its cell index (shells, rings of a shell,
    /// sectors of a ring): the largest n with n^3 below 2^63, so that every
    /// cell row fits an int64 index.
    constexpr std::size_t kMaxCellsAxis = 2097151;

    /**
     *  @brief Build the box grid of a search.
     *
     *  The extent is search.box when set, otherwise the bounding box of the
     *  tracers read and of the endpoints of their valid displacements. It is
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

    /**
     *  @brief A field on a box grid at a point, as BasinCatalog::fieldCentre
     *  describes: trilinear between the cell centres, the point first
     *  clamped per axis to the range of the cell centres.
     *
     *  The cell centres are those of box_centres. On an axis with more
     *  than one cell, the point falls between centres i and i + 1, with
     *  weights 1 - f and f, f = (x - centre_i) / cellSize, and f = 1 at the
     *  last centre, so that at a cell centre the result is that cell's
     *  value exactly. The eight terms are summed in ascending row order.
     *
     *  @param field nCells values, [cell]
     *  @param x the point, x, y, z
     */
    double box_field_at (const DivergenceGrid& grid, const std::vector<double>& field, const double* x);

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

    /// The base-2 logarithm of the largest count of a cell the scale of a
    /// fixed-point accumulation of the divergence assumes; a larger count
    /// reruns the accumulation at the scale of the real count.
    constexpr int kAssumedCountBits = 20;

    /**
     *  @brief The scale exponent q = 61 - e - L of a fixed-point
     *  accumulation of the divergence, with B < 2^e, capped at 1000 so that
     *  2^-q times a sum stays a normal double; 0 when B is 0.
     *  @param B the largest absolute value of a term
     *  @param L the base-2 logarithm of the largest count of a cell
     */
    int scale_exponent (double B, int L);

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
     *  When none is, because x or x + d lies on a face, edge or vertex of
     *  the cell, the segment runs along a face, or no component of d
     *  exceeds 1e-9, the face follows the cells of box_cell, [lower, upper)
     *  on each axis: among the axes whose component of d is not zero, the
     *  plane of the cell the segment reaches along each, t = (plane - x) /
     *  d; the start cell is left through the plane of smallest t, the end
     *  cell entered through the plane of largest t. Equal t are decided as
     *  for x moved by an equal, vanishing step along +x, +y and +z: by the
     *  smaller -1 / d for the start cell, the larger for the end cell, then
     *  the lower axis. So a tracer on a face gives the cells and the terms
     *  of the same tracer moved one ulp into its cell, and every segment
     *  that leaves its cell has a face in each of its two cells. Axis -1
     *  only when d is zero.
     *
     *  @param cell the cell, (i, j, k)
     *  @param x the start of the segment, x, y, z
     *  @param d the displacement, x, y, z
     *  @param start true for the cell of x, which the segment leaves;
     *  false for the cell of x + d, which it enters
     */
    BoxFace box_face (const DivergenceGrid& grid, const std::array<std::size_t, 3>& cell,
                      const double* x, const double* d, bool start);

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
     *  it, an internal error that box_face rules out; then if the fill
     *  cannot reach every empty cell
     */
    BoxDivergence box_divergence (const BoxSearch& search, DivergenceGrid& grid);

    /**
     *  @brief Fill the divergence of the empty cells (nDisplacements 0) of
     *  a box grid, as DivergenceGrid::weight describes; their weight is
     *  left unchanged.
     *
     *  Reads only nDisplacements, divergence and the fields that define
     *  the cells.
     *
     *  @return the number of passes run, 0 with no empty cell
     *  @exception Error if a pass fills no cell while some are empty
     */
    unsigned box_fill (DivergenceGrid& grid);

    /**
     *  @brief The cells of the tracers a box search reads, without
     *  displacements: what box_divergence finds of them.
     *
     *  Sets DivergenceGrid::nTracers, counted in integers so that the
     *  result does not depend on the number of threads, and the cell of
     *  each tracer read, -1 for the others. used is left empty, nUsed 0;
     *  scaleExponent and fillPasses are 0. Reads only the fields of grid
     *  that define the cells.
     *
     *  @param search a checked search whose tracers read all lie inside a
     *  cell of grid (check_box_from_divergence)
     */
    BoxDivergence box_tracer_cells (const BoxSearch& search, DivergenceGrid& grid);

    //==========================================================================

    /**
     *  @brief Smooth the divergence of a box grid into divergenceSmoothed,
     *  as Smoothing describes for a box, with weight as the tracer weight.
     *
     *  The mean over the cube of cells around a cell is computed as three
     *  1-D convolutions of weight times divergence and of weight, along k,
     *  j, then i, each sum taken from 0. in ascending order; the result
     *  does not depend on the number of threads. The half-width of the cube
     *  is floor(3.5 sigmaLength / cellSize), or the nearest integer when
     *  3.5 sigmaLength / cellSize lies within 1e-9 relative of it, capped
     *  per axis at nCellsAxis - 1. A cell whose sum of weights is at most
     *  threshold takes 0. Any nCellsAxis of at least 1 per axis.
     *
     *  With sigmaLength > 0, divergenceSmoothed receives the smoothed
     *  field, and weightSmoothed the sum of weight over the cube from the
     *  same convolutions, divided by the product S_i S_j S_k, S_a the sum
     *  of the 1-D kernel over the cells of the axis within the half-width,
     *  cut at the faces, each taken from 0. in ascending order. With
     *  sigmaLength 0 divergenceSmoothed and weightSmoothed are emptied.
     *  divergence and weight are never changed.
     *
     *  Reads only divergence, weight and the fields that define the cells.
     *
     *  @param sigmaLength the sigma in coordinate units: 0 for no
     *  smoothing, and when positive with (2 sigmaLength) sigmaLength
     *  positive and finite (check_box checks it)
     *  @param threshold the sum of weights at or below which a cell takes
     *  0: the searches pass 1e-12 times the largest weight of the grid
     */
    void box_smooth (DivergenceGrid& grid, double sigmaLength, double threshold);

    //==========================================================================

    /// The basins of one kind found by the watershed, numbered by ascending
    /// seed cell.
    struct BasinSet {

      /// The row of the seed cell of each basin, ascending.
      std::vector<std::int64_t> seed;

      /// The cells of each basin: the extended basin.
      std::vector<std::uint64_t> nCells;

      /// The cells of the core of each basin: divergence below 0 for a
      /// void, at or above 0 for an overdense basin.
      std::vector<std::uint64_t> nCoreCells;

      /// The sum of DivergenceGrid::nTracers over the core cells of each
      /// basin; empty when the grid has no nTracers.
      std::vector<std::uint64_t> nCoreTracers;
    };

    /// The voids and the overdense basins of a box grid.
    struct BoxBasins {
      BasinSet voids;
      BasinSet overdense;
    };

    /// The cells a walk of the watershed records per thread. A longer walk
    /// labels its first kPathCells cells; the others are labelled when the
    /// loop over the cells reaches them, with the same result.
    constexpr std::size_t kPathCells = 4096;

    /// The cells of a block of the seed search: blocks of fixed size, so that
    /// the seeds are gathered in the same order with any number of threads.
    constexpr std::int64_t kSeedBlockRows = std::int64_t(1) << 16;

    /**
     *  @brief Split a box grid into basins by steepest descent of its
     *  divergence, as BasinCatalog describes.
     *
     *  Each cell flows to the first neighbour, among the 26 cut at the
     *  faces of the grid and taken in ascending row order, whose
     *  divergence is strictly below that of the cell and of every
     *  neighbour before it; a cell with no such neighbour is a seed. A
     *  seed below 0 starts a void, one at or above 0 an overdense basin;
     *  each kind is numbered by ascending seed row. Every cell belongs to
     *  the basin its flow reaches. Any nCellsAxis of at least 1 per axis;
     *  with no cell below 0 there is no void, and no error. The result
     *  does not depend on the number of threads.
     *
     *  The field split is divergenceSmoothed when it is not empty,
     *  otherwise divergence: "divergence" below means that field.
     *
     *  Reads only divergence, divergenceSmoothed, nCells, nCellsAxis and,
     *  when not empty, nTracers. Writes voidIndex and overdenseIndex, one of which is -1
     *  in every cell, and every field of basins (nCoreTracers left empty
     *  without nTracers).
     */
    void box_watershed (DivergenceGrid& grid, BoxBasins& basins);

    //==========================================================================

    /**
     *  @brief Fill the catalogues of a box grid split by box_watershed, as
     *  BasinCatalog describes for a box.
     *
     *  The centre is computed from the divergence of the seed's neighbours,
     *  with the shift's sign reversed for overdense basins, cut at the faces
     *  of the grid and taken in ascending row order, with a
     *  sum of 0 over |div_n| read as 1. It is then clamped to the seed cell,
     *  per axis, and moved by the least steps of one ulp that make
     *  box_cell return the seed for it. The radii and deltaAv come from the
     *  counts of basins, edgeFlag and contamination from the faces of the
     *  cells counted in integers, delta2Mps from the tracers read (cell of
     *  div at or above 0) within 2 mps of the centre, the sphere closed.
     *  The result does not depend on the number of threads. Any
     *  nCellsAxis of at least 1 per axis.
     *
     *  The divergence is the field box_watershed split: divergenceSmoothed
     *  when it is not empty, otherwise divergence.
     *
     *  Reads divergence, divergenceSmoothed, voidIndex, overdenseIndex and
     *  the fields that define the cells of grid; with search, also mps of grid, cell of
     *  div and the tracers of search. Writes every field of both
     *  catalogues; centreSky is left empty, and without search so are
     *  deltaAv, delta2Mps and nTracers.
     *
     *  @param search the search of the tracers, or null for none; basins
     *  then needs no nCoreTracers
     *  @param div the divergence of search; read only with search
     *  @exception Error internal, if the search of the tracers around the
     *  centres fails, with the meshsearch error nested when there is one
     */
    void box_basins (const BoxSearch* search, const BoxDivergence* div, const DivergenceGrid& grid,
                     const BoxBasins& basins, BasinCatalog& voids, BasinCatalog& overdense);

    /**
     *  @brief Where each input tracer of a box search ended up, from the
     *  catalogues box_basins filled.
     *
     *  A tracer read gets its cell, the void or the overdense basin of that
     *  cell, and its distance from the centre of that basin; a tracer not
     *  read gets -1, -1 and NaN. inMask is all 1, inBox and used are those
     *  of search and div; used is empty when div.used is.
     */
    TracerAssignment box_assignment (const BoxSearch& search, const BoxDivergence& div, const DivergenceGrid& grid,
                                     const BasinCatalog& voids, const BasinCatalog& overdense);

    /**
     *  @brief The extent of the objects a new lightcone grid must hold.
     *
     *  Redshifts, distances and declinations are the smallest and largest
     *  over the objects; each distance is DistanceTable::distanceAt of the
     *  object's redshift, so distanceMin and distanceMax are those of
     *  redshiftMin and redshiftMax. The right ascensions are cut at their
     *  largest gap on the circle: raMin is the first value after the gap,
     *  in [0, 2 pi), and raMax the last before it, plus 2 pi when the gap
     *  does not contain RA = 0, so raMax may exceed 2 pi.
     */
    struct LightconeExtent {
      double redshiftMin, redshiftMax;
      double distanceMin, distanceMax;
      double raMin, raMax;
      double decMin, decMax;
    };

    /**
     *  @brief Measure the extent of the kept tracers and the randoms.
     *
     *  The largest gap of the right ascensions is found without sorting:
     *  the circle is split into k = min(n, 2^20) equal bins, n the number
     *  of objects, each keeping its smallest and largest value; the gaps
     *  between consecutive non-empty bins, the wrap-around one included,
     *  are exact, and the largest is the cut. The wrap-around gap wins a
     *  tie, then the first gap from RA = 0. When the largest gap is
     *  narrower than a bin, the objects cover the whole circle and every
     *  ring of the grid is full whatever the cut. The result does not
     *  depend on the number of threads.
     *
     *  @param tracersSky right ascension in [0, 2 pi), declination and
     *  redshift of each tracer, radians
     *  @param inMask,inCut a tracer is read when both are 1; empty for all
     *  @param randomsSky as tracersSky, every random read
     *  @param distances the table of the call; every redshift read must lie
     *  inside it
     */
    LightconeExtent lightcone_extent (const std::vector<double>& tracersSky, const std::vector<std::uint8_t>& inMask,
                                      const std::vector<std::uint8_t>& inCut, const std::vector<double>& randomsSky,
                                      const otswap::DistanceTable& distances);

    /**
     *  @brief A lightcone grid as built, with the layout of its cells.
     *
     *  Cells are numbered by shell, then by declination ring within the
     *  shell, then by right-ascension sector within the ring, all from 0;
     *  every index is an int64. Shell 0 is the central sphere when there is
     *  one: a single ring and a single sector over the full sky.
     */
    struct LightconeGrid {

      /// The grid: the record fields, mpsProfile, nCells and the per-cell
      /// limits, centre, centreSky, volume and masked.
      DivergenceGrid cells;

      /// Whether shell 0 is the central sphere.
      bool centralSphere = false;

      /// The distance edges of the shells, increasing, [shell + 1].
      std::vector<double> shellEdge;

      /// The distance of each shell centre, the mean of its edges; 0 for
      /// the central sphere. [shell]
      std::vector<double> shellCentre;

      /// The redshift of each shell centre, through the distance table.
      /// [shell]
      std::vector<double> shellRedshift;

      /// mps(z) at shellRedshift, in coordinate units. [shell]
      std::vector<double> shellMps;

      /// The first ring of each shell, and the number of rings last.
      /// [shell + 1]
      std::vector<std::int64_t> shellFirstRing;

      /// The declination edges of each ring, radians. [ring]
      std::vector<double> ringDecMin, ringDecMax;

      /// 1 for a ring whose centre lies within 1e-5 radians of a pole,
      /// polar caps included: one sector over the full circle, with no
      /// right-ascension faces. [ring]
      std::vector<std::uint8_t> ringPole;

      /// 1 for a ring whose sectors cover the full circle from RA = 0 to
      /// 2 pi, polar caps and the central sphere included. [ring]
      std::vector<std::uint8_t> ringFull;

      /// The first cell of each ring, and the number of cells last.
      /// [ring + 1]
      std::vector<std::int64_t> ringFirstCell;

      /// The redshifts of the first and the last shell edge (0 with the
      /// central sphere, through the table), and the smallest mps(z)
      /// between them, at the profile's nodes inside the range and at its
      /// ends, in coordinate units.
      double redshiftFirstEdge = 0., redshiftLastEdge = 0., mpsMin = 0.;
    };

    /**
     *  @brief Build a lightcone grid from its record.
     *
     *  The record is the distance, right-ascension and declination ranges,
     *  cellSizeMps and mpsProfile of @p record; nothing else of it is read.
     *  The same record, table and mask give the same grid, bit for bit.
     *
     *  Shells: with s the cell size in mps, zMin and zMax the redshifts of
     *  distanceMin and distanceMax, the first edge is distanceMin - 0.1 s
     *  mps(zMin); when distanceMin does not exceed (3 / (4 pi))^(1/3) s
     *  mps(zMin), a central sphere of that radius replaces it. Each shell
     *  then has width s mps(z), z the redshift at the previous edge plus
     *  half the previous width (for the first shell, a width of s mps at
     *  the first edge), until an edge reaches distanceMax + 0.1 s mps(zMax).
     *
     *  Rings: the declination range padded by 0.1 s mps / R on each side
     *  (R the shell centre, mps at its redshift), with a polar cap of
     *  half-width s mps / (sqrt(pi) R) where the range comes that close to
     *  a pole; the rest is split into round(...) rings of equal width.
     *  Edges within 1e-12 of 0 are set to 0.
     *
     *  Sectors: a ring is a full circle when its padded right-ascension arc
     *  (0.1 s mps / (R cos dec) on each side) comes within two cells of the
     *  circumference, then split into round(R cos dec 2 pi / (s mps))
     *  sectors from RA = 0; otherwise into ceil(arc / (s mps)) sectors
     *  centred on the padded range.
     *
     *  masked is 1 where the cell centre, right ascension folded into
     *  [0, 2 pi), falls on an unobserved pixel of @p mask.
     *
     *  @param reused whether the record is a grid given for reuse; the
     *  messages then name "the grid to reuse"
     *  @exception Error, raised serially before any cell is filled: if an
     *  edge, a shell centre or a point where mps(z) is evaluated lies
     *  outside the distance table; if mps(z) is not positive inside the
     *  grid's range; if a shell would hold no declination ring, or a full
     *  ring no sector; if the shells, the rings of a shell or the sectors
     *  of a ring would exceed kMaxCellsAxis
     */
    LightconeGrid lightcone_grid (const DivergenceGrid& record, const otswap::DistanceTable& distances,
                                  const otswap::Mask& mask, bool reused = false);

    /**
     *  @brief The cell of a lightcone grid that holds a point.
     *
     *  Every lookup of a point in a lightcone grid goes through this
     *  function. A point on an edge belongs to the cell above it, in
     *  distance, declination and right ascension: a point on the first edge
     *  of the grid is inside, on the last outside, but for a declination of
     *  exactly pi/2 inside a north polar cap. In a full ring, the right
     *  ascension is folded into [0, 2 pi), so RA = 0 is in the first
     *  sector; in a partial ring, it is moved by a multiple of 2 pi next to
     *  the ring's centre.
     *
     *  @param distance the comoving distance, in coordinate units
     *  @param ra,dec radians; any finite right ascension
     *  @return the cell, or -1 outside every cell
     */
    std::int64_t lightcone_cell (const LightconeGrid& grid, double distance, double ra, double dec);

    /// Where a point lies with respect to one cell of a lightcone grid, per
    /// coordinate: -1 below the cell's lower edge, +1 on or above its upper
    /// edge, 0 inside its range. [0] distance, [1] declination, [2] right
    /// ascension.
    struct LightconeOutside {
      std::array<int, 3> side{};
    };

    /**
     *  @brief Where a point lies with respect to one cell of a lightcone
     *  grid, with the comparisons of lightcone_cell.
     *
     *  The three sides are all 0 exactly when lightcone_cell returns @p cell
     *  for the same point. Each coordinate is compared with the cell's own
     *  range, [lower, upper), whatever the other two: the distance with the
     *  edges of the cell's shell, the declination with those of its ring
     *  (a declination of exactly pi/2 is inside a ring whose upper edge is
     *  pi/2), the right ascension with those of its sector. Whether the
     *  right ascension is inside is decided as in lightcone_cell; when it is
     *  not, its side is taken with the right ascension moved by a multiple
     *  of 2 pi to within pi of the sector's centre, so that a point just
     *  below RA = 0 in a full ring is below the first sector.
     *
     *  @param cell a cell of the grid, 0 <= cell < nCells
     *  @param distance the comoving distance, in coordinate units
     *  @param ra,dec radians; any finite right ascension
     */
    LightconeOutside lightcone_outside (const LightconeGrid& grid, std::int64_t cell, double distance, double ra,
                                        double dec);

    /**
     *  @brief The sky coordinates of a Cartesian point, observer at the
     *  origin, exactly as otswap::toSky computes them: distance
     *  sqrt(x^2 + y^2 + z^2), right ascension atan2(y, x) folded into
     *  [0, 2 pi), declination asin(z / distance) with the ratio clamped to
     *  [-1, 1], with otswap's maths; at the origin both angles are 0.
     *
     *  @param x the point, x, y, z, in coordinate units; finite
     */
    void lightcone_sky (const double* x, double& distance, double& ra, double& dec);

    /**
     *  @brief The cell of a lightcone grid that holds a Cartesian point,
     *  observer at the origin.
     *
     *  The point is converted by lightcone_sky, as otswap::toSky does, and
     *  the cell is lightcone_cell's for those values. Every lookup of a
     *  Cartesian point goes through this function: endpoints of
     *  displacements, points along them and the clamp of the centres.
     *
     *  @param x the point, x, y, z, in coordinate units; finite
     *  @return the cell, or -1 outside every cell
     */
    std::int64_t lightcone_cell_at (const LightconeGrid& grid, const double* x);

    /// The face of a lightcone cell through which a displacement leaves or
    /// enters it.
    struct LightconeFace {

      /// 0 distance (a sphere), 1 declination (a cone, or the plane z = 0
      /// for an edge at 0), 2 right ascension (a plane through the z axis);
      /// -1 when no face is found.
      int coordinate = -1;

      /// -1 the cell's lower edge in that coordinate, +1 its upper edge.
      int side = 0;

      /// The crossing, x + t d, 0 <= t <= 1.
      double t = 0.;

      /// The crossing point x + t d, componentwise; x itself at t = 0 and
      /// x + d at t = 1.
      std::array<double, 3> point{};

      /// The unit normal of the face at point, outward from the cell (unit
      /// up to rounding).
      std::array<double, 3> normal{};
    };

    /**
     *  @brief The face of a lightcone cell through which the segment from
     *  x to x + d first leaves it (start) or last enters it (end).
     *
     *  Roots. Every root in (0, 1) of the surfaces bounding the cell is
     *  found, with no tolerance:
     *  - the spheres of the shell's two edges (only the outer one for the
     *    central sphere);
     *  - the cones of the ring's two declination edges, both nappes,
     *    (x^2 + y^2) sin^2 dec - z^2 cos^2 dec = 0; an edge at 0 is the
     *    plane z = 0; an edge at +-pi/2 is not searched;
     *  - the planes through the z axis of the sector's two right-ascension
     *    edges, whole planes, except in a ring with a single sector over
     *    the full circle (polar caps and the central sphere included). The
     *    edge 2 pi of a full ring is taken as 0, so that both sectors at
     *    RA = 0 see the same plane, bit for bit.
     *  A quadratic A t^2 + B t + C has roots q / A and C / q, with
     *  q = -(B / 2 + sign(B) sqrt(D)), sign(0) = +1, and D = (B / 2)^2 - AC
     *  written from L = x cross d: R^2 |d|^2 - |L|^2 for a sphere of radius
     *  R, with C = (|x| - R)(|x| + R); sin^2 dec (cos^2 dec (L_x^2 + L_y^2)
     *  - sin^2 dec L_z^2) for a cone. A surface gives roots only when D > 0;
     *  a plane only when its denominator is not 0; a root that is not
     *  finite, or not strictly inside (0, 1), is dropped; equal roots count
     *  once. The sines and cosines of the edges come from otswap's maths,
     *  computed in each call. A root that should be 0 or 1, for an endpoint
     *  on a face, may round just outside (0, 1) and be dropped; the
     *  intervals below give the same crossing either way.
     *
     *  Intervals. With the roots sorted, t = 0 and t = 1 bound the list;
     *  the cell cannot change inside an interval, so each is classified by
     *  lightcone_outside at its midpoint (a + b) / 2, the point x + t d
     *  converted by lightcone_sky.
     *  - start: when x itself lies outside the cell (by rounding, when the
     *    cell is that of the tracer's sky coordinates), the crossing is
     *    t = 0 and the point classified is x. Otherwise, from t = 0
     *    forwards, the first interval whose midpoint is outside the cell;
     *    the crossing is the root opening it. When none is, the crossing is
     *    t = 1 and the point classified is x + d.
     *  - end: from t = 1 backwards, the first interval whose midpoint is
     *    outside the cell; the crossing is the root closing it (t = 1 for
     *    the last interval). When none is, the crossing is t = 0 and the
     *    point classified is x; when x lies inside the cell and @p sky is
     *    given, the point classified is the tracer's sky position instead.
     *    That happens only when x lies within rounding of a face of the
     *    cell while its sky position lies in the cell across that face:
     *    the segment then enters the cell at x, through the face its sky
     *    position lies beyond, as the start search leaves the other cell
     *    there.
     *  The face is the coordinate the classified point violates, at the
     *  crossing, whichever surface gave the root. When it violates more than
     *  one, the first of distance-upper, declination-upper, right-ascension-
     *  upper, right-ascension-lower, declination-lower, distance-lower is
     *  taken (a vanishing step along +r, +dec, +ra, in that order of size).
     *  coordinate is -1 only when the classified point violates none: then
     *  every point examined lies in the cell.
     *
     *  Normal at the crossing point p, outward: distance, side * p / |p|;
     *  declination, side * (-sin dec_e p_x / rho, -sin dec_e p_y / rho,
     *  cos dec_e), dec_e the edge, rho = sqrt(p_x^2 + p_y^2) ((0, 0, side)
     *  for an edge at 0); right ascension, side * (-sin ra_e, cos ra_e, 0),
     *  ra_e the edge, as taken for the roots. These hold also when p is not on the face. At rho = 0
     *  (p on the z axis) the azimuth is 0, as lightcone_sky takes it.
     *
     *  Never throws: it is called inside parallel regions.
     *
     *  @param cell the cell searched: for start, the cell of the tracer's
     *  sky coordinates; for end, lightcone_cell_at(x + d). 0 <= cell < nCells
     *  @param x the start of the segment, the tracer's Cartesian position;
     *  finite
     *  @param d the displacement; finite
     *  @param start true for the cell left, false for the cell entered
     *  @param sky for the end search, the tracer's distance, right
     *  ascension and declination, from which its start cell was found; null
     *  for none. Not read by the start search.
     */
    LightconeFace lightcone_face (const LightconeGrid& grid, std::int64_t cell, const double* x, const double* d,
                                  bool start, const double* sky = nullptr);

    /// The tracers and displacements a lightcone divergence reads, after
    /// the entry checks. The pointers refer to arrays that outlive the
    /// call.
    struct LightconeSweep {

      std::size_t nTracers = 0;
      unsigned nRealizations = 0;

      /// [tracer][xyz], the Cartesian positions the displacements start
      /// from; rows of tracers not kept are not read.
      const std::vector<double>* tracers = nullptr;

      /// [tracer][3]: right ascension in [0, 2 pi), declination, redshift.
      const std::vector<double>* tracersSky = nullptr;

      /// [tracer]: a tracer is kept when both are 1.
      const std::vector<std::uint8_t>* inMask = nullptr;
      const std::vector<std::uint8_t>* inCut = nullptr;

      /// [realization][tracer][xyz] and [realization][tracer]; valid empty
      /// means all valid. With a new grid, every valid displacement of a
      /// kept tracer is finite; with a reused grid, those of the tracers
      /// read are checked here.
      const std::vector<double>* displacement = nullptr;
      const std::vector<std::uint8_t>* valid = nullptr;

      /// The distance table of the call; every kept redshift lies inside
      /// it.
      const otswap::DistanceTable* distances = nullptr;

      /// Whether the grid is a grid given for reuse, which kept tracers may
      /// lie outside; a new grid holds every kept tracer.
      bool reusedGrid = false;
    };

    /// What the divergence of a lightcone search finds besides the grid
    /// fields.
    struct LightconeDivergence {

      /// [tracer]: the cell of the tracer's sky position, distance
      /// DistanceTable::distanceAt of its redshift; -1 for a tracer not
      /// kept or outside the grid.
      std::vector<std::int64_t> cell;

      /// [tracer]: 1 for a kept tracer outside every cell of a reused grid,
      /// whose displacements are not read; empty with a new grid.
      std::vector<std::uint8_t> outsideGrid;

      /// [tracer]: 1 when a valid displacement of the tracer entered the
      /// divergence.
      std::vector<std::uint8_t> used;

      /// The tracers with used 1, those with outsideGrid 1, and the valid
      /// displacements of the tracers read that were rejected.
      std::size_t nUsed = 0, nOutsideGrid = 0, nRejected = 0;

      /// The largest norm of a valid displacement of a tracer read, in
      /// coordinate units.
      double maxDisplacement = 0.;

      /// The exponent q of the fixed-point accumulation: each term is
      /// rounded to a multiple of 2^-q.
      int scaleExponent = 0;
    };

    /**
     *  @brief The divergence of the displacement field on a lightcone grid,
     *  before the fill and the smoothing.
     *
     *  A tracer is read when it is kept and its sky position lies in a cell;
     *  DivergenceGrid::nTracers counts the tracers read in each cell, masked
     *  cells included. A valid displacement d of a tracer read, from its
     *  Cartesian position x, is rejected, and counts nowhere, when its start
     *  cell or its end cell (lightcone_cell_at(x + d)) is masked, when x + d
     *  lies outside every cell, or, without a central sphere, when the
     *  segment passes closer to the observer than the grid's first edge
     *  (|x cross d|^2 < r^2 |d|^2 with the closest point strictly inside the
     *  segment). An accepted displacement that ends in its start cell counts
     *  2 there and adds nothing. One that ends in another cell counts 1 in
     *  each and adds n . d / (s mps(z)) to each, n the outward normal at the
     *  crossing point of lightcone_face (start search in the start cell,
     *  end search in the end cell with the tracer's sky position), s
     *  cellSizeMps, z the redshift of the crossing point's distance clamped
     *  to the cell's distance range. The divergence of a cell is its sum
     *  over its count, 0 for a count of 0; its weight is the count over the
     *  number of realizations.
     *
     *  The sums are accumulated in fixed point as in a box, with the term
     *  bound maxDisplacement / (s mpsMin), so that the result is the same
     *  bit for bit with any number of threads.
     *
     *  @param grid the grid; its cells' nTracers, nDisplacements,
     *  divergence and weight are set
     *  @exception Error, in this order: an internal error if a kept
     *  tracer's distance cannot be read or the tracer lies outside a new
     *  grid; with a reused grid, if a valid displacement of a tracer read
     *  is not finite, or if no valid displacement belongs to a tracer read;
     *  after the sweep, an internal error if no face is found
     *  for an accepted displacement, or if a redshift or mps lookup fails at
     *  a crossing; if every valid displacement read was rejected
     */
    LightconeDivergence lightcone_divergence (const LightconeSweep& sweep, LightconeGrid& grid);

    /// The two neighbour definitions of a lightcone grid, each in compressed
    /// rows: the neighbours of cell c are contactNeighbour[contactFirst[c]]
    /// to contactNeighbour[contactFirst[c + 1] - 1], and likewise for the
    /// radius lists, in ascending order; a masked cell has no list in
    /// either.
    struct LightconeAdjacency {

      /// The contact lists: [cell + 1], and the neighbours of every cell,
      /// one list after the other. They may hold masked cells.
      std::vector<std::int64_t> contactFirst, contactNeighbour;

      /// The radius lists: [cell + 1], and the neighbours of every cell,
      /// one list after the other. They hold unmasked cells only.
      std::vector<std::int64_t> radiusFirst, radiusNeighbour;
    };

    /**
     *  @brief The contact lists and the radius lists of a lightcone grid.
     *
     *  Contact lists. The window of a cell is its own shell and the shells
     *  next to it, its declination range widened by 1 % of its width on
     *  each side, and its right-ascension range widened by 1 % of its
     *  width on each side; on a full ring the window is also shifted by
     *  2 pi, or by -2 pi, where it passes below the ring's first edge or
     *  above its last. A cell lies in the window when both its ranges
     *  overlap the window's, touching included. Two distinct cells are
     *  adjacent when either lies in the window of the other, so the
     *  relation is symmetric. Every pair is evaluated, masked cells
     *  included; the lists of masked cells are left empty, and the lists of
     *  the others may hold masked cells.
     *
     *  Radius lists. With side = cbrt(volume) and r = (sqrt(3) + 2) / 2,
     *  two distinct cells that are not masked are neighbours when the
     *  distance between their Cartesian centres is below
     *  r (side_i + side_j) / 2; the relation is symmetric, and the lists of
     *  masked cells are empty.
     *
     *  Reads shellFirstRing, ringDecMin, ringDecMax, ringFull,
     *  ringFirstCell and, of the cells, nCells, cellRaMin, cellRaMax,
     *  masked, centre and volume. The result is the same with any number of
     *  threads.
     *
     *  @exception Error if the search of the radius lists fails, an
     *  internal error with the meshsearch error nested
     */
    LightconeAdjacency lightcone_adjacency (const LightconeGrid& grid);

    /// What the fill of a lightcone grid did.
    struct LightconeFill {

      /// The cells filled, and the empty cells no pass could fill.
      std::size_t filled = 0, leftEmpty = 0;

      /// The passes that filled at least one cell; 0 with no empty cell.
      unsigned passes = 0;
    };

    /**
     *  @brief Fill the divergence of the empty cells of a lightcone grid,
     *  as DivergenceGrid::weight describes; the weights are left unchanged.
     *
     *  A cell is empty when its nDisplacements is 0 and it is not masked.
     *  In each pass, every empty cell whose radius neighbours include at least
     *  three that hold a value, a contribution or a value filled in an
     *  earlier pass, takes the sum of their divergence
     *  over r^2 divided by the sum of 1 / r^2, r the distance between the
     *  Cartesian centres, summed in ascending cell index; a pass reads only
     *  the values left by the previous one. The passes repeat until one
     *  fills no cell; the cells then left empty get a divergence of 0. The
     *  result is the same with any number of threads.
     *
     *  Reads nDisplacements, masked, centre and divergence of @p grid.
     *
     *  @param adjacency lightcone_adjacency of the grid; only its radius
     *  lists are read
     */
    LightconeFill lightcone_fill (const LightconeAdjacency& adjacency, DivergenceGrid& grid);

    /**
     *  @brief The smoothing sigma of each shell of a lightcone grid, in
     *  coordinate units.
     *
     *  In Mps, sigma times shellMps of the shell; in coordinate units,
     *  sigma for every shell. Checked serially, shells in ascending order:
     *  2 sigma_i^2 must be non-zero and finite, and 16 sigma_i^2 finite.
     *
     *  @return one sigma per shell; empty when sigma is 0
     *  @exception Error if a shell's sigma is too small or too large to
     *  build the kernel; in Mps the message names the first such shell and
     *  its redshift
     */
    std::vector<double> lightcone_sigma (const Smoothing& smoothing, const LightconeGrid& grid);

    /**
     *  @brief Smooth the divergence of a lightcone grid into
     *  divergenceSmoothed and weightSmoothed, as Smoothing describes for a
     *  lightcone.
     *
     *  For a cell m that is not masked, of shell i, the sums run over the
     *  cells n that are not masked and whose Cartesian centre lies at d^2
     *  < 16 sigma_i^2 from m's, m included, in ascending cell index: K =
     *  exp(-d^2 / (2 sigma_i^2)), num += K (weight_n divergence_n), den +=
     *  K weight_n, sum += K. divergenceSmoothed is num / den when den
     *  exceeds @p threshold, otherwise 0; weightSmoothed is den / sum. A
     *  masked cell takes 0 in both. With @p sigmaShell empty both fields
     *  are emptied. The result is the same with any number of threads.
     *
     *  Reads shellCentre, shellFirstRing, ringDecMin, ringDecMax,
     *  ringFirstCell and, of the cells, nCells, cellRaMin, cellRaMax,
     *  centre, centreSky, masked, divergence and weight.
     *
     *  @param sigmaShell lightcone_sigma of the grid
     *  @param threshold 1e-12 times the largest weight of the grid
     */
    void lightcone_smooth (LightconeGrid& grid, const std::vector<double>& sigmaShell, double threshold);

    /// The voids and the overdense basins of a lightcone grid.
    /// lightcone_watershed sets the seeds; lightcone_basins the counts and
    /// the volumes.
    struct LightconeBasins {
      BasinSet voids;
      BasinSet overdense;

      /// The volume of the core and of the extended basin of each void, in
      /// cubed coordinate units.
      std::vector<double> voidCoreVolume, voidVolume;

      /// The volume of the core and of the extended basin of each
      /// overdense basin, in cubed coordinate units.
      std::vector<double> overdenseCoreVolume, overdenseVolume;
    };

    /**
     *  @brief Split a lightcone grid into basins by steepest descent of its
     *  divergence.
     *
     *  Each cell that is not masked flows to the first of its radius
     *  neighbours, in ascending cell index, whose divergence is strictly
     *  below that of the cell and of every neighbour before it; a cell with
     *  no such neighbour is a seed. A seed below 0 starts a void, one at or
     *  above 0 an overdense basin; each kind is numbered by ascending seed
     *  cell index. Every cell that is not masked belongs to the basin its
     *  flow reaches; a masked cell is in none and takes -1 in voidIndex and
     *  overdenseIndex, as does a cell in a basin of the other kind. With no
     *  cell below 0 there is no void, and no error. The result does not
     *  depend on the number of threads.
     *
     *  The field split is divergenceSmoothed when it is not empty,
     *  otherwise divergence.
     *
     *  Reads masked, nCells and those fields, and the radius lists. Writes
     *  voidIndex and overdenseIndex; of the result only the seeds are set,
     *  the counts of each BasinSet and the volumes are left empty.
     *
     *  @param adjacency lightcone_adjacency of the grid; only its radius
     *  lists are read
     */
    LightconeBasins lightcone_watershed (const LightconeAdjacency& adjacency, DivergenceGrid& grid);

    /**
     *  @brief The area of the face two cells of a lightcone grid share.
     *
     *  The ranges of the two cells in distance, declination and right
     *  ascension (the second cell's moved by the multiple of 2 pi that
     *  brings its centre nearest to the first's) overlap with lengths Lr,
     *  Ld, La. With each tolerance 1e-9 times the first cell's own extent
     *  on that axis, the cells share a face when exactly one length is
     *  within its tolerance of 0 and the other two exceed theirs: a sphere
     *  of radius R, R^2 (sin dec_hi - sin dec_lo) La; a cone at
     *  declination D, (r_hi^2 - r_lo^2) / 2 cos D La; a half-plane of
     *  constant right ascension, (r_hi^2 - r_lo^2) / 2 Ld, over the
     *  overlapping ranges. Otherwise, a gap beyond a tolerance included,
     *  the area is 0.
     *
     *  Reads cellDistanceMin, cellDistanceMax, cellDecMin, cellDecMax,
     *  cellRaMin and cellRaMax.
     */
    double lightcone_shared_area (const DivergenceGrid& grid, std::size_t c, std::size_t n);

    /**
     *  @brief Count, measure and flag the basins lightcone_watershed found.
     *
     *  The cells of each basin are taken in ascending cell index and summed
     *  serially in that order, basins in parallel, so the result is the
     *  same with any number of threads. The core is the cells below 0 for
     *  a void, at or above 0 for an overdense basin, in the field the
     *  watershed split.
     *
     *  Per basin: the cells, the core cells and, when the grid has
     *  nTracers, the tracers of the core; the volume of the core and of
     *  the extended basin; radius and radiusExt, the radii of spheres of
     *  those volumes.
     *
     *  edgeFlag: a cell touches the boundary when it is on the border of
     *  the grid, or when a masked cell of its contact list shares a face
     *  with it (lightcone_shared_area above 0). The border is shell 0
     *  unless it is the central sphere, the last shell, the first and the
     *  last ring of a shell unless its edge is on the pole, and the first
     *  and the last sector of a ring that is not a full circle. 2 when a
     *  core cell touches it, 1 when only another cell does, 0 otherwise.
     *
     *  contamination: for each cell, its surface (the faces with area: no
     *  face on a pole or at the centre, no right-ascension face in a ring
     *  of one sector over the full circle) is split by its contact list
     *  into the area shared with its own basin, with masked cells, and with
     *  cells of other basins; the rest is on the border. Summed over the
     *  basin, the area on masked cells and on the border over that area
     *  plus the area shared with other basins; a small negative sum, above
     *  -1024 epsilon times the surface, is read as 0; 0 when the
     *  denominator is not positive.
     *
     *  Reads, of the grid, shellCentre, shellFirstRing, ringFull,
     *  ringFirstCell, centralSphere and, of the cells, nCells,
     *  masked, the per-cell limits, volume, voidIndex, overdenseIndex,
     *  divergence, divergenceSmoothed and nTracers; of the adjacency, the
     *  contact lists. Writes the counts and volumes of @p basins and, of
     *  each catalogue, nBasins, radius, radiusExt, edgeFlag, contamination
     *  and nTracers (empty without nTracers); centre, centreSky, deltaAv
     *  and delta2Mps are left empty.
     */
    void lightcone_basins (const LightconeGrid& grid, const LightconeAdjacency& adjacency, LightconeBasins& basins,
                           BasinCatalog& voids, BasinCatalog& overdense);

    /**
     *  @brief The centres of the basins lightcone_watershed found, Cartesian
     *  and on the sky.
     *
     *  For the seed cell of each basin, of centre x0: x0 + s sum_n div_n
     *  (x_n - x0) / sum_n |div_n| over its radius neighbours n, summed in
     *  ascending cell index, with x_n their centres, div the field the
     *  watershed split, a sum of 0 over |div_n| read as 1, and s = -1 for a
     *  void, +1 for an overdense basin, so that a void's centre moves
     *  towards lower divergence and an overdense basin's towards higher.
     *
     *  A point that lightcone_cell_at places in the seed cell is kept as
     *  it is. Otherwise its distance, declination and right ascension (the
     *  last moved by a multiple of 2 pi next to the cell's centre) are
     *  clamped to the cell's ranges, the declination alone in a polar cap
     *  and the distance alone in the central sphere, and converted back to
     *  Cartesian with otswap's maths; while lightcone_cell_at does not
     *  return the seed cell, every coordinate lightcone_outside finds
     *  outside moves inward by one std::nextafter step and the point is
     *  converted again. Every centre is then in its seed cell.
     *
     *  centreSky is lightcone_sky of the centre, so its right ascension is
     *  in [0, 2 pi), with the redshift of its distance from
     *  DistanceTable::redshiftAt. The result is the same with any number of
     *  threads.
     *
     *  Reads, of the grid, centralSphere, ringPole, ringFirstCell,
     *  shellFirstRing, the shell and ring edges and, of the cells, centre,
     *  the per-cell limits, divergence and divergenceSmoothed; of the
     *  adjacency, the radius lists; of @p basins, the seeds. Writes centre
     *  and centreSky of both catalogues.
     *
     *  @exception Error, an internal error raised after the loop, if a
     *  centre cannot be placed in its seed cell, or if the redshift lookup
     *  fails at a centre, the otswap error nested
     */
    void lightcone_centres (const LightconeGrid& grid, const LightconeAdjacency& adjacency,
                            const LightconeBasins& basins, const otswap::DistanceTable& distances,
                            BasinCatalog& voids, BasinCatalog& overdense);

    /**
     *  @brief The density contrasts of the basins of a lightcone search,
     *  deltaAv and delta2Mps.
     *
     *  With m = mps(z_c) from the grid's mpsProfile at the redshift of each
     *  centre: deltaAv = (nTracers / core volume) / m^-3 - 1; delta2Mps =
     *  n_t / (alpha n_r) - 1, with n_t and n_r the tracers and randoms
     *  within 2 m of the centre, the sphere closed, and alpha = N_t / N_r,
     *  computed once from the sizes of the two sets; NaN when n_r is 0.
     *  The counts use one meshsearch grid per set, built in turn, of cell
     *  4 mpsMin and limits the Cartesian bounding box of the grid's cells
     *  padded by the largest sphere radius. Basins are taken in parallel
     *  with integer counts, so the result is the same with any number of
     *  threads.
     *
     *  Reads, of the grid, cells (mpsProfile, nCells and the per-cell
     *  limits) and mpsMin; of @p basins, the core volumes; of each
     *  catalogue, nBasins, centre, centreSky and nTracers. Writes deltaAv
     *  and delta2Mps of both catalogues.
     *
     *  @param tracers the Cartesian positions of the tracers the search
     *  reads, [tracer][xyz]; their number is N_t
     *  @param randoms the Cartesian positions of the randoms the search
     *  reads, [random][xyz]; their number is N_r, at least 1
     *  @exception Error, an internal error raised after the loop, if the
     *  mps lookup fails at a centre, with the otswap error nested, or if
     *  meshsearch fails, with its error nested
     */
    void lightcone_densities (const LightconeGrid& grid, const LightconeBasins& basins,
                              const std::vector<double>& tracers, const std::vector<double>& randoms,
                              BasinCatalog& voids, BasinCatalog& overdense);

    /**
     *  @brief The lightcone search every lightcone entry converges into:
     *  the entry checks, the selection and the conversion of what the
     *  search reads.
     *
     *  @param tracers the Cartesian positions of the tracers, or null to
     *  convert the kept rows of tracersSky with otswap::toCartesian
     *  @param tracersSky the sky positions of the tracers; ignored, and
     *  computed from tracers with otswap::toSky, when randoms is not null
     *  @param randoms the Cartesian positions of the randoms, or null to
     *  convert the kept rows of randomsSky with otswap::toCartesian; when
     *  not null, the sky positions of tracers and randoms are computed
     *  from the Cartesian ones and randomsSky is ignored
     *  @param reconstruction the result the other arrays come from, already
     *  checked for consistency, or null for the array entries; its
     *  outsideMask and outsideRedshiftCut flags then leave tracers out, and
     *  its profile is not read with a grid
     *  @param grid a grid to reuse, or null to build one from the kept
     *  tracers and randoms. Its distance, right-ascension and declination
     *  ranges, cellSizeMps and mpsProfile define the cells, rebuilt with
     *  the call's distance table and mask; its nCells must equal the
     *  rebuilt count, and in the array entries its profile must equal
     *  mps in redshift and mps, bit for bit. Every other field is ignored.
     *  @exception Error on every case of the entry of findVoidsLightcone,
     *  of the grid and of the divergence; after the divergence, always,
     *  since the rest of the search is not written
     */
    Result runLightcone (const std::vector<double>* tracers, const std::vector<double>& tracersSky,
                         const std::vector<double>& displacement, const std::vector<std::uint8_t>& valid,
                         const std::vector<double>* randoms, const std::vector<double>& randomsSky,
                         const otswap::Mask& mask, const otswap::DistanceTable& distances,
                         const otswap::MpsProfile& mps, const Config& config,
                         const otswap::RedshiftCut& redshiftCut, const std::vector<std::int64_t>& ids,
                         const DivergenceGrid* grid, const otswap::Result* reconstruction);

  }

  /**
   *  @brief The grid findVoidsLightcone builds from this reconstruction and
   *  these randoms, bit for bit: on the tracers and randoms the search
   *  keeps, with the reconstruction's profile and distance table.
   *  Declared here until the pipeline runs end to end; it then moves to
   *  VoidFinder.h.
   *
   *  Only the fields that define the cells are set: geometry, nCells,
   *  distanceMin, distanceMax, raMin, raMax, decMin, decMax, cellSizeMps
   *  and mpsProfile; every per-cell field is left empty. The grid passes
   *  the checks of a grid to reuse, so that several searches can share
   *  it.
   *
   *  @param divergenceCellSize the cell size in units of mps, as
   *  Config::divergenceCellSize; without it, the default for the
   *  reconstruction's realizations, min(1, (20 / nRealizations)^(1/3))
   *  mps
   *  @exception Error on every case of the entry of the same overload of
   *  findVoidsLightcone, the checks of the displacements included,
   *  and of the grid's construction
   */
  DivergenceGrid lightconeGrid (const otswap::Result& reconstruction, const std::vector<double>& randomsSky,
                                const otswap::Mask& mask,
                                std::optional<double> divergenceCellSize = std::nullopt);

  /**
   *  @brief The grid the sky-position overload of findVoidsLightcone
   *  builds, bit for bit: the same selection of the tracers and randoms,
   *  as the overload above, without displacements.
   *
   *  @param divergenceCellSize the cell size in units of mps; finite and
   *  positive. A reasonable starting value is min(1, (20 / R)^(1/3)), R
   *  the number of realizations the search will read.
   *  @exception Error on every case of the entry of findVoidsLightcone
   *  that does not concern displacements, ids or the configuration, and
   *  of the grid's construction
   */
  DivergenceGrid lightconeGrid (const std::vector<double>& tracersSky, const std::vector<double>& randomsSky,
                                const otswap::Mask& mask, const otswap::DistanceTable& distances,
                                const otswap::MpsProfile& mps, double divergenceCellSize,
                                const otswap::RedshiftCut& redshiftCut = {});

  /**
   *  @brief Find the voids of a lightcone from sky positions.
   *  Declared here until the pipeline runs end to end; it then moves to
   *  VoidFinder.h.
   */
  Result findVoidsLightcone (const std::vector<double>& tracersSky,
                             const std::vector<double>& displacement,
                             const std::vector<std::uint8_t>& valid,
                             const std::vector<double>& randomsSky,
                             const otswap::Mask& mask,
                             const otswap::DistanceTable& distances,
                             const otswap::MpsProfile& mps,
                             const Config& config,
                             const otswap::RedshiftCut& redshiftCut = {},
                             const std::vector<std::int64_t>& ids = {},
                             const DivergenceGrid* grid = nullptr);

  /**
   *  @brief As findVoidsLightcone from sky positions, with the tracers and
   *  randoms given as Cartesian positions, the observer at the origin.
   *  Declared here until the pipeline runs end to end; it then moves to
   *  VoidFinder.h.
   *
   *  The sky coordinates of every row are computed first with
   *  otswap::toSky and the distance table given: the right ascension
   *  folded into [0, 2 pi), the declination, and the redshift at the
   *  row's distance from the origin. They alone decide the mask, the
   *  redshift cut and the cell of every tracer and random; the Cartesian
   *  positions are those the displacements start from. A position and its
   *  sky coordinates can then disagree on a cell only by rounding, for a
   *  point within rounding of a cell edge.
   *
   *  Every row must be finite, must not lie at the origin, and must lie
   *  within the distance table, rows on unobserved pixels or outside the
   *  cut included, since its sky coordinates decide whether it is kept. A
   *  table covering a wider redshift range than the catalogue's is a safe
   *  choice.
   *
   *  @param tracers 3 per tracer, x, y, z, in the unit of the distance
   *  table
   *  @param randoms as tracers, for the randoms; required
   *  @exception Error as findVoidsLightcone from sky positions, and if a
   *  row of tracers or randoms is not finite, lies at the origin, or lies
   *  at a distance outside the distance table; the message names the
   *  object
   */
  Result findVoidsLightconeFromCartesian (const std::vector<double>& tracers,
                                          const std::vector<double>& displacement,
                                          const std::vector<std::uint8_t>& valid,
                                          const std::vector<double>& randoms,
                                          const otswap::Mask& mask,
                                          const otswap::DistanceTable& distances,
                                          const otswap::MpsProfile& mps,
                                          const Config& config,
                                          const otswap::RedshiftCut& redshiftCut = {},
                                          const std::vector<std::int64_t>& ids = {},
                                          const DivergenceGrid* grid = nullptr);

  /**
   *  @brief As the functions above, from a lightcone reconstruction: its
   *  tracers, tracersSky, displacement, valid, ids, distance table and,
   *  without a grid, its mps profile.
   *
   *  A tracer is read only if neither outsideMask nor outsideRedshiftCut
   *  flags it and the mask given allows it; the mask must not be more
   *  permissive than the reconstruction's. The randoms are converted with
   *  the reconstruction's distance table and cut by its redshift cut
   *  (selection.redshiftCut). A mean read with
   *  otswap::io::readDisplacementMean is accepted.
   */
  Result findVoidsLightcone (const otswap::Result& reconstruction,
                             const std::vector<double>& randomsSky,
                             const otswap::Mask& mask,
                             const Config& config,
                             const DivergenceGrid* grid = nullptr);

}

#endif
