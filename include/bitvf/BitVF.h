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
 *  @file include/bitvf/BitVF.h
 *
 *  @brief The BitVF pipeline in one call, in a box or a lightcone: the
 *  optimal transport reconstruction, optionally the correction of
 *  redshift-space distortions and a second reconstruction in real space,
 *  and the void search.
 *
 *  The pipeline calls otswap and the void finder of bitvf/VoidFinder.h
 *  through their public functions only; each step can also be run by
 *  hand with them.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef BITVF_BITVF_H
#define BITVF_BITVF_H

#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

#include "otswap/OT.h"
#include "otswap/RSD.h"
#include "bitvf/VoidFinder.h"

namespace bitvf {

  // ==========================================================================
  // Configuration
  // ==========================================================================

  /**
   *  @brief Correction of redshift-space distortions in a box, as
   *  otswap::realSpaceBox applies it.
   *
   *  axis, redshift, bias and sigma have no default: each must be set.
   */
  struct BoxRsd {

    /// The line of sight: 0, 1 or 2, for x, y or z. Empty is an error.
    std::optional<unsigned> axis;

    /// The redshift of the box, at which the growth rate f is read from
    /// the distance table; inside the table's range.
    double redshift = std::numeric_limits<double>::quiet_NaN();

    /// The linear bias of the tracers, constant; finite and positive.
    double bias = std::numeric_limits<double>::quiet_NaN();

    /// Width of the neighbour average of the displacements, in coordinate
    /// units; 0 for none, finite and non-negative. For coordinates in
    /// Mpc/h, 10 Mpc/h is a reasonable starting value; the best value
    /// depends on the sample, and should be checked in each analysis.
    double sigma = std::numeric_limits<double>::quiet_NaN();

    /// The weighting of the neighbour average, and what the correction
    /// reports; passed to otswap::realSpaceBox unchanged.
    otswap::CorrectionConfig correction;

    /// Also search voids on the first, redshift-space, reconstruction
    /// (BitVFResult::zspaceVoids).
    bool zspaceVoids = false;
  };

  /**
   *  @brief The parameters of runBox.
   */
  struct BoxConfig {

    /// The reconstruction, used for both reconstructions of a run with
    /// correction, with the same seed. A seed of 0 is replaced by one
    /// drawn from std::random_device before the first reconstruction, and
    /// both use it; otswap::Result::config records it.
    otswap::Config ot;

    /// The box of the run, given to otswap::reconstructBox as limits on
    /// both reconstructions. Empty: the bounding box of the tracers, or of
    /// the tracers and the randoms when randoms are given. A point on a
    /// face is inside.
    std::optional<otswap::BoxLimits> limits;

    /// The mean particle separation, in coordinate units, given to both
    /// reconstructions. Empty: otswap computes it on each from the volume
    /// of the box and the tracers inside it (otswap::Result::mps).
    std::optional<double> mps;

    /// The correction of redshift-space distortions. Empty: one
    /// reconstruction and one search, no correction.
    std::optional<BoxRsd> rsd;

    /// The void search, on every pass. Its verbosity and log also govern
    /// the lines of the pipeline itself (runBox).
    Config voids;

    /// Keep the otswap results and the corrected catalogue in
    /// BitVFResult.
    bool keepReconstructions = false;
  };

  /// The randoms of the second, real-space, reconstruction of a lightcone
  /// run (runLightcone).
  enum class SecondPassRandoms {
    None,      ///< the first reconstruction's randoms, unchanged
    Resample,  ///< the first reconstruction's randoms, with redshifts drawn from the corrected tracers
    External   ///< randoms given to runLightcone for the second reconstruction
  };

  /**
   *  @brief Correction of redshift-space distortions in a lightcone, as
   *  otswap::realSpaceLightcone applies it.
   *
   *  bias and sigma have no default: each must be set.
   */
  struct LightconeRsd {

    /// The linear bias of the tracers as a function of redshift; it must
    /// give a positive bias at the redshift of every tracer kept.
    otswap::BiasTable bias;

    /// Width of the neighbour average of the displacements, in coordinate
    /// units; 0 for none, finite and non-negative. For coordinates in
    /// Mpc/h, 10 Mpc/h is a reasonable starting value; the best value
    /// depends on the sample, and should be checked in each analysis.
    double sigma = std::numeric_limits<double>::quiet_NaN();

    /// The weighting of the neighbour average, and what the correction
    /// reports; passed to otswap::realSpaceLightcone unchanged.
    otswap::CorrectionConfig correction;

    /// The randoms of the second reconstruction. Resample, the default,
    /// gives them the redshift distribution of the corrected tracers; None
    /// keeps the observed one; External takes randoms built by the user
    /// for the real-space sample.
    SecondPassRandoms secondPassRandoms = SecondPassRandoms::Resample;

    /// Also search voids on the first, redshift-space, reconstruction
    /// (BitVFResult::zspaceVoids).
    bool zspaceVoids = false;
  };

  /**
   *  @brief The parameters of runLightcone.
   */
  struct LightconeConfig {

    /// The reconstruction, used for both reconstructions of a run with
    /// correction, with the same seed. A seed of 0 is replaced by one
    /// drawn from std::random_device before the first reconstruction, and
    /// both use it; otswap::Result::config records it.
    otswap::Config ot;

    /// The redshift bins of the mps(z) profile of each reconstruction
    /// (otswap::reconstructLightcone); at least 1. 50 is a reasonable
    /// starting value.
    unsigned nBins = 0;

    /// The closed range of observed redshift of the run, given to both
    /// reconstructions and to every search: the tracers outside keep their
    /// row and are skipped, the randoms outside are removed. The default
    /// cuts nothing.
    otswap::RedshiftCut cut;

    /// The correction of redshift-space distortions. Empty: one
    /// reconstruction and one search, no correction.
    std::optional<LightconeRsd> rsd;

    /// The void search, on every pass. Its verbosity and log also govern
    /// the lines of the pipeline itself (runLightcone).
    Config voids;

    /// Keep the otswap results and the corrected catalogue in
    /// BitVFResult.
    bool keepReconstructions = false;
  };

  // ==========================================================================
  // Result
  // ==========================================================================

  /**
   *  @brief Output of runBox and runLightcone.
   *
   *  Given the same inputs and a seed other than 0, every field but the
   *  elapsedSeconds of the searches, of the kept otswap results and of the
   *  kept catalogue is the same bit for bit whatever the number
   *  of threads.
   */
  struct BitVFResult {

    /// The voids of the run: of the real-space pass with a correction, of
    /// the only pass without. Its per-tracer columns have one row per input tracer, in input order.
    Result voids;

    /// The voids of the redshift-space pass, when BoxRsd::zspaceVoids or
    /// LightconeRsd::zspaceVoids is set: the search of zspaceReconstruction
    /// (findVoidsBox, findVoidsLightcone).
    std::optional<Result> zspaceVoids;

    /// With keepReconstructions: the reconstruction the voids were found
    /// on. With a correction its rows are the second reconstruction's
    /// tracers only (runBox, runLightcone), and its ids are their identity
    /// among the input's.
    std::optional<otswap::Result> reconstruction;

    /// With keepReconstructions and a correction: the first,
    /// redshift-space, reconstruction, one row per input tracer.
    std::optional<otswap::Result> zspaceReconstruction;

    /// With keepReconstructions and a correction: the corrected
    /// catalogue, one row per input tracer.
    std::optional<otswap::RealSpaceCatalog> realSpace;
  };

  // ==========================================================================
  // Pipeline
  // ==========================================================================

  /**
   *  @brief Reconstruct a box and find its voids, without correction of
   *  redshift-space distortions.
   *
   *  The run:
   *
   *  1. Checks every argument and every field of config that can be
   *     checked before any work, so that an invalid run fails before the
   *     reconstruction.
   *  2. Fixes the box: config.limits, or the bounding box of the tracers,
   *     or of the tracers and the randoms when randoms are given.
   *  3. Without randoms, draws config.ot.nRealizations * N of them,
   *     uniformly in the box, N the tracers inside it (Result::randomSeed
   *     describes the draw).
   *  4. Reconstructs with otswap::reconstructBox: the randoms, the box as
   *     limits, config.mps when given.
   *  5. Searches with findVoidsBox on that reconstruction. The search reads
   *     otswap::Result::mps, so the reconstruction and the search use one
   *     mean particle separation.
   *
   *  The box always reaches otswap as limits, so Result::limits holds it,
   *  and the files record it, also when config.limits is empty. Its mean
   *  particle separation is then otswap's without limits, bit for bit,
   *  unless the randoms given reach beyond the tracers, which enlarges the
   *  box and so the mps.
   *
   *  @par Randoms drawn
   *  The seed is s = config.ot.seed + 1299709, modulo 2^32, a seed of 0
   *  counting as the one drawn to replace it (BoxConfig::ot); s is
   *  recorded in Result::randomSeed of every pass and written as the keyword RANDSEED
   *  of bitvf's files. One std::mt19937 seeded with s draws, for each
   *  random in turn, x, y and z, each as lo + u (hi - lo) over the box's
   *  range on that axis, with u = (a 2^26 + b) / 2^53, a the first output
   *  shifted right by 5 bits and b the second shifted right by 6; a value
   *  that rounds to hi or above is replaced by the largest double below
   *  hi, std::nextafter(hi, lo), with no further draw, so every random lies
   *  inside the box and the sequence of outputs is unchanged. s alone
   *  reproduces the randoms. They are given to otswap as given randoms;
   *  otswap's displacement files record config.ot.seed in SEED.
   *
   *  @par Identity
   *  ids are the tracers' identity: empty means the rows received. They go
   *  to otswap::Result::ids and Result::ids.
   *
   *  @par Reports
   *  Each step reports at its own configuration's level:
   *  config.ot.verbosity for the reconstruction, config.voids.verbosity
   *  for the search (Config::verbosity). At config.voids.verbosity
   *  Detailed, runBox adds, after drawing the randoms:
   *
   *      bitvf: drew <n> randoms in <box> with seed <s> in <t> s
   *
   *  and at Normal, last:
   *
   *      bitvf: runBox: <V> voids and <O> overdense basins from <N> tracers in <T> s
   *
   *  T the time of the whole call. Each line is written as soon as its step
   *  finishes, from the calling thread; a call that throws never writes the
   *  call line.
   *
   *  @param tracers Cartesian coordinates, in coordinate units, 3 * N
   *  entries.
   *  @param randoms Cartesian coordinates of the randoms, or empty to have
   *  them drawn. Given, they must number at least config.ot.nRealizations
   *  * N inside the box; those outside it are dropped by otswap. They go to
   *  the reconstruction only: the box search takes no randoms.
   *  @param config the parameters; config.rsd must be empty (the
   *  correction needs a distance table: use the overload below).
   *  @param ids the tracers' identity, or empty; one per tracer, no value
   *  twice.
   *
   *  @return the voids, and with config.keepReconstructions the
   *  reconstruction.
   *
   *  @exception Error if an argument or a field of config is invalid, as
   *  findVoidsBox and otswap::reconstructBox require, if config.rsd is
   *  set, if the bounding box spans no volume or has an infinite width on an axis, if randoms are to be drawn
   *  and no tracer lies inside the box, from otswap (prefixed "otswap: ",
   *  the original nested), and from the search, as findVoidsBox.
   */
  BitVFResult runBox (const std::vector<double>& tracers,
                      const std::vector<double>& randoms,
                      const BoxConfig& config,
                      const std::vector<std::int64_t>& ids = {});

  /**
   *  @brief Reconstruct a box, correct its redshift-space distortions,
   *  reconstruct it again in real space and find its voids.
   *
   *  Without config.rsd, as the overload above, and distances is not read.
   *  With it, the first four steps above, then:
   *
   *  - With BoxRsd::zspaceVoids, searches the first reconstruction
   *     (BitVFResult::zspaceVoids).
   *  - Builds the grid of the first reconstruction, boxGrid(first
   *     reconstruction, config.voids.divergenceCellSize): the grid of
   *     every search of the run, the one the redshift-space search builds.
   *  - Corrects the tracers with otswap::realSpaceBox (BoxRsd, distances).
   *  - Selects, in input order, the tracers corrected (status Corrected or
   *     MovedByNeighbours) whose corrected position lies in a cell of that
   *     grid. The others are excluded, with their row kept: those the
   *     correction left unmoved (NoValidNeighbour, and LeftOut, which
   *     includes the tracers outside the box) are flagged in
   *     Result::uncorrected, those moved outside every cell in
   *     Result::outsideGrid.
   *  - Reconstructs the selected tracers at their corrected positions,
   *     with the randoms of the first reconstruction (given or drawn), the
   *     same box as limits, config.mps when given, and the same seed.
   *     Their ids are the input's ids of the rows selected, or the rows
   *     themselves when no ids were given.
   *  - Searches that reconstruction on the grid built above, with
   *     findVoidsBox's grid argument: the same cells, and the smoothing and
   *     the delta2Mps radius from the first reconstruction's mean particle
   *     separation, which is Result::mps. The second reconstruction's own
   *     mps, recomputed by otswap from the fewer tracers when config.mps is
   *     empty, serves its reconstruction only. A tracer moved outside the
   *     box but into a cell is outside the box for the second
   *     reconstruction and is skipped by the search
   *     (TracerAssignment::inBox).
   *  - Returns the search with one row per input tracer.
   *
   *  @par Rows of Result::voids
   *  The selected rows carry the search's values. An excluded row carries
   *  -1 in every index, NaN in TracerAssignment::distance, 0 in
   *  TracerAssignment::used, and the first reconstruction's value in
   *  TracerAssignment::inBox (0 for a tracer outside BoxConfig::limits, 1
   *  otherwise); and one flag, uncorrected before outsideGrid. A tracer
   *  outside the box of the first reconstruction is uncorrected (LeftOut),
   *  and uncorrected means only that the correction did not move it. Result::ids are
   *  the ids given, or empty when none. nTracers is N; Result::counts
   *  counts the input rows: tracers N; tracersUncorrected and
   *  tracersOutsideGrid the rows flagged; tracersOutsideBox the rows outside
   *  the limits in the first reconstruction plus those the second search
   *  found outside among the selected tracers; tracersUsed and displacements
   *  as the second search counted them among the selected tracers.
   *
   *  @par Reports
   *  As the overload above, the correction at
   *  BoxRsd::correction::verbosity, and at config.voids.verbosity Detailed,
   *  after the selection of the corrected tracers:
   *
   *      bitvf: second pass on <k> of <N> tracers: <u> uncorrected, <g> outside the first pass's grid in <t> s
   *
   *  The second search's own lines count the selected tracers.
   *
   *  @param tracers as above, in redshift space.
   *  @param randoms as above; the first reconstruction's randoms are
   *  enough for the second, which has fewer tracers.
   *  @param distances a table carrying the growth rate at BoxRsd::redshift
   *  (otswap::rsdFactorBox).
   *  @param config as above.
   *  @param ids as above.
   *
   *  @return the voids, and as config asks the redshift-space voids, the
   *  two reconstructions and the corrected catalogue.
   *
   *  @exception Error as the overload above, but for config.rsd; if a
   *  field of config.rsd is invalid or the distance table cannot give the
   *  factor (otswap::rsdFactorBox), checked before the first
   *  reconstruction; if no tracer is selected for the second
   *  reconstruction; from otswap (prefixed "otswap: ", the original
   *  nested), and from each search, as findVoidsBox.
   */
  BitVFResult runBox (const std::vector<double>& tracers,
                      const std::vector<double>& randoms,
                      const otswap::DistanceTable& distances,
                      const BoxConfig& config,
                      const std::vector<std::int64_t>& ids = {});

  /**
   *  @brief Reconstruct a lightcone and find its voids, optionally
   *  correcting its redshift-space distortions and reconstructing it
   *  again in real space.
   *
   *  The run:
   *
   *  1. Checks every argument and every field of config that can be
   *     checked before any work, so that an invalid run fails before the
   *     reconstruction. With a correction, the factor of the correction
   *     is computed once at the observed redshifts of the tracers kept, so
   *     that a table without growth rate or an unusable bias fails here.
   *  2. Reconstructs with otswap::reconstructLightcone: the randoms, the
   *     mask, config.nBins, the distance table, config.ot and config.cut.
   *     The tracers on unobserved pixels or outside the cut take no part.
   *  3. Without config.rsd, searches that reconstruction with
   *     findVoidsLightcone (reconstruction, randomsSky, mask, config.voids)
   *     and returns.
   *
   *  With config.rsd, after the first reconstruction:
   *
   *  - Builds the grid of the first reconstruction, with
   *     config.voids.divergenceCellSize: the grid of every search of the
   *     run, the one the redshift-space search would build.
   *  - With LightconeRsd::zspaceVoids, searches the first reconstruction
   *     on it (BitVFResult::zspaceVoids).
   *  - Corrects the tracers with otswap::realSpaceLightcone (LightconeRsd,
   *     distances). The correction moves a tracer in redshift only.
   *  - Selects, in input order, the tracers corrected (status Corrected or
   *     MovedByNeighbours) whose corrected redshift lies inside config.cut
   *     and whose corrected position lies in a cell of the grid. The others
   *     are excluded, with their row kept, and flagged: those the
   *     correction left unmoved (NoValidNeighbour, and LeftOut, which
   *     includes the tracers on unobserved pixels or outside the cut) in
   *     Result::uncorrected, those moved outside the cut in
   *     Result::outsideCut, those moved outside every cell in
   *     Result::outsideGrid.
   *  - Reconstructs the selected tracers at their corrected redshifts, with
   *     the randoms of LightconeRsd::secondPassRandoms (below), the same
   *     mask, cut, bins and seed. Their ids are the input's ids of the rows
   *     selected, or the rows themselves when no ids were given.
   *  - Searches that reconstruction on the grid built above: the same cells,
   *     and the first reconstruction's mps(z) profile for the smoothing and
   *     the densities. The second reconstruction's own profile serves its
   *     reconstruction only.
   *  - Returns the search with one row per input tracer.
   *
   *  @par Rows of Result::voids
   *  The selected rows carry the search's values. An excluded row carries
   *  -1 in every index, NaN in TracerAssignment::distance and 0 in
   *  TracerAssignment::used, and exactly one flag, uncorrected before
   *  outsideCut before outsideGrid. TracerAssignment::inMask and inCut are
   *  those of the first reconstruction for every row: a tracer on an
   *  unobserved pixel or outside the cut is uncorrected (LeftOut) as well,
   *  and uncorrected means only that the correction did not move it.
   *  Result::ids are the ids given, or empty when none. nTracers is N;
   *  Result::counts counts the input rows: tracers N; tracersUncorrected,
   *  tracersOutsideCut and tracersOutsideGrid the rows flagged;
   *  tracersOutsideMask and tracersOutsideRedshiftCut those of the first
   *  reconstruction; tracersUsed, tracersInMaskedCells, displacements,
   *  displacementsRejected and the counts of the randoms as the second
   *  search counted them.
   *
   *  @par Second-pass randoms
   *  - None: the randoms given, unchanged; both reconstructions remove the
   *    same ones by the mask and the cut.
   *  - Resample: one random for each of the K randoms given that lie on
   *    observed pixels inside the cut, in row order, with the same right
   *    ascension and declination, bit for bit. Their redshifts are drawn
   *    from the k tracers selected after the correction: the seed is s =
   *    config.ot.seed + 104729, modulo 2^32, a seed of 0 counting as the
   *    one drawn to replace it; one std::mt19937 seeded with s gives, for
   *    each random in turn, otswap::uniformInt(engine, k) (Lemire's
   *    multiply-and-reject), and the random takes the corrected redshift of
   *    the selected tracer of that index, counted in row order. The draw is
   *    serial and takes nothing else from the engine, so s alone
   *    reproduces the redshifts. s is recorded in Result::randomSeed of
   *    BitVFResult::voids only. The search removes and counts the
   *    randoms that fall outside the grid.
   *  - External: @p randomsSkySecondPass, checked before the first
   *    reconstruction as @p randomsSky is (those on unobserved pixels or
   *    outside the cut are not read further). After the selection, those outside
   *    the grid are removed from the array the second reconstruction
   *    takes; the search takes the whole array and removes and counts them
   *    itself. At least config.ot.nRealizations times the k tracers
   *    selected must remain.
   *
   *  @par Identity
   *  ids are the tracers' identity: empty means the rows received. They go
   *  to otswap::Result::ids and Result::ids.
   *
   *  @par Reports
   *  Each step reports at its own configuration's level:
   *  config.ot.verbosity for the reconstructions,
   *  LightconeRsd::correction::verbosity for the correction,
   *  config.voids.verbosity for the searches (Config::verbosity). At
   *  config.voids.verbosity Detailed, runLightcone adds, after the selection of the corrected tracers:
   *
   *      bitvf: second pass on <k> of <N> tracers: <u> uncorrected, <c> outside the redshift cut, <g> outside the first pass's grid in <t> s
   *
   *  followed, with None, by
   *
   *      bitvf: second-pass randoms: the first pass's <K>
   *
   *  K the randoms on observed pixels inside the cut, or, with Resample, by
   *
   *      bitvf: resampled the redshifts of <K> randoms from <k> tracers with seed <s>
   *
   *  With External, at Normal, in place of that line:
   *
   *      bitvf: external randoms: <E> given, <m> on unobserved pixels, <c> outside the redshift cut, <g> outside the first pass's grid, <K> kept of <n> needed
   *
   *  and at Normal, last:
   *
   *      bitvf: runLightcone: <V> voids and <O> overdense basins from <N> tracers in <T> s
   *
   *  T the time of the whole call. Each line is written as soon as its step
   *  finishes, from the calling thread; a call that throws never writes the
   *  call line. The second search's own lines count the selected tracers.
   *
   *  @param tracersSky sky positions of the tracers, [tracer][ra, dec, z],
   *  radians, z the observed redshift, as findVoidsLightcone takes them.
   *  @param randomsSky sky positions of the randoms, as tracersSky;
   *  required. They describe the survey: every systematic belongs in them.
   *  The reconstruction needs at least config.ot.nRealizations times the
   *  tracers kept on observed pixels inside the cut.
   *  @param mask the observed pixels of the survey, given to both
   *  reconstructions and to every search.
   *  @param distances the redshift-distance table, carrying the growth
   *  rate when config.rsd is set.
   *  @param config the parameters.
   *  @param randomsSkySecondPass the randoms of the second reconstruction
   *  with SecondPassRandoms::External, as randomsSky; empty otherwise.
   *  @param ids the tracers' identity, or empty; one per tracer, no value
   *  twice.
   *
   *  @return the voids, and as config asks the redshift-space voids, the
   *  reconstructions and the corrected catalogue.
   *
   *  @exception Error
   *  - before the first reconstruction: on tracersSky, randomsSky, ids and
   *    config.cut as findVoidsLightcone; if no random is left on observed
   *    pixels inside the cut; if config.nBins or config.ot.nRealizations is
   *    0; on config.voids as findVoidsLightcone; if randomsSkySecondPass is
   *    given without a correction or with None or Resample, or is empty
   *    with External; if LightconeRsd::sigma is not finite and
   *    non-negative; with External, if randomsSkySecondPass is malformed,
   *    has a kept row with a redshift outside the distance table, or keeps
   *    no row; if the factor of the correction cannot be computed (the
   *    table has no growth rate, the bias is not positive), from otswap;
   *  - during the run: from otswap (prefixed "otswap: ", the original
   *    nested), a corrected distance outside the table included; from each
   *    search and from the grid of the run, as findVoidsLightcone; if no
   *    tracer is selected, giving the three counts; with
   *    External, if fewer randoms remain than needed, giving the counts.
   */
  BitVFResult runLightcone (const std::vector<double>& tracersSky,
                            const std::vector<double>& randomsSky,
                            const otswap::Mask& mask,
                            const otswap::DistanceTable& distances,
                            const LightconeConfig& config,
                            const std::vector<double>& randomsSkySecondPass = {},
                            const std::vector<std::int64_t>& ids = {});

}

#endif
