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
 *  @file src/BitVFBox.cpp
 *
 *  @brief The box pipeline, runBox.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "bitvf/BitVF.h"
#include "internal.h"

namespace {

  /// The offset from the reconstruction's seed to the seed of the randoms
  /// runBox draws.
  constexpr unsigned kRandomSeedOffset = 1299709u;

  /**
   *  @brief Check every argument of runBox that can be checked before any
   *  work, in the order of the table of its errors, and fix the box.
   *  @param distances the distance table, or null for the overload without
   *  @return the box of the run
   */
  otswap::BoxLimits check_run (const std::vector<double>& tracers, const std::vector<double>& randoms,
                               const otswap::DistanceTable* distances, const bitvf::BoxConfig& config,
                               const std::vector<std::int64_t>& ids)
  {
    using bitvf::Error;

    bitvf::internal::check_coordinates(tracers, "tracers", false);
    bitvf::internal::check_coordinates(randoms, "the random array", true);
    bitvf::internal::check_ids(ids, tracers.size() / 3, "ids");
    if (config.limits) bitvf::internal::check_limits(*config.limits, "the box limits");
    if (config.mps && (!std::isfinite(*config.mps) || *config.mps <= 0.))
      throw Error("mps is " + std::to_string(*config.mps) + " in coordinate units; it must be finite and positive");
    if (config.ot.nRealizations == 0)
      throw Error("the number of realizations is 0; it must be at least 1");
    bitvf::internal::check_config(config.voids, false);

    if (config.rsd) {
      if (distances == nullptr)
        throw Error("the correction of redshift-space distortions (config.rsd) needs a distance table carrying "
                    "the growth rate; call runBox with distances");
      const bitvf::BoxRsd& rsd = *config.rsd;
      if (!rsd.axis)
        throw Error("the line-of-sight axis (BoxRsd::axis) is not set; it must be 0, 1 or 2");
      if (*rsd.axis > 2)
        throw Error("the line-of-sight axis is " + std::to_string(*rsd.axis) + "; it must be 0, 1 or 2");
      if (!std::isfinite(rsd.sigma) || rsd.sigma < 0.)
        throw Error("the smoothing width is " + std::to_string(rsd.sigma) +
                    " in coordinate units; it must be finite and non-negative");
      bitvf::internal::call_otswap([&] { return otswap::rsdFactorBox(rsd.redshift, *distances, rsd.bias); });
    }

    if (config.limits) return *config.limits;

    otswap::BoxLimits box;
    box.lower.fill(std::numeric_limits<double>::infinity());
    box.upper.fill(-std::numeric_limits<double>::infinity());
    for (const std::vector<double>* a : {&tracers, &randoms})
      for (std::size_t j = 0; j < a->size(); ++j) {
        box.lower[j%3] = std::min(box.lower[j%3], (*a)[j]);
        box.upper[j%3] = std::max(box.upper[j%3], (*a)[j]);
      }
    for (std::size_t k = 0; k < 3; ++k)
      if (!(box.lower[k] < box.upper[k]))
        throw Error(std::string(randoms.empty() ? "the tracers span" : "the tracers and the given randoms span") +
                    " no volume (all of them in one plane, on one line or at one point), so they fix no box; "
                    "give the box limits (BoxConfig::limits)");
    return box;
  }

  /**
   *  @brief Draw the randoms of a run given none, as runBox describes.
   *  @param seed the seed of the draw, the reconstruction's plus the offset
   *  @exception Error if no tracer lies inside the box, or if the randoms
   *  would not fit in a std::vector
   */
  std::vector<double> draw_randoms (const std::vector<double>& tracers, const otswap::BoxLimits& box,
                                    const unsigned nRealizations, const unsigned seed)
  {
    const std::size_t N = tracers.size() / 3;
    std::size_t inside = 0;
    for (std::size_t i = 0; i < N; ++i) {
      bool in = true;
      for (std::size_t k = 0; k < 3; ++k)
        in = in && tracers[3*i+k] >= box.lower[k] && tracers[3*i+k] <= box.upper[k];
      inside += in;
    }

    if (inside == 0) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "the box ";
      bitvf::internal::Report::box(out, box.lower, box.upper);
      out << " keeps 0 of " << N << (N == 1 ? " tracer" : " tracers")
          << "; at least one is needed to draw the randoms";
      throw bitvf::Error(out.str());
    }

    std::vector<double> randoms;
    if (inside > randoms.max_size() / 3 / nRealizations)
      throw bitvf::Error("the randoms to draw, " + std::to_string(nRealizations) + " x " + std::to_string(inside) +
                         ", do not fit in memory addressing; give fewer realizations or tracers");

    randoms.resize(3 * (std::size_t)nRealizations * inside);
    std::mt19937 rng(seed);
    for (std::size_t j = 0; j < randoms.size(); ++j) {
      const std::uint64_t a = (std::uint64_t)(rng() >> 5);
      const std::uint64_t b = (std::uint64_t)(rng() >> 6);
      randoms[j] = bitvf::internal::box_random(box.lower[j%3], box.upper[j%3],
                                               (double)((a << 26) + b) / 9007199254740992.);
    }
    return randoms;
  }

  /**
   *  @brief Spread the search of the selected rows over every input row,
   *  as runBox describes the rows of its result.
   *  @param rows the input row of each selected tracer, in order
   */
  void map_back (bitvf::Result& result, const std::vector<std::size_t>& rows, const std::size_t N,
                 const std::vector<std::int64_t>& ids, std::vector<std::uint8_t>&& uncorrected,
                 std::vector<std::uint8_t>&& outsideGrid, const std::size_t nUncorrected,
                 const std::size_t nOutsideGrid)
  {
    result.nTracers = N;
    result.ids = ids;
    result.uncorrected = std::move(uncorrected);
    result.outsideGrid = std::move(outsideGrid);
    result.counts.tracers = N;
    result.counts.tracersUncorrected = nUncorrected;
    result.counts.tracersOutsideGrid = nOutsideGrid;

    if (!result.assignment) return;
    const bitvf::TracerAssignment& from = *result.assignment;
    bitvf::TracerAssignment to;
    to.cellIndex.assign(N, -1);
    to.voidIndex.assign(N, -1);
    to.overdenseIndex.assign(N, -1);
    to.distance.assign(N, std::numeric_limits<double>::quiet_NaN());
    to.inMask.assign(N, 1);
    to.inBox.assign(N, 1);
    if (!from.used.empty()) to.used.assign(N, 0);
    for (std::size_t j = 0; j < rows.size(); ++j) {
      const std::size_t i = rows[j];
      to.cellIndex[i] = from.cellIndex[j];
      to.voidIndex[i] = from.voidIndex[j];
      to.overdenseIndex[i] = from.overdenseIndex[j];
      to.distance[i] = from.distance[j];
      to.inMask[i] = from.inMask[j];
      to.inBox[i] = from.inBox[j];
      if (!from.used.empty()) to.used[i] = from.used[j];
    }
    result.assignment = std::move(to);
  }

  /**
   *  @brief Both overloads of runBox.
   *  @param distances the distance table, or null for the overload without
   */
  bitvf::BitVFResult run_box (const std::vector<double>& tracers, const std::vector<double>& randoms,
                              const otswap::DistanceTable* distances, const bitvf::BoxConfig& config,
                              const std::vector<std::int64_t>& ids)
  {
    using namespace bitvf;

    internal::Report report(config.voids.verbosity, config.voids.log);
    const otswap::BoxLimits box = check_run(tracers, randoms, distances, config, ids);
    const std::size_t N = tracers.size() / 3;

    otswap::Config ot = config.ot;
    if (ot.seed == 0) {
      std::random_device device;
      unsigned drawn = device();
      while (drawn == 0) drawn = device();
      ot.seed = drawn;
    }

    std::optional<unsigned> randomSeed;
    std::vector<double> drawn;
    if (randoms.empty()) {
      randomSeed = ot.seed + kRandomSeedOffset;
      report.lap();
      drawn = draw_randoms(tracers, box, ot.nRealizations, *randomSeed);
      if (report.detailed()) {
        const std::size_t n = drawn.size() / 3;
        std::ostringstream out = internal::Report::line();
        out << "drew " << n << (n == 1 ? " random in " : " randoms in ");
        internal::Report::box(out, box.lower, box.upper);
        out << " with seed " << *randomSeed;
        report.step(out);
      }
    }
    const std::vector<double>& otRandoms = randoms.empty() ? drawn : randoms;

    auto reconstruct = [&] (const std::vector<double>& t) {
      return internal::call_otswap([&] {
        return config.mps ? otswap::reconstructBox(t, otRandoms, *config.mps, ot, box)
                          : otswap::reconstructBox(t, otRandoms, ot, box);
      });
    };

    otswap::Result first = reconstruct(tracers);
    first.ids = ids;

    BitVFResult out;
    if (!config.rsd) {
      out.voids = findVoidsBox(first, config.voids);
      out.voids.randomSeed = randomSeed;
      if (config.keepReconstructions) out.reconstruction = std::move(first);
    }
    else {
      const BoxRsd& rsd = *config.rsd;
      if (rsd.zspaceVoids) {
        out.zspaceVoids = findVoidsBox(first, config.voids);
        out.zspaceVoids->randomSeed = randomSeed;
      }

      const DivergenceGrid grid = boxGrid(first, config.voids.divergenceCellSize);
      otswap::RealSpaceCatalog real = internal::call_otswap([&] {
        return otswap::realSpaceBox(first, *rsd.axis, rsd.redshift, *distances, rsd.bias, rsd.sigma,
                                    rsd.correction);
      });

      report.lap();
      std::vector<std::uint8_t> uncorrected(N, 0), outsideGrid(N, 0);
      std::vector<std::size_t> rows;
      std::size_t nUncorrected = 0, nOutsideGrid = 0;
      for (std::size_t i = 0; i < N; ++i) {
        const otswap::CorrectionStatus status = real.status[i];
        if (status != otswap::CorrectionStatus::Corrected && status != otswap::CorrectionStatus::MovedByNeighbours) {
          uncorrected[i] = 1;
          ++nUncorrected;
        }
        else if (internal::box_cell(grid, &real.cartesian[3*i]) < 0) {
          outsideGrid[i] = 1;
          ++nOutsideGrid;
        }
        else rows.push_back(i);
      }
      if (rows.empty())
        throw Error("no tracer is left for the real-space reconstruction: of " + std::to_string(N) +
                    " tracers, " + std::to_string(nUncorrected) + " were left uncorrected and " +
                    std::to_string(nOutsideGrid) + " were moved outside the first pass's grid");

      const std::size_t k = rows.size();
      std::vector<double> moved(3 * k);
      std::vector<std::int64_t> movedIds(k);
      for (std::size_t j = 0; j < k; ++j) {
        for (std::size_t c = 0; c < 3; ++c) moved[3*j+c] = real.cartesian[3*rows[j]+c];
        movedIds[j] = ids.empty() ? (std::int64_t)rows[j] : ids[rows[j]];
      }
      if (report.detailed()) {
        std::ostringstream line = internal::Report::line();
        line << "second pass on " << k << " of " << N << (N == 1 ? " tracer: " : " tracers: ")
             << nUncorrected << " uncorrected, " << nOutsideGrid << " outside the first pass's grid";
        report.step(line);
      }

      otswap::Result second = reconstruct(moved);
      second.ids = std::move(movedIds);

      Config voids = config.voids;
      voids.divergenceCellSize.reset();
      out.voids = findVoidsBox(second, voids, &grid);
      out.voids.randomSeed = randomSeed;
      map_back(out.voids, rows, N, ids, std::move(uncorrected), std::move(outsideGrid), nUncorrected,
               nOutsideGrid);

      if (config.keepReconstructions) {
        out.reconstruction = std::move(second);
        out.zspaceReconstruction = std::move(first);
        out.realSpace = std::move(real);
      }
    }

    if (report.normal()) {
      std::ostringstream line = internal::Report::line();
      line << "runBox: ";
      internal::Report::basins(line, out.voids.voids.nBasins, out.voids.overdense.nBasins);
      line << " from " << N << (N == 1 ? " tracer " : " tracers ");
      internal::Report::time(line, report.elapsed());
      report.write(line.str());
    }
    return out;
  }

}

// ============================================================================


double bitvf::internal::box_random (const double lo, const double hi, const double u)
{
  const double x = lo + u * (hi - lo);
  return x >= hi ? std::nextafter(hi, lo) : x;
}


// ============================================================================


bitvf::BitVFResult bitvf::runBox (const std::vector<double>& tracers,
                                  const std::vector<double>& randoms,
                                  const BoxConfig& config,
                                  const std::vector<std::int64_t>& ids)
{
  return run_box(tracers, randoms, nullptr, config, ids);
}


// ============================================================================


bitvf::BitVFResult bitvf::runBox (const std::vector<double>& tracers,
                                  const std::vector<double>& randoms,
                                  const otswap::DistanceTable& distances,
                                  const BoxConfig& config,
                                  const std::vector<std::int64_t>& ids)
{
  return run_box(tracers, randoms, &distances, config, ids);
}
