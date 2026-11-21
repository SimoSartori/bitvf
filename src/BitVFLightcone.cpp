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
 *  @file src/BitVFLightcone.cpp
 *
 *  @brief The lightcone pipeline, runLightcone.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <cmath>
#include <cstdint>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "otswap/Random.h"
#include "internal.h"


bitvf::internal::LightconeSecondPass bitvf::internal::lightcone_second_pass (const otswap::RealSpaceCatalog& real,
                                                                             const LightconeGrid& grid,
                                                                             const otswap::DistanceTable& distances,
                                                                             const otswap::RedshiftCut& cut,
                                                                             const std::vector<std::int64_t>& ids)
{
  const std::size_t N = real.nObjects;
  LightconeSecondPass pass;
  pass.uncorrected.assign(N, 0);
  pass.outsideCut.assign(N, 0);
  pass.outsideGrid.assign(N, 0);
  for (std::size_t i = 0; i < N; ++i) {
    const otswap::CorrectionStatus status = real.status[i];
    if (status != otswap::CorrectionStatus::Corrected && status != otswap::CorrectionStatus::MovedByNeighbours) {
      pass.uncorrected[i] = 1;
      ++pass.nUncorrected;
      continue;
    }
    const double* sky = &real.sky[3*i];
    if (!(sky[2] >= cut.min && sky[2] <= cut.max)) {
      pass.outsideCut[i] = 1;
      ++pass.nOutsideCut;
      continue;
    }
    const double distance = call_otswap([&] { return distances.distanceAt(sky[2]); });
    if (lightcone_cell(grid, distance, sky[0], sky[1]) < 0) {
      pass.outsideGrid[i] = 1;
      ++pass.nOutsideGrid;
      continue;
    }
    pass.rows.push_back(i);
    pass.tracersSky.insert(pass.tracersSky.end(), sky, sky + 3);
    pass.ids.push_back(ids.empty() ? (std::int64_t)i : ids[i]);
  }
  if (pass.rows.empty())
    throw Error("no tracer is left for the real-space reconstruction: of " + std::to_string(N) + " tracers, " +
                std::to_string(pass.nUncorrected) + " were left uncorrected, " + std::to_string(pass.nOutsideCut) +
                " were moved outside the redshift cut and " + std::to_string(pass.nOutsideGrid) +
                " were moved outside the first pass's grid");
  return pass;
}


// ============================================================================


std::vector<double> bitvf::internal::resample_randoms (const std::vector<double>& randomsSky,
                                                       const std::vector<std::size_t>& rows,
                                                       const std::vector<double>& tracersSky, const unsigned seed)
{
  const unsigned k = (unsigned)(tracersSky.size() / 3);
  std::mt19937 engine(seed);
  std::vector<double> out;
  out.reserve(3 * rows.size());
  for (const std::size_t row : rows) {
    const unsigned t = otswap::uniformInt(engine, k);
    out.push_back(randomsSky[3*row]);
    out.push_back(randomsSky[3*row+1]);
    out.push_back(tracersSky[3*t+2]);
  }
  return out;
}


// ============================================================================


namespace {

  /// The offset from the reconstruction's seed to the seed of the resampled randoms
  constexpr unsigned kRandomSeedOffset = 104729u;

  /**
   *  @brief Check every argument of runLightcone that can be checked before
   *  any work, in the order of the table of its errors, and select the rows
   *  the first reconstruction keeps.
   *  @return the redshifts of the kept tracers and the number of randoms kept
   */
  bitvf::internal::PipelineSelection check_run (const std::vector<double>& tracersSky,
                                                const std::vector<double>& randomsSky, const otswap::Mask& mask,
                                                const otswap::DistanceTable& distances,
                                                const bitvf::LightconeConfig& config,
                                                const std::vector<double>& randomsSkySecondPass,
                                                const std::vector<std::int64_t>& ids, const std::string& fn,
                                                bitvf::internal::ExternalSelection& external)
  {
    using bitvf::Error;

    const bitvf::internal::PipelineSelection selected =
      bitvf::internal::check_pipeline_sky(tracersSky, randomsSky, mask, distances, config.cut, ids, fn);
    if (config.nBins == 0) throw Error("the number of bins is 0; it must be at least 1");
    if (config.ot.nRealizations == 0) throw Error("the number of realizations is 0; it must be at least 1");
    bitvf::internal::check_config(config.voids, false);

    if (!randomsSkySecondPass.empty() && !config.rsd)
      throw Error("second-pass randoms were given, but no correction of redshift-space distortions is requested");
    if (config.rsd) {
      const bitvf::LightconeRsd& rsd = *config.rsd;
      if (rsd.secondPassRandoms == bitvf::SecondPassRandoms::External && randomsSkySecondPass.empty())
        throw Error("the second-pass randoms are External but none were given");
      if (rsd.secondPassRandoms != bitvf::SecondPassRandoms::External && !randomsSkySecondPass.empty())
        throw Error("second-pass randoms were given, but their mode is None or Resample, which draws its own");
      if (!std::isfinite(rsd.sigma) || rsd.sigma < 0.)
        throw Error("the smoothing width is " + std::to_string(rsd.sigma) +
                         " in coordinate units; it must be finite and non-negative");
      if (rsd.secondPassRandoms == bitvf::SecondPassRandoms::External)
        external = bitvf::internal::select_external_randoms(randomsSkySecondPass, mask, distances, config.cut, fn);
      bitvf::internal::call_otswap([&] { return otswap::rsdFactor(selected.redshift, distances, rsd.bias); });
    }

    return selected;
  }

}


// ============================================================================


bitvf::BitVFResult bitvf::runLightcone (const std::vector<double>& tracersSky, const std::vector<double>& randomsSky,
                                        const otswap::Mask& mask, const otswap::DistanceTable& distances,
                                        const LightconeConfig& config,
                                        const std::vector<double>& randomsSkySecondPass,
                                        const std::vector<std::int64_t>& ids)
{
  const std::string fn = "runLightcone";
  internal::Report report(config.voids.verbosity, config.voids.log);

  internal::ExternalSelection external;
  const internal::PipelineSelection selected =
    check_run(tracersSky, randomsSky, mask, distances, config, randomsSkySecondPass, ids, fn, external);

  otswap::Config ot = config.ot;
  if (ot.seed == 0) {
    std::random_device device;
    unsigned drawn = device();
    while (drawn == 0) drawn = device();
    ot.seed = drawn;
  }

  auto reconstruct = [&] (const std::vector<double>& sky, const std::vector<double>& randoms) {
    return internal::call_otswap([&] {
      return otswap::reconstructLightcone(sky, randoms, mask, config.nBins, distances, ot, config.cut);
    });
  };

  otswap::Result first = reconstruct(tracersSky, randomsSky);
  first.ids = ids;

  BitVFResult out;
  if (!config.rsd) {
    out.voids = findVoidsLightcone(first, randomsSky, mask, config.voids);
    if (config.keepReconstructions) out.reconstruction = std::move(first);
  }
  else {
    const LightconeRsd& rsd = *config.rsd;
    const DivergenceGrid grid = lightconeGrid(first, randomsSky, mask, config.voids.divergenceCellSize);
    Config voids = config.voids;
    voids.divergenceCellSize.reset();
    if (rsd.zspaceVoids) out.zspaceVoids = findVoidsLightcone(first, randomsSky, mask, voids, &grid);

    otswap::RealSpaceCatalog real = internal::call_otswap([&] {
      return otswap::realSpaceLightcone(first, distances, rsd.bias, rsd.sigma, rsd.correction);
    });

    report.lap();
    const internal::LightconeGrid lc = internal::lightcone_grid(grid, distances, mask, true);
    internal::LightconeSecondPass pass = internal::lightcone_second_pass(real, lc, distances, config.cut, ids);
    const std::size_t N = tracersSky.size() / 3;
    if (report.detailed()) {
      std::ostringstream line = internal::Report::line();
      line << "second pass on " << pass.rows.size() << " of " << N << (N == 1 ? " tracer: " : " tracers: ")
           << pass.nUncorrected << " uncorrected, " << pass.nOutsideCut << " outside the redshift cut, "
           << pass.nOutsideGrid << " outside the first pass's grid";
      report.step(line);
    }

    std::vector<double> resampled;
    std::optional<unsigned> randomSeed;
    if (rsd.secondPassRandoms == SecondPassRandoms::Resample) {
      randomSeed = ot.seed + kRandomSeedOffset;
      resampled = internal::resample_randoms(randomsSky, selected.randomRows, pass.tracersSky, *randomSeed);
      if (report.detailed()) {
        const std::size_t k = pass.rows.size();
        std::ostringstream line = internal::Report::line();
        line << "resampled the redshifts of " << selected.nRandomsKept << " randoms from " << k
             << (k == 1 ? " tracer" : " tracers") << " with seed " << *randomSeed;
        report.write(line.str());
      }
    }
    std::vector<double> reduced;
    if (rsd.secondPassRandoms == SecondPassRandoms::External) {
      const std::vector<std::size_t> inside = internal::external_in_grid(randomsSkySecondPass, external.rows, lc, distances);
      const std::size_t needed = (std::size_t)config.ot.nRealizations * pass.rows.size();
      const std::size_t outsideGrid = external.rows.size() - inside.size();
      if (report.normal()) {
        std::ostringstream line = internal::Report::line();
        line << "external randoms: " << external.nGiven << " given, " << external.outsideMask
             << " on unobserved pixels, " << external.outsideCut << " outside the redshift cut, " << outsideGrid
             << " outside the first pass's grid, " << inside.size() << " kept of " << needed << " needed";
        report.write(line.str());
      }
      if (inside.size() < needed)
        throw Error("the external randoms leave " + std::to_string(inside.size()) + " of the " +
                    std::to_string(needed) + " needed (the number of realizations times the " +
                    std::to_string(pass.rows.size()) + " tracers the mask keeps): of " +
                    std::to_string(external.nGiven) + " given, " + std::to_string(external.outsideMask) +
                    " lie on unobserved pixels and " + std::to_string(external.outsideCut) +
                    " outside the redshift cut (removed by the reconstruction), and " +
                    std::to_string(outsideGrid) + " outside the first pass's grid");
      reduced.reserve(3 * inside.size());
      for (const std::size_t r : inside)
        reduced.insert(reduced.end(), randomsSkySecondPass.begin() + 3*r, randomsSkySecondPass.begin() + 3*r + 3);
    }
    else if (!randomSeed && report.detailed()) {
      std::ostringstream line = internal::Report::line();
      line << "second-pass randoms: the first pass's " << selected.nRandomsKept;
      report.write(line.str());
    }
    const bool external2 = rsd.secondPassRandoms == SecondPassRandoms::External;
    const std::vector<double>& randoms2 = external2 ? randomsSkySecondPass : randomSeed ? resampled : randomsSky;
    const std::vector<double>& otRandoms2 = external2 ? reduced : randoms2;

    otswap::Result second = reconstruct(pass.tracersSky, otRandoms2);
    second.ids = std::move(pass.ids);
    out.voids = findVoidsLightcone(second, randoms2, mask, voids, &grid);
    out.voids.randomSeed = randomSeed;

    if (config.keepReconstructions) {
      out.reconstruction = std::move(second);
      out.zspaceReconstruction = std::move(first);
      out.realSpace = std::move(real);
    }
  }

  throw Error(fn + ": the result of the run is not implemented yet");
}
