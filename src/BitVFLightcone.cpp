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
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

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


namespace {

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
                                                const std::vector<std::int64_t>& ids, const std::string& fn)
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
        throw Error(fn + ": external second-pass randoms are not implemented yet");
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

  const internal::PipelineSelection selected =
    check_run(tracersSky, randomsSky, mask, distances, config, randomsSkySecondPass, ids, fn);

  otswap::Config ot = config.ot;
  if (ot.seed == 0) {
    std::random_device device;
    unsigned drawn = device();
    while (drawn == 0) drawn = device();
    ot.seed = drawn;
  }

  auto reconstruct = [&] (const std::vector<double>& sky) {
    return internal::call_otswap([&] {
      return otswap::reconstructLightcone(sky, randomsSky, mask, config.nBins, distances, ot, config.cut);
    });
  };

  otswap::Result first = reconstruct(tracersSky);
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

    if (rsd.secondPassRandoms == SecondPassRandoms::Resample)
      throw Error(fn + ": resampled second-pass randoms are not implemented yet");
    if (report.detailed()) {
      std::ostringstream line = internal::Report::line();
      line << "second-pass randoms: the first pass's " << selected.nRandomsKept;
      report.write(line.str());
    }

    otswap::Result second = reconstruct(pass.tracersSky);
    second.ids = std::move(pass.ids);
    out.voids = findVoidsLightcone(second, randomsSky, mask, voids, &grid);

    if (config.keepReconstructions) {
      out.reconstruction = std::move(second);
      out.zspaceReconstruction = std::move(first);
      out.realSpace = std::move(real);
    }
  }

  throw Error(fn + ": the result of the run is not implemented yet");
}
