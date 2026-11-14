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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;
  constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

  std::string num (const double x)
  {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out.precision(17);
    out << x;
    return out.str();
  }

  struct Selection {
    std::vector<std::uint8_t> inMask;
    std::vector<std::uint8_t> inCut;
    std::size_t outsideMask = 0;
    std::size_t outsideCut = 0;
  };

  struct Entry {
    std::size_t nTracers = 0;
    unsigned nRealizations = 0;
    std::vector<double> tracersSky;
    std::vector<double> tracers;
    std::vector<double> randomsSky;
    std::vector<double> randoms;
    Selection tracerSelection;
    Selection randomSelection;
    std::size_t nRandoms = 0;
    double maxDisplacement = 0.;
    bitvf::Counts counts;
  };

  std::size_t check_sky (const std::vector<double>& sky, const std::string& name, const std::string& fn)
  {
    if (sky.empty()) throw bitvf::Error(fn + ": " + name + " is empty");
    if (sky.size() % 3 != 0)
      throw bitvf::Error(fn + ": " + name + " holds " + std::to_string(sky.size()) +
                         " entries; it must hold 3 per object");
    return sky.size() / 3;
  }

  void check_sky_rows (const std::vector<double>& sky, const std::string& noun, const std::string& plural)
  {
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      for (std::size_t k = 0; k < 3; ++k)
        if (!std::isfinite(sky[3*i+k]))
          throw bitvf::Error("the sky positions of the " + plural + " hold a non-finite value at entry " +
                             std::to_string(3*i+k) + " (object " + std::to_string(i) + ", component " +
                             std::to_string(k) + ")");
      if (!(sky[3*i+1] >= -kPi/2 && sky[3*i+1] <= kPi/2))
        throw bitvf::Error("the declination of " + noun + " " + std::to_string(i) + " is " + num(sky[3*i+1]) +
                           " radians; it must lie in [-pi/2, pi/2]");
    }
  }

  void check_profile (const otswap::MpsProfile& profile, const std::string& fn)
  {
    const std::vector<double>& z = profile.redshift;
    const std::vector<double>& m = profile.mps;
    if (z.empty() && m.empty())
      throw bitvf::Error("the mps profile is empty; it must hold at least one entry");
    if (z.size() != m.size())
      throw bitvf::Error(fn + ": the redshift and mps lists of the profile hold " + std::to_string(z.size()) +
                         " and " + std::to_string(m.size()) + " entries");
    for (std::size_t i = 0; i < z.size(); ++i) {
      if (!std::isfinite(z[i]))
        throw bitvf::Error("the mps profile holds a redshift of " + num(z[i]) + " at entry " +
                           std::to_string(i) + "; it must be finite");
      if (!std::isfinite(m[i]))
        throw bitvf::Error("the mps profile holds a mean particle separation of " + num(m[i]) + " at entry " +
                           std::to_string(i) + "; it must be finite");
    }
    for (std::size_t i = 1; i < z.size(); ++i)
      if (!(z[i] > z[i-1]))
        throw bitvf::Error("the mps profile holds a redshift of " + num(z[i]) + " at entry " + std::to_string(i) +
                           ", not above the " + num(z[i-1]) + " before it; the redshifts must increase strictly");
    for (std::size_t i = 0; i < m.size(); ++i)
      if (!(m[i] > 0.))
        throw bitvf::Error("the mps profile holds a mean particle separation of " + num(m[i]) + " at entry " +
                           std::to_string(i) + "; it must be positive");
  }

  void fold_right_ascension (std::vector<double>& sky)
  {
    for (std::size_t i = 0; i + 2 < sky.size(); i += 3) {
      double ra = std::fmod(sky[i], kTwoPi);
      if (ra < 0.) ra += kTwoPi;
      if (ra >= kTwoPi) ra = 0.;
      sky[i] = ra;
    }
  }

  Selection select_rows (const std::vector<double>& sky, const otswap::Mask& mask,
                         const otswap::RedshiftCut& cut)
  {
    const std::int64_t n = (std::int64_t)(sky.size() / 3);
    Selection selection;
    selection.inMask.assign((std::size_t)n, 1);
    selection.inCut.assign((std::size_t)n, 1);
    std::size_t outsideMask = 0, outsideCut = 0;
    std::exception_ptr failure;
#pragma omp parallel for schedule(static) reduction(+:outsideMask, outsideCut)
    for (std::int64_t i = 0; i < n; ++i) {
      const std::size_t r = (std::size_t)i;
      try {
        const bool masked = !mask.allows(sky[3*r], sky[3*r+1]);
        const bool cutOut = !(sky[3*r+2] >= cut.min && sky[3*r+2] <= cut.max);
        if (masked) {
          selection.inMask[r] = 0;
          ++outsideMask;
        }
        if (cutOut) {
          selection.inCut[r] = 0;
          if (!masked) ++outsideCut;
        }
      }
      catch (...) {
#pragma omp critical(bitvf_lightcone_select)
        if (!failure) failure = std::current_exception();
      }
    }
    if (failure) {
      try {
        std::rethrow_exception(failure);
      }
      catch (const otswap::Error& e) {
        std::throw_with_nested(bitvf::Error(std::string("otswap: ") + e.what()));
      }
    }
    selection.outsideMask = outsideMask;
    selection.outsideCut = outsideCut;
    return selection;
  }

  void check_usable (const std::vector<double>& sky, const std::vector<double>* cartesian, const Selection& selection,
                     const otswap::DistanceTable& distances, const std::string& noun, const std::string& plural)
  {
    const double lo = distances.minRedshift(), hi = distances.maxRedshift();
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      if (!selection.inMask[i] || !selection.inCut[i]) continue;
      const double z = sky[3*i+2];
      if (!(z >= lo && z <= hi))
        throw bitvf::Error("the redshift of " + noun + " " + std::to_string(i) + " is " + num(z) +
                           ", outside the distance table [" + num(lo) + ", " + num(hi) + "]");
    }
    if (cartesian == nullptr) return;
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      if (!selection.inMask[i] || !selection.inCut[i]) continue;
      for (std::size_t k = 0; k < 3; ++k)
        if (!std::isfinite((*cartesian)[3*i+k]))
          throw bitvf::Error("the " + plural + " hold a non-finite value at entry " + std::to_string(3*i+k) +
                             " (object " + std::to_string(i) + ", component " + std::to_string(k) + ")");
    }
  }

  void keep_randoms (Entry& entry, const std::vector<double>& sky, const std::vector<double>* cartesian)
  {
    const Selection& selection = entry.randomSelection;
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      if (!selection.inMask[i] || !selection.inCut[i]) continue;
      entry.randomsSky.insert(entry.randomsSky.end(), sky.begin() + 3*i, sky.begin() + 3*i + 3);
      if (cartesian != nullptr)
        entry.randoms.insert(entry.randoms.end(), cartesian->begin() + 3*i, cartesian->begin() + 3*i + 3);
    }
    entry.nRandoms = entry.randomsSky.size() / 3;
  }

  std::string left_out (const std::size_t outsideMask, const std::size_t outsideCut)
  {
    return std::to_string(outsideMask) + " lie on unobserved pixels and " + std::to_string(outsideCut) +
           " outside the redshift cut";
  }

  void report_selection (const bitvf::internal::Report& report, const std::string& plural, const std::size_t total,
                         const Selection& selection, const bool boundedCut)
  {
    if (!report.normal()) return;
    const std::size_t kept = total - selection.outsideMask - selection.outsideCut;
    std::ostringstream out = bitvf::internal::Report::line();
    out << "kept " << kept << " of " << total << (total == 1 ? " " + plural.substr(0, plural.size() - 1) + ": " : " " + plural + ": ")
        << selection.outsideMask << " on unobserved pixels";
    if (boundedCut) out << ", " << selection.outsideCut << " outside the redshift cut";
    report.write(out.str());
  }

  Entry check_entry (const std::vector<double>* tracers, const std::vector<double>& tracersSky,
                     const std::vector<double>& displacement, const std::vector<std::uint8_t>& valid,
                     const std::vector<double>* randoms, const std::vector<double>& randomsSky,
                     const otswap::Mask& mask, const otswap::DistanceTable& distances,
                     const otswap::MpsProfile& mps, const bitvf::Config& config,
                     const otswap::RedshiftCut& cut, const std::vector<std::int64_t>& ids,
                     const bitvf::DivergenceGrid* grid, const bitvf::internal::Report& report)
  {
    using bitvf::Error;
    const std::string fn = "findVoidsLightcone";

    const std::size_t N = check_sky(tracersSky, "tracersSky", fn);
    const std::size_t M = check_sky(randomsSky, "randomsSky", fn);

    if (tracers != nullptr) {
      if (tracers->size() != tracersSky.size())
        throw Error(fn + ": tracers holds " + std::to_string(tracers->size()) + " entries for " +
                    std::to_string(N) + " objects in tracersSky");
      if (randoms->size() != randomsSky.size())
        throw Error(fn + ": randoms holds " + std::to_string(randoms->size()) + " entries for " +
                    std::to_string(M) + " objects in randomsSky");
    }

    const std::size_t R = bitvf::internal::check_displacement_shape(displacement, valid, N, fn);

    if (std::isnan(cut.min) || std::isnan(cut.max) || cut.min > cut.max)
      throw Error("the redshift cut is [" + num(cut.min) + ", " + num(cut.max) +
                  "]; neither bound may be NaN and the lower bound must not exceed the upper one");

    check_profile(mps, fn);
    bitvf::internal::check_ids(ids, N, fn, false);
    bitvf::internal::check_config(config, grid != nullptr);

    check_sky_rows(tracersSky, "tracer", "tracers");
    check_sky_rows(randomsSky, "random", "randoms");

    Entry entry;
    entry.nTracers = N;
    entry.nRealizations = (unsigned)R;
    entry.tracersSky = tracersSky;
    fold_right_ascension(entry.tracersSky);
    std::vector<double> randomsFolded = randomsSky;
    fold_right_ascension(randomsFolded);

    entry.tracerSelection = select_rows(entry.tracersSky, mask, cut);
    entry.randomSelection = select_rows(randomsFolded, mask, cut);

    check_usable(entry.tracersSky, tracers, entry.tracerSelection, distances, "tracer", "tracers");
    check_usable(randomsFolded, randoms, entry.randomSelection, distances, "random", "randoms");

    const Selection& ts = entry.tracerSelection;
    const std::int64_t rows = (std::int64_t)(R * N);
    std::size_t nValid = 0, nRead = 0, firstNonFinite = kNone;
    double maxDisplacement = 0.;
#pragma omp parallel for schedule(static) reduction(+:nValid, nRead) \
  reduction(min:firstNonFinite) reduction(max:maxDisplacement)
    for (std::int64_t row = 0; row < rows; ++row) {
      const std::size_t r = (std::size_t)row;
      const std::size_t i = r % N;
      if (!valid.empty() && !valid[r]) continue;
      ++nValid;
      if (!ts.inMask[i] || !ts.inCut[i]) continue;
      bool finite = true;
      for (std::size_t k = 0; k < 3; ++k)
        if (!std::isfinite(displacement[3*r+k])) {
          firstNonFinite = std::min(firstNonFinite, 3*r+k);
          finite = false;
        }
      if (!finite) continue;
      const double norm = std::sqrt(displacement[3*r]*displacement[3*r] + displacement[3*r+1]*displacement[3*r+1] +
                                    displacement[3*r+2]*displacement[3*r+2]);
      maxDisplacement = std::max(maxDisplacement, norm);
      ++nRead;
    }
    if (firstNonFinite != kNone) {
      const std::size_t row = firstNonFinite / 3;
      throw Error("the displacements hold a non-finite value at entry " + std::to_string(firstNonFinite) +
                  " (realization " + std::to_string(row / N) + ", tracer " + std::to_string(row % N) +
                  ", component " + std::to_string(firstNonFinite % 3) +
                  "), a valid displacement of a tracer the search reads");
    }

    const std::size_t keptTracers = N - ts.outsideMask - ts.outsideCut;
    if (keptTracers == 0)
      throw Error("no tracer is left to read: of " + std::to_string(N) + " tracers, " +
                  left_out(ts.outsideMask, ts.outsideCut));
    if (nRead == 0)
      throw Error("no valid displacement is left to read: of " + std::to_string(nValid) +
                  " valid displacements, none belongs to a tracer on an observed pixel inside the redshift cut (" +
                  left_out(ts.outsideMask, ts.outsideCut) + ")");

    const Selection& rs = entry.randomSelection;
    if (M - rs.outsideMask - rs.outsideCut == 0)
      throw Error("no random is left: of " + std::to_string(M) + " randoms, " +
                  left_out(rs.outsideMask, rs.outsideCut));

    keep_randoms(entry, randomsFolded, randoms);

    if (keptTracers > std::numeric_limits<unsigned>::max())
      throw Error("the search reads " + std::to_string(keptTracers) + " tracers, more than " +
                  std::to_string(std::numeric_limits<unsigned>::max()) +
                  ", the most that the density contrast within 2 mean particle separations of a centre supports");
    if (entry.nRandoms > std::numeric_limits<unsigned>::max())
      throw Error("the search reads " + std::to_string(entry.nRandoms) + " randoms, more than " +
                  std::to_string(std::numeric_limits<unsigned>::max()) +
                  ", the most that the density contrast within 2 mean particle separations of a centre supports");

    if (tracers != nullptr) entry.tracers = *tracers;
    else {
      std::vector<double> kept;
      kept.reserve(3 * keptTracers);
      for (std::size_t i = 0; i < N; ++i)
        if (ts.inMask[i] && ts.inCut[i]) kept.insert(kept.end(), entry.tracersSky.begin() + 3*i, entry.tracersSky.begin() + 3*i + 3);
      const std::vector<double> converted = bitvf::internal::call_otswap([&] { return otswap::toCartesian(kept, distances); });
      entry.tracers.assign(3 * N, std::numeric_limits<double>::quiet_NaN());
      std::size_t j = 0;
      for (std::size_t i = 0; i < N; ++i)
        if (ts.inMask[i] && ts.inCut[i]) {
          for (std::size_t k = 0; k < 3; ++k) entry.tracers[3*i+k] = converted[3*j+k];
          ++j;
        }
    }
    if (randoms == nullptr)
      entry.randoms = bitvf::internal::call_otswap([&] { return otswap::toCartesian(entry.randomsSky, distances); });

    entry.maxDisplacement = maxDisplacement;
    entry.counts.tracers = N;
    entry.counts.tracersOutsideMask = ts.outsideMask;
    entry.counts.tracersOutsideRedshiftCut = ts.outsideCut;
    entry.counts.randoms = entry.nRandoms;
    entry.counts.randomsOutsideMask = rs.outsideMask;
    entry.counts.randomsOutsideRedshiftCut = rs.outsideCut;
    entry.counts.displacements = nValid;

    const bool boundedCut = std::isfinite(cut.min) || std::isfinite(cut.max);
    report_selection(report, "tracers", N, ts, boundedCut);
    report_selection(report, "randoms", M, rs, boundedCut);
    if (report.detailed()) {
      std::ostringstream out = bitvf::internal::Report::line();
      out << "checked " << N << (N == 1 ? " tracer, " : " tracers, ") << M << (M == 1 ? " random and " : " randoms and ")
          << nValid << (nValid == 1 ? " valid displacement" : " valid displacements");
      report.write(out.str());
      std::ostringstream profile = bitvf::internal::Report::line();
      profile << "mps profile of " << mps.redshift.size() << (mps.redshift.size() == 1 ? " entry" : " entries")
              << " over redshifts [" << mps.redshift.front() << ", " << mps.redshift.back() << "] (given)";
      report.write(profile.str());
    }
    return entry;
  }

}

bitvf::Result bitvf::internal::runLightcone (const std::vector<double>* tracers, const std::vector<double>& tracersSky,
                                             const std::vector<double>& displacement,
                                             const std::vector<std::uint8_t>& valid,
                                             const std::vector<double>* randoms, const std::vector<double>& randomsSky,
                                             const otswap::Mask& mask, const otswap::DistanceTable& distances,
                                             const otswap::MpsProfile& mps, const Config& config,
                                             const otswap::RedshiftCut& redshiftCut,
                                             const std::vector<std::int64_t>& ids, const DivergenceGrid* grid)
{
  Report report(config.verbosity, config.log);
  const Entry entry = check_entry(tracers, tracersSky, displacement, valid, randoms, randomsSky, mask, distances,
                                  mps, config, redshiftCut, ids, grid, report);
  throw Error("findVoidsLightcone: the search after the entry is not implemented yet");
}

bitvf::Result bitvf::findVoidsLightcone (const std::vector<double>& tracersSky,
                                         const std::vector<double>& displacement,
                                         const std::vector<std::uint8_t>& valid,
                                         const std::vector<double>& randomsSky, const otswap::Mask& mask,
                                         const otswap::DistanceTable& distances, const otswap::MpsProfile& mps,
                                         const Config& config, const otswap::RedshiftCut& redshiftCut,
                                         const std::vector<std::int64_t>& ids, const DivergenceGrid* grid)
{
  return internal::runLightcone(nullptr, tracersSky, displacement, valid, nullptr, randomsSky, mask, distances,
                                mps, config, redshiftCut, ids, grid);
}

bitvf::Result bitvf::findVoidsLightcone (const std::vector<double>& tracers, const std::vector<double>& tracersSky,
                                         const std::vector<double>& displacement,
                                         const std::vector<std::uint8_t>& valid, const std::vector<double>& randoms,
                                         const std::vector<double>& randomsSky, const otswap::Mask& mask,
                                         const otswap::DistanceTable& distances, const otswap::MpsProfile& mps,
                                         const Config& config, const otswap::RedshiftCut& redshiftCut,
                                         const std::vector<std::int64_t>& ids, const DivergenceGrid* grid)
{
  return internal::runLightcone(&tracers, tracersSky, displacement, valid, &randoms, randomsSky, mask, distances,
                                mps, config, redshiftCut, ids, grid);
}
