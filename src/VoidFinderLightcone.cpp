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
#include <cstring>
#include <cstdint>
#include <exception>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "otswap/Maths.h"

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
    std::size_t nRandomsGiven = 0;
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

  void check_sky_rows (const std::vector<double>& sky, const std::string& plural,
                       const std::string& fromResult)
  {
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      for (std::size_t k = 0; k < 3; ++k)
        if (!std::isfinite(sky[3*i+k]))
          throw bitvf::Error(fromResult.empty()
                             ? "the sky positions of the " + plural + " hold a non-finite value for object " +
                               std::to_string(i) + " (from 0), component " + std::to_string(k)
                             : fromResult + ": the reconstruction's tracersSky holds a non-finite value for object " +
                               std::to_string(i) + " (from 0), component " + std::to_string(k));
      if (!(sky[3*i+1] >= -kPi/2 && sky[3*i+1] <= kPi/2))
        throw bitvf::Error(fromResult.empty()
                           ? "the sky positions of the " + plural + " hold a declination of " + num(sky[3*i+1]) +
                             " radians for object " + std::to_string(i) + " (from 0); it must lie in [-pi/2, pi/2]"
                           : fromResult + ": the reconstruction's tracersSky holds a declination of " +
                             num(sky[3*i+1]) + " radians for object " + std::to_string(i) +
                             " (from 0); it must lie in [-pi/2, pi/2]");
    }
  }

  void check_finite_rows (const std::vector<double>& xyz, const std::string& plural)
  {
    for (std::size_t i = 0; i < xyz.size() / 3; ++i)
      for (std::size_t k = 0; k < 3; ++k)
        if (!std::isfinite(xyz[3*i+k]))
          throw bitvf::Error("the " + plural + " hold a non-finite value for object " + std::to_string(i) +
                             " (from 0), component " + std::to_string(k));
  }

  std::vector<double> sky_of (const std::vector<double>& xyz, const std::string& plural,
                              const otswap::DistanceTable& distances)
  {
    for (std::size_t i = 0; i < xyz.size() / 3; ++i)
      if (xyz[3*i] == 0. && xyz[3*i+1] == 0. && xyz[3*i+2] == 0.)
        throw bitvf::Error("the " + plural + " hold object " + std::to_string(i) +
                           " (from 0) at the origin, the observer's position, which has no sky position");
    try {
      return otswap::toSky(xyz, distances);
    }
    catch (const otswap::Error& e) {
      const std::vector<double>& table = distances.distance();
      for (std::size_t i = 0; i < xyz.size() / 3; ++i) {
        const std::vector<double> row {xyz[3*i], xyz[3*i+1], xyz[3*i+2]};
        bool outside = false;
        try {
          otswap::toSky(row, distances);
        }
        catch (const otswap::Error&) {
          outside = true;
        }
        if (outside)
          throw bitvf::Error("the " + plural + " hold a distance of " +
                             num(std::sqrt(row[0]*row[0] + row[1]*row[1] + row[2]*row[2])) + " for object " +
                             std::to_string(i) + " (from 0), outside the distance table [" + num(table.front()) +
                             ", " + num(table.back()) + "]");
      }
      std::throw_with_nested(bitvf::Error(std::string("otswap: ") + e.what()));
    }
  }

  void check_profile (const otswap::MpsProfile& profile, const std::string& fn, const bool ofGrid)
  {
    const std::vector<double>& z = profile.redshift;
    const std::vector<double>& m = profile.mps;
    const std::string subject = ofGrid ? fn + ": the mps profile of the grid to reuse" : "the mps profile";
    if (z.empty() && m.empty())
      throw bitvf::Error(subject + " is empty; it must hold at least one entry");
    if (z.size() != m.size())
      throw bitvf::Error(fn + ": the redshift and mps lists of " +
                         (ofGrid ? "the profile of the grid to reuse" : "the profile") + " hold " +
                         std::to_string(z.size()) + " and " + std::to_string(m.size()) + " entries");
    for (std::size_t i = 0; i < z.size(); ++i) {
      if (!std::isfinite(z[i]))
        throw bitvf::Error(subject + " holds a redshift of " + num(z[i]) + " at entry " +
                           std::to_string(i) + " (from 0); it must be finite");
      if (!std::isfinite(m[i]))
        throw bitvf::Error(subject + " holds a mean particle separation of " + num(m[i]) + " at entry " +
                           std::to_string(i) + " (from 0); it must be finite");
    }
    for (std::size_t i = 1; i < z.size(); ++i)
      if (!(z[i] > z[i-1]))
        throw bitvf::Error(subject + " holds a redshift of " + num(z[i]) + " at entry " + std::to_string(i) +
                           " (from 0), not above the " + num(z[i-1]) + " before it; the redshifts must increase strictly");
    for (std::size_t i = 0; i < m.size(); ++i)
      if (!(m[i] > 0.))
        throw bitvf::Error(subject + " holds a mean particle separation of " + num(m[i]) + " at entry " +
                           std::to_string(i) + " (from 0); it must be positive");
  }

  void check_grid (const bitvf::DivergenceGrid& grid, const std::string& fn)
  {
    using bitvf::Error;
    const std::string bad = fn + ": the grid to reuse ";
    if (grid.geometry != otswap::Geometry::Lightcone) throw Error(bad + "is not a lightcone grid");
    if (!std::isfinite(grid.distanceMin) || !std::isfinite(grid.distanceMax) ||
        !(grid.distanceMin >= 0. && grid.distanceMin <= grid.distanceMax))
      throw Error(bad + "has a distance range [" + num(grid.distanceMin) + ", " + num(grid.distanceMax) +
                  "] in coordinate units; both bounds must be finite, with 0 <= min <= max");
    if (!std::isfinite(grid.raMin) || !std::isfinite(grid.raMax) ||
        !(grid.raMin <= grid.raMax && grid.raMax - grid.raMin <= kTwoPi))
      throw Error(bad + "has an RA range [" + num(grid.raMin) + ", " + num(grid.raMax) +
                  "] radians; both bounds must be finite, with min <= max and a width of at most 2 pi");
    if (!(grid.decMin >= -kPi/2 && grid.decMin <= kPi/2 && grid.decMax >= -kPi/2 && grid.decMax <= kPi/2 &&
          grid.decMin <= grid.decMax))
      throw Error(bad + "has a declination range [" + num(grid.decMin) + ", " + num(grid.decMax) +
                  "] radians; it must lie in [-pi/2, pi/2], with min <= max");
    if (!std::isfinite(grid.cellSizeMps) || grid.cellSizeMps <= 0.)
      throw Error(bad + "has a cell size of " + num(grid.cellSizeMps) +
                  " mps (cellSizeMps); it must be finite and positive");
    check_profile(grid.mpsProfile, fn, true);
    std::string boxField;
    if (!std::isnan(grid.cellSize)) boxField = "cellSize";
    else if (!std::isnan(grid.lower[0]) || !std::isnan(grid.lower[1]) || !std::isnan(grid.lower[2]))
      boxField = "lower";
    else if (!std::isnan(grid.upper[0]) || !std::isnan(grid.upper[1]) || !std::isnan(grid.upper[2]))
      boxField = "upper";
    else if (!std::isnan(grid.mps)) boxField = "mps";
    else if (grid.nCellsAxis[0] != 0 || grid.nCellsAxis[1] != 0 || grid.nCellsAxis[2] != 0)
      boxField = "nCellsAxis";
    if (!boxField.empty()) throw Error(bad + "is a lightcone grid but holds the box field " + boxField);
  }

  void compare_profiles (const otswap::MpsProfile& given, const otswap::MpsProfile& grid, const std::string& fn)
  {
    using bitvf::Error;
    if (given.redshift.size() != grid.redshift.size())
      throw Error(fn + ": the mps profile given holds " + std::to_string(given.redshift.size()) +
                  " entries and the grid's " + std::to_string(grid.redshift.size()) + "; pass the grid's profile");
    for (std::size_t i = 0; i < given.redshift.size(); ++i)
      for (const bool redshift : {true, false}) {
        const double a = redshift ? given.redshift[i] : given.mps[i];
        const double b = redshift ? grid.redshift[i] : grid.mps[i];
        if (std::memcmp(&a, &b, sizeof a) != 0)
          throw Error(fn + ": the mps profile given differs from the grid's at entry " + std::to_string(i) +
                      " (from 0): " + (redshift ? "redshift " : "mps ") + num(a) + " against " + num(b) +
                      "; pass the grid's profile");
      }
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

  void check_usable (const std::vector<double>& sky, const Selection& selection,
                     const otswap::DistanceTable& distances, const std::string& plural)
  {
    const double lo = distances.minRedshift(), hi = distances.maxRedshift();
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      if (!selection.inMask[i] || !selection.inCut[i]) continue;
      const double z = sky[3*i+2];
      if (!(z >= lo && z <= hi))
        throw bitvf::Error("the sky positions of the " + plural + " hold a redshift of " + num(z) + " for object " +
                           std::to_string(i) + " (from 0), outside the distance table [" + num(lo) + ", " + num(hi) + "]");
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

  otswap::RedshiftCut check_reconstruction (const otswap::Result& rec, const otswap::Mask& mask,
                                             const bitvf::DivergenceGrid* grid, const bool builder)
  {
    using bitvf::Error;
    const std::string fn = builder ? "lightconeGrid" : "findVoidsLightcone";

    if (rec.geometry != otswap::Geometry::Lightcone)
      throw Error(fn + ": the reconstruction is not a lightcone result; use " + (builder ? "boxGrid" : "findVoidsBox"));

    std::vector<std::string> missing;
    if (rec.tracersSky.empty()) missing.push_back("tracersSky");
    if (!rec.distances) missing.push_back("distances");
    if (grid == nullptr && rec.mpsProfile.redshift.empty() && rec.mpsProfile.mps.empty())
      missing.push_back("mpsProfile");
    if (!missing.empty()) {
      std::string list = missing[0];
      for (std::size_t i = 1; i < missing.size(); ++i)
        list += (i + 1 == missing.size() ? " and " : ", ") + missing[i];
      throw Error(fn + ": the reconstruction has no " + list + ", which a lightcone search needs" +
                  (!builder && missing.size() == 1 && missing[0] == "mpsProfile"
                   ? "; give a grid to search without mpsProfile" : ""));
    }

    const std::size_t N = rec.nObjects;
    const std::size_t R = rec.nRealizations;
    const std::string malformed = fn + ": the reconstruction is malformed: ";
    bitvf::internal::check_reconstruction_shape(rec, fn);
    auto check_size = [&] (const std::size_t size, const std::size_t expected,
                           const std::string& name, const bool allowEmpty) {
      if (size != expected && !(allowEmpty && size == 0))
        throw Error(malformed + name + " holds " + std::to_string(size) + " entries for " +
                    std::to_string(R) + " realizations of " + std::to_string(N) + " objects");
    };
    check_size(rec.tracersSky.size(), 3 * N, "tracersSky", false);
    check_size(rec.outsideMask.size(), N, "outsideMask", true);
    check_size(rec.outsideRedshiftCut.size(), N, "outsideRedshiftCut", true);

    if (rec.filteredNside != 0 && rec.filteredNside != mask.nside())
      throw Error(fn + ": the reconstruction was filtered against a mask of NSIDE " +
                  std::to_string(rec.filteredNside) +
                  " (otswap::Result::filteredNside), and the mask given has NSIDE " + std::to_string(mask.nside()) +
                  "; the two must agree");

    otswap::RedshiftCut cut;
    if (rec.selection.redshiftCut) {
      cut = *rec.selection.redshiftCut;
      if (std::isnan(cut.min) || std::isnan(cut.max) || cut.min > cut.max)
        throw Error(malformed + "its redshift cut (otswap::Result::selection) is [" + num(cut.min) + ", " +
                    num(cut.max) + "]");
    }
    return cut;
  }

  Entry check_entry (const std::vector<double>* tracers, const std::vector<double>& tracersSky,
                     const std::vector<double>& displacement, const std::vector<std::uint8_t>& valid,
                     const std::vector<double>* randoms, const std::vector<double>& randomsSky,
                     const otswap::Mask& mask, const otswap::DistanceTable& distances,
                     const otswap::MpsProfile& mps, const bitvf::Config& config,
                     const otswap::RedshiftCut& cut, const std::vector<std::int64_t>& ids,
                     const bitvf::DivergenceGrid* grid, const otswap::Result* rec, const std::string& fn,
                     const bool withDisplacement, const bitvf::internal::Report& report)
  {
    using bitvf::Error;

    const bool cartesian = randoms != nullptr;
    std::vector<double> tracersConverted, randomsConverted;
    if (cartesian) {
      check_sky(*tracers, "tracers", fn);
      check_sky(*randoms, "randoms", fn);
      check_finite_rows(*tracers, "tracers");
      check_finite_rows(*randoms, "randoms");
      tracersConverted = sky_of(*tracers, "tracers", distances);
      randomsConverted = sky_of(*randoms, "randoms", distances);
    }
    const std::vector<double>& tSky = cartesian ? tracersConverted : tracersSky;
    const std::vector<double>& rSky = cartesian ? randomsConverted : randomsSky;

    const std::size_t N = check_sky(tSky, "tracersSky", fn);
    const std::size_t M = check_sky(rSky, "randomsSky", fn);

    const std::size_t R = withDisplacement ? bitvf::internal::check_displacement_shape(displacement, valid, N, fn) : 0;

    if (std::isnan(cut.min) || std::isnan(cut.max) || cut.min > cut.max)
      throw Error("the redshift cut is [" + num(cut.min) + ", " + num(cut.max) +
                  "]; neither bound may be NaN and the lower bound must not exceed the upper one");

    const bool readProfile = rec == nullptr || grid == nullptr;
    if (readProfile) check_profile(mps, fn, false);
    if (rec == nullptr) bitvf::internal::check_ids(ids, N, fn, false);
    bitvf::internal::check_config(config, grid != nullptr);
    if (grid != nullptr) {
      check_grid(*grid, fn);
      if (rec == nullptr) compare_profiles(mps, grid->mpsProfile, fn);
    }

    check_sky_rows(tSky, "tracers", rec == nullptr ? "" : fn);
    check_sky_rows(rSky, "randoms", "");

    Entry entry;
    entry.nTracers = N;
    entry.nRealizations = (unsigned)R;
    entry.nRandomsGiven = M;
    entry.tracersSky = tSky;
    fold_right_ascension(entry.tracersSky);
    std::vector<double> randomsFolded = rSky;
    fold_right_ascension(randomsFolded);

    if (rec == nullptr) entry.tracerSelection = select_rows(entry.tracersSky, mask, cut);
    else {
      Selection& selection = entry.tracerSelection;
      selection = select_rows(entry.tracersSky, mask, otswap::RedshiftCut{});
      for (std::size_t i = 0; i < rec->outsideMask.size(); ++i)
        if (rec->outsideMask[i] && selection.inMask[i])
          throw Error(fn + ": the mask given allows object " + std::to_string(i) +
                      " (from 0) of tracers, which the reconstruction excluded as outside its mask (Result::outsideMask); "
                      "the mask must not be more permissive than the reconstruction's");
      for (std::size_t i = 0; i < rec->outsideRedshiftCut.size(); ++i)
        if (rec->outsideRedshiftCut[i]) {
          selection.inCut[i] = 0;
          if (selection.inMask[i]) ++selection.outsideCut;
        }
    }
    entry.randomSelection = select_rows(randomsFolded, mask, cut);

    const Selection& ts = entry.tracerSelection;
    check_usable(entry.tracersSky, ts, distances, "tracers");
    if (rec != nullptr)
      for (std::size_t i = 0; i < N; ++i) {
        if (!ts.inMask[i] || !ts.inCut[i]) continue;
        for (std::size_t k = 0; k < 3; ++k)
          if (!std::isfinite((*tracers)[3*i+k]))
            throw Error(fn + ": the reconstruction's tracers hold a non-finite value for object " +
                        std::to_string(i) + " (from 0), component " + std::to_string(k));
      }
    check_usable(randomsFolded, entry.randomSelection, distances, "randoms");

    const std::int64_t rows = (std::int64_t)(R * N);
    std::size_t nValid = 0, nRead = 0, firstNonFinite = kNone;
#pragma omp parallel for schedule(static) reduction(+:nValid, nRead) reduction(min:firstNonFinite)
    for (std::int64_t row = 0; row < rows; ++row) {
      const std::size_t r = (std::size_t)row;
      const std::size_t i = r % N;
      if (!valid.empty() && !valid[r]) continue;
      ++nValid;
      if (!ts.inMask[i] || !ts.inCut[i]) continue;
      if (grid != nullptr) {
        ++nRead;
        continue;
      }
      bool finite = true;
      for (std::size_t k = 0; k < 3; ++k)
        if (!std::isfinite(displacement[3*r+k])) {
          firstNonFinite = std::min(firstNonFinite, 3*r+k);
          finite = false;
        }
      if (finite) ++nRead;
    }
    if (firstNonFinite != kNone) {
      const std::size_t row = firstNonFinite / 3;
      throw Error("the displacements hold a non-finite value for realization " + std::to_string(row / N) +
                  " (from 0), object " + std::to_string(row % N) + " (from 0), component " +
                  std::to_string(firstNonFinite % 3) + ", a valid displacement of a tracer the search reads");
    }

    const std::size_t keptTracers = N - ts.outsideMask - ts.outsideCut;
    if (keptTracers == 0)
      throw Error("no tracer is left to read: of " + std::to_string(N) + " tracers, " +
                  left_out(ts.outsideMask, ts.outsideCut));
    if (withDisplacement && nRead == 0)
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
    }
    if (report.detailed() && readProfile) {
      std::ostringstream profile = bitvf::internal::Report::line();
      profile << "mps profile of " << mps.redshift.size() << (mps.redshift.size() == 1 ? " entry" : " entries")
              << " over redshifts [" << mps.redshift.front() << ", " << mps.redshift.back()
              << (rec == nullptr ? "] (given)" : "] (from the reconstruction)");
      report.write(profile.str());
    }
    return entry;
  }

  bitvf::internal::LightconeGrid entry_grid (Entry& entry, const otswap::Mask& mask,
                                             const otswap::DistanceTable& distances, const otswap::MpsProfile& mps,
                                             const std::optional<double>& divergenceCellSize,
                                             const bitvf::DivergenceGrid* grid, const std::string& fn,
                                             bitvf::internal::Report& report)
  {
    using bitvf::Error;
    using bitvf::internal::LightconeGrid;
    using bitvf::internal::Report;

    bitvf::DivergenceGrid record;
    if (grid == nullptr) {
      const bitvf::internal::LightconeExtent extent =
        bitvf::internal::lightcone_extent(entry.tracersSky, entry.tracerSelection.inMask, entry.tracerSelection.inCut,
                                          entry.randomsSky, distances);
      for (const double z : {extent.redshiftMin, extent.redshiftMax}) {
        try {
          mps.at(z);
        }
        catch (const otswap::Error&) {
          throw Error("the mps profile is not positive at redshift " + num(z) + " (tracers and randoms span [" +
                      num(extent.redshiftMin) + ", " + num(extent.redshiftMax) + "])");
        }
      }
      record.geometry = otswap::Geometry::Lightcone;
      record.distanceMin = extent.distanceMin;
      record.distanceMax = extent.distanceMax;
      record.raMin = extent.raMin;
      record.raMax = extent.raMax;
      record.decMin = extent.decMin;
      record.decMax = extent.decMax;
      record.cellSizeMps = divergenceCellSize ? *divergenceCellSize
                                              : std::min(1., otswap::maths::cbrt(20. / entry.nRealizations));
      record.mpsProfile = mps;
    }
    else {
      record.geometry = otswap::Geometry::Lightcone;
      record.distanceMin = grid->distanceMin;
      record.distanceMax = grid->distanceMax;
      record.raMin = grid->raMin;
      record.raMax = grid->raMax;
      record.decMin = grid->decMin;
      record.decMax = grid->decMax;
      record.cellSizeMps = grid->cellSizeMps;
      record.mpsProfile = grid->mpsProfile;
    }
    const LightconeGrid lc = bitvf::internal::lightcone_grid(record, distances, mask, grid != nullptr);
    if (grid != nullptr && grid->nCells != lc.cells.nCells)
      throw Error(fn + ": the grid to reuse holds " + std::to_string(grid->nCells) + " cells but the rebuild gives " +
                  std::to_string(lc.cells.nCells) + ": the record was modified or the distance table differs");

    const std::int64_t nRandoms = (std::int64_t)entry.nRandoms;
    std::vector<std::uint8_t> inGrid((std::size_t)nRandoms, 1);
    std::size_t outsideGrid = 0;
    std::exception_ptr failure;
#pragma omp parallel for schedule(static) reduction(+:outsideGrid)
    for (std::int64_t i = 0; i < nRandoms; ++i) {
      const std::size_t r = (std::size_t)i;
      try {
        const double d = distances.distanceAt(entry.randomsSky[3*r+2]);
        if (bitvf::internal::lightcone_cell(lc, d, entry.randomsSky[3*r], entry.randomsSky[3*r+1]) < 0) {
          inGrid[r] = 0;
          ++outsideGrid;
        }
      }
      catch (...) {
#pragma omp critical(bitvf_lightcone_randoms_grid)
        if (!failure) failure = std::current_exception();
      }
    }
    if (failure) {
      try {
        std::rethrow_exception(failure);
      }
      catch (const otswap::Error& e) {
        std::throw_with_nested(Error(std::string("otswap: ") + e.what()));
      }
    }
    if (outsideGrid > 0) {
      std::size_t j = 0;
      for (std::size_t i = 0; i < entry.nRandoms; ++i) {
        if (!inGrid[i]) continue;
        for (std::size_t k = 0; k < 3; ++k) {
          entry.randomsSky[3*j+k] = entry.randomsSky[3*i+k];
          entry.randoms[3*j+k] = entry.randoms[3*i+k];
        }
        ++j;
      }
      entry.nRandoms = j;
      entry.randomsSky.resize(3 * j);
      entry.randoms.resize(3 * j);
    }
    entry.counts.randoms = entry.nRandoms;
    entry.counts.randomsOutsideGrid = outsideGrid;
    const Selection& rs = entry.randomSelection;
    if (entry.nRandoms == 0)
      throw Error("no random is left: of " + std::to_string(entry.nRandomsGiven) + " randoms, " +
                  left_out(rs.outsideMask, rs.outsideCut) + ", " + std::to_string(outsideGrid) +
                  (grid == nullptr ? " outside the grid" : " outside the reused grid"));
    if (entry.nRandoms > std::numeric_limits<unsigned>::max())
      throw Error("the search reads " + std::to_string(entry.nRandoms) + " randoms, more than " +
                  std::to_string(std::numeric_limits<unsigned>::max()) +
                  ", the most that the density contrast within 2 mean particle separations of a centre supports");

    if (report.normal()) {
      constexpr double kDegrees = 180. / kPi;
      const std::size_t nShells = lc.shellCentre.size();
      std::ostringstream out = Report::line();
      out << "grid of " << lc.cells.nCells << (lc.cells.nCells == 1 ? " cell in " : " cells in ") << nShells
          << (nShells == 1 ? " shell" : " shells") << (lc.centralSphere ? ", central sphere included" : "")
          << ", cell size " << record.cellSizeMps << " mps ("
          << (grid != nullptr ? "reused" : divergenceCellSize ? "given" : "default") << "), distances ["
          << record.distanceMin << ", " << record.distanceMax << "], RA [" << record.raMin * kDegrees << ", "
          << record.raMax * kDegrees << "] deg, Dec [" << record.decMin * kDegrees << ", " << record.decMax * kDegrees
          << (grid != nullptr ? "] deg (from the reused grid)" : "] deg (from the tracers and randoms)");
      report.step(out);
    }
    if (report.detailed()) {
      std::size_t masked = 0;
      for (const std::uint8_t m : lc.cells.masked) masked += m;
      std::ostringstream out = Report::line();
      out << "grid: " << masked << " of " << lc.cells.nCells << " cells masked, redshifts [" << lc.redshiftFirstEdge
          << ", " << lc.redshiftLastEdge << "] from the first to the last edge, smallest mps " << lc.mpsMin;
      report.write(out.str());
    }
    return lc;
  }

  bitvf::DivergenceGrid defining_fields (const bitvf::DivergenceGrid& cells)
  {
    bitvf::DivergenceGrid grid;
    grid.geometry = otswap::Geometry::Lightcone;
    grid.nCells = cells.nCells;
    grid.distanceMin = cells.distanceMin;
    grid.distanceMax = cells.distanceMax;
    grid.raMin = cells.raMin;
    grid.raMax = cells.raMax;
    grid.decMin = cells.decMin;
    grid.decMax = cells.decMax;
    grid.cellSizeMps = cells.cellSizeMps;
    grid.mpsProfile = cells.mpsProfile;
    return grid;
  }

}

bitvf::Result bitvf::internal::runLightcone (const std::vector<double>* tracers, const std::vector<double>& tracersSky,
                                             const std::vector<double>& displacement,
                                             const std::vector<std::uint8_t>& valid,
                                             const std::vector<double>* randoms, const std::vector<double>& randomsSky,
                                             const otswap::Mask& mask, const otswap::DistanceTable& distances,
                                             const otswap::MpsProfile& mps, const Config& config,
                                             const otswap::RedshiftCut& redshiftCut,
                                             const std::vector<std::int64_t>& ids, const DivergenceGrid* grid,
                                             const otswap::Result* reconstruction)
{
  const std::string fn = "findVoidsLightcone";
  Report report(config.verbosity, config.log);
  Entry entry = check_entry(tracers, tracersSky, displacement, valid, randoms, randomsSky, mask, distances,
                            mps, config, redshiftCut, ids, grid, reconstruction, fn, true, report);
  LightconeGrid lc = entry_grid(entry, mask, distances, mps, config.divergenceCellSize, grid, fn, report);
  const std::vector<double> sigmaShell = lightcone_sigma(config.smoothing, lc);

  LightconeSweep sweep;
  sweep.nTracers = entry.nTracers;
  sweep.nRealizations = entry.nRealizations;
  sweep.tracers = &entry.tracers;
  sweep.tracersSky = &entry.tracersSky;
  sweep.inMask = &entry.tracerSelection.inMask;
  sweep.inCut = &entry.tracerSelection.inCut;
  sweep.displacement = &displacement;
  sweep.valid = &valid;
  sweep.distances = &distances;
  sweep.reusedGrid = grid != nullptr;
  const LightconeDivergence div = lightcone_divergence(sweep, lc);
  entry.counts.tracersOutsideGrid = div.nOutsideGrid;
  entry.counts.tracersUsed = div.nUsed;
  entry.counts.displacementsRejected = div.nRejected;

  if (grid != nullptr && report.normal()) {
    const std::size_t kept = entry.nTracers - entry.counts.tracersOutsideMask - entry.counts.tracersOutsideRedshiftCut;
    std::ostringstream out = Report::line();
    out << "kept " << kept - div.nOutsideGrid << " of " << kept << (kept == 1 ? " tracer: " : " tracers: ")
        << div.nOutsideGrid << " outside the reused grid";
    report.write(out.str());
  }
  if (report.detailed() || (report.normal() && div.nRejected > 0)) {
    std::ostringstream out = Report::line();
    out << "divergence of the valid displacements of " << div.nUsed << (div.nUsed == 1 ? " tracer, " : " tracers, ")
        << div.nRejected << (div.nRejected == 1 ? " displacement" : " displacements")
        << " rejected (masked cells, outside or below the grid)";
    report.step(out);
  }

  const LightconeAdjacency adjacency = lightcone_adjacency(lc);
  const LightconeFill fill = lightcone_fill(adjacency, lc.cells);
  if (report.detailed() || (report.normal() && fill.leftEmpty > 0)) {
    std::ostringstream out = Report::line();
    out << "filled " << fill.filled << (fill.filled == 1 ? " cell in " : " cells in ") << fill.passes
        << (fill.passes == 1 ? " pass, " : " passes, ") << fill.leftEmpty << " left empty";
    report.step(out);
  }

  double wMax = 0.;
  if (!sigmaShell.empty()) {
    const std::int64_t C = (std::int64_t)lc.cells.nCells;
#pragma omp parallel for schedule(static) reduction(max:wMax)
    for (std::int64_t c = 0; c < C; ++c) wMax = std::max(wMax, lc.cells.weight[(std::size_t)c]);
  }
  lightcone_smooth(lc, sigmaShell, 1e-12 * wMax);
  if (report.detailed()) {
    if (!sigmaShell.empty()) {
      std::ostringstream out = Report::line();
      if (config.smoothing.unit == LengthUnit::Mps)
        out << "smoothed with sigma " << config.smoothing.sigma << " mps ("
            << *std::min_element(sigmaShell.begin(), sigmaShell.end()) << " to "
            << *std::max_element(sigmaShell.begin(), sigmaShell.end()) << " in coordinate units)";
      else
        out << "smoothed with sigma " << sigmaShell.front() << " in coordinate units";
      report.step(out);
    }
    else {
      report.lap();
      report.write("not smoothed (sigma 0)");
    }
  }

  LightconeBasins basins = lightcone_watershed(adjacency, lc.cells);
  if (report.detailed()) {
    std::ostringstream out = Report::line();
    out << "watershed into ";
    Report::basins(out, basins.voids.seed.size(), basins.overdense.seed.size());
    report.step(out);
  }
  if (basins.voids.seed.empty())
    throw Error("no void found: the divergence is 0 or positive in every unmasked cell");
  BasinCatalog voids, overdense;
  lightcone_basins(lc, adjacency, basins, voids, overdense);

  throw Error("findVoidsLightcone: the search after the divergence is not implemented yet");
}


// ============================================================================


bitvf::DivergenceGrid bitvf::lightconeGrid (const otswap::Result& reconstruction, const std::vector<double>& randomsSky,
                                            const otswap::Mask& mask, const std::optional<double> divergenceCellSize)
{
  const std::string fn = "lightconeGrid";
  const otswap::RedshiftCut cut = check_reconstruction(reconstruction, mask, nullptr, true);
  const otswap::Result& rec = reconstruction;
  Config config;
  config.verbosity = otswap::Verbosity::Silent;
  config.divergenceCellSize = divergenceCellSize;
  internal::Report report(config.verbosity, nullptr);
  Entry entry = check_entry(&rec.tracers, rec.tracersSky, rec.displacement, rec.valid, nullptr, randomsSky, mask,
                            *rec.distances, rec.mpsProfile, config, cut, rec.ids, nullptr, &rec, fn, true, report);
  const internal::LightconeGrid lc = entry_grid(entry, mask, *rec.distances, rec.mpsProfile,
                                                divergenceCellSize, nullptr, fn, report);
  return defining_fields(lc.cells);
}


// ============================================================================


bitvf::DivergenceGrid bitvf::lightconeGrid (const std::vector<double>& tracersSky, const std::vector<double>& randomsSky,
                                            const otswap::Mask& mask, const otswap::DistanceTable& distances,
                                            const otswap::MpsProfile& mps, const double divergenceCellSize,
                                            const otswap::RedshiftCut& redshiftCut)
{
  const std::string fn = "lightconeGrid";
  Config config;
  config.verbosity = otswap::Verbosity::Silent;
  config.divergenceCellSize = divergenceCellSize;
  internal::Report report(config.verbosity, nullptr);
  Entry entry = check_entry(nullptr, tracersSky, {}, {}, nullptr, randomsSky, mask, distances, mps, config,
                            redshiftCut, {}, nullptr, nullptr, fn, false, report);
  const internal::LightconeGrid lc = entry_grid(entry, mask, distances, mps,
                                                config.divergenceCellSize, nullptr, fn, report);
  return defining_fields(lc.cells);
}


// ============================================================================


bitvf::Result bitvf::findVoidsLightcone (const std::vector<double>& tracersSky,
                                         const std::vector<double>& displacement,
                                         const std::vector<std::uint8_t>& valid,
                                         const std::vector<double>& randomsSky, const otswap::Mask& mask,
                                         const otswap::DistanceTable& distances, const otswap::MpsProfile& mps,
                                         const Config& config, const otswap::RedshiftCut& redshiftCut,
                                         const std::vector<std::int64_t>& ids, const DivergenceGrid* grid)
{
  return internal::runLightcone(nullptr, tracersSky, displacement, valid, nullptr, randomsSky, mask, distances,
                                mps, config, redshiftCut, ids, grid, nullptr);
}

bitvf::Result bitvf::findVoidsLightconeFromCartesian (const std::vector<double>& tracers,
                                                      const std::vector<double>& displacement,
                                                      const std::vector<std::uint8_t>& valid,
                                                      const std::vector<double>& randoms, const otswap::Mask& mask,
                                                      const otswap::DistanceTable& distances,
                                                      const otswap::MpsProfile& mps, const Config& config,
                                                      const otswap::RedshiftCut& redshiftCut,
                                                      const std::vector<std::int64_t>& ids, const DivergenceGrid* grid)
{
  return internal::runLightcone(&tracers, {}, displacement, valid, &randoms, {}, mask, distances, mps, config,
                                redshiftCut, ids, grid, nullptr);
}

bitvf::Result bitvf::findVoidsLightcone (const otswap::Result& reconstruction, const std::vector<double>& randomsSky,
                                         const otswap::Mask& mask, const Config& config, const DivergenceGrid* grid)
{
  const otswap::RedshiftCut cut = check_reconstruction(reconstruction, mask, grid, false);
  const otswap::Result& rec = reconstruction;
  return internal::runLightcone(&rec.tracers, rec.tracersSky, rec.displacement, rec.valid, nullptr, randomsSky, mask,
                                *rec.distances, rec.mpsProfile, config, cut, rec.ids, grid, &rec);
}
