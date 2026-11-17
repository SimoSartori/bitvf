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
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;
  constexpr double kHalfPi = kPi / 2.;
  constexpr std::int64_t kMaxBins = std::int64_t(1) << 20;
  constexpr std::uint64_t kEmpty = std::numeric_limits<std::uint64_t>::max();

  std::string num (const double x)
  {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out.precision(17);
    out << x;
    return out.str();
  }

  double fold (const double ra)
  {
    double folded = std::fmod(ra, kTwoPi);
    if (folded < 0.) folded += kTwoPi;
    if (folded >= kTwoPi) folded = 0.;
    return folded + 0.;
  }

  std::uint64_t bits (const double x)
  {
    std::uint64_t b;
    std::memcpy(&b, &x, sizeof b);
    return b;
  }

  double value (const std::uint64_t b)
  {
    double x;
    std::memcpy(&x, &b, sizeof x);
    return x;
  }

  void atomic_min (std::atomic<std::uint64_t>& a, const std::uint64_t v)
  {
    std::uint64_t current = a.load(std::memory_order_relaxed);
    while (v < current && !a.compare_exchange_weak(current, v, std::memory_order_relaxed)) {}
  }

  void atomic_max (std::atomic<std::uint64_t>& a, const std::uint64_t v)
  {
    std::uint64_t current = a.load(std::memory_order_relaxed);
    while (v > current && !a.compare_exchange_weak(current, v, std::memory_order_relaxed)) {}
  }

  std::string limit_message (const std::string& grid, const std::string& what)
  {
    return grid + " would need more than " + std::to_string(bitvf::internal::kMaxCellsAxis) + " " + what +
           "; increase the divergence cell size";
  }

}


// ============================================================================


bitvf::internal::LightconeExtent bitvf::internal::lightcone_extent (const std::vector<double>& tracersSky,
                                                                    const std::vector<std::uint8_t>& inMask,
                                                                    const std::vector<std::uint8_t>& inCut,
                                                                    const std::vector<double>& randomsSky,
                                                                    const otswap::DistanceTable& distances)
{
  const std::int64_t nTracers = (std::int64_t)(tracersSky.size() / 3);
  const std::int64_t nRandoms = (std::int64_t)(randomsSky.size() / 3);
  auto kept = [&] (const std::int64_t i) {
    return (inMask.empty() || inMask[(std::size_t)i]) && (inCut.empty() || inCut[(std::size_t)i]);
  };

  std::int64_t n = 0;
  double zMin = std::numeric_limits<double>::infinity(), zMax = -zMin;
  double decMin = zMin, decMax = zMax;
#pragma omp parallel for schedule(static) reduction(+:n) reduction(min:zMin, decMin) reduction(max:zMax, decMax)
  for (std::int64_t i = 0; i < nTracers + nRandoms; ++i) {
    const double* p = i < nTracers ? &tracersSky[3*(std::size_t)i] : &randomsSky[3*(std::size_t)(i - nTracers)];
    if (i < nTracers && !kept(i)) continue;
    ++n;
    zMin = std::min(zMin, p[2]);
    zMax = std::max(zMax, p[2]);
    decMin = std::min(decMin, p[1]);
    decMax = std::max(decMax, p[1]);
  }
  if (n == 0) throw Error("lightcone_extent: no object to measure");

  const std::int64_t k = std::min(n, kMaxBins);
  const std::unique_ptr<std::atomic<std::uint64_t>[]> lo(new std::atomic<std::uint64_t>[(std::size_t)k]);
  const std::unique_ptr<std::atomic<std::uint64_t>[]> hi(new std::atomic<std::uint64_t>[(std::size_t)k]);
  for (std::int64_t b = 0; b < k; ++b) {
    lo[(std::size_t)b].store(kEmpty, std::memory_order_relaxed);
    hi[(std::size_t)b].store(0, std::memory_order_relaxed);
  }
  const double scale = (double)k / kTwoPi;
#pragma omp parallel for schedule(static)
  for (std::int64_t i = 0; i < nTracers + nRandoms; ++i) {
    const double* p = i < nTracers ? &tracersSky[3*(std::size_t)i] : &randomsSky[3*(std::size_t)(i - nTracers)];
    if (i < nTracers && !kept(i)) continue;
    const double ra = p[0] + 0.;
    const std::int64_t b = std::min(k - 1, (std::int64_t)(ra * scale));
    atomic_min(lo[(std::size_t)b], bits(ra));
    atomic_max(hi[(std::size_t)b], bits(ra));
  }

  std::int64_t first = -1, last = -1;
  for (std::int64_t b = 0; b < k; ++b)
    if (lo[(std::size_t)b].load(std::memory_order_relaxed) != kEmpty) {
      if (first < 0) first = b;
      last = b;
    }
  const double firstValue = value(lo[(std::size_t)first].load(std::memory_order_relaxed));
  const double lastValue = value(hi[(std::size_t)last].load(std::memory_order_relaxed));

  LightconeExtent extent;
  extent.raMin = firstValue;
  extent.raMax = lastValue;
  double largest = firstValue + kTwoPi - lastValue;
  std::int64_t previous = first;
  for (std::int64_t b = first + 1; b <= last; ++b) {
    const std::uint64_t low = lo[(std::size_t)b].load(std::memory_order_relaxed);
    if (low == kEmpty) continue;
    const double before = value(hi[(std::size_t)previous].load(std::memory_order_relaxed));
    const double after = value(low);
    if (after - before > largest) {
      largest = after - before;
      extent.raMin = after;
      extent.raMax = before + kTwoPi;
    }
    previous = b;
  }

  extent.redshiftMin = zMin;
  extent.redshiftMax = zMax;
  extent.decMin = decMin;
  extent.decMax = decMax;
  extent.distanceMin = call_otswap([&] { return distances.distanceAt(zMin); });
  extent.distanceMax = call_otswap([&] { return distances.distanceAt(zMax); });
  return extent;
}


// ============================================================================


bitvf::internal::LightconeGrid bitvf::internal::lightcone_grid (const DivergenceGrid& record,
                                                                const otswap::DistanceTable& distances,
                                                                const otswap::Mask& mask, const bool reused)
{
  const std::string name = reused ? "the grid to reuse" : "the grid";
  const std::string of = reused ? " of the grid to reuse" : "";
  const otswap::MpsProfile& profile = record.mpsProfile;
  const double s = record.cellSizeMps;
  const double d0 = distances.distance().front(), d1 = distances.distance().back();

  auto redshift = [&] (const double d) {
    if (!(d >= d0 && d <= d1))
      throw Error(name + " reaches distance " + num(d) + ", outside the distance table [" + num(d0) + ", " +
                  num(d1) + "]; extend the table");
    return call_otswap([&] { return distances.redshiftAt(d); });
  };

  double rangeLo = 0., rangeHi = 0.;
  auto mps = [&] (const double z) {
    try {
      return profile.at(z);
    }
    catch (const otswap::Error&) {
      throw Error("the mps profile" + of + " is not positive at redshift " + num(z) + ", inside the grid's range [" +
                  num(std::min(rangeLo, z)) + ", " + num(std::max(rangeHi, z)) + "]");
    }
  };

  LightconeGrid g;
  const double zMin = redshift(record.distanceMin);
  const double zMax = redshift(record.distanceMax);
  rangeLo = zMin;
  rangeHi = zMax;
  const double mpsIn = mps(zMin);
  const double mpsOut = mps(zMax);
  const double marginIn = 0.1 * s * mpsIn;
  const double marginOut = 0.1 * s * mpsOut;
  const double sphere = otswap::maths::cbrt(3. / (4. * kPi)) * s * mpsIn;

  std::vector<double>& edge = g.shellEdge;
  if (record.distanceMin > sphere) edge.push_back(record.distanceMin - marginIn);
  else {
    g.centralSphere = true;
    edge.push_back(0.);
    g.shellCentre.push_back(0.);
    g.shellRedshift.push_back(redshift(0.));
    edge.push_back(sphere);
  }
  g.redshiftFirstEdge = redshift(edge.front());
  rangeLo = std::min(rangeLo, g.redshiftFirstEdge);

  double width = s * mps(redshift(edge.back()));
  const double outer = record.distanceMax + marginOut;
  while (edge.back() < outer) {
    if (g.shellCentre.size() == kMaxCellsAxis) throw Error(limit_message(name, "shells"));
    const double previous = edge.back();
    const double zStep = redshift(previous + width / 2.);
    rangeHi = std::max(rangeHi, zStep);
    const double next = previous + s * mps(zStep);
    redshift(next);
    width = next - previous;
    const double centre = (next + previous) / 2.;
    edge.push_back(next);
    g.shellCentre.push_back(centre);
    g.shellRedshift.push_back(redshift(centre));
  }
  g.redshiftLastEdge = redshift(edge.back());

  rangeLo = g.redshiftFirstEdge;
  rangeHi = g.redshiftLastEdge;
  g.mpsMin = std::min(mps(g.redshiftFirstEdge), mps(g.redshiftLastEdge));
  for (std::size_t i = 0; i < profile.redshift.size(); ++i)
    if (profile.redshift[i] > g.redshiftFirstEdge && profile.redshift[i] < g.redshiftLastEdge)
      g.mpsMin = std::min(g.mpsMin, mps(profile.redshift[i]));

  const std::size_t nShells = g.shellCentre.size();
  for (std::size_t i = 0; i < nShells; ++i) g.shellMps.push_back(mps(g.shellRedshift[i]));

  std::vector<double> ringCentre, ringRaStart, ringRaStep;
  std::vector<std::int64_t> ringSectors;
  auto add_ring = [&] (const double decLo, const double decHi, const double centre, const bool pole, const bool full,
                       const double raStart, const double raStep, const std::int64_t sectors) {
    g.ringDecMin.push_back(decLo);
    ringCentre.push_back(centre);
    g.ringDecMax.push_back(decHi);
    g.ringPole.push_back(pole ? 1 : 0);
    g.ringFull.push_back(full ? 1 : 0);
    ringRaStart.push_back(raStart);
    ringRaStep.push_back(raStep);
    ringSectors.push_back(sectors);
  };

  const double decMin = record.decMin, decMax = record.decMax;
  const double raMin = record.raMin, raMax = record.raMax;
  for (std::size_t i = 0; i < nShells; ++i) {
    g.shellFirstRing.push_back((std::int64_t)g.ringDecMin.size());
    if (i == 0 && g.centralSphere) {
      add_ring(-kHalfPi, kHalfPi, 0., false, true, 0., kTwoPi, 1);
      continue;
    }
    const double cell = s * g.shellMps[i];
    const double R = g.shellCentre[i];
    const double axisCell = cell / std::sqrt(kPi);
    const double capWidth = axisCell / R;
    double decLo = decMin - 0.1 * cell / R;
    double decHi = decMax + 0.1 * cell / R;
    const bool down = decMin < -kHalfPi + capWidth;
    const bool up = decMax > kHalfPi - capWidth;
    if (down) decLo = -kHalfPi;
    if (up) decHi = kHalfPi;
    const int caps = (down ? 1 : 0) + (up ? 1 : 0);
    const double middle = std::round(((decHi - decLo) * R - caps * axisCell) / cell);
    if (!(middle >= 0.) || middle + caps == 0.)
      throw Error("shell " + std::to_string(i) + " (from 0)" + of +
                  " holds no declination ring; use a smaller cell or a central sphere");
    if (middle + caps > (double)kMaxCellsAxis)
      throw Error(limit_message(name, "rings in shell " + std::to_string(i) + " (from 0)"));
    const std::int64_t nMiddle = (std::int64_t)middle;
    const std::int64_t nRings = nMiddle + caps;
    const double step = ((decHi - decLo) - caps * capWidth) / (double)std::max<std::int64_t>(1, nMiddle);

    double lower = std::fabs(decLo) < 1e-12 ? 0. : decLo;
    for (std::int64_t j = 0; j < nRings; ++j) {
      const bool capDown = j == 0 && down;
      const bool capUp = j == nRings - 1 && up;
      double upper = j == nRings - 1 ? decHi : lower + (capDown || capUp ? capWidth : step);
      if (std::fabs(upper) < 1e-12) upper = 0.;
      double centre = (upper + lower) / 2.;
      if (capDown) centre = -kHalfPi;
      if (capUp) centre = kHalfPi;

      if (std::fabs(centre) >= kHalfPi - 1e-5)
        add_ring(lower, upper, centre, true, true, 0., kTwoPi, 1);
      else {
        const double RcosDec = R * otswap::maths::cos(centre);
        const double pad = 0.1 * cell / RcosDec;
        const double arc = (raMax - raMin) + 2. * pad;
        if (RcosDec * arc > RcosDec * 2. * kPi - 2. * cell) {
          const double sectors = std::round((RcosDec * kTwoPi) / cell);
          if (sectors == 0.)
            throw Error("ring " + std::to_string(j) + " (from 0) of shell " + std::to_string(i) +
                        " (from 0)" + of + " is a full circle with no right-ascension sector; use a smaller cell "
                        "or a central sphere");
          if (sectors > (double)kMaxCellsAxis)
            throw Error(limit_message(name, "sectors in ring " + std::to_string(j) + " (from 0) of shell " +
                                      std::to_string(i) + " (from 0)"));
          add_ring(lower, upper, centre, false, true, 0., kTwoPi / sectors, (std::int64_t)sectors);
        }
        else {
          const double ratio = (RcosDec * arc) / cell;
          const double sectors = std::ceil(ratio);
          if (sectors > (double)kMaxCellsAxis)
            throw Error(limit_message(name, "sectors in ring " + std::to_string(j) + " (from 0) of shell " +
                                      std::to_string(i) + " (from 0)"));
          const double raStep = arc / sectors;
          add_ring(lower, upper, centre, false, false, (raMin - pad) - (sectors * raStep - arc) / 2., raStep,
                   (std::int64_t)sectors);
        }
      }
      lower = upper;
    }
  }
  g.shellFirstRing.push_back((std::int64_t)g.ringDecMin.size());

  const std::int64_t nRings = (std::int64_t)g.ringDecMin.size();
  g.ringFirstCell.assign((std::size_t)nRings + 1, 0);
  for (std::int64_t r = 0; r < nRings; ++r)
    g.ringFirstCell[(std::size_t)r + 1] = g.ringFirstCell[(std::size_t)r] + ringSectors[(std::size_t)r];
  const std::int64_t nCells = g.ringFirstCell.back();

  std::vector<std::int64_t> ringShell((std::size_t)nRings);
  for (std::size_t i = 0; i < nShells; ++i)
    for (std::int64_t r = g.shellFirstRing[i]; r < g.shellFirstRing[i + 1]; ++r) ringShell[(std::size_t)r] = (std::int64_t)i;

  DivergenceGrid& cells = g.cells;
  cells.geometry = otswap::Geometry::Lightcone;
  cells.nCells = (std::size_t)nCells;
  cells.distanceMin = record.distanceMin;
  cells.distanceMax = record.distanceMax;
  cells.raMin = record.raMin;
  cells.raMax = record.raMax;
  cells.decMin = record.decMin;
  cells.decMax = record.decMax;
  cells.cellSizeMps = record.cellSizeMps;
  cells.mpsProfile = record.mpsProfile;
  const std::size_t n = (std::size_t)nCells;
  cells.cellRaMin.resize(n);
  cells.cellRaMax.resize(n);
  cells.cellDecMin.resize(n);
  cells.cellDecMax.resize(n);
  cells.cellDistanceMin.resize(n);
  cells.cellDistanceMax.resize(n);
  cells.centre.resize(3 * n);
  cells.centreSky.resize(3 * n);
  cells.volume.resize(n);
  cells.masked.resize(n);

  std::exception_ptr failure;
#pragma omp parallel for schedule(dynamic, 16)
  for (std::int64_t r = 0; r < nRings; ++r) {
    try {
      const std::size_t ring = (std::size_t)r;
      const std::size_t i = (std::size_t)ringShell[ring];
      const double rIn = g.shellEdge[i], rOut = g.shellEdge[i + 1];
      const double decLo = g.ringDecMin[ring], decHi = g.ringDecMax[ring];
      const bool sphereCell = i == 0 && g.centralSphere;
      const double shellVolume = sphereCell ? 4. * kPi * rOut * rOut * rOut / 3.
                                            : 4. * kPi * (rOut * rOut * rOut - rIn * rIn * rIn) / 3.;
      const double ringFraction = sphereCell ? 1.
        : (otswap::maths::sin(decHi) - otswap::maths::sin(decLo)) / 2.;
      const double decCentre = ringCentre[ring];
      const double distance = g.shellCentre[i];
      const double cosDec = otswap::maths::cos(decCentre), sinDec = otswap::maths::sin(decCentre);
      const std::int64_t sectors = ringSectors[ring];
      for (std::int64_t k = 0; k < sectors; ++k) {
        const std::size_t c = (std::size_t)(g.ringFirstCell[ring] + k);
        const double raLo = ringRaStart[ring] + (double)k * ringRaStep[ring];
        const double raHi = g.ringFull[ring] && k == sectors - 1 ? kTwoPi
                                                                  : ringRaStart[ring] + (double)(k + 1) * ringRaStep[ring];
        cells.cellRaMin[c] = raLo;
        cells.cellRaMax[c] = raHi;
        cells.cellDecMin[c] = decLo;
        cells.cellDecMax[c] = decHi;
        cells.cellDistanceMin[c] = rIn;
        cells.cellDistanceMax[c] = rOut;
        cells.volume[c] = shellVolume * ringFraction * (ringRaStep[ring] / kTwoPi);
        const double ra = fold((raLo + raHi) / 2.);
        cells.centreSky[3*c] = ra;
        cells.centreSky[3*c+1] = decCentre;
        cells.centreSky[3*c+2] = g.shellRedshift[i];
        cells.centre[3*c] = distance * cosDec * otswap::maths::cos(ra);
        cells.centre[3*c+1] = distance * cosDec * otswap::maths::sin(ra);
        cells.centre[3*c+2] = distance * sinDec;
        cells.masked[c] = mask.allows(ra, decCentre) ? 0 : 1;
      }
    }
    catch (...) {
#pragma omp critical(bitvf_lightcone_grid)
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
  return g;
}


// ============================================================================


std::int64_t bitvf::internal::lightcone_cell (const LightconeGrid& grid, const double distance, const double ra,
                                              const double dec)
{
  const std::vector<double>& edge = grid.shellEdge;
  const std::int64_t shell = (std::int64_t)(std::upper_bound(edge.begin(), edge.end(), distance) - edge.begin()) - 1;
  if (shell < 0 || shell >= (std::int64_t)edge.size() - 1) return -1;

  const std::int64_t firstRing = grid.shellFirstRing[(std::size_t)shell];
  const std::int64_t endRing = grid.shellFirstRing[(std::size_t)shell + 1];
  const auto decBegin = grid.ringDecMin.begin() + firstRing;
  const std::int64_t ring = firstRing + (std::int64_t)(std::upper_bound(decBegin, grid.ringDecMin.begin() + endRing, dec) -
                                                 decBegin) - 1;
  if (ring < firstRing) return -1;
  if (!(dec < grid.ringDecMax[(std::size_t)ring]) && !(dec == kHalfPi && grid.ringDecMax[(std::size_t)ring] == kHalfPi))
    return -1;

  const std::int64_t firstCell = grid.ringFirstCell[(std::size_t)ring];
  const std::int64_t endCell = grid.ringFirstCell[(std::size_t)ring + 1];
  const std::vector<double>& raLo = grid.cells.cellRaMin;
  double a = fold(ra);
  if (!grid.ringFull[(std::size_t)ring]) {
    const double centre = (raLo[(std::size_t)firstCell] + grid.cells.cellRaMax[(std::size_t)endCell - 1]) / 2.;
    a += kTwoPi * std::round((centre - a) / kTwoPi);
  }
  const auto raBegin = raLo.begin() + firstCell;
  const std::int64_t cell = firstCell + (std::int64_t)(std::upper_bound(raBegin, raLo.begin() + endCell, a) -
                                                       raBegin) - 1;
  if (cell < firstCell || !(a < grid.cells.cellRaMax[(std::size_t)cell])) return -1;
  return cell;
}
