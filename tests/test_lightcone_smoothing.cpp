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
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bitvf/VoidFinder.h"
#include "otswap/Maths.h"

#include "check.h"
#include "internal.h"

using namespace bitvf;
using internal::LightconeGrid;

namespace {

  constexpr std::uint64_t kRecordedHash = 0xa4340a829441122dull;
  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  /// One shell of radius R holding one full ring at the equator, split into n equal sectors from RA = 0.
  LightconeGrid ring (const double R, const std::size_t n)
  {
    LightconeGrid g;
    DivergenceGrid& cells = g.cells;
    cells.geometry = otswap::Geometry::Lightcone;
    g.shellCentre = {R};
    g.shellFirstRing = {0, 1};
    g.ringDecMin = {-0.1};
    g.ringDecMax = {0.1};
    g.ringFull = {1};
    g.ringPole = {0};
    g.ringFirstCell = {0, (std::int64_t)n};
    cells.nCells = n;
    for (std::size_t k = 0; k < n; ++k) {
      const double lo = kTwoPi * (double)k / (double)n, hi = k + 1 == n ? kTwoPi : kTwoPi * (double)(k + 1) / (double)n;
      const double ra = (lo + hi) / 2.;
      cells.cellRaMin.push_back(lo);
      cells.cellRaMax.push_back(hi);
      cells.centreSky.insert(cells.centreSky.end(), {ra, 0., 0.1});
      cells.centre.insert(cells.centre.end(), {R * otswap::maths::cos(ra), R * otswap::maths::sin(ra), 0.});
    }
    cells.masked.assign(n, 0);
    cells.divergence.assign(n, 0.);
    cells.weight.assign(n, 1.);
    return g;
  }

  double kernel (const DivergenceGrid& cells, const std::size_t m, const std::size_t n, const double sigma)
  {
    const double dx = cells.centre[3*m] - cells.centre[3*n], dy = cells.centre[3*m+1] - cells.centre[3*n+1],
                 dz = cells.centre[3*m+2] - cells.centre[3*n+2];
    const double d2 = dx * dx + dy * dy + dz * dz;
    return d2 < 16. * sigma * sigma ? otswap::maths::exp(-d2 / (2. * sigma * sigma)) : 0.;
  }

  /// The smoothing by its definition: every pair of cells, in ascending index.
  void reference (const LightconeGrid& g, const std::vector<double>& sigmaShell, const double threshold,
                  std::vector<double>& smoothed, std::vector<double>& weight)
  {
    const DivergenceGrid& cells = g.cells;
    smoothed.assign(cells.nCells, 0.);
    weight.assign(cells.nCells, 0.);
    std::vector<std::size_t> shell(cells.nCells);
    for (std::size_t i = 0; i + 1 < g.shellFirstRing.size(); ++i)
      for (std::int64_t r = g.shellFirstRing[i]; r < g.shellFirstRing[i + 1]; ++r)
        for (std::int64_t c = g.ringFirstCell[(std::size_t)r]; c < g.ringFirstCell[(std::size_t)r + 1]; ++c)
          shell[(std::size_t)c] = i;
    for (std::size_t m = 0; m < cells.nCells; ++m) {
      if (cells.masked[m]) continue;
      const double sigma = sigmaShell[shell[m]];
      double num = 0., den = 0., sum = 0.;
      for (std::size_t n = 0; n < cells.nCells; ++n) {
        if (cells.masked[n]) continue;
        const double dx = cells.centre[3*m] - cells.centre[3*n], dy = cells.centre[3*m+1] - cells.centre[3*n+1],
                     dz = cells.centre[3*m+2] - cells.centre[3*n+2];
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (!(d2 < 16. * sigma * sigma)) continue;
        const double K = otswap::maths::exp(-d2 / (2. * sigma * sigma));
        num += K * (cells.weight[n] * cells.divergence[n]);
        den += K * cells.weight[n];
        sum += K;
      }
      smoothed[m] = den > threshold ? num / den : 0.;
      weight[m] = den / sum;
    }
  }

  /// A grid built from objects spread uniformly over a redshift range and a patch of sky.
  LightconeGrid built_grid (const double raLo, const double raHi, const double decLo, const double decHi,
                            const double zLo, const double zHi, const double cellSizeMps)
  {
    std::mt19937 rng(1);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < 4000; ++i) {
      sky.push_back(uniform(raLo, raHi));
      sky.push_back(otswap::maths::asin(uniform(otswap::maths::sin(decLo), otswap::maths::sin(decHi))));
      sky.push_back(uniform(zLo, zHi));
    }
    const internal::LightconeExtent e = internal::lightcone_extent(sky, {}, {}, {}, table());
    DivergenceGrid record;
    record.geometry = otswap::Geometry::Lightcone;
    record.distanceMin = e.distanceMin;
    record.distanceMax = e.distanceMax;
    record.raMin = e.raMin;
    record.raMax = e.raMax;
    record.decMin = e.decMin;
    record.decMax = e.decMax;
    record.cellSizeMps = cellSizeMps;
    record.mpsProfile.redshift = {0., 0.5};
    record.mpsProfile.mps = {5., 6.};
    const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    return internal::lightcone_grid(record, table(), mask);
  }

  /// Random divergence and weights, some weights 0, every 23rd cell masked.
  void fill_values (LightconeGrid& g, const unsigned seed)
  {
    DivergenceGrid& cells = g.cells;
    std::mt19937 rng(seed);
    cells.divergence.assign(cells.nCells, 0.);
    cells.weight.assign(cells.nCells, 0.);
    for (std::size_t c = 0; c < cells.nCells; ++c) {
      if (c % 23 == 0) cells.masked[c] = 1;
      cells.divergence[c] = -1. + 2. * (rng() / 4294967296.);
      cells.weight[c] = (double)(rng() % 4);
    }
  }

  std::vector<double> sigma_mps (const LightconeGrid& g, const double sigma)
  {
    Smoothing s;
    s.sigma = sigma;
    return internal::lightcone_sigma(s, g);
  }

  bool matches_reference (LightconeGrid& g, const double sigma)
  {
    const std::vector<double> sigmaShell = sigma_mps(g, sigma);
    internal::lightcone_smooth(g, sigmaShell, 1e-12 * 3.);
    std::vector<double> smoothed, weight;
    reference(g, sigmaShell, 1e-12 * 3., smoothed, weight);
    return g.cells.divergenceSmoothed == smoothed && g.cells.weightSmoothed == weight;
  }

  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  std::uint64_t smoothing_hash ()
  {
    LightconeGrid g = built_grid(0.5, 0.8, -0.1, 0.1, 0.25, 0.28, 2.);
    fill_values(g, 7);
    internal::lightcone_smooth(g, sigma_mps(g, 1.), 1e-12 * 3.);
    std::uint64_t h = 0xcbf29ce484222325ull;
    hash_bytes(h, g.cells.divergenceSmoothed.data(), g.cells.divergenceSmoothed.size() * sizeof(double));
    hash_bytes(h, g.cells.weightSmoothed.data(), g.cells.weightSmoothed.size() * sizeof(double));
    return h;
  }

  std::string run_log (const double sigma, const LengthUnit unit)
  {
    const std::vector<double> tracersSky {1., 0.2, 0.3,  2., -0.2, 0.35,  3., 0.1, 0.4,  4., 0.3, 0.25};
    const std::vector<double> displacement(12, 0.1);
    const std::vector<double> randomsSky {1.5, 0.1, 0.3,  2.5, 0., 0.3,  3.5, 0.2, 0.35,  4.5, -0.1, 0.4,
                                          5., 0.2, 0.3,  6., 0.1, 0.25,  0.5, 0.3, 0.3,  0.2, -0.3, 0.35};
    const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    otswap::MpsProfile profile;
    profile.redshift = {0.1, 0.5};
    profile.mps = {5., 6.};
    Config config;
    std::ostringstream log;
    config.smoothing.sigma = sigma;
    config.smoothing.unit = unit;
    config.divergenceCellSize = 15.;
    config.verbosity = otswap::Verbosity::Detailed;
    config.log = &log;
    throws_naming([&] { findVoidsLightcone(tracersSky, displacement, {}, randomsSky, mask, table(), profile, config); },
                  "not implemented yet", "run with a log");
    return log.str();
  }

}

int main ()
{
  // A uniform divergence stays uniform after the smoothing, and a uniform weight gives a uniform smoothed weight.
  group("a uniform field stays uniform");
  {
    LightconeGrid g = built_grid(0.5, 0.8, -0.1, 0.1, 0.25, 0.28, 2.);
    g.cells.divergence.assign(g.cells.nCells, 0.7);
    g.cells.weight.assign(g.cells.nCells, 2.);
    internal::lightcone_smooth(g, sigma_mps(g, 1.), 1e-12 * 2.);
    bool uniform = true;
    for (std::size_t c = 0; c < g.cells.nCells; ++c)
      uniform = uniform && std::abs(g.cells.divergenceSmoothed[c] - 0.7) <= 1e-14 && g.cells.weightSmoothed[c] == 2.;
    check(uniform, "every cell keeps 0.7 and a smoothed weight of 2");
  }

  // A single non-zero cell spreads to the cells within 4 sigma as its kernel over the kernel sum, and nowhere beyond.
  group("a single non-zero cell spreads as the kernel says");
  {
    LightconeGrid g = built_grid(0.5, 0.8, -0.1, 0.1, 0.25, 0.28, 2.);
    DivergenceGrid& cells = g.cells;
    const std::size_t c0 = cells.nCells / 2;
    cells.divergence.assign(cells.nCells, 0.);
    cells.weight.assign(cells.nCells, 1.);
    cells.divergence[c0] = 3.;
    const std::vector<double> sigmaShell = sigma_mps(g, 1.5);
    internal::lightcone_smooth(g, sigmaShell, 1e-12);
    const auto sigma_of = [&] (const std::size_t c) {
      return sigmaShell[(std::size_t)std::distance(g.shellFirstRing.begin(),
        std::upper_bound(g.shellFirstRing.begin(), g.shellFirstRing.end(),
          std::distance(g.ringFirstCell.begin(),
            std::upper_bound(g.ringFirstCell.begin(), g.ringFirstCell.end(), (std::int64_t)c)) - 1)) - 1];
    };
    const double sigma = sigma_of(c0);
    bool spread = true;
    std::size_t reached = 0;
    for (std::size_t m = 0; m < cells.nCells; ++m) {
      if (kernel(cells, c0, m, sigma) > 0.) ++reached;
      if (kernel(cells, m, c0, sigma_of(m)) == 0. && cells.divergenceSmoothed[m] != 0.) spread = false;
      if (cells.weightSmoothed[m] != 1.) spread = false;
    }
    double sum = 0.;
    for (std::size_t n = 0; n < cells.nCells; ++n) sum += kernel(cells, c0, n, sigma);
    check(spread, "zero beyond 4 sigma, and a smoothed weight of 1 everywhere");
    check(reached > 7, "the cell reaches its neighbours");
    check(cells.divergenceSmoothed[c0] == 1. * (1. * 3.) / sum, "the cell itself keeps its value over the kernel sum");
  }

  // On a full ring the cells on both sides of RA = 0 enter each other's mean, and each cell is counted once.
  group("a full ring wraps across RA = 0");
  {
    LightconeGrid g = ring(100., 12);
    g.cells.divergence[11] = 1.;
    const std::vector<double> sigmaShell {15.};
    internal::lightcone_smooth(g, sigmaShell, 1e-12);
    const DivergenceGrid& cells = g.cells;
    const double k0 = kernel(cells, 0, 0, 15.), k1 = kernel(cells, 0, 1, 15.), k11 = kernel(cells, 0, 11, 15.);
    check(k1 > 0. && k11 > 0. && kernel(cells, 0, 2, 15.) == 0., "the window holds one sector on each side");
    check(cells.divergenceSmoothed[0] == k11 * (1. * 1.) / (k0 + k1 + k11), "the cell after RA = 0 sees the last one");
    check(cells.divergenceSmoothed[10] > 0. && cells.divergenceSmoothed[1] == 0. && cells.divergenceSmoothed[6] == 0.,
          "only the neighbours of the last cell take a value");
    check(cells.weightSmoothed[0] == 1., "each cell counted once");
  }

  // Masked cells take 0 and enter no mean, whatever they hold.
  group("masked cells take 0 and enter no sum");
  {
    LightconeGrid g = ring(100., 12);
    g.cells.divergence = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
    g.cells.masked[1] = 1;
    internal::lightcone_smooth(g, {15.}, 1e-12);
    const std::vector<double> first = g.cells.divergenceSmoothed, firstWeight = g.cells.weightSmoothed;
    g.cells.divergence[1] = 1e6;
    g.cells.weight[1] = 50.;
    internal::lightcone_smooth(g, {15.}, 1e-12);
    std::vector<double> smoothed, weight;
    reference(g, {15.}, 1e-12, smoothed, weight);
    check(g.cells.divergenceSmoothed[1] == 0. && g.cells.weightSmoothed[1] == 0., "a masked cell takes 0 and 0");
    check(g.cells.divergenceSmoothed == first && g.cells.weightSmoothed == firstWeight,
          "the result does not depend on the masked cell's values");
    check(g.cells.divergenceSmoothed == smoothed && g.cells.weightSmoothed == weight, "the sums skip the masked cell");
  }

  // A cell whose weighted kernel sum is at most the threshold takes 0, while its smoothed weight is kept.
  group("below the threshold the smoothed divergence is 0 and the smoothed weight is kept");
  {
    LightconeGrid g = ring(100., 12);
    g.cells.divergence.assign(12, 5.);
    g.cells.weight.assign(12, 0.);
    g.cells.weight[0] = 1e-13;
    g.cells.weight[6] = 1.;
    internal::lightcone_smooth(g, {15.}, 1e-12 * 1.);
    check(g.cells.divergenceSmoothed[0] == 0. && g.cells.weightSmoothed[0] > 0., "below the threshold");
    check(g.cells.divergenceSmoothed[6] == 5., "above it");
    check(g.cells.divergenceSmoothed[3] == 0. && g.cells.weightSmoothed[3] == 0., "no weight around: 0 and 0");
  }

  // The smoothing equals its definition, every pair of cells within 4 sigma, on partial, polar and full-sky grids.
  group("the smoothing equals the sum over every pair of cells within 4 sigma");
  {
    LightconeGrid partial = built_grid(0.5, 1., -0.15, 0.15, 0.25, 0.3, 3.);
    fill_values(partial, 3);
    check(matches_reference(partial, 1.), "a partial grid, sigma 1 mps");
    fill_values(partial, 4);
    check(matches_reference(partial, 2.5), "a partial grid, sigma 2.5 mps");
    LightconeGrid polar = built_grid(0., kTwoPi, 1., kPi / 2., 0.25, 0.3, 6.);
    fill_values(polar, 5);
    check(matches_reference(polar, 4.), "full rings and a polar cap");
    LightconeGrid sky = built_grid(0., kTwoPi, -kPi / 2., kPi / 2., 0.001, 0.03, 3.);
    fill_values(sky, 6);
    check(sky.centralSphere, "the full-sky grid has a central sphere");
    check(matches_reference(sky, 2.), "the full sky with a central sphere");
  }

  // With sigma 0 the smoothed fields are empty.
  group("sigma 0 gives empty smoothed fields");
  {
    LightconeGrid g = ring(100., 12);
    internal::lightcone_smooth(g, {15.}, 1e-12);
    const std::vector<double> sigmaShell = sigma_mps(g, 0.);
    check(sigmaShell.empty(), "no sigma per shell");
    internal::lightcone_smooth(g, sigmaShell, 0.);
    check(g.cells.divergenceSmoothed.empty() && g.cells.weightSmoothed.empty(), "both fields emptied");
  }

  // A sigma too small or too large to build the kernel is refused, in mps naming the first such shell and its redshift.
  group("a sigma too small or too large for the kernel is refused");
  {
    LightconeGrid g = ring(100., 2);
    g.shellCentre = {100., 200.};
    g.shellFirstRing = {0, 1, 1};
    g.shellRedshift = {0.1, 0.2};
    g.shellMps = {1., 1e200};
    Smoothing s;
    s.sigma = 1e150;
    throws_naming([&] { internal::lightcone_sigma(s, g); },
                  "the smoothing sigma is 9.9999999999999998e+149 mps, which is inf in coordinate units", "too large in mps");
    s.sigma = 1e60;
    throws_naming([&] { internal::lightcone_sigma(s, g); }, "at redshift 0.20000000000000001 of shell 1 (from 0), too large",
                  "the overflowing shell");
    g.shellMps = {1e-200, 1.};
    s.sigma = 1e-150;
    throws_naming([&] { internal::lightcone_sigma(s, g); }, "of shell 0 (from 0), too small to build the kernel",
                  "too small in mps");
    s.unit = LengthUnit::Coordinate;
    s.sigma = 5e153;
    throws_naming([&] { internal::lightcone_sigma(s, g); },
                  "the smoothing sigma is 5.0000000000000002e+153 in coordinate units", "16 sigma^2 overflows");
    s.sigma = 1e-170;
    throws_naming([&] { internal::lightcone_sigma(s, g); }, "too small to build the kernel", "too small in coordinate units");
    s.sigma = 3.;
    const std::vector<double> sigmaShell = internal::lightcone_sigma(s, g);
    check(sigmaShell == std::vector<double>{3., 3.}, "in coordinate units every shell takes sigma");
  }

  // The search reports the smoothing at Detailed, with the sigma range in coordinate units when sigma is in mps.
  group("the smoothing line is written at Detailed");
  {
    const std::string mps = run_log(1., LengthUnit::Mps);
    check(mps.find("bitvf: smoothed with sigma 1 mps (") != std::string::npos &&
          mps.find(" in coordinate units) in ") != std::string::npos, "in mps: " + mps);
    const std::string coordinate = run_log(3., LengthUnit::Coordinate);
    check(coordinate.find("bitvf: smoothed with sigma 3 in coordinate units in ") != std::string::npos,
          "in coordinate units: " + coordinate);
    const std::string none = run_log(0., LengthUnit::Mps);
    check(none.find("bitvf: not smoothed (sigma 0)") != std::string::npos, "sigma 0: " + none);
  }

  // The smoothed fields of a fixed grid are identical bit for bit with any number of threads and match the recorded hash.
  group("the smoothing is the same with 1 and 4 threads, and matches the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = smoothing_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = smoothing_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "smoothing hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_smoothing");
}
