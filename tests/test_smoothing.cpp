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
 *  @file tests/test_smoothing.cpp
 *
 *  @brief The Gaussian smoothing of a box search: no smoothing at sigma
 *  0, the unsmoothed field kept, agreement with the direct weighted mean
 *  over the cube within the rounding bound, the 1e-12 threshold, the
 *  half-window clamp, and its independence of the number of threads.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "otswap/Maths.h"

#include "internal.h"
#include "check.h"

using namespace bitvf;

namespace {

  // The hashes of the smoothed divergence of the generated input, cells of
  // 1 and 0.5 mps, sigma 1 mps, recorded once. A different value on some
  // platform means results are not portable there.
  constexpr std::uint64_t kRecordedHash1 = 0x5d8a2916633921e9ull;
  constexpr std::uint64_t kRecordedHash05 = 0x5013e745876ba35dull;

  constexpr double kU = 0x1p-53;

  void threads (const int n)
  {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
  }

  /// A box grid with lower corner 0, nCellsAxis cells of side s, mps 1;
  /// the per-cell arrays empty.
  DivergenceGrid box (const std::array<std::size_t, 3>& nCellsAxis, const double s)
  {
    DivergenceGrid g;
    g.geometry = otswap::Geometry::Box;
    g.mps = 1.;
    g.cellSize = s;
    g.lower = {0., 0., 0.};
    g.upper = {nCellsAxis[0] * s, nCellsAxis[1] * s, nCellsAxis[2] * s};
    g.nCellsAxis = nCellsAxis;
    g.nCells = nCellsAxis[0] * nCellsAxis[1] * nCellsAxis[2];
    return g;
  }

  std::size_t row (const std::size_t n, const std::size_t i, const std::size_t j,
                   const std::size_t k)
  {
    return (i * n + j) * n + k;
  }

  /// A grid of random divergence in [-1, 1) and weights in [1/4, 2); with
  /// wide, about one cell in 16 has its weight multiplied by 1e4 and one in
  /// 16 its divergence by 1e6.
  DivergenceGrid random_grid (const std::array<std::size_t, 3>& nCellsAxis, const double s, std::mt19937& rng,
                              const bool wide)
  {
    auto uniform = [&rng] (const double lo, const double hi) {
      return lo + (hi - lo) * (rng() / 4294967296.);
    };
    DivergenceGrid g = box(nCellsAxis, s);
    g.divergence.resize(g.nCells);
    g.weight.resize(g.nCells);
    for (std::size_t c = 0; c < g.nCells; ++c) {
      g.divergence[c] = uniform(-1., 1.);
      g.weight[c] = uniform(0.25, 2.);
      if (wide && rng() % 16 == 0) g.weight[c] *= 1e4;
      if (wide && rng() % 16 == 0) g.divergence[c] *= 1e6;
    }
    return g;
  }

  DivergenceGrid random_grid (const std::size_t n, const double s, std::mt19937& rng, const bool wide)
  {
    return random_grid({n, n, n}, s, rng, wide);
  }

  /// The direct form, per cell: the weighted mean over the cube of half-width
  /// h cut at the faces, summed in ascending (i, j, k) order with total
  /// weight ((K_i K_j) K_k) w; the sum of the weights; the weighted mean of
  /// |D|; m_x + m_y + m_z and m_x m_y m_z, the window's cells. h is
  /// halfWidth when given, otherwise from 3.5 sigma / s as the header states.
  struct Direct {
    std::vector<double> out, den, mean;
    std::vector<std::size_t> axes, cells;
  };

  Direct direct (const DivergenceGrid& g, const double sigmaLength,
                 const std::optional<double> halfWidth = std::nullopt)
  {
    const std::int64_t nx = (std::int64_t)g.nCellsAxis[0];
    const std::int64_t ny = (std::int64_t)g.nCellsAxis[1];
    const std::int64_t nz = (std::int64_t)g.nCellsAxis[2];
    const double s = g.cellSize;
    double hd = 3.5 * sigmaLength / s;
    if (std::fabs(hd - std::round(hd)) <= 1e-9 * std::round(hd)) hd = std::round(hd);
    if (halfWidth) hd = *halfWidth;
    const std::int64_t hx = (std::int64_t)std::min(hd, (double)(nx - 1));
    const std::int64_t hy = (std::int64_t)std::min(hd, (double)(ny - 1));
    const std::int64_t hz = (std::int64_t)std::min(hd, (double)(nz - 1));
    const std::int64_t h = std::max(hx, std::max(hy, hz));
    const double two = 2 * sigmaLength * sigmaLength;
    std::vector<double> K(2 * h + 1);
    for (std::int64_t i = 0; i < 2 * h + 1; ++i) {
      const double dist = (double)(i - h) * s;
      K[i] = otswap::maths::exp(-(dist * dist) / two);
    }

    Direct d;
    d.out.resize(g.nCells);
    d.den.resize(g.nCells);
    d.mean.resize(g.nCells);
    d.axes.resize(g.nCells);
    d.cells.resize(g.nCells);
    for (std::int64_t i = 0; i < nx; ++i)
      for (std::int64_t j = 0; j < ny; ++j)
        for (std::int64_t k = 0; k < nz; ++k) {
          const std::int64_t i0 = std::max<std::int64_t>(0, i - hx), i1 = std::min(nx - 1, i + hx);
          const std::int64_t j0 = std::max<std::int64_t>(0, j - hy), j1 = std::min(ny - 1, j + hy);
          const std::int64_t k0 = std::max<std::int64_t>(0, k - hz), k1 = std::min(nz - 1, k + hz);
          double num = 0., den = 0., abs = 0.;
          for (std::int64_t ii = i0; ii <= i1; ++ii)
            for (std::int64_t jj = j0; jj <= j1; ++jj) {
              const double wxy = K[ii - i + h] * K[jj - j + h];
              for (std::int64_t kk = k0; kk <= k1; ++kk) {
                const std::size_t cc = (std::size_t)((ii * ny + jj) * nz + kk);
                const double total = (wxy * K[kk - k + h]) * g.weight[cc];
                num += g.divergence[cc] * total;
                den += total;
                abs += std::fabs(g.divergence[cc]) * total;
              }
            }
          const std::size_t c = (std::size_t)((i * ny + j) * nz + k);
          d.out[c] = den > 1e-12 ? num / den : 0.;
          d.den[c] = den;
          d.mean[c] = den > 0. ? abs / den : 0.;
          d.axes[c] = (std::size_t)((i1 - i0 + 1) + (j1 - j0 + 1) + (k1 - k0 + 1));
          d.cells[c] = (std::size_t)((i1 - i0 + 1) * (j1 - j0 + 1) * (k1 - k0 + 1));
        }
    return d;
  }

  double gamma (const std::size_t k)
  {
    return k * kU / (1 - k * kU);
  }

  /// The bound beta_c between the separable output sep of cell c and a value
  /// other of the direct form, widened by 2^-20 for the rounding of the mean
  /// of |D|.
  double bound (const Direct& d, const std::size_t c, const double sep, const double other)
  {
    const double as = gamma(d.axes[c] + 1), ad = gamma(d.cells[c] + 3);
    const double mu = (3. * d.axes[c] + 2. * d.cells[c] + 8.) * 0x1p-1074 / d.den[c];
    return (2 * (as / (1 - as) + ad / (1 - ad)) * d.mean[c] + kU * (std::fabs(sep) + std::fabs(other)) +
            mu) * (1 + 0x1p-20);
  }

  /// Every cell of out against the direct form: 0 on both sides where the
  /// weights sum to 1e-12 or less, within beta_c elsewhere; the cells equal
  /// bit for bit counted in equal.
  bool agrees (const std::vector<double>& out, const Direct& d, std::size_t& equal)
  {
    bool ok = true;
    equal = 0;
    for (std::size_t c = 0; c < out.size(); ++c) {
      if (out[c] == d.out[c]) ++equal;
      if (d.den[c] <= 1e-12) {
        ok = ok && out[c] == 0. && d.out[c] == 0.;
        continue;
      }
      ok = ok && std::fabs(out[c] - d.out[c]) <= bound(d, c, out[c], d.out[c]);
    }
    return ok;
  }

  /// The smoothed field of a copy of g.
  std::vector<double> smoothed (const DivergenceGrid& g, const double sigmaLength)
  {
    DivergenceGrid copy = g;
    internal::box_smooth(copy, sigmaLength, nullptr);
    return copy.divergence;
  }

  /// FNV-1a over raw bytes.
  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  /// Generated input: tracers uniform in [0, 100)^3, R = 4, displacements
  /// uniform within 6 of zero, about 5% invalid; limits the bounding box
  /// of the tracers and of randoms uniform in [-8, 108)^3, which covers
  /// every endpoint.
  struct Generated {
    static constexpr std::size_t N = 20000;
    static constexpr unsigned R = 4;
    std::vector<double> tracers, displacement, randoms;
    std::vector<std::uint8_t> valid;
    std::vector<std::int64_t> ids;
    otswap::BoxLimits limits;
    Config config;

    Generated ()
    {
      std::mt19937 rng(11);
      auto uniform = [&rng] (const double lo, const double hi) {
        return lo + (hi - lo) * (rng() / 4294967296.);
      };
      for (std::size_t i = 0; i < 3 * N; ++i) tracers.push_back(uniform(0., 100.));
      for (std::size_t i = 0; i < 3 * R * N; ++i) displacement.push_back(uniform(-6., 6.));
      for (std::size_t i = 0; i < R * N; ++i) valid.push_back(rng() % 20 != 0);
      for (std::size_t i = 0; i < 3 * N; ++i) randoms.push_back(uniform(-8., 108.));
      limits.lower = {tracers[0], tracers[1], tracers[2]};
      limits.upper = limits.lower;
      for (const std::vector<double>* points : {&tracers, &randoms})
        for (std::size_t i = 0; i < points->size(); ++i) {
          limits.lower[i%3] = std::min(limits.lower[i%3], (*points)[i]);
          limits.upper[i%3] = std::max(limits.upper[i%3], (*points)[i]);
        }
      config.log = nullptr;
    }

    /// The hash of the smoothed divergence on cells of cellSizeMps, sigma 1
    /// mps; the search's sigmaLength in sigmaLength.
    std::uint64_t hash (const double cellSizeMps, double& sigmaLength) const
    {
      const internal::BoxSearch s = internal::check_box(tracers, displacement, valid, 4., config,
                                                        limits, ids, nullptr);
      DivergenceGrid g = internal::box_grid(s, cellSizeMps);
      internal::box_divergence(s, g);
      internal::box_smooth(g, s.sigmaLength, nullptr);
      sigmaLength = s.sigmaLength;

      std::uint64_t h = 0xcbf29ce484222325ull;
      hash_bytes(h, g.divergence.data(), g.divergence.size() * sizeof(double));
      return h;
    }
  };

}

int main ()
{
  threads(4);
  std::mt19937 rng(2026);

  group("sigma 0: the field unchanged, unsmoothed equal to it, nullptr accepted");
  {
    const DivergenceGrid base = random_grid(9, 1., rng, true);
    DivergenceGrid g = base;
    std::vector<double> unsmoothed;
    internal::box_smooth(g, 0., &unsmoothed);
    check(g.divergence == base.divergence, "the field is unchanged");
    check(unsmoothed == base.divergence, "unsmoothed is the field");
    internal::box_smooth(g, 0., nullptr);
    check(g.divergence == base.divergence, "the field is unchanged with nullptr");
  }

  group("sigma > 0: unsmoothed receives the input bit for bit, the output is the same as without it");
  {
    const DivergenceGrid base = random_grid(9, 1., rng, true);
    DivergenceGrid g = base;
    std::vector<double> unsmoothed;
    internal::box_smooth(g, 1., &unsmoothed);
    check(unsmoothed == base.divergence, "unsmoothed is the input");
    check(g.divergence == smoothed(base, 1.), "the output does not depend on unsmoothed");
  }

  group("half-width 0 (sigma below s / 3.5): equal to the direct form bit for bit");
  for (const bool wide : {false, true}) {
    const DivergenceGrid base = random_grid(9, 1., rng, wide);
    const double sigma = 0.2;
    check(smoothed(base, sigma) == direct(base, sigma).out,
          wide ? "wide weights and values" : "weights in [1/4, 2)");
  }

  group("half-width: 3.5 sigma / s within 1e-9 of an integer takes that integer, s 1 ulp larger included; "
        "farther below it, the integer below");
  {
    const std::size_t n = 17;
    const double sigma = 1.;
    const double exact = 0.5, above = std::nextafter(exact, 1.), farther = 3.5 / (7. * (1. - 1e-8));
    check(3.5 * sigma / exact == 7. && 3.5 * sigma / above < 7. && 3.5 * sigma / farther < 7.,
          "3.5 sigma / s is 7, just below 7, and 7 (1 - 1e-8)");
    for (const double s : {exact, above, farther})
      for (const bool wide : {false, true}) {
        const DivergenceGrid base = random_grid(n, s, rng, wide);
        const std::vector<double> out = smoothed(base, sigma);
        const double h = s == farther ? 6. : 7., other = s == farther ? 7. : 6.;
        std::size_t equal = 0;
        std::ostringstream what;
        what << "3.5 sigma / s = " << std::setprecision(17) << 3.5 * sigma / s << (wide ? ", wide" : "");
        check(agrees(out, direct(base, sigma, h), equal), what.str() + ": half-width " + std::to_string((int)h));
        check(!agrees(out, direct(base, sigma, other), equal),
              what.str() + ": not half-width " + std::to_string((int)other));
      }
  }

  group("separable against direct: every cell within beta_c, faces, edges and corners included");
  for (const std::size_t n : {5, 9, 12})
    for (const std::size_t h : {1, 2, 4, 1000})
      for (const bool wide : {false, true}) {
        const double s = 0.75;
        const double sigma = h == 1000 ? 10. * n * s : (h + 0.5) * s / 3.5;
        const DivergenceGrid base = random_grid(n, s, rng, wide);
        std::size_t equal = 0;
        const std::string what = "n " + std::to_string(n) + ", h " +
                                 (h == 1000 ? std::string("clamped") : std::to_string(h)) +
                                 (wide ? ", wide" : "");
        check(agrees(smoothed(base, sigma), direct(base, sigma), equal), what);
        std::cout << "    " << what << ": " << equal << " of " << base.nCells
                  << " cells equal bit for bit" << std::endl;
      }

  group("threshold: weights summing to 1e-12 or less give 0, just above give num / den, whatever the "
        "largest weight");
  {
    const std::size_t n = 12;
    DivergenceGrid base = random_grid(n, 1., rng, false);
    for (std::size_t i = 3; i <= 8; ++i)
      for (std::size_t j = 3; j <= 8; ++j)
        for (std::size_t k = 3; k <= 8; ++k) base.weight[row(n, i, j, k)] = 0.;
    const double sigma = 2.5 / 3.5;
    const std::vector<double> out = smoothed(base, sigma);
    const Direct d = direct(base, sigma);
    bool zero = true;
    for (std::size_t i = 5; i <= 6; ++i)
      for (std::size_t j = 5; j <= 6; ++j)
        for (std::size_t k = 5; k <= 6; ++k)
          zero = zero && out[row(n, i, j, k)] == 0. && d.den[row(n, i, j, k)] == 0.;
    check(zero, "cells whose window has weight 0 take 0");
    std::size_t equal = 0;
    check(agrees(out, d, equal), "the other cells within beta_c");

    for (const double wMax : {1., 1e20}) {
      DivergenceGrid tiny = random_grid(5, 1., rng, false);
      for (double& w : tiny.weight) w *= wMax / 2.;
      tiny.weight[2] = wMax;
      tiny.weight[0] = 1e-12;
      tiny.weight[1] = std::nextafter(1e-12, 1.);
      const std::vector<double> single = smoothed(tiny, 0.2);
      std::ostringstream what;
      what << "largest weight " << wMax;
      check(single[0] == 0., what.str() + ": a weight sum of exactly 1e-12 gives 0");
      check(single[1] == tiny.weight[1] * tiny.divergence[1] / tiny.weight[1],
            what.str() + ": a weight sum just above gives num / den");
    }
  }

  group("any cell counts: non-cubic grids, an axis of 1 cell included, within beta_c of the direct form");
  for (const std::array<std::size_t, 3> n : {std::array<std::size_t, 3>{5, 9, 13}, {13, 5, 9}, {9, 13, 5},
                                             {1, 7, 11}, {6, 1, 1}, {1, 1, 9}, {2, 3, 4}})
    for (const std::size_t h : {1, 3, 1000})
      for (const bool wide : {false, true}) {
        const double s = 0.75;
        const double sigma = h == 1000 ? 10. * 13 * s : (h + 0.5) * s / 3.5;
        const DivergenceGrid base = random_grid(n, s, rng, wide);
        std::size_t equal = 0;
        std::ostringstream what;
        what << "(" << n[0] << ", " << n[1] << ", " << n[2] << "), h "
             << (h == 1000 ? std::string("clamped") : std::to_string(h)) << (wide ? ", wide" : "");
        check(agrees(smoothed(base, sigma), direct(base, sigma), equal), what.str());
      }

  group("any cell counts: a field constant along one axis gives on each plane the smoothing of one plane, "
        "within the sum of both beta_c");
  for (const int axis : {0, 1, 2}) {
    std::array<std::size_t, 3> n = {5, 9, 13};
    std::array<std::size_t, 3> plane = n;
    plane[(std::size_t)axis] = 1;
    const DivergenceGrid one = random_grid(plane, 1., rng, false);
    DivergenceGrid g = box(n, 1.);
    g.divergence.resize(g.nCells);
    g.weight.resize(g.nCells);
    for (std::size_t i = 0; i < n[0]; ++i)
      for (std::size_t j = 0; j < n[1]; ++j)
        for (std::size_t k = 0; k < n[2]; ++k) {
          std::size_t idx[3] = {i, j, k};
          idx[axis] = 0;
          const std::size_t p = (idx[0] * plane[1] + idx[1]) * plane[2] + idx[2];
          g.divergence[(i * n[1] + j) * n[2] + k] = one.divergence[p];
          g.weight[(i * n[1] + j) * n[2] + k] = one.weight[p];
        }
    const double sigma = 2.5 / 3.5;
    const std::vector<double> flat = smoothed(one, sigma), full = smoothed(g, sigma);
    const Direct d = direct(one, sigma), dFull = direct(g, sigma);
    bool ok = true;
    for (std::size_t i = 0; i < n[0]; ++i)
      for (std::size_t j = 0; j < n[1]; ++j)
        for (std::size_t k = 0; k < n[2]; ++k) {
          std::size_t idx[3] = {i, j, k};
          idx[axis] = 0;
          const std::size_t p = (idx[0] * plane[1] + idx[1]) * plane[2] + idx[2];
          const std::size_t c = (i * n[1] + j) * n[2] + k;
          const double v = full[c];
          ok = ok && std::fabs(v - flat[p]) <= bound(d, p, flat[p], v) + bound(dFull, c, v, flat[p]);
        }
    check(ok, "constant along axis " + std::to_string(axis));
  }

  group("clamp: a large sigma caps the half-width, every cell within beta_c of the global "
        "weighted mean");
  for (const double sigma : {1e100, 9e153}) {
    const std::size_t n = 7;
    const DivergenceGrid base = random_grid(n, 1., rng, true);
    const std::vector<double> out = smoothed(base, sigma);
    const Direct d = direct(base, sigma);
    double num = 0., den = 0.;
    for (std::size_t c = 0; c < base.nCells; ++c) {
      num += base.divergence[c] * base.weight[c];
      den += base.weight[c];
    }
    const double mean = num / den;
    bool within = true;
    for (std::size_t c = 0; c < base.nCells; ++c)
      within = within && d.cells[c] == base.nCells &&
               std::fabs(out[c] - mean) <= bound(d, c, out[c], mean);
    check(within, "sigma " + std::to_string(sigma) + ": the global weighted mean");
  }

  group("invariance: a constant field stays within beta_c of the constant");
  for (const double sigma : {2.5 / 3.5, 4.5 / 3.5}) {
    DivergenceGrid base = random_grid(9, 1., rng, true);
    const double value = -0.37;
    base.divergence.assign(base.nCells, value);
    const std::vector<double> out = smoothed(base, sigma);
    const Direct d = direct(base, sigma);
    bool within = true;
    for (std::size_t c = 0; c < base.nCells; ++c)
      within = within && std::fabs(out[c] - value) <= bound(d, c, out[c], value);
    check(within, "sigma " + std::to_string(sigma));
  }

  group("the same smoothed divergence with 1, 3 and 4 threads, and the recorded hashes");
  {
    const Generated in;
    for (const double cellSizeMps : {1., 0.5}) {
      const std::string what = cellSizeMps == 1. ? "cell 1 mps: " : "cell 0.5 mps: ";
      double sigmaLength = 0.;
      threads(1);
      const std::uint64_t one = in.hash(cellSizeMps, sigmaLength);
      threads(3);
      const std::uint64_t three = in.hash(cellSizeMps, sigmaLength);
      threads(4);
      const std::uint64_t four = in.hash(cellSizeMps, sigmaLength);
      check(sigmaLength == 4., what + "sigmaLength is 1 mps");
      check(one == three && one == four, what + "1, 3 and 4 threads give the same hash");

      const std::uint64_t recorded = cellSizeMps == 1. ? kRecordedHash1 : kRecordedHash05;
      if (one != recorded)
        std::cerr << what << "hash 0x" << std::hex << one << std::dec << std::endl;
      check(one == recorded, what + "the hash equals the recorded one");
    }
  }

  return report("test_smoothing");
}
