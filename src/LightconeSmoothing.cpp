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
#include <limits>
#include <vector>

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;
  constexpr double kMargin = 1e-9;
  constexpr double kWhole = std::numeric_limits<double>::infinity();

  struct RingWindow {
    std::int64_t ring;
    double halfRa;
  };

}


// ============================================================================


void bitvf::internal::lightcone_smooth (LightconeGrid& grid, const std::vector<double>& sigmaShell,
                                        const double threshold)
{
  DivergenceGrid& cells = grid.cells;
  if (sigmaShell.empty()) {
    cells.divergenceSmoothed.clear();
    cells.weightSmoothed.clear();
    return;
  }

  const std::int64_t nShells = (std::int64_t)grid.shellCentre.size();
  const std::int64_t nRings = grid.shellFirstRing.back();
  const std::vector<double>& R = grid.shellCentre;

  std::vector<std::int64_t> ringShell((std::size_t)nRings);
  for (std::int64_t i = 0; i < nShells; ++i)
    for (std::int64_t r = grid.shellFirstRing[(std::size_t)i]; r < grid.shellFirstRing[(std::size_t)i + 1]; ++r)
      ringShell[(std::size_t)r] = i;

  std::vector<std::int64_t> shellLo((std::size_t)nShells), shellHi((std::size_t)nShells);
  std::vector<std::int64_t> pairFirst((std::size_t)nShells + 1, 0);
  std::vector<double> pairSin, pairDec;
  for (std::int64_t i = 0; i < nShells; ++i) {
    const std::size_t ii = (std::size_t)i;
    const double reach = 4. * sigmaShell[ii] * (1. + kMargin);
    shellLo[ii] = std::lower_bound(R.begin(), R.end(), R[ii] - reach) - R.begin();
    shellHi[ii] = (std::upper_bound(R.begin(), R.end(), R[ii] + reach) - R.begin()) - 1;
    for (std::int64_t n = shellLo[ii]; n <= shellHi[ii]; ++n) {
      const double rMin = std::min(R[ii], R[(std::size_t)n]);
      const double s = rMin > 0. ? reach / (2. * rMin) : kWhole;
      pairSin.push_back(s);
      pairDec.push_back(s < 1. ? 2. * otswap::maths::asin(s) * (1. + kMargin) : kWhole);
    }
    pairFirst[ii + 1] = (std::int64_t)pairSin.size();
  }

  const std::size_t C = cells.nCells;
  cells.divergenceSmoothed.assign(C, 0.);
  cells.weightSmoothed.assign(C, 0.);

  const std::vector<double>& raMin = cells.cellRaMin;
  const std::vector<double>& raMax = cells.cellRaMax;
  const std::vector<double>& centre = cells.centre;
  const std::vector<std::uint8_t>& masked = cells.masked;
  const std::vector<double>& divergence = cells.divergence;
  const std::vector<double>& weight = cells.weight;

#pragma omp parallel
  {
    std::vector<RingWindow> window;

#pragma omp for schedule(dynamic)
    for (std::int64_t r = 0; r < nRings; ++r) {
      const std::size_t ring = (std::size_t)r;
      const std::int64_t first = grid.ringFirstCell[ring], last = grid.ringFirstCell[ring + 1];
      if (first == last) continue;
      const std::size_t i = (std::size_t)ringShell[ring];
      const double decT = cells.centreSky[3 * (std::size_t)first + 1];
      const double cosT = std::abs(otswap::maths::cos(decT));
      const double sigma = sigmaShell[i];
      const double two = 2. * sigma * sigma;
      const double reach2 = 16. * sigma * sigma;

      window.clear();
      for (std::int64_t n = shellLo[i]; n <= shellHi[i]; ++n) {
        const std::size_t p = (std::size_t)(pairFirst[i] + n - shellLo[i]);
        const double s = pairSin[p], dDec = pairDec[p];
        const auto begin = grid.ringDecMax.begin() + grid.shellFirstRing[(std::size_t)n];
        const std::int64_t end = grid.shellFirstRing[(std::size_t)n + 1];
        for (std::int64_t nr = std::lower_bound(begin, grid.ringDecMax.begin() + end, decT - dDec) - grid.ringDecMax.begin();
             nr < end && grid.ringDecMin[(std::size_t)nr] <= decT + dDec; ++nr) {
          const std::int64_t nFirst = grid.ringFirstCell[(std::size_t)nr];
          if (nFirst == grid.ringFirstCell[(std::size_t)nr + 1]) continue;
          const double cosMin = std::min(cosT, std::abs(otswap::maths::cos(cells.centreSky[3 * (std::size_t)nFirst + 1])));
          double half = kWhole;
          if (s < 1. && cosMin > 1e-12 && s / cosMin < 1.) {
            half = 2. * otswap::maths::asin(s / cosMin) * (1. + kMargin);
            if (2. * half >= kTwoPi) half = kWhole;
          }
          window.push_back({nr, half});
        }
      }

      for (std::int64_t c = first; c < last; ++c) {
        const std::size_t cc = (std::size_t)c;
        if (masked[cc]) continue;
        const double xt = centre[3*cc], yt = centre[3*cc+1], zt = centre[3*cc+2];
        const double ra = cells.centreSky[3*cc];
        double num = 0., den = 0., sum = 0.;

        auto gather = [&] (const std::int64_t b, const std::int64_t e) {
          for (std::int64_t n = b; n < e; ++n) {
            const std::size_t nn = (std::size_t)n;
            if (masked[nn]) continue;
            const double dx = xt - centre[3*nn], dy = yt - centre[3*nn+1], dz = zt - centre[3*nn+2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (!(d2 < reach2)) continue;
            const double K = otswap::maths::exp(-d2 / two);
            num += K * (weight[nn] * divergence[nn]);
            den += K * weight[nn];
            sum += K;
          }
        };

        for (const RingWindow& w : window) {
          const std::int64_t rf = grid.ringFirstCell[(std::size_t)w.ring];
          const std::int64_t rl = grid.ringFirstCell[(std::size_t)w.ring + 1];
          if (w.halfRa == kWhole) {
            gather(rf, rl);
            continue;
          }
          const double a0 = raMin[(std::size_t)rf];
          double offset = ra - w.halfRa - a0;
          while (offset < 0.) offset += kTwoPi;
          while (offset >= kTwoPi) offset -= kTwoPi;
          const double lo = a0 + offset, hi = lo + 2. * w.halfRa;
          const auto maxBegin = raMax.begin() + rf, maxEnd = raMax.begin() + rl;
          const auto minBegin = raMin.begin() + rf, minEnd = raMin.begin() + rl;
          const std::int64_t b1 = rf + (std::lower_bound(maxBegin, maxEnd, lo) - maxBegin);
          const std::int64_t e1 = rf + (std::upper_bound(minBegin, minEnd, hi) - minBegin);
          const std::int64_t b2 = rf + (std::lower_bound(maxBegin, maxEnd, lo - kTwoPi) - maxBegin);
          const std::int64_t e2 = rf + (std::upper_bound(minBegin, minEnd, hi - kTwoPi) - minBegin);
          if (b2 >= e2) gather(b1, e1);
          else if (b1 >= e1) gather(b2, e2);
          else if (e2 >= b1) gather(std::min(b1, b2), std::max(e1, e2));
          else {
            gather(b2, e2);
            gather(b1, e1);
          }
        }

        cells.divergenceSmoothed[cc] = den > threshold ? num / den : 0.;
        cells.weightSmoothed[cc] = den / sum;
      }
    }
  }
}
