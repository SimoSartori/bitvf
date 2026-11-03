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
 *  @file tests/fixtures.h
 *
 *  @brief Inputs shared by the tests of the box entry points: a search
 *  that finds both kinds of basin, and an analytic field on a non-cubic
 *  lattice.
 */

#ifndef BITVF_TESTS_FIXTURES_H
#define BITVF_TESTS_FIXTURES_H

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <optional>
#include <random>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "otswap/Maths.h"

#include "bitvf/VoidFinder.h"

namespace fixtures {

  inline void set_threads (const int n)
  {
#ifdef _OPENMP
    omp_set_num_threads(n);
#else
    (void)n;
#endif
  }

  inline bool same_bits (const std::vector<double>& a, const std::vector<double>& b)
  {
    return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
  }

  /// Every column of two catalogues equal, bit for bit.
  inline bool same_catalog (const bitvf::BasinCatalog& a, const bitvf::BasinCatalog& b)
  {
    return a.nBasins == b.nBasins && same_bits(a.centre, b.centre) && same_bits(a.radius, b.radius) &&
           same_bits(a.radiusExt, b.radiusExt) && same_bits(a.deltaAv, b.deltaAv) &&
           same_bits(a.delta2Mps, b.delta2Mps) && a.edgeFlag == b.edgeFlag &&
           same_bits(a.contamination, b.contamination) && a.nTracers == b.nTracers &&
           same_bits(a.fieldCentre, b.fieldCentre);
  }

  /// FNV-1a over raw bytes.
  inline void fnv (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  template <typename T>
  inline void fnv_vector (std::uint64_t& h, const std::vector<T>& v)
  {
    fnv(h, v.data(), v.size() * sizeof(T));
  }

  /// A search on 3000 tracers uniform in [0, 100)^3 and kOutside more
  /// beyond its limits, over two realizations: the tracers with x below
  /// 50 move by about 30 along x, the others by at most 2, so that the
  /// cells they leave form a positive region whose dips seed overdense
  /// basins, and the cells they reach seed voids. Sigma 0.5 mps, every
  /// output kept, IDs given, three realizations averaged on record.
  struct Translated {

    static constexpr std::size_t kN = 3000;
    static constexpr std::size_t kOutside = 4;
    static constexpr double kMps = 7.;

    std::vector<double> tracers, displacement;
    std::vector<std::uint8_t> valid;
    std::optional<otswap::BoxLimits> limits;
    std::vector<std::int64_t> ids;
    bitvf::Config config;

    Translated ()
    {
      std::mt19937 rng(2026);
      auto uniform = [&rng] (const double lo, const double hi) {
        return lo + (hi - lo) * (rng() / 4294967296.);
      };
      const std::size_t n = kN + kOutside;
      for (std::size_t i = 0; i < 3 * kN; ++i) tracers.push_back(uniform(0., 100.));
      for (std::size_t i = 0; i < 3 * kOutside; ++i) tracers.push_back(150. + (double)i);
      for (std::size_t r = 0; r < 2; ++r)
        for (std::size_t i = 0; i < n; ++i) {
          const bool moved = tracers[3*i] < 50.;
          displacement.push_back(moved ? 30. + uniform(-3., 3.) : uniform(-2., 2.));
          displacement.push_back(moved ? uniform(-3., 3.) : uniform(-2., 2.));
          displacement.push_back(moved ? uniform(-3., 3.) : uniform(-2., 2.));
        }
      otswap::BoxLimits box;
      box.lower = {-5., -5., -5.};
      box.upper = {135., 105., 105.};
      limits = box;
      for (std::size_t i = 0; i < n; ++i) ids.push_back(i % 2 == 0 ? (std::int64_t)(500 + i) : -(std::int64_t)(500 + i));
      config.log = nullptr;
      config.smoothing.sigma = 0.5;
      config.keepDivergence = true;
      config.keepAssignment = true;
      config.nRealizationsMean = 3;
    }

    bitvf::Result run (const bitvf::DivergenceGrid* grid = nullptr) const
    {
      return bitvf::findVoidsBox(tracers, displacement, valid, kMps, config, limits, ids, grid);
    }

    bitvf::Result run (const bitvf::Config& c) const
    {
      return bitvf::findVoidsBox(tracers, displacement, valid, kMps, c, limits, ids);
    }
  };

  /// A field at the centres of an (nx, ny, nz) lattice of side s, lowest
  /// centre at origin: cos(k x) + cos(k y) + cos(k z) with a period of six
  /// cells, plus a tilt along x that lifts the minima from about -3 to
  /// about +3 across the lattice, so that both kinds of basin appear, plus
  /// a noise of 0.01 that breaks the symmetric ties. Rows in lattice order;
  /// weights in [0.5, 2].
  struct LatticeField {

    std::array<std::size_t, 3> n;
    double s;
    std::array<double, 3> origin;
    std::vector<double> centres, values, weights;

    LatticeField (const std::array<std::size_t, 3>& counts = {12, 17, 9}, const double side = 1.5,
                  const std::array<double, 3>& lowest = {-3.25, 10.5, 0.75})
      : n(counts), s(side), origin(lowest)
    {
      std::mt19937 rng(7);
      auto uniform = [&rng] (const double lo, const double hi) {
        return lo + (hi - lo) * (rng() / 4294967296.);
      };
      const double k = 2. * 3.141592653589793 / (6. * s);
      const double tilt = 6. / ((double)n[0] * s);
      for (std::size_t i = 0; i < n[0]; ++i)
        for (std::size_t j = 0; j < n[1]; ++j)
          for (std::size_t l = 0; l < n[2]; ++l) {
            const double x[3] = {origin[0] + s * (double)i, origin[1] + s * (double)j, origin[2] + s * (double)l};
            centres.insert(centres.end(), x, x + 3);
            values.push_back(otswap::maths::cos(k * (double)i * s) + otswap::maths::cos(k * (double)j * s) +
                             otswap::maths::cos(k * (double)l * s) + tilt * (double)i * s + uniform(-0.01, 0.01));
            weights.push_back(uniform(0.5, 2.));
          }
    }
  };

}

#endif
