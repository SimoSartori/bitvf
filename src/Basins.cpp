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
 *  @file src/Basins.cpp
 *
 *  @brief The catalogues of the voids and overdense basins of a box
 *  search, and the assignment of its tracers.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <initializer_list>
#include <limits>
#include <memory>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "meshsearch/MeshGrid.h"
#include "otswap/Maths.h"

#include "internal.h"

namespace {

  constexpr double kPi = 3.141592653589793;

  /// The position of x in the order of the finite doubles: consecutive
  /// doubles get consecutive integers, and both zeros get 0.
  std::int64_t ordinal (const double x)
  {
    std::int64_t bits;
    std::memcpy(&bits, &x, sizeof bits);
    return bits >= 0 ? bits : -(bits & std::numeric_limits<std::int64_t>::max());
  }

  /// The double at position k of ordinal; +0. for 0.
  double from_ordinal (const std::int64_t k)
  {
    const std::uint64_t bits = k >= 0 ? (std::uint64_t)k : ((std::uint64_t)0 - (std::uint64_t)k) |
                                                           ((std::uint64_t)1 << 63);
    double x;
    std::memcpy(&x, &bits, sizeof x);
    return x;
  }

  /**
   *  @brief The double reached from x by steps of std::nextafter towards
   *  -infinity (down) or +infinity (up) while inside(x) is false, in
   *  O(log steps) evaluations of inside.
   *
   *  inside must be false at x, and monotone along the direction: once true
   *  at some double, true at every double beyond it. A zero reached by
   *  steps is +0. going down and -0. going up, as std::nextafter gives.
   */
  template <typename Inside>
  double step_until (const double x, const bool down, const Inside& inside)
  {
    const std::int64_t start = ordinal(x);
    const std::int64_t sign = down ? -1 : 1;
    std::uint64_t fail = 0, pass = 1;
    while (!inside(from_ordinal(start + sign * (std::int64_t)pass))) {
      fail = pass;
      pass *= 2;
    }
    while (pass - fail > 1) {
      const std::uint64_t mid = fail + (pass - fail) / 2;
      if (inside(from_ordinal(start + sign * (std::int64_t)mid))) pass = mid;
      else fail = mid;
    }
    const std::int64_t end = start + sign * (std::int64_t)pass;
    return end == 0 && !down ? -0. : from_ordinal(end);
  }

}


// ============================================================================


void bitvf::internal::box_basins (const BoxSearch* search, const BoxDivergence* div, const DivergenceGrid& grid,
                                  const BoxBasins& basins, BasinCatalog& voids, BasinCatalog& overdense)
{
  const bool withTracers = search != nullptr;
  const std::int64_t nx = (std::int64_t)grid.nCellsAxis[0];
  const std::int64_t ny = (std::int64_t)grid.nCellsAxis[1];
  const std::int64_t nz = (std::int64_t)grid.nCellsAxis[2];
  const std::int64_t n[3] = {nx, ny, nz};
  const std::int64_t C = (std::int64_t)grid.nCells;
  const double s = grid.cellSize;
  const double m = withTracers ? grid.mps : 0.;
  const double inverse = 1. / s;
  const double avDens = withTracers ? 1. / (m * m * m) : 0.;
  const double* D = grid.divergence.data();
  const std::array<double, 3> lower = grid.lower;

  const BasinSet* sets[2] = {&basins.voids, &basins.overdense};
  BasinCatalog* catalogues[2] = {&voids, &overdense};

  for (int kind = 0; kind < 2; ++kind) {
    const BasinSet& set = *sets[kind];
    BasinCatalog& cat = *catalogues[kind];
    const std::int64_t B = (std::int64_t)set.seed.size();

    cat.nBasins = (std::size_t)B;
    cat.centre.resize(3 * (std::size_t)B);
    cat.centreSky.clear();
    cat.radius.resize((std::size_t)B);
    cat.radiusExt.resize((std::size_t)B);
    cat.edgeFlag.resize((std::size_t)B);
    cat.contamination.resize((std::size_t)B);
    if (withTracers) {
      cat.deltaAv.resize((std::size_t)B);
      cat.delta2Mps.resize((std::size_t)B);
      cat.nTracers.resize((std::size_t)B);
    }
    else {
      cat.deltaAv.clear();
      cat.delta2Mps.clear();
      cat.nTracers.clear();
    }

#pragma omp parallel for schedule(static)
    for (std::int64_t b = 0; b < B; ++b) {
      const std::int64_t c = set.seed[(std::size_t)b];
      const std::int64_t idx[3] = {c / (ny * nz), (c / nz) % ny, c % nz};
      std::int64_t lo[3], hi[3];
      for (int a = 0; a < 3; ++a) {
        lo[a] = std::max<std::int64_t>(0, idx[a] - 1);
        hi[a] = std::min(n[a] - 1, idx[a] + 1);
      }

      double num[3] = {0., 0., 0.};
      double den = 0.;
      for (std::int64_t ii = lo[0]; ii <= hi[0]; ++ii)
        for (std::int64_t jj = lo[1]; jj <= hi[1]; ++jj)
          for (std::int64_t kk = lo[2]; kk <= hi[2]; ++kk) {
            const std::int64_t nb = (ii * ny + jj) * nz + kk;
            if (nb == c) continue;
            num[0] += D[nb] * s * (double)(idx[0] - ii);
            num[1] += D[nb] * s * (double)(idx[1] - jj);
            num[2] += D[nb] * s * (double)(idx[2] - kk);
            den += std::abs(D[nb]);
          }
      if (den == 0.) den = 1.;

      for (int a = 0; a < 3; ++a) {
        const double i = (double)idx[a];
        double x = lower[a] + s * (i + 0.5) + num[a] / den;
        const double faceLo = lower[a] + s * i;
        const double faceHi = lower[a] + s * (i + 1.);
        x = std::max(faceLo, std::min(faceHi, x));
        const double low = lower[a];
        if ((x - low) * inverse >= i + 1.)
          x = step_until(x, true, [low, inverse, i] (const double y) { return (y - low) * inverse < i + 1.; });
        if ((x - low) * inverse < i)
          x = step_until(x, false, [low, inverse, i] (const double y) { return (y - low) * inverse >= i; });
        cat.centre[3 * (std::size_t)b + (std::size_t)a] = x;
      }

      const std::uint64_t nCore = set.nCoreCells[(std::size_t)b];
      cat.radius[(std::size_t)b] = s * otswap::maths::cbrt(3.0 * (double)nCore / (4.0 * kPi));
      cat.radiusExt[(std::size_t)b] =
        s * otswap::maths::cbrt(3.0 * (double)set.nCells[(std::size_t)b] / (4.0 * kPi));
      if (withTracers) {
        const std::uint64_t nCoreTracers = set.nCoreTracers[(std::size_t)b];
        cat.deltaAv[(std::size_t)b] = ((double)nCoreTracers / ((double)nCore * s * s * s)) / avDens - 1.;
        cat.nTracers[(std::size_t)b] = (std::size_t)nCoreTracers;
      }
    }
  }

  const std::size_t nVoid = voids.nBasins;
  const std::size_t nAll = nVoid + overdense.nBasins;
  const std::int64_t* voidIndex = grid.voidIndex.data();
  const std::int64_t* overdenseIndex = grid.overdenseIndex.data();

#ifdef _OPENMP
  const std::size_t nThreads = (std::size_t)omp_get_max_threads();
#else
  const std::size_t nThreads = 1;
#endif
  const std::size_t stride = 3 * nAll;
  std::unique_ptr<std::uint64_t[]> faces(new std::uint64_t[nThreads * stride]);
  std::size_t nTeam = 1;

#pragma omp parallel
  {
#ifdef _OPENMP
    std::uint64_t* out = &faces[(std::size_t)omp_get_thread_num() * stride];
#pragma omp single nowait
    nTeam = (std::size_t)omp_get_num_threads();
#else
    std::uint64_t* out = faces.get();
#endif
    std::fill(out, out + stride, std::uint64_t(0));
    std::uint64_t* outCore = out + nAll;
    std::uint64_t* phy = outCore + nAll;

#pragma omp for schedule(static)
    for (std::int64_t c = 0; c < C; ++c) {
      const bool isVoid = voidIndex[c] >= 0;
      const std::int64_t g = isVoid ? voidIndex[c] : (std::int64_t)nVoid + overdenseIndex[c];
      const bool core = isVoid ? D[c] < 0. : D[c] >= 0.;
      const std::int64_t idx[3] = {c / (ny * nz), (c / nz) % ny, c % nz};
      const std::int64_t step[3] = {ny * nz, nz, 1};
      for (int a = 0; a < 3; ++a)
        for (const std::int64_t d : {-1, 1}) {
          const std::int64_t to = idx[a] + d;
          if (to < 0 || to >= n[a]) {
            out[g] += 1;
            if (core) outCore[g] += 1;
            continue;
          }
          const std::int64_t nb = c + d * step[a];
          const std::int64_t h = voidIndex[nb] >= 0 ? voidIndex[nb] : (std::int64_t)nVoid + overdenseIndex[nb];
          if (h != g) phy[g] += 1;
        }
    }
  }

#pragma omp parallel for schedule(static)
  for (std::int64_t g = 0; g < (std::int64_t)nAll; ++g) {
    std::uint64_t out = 0, outCore = 0, phy = 0;
    for (std::size_t t = 0; t < nTeam; ++t) {
      out += faces[t * stride + (std::size_t)g];
      outCore += faces[t * stride + nAll + (std::size_t)g];
      phy += faces[t * stride + 2 * nAll + (std::size_t)g];
    }
    BasinCatalog& cat = (std::size_t)g < nVoid ? voids : overdense;
    const std::size_t b = (std::size_t)g < nVoid ? (std::size_t)g : (std::size_t)g - nVoid;
    cat.edgeFlag[b] = outCore > 0 ? 2 : out > 0 ? 1 : 0;
    cat.contamination[b] = out + phy > 0 ? (double)out / (double)(out + phy) : 0.;
  }
  faces.reset();

  if (!withTracers) return;

  const std::size_t N = search->nTracers;
  const std::vector<double>& tracers = *search->tracers;
  const std::size_t nRead = N - search->counts.tracersOutsideBox - search->counts.tracersOutsideGrid;
  std::vector<double> X, Y, Z;
  X.reserve(nRead);
  Y.reserve(nRead);
  Z.reserve(nRead);
  for (std::size_t i = 0; i < N; ++i)
    if (div->cell[i] >= 0) {
      X.push_back(tracers[3*i]);
      Y.push_back(tracers[3*i+1]);
      Z.push_back(tracers[3*i+2]);
    }

  std::vector<std::vector<double>> limits(3);
  for (std::size_t a = 0; a < 3; ++a)
    limits[a] = {lower[a] - s, lower[a] + (double)(n[a] + 1) * s};

  std::unique_ptr<meshsearch::MeshGrid> mesh;
  try {
    mesh.reset(new meshsearch::MeshGrid(X, Y, Z, 4 * m, limits));
  }
  catch (const meshsearch::Error& e) {
    std::throw_with_nested(Error(std::string("internal error: meshsearch: ") + e.what()));
  }
  std::vector<double>().swap(X);
  std::vector<double>().swap(Y);
  std::vector<double>().swap(Z);

  const double radius2Mps = 2 * m;
  const double sphere = 4 * kPi * 8 * m * m * m / 3.;
  constexpr std::int64_t none = std::numeric_limits<std::int64_t>::max();
  std::int64_t failed = none;

#pragma omp parallel for schedule(dynamic) reduction(min:failed)
  for (std::int64_t g = 0; g < (std::int64_t)nAll; ++g) {
    BasinCatalog& cat = (std::size_t)g < nVoid ? voids : overdense;
    const std::size_t b = (std::size_t)g < nVoid ? (std::size_t)g : (std::size_t)g - nVoid;
    const double* x = &cat.centre[3 * b];
    try {
      const std::size_t nt = mesh->closeObjects(x[0], x[1], x[2], radius2Mps).size();
      cat.delta2Mps[b] = ((double)nt / sphere) / avDens - 1.;
    }
    catch (...) {
      failed = std::min(failed, g);
    }
  }

  if (failed != none) {
    const bool isVoid = (std::size_t)failed < nVoid;
    const BasinCatalog& cat = isVoid ? voids : overdense;
    const std::size_t b = isVoid ? (std::size_t)failed : (std::size_t)failed - nVoid;
    const double* x = &cat.centre[3 * b];
    const std::string what = std::string("internal error: counting the tracers within 2 mps of the centre of ") +
      (isVoid ? "void " : "overdense basin ") + std::to_string(b);
    try {
      mesh->closeObjects(x[0], x[1], x[2], radius2Mps);
    }
    catch (const meshsearch::Error& e) {
      std::throw_with_nested(Error(what + ": meshsearch: " + e.what()));
    }
    throw Error(what + " failed inside the parallel region");
  }
}


// ============================================================================


bitvf::TracerAssignment bitvf::internal::box_assignment (const BoxSearch& search, const BoxDivergence& div,
                                                         const DivergenceGrid& grid, const BasinCatalog& voids,
                                                         const BasinCatalog& overdense)
{
  const std::int64_t N = (std::int64_t)search.nTracers;
  const std::vector<double>& tracers = *search.tracers;

  TracerAssignment result;
  result.cellIndex.resize((std::size_t)N);
  result.voidIndex.resize((std::size_t)N);
  result.overdenseIndex.resize((std::size_t)N);
  result.distance.resize((std::size_t)N);
  result.inMask.assign((std::size_t)N, 1);
  result.inBox = search.inBox;
  result.used = div.used;

#pragma omp parallel for schedule(static)
  for (std::int64_t i = 0; i < N; ++i) {
    const std::int64_t c = div.cell[(std::size_t)i];
    result.cellIndex[(std::size_t)i] = c;
    if (c < 0) {
      result.voidIndex[(std::size_t)i] = -1;
      result.overdenseIndex[(std::size_t)i] = -1;
      result.distance[(std::size_t)i] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    const std::int64_t v = grid.voidIndex[(std::size_t)c];
    const std::int64_t o = grid.overdenseIndex[(std::size_t)c];
    result.voidIndex[(std::size_t)i] = v;
    result.overdenseIndex[(std::size_t)i] = o;
    const double* centre = v >= 0 ? &voids.centre[3 * (std::size_t)v] : &overdense.centre[3 * (std::size_t)o];
    const double* x = &tracers[3 * (std::size_t)i];
    const double dx = x[0] - centre[0];
    const double dy = x[1] - centre[1];
    const double dz = x[2] - centre[2];
    result.distance[(std::size_t)i] = std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  return result;
}
