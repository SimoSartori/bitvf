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
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  constexpr double kHalfPi = 3.14159265358979323846 / 2.;
  constexpr double kTwoPi = 2. * 3.14159265358979323846;

  struct Roots {
    std::array<double, 12> t{};
    std::size_t n = 0;

    void add (const double root)
    {
      if (std::isfinite(root) && root > 0. && root < 1.) t[n++] = root;
    }

    void quadratic (const double a, const double halfB, const double c, const double discriminant)
    {
      if (!(discriminant > 0.)) return;
      const double q = -(halfB + (halfB >= 0. ? 1. : -1.) * std::sqrt(discriminant));
      add(q / a);
      add(c / q);
    }
  };

  std::array<double, 3> point_at (const double* x, const double* d, const double t)
  {
    return {x[0] + t * d[0], x[1] + t * d[1], x[2] + t * d[2]};
  }

  bitvf::internal::LightconeOutside outside_at (const bitvf::internal::LightconeGrid& grid, const std::int64_t cell,
                                                 const std::array<double, 3>& p)
  {
    double distance = 0., ra = 0., dec = 0.;
    bitvf::internal::lightcone_sky(p.data(), distance, ra, dec);
    return bitvf::internal::lightcone_outside(grid, cell, distance, ra, dec);
  }

  bool outside (const bitvf::internal::LightconeOutside& o)
  {
    return o.side[0] != 0 || o.side[1] != 0 || o.side[2] != 0;
  }

  constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

  struct SweepErrors {
    std::size_t noFace = kNone;
    std::int64_t noFaceCell = -1;
    bool noFaceStart = true;
    std::size_t lookup = kNone;
    double lookupDistance = 0.;
    std::exception_ptr lookupFailure;
    std::size_t rejected = 0;
  };

  SweepErrors sweep (const bitvf::internal::LightconeSweep& in, bitvf::internal::LightconeGrid& grid,
                     bitvf::internal::LightconeDivergence& result, const std::vector<double>& tracerDistance,
                     std::vector<std::uint64_t>& numerator, const int q, const bool count)
  {
    using bitvf::internal::LightconeFace;
    const std::size_t N = in.nTracers;
    const std::size_t R = in.nRealizations;
    const std::vector<double>& tracers = *in.tracers;
    const std::vector<double>& tracersSky = *in.tracersSky;
    const std::vector<double>& displacement = *in.displacement;
    const std::vector<std::uint8_t>& valid = *in.valid;
    const otswap::DistanceTable& distances = *in.distances;
    bitvf::DivergenceGrid& cells = grid.cells;
    const otswap::MpsProfile& profile = cells.mpsProfile;
    const double s = cells.cellSizeMps;
    const double firstEdge = grid.shellEdge.front();

    SweepErrors errors;
    std::size_t rejected = 0;

#pragma omp parallel for schedule(static) reduction(+:rejected)
    for (std::size_t i = 0; i < N; ++i) {
      const std::int64_t c0 = result.cell[i];
      if (c0 < 0) continue;
      const double* x = &tracers[3*i];
      const double sky[3] = {tracerDistance[i], tracersSky[3*i], tracersSky[3*i+1]};

      for (std::size_t r = 0; r < R; ++r) {
        const std::size_t row = r * N + i;
        if (!valid.empty() && !valid[row]) continue;
        const double* d = &displacement[3*row];
        const double end[3] = {x[0] + d[0], x[1] + d[1], x[2] + d[2]};
        const std::int64_t c1 = bitvf::internal::lightcone_cell_at(grid, end);

        bool reject = c1 < 0 || cells.masked[(std::size_t)c0] || cells.masked[(std::size_t)c1];
        if (!reject && !grid.centralSphere) {
          const double xd = x[0] * d[0] + x[1] * d[1] + x[2] * d[2];
          const double dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
          if (xd < 0. && -xd < dd) {
            const double lx = x[1] * d[2] - x[2] * d[1];
            const double ly = x[2] * d[0] - x[0] * d[2];
            const double lz = x[0] * d[1] - x[1] * d[0];
            reject = lx * lx + ly * ly + lz * lz < firstEdge * firstEdge * dd;
          }
        }
        if (reject) {
          ++rejected;
          continue;
        }
        if (count) result.used[i] = 1;

        if (c1 == c0) {
          if (count) {
#pragma omp atomic
            cells.nDisplacements[(std::size_t)c0] += 2;
          }
          continue;
        }

        const std::size_t index = i * R + r;
        const LightconeFace faces[2] = {bitvf::internal::lightcone_face(grid, c0, x, d, true),
                                        bitvf::internal::lightcone_face(grid, c1, x, d, false, sky)};
        for (std::size_t k = 0; k < 2; ++k) {
          const std::int64_t c = k == 0 ? c0 : c1;
          const std::size_t u = (std::size_t)c;
          const LightconeFace& face = faces[k];
          if (face.coordinate < 0) {
#pragma omp critical(bitvf_lightcone_sweep)
            if (index < errors.noFace) {
              errors.noFace = index;
              errors.noFaceCell = c;
              errors.noFaceStart = k == 0;
            }
            continue;
          }
          const std::array<double, 3>& p = face.point;
          const double distance = std::min(std::max(std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]),
                                                    cells.cellDistanceMin[u]), cells.cellDistanceMax[u]);
          double mps = 0.;
          try {
            mps = profile.at(distances.redshiftAt(distance));
          }
          catch (...) {
#pragma omp critical(bitvf_lightcone_sweep)
            if (index < errors.lookup) {
              errors.lookup = index;
              errors.lookupDistance = distance;
              errors.lookupFailure = std::current_exception();
            }
            continue;
          }
          if (count) {
#pragma omp atomic
            cells.nDisplacements[u] += 1;
          }
          const double term = (face.normal[0] * d[0] + face.normal[1] * d[1] + face.normal[2] * d[2]) / (s * mps);
          const std::uint64_t fixed = (std::uint64_t)(std::int64_t)std::llrint(std::ldexp(term, q));
#pragma omp atomic
          numerator[u] += fixed;
        }
      }
    }

    errors.rejected = rejected;
    return errors;
  }

}


// ============================================================================


bitvf::internal::LightconeFace bitvf::internal::lightcone_face (const LightconeGrid& grid, const std::int64_t cell,
                                                                 const double* x, const double* d, const bool start,
                                                                 const double* sky)
{
  const std::vector<std::int64_t>& firstCell = grid.ringFirstCell;
  const std::size_t ring = (std::size_t)(std::upper_bound(firstCell.begin(), firstCell.end(), cell) -
                                         firstCell.begin()) - 1;
  const std::vector<std::int64_t>& firstRing = grid.shellFirstRing;
  const std::size_t shell = (std::size_t)(std::upper_bound(firstRing.begin(), firstRing.end(), (std::int64_t)ring) -
                                          firstRing.begin()) - 1;
  const std::size_t c = (std::size_t)cell;
  const double raLo = grid.cells.cellRaMin[c];
  const double raHi = grid.ringFull[ring] && grid.cells.cellRaMax[c] == kTwoPi ? 0. : grid.cells.cellRaMax[c];

  const double lx = x[1] * d[2] - x[2] * d[1];
  const double ly = x[2] * d[0] - x[0] * d[2];
  const double lz = x[0] * d[1] - x[1] * d[0];
  const double dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  const double xd = x[0] * d[0] + x[1] * d[1] + x[2] * d[2];
  const double norm = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  const double ll = lx * lx + ly * ly + lz * lz;

  Roots roots;
  for (const bool upper : {false, true}) {
    if (!upper && grid.centralSphere && shell == 0) continue;
    const double r = grid.shellEdge[upper ? shell + 1 : shell];
    roots.quadratic(dd, xd, (norm - r) * (norm + r), r * r * dd - ll);
  }
  for (const double edge : {grid.ringDecMin[ring], grid.ringDecMax[ring]}) {
    if (edge == kHalfPi || edge == -kHalfPi) continue;
    if (edge == 0.) {
      if (d[2] != 0.) roots.add(-x[2] / d[2]);
      continue;
    }
    const double s = otswap::maths::sin(edge), co = otswap::maths::cos(edge);
    const double s2 = s * s, c2 = co * co;
    roots.quadratic((d[0] * d[0] + d[1] * d[1]) * s2 - d[2] * d[2] * c2,
                    (x[0] * d[0] + x[1] * d[1]) * s2 - x[2] * d[2] * c2,
                    (x[0] * x[0] + x[1] * x[1]) * s2 - x[2] * x[2] * c2,
                    s2 * (c2 * (lx * lx + ly * ly) - s2 * lz * lz));
  }
  if (!(grid.ringFull[ring] && firstCell[ring + 1] - firstCell[ring] == 1))
    for (const double edge : {raLo, raHi}) {
      const double nx = -otswap::maths::sin(edge), ny = otswap::maths::cos(edge);
      const double den = nx * d[0] + ny * d[1];
      if (den != 0.) roots.add(-(nx * x[0] + ny * x[1]) / den);
    }

  std::array<double, 14> bound{};
  std::sort(roots.t.begin(), roots.t.begin() + (std::ptrdiff_t)roots.n);
  std::size_t nBounds = 0;
  bound[nBounds++] = 0.;
  for (std::size_t i = 0; i < roots.n; ++i)
    if (roots.t[i] != bound[nBounds - 1]) bound[nBounds++] = roots.t[i];
  bound[nBounds++] = 1.;

  double t = 0.;
  LightconeOutside o;
  bool found = false;
  if (start) {
    o = outside_at(grid, cell, {x[0], x[1], x[2]});
    found = outside(o);
    for (std::size_t i = 0; !found && i + 1 < nBounds; ++i) {
      o = outside_at(grid, cell, point_at(x, d, (bound[i] + bound[i + 1]) / 2.));
      if (outside(o)) {
        t = bound[i];
        found = true;
      }
    }
    if (!found) {
      t = 1.;
      o = outside_at(grid, cell, point_at(x, d, 1.));
    }
  }
  else {
    for (std::size_t i = nBounds - 1; !found && i > 0; --i) {
      o = outside_at(grid, cell, point_at(x, d, (bound[i - 1] + bound[i]) / 2.));
      if (outside(o)) {
        t = bound[i];
        found = true;
      }
    }
    if (!found) {
      t = 0.;
      o = outside_at(grid, cell, {x[0], x[1], x[2]});
      if (!outside(o) && sky != nullptr) o = lightcone_outside(grid, cell, sky[0], sky[1], sky[2]);
    }
  }

  LightconeFace face;
  face.t = t;
  face.point = point_at(x, d, t);
  const std::array<double, 3>& p = face.point;
  constexpr int order[6][2] = {{0, 1}, {1, 1}, {2, 1}, {2, -1}, {1, -1}, {0, -1}};
  for (const auto& face_side : order)
    if (o.side[face_side[0]] == face_side[1]) {
      face.coordinate = face_side[0];
      face.side = face_side[1];
      break;
    }
  if (face.coordinate < 0) return face;

  const double side = face.side;
  if (face.coordinate == 0) {
    const double r = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    face.normal = {side * p[0] / r, side * p[1] / r, side * p[2] / r};
  }
  else if (face.coordinate == 1) {
    const double edge = face.side < 0 ? grid.ringDecMin[ring] : grid.ringDecMax[ring];
    if (edge == 0.) face.normal = {0., 0., side};
    else {
      const double s = otswap::maths::sin(edge), co = otswap::maths::cos(edge);
      const double rho = std::sqrt(p[0] * p[0] + p[1] * p[1]);
      const double cosAz = rho == 0. ? 1. : p[0] / rho, sinAz = rho == 0. ? 0. : p[1] / rho;
      face.normal = {-side * s * cosAz, -side * s * sinAz, side * co};
    }
  }
  else {
    const double edge = face.side < 0 ? raLo : raHi;
    face.normal = {-side * otswap::maths::sin(edge), side * otswap::maths::cos(edge), 0.};
  }
  return face;
}


// ============================================================================


bitvf::internal::LightconeDivergence bitvf::internal::lightcone_divergence (const LightconeSweep& in,
                                                                             LightconeGrid& grid)
{
  const std::size_t N = in.nTracers;
  const std::size_t R = in.nRealizations;
  const std::vector<double>& tracersSky = *in.tracersSky;
  const std::vector<std::uint8_t>& inMask = *in.inMask;
  const std::vector<std::uint8_t>& inCut = *in.inCut;
  const std::vector<double>& displacement = *in.displacement;
  const std::vector<std::uint8_t>& valid = *in.valid;
  const otswap::DistanceTable& distances = *in.distances;
  DivergenceGrid& cells = grid.cells;
  const std::size_t C = cells.nCells;

  LightconeDivergence result;
  result.cell.assign(N, -1);
  result.used.assign(N, 0);
  if (in.reusedGrid) result.outsideGrid.assign(N, 0);
  cells.nTracers.assign(C, 0);
  cells.nDisplacements.assign(C, 0);
  std::vector<double> tracerDistance(N, 0.);

  std::size_t nOutside = 0, nOutsideValid = 0, nRead = 0, outsideNew = kNone, failed = kNone, nonFinite = kNone;
  double maxDisplacement = 0.;
  std::exception_ptr failure;
#pragma omp parallel for schedule(static) reduction(+:nOutside, nOutsideValid, nRead) \
  reduction(min:outsideNew, nonFinite) reduction(max:maxDisplacement)
  for (std::size_t i = 0; i < N; ++i) {
    if (!inMask[i] || !inCut[i]) continue;
    double distance = 0.;
    try {
      distance = distances.distanceAt(tracersSky[3*i+2]);
    }
    catch (...) {
#pragma omp critical(bitvf_lightcone_tracers)
      if (i < failed) {
        failed = i;
        failure = std::current_exception();
      }
      continue;
    }
    const std::int64_t c = lightcone_cell(grid, distance, tracersSky[3*i], tracersSky[3*i+1]);
    if (c < 0) {
      if (in.reusedGrid) {
        result.outsideGrid[i] = 1;
        ++nOutside;
        for (std::size_t r = 0; r < R; ++r)
          if (valid.empty() || valid[r * N + i]) {
            ++nOutsideValid;
            break;
          }
      }
      else outsideNew = std::min(outsideNew, i);
      continue;
    }
    tracerDistance[i] = distance;
    result.cell[i] = c;
#pragma omp atomic
    cells.nTracers[(std::size_t)c] += 1;
    for (std::size_t r = 0; r < R; ++r) {
      const std::size_t row = r * N + i;
      if (!valid.empty() && !valid[row]) continue;
      const double* d = &displacement[3*row];
      if (!std::isfinite(d[0]) || !std::isfinite(d[1]) || !std::isfinite(d[2])) {
        for (std::size_t k = 0; k < 3; ++k)
          if (!std::isfinite(d[k])) nonFinite = std::min(nonFinite, 3*row+k);
        continue;
      }
      maxDisplacement = std::max(maxDisplacement, std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]));
      ++nRead;
    }
  }
  if (failure) {
    try {
      std::rethrow_exception(failure);
    }
    catch (...) {
      std::throw_with_nested(Error("internal error: the distance of the redshift of tracer " + std::to_string(failed) +
                                   " (from 0) cannot be read from the distance table; please report"));
    }
  }
  if (outsideNew != kNone)
    throw Error("internal error: tracer " + std::to_string(outsideNew) +
                " (from 0) lies outside every cell of the grid built on the tracers; please report");
  if (nonFinite != kNone) {
    const std::size_t row = nonFinite / 3;
    throw Error("the displacements hold a non-finite value for realization " + std::to_string(row / N) +
                " (from 0), object " + std::to_string(row % N) + " (from 0), component " +
                std::to_string(nonFinite % 3) + ", a valid displacement of a tracer the search reads");
  }
  if (nRead == 0)
    throw Error("all " + std::to_string(nOutsideValid) + " kept tracers with a valid displacement lie outside the "
                "reused grid; build the grid from these data or pass one that covers them");
  result.nOutsideGrid = nOutside;
  result.maxDisplacement = maxDisplacement;

  const double B = maxDisplacement / (cells.cellSizeMps * grid.mpsMin);
  int q = scale_exponent(B, kAssumedCountBits);
  std::vector<std::uint64_t> numerator(C, 0);
  const SweepErrors errors = sweep(in, grid, result, tracerDistance, numerator, q, true);

  if (errors.noFace != kNone)
    throw Error("internal error: no face found for the displacement of object " + std::to_string(errors.noFace / R) +
                " (from 0) in realization " + std::to_string(errors.noFace % R) + " (from 0), " +
                (errors.noFaceStart ? "start" : "end") + " cell " + std::to_string(errors.noFaceCell) +
                " (from 0); please report");
  if (errors.lookupFailure) {
    std::ostringstream at;
    at.imbue(std::locale::classic());
    at.precision(17);
    at << errors.lookupDistance;
    try {
      std::rethrow_exception(errors.lookupFailure);
    }
    catch (...) {
      std::throw_with_nested(Error("internal error: the redshift or mps lookup failed at distance " + at.str() +
                                   " for the displacement of object " + std::to_string(errors.lookup / R) +
                                   " (from 0) in realization " + std::to_string(errors.lookup % R) +
                                   " (from 0); please report"));
    }
  }
  if (errors.rejected == nRead)
    throw Error(nRead == 1 ? std::string("the 1 valid displacement was rejected (masked cells, outside or below the grid)")
                           : "all " + std::to_string(nRead) +
                             " valid displacements were rejected (masked cells, outside or below the grid)");
  result.nRejected = errors.rejected;

  std::uint64_t K = 0;
#pragma omp parallel for schedule(static) reduction(max:K)
  for (std::size_t c = 0; c < C; ++c)
    K = std::max(K, cells.nDisplacements[c]);
  if (K > (std::uint64_t{1} << kAssumedCountBits)) {
    int L = kAssumedCountBits;
    while ((std::uint64_t{1} << L) < K) ++L;
    q = scale_exponent(B, L);
    std::fill(numerator.begin(), numerator.end(), 0);
    sweep(in, grid, result, tracerDistance, numerator, q, false);
  }
  result.scaleExponent = q;

  cells.divergence.assign(C, 0.);
  cells.weight.assign(C, 0.);
#pragma omp parallel for schedule(static)
  for (std::size_t c = 0; c < C; ++c) {
    const std::uint64_t n = cells.nDisplacements[c];
    cells.divergence[c] = n > 0 ? std::ldexp((double)(std::int64_t)numerator[c], -q) / (double)n : 0.;
    cells.weight[c] = (double)n / (double)R;
  }

  std::size_t nUsed = 0;
#pragma omp parallel for schedule(static) reduction(+:nUsed)
  for (std::size_t i = 0; i < N; ++i)
    nUsed += result.used[i];
  result.nUsed = nUsed;

  return result;
}


// ============================================================================


bitvf::internal::LightconeFill bitvf::internal::lightcone_fill (const LightconeAdjacency& adjacency, DivergenceGrid& grid)
{
  LightconeFill result;
  std::vector<std::int64_t> empty;
  for (std::size_t c = 0; c < grid.nCells; ++c)
    if (grid.nDisplacements[c] == 0 && !grid.masked[c]) {
      empty.push_back((std::int64_t)c);
      grid.divergence[c] = std::numeric_limits<double>::quiet_NaN();
    }

  std::vector<double> value(empty.size());

  while (!empty.empty()) {
    const std::size_t E = empty.size();
    std::size_t nFilled = 0;

#pragma omp parallel for schedule(dynamic, 64) reduction(+:nFilled)
    for (std::size_t e = 0; e < E; ++e) {
      const std::size_t c = (std::size_t)empty[e];
      const double* x = &grid.centre[3*c];
      double num = 0., den = 0.;
      unsigned sources = 0;
      for (std::int64_t k = adjacency.first[c]; k < adjacency.first[c + 1]; ++k) {
        const std::size_t n = (std::size_t)adjacency.neighbour[(std::size_t)k];
        if (grid.masked[n] || std::isnan(grid.divergence[n])) continue;
        ++sources;
        const double dx = x[0] - grid.centre[3*n], dy = x[1] - grid.centre[3*n+1], dz = x[2] - grid.centre[3*n+2];
        const double factor = 1. / (dx*dx + dy*dy + dz*dz);
        num += grid.divergence[n] * factor;
        den += factor;
      }
      if (sources >= 3) {
        value[e] = num / den;
        ++nFilled;
      }
      else
        value[e] = std::numeric_limits<double>::quiet_NaN();
    }

    if (nFilled == 0) break;
    ++result.passes;
    result.filled += nFilled;

    std::size_t kept = 0;
    for (std::size_t e = 0; e < E; ++e)
      if (!std::isnan(value[e]))
        grid.divergence[(std::size_t)empty[e]] = value[e];
      else
        empty[kept++] = empty[e];
    empty.resize(kept);
  }

  for (const std::int64_t c : empty) grid.divergence[(std::size_t)c] = 0.;
  result.leftEmpty = empty.size();
  return result;
}
