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

}


// ============================================================================


bitvf::internal::LightconeFace bitvf::internal::lightcone_face (const LightconeGrid& grid, const std::int64_t cell,
                                                                 const double* x, const double* d, const bool start)
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
