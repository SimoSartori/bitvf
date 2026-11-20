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

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bitvf/VoidFinder.h"
#include "otswap/Maths.h"
#include "otswap/OT.h"

#include "check.h"
#include "internal.h"

using namespace bitvf;
using internal::LightconeFace;
using internal::LightconeGrid;

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;
  constexpr double kHalfPi = kPi / 2.;

  constexpr std::uint64_t kRecordedHash = 0x8380b08b0325a45cull;

  using Point = std::array<double, 3>;

  struct RingSpec {
    double decMin, decMax;
    std::vector<double> raEdges;
    bool full, pole;
  };

  LightconeGrid hand_grid (const bool centralSphere, const std::vector<double>& shellEdge,
                           const std::vector<std::vector<RingSpec>>& shells)
  {
    LightconeGrid g;
    g.centralSphere = centralSphere;
    g.shellEdge = shellEdge;
    std::int64_t ring = 0, cell = 0;
    for (const std::vector<RingSpec>& rings : shells) {
      g.shellFirstRing.push_back(ring);
      for (const RingSpec& r : rings) {
        g.ringDecMin.push_back(r.decMin);
        g.ringDecMax.push_back(r.decMax);
        g.ringPole.push_back(r.pole ? 1 : 0);
        g.ringFull.push_back(r.full ? 1 : 0);
        g.ringFirstCell.push_back(cell);
        for (std::size_t k = 0; k + 1 < r.raEdges.size(); ++k) {
          g.cells.cellRaMin.push_back(r.raEdges[k]);
          g.cells.cellRaMax.push_back(r.raEdges[k + 1]);
          ++cell;
        }
        ++ring;
      }
    }
    g.shellFirstRing.push_back(ring);
    g.ringFirstCell.push_back(cell);
    g.cells.nCells = (std::size_t)cell;
    return g;
  }

  double meridian ()
  {
    return otswap::maths::atan2(1., 0.);
  }

  const LightconeGrid& toy ()
  {
    static const LightconeGrid g = [] {
      const std::vector<double> ra {1.3, meridian(), 1.85};
      std::vector<RingSpec> rings;
      const double dec[] = {-0.3, -0.15, 0., 0.15, 0.3};
      for (int r = 0; r < 4; ++r) rings.push_back({dec[r], dec[r + 1], ra, false, false});
      return hand_grid(false, {100., 110., 120.}, {rings, rings});
    }();
    return g;
  }

  std::int64_t toy_cell (const int shell, const int ring, const int sector)
  {
    return (shell * 4 + ring) * 2 + sector;
  }

  const LightconeGrid& polar ()
  {
    static const LightconeGrid g = [] {
      const std::vector<double> quarters {0., kHalfPi, kPi, 3. * kHalfPi, kTwoPi};
      const std::vector<RingSpec> shell {{-kHalfPi, -1.2, {0., kTwoPi}, true, true},
                                         {-1.2, 0., quarters, true, false},
                                         {0., 1.2, quarters, true, false},
                                         {1.2, kHalfPi, {0., kTwoPi}, true, true}};
      return hand_grid(true, {0., 50., 60., 70.}, {{{-kHalfPi, kHalfPi, {0., kTwoPi}, true, false}}, shell, shell});
    }();
    return g;
  }

  Point sky_point (const double r, const double dec, const double ra)
  {
    const double c = otswap::maths::cos(dec);
    return {r * c * otswap::maths::cos(ra), r * c * otswap::maths::sin(ra), r * otswap::maths::sin(dec)};
  }

  Point minus (const Point& a, const Point& b)
  {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  }

  Point end_of (const Point& x, const Point& d)
  {
    return {x[0] + d[0], x[1] + d[1], x[2] + d[2]};
  }

  LightconeFace leave (const LightconeGrid& g, const std::int64_t cell, const Point& x, const Point& d)
  {
    return internal::lightcone_face(g, cell, x.data(), d.data(), true);
  }

  LightconeFace enter (const LightconeGrid& g, const Point& x, const Point& d)
  {
    const Point e = end_of(x, d);
    return internal::lightcone_face(g, internal::lightcone_cell_at(g, e.data()), x.data(), d.data(), false);
  }

  bool is_face (const LightconeFace& f, const int coordinate, const int side)
  {
    return f.coordinate == coordinate && f.side == side;
  }

  double norm (const Point& p)
  {
    return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
  }

  bool close (const Point& a, const Point& b, const double tolerance)
  {
    return std::fabs(a[0] - b[0]) <= tolerance && std::fabs(a[1] - b[1]) <= tolerance &&
           std::fabs(a[2] - b[2]) <= tolerance;
  }

  std::size_t ring_of (const LightconeGrid& g, const std::int64_t cell)
  {
    std::size_t ring = 0;
    while (g.ringFirstCell[ring + 1] <= cell) ++ring;
    return ring;
  }

  std::size_t shell_of (const LightconeGrid& g, const std::size_t ring)
  {
    std::size_t shell = 0;
    while (g.shellFirstRing[shell + 1] <= (std::int64_t)ring) ++shell;
    return shell;
  }

  double angle_gap (const double a, const double b)
  {
    const double gap = std::fmod(std::fabs(a - b), kTwoPi);
    return std::min(gap, kTwoPi - gap);
  }

  bool on_face (const LightconeGrid& g, const std::int64_t cell, const LightconeFace& f)
  {
    if (f.coordinate < 0) return false;
    const std::size_t ring = ring_of(g, cell), shell = shell_of(g, ring);
    const Point& p = f.point;
    const double r = norm(p), dec = std::asin(p[2] / r), ra = std::atan2(p[1], p[0]);
    const double rLo = g.shellEdge[shell], rHi = g.shellEdge[shell + 1];
    const double decLo = g.ringDecMin[ring], decHi = g.ringDecMax[ring];
    const double raLo = g.cells.cellRaMin[(std::size_t)cell], raHi = g.cells.cellRaMax[(std::size_t)cell];
    const double tol = 1e-11;
    bool ok = r >= rLo - tol * rHi && r <= rHi + tol * rHi && dec >= decLo - tol && dec <= decHi + tol;
    if (f.coordinate == 0) ok = ok && std::fabs(r - (f.side < 0 ? rLo : rHi)) <= tol * rHi;
    if (f.coordinate == 1) ok = ok && std::fabs(dec - (f.side < 0 ? decLo : decHi)) <= tol;
    if (f.coordinate == 2) ok = ok && angle_gap(ra, f.side < 0 ? raLo : raHi) <= tol;
    return ok;
  }

  Point expected_normal (const LightconeGrid& g, const std::int64_t cell, const LightconeFace& f)
  {
    const std::size_t ring = ring_of(g, cell);
    const Point& p = f.point;
    const double s = f.side;
    if (f.coordinate == 0) {
      const double r = norm(p);
      return {s * p[0] / r, s * p[1] / r, s * p[2] / r};
    }
    if (f.coordinate == 1) {
      const double e = s < 0 ? g.ringDecMin[ring] : g.ringDecMax[ring];
      const double az = std::atan2(p[1], p[0]);
      return {-s * std::sin(e) * std::cos(az), -s * std::sin(e) * std::sin(az), s * std::cos(e)};
    }
    const double e = s < 0 ? g.cells.cellRaMin[(std::size_t)cell] : g.cells.cellRaMax[(std::size_t)cell];
    return {-s * std::sin(e), s * std::cos(e), 0.};
  }

  void check_face (const LightconeGrid& g, const std::int64_t cell, const LightconeFace& f, const int coordinate,
                   const int side, const std::string& what)
  {
    check(is_face(f, coordinate, side), what + ": face " + std::to_string(f.coordinate) + " side " +
                                        std::to_string(f.side));
    check(f.t > 0. && f.t < 1., what + ": crossing inside the segment, t = " + std::to_string(f.t));
    check(on_face(g, cell, f), what + ": the crossing lies on the face");
    check(close(f.normal, expected_normal(g, cell, f), 1e-12), what + ": the outward normal");
  }

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  std::vector<double> objects (const std::size_t n, const unsigned seed)
  {
    std::mt19937 rng(seed);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < n; ++i) {
      sky.push_back(uniform(0.5, 2.));
      sky.push_back(otswap::maths::asin(uniform(otswap::maths::sin(-0.3), otswap::maths::sin(0.3))));
      sky.push_back(uniform(0.25, 0.4));
    }
    return sky;
  }

  const LightconeGrid& built ()
  {
    static const LightconeGrid g = [] {
      const std::vector<double> sky = objects(4000, 1);
      const internal::LightconeExtent e = internal::lightcone_extent(sky, {}, {}, {}, table());
      DivergenceGrid record;
      record.geometry = otswap::Geometry::Lightcone;
      record.distanceMin = e.distanceMin;
      record.distanceMax = e.distanceMax;
      record.raMin = e.raMin;
      record.raMax = e.raMax;
      record.decMin = e.decMin;
      record.decMax = e.decMax;
      record.cellSizeMps = 2.;
      record.mpsProfile.redshift = {0.1, 0.5};
      record.mpsProfile.mps = {5., 6.};
      const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
      return internal::lightcone_grid(record, table(), mask);
    }();
    return g;
  }

  struct Segment {
    Point x, d;
    std::int64_t startCell, endCell;
  };

  std::vector<Segment> segments (const std::size_t n, const double length, const unsigned seed)
  {
    const std::vector<double> sky = objects(n, seed);
    const std::vector<double> xyz = otswap::toCartesian(sky, table());
    std::mt19937 rng(seed + 100);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<Segment> out;
    for (std::size_t i = 0; i < n; ++i) {
      Segment s;
      s.x = {xyz[3*i], xyz[3*i+1], xyz[3*i+2]};
      s.d = {uniform(-length, length), uniform(-length, length), uniform(-length, length)};
      s.startCell = internal::lightcone_cell(built(), table().distanceAt(sky[3*i+2]), sky[3*i], sky[3*i+1]);
      const Point e = end_of(s.x, s.d);
      s.endCell = internal::lightcone_cell_at(built(), e.data());
      if (s.startCell >= 0 && s.endCell >= 0) out.push_back(s);
    }
    return out;
  }

  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  std::uint64_t faces_hash ()
  {
    const std::vector<Segment> s = segments(3000, 20., 7);
    const std::int64_t n = (std::int64_t)s.size();
    std::vector<LightconeFace> exits((std::size_t)n), entries((std::size_t)n);
    std::vector<std::int64_t> cells((std::size_t)n);
#pragma omp parallel for schedule(dynamic, 7)
    for (std::int64_t i = 0; i < n; ++i) {
      const Segment& g = s[(std::size_t)i];
      cells[(std::size_t)i] = internal::lightcone_cell_at(built(), g.x.data());
      exits[(std::size_t)i] = internal::lightcone_face(built(), g.startCell, g.x.data(), g.d.data(), true);
      entries[(std::size_t)i] = internal::lightcone_face(built(), g.endCell, g.x.data(), g.d.data(), false);
    }
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (std::size_t i = 0; i < (std::size_t)n; ++i)
      for (const LightconeFace* f : {&exits[i], &entries[i]}) {
        hash_bytes(h, &f->coordinate, sizeof f->coordinate);
        hash_bytes(h, &f->side, sizeof f->side);
        hash_bytes(h, &f->t, sizeof f->t);
        hash_bytes(h, f->point.data(), sizeof f->point);
        hash_bytes(h, f->normal.data(), sizeof f->normal);
      }
    hash_bytes(h, cells.data(), cells.size() * sizeof(std::int64_t));
    return h;
  }

}

int main ()
{
  const LightconeGrid& t = toy();
  const LightconeGrid& p = polar();

  // A Cartesian point is found in the cell of its sky coordinates; a point on an edge belongs to the cell above it, the north pole to its cap, RA = 0 to the first sector of a full ring, the origin to the central sphere, and a point beyond the grid to no cell.
  group("a Cartesian point is in the cell of its sky coordinates, edges belonging to the cell above");
  {
    const std::vector<double> sky = objects(2000, 3);
    const std::vector<double> xyz = otswap::toCartesian(sky, table());
    const std::vector<double> back = otswap::toSky(xyz, table());
    bool same = true, sameSky = true;
    std::size_t inside = 0;
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      double r = 0., ra = 0., dec = 0.;
      internal::lightcone_sky(&xyz[3*i], r, ra, dec);
      sameSky = sameSky && ra == back[3*i] && dec == back[3*i+1] && r == norm({xyz[3*i], xyz[3*i+1], xyz[3*i+2]});
      const std::int64_t c = internal::lightcone_cell_at(built(), &xyz[3*i]);
      same = same && c == internal::lightcone_cell(built(), r, back[3*i], back[3*i+1]);
      inside += c >= 0;
    }
    check(sameSky, "the sky coordinates are otswap::toSky's, bit for bit");
    check(same && inside > 1900, "the cell is that of the sky coordinates");
    const Point vertex {0., 110., 0.};
    check(internal::lightcone_cell_at(t, vertex.data()) == toy_cell(1, 2, 1),
          "a point on a distance, declination and right-ascension edge belongs to the cell above in all three");
    const Point pole {0., 0., 55.};
    check(internal::lightcone_cell_at(p, pole.data()) == 10, "the north pole belongs to its cap");
    const Point seam {55., 0., 0.};
    check(internal::lightcone_cell_at(p, seam.data()) == 6, "RA = 0 belongs to the first sector of a full ring");
    const Point below {55., -1e-9, 0.};
    check(internal::lightcone_cell_at(p, below.data()) == 9, "just below RA = 0 belongs to the last sector");
    const Point origin {0., 0., 0.};
    check(internal::lightcone_cell_at(p, origin.data()) == 0, "the origin belongs to the central sphere");
    const Point far {0., 0., 100.};
    check(internal::lightcone_cell_at(p, far.data()) == -1, "a point beyond the last edge is in no cell");
    const Point aside = sky_point(105., 0.1, 2.);
    check(internal::lightcone_cell_at(t, aside.data()) == -1, "a point beside the sectors of a partial ring is in no cell");
  }

  // A point is inside a cell exactly when the lookup returns that cell, and otherwise each coordinate reports the side of the cell it lies on, a point just below RA = 0 being below the first sector of a full ring.
  group("the sides of a point agree with the lookup, coordinate by coordinate");
  {
    std::mt19937 rng(5);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    bool agree = true;
    for (const LightconeGrid* g : {&t, &p})
      for (int i = 0; i < 400; ++i) {
        const double r = g == &t ? uniform(95., 125.) : uniform(0., 75.);
        const Point x = sky_point(r, uniform(-1.5, 1.5), uniform(-1., 7.));
        double d = 0., ra = 0., dec = 0.;
        internal::lightcone_sky(x.data(), d, ra, dec);
        const std::int64_t found = internal::lightcone_cell_at(*g, x.data());
        for (std::int64_t c = 0; c < (std::int64_t)g->cells.nCells; ++c) {
          const internal::LightconeOutside o = internal::lightcone_outside(*g, c, d, ra, dec);
          agree = agree && ((o.side[0] == 0 && o.side[1] == 0 && o.side[2] == 0) == (found == c));
        }
      }
    for (int i = 0; i < 2000; ++i) {
      const Point x = sky_point(uniform(700., 1100.), uniform(-0.3, 0.3), uniform(0.5, 2.));
      double d = 0., ra = 0., dec = 0.;
      internal::lightcone_sky(x.data(), d, ra, dec);
      const std::int64_t found = internal::lightcone_cell_at(built(), x.data());
      if (found < 0) continue;
      for (std::int64_t c = std::max<std::int64_t>(0, found - 2);
           c <= std::min<std::int64_t>((std::int64_t)built().cells.nCells - 1, found + 2); ++c) {
        const internal::LightconeOutside o = internal::lightcone_outside(built(), c, d, ra, dec);
        agree = agree && ((o.side[0] == 0 && o.side[1] == 0 && o.side[2] == 0) == (found == c));
      }
    }
    check(agree, "all sides are 0 exactly in the cell the lookup returns");
    const internal::LightconeOutside below = internal::lightcone_outside(p, 6, 55., kTwoPi - 1e-9, 0.5);
    check(below.side[0] == 0 && below.side[1] == 0 && below.side[2] == -1,
          "just below RA = 0 is below the first sector of a full ring");
    const internal::LightconeOutside corner = internal::lightcone_outside(t, toy_cell(0, 2, 0), 115., 1.9, -0.1);
    check(corner.side[0] == 1 && corner.side[1] == -1 && corner.side[2] == 1, "each coordinate reports its own side");
    const internal::LightconeOutside low = internal::lightcone_outside(t, toy_cell(1, 3, 1), 105., 1.0, 0.4);
    check(low.side[0] == -1 && low.side[1] == 1 && low.side[2] == -1, "below in distance and right ascension");
  }

  // A segment leaving a cell through one face gives that face, the crossing on it and the outward normal.
  group("a segment leaving through one face gives that face, its crossing and its outward normal");
  {
    const std::int64_t north = toy_cell(0, 3, 0), south = toy_cell(0, 0, 0), equator = toy_cell(0, 2, 0);
    const Point n = sky_point(105., 0.225, 1.43), s = sky_point(105., -0.225, 1.43), e = sky_point(105., 0.075, 1.43);
    check_face(t, north, leave(t, north, n, minus(sky_point(115., 0.225, 1.43), n)), 0, 1, "outer sphere");
    check_face(t, north, leave(t, north, n, minus(sky_point(95., 0.225, 1.43), n)), 0, -1, "inner sphere");
    check_face(t, north, leave(t, north, n, minus(sky_point(105., 0.4, 1.43), n)), 1, 1, "upper cone, north");
    check_face(t, north, leave(t, north, n, minus(sky_point(105., 0.05, 1.43), n)), 1, -1, "lower cone, north");
    check_face(t, south, leave(t, south, s, minus(sky_point(105., -0.1, 1.43), s)), 1, 1, "upper cone, south");
    check_face(t, south, leave(t, south, s, minus(sky_point(105., -0.4, 1.43), s)), 1, -1, "lower cone, south");
    const LightconeFace plane = leave(t, equator, e, minus(sky_point(105., -0.05, 1.43), e));
    check_face(t, equator, plane, 1, -1, "plane of a declination edge at 0");
    check(plane.normal == Point{0., 0., -1.}, "the plane's normal is exactly -z");
    check_face(t, north, leave(t, north, n, minus(sky_point(105., 0.225, 1.2), n)), 2, -1, "lower right ascension");
    check_face(t, north, leave(t, north, n, minus(sky_point(105., 0.225, 1.6), n)), 2, 1, "upper right ascension");
  }

  // Two adjacent cells crossed directly share the crossing, with opposite normals.
  group("adjacent cells crossed directly share the crossing, with opposite normals");
  {
    struct Pair { const LightconeGrid* g; std::int64_t from; Point x, e; std::string what; };
    const std::vector<Pair> pairs {
      {&t, toy_cell(0, 3, 0), sky_point(105., 0.225, 1.43), sky_point(113., 0.24, 1.45), "across a sphere"},
      {&t, toy_cell(0, 2, 0), sky_point(105., 0.075, 1.43), sky_point(106., 0.2, 1.44), "across a cone"},
      {&t, toy_cell(0, 2, 0), sky_point(105., 0.075, 1.43), sky_point(104., -0.05, 1.42), "across the plane z = 0"},
      {&t, toy_cell(0, 3, 0), sky_point(105., 0.225, 1.43), sky_point(105., 0.23, 1.62), "across a right-ascension plane"},
      {&p, 6, sky_point(55., 0.5, 0.3), sky_point(55., 0.5, -0.3), "across RA = 0 in a full ring"}};
    for (const Pair& q : pairs) {
      const Point d = minus(q.e, q.x);
      const LightconeFace out = leave(*q.g, q.from, q.x, d), in = enter(*q.g, q.x, d);
      check(internal::lightcone_cell_at(*q.g, q.x.data()) == q.from &&
            internal::lightcone_cell_at(*q.g, q.e.data()) != q.from, q.what + ": the segment changes cell");
      check(out.coordinate >= 0 && out.coordinate == in.coordinate && out.side == -in.side, q.what + ": the same face");
      check(out.t == in.t && out.point == in.point, q.what + ": the same crossing, bit for bit");
      check(out.normal == Point{-in.normal[0], -in.normal[1], -in.normal[2]}, q.what + ": opposite normals");
    }
  }

  // A segment that leaves its cell and comes back is found leaving at its first exit and entering at its last entry.
  group("a segment that leaves and re-enters its cell: the first exit and the last entry");
  {
    const std::int64_t c = toy_cell(0, 2, 0);
    const Point x = sky_point(100.2, 0.07, 1.33), e = sky_point(100.2, 0.07, 1.53);
    const Point d = minus(e, x);
    check(internal::lightcone_cell_at(t, x.data()) == c && internal::lightcone_cell_at(t, e.data()) == c,
          "both ends in the same cell");
    const LightconeFace out = leave(t, c, x, d), in = enter(t, x, d);
    check_face(t, c, out, 0, -1, "first exit");
    check_face(t, c, in, 0, -1, "last entry");
    check(out.t < 0.5 && in.t > 0.5, "the exit comes before the dip, the entry after it");
  }

  // Crossing the other nappe of a declination cone or the other half of a right-ascension plane is not taken as an exit.
  group("the other nappe of a cone and the other half of a right-ascension plane are not faces");
  {
    const std::int64_t north = toy_cell(0, 3, 0);
    const Point n = sky_point(105., 0.225, 1.43), s = sky_point(105., -0.225, 1.43);
    const LightconeFace out = leave(t, north, n, minus(s, n));
    check_face(t, north, out, 1, -1, "leaving toward the other nappe");
    const LightconeFace in = enter(t, s, minus(n, s));
    check_face(t, north, in, 1, -1, "entering after the other nappe");
    const Point far = sky_point(105., 0.225, 1.43 + kPi + 0.3);
    const LightconeFace across = enter(t, far, minus(n, far));
    check(across.coordinate >= 0 && on_face(t, north, across), "entering after the other half of a plane");
    check(across.t > 0.5, "the entry is near the cell, not at the other half");
    const LightconeFace back = leave(t, north, n, minus(far, n));
    check(back.coordinate >= 0 && on_face(t, north, back) && back.t < 0.5, "leaving toward the other half");
  }

  // A start or an end lying on a face gives the crossing at t = 0 or t = 1, also when the root rounds just outside [0, 1].
  group("a start or an end on a face gives the crossing at t = 0 or t = 1");
  {
    const std::int64_t equator = toy_cell(0, 2, 0);
    const Point x {105. * std::cos(1.43), 105. * std::sin(1.43), 0.};
    check(internal::lightcone_cell_at(t, x.data()) == equator, "a start on the plane z = 0 is in the ring above");
    const LightconeFace down = leave(t, equator, x, {0.1, 0.1, -5.});
    check(is_face(down, 1, -1) && down.t == 0. && down.point == x, "leaving downward at t = 0");
    const LightconeFace up = leave(t, equator, x, {0.1, 0.1, 30.});
    check(up.coordinate >= 0 && up.t > 0., "leaving upward later");
    const Point below {x[0], x[1], -3.};
    const Point d {0.5, -0.2, 3.};
    check(end_of(below, d)[2] == 0., "an end exactly on the plane");
    const LightconeFace onPlane = enter(t, below, d);
    check(is_face(onPlane, 1, -1) && onPlane.t == 1. && onPlane.point == end_of(below, d), "entering at t = 1");
    const Point sphere {20., 108., 6.};
    const std::int64_t outer = toy_cell(1, 2, 0);
    check(internal::lightcone_cell_at(t, sphere.data()) == outer, "a point on the sphere is in the shell above");
    bool ends = true, starts = true;
    for (const Point& dir : {Point{1., 5., 0.25}, Point{0.5, 4., -0.125}, Point{2., 6., 0.5},
                             Point{-0.5, 3., -0.25}, Point{0.125, 1., 0.0625}}) {
      const LightconeFace in = enter(t, minus(sphere, dir), dir);
      ends = ends && is_face(in, 0, -1) && in.t >= 1. - 1e-13 && in.t <= 1.;
      const LightconeFace out = leave(t, outer, sphere, {-dir[0], -dir[1], -dir[2]});
      starts = starts && is_face(out, 0, -1) && out.t >= 0. && out.t <= 1e-13;
    }
    check(ends, "an end exactly on a sphere enters at t within rounding of 1");
    check(starts, "a start on a sphere leaves at t within rounding of 0");
  }

  // A tracer whose Cartesian position lies just outside the cell of its sky coordinates leaves that cell at t = 0.
  group("a start just outside its sky cell leaves it at t = 0 through the face it violates");
  {
    const std::int64_t below = toy_cell(0, 1, 0);
    const Point x {105. * std::cos(1.43), 105. * std::sin(1.43), 0.};
    for (const Point& d : {Point{0.1, 0.1, -2.}, Point{0.1, 0.1, 2.}, Point{0., 0., -0.5}}) {
      const LightconeFace f = leave(t, below, x, d);
      check(is_face(f, 1, 1) && f.t == 0. && f.point == x && f.normal == Point{0., 0., 1.},
            "the start on the face above leaves at t = 0, whatever the direction");
    }
    std::size_t moved = 0;
    bool atZero = true;
    for (const Segment& s : segments(4000, 15., 11)) {
      if (internal::lightcone_cell_at(built(), s.x.data()) == s.startCell) continue;
      ++moved;
      const LightconeFace f = internal::lightcone_face(built(), s.startCell, s.x.data(), s.d.data(), true);
      atZero = atZero && f.t == 0. && f.coordinate >= 0;
    }
    check(atZero, "every tracer outside its sky cell (" + std::to_string(moved) + " found) leaves at t = 0");
  }

  // Through an edge or a vertex the face follows the fixed order of coordinates and sides.
  group("a point violating several sides gives the face of the fixed order");
  {
    struct Tie { std::int64_t cell; Point x; int coordinate, side; std::string what; };
    const std::vector<Tie> ties {
      {5, {60., 0., 0.}, 0, 1, "distance above before declination and right ascension above"},
      {2, {60., 0., 0.}, 0, 1, "distance above before declination above"},
      {15, {60., 0., 0.}, 1, 1, "declination above before right ascension above"},
      {9, {60., 0., 0.}, 0, 1, "distance above before right ascension above"},
      {16, {59.9, -0.1, -0.1}, 2, -1, "right ascension below before declination and distance below"},
      {19, {59.9, -0.1, -0.1}, 1, -1, "declination below before distance below"},
      {6, {60.5, 0., -0.1}, 0, 1, "distance above before declination below"},
      {9, {55., 0.1, -0.1}, 2, 1, "right ascension above before declination below"},
      {2, {55., -0.1, 0.1}, 1, 1, "declination above before right ascension below"},
      {6, {49., 0.1, -0.1}, 1, -1, "declination below before distance below"}};
    for (const Tie& q : ties) {
      const LightconeFace f = leave(p, q.cell, q.x, {0.01, 0.02, 0.03});
      check(is_face(f, q.coordinate, q.side) && f.t == 0., q.what);
    }
  }

  // Faces without area are never reported: the inner surface of the central sphere, the pole of a cap, the right-ascension planes of a ring with a single sector over the full circle.
  group("faces without area are never reported");
  {
    std::mt19937 rng(9);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    bool sphere = true;
    for (int i = 0; i < 200; ++i) {
      const Point x = sky_point(uniform(0., 40.), uniform(-1.5, 1.5), uniform(0., kTwoPi));
      const Point e = sky_point(uniform(52., 58.), uniform(-1.5, 1.5), uniform(0., kTwoPi));
      sphere = sphere && is_face(leave(p, 0, x, minus(e, x)), 0, 1);
    }
    check(sphere, "the central sphere is left only through its outer surface");
    check(leave(p, 0, {10., 0., 0.}, {-20., 0., 0.}).coordinate == -1, "through the observer, nothing");
    check(leave(p, 0, {10., -5., 1.}, {0., 10., 0.}).coordinate == -1, "across RA = 0 in a single sector, nothing");
    const Point a = sky_point(55., 1.4, 0.5);
    check(leave(p, 10, a, minus(sky_point(55., 1.4, 0.5 + kPi), a)).coordinate == -1, "across the pole, nothing");
    const LightconeFace cap = leave(p, 10, a, minus(sky_point(55., 1.0, 0.5 + kPi), a));
    check_face(p, 10, cap, 1, -1, "a cap is left through its cone only");
  }

  // A displacement that ends on the lower declination face of its cell close to the equator, from a tracer two rings below, enters through that face at t in [1 - 1e-15, 1], with the cone's outward normal.
  group("a displacement ending just above a declination edge near the equator enters through it");
  {
    const double step = 0.006764048031421698, edge = 3.4906585037394744e-08;
    const std::vector<double> sector {6.3172599080748038, 6.3239145784745503};
    const LightconeGrid g = hand_grid(false, {1806.3833196237779, 1818.4826654701524},
                                      {{{edge - 2. * step, edge - step, sector, false, false},
                                        {edge - step, edge, sector, false, false},
                                        {edge, 0.006764082938006735, sector, false, false}}});
    const Point x {1814.0627046203342, 72.632939349913272, -17.099523433542164};
    const Point d {-0.75758173299800546, -0.54578232868746568, 17.099586779829426};
    const Point e = end_of(x, d);
    check(internal::lightcone_cell_at(g, e.data()) == 2, "the end lies in the cell above the edge");
    const LightconeFace in = enter(g, x, d);
    check(is_face(in, 1, -1) && in.t >= 1. - 1e-15 && in.t <= 1., "the entry is the lower declination face at t ~ 1");
    check(close(in.normal, expected_normal(g, 2, in), 1e-12) && in.normal[2] < -0.999, "the cone's outward normal");
    const std::int64_t start = internal::lightcone_cell_at(g, x.data());
    check(start == 0, "the tracer lies two rings below");
    check(leave(g, start, x, d).coordinate >= 0, "the tracer's cell gives a face");
  }

  // No face is reported only when every point examined lies in the cell.
  group("no face only when the segment never leaves the cell");
  {
    const std::int64_t c = toy_cell(0, 3, 0);
    const Point x = sky_point(105., 0.225, 1.43);
    const Point d = minus(sky_point(106., 0.23, 1.44), x);
    const LightconeFace out = leave(t, c, x, d), in = internal::lightcone_face(t, c, x.data(), d.data(), false);
    check(out.coordinate == -1 && out.t == 1. && in.coordinate == -1 && in.t == 0., "a segment inside its cell");
    bool always = true;
    std::size_t crossing = 0;
    for (const Segment& s : segments(4000, 15., 13)) {
      if (s.startCell == s.endCell) continue;
      ++crossing;
      always = always && internal::lightcone_face(built(), s.startCell, s.x.data(), s.d.data(), true).coordinate >= 0;
    }
    check(crossing > 1000 && always, "a segment ending in another cell always leaves through a face");
  }

  // Faces and lookups of a fixed set of segments on a built grid are identical bit for bit with any number of threads and match the recorded hash.
  group("the faces and lookups are the same with 1 and 4 threads, and match the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = faces_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = faces_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "faces hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_faces");
}
