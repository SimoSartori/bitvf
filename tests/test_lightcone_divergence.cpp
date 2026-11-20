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
#include <limits>
#include <random>
#include <sstream>
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
using internal::LightconeDivergence;
using internal::LightconeFace;
using internal::LightconeGrid;

namespace {

  constexpr std::uint64_t kRecordedHash = 0x84aa9dcad44025c8ull;

  using Point = std::array<double, 3>;

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

  LightconeGrid make_grid ()
  {
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
  }

  struct Sweep {

    LightconeGrid grid = make_grid();
    std::size_t N = 0;
    unsigned R = 1;
    std::vector<double> tracers, sky, displacement;
    std::vector<std::uint8_t> valid, inMask, inCut;
    bool reused = false;

    std::size_t add (const double ra, const double dec, const double z, const Point& d)
    {
      const std::vector<double> s {ra, dec, z};
      const std::vector<double> x = otswap::toCartesian(s, table());
      return add_at({x[0], x[1], x[2]}, s, d);
    }

    std::size_t add_centre (const Point& d)
    {
      const std::size_t c = grid.cells.nCells / 2;
      const std::vector<double>& centre = grid.cells.centreSky;
      return add(centre[3*c], centre[3*c+1], centre[3*c+2], d);
    }

    std::size_t add_at (const Point& x, const std::vector<double>& s, const Point& d)
    {
      tracers.insert(tracers.end(), x.begin(), x.end());
      sky.insert(sky.end(), s.begin(), s.end());
      displacement.insert(displacement.end(), d.begin(), d.end());
      inMask.push_back(1);
      inCut.push_back(1);
      return N++;
    }

    LightconeDivergence run ()
    {
      internal::LightconeSweep s;
      s.nTracers = N;
      s.nRealizations = R;
      s.tracers = &tracers;
      s.tracersSky = &sky;
      s.inMask = &inMask;
      s.inCut = &inCut;
      s.displacement = &displacement;
      s.valid = &valid;
      s.distances = &table();
      s.reusedGrid = reused;
      return internal::lightcone_divergence(s, grid);
    }

    Point x (const std::size_t i) const
    {
      return {tracers[3*i], tracers[3*i+1], tracers[3*i+2]};
    }

    std::int64_t start_cell (const std::size_t i) const
    {
      return internal::lightcone_cell(grid, table().distanceAt(sky[3*i+2]), sky[3*i], sky[3*i+1]);
    }

    std::int64_t end_cell (const std::size_t i) const
    {
      const double e[3] = {tracers[3*i] + displacement[3*i], tracers[3*i+1] + displacement[3*i+1],
                           tracers[3*i+2] + displacement[3*i+2]};
      return internal::lightcone_cell_at(grid, e);
    }
  };

  std::uint64_t total (const std::vector<std::uint64_t>& v)
  {
    std::uint64_t n = 0;
    for (const std::uint64_t k : v) n += k;
    return n;
  }

  Point radial (const double ra, const double dec, const double length)
  {
    const double c = otswap::maths::cos(dec);
    return {length * c * otswap::maths::cos(ra), length * c * otswap::maths::sin(ra), length * otswap::maths::sin(dec)};
  }

  Point to_sky_point (const double ra, const double dec, const double z, const Point& from)
  {
    const std::vector<double> e = otswap::toCartesian({ra, dec, z}, table());
    return {e[0] - from[0], e[1] - from[1], e[2] - from[2]};
  }

  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  std::uint64_t divergence_hash ()
  {
    Sweep s;
    s.R = 3;
    const std::vector<double> sky = objects(3000, 11);
    std::mt19937 rng(5);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    const std::size_t n = sky.size() / 3;
    const std::vector<double> x = otswap::toCartesian(sky, table());
    s.tracers = x;
    s.sky = sky;
    s.N = n;
    s.inMask.assign(n, 1);
    s.inCut.assign(n, 1);
    s.displacement.resize(3 * n * s.R);
    for (double& d : s.displacement) d = uniform(-15., 15.);
    s.valid.assign(n * s.R, 1);
    for (std::size_t k = 0; k < s.valid.size(); k += 17) s.valid[k] = 0;
    for (std::size_t c = 0; c < s.grid.cells.nCells; c += 23) s.grid.cells.masked[c] = 1;
    const LightconeDivergence div = s.run();
    std::uint64_t h = 0xcbf29ce484222325ull;
    const DivergenceGrid& g = s.grid.cells;
    hash_bytes(h, g.divergence.data(), g.divergence.size() * sizeof(double));
    hash_bytes(h, g.weight.data(), g.weight.size() * sizeof(double));
    hash_bytes(h, g.nDisplacements.data(), g.nDisplacements.size() * sizeof(std::uint64_t));
    hash_bytes(h, g.nTracers.data(), g.nTracers.size() * sizeof(std::uint64_t));
    hash_bytes(h, div.cell.data(), div.cell.size() * sizeof(std::int64_t));
    hash_bytes(h, div.used.data(), div.used.size());
    hash_bytes(h, &div.nRejected, sizeof div.nRejected);
    hash_bytes(h, &div.scaleExponent, sizeof div.scaleExponent);
    return h;
  }

  struct Case {

    std::vector<double> tracersSky {1., 0.2, 0.3,  2., -0.2, 0.35,  3., 0.1, 0.4,  4., 0.3, 0.25};
    std::vector<double> displacement = std::vector<double>(12, 0.1);
    std::vector<double> randomsSky {1.5, 0.1, 0.3,  2.5, 0., 0.3,  3.5, 0.2, 0.35,  4.5, -0.1, 0.4,
                                    5., 0.2, 0.3,  6., 0.1, 0.25,  0.5, 0.3, 0.3,  0.2, -0.3, 0.35};
    otswap::Mask mask {std::vector<double>(12, 1.), otswap::PixelOrdering::Ring};
    otswap::MpsProfile profile;
    Config config;

    Case ()
    {
      profile.redshift = {0.1, 0.5};
      profile.mps = {5., 6.};
      config.smoothing.sigma = 0.;
      config.verbosity = otswap::Verbosity::Silent;
    }

    std::string run (const otswap::Verbosity verbosity, const std::string& text)
    {
      std::ostringstream log;
      config.verbosity = verbosity;
      config.log = &log;
      throws_naming([&] { findVoidsLightcone(tracersSky, displacement, {}, randomsSky, mask, table(), profile, config); },
                    text, "run with a log");
      return log.str();
    }
  };

}

int main ()
{
  // A displacement whose end lies outside the grid is rejected and counts in no cell, while its tracer is still counted.
  group("an end outside the grid is rejected");
  {
    Sweep s;
    const std::size_t kept = s.add_centre({0.1, 0.1, 0.1});
    const std::size_t out = s.add(1.0, 0.0, 0.32, radial(1.0, 0.0, 2000.));
    check(s.end_cell(out) < 0 && s.end_cell(kept) == s.start_cell(kept), "the ends lie where intended ");
    const LightconeDivergence div = s.run();
    check(div.nRejected == 1 && div.used[out] == 0 && div.used[kept] == 1 && div.nUsed == 1, "one rejected");
    check(total(s.grid.cells.nDisplacements) == 2 && s.grid.cells.nDisplacements[(std::size_t)div.cell[kept]] == 2,
          "only the accepted displacement counts, twice in its cell");
    check(total(s.grid.cells.nTracers) == 2 && s.grid.cells.nTracers[(std::size_t)div.cell[out]] == 1,
          "both tracers are counted in their cells");
  }

  // A displacement that ends in a masked cell, or starts in one, is rejected; the tracer in a masked cell is still counted.
  group("an end or a start in a masked cell is rejected");
  {
    Sweep s;
    const std::size_t kept = s.add_centre({0.1, 0.1, 0.1});
    const std::size_t toMasked = s.add(1.5, 0.0, 0.32, radial(1.5, 0.0, 40.));
    const std::size_t fromMasked = s.add(0.8, -0.1, 0.3, radial(0.8, -0.1, 40.));
    const std::int64_t endCell = s.end_cell(toMasked), startCell = s.start_cell(fromMasked);
    check(endCell >= 0 && endCell != s.start_cell(toMasked) && s.end_cell(fromMasked) != startCell,
          "both displacements cross into another cell");
    s.grid.cells.masked[(std::size_t)endCell] = 1;
    s.grid.cells.masked[(std::size_t)startCell] = 1;
    const LightconeDivergence div = s.run();
    check(div.nRejected == 2 && div.used[toMasked] == 0 && div.used[fromMasked] == 0 && div.used[kept] == 1,
          "both rejected");
    check(s.grid.cells.nDisplacements[(std::size_t)endCell] == 0 && s.grid.cells.nDisplacements[(std::size_t)startCell] == 0 &&
          total(s.grid.cells.nDisplacements) == 2, "the masked cells receive nothing");
    check(s.grid.cells.nTracers[(std::size_t)startCell] == 1, "the tracer in a masked cell is counted");
  }

  // Without a central sphere, a path passing closer to the observer than the grid's first edge is rejected; one whose closest point is an end, or that stays above the edge, is not.
  group("a path below the first edge is rejected; an inward step and a path above the edge are not");
  {
    Sweep s;
    check(!s.grid.centralSphere, "the grid has no central sphere");
    const double z = 0.2505;
    const std::size_t chord = s.add(0.6, 0.0, z, {0., 0., 0.});
    const Point d = to_sky_point(1.9, 0.0, z, s.x(chord));
    for (std::size_t k = 0; k < 3; ++k) s.displacement[3*chord+k] = d[k];
    const std::size_t inward = s.add(1.2, 0.1, z, radial(1.2, 0.1, -1.));
    const std::size_t shortChord = s.add(1.0, 0.0, 0.39, {0., 0., 0.});
    const Point e = to_sky_point(1.1, 0.0, 0.39, s.x(shortChord));
    for (std::size_t k = 0; k < 3; ++k) s.displacement[3*shortChord+k] = e[k];
    check(s.end_cell(chord) >= 0 && s.end_cell(inward) >= 0 && s.end_cell(shortChord) >= 0, "every end is in the grid");
    const LightconeDivergence div = s.run();
    check(div.nRejected == 1 && div.used[chord] == 0, "the long chord is rejected");
    check(div.used[inward] == 1 && div.used[shortChord] == 1, "the inward step and the short chord are accepted");
  }

  // A reused grid leaves out the kept tracers outside it, as masked and cut rows: their displacements are not read, not even for finiteness or the largest norm.
  group("tracers outside a reused grid are not read; those inside must have finite displacements");
  {
    Sweep s;
    s.reused = true;
    const std::size_t kept = s.add(1.2, 0.05, 0.32, {3., 4., 0.});
    const std::size_t outside = s.add(3.0, 0.0, 0.32, {1e6, 0., 0.});
    const std::size_t nan = s.add(3.1, 0.0, 0.32, {std::numeric_limits<double>::quiet_NaN(), 0., 0.});
    const LightconeDivergence div = s.run();
    check(div.outsideGrid.size() == 3 && div.outsideGrid[outside] == 1 && div.outsideGrid[nan] == 1 &&
          div.outsideGrid[kept] == 0 && div.nOutsideGrid == 2, "two outside");
    check(div.cell[outside] == -1 && div.cell[nan] == -1 && total(s.grid.cells.nTracers) == 1, "only one counted");
    check(div.maxDisplacement == 5., "the largest norm is that of the tracer read");
    Sweep bad = s;
    bad.displacement[3*kept] = std::numeric_limits<double>::quiet_NaN();
    throws_naming([&] { bad.run(); }, "the displacements hold a non-finite value for realization 0 (from 0), object 0 "
                  "(from 0), component 0", "a non-finite displacement of a tracer read is refused");
    Case c;
    c.randomsSky = {1.5, 0.1, 0.3,  2.5, 0., 0.3,  0.5, 0.3, 0.3,  1.2, -0.3, 0.35};
    const std::vector<double> firstTwo(c.tracersSky.begin(), c.tracersSky.begin() + 6);
    const DivergenceGrid g = lightconeGrid(firstTwo, c.randomsSky, c.mask, table(), c.profile, 4.);
    c.displacement[9] = std::numeric_limits<double>::quiet_NaN();
    throws_naming([&] { findVoidsLightcone(c.tracersSky, c.displacement, {}, c.randomsSky, c.mask, table(), c.profile,
                                           c.config, {}, {}, &g); },
                  "not implemented yet", "a non-finite displacement outside the reused grid is not read");
    c.displacement[0] = std::numeric_limits<double>::quiet_NaN();
    throws_naming([&] { findVoidsLightcone(c.tracersSky, c.displacement, {}, c.randomsSky, c.mask, table(), c.profile,
                                           c.config, {}, {}, &g); },
                  "non-finite value for realization 0 (from 0), object 0 (from 0)",
                  "a non-finite displacement inside the reused grid is refused");
    Sweep fresh;
    fresh.add(1.2, 0.05, 0.32, {3., 4., 0.});
    check(fresh.run().outsideGrid.empty(), "a new grid gives no outside flags");
  }

  // Tracers not kept, or displacements not valid, are not read.
  group("tracers not kept and invalid displacements are not read");
  {
    Sweep s;
    const std::size_t kept = s.add_centre({0.1, 0.1, 0.1});
    const std::size_t masked = s.add(1.3, 0.05, 0.32, radial(1.3, 0.05, 2000.));
    const std::size_t invalid = s.add(1.4, 0.05, 0.32, radial(1.4, 0.05, 2000.));
    s.inMask[masked] = 0;
    s.valid = {1, 1, 0};
    const LightconeDivergence div = s.run();
    check(div.nRejected == 0 && div.cell[masked] == -1 && div.cell[invalid] >= 0 && div.used[invalid] == 0 &&
          div.used[kept] == 1 && total(s.grid.cells.nTracers) == 2, "neither is read nor rejected");
  }

  // A tracer whose position lies within rounding of a face, inside the cell across it from its sky position, enters the end cell at its own position through that face, and the displacement contributes normally.
  group("a tracer on a face enters the end cell at t = 0 through the face its sky position lies beyond");
  {
    Sweep s;
    const double dec = 0.02, z = 0.33;
    const std::int64_t below = internal::lightcone_cell(s.grid, table().distanceAt(z), 1.2, dec);
    const std::int64_t above = below + 1;
    const double edge = s.grid.cells.cellRaMin[(std::size_t)above];
    check(s.grid.cells.cellRaMax[(std::size_t)below] == edge, "two sectors share the edge");
    const std::vector<double> sky {std::nextafter(edge, 0.), dec, z};
    const std::vector<double> onEdge = otswap::toCartesian({edge + 1e-13, dec, z}, table());
    const Point x {onEdge[0], onEdge[1], onEdge[2]};
    const Point step {-otswap::maths::sin(edge) * 0.5, otswap::maths::cos(edge) * 0.5, 0.01};
    const std::size_t i = s.add_at(x, sky, step);
    check(s.start_cell(i) == below && internal::lightcone_cell_at(s.grid, x.data()) == above &&
          s.end_cell(i) == above, "the sky position lies below the edge, the position and the end above it");
    const double tracerSky[3] = {table().distanceAt(z), sky[0], sky[1]};
    const LightconeFace in = internal::lightcone_face(s.grid, above, x.data(), step.data(), false, tracerSky);
    const LightconeFace none = internal::lightcone_face(s.grid, above, x.data(), step.data(), false);
    check(none.coordinate == -1, "without the sky position no face is found");
    check(in.coordinate == 2 && in.side == -1 && in.t == 0. && in.point == x, "the entry is the lower RA face at x");
    double r = 0., ra = 0., de = 0.;
    internal::lightcone_sky(x.data(), r, ra, de);
    check(std::fabs(ra - edge) < 1e-12, "x lies on that face");
    const Point before {x[0] - 1e-6 * step[0], x[1] - 1e-6 * step[1], x[2] - 1e-6 * step[2]};
    check(internal::lightcone_cell_at(s.grid, before.data()) == below,
          "the cell across the face is the cell just before x along the segment");
    const LightconeFace out = internal::lightcone_face(s.grid, below, x.data(), step.data(), true);
    check(out.coordinate == 2 && out.side == 1 && out.t == 0. && out.normal[0] == -in.normal[0] &&
          out.normal[1] == -in.normal[1] && out.normal[2] == -in.normal[2], "the exit from the sky cell mirrors it");
    const LightconeDivergence div = s.run();
    const DivergenceGrid& g = s.grid.cells;
    check(div.nRejected == 0 && div.used[i] == 1 && g.nDisplacements[(std::size_t)below] == 1 &&
          g.nDisplacements[(std::size_t)above] == 1, "the displacement counts once in each cell");
    check(g.divergence[(std::size_t)below] > 0. && g.divergence[(std::size_t)below] == -g.divergence[(std::size_t)above],
          "it leaves one cell and enters the other with opposite terms");
  }

  // Every valid displacement read being rejected is an error, as is, with a reused grid, no valid displacement left to read.
  group("all displacements rejected, or none left inside a reused grid, is refused");
  {
    Sweep s;
    s.add(1.0, 0.0, 0.32, radial(1.0, 0.0, 2000.));
    s.add(1.1, 0.0, 0.32, radial(1.1, 0.0, 2000.));
    throws_naming([&] { s.run(); },
                  "all 2 valid displacements were rejected (masked cells, outside or below the grid)", "all rejected");
    Sweep one;
    one.add(1.0, 0.0, 0.32, radial(1.0, 0.0, 2000.));
    throws_naming([&] { one.run(); }, "the 1 valid displacement was rejected", "the only one rejected");
    Sweep r;
    r.reused = true;
    r.add(3.0, 0.0, 0.32, {0.1, 0.1, 0.1});
    r.add(3.1, 0.0, 0.32, {0.1, 0.1, 0.1});
    r.add(3.2, 0.0, 0.32, {0.1, 0.1, 0.1});
    r.valid = {1, 1, 0};
    throws_naming([&] { r.run(); },
                  "all 2 kept tracers with a valid displacement lie outside the reused grid; build the grid from "
                  "these data or pass one that covers them", "none inside the reused grid");
    Case c;
    c.displacement = std::vector<double>(12, 1e4);
    throws_naming([&] { findVoidsLightcone(c.tracersSky, c.displacement, {}, c.randomsSky, c.mask, table(), c.profile,
                                           c.config); },
                  "all 4 valid displacements were rejected", "through the public entry");
  }

  // The guards that cannot fire on consistent inputs are raised after the loops as internal errors.
  group("internal guards are raised after the loops");
  {
    Sweep s;
    s.add(3.0, 0.0, 0.32, {0.1, 0.1, 0.1});
    throws_naming([&] { s.run(); }, "internal error: tracer 0 (from 0) lies outside every cell of the grid built",
                  "a tracer outside a new grid");
    Sweep z;
    z.add_centre({0.1, 0.1, 0.1});
    z.sky[2] = 5.;
    throws_naming([&] { z.run(); }, "internal error: the distance of the redshift of tracer 0 (from 0)",
                  "a redshift outside the table");
    Sweep m;
    const std::size_t i = m.add(1.5, 0.0, 0.32, radial(1.5, 0.0, 40.));
    check(m.end_cell(i) != m.start_cell(i), "the displacement crosses a face");
    m.grid.cells.mpsProfile.mps = {-1., -1.};
    throws_naming([&] { m.run(); }, "internal error: the redshift or mps lookup failed at distance",
                  "an mps lookup that fails");
  }

  // The report gives the divergence line at Detailed, and at Normal only when displacements were rejected.
  group("the divergence line is written at Detailed, and at Normal with rejections");
  {
    Case c;
    const std::string detailed = c.run(otswap::Verbosity::Detailed, "not implemented yet");
    check(detailed.find("bitvf: divergence of the valid displacements of 4 tracers, 0 displacements rejected "
                        "(masked cells, outside or below the grid) in ") != std::string::npos, "Detailed: " + detailed);
    const std::string normal = c.run(otswap::Verbosity::Normal, "not implemented yet");
    check(normal.find("bitvf: divergence") == std::string::npos, "no line at Normal without rejections");
    c.displacement[0] = 1e4;
    const std::string rejected = c.run(otswap::Verbosity::Normal, "not implemented yet");
    check(rejected.find("bitvf: divergence of the valid displacements of 3 tracers, 1 displacement rejected "
                        "(masked cells, outside or below the grid) in ") != std::string::npos, "Normal: " + rejected);
  }

  // The divergence of a fixed set of displacements is identical bit for bit with any number of threads and matches the recorded hash.
  group("the divergence is the same with 1 and 4 threads, and matches the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = divergence_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = divergence_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "divergence hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_divergence");
}
