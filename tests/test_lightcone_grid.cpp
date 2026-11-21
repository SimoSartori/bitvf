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
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
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

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;

  constexpr std::uint64_t kRecordedHash = 0x7e64a2eb37c2bc70ull;
  constexpr double kCaseCellSize = 12.;

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  otswap::MpsProfile flat_profile ()
  {
    otswap::MpsProfile p;
    p.redshift = {0.1, 0.5};
    p.mps = {5., 6.};
    return p;
  }

  std::vector<double> objects (const std::size_t n, const double ra0, const double ra1, const double dec0,
                               const double dec1, const double z0, const double z1, const unsigned seed)
  {
    std::mt19937 rng(seed);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < n; ++i) {
      double ra = std::fmod(uniform(ra0, ra1), kTwoPi);
      if (ra < 0.) ra += kTwoPi;
      sky.push_back(ra);
      sky.push_back(otswap::maths::asin(uniform(otswap::maths::sin(dec0), otswap::maths::sin(dec1))));
      sky.push_back(uniform(z0, z1));
    }
    return sky;
  }

  struct Built {
    internal::LightconeExtent extent;
    internal::LightconeGrid grid;
  };

  Built build (const std::vector<double>& tracers, const std::vector<double>& randoms, const double cellSizeMps,
               const otswap::Mask& mask, const otswap::MpsProfile& profile = flat_profile())
  {
    Built b;
    b.extent = internal::lightcone_extent(tracers, {}, {}, randoms, table());
    DivergenceGrid record;
    record.geometry = otswap::Geometry::Lightcone;
    record.distanceMin = b.extent.distanceMin;
    record.distanceMax = b.extent.distanceMax;
    record.raMin = b.extent.raMin;
    record.raMax = b.extent.raMax;
    record.decMin = b.extent.decMin;
    record.decMax = b.extent.decMax;
    record.cellSizeMps = cellSizeMps;
    record.mpsProfile = profile;
    b.grid = internal::lightcone_grid(record, table(), mask);
    return b;
  }

  const otswap::Mask& full_mask ()
  {
    static const otswap::Mask m(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    return m;
  }

  struct Case {
    std::vector<double> tracersSky = objects(40, 0.5, 2., -0.3, 0.3, 0.25, 0.4, 1);
    std::vector<double> displacement = std::vector<double>(120, 0.1);
    std::vector<double> randomsSky = objects(200, 0.5, 2., -0.3, 0.3, 0.25, 0.4, 2);
    otswap::DistanceTable distances = table();
    otswap::MpsProfile profile = flat_profile();
    const otswap::Mask* mask = &full_mask();
    Config config;
    const DivergenceGrid* grid = nullptr;
    Case ()
    {
      config.verbosity = otswap::Verbosity::Silent;
      config.divergenceCellSize = kCaseCellSize;
    }
    void run () const
    {
      findVoidsLightcone(tracersSky, displacement, {}, randomsSky, *mask, distances, profile, config, {}, {}, grid);
    }
    DivergenceGrid built () const
    {
      return lightconeGrid(tracersSky, randomsSky, *mask, distances, profile, *config.divergenceCellSize);
    }
  };

  Case small_footprint ()
  {
    Case c;
    c.tracersSky = objects(40, 1., 1.1, 0., 0.1, 0.3, 0.32, 30);
    c.randomsSky = objects(200, 1., 1.1, 0., 0.1, 0.3, 0.32, 31);
    return c;
  }

  otswap::Result reconstruction (const Case& c)
  {
    otswap::Result rec;
    const std::size_t n = c.tracersSky.size() / 3;
    rec.geometry = otswap::Geometry::Lightcone;
    rec.nObjects = n;
    rec.nRealizations = 1;
    rec.config.nRealizations = 1;
    rec.tracersSky = c.tracersSky;
    rec.tracers = otswap::toCartesian(c.tracersSky, c.distances);
    rec.displacement = c.displacement;
    rec.valid.assign(n, 1);
    rec.outsideMask.assign(n, 0);
    rec.outsideRedshiftCut.assign(n, 0);
    rec.mpsProfile = c.profile;
    rec.mpsProfile.count = {2, 2};
    rec.distances = std::make_shared<const otswap::DistanceTable>(c.distances);
    otswap::recomputeMeans(rec);
    return rec;
  }

  std::string grid_lines (const std::string& log, const bool built)
  {
    static const std::regex time(" in [0-9]+\\.[0-9]+ s$");
    std::istringstream in(log);
    std::string line, kept;
    while (std::getline(in, line)) {
      if (line.find("bitvf: kept ") != 0 && line.find("bitvf: grid") != 0) continue;
      if (line.find(" outside the reused grid") != std::string::npos) continue;
      line = std::regex_replace(line, time, "");
      if (built) {
        for (const auto& swap : {std::pair<std::string, std::string>{" mps (given)", " mps (reused)"},
                                 std::pair<std::string, std::string>{"(from the tracers and randoms)",
                                                                     "(from the reused grid)"}}) {
          const std::size_t at = line.find(swap.first);
          if (at != std::string::npos) line.replace(at, swap.first.size(), swap.second);
        }
      }
      kept += line + "\n";
    }
    return kept;
  }

  /// @p f must get through the grid: it returns, or finds no void.
  template <typename F>
  void passes_grid (F&& f, const std::string& what)
  {
    try {
      f();
    }
    catch (const Error& e) {
      const std::string m = e.what();
      check(m.find("no void found") != std::string::npos, what + " (the message is \"" + m + "\")");
      return;
    }
    catch (const std::exception& e) {
      check(false, what + " (threw the wrong type: " + e.what() + ")");
      return;
    }
    check(true, what);
  }

  std::string run_log (const Case& c)
  {
    Case d = c;
    std::ostringstream log;
    d.config.verbosity = otswap::Verbosity::Detailed;
    d.config.log = &log;
    passes_grid([&] { d.run(); }, "run to the grid");
    return log.str();
  }

  bool same_record (const DivergenceGrid& a, const DivergenceGrid& b)
  {
    return a.geometry == b.geometry && a.nCells == b.nCells && a.distanceMin == b.distanceMin &&
           a.distanceMax == b.distanceMax && a.raMin == b.raMin && a.raMax == b.raMax && a.decMin == b.decMin &&
           a.decMax == b.decMax && a.cellSizeMps == b.cellSizeMps && a.mpsProfile.redshift == b.mpsProfile.redshift &&
           a.mpsProfile.mps == b.mpsProfile.mps;
  }

  void fails (const Case& c, const std::string& text, const std::string& what)
  {
    throws_naming([&] { c.run(); }, text, what);
  }

  void check_objects (const internal::LightconeGrid& g, const std::vector<double>& sky, const std::string& what)
  {
    const DivergenceGrid& c = g.cells;
    bool inside = true, held = true;
    for (std::size_t i = 0; i < sky.size() / 3; ++i) {
      const double d = table().distanceAt(sky[3*i+2]);
      const std::int64_t k = internal::lightcone_cell(g, d, sky[3*i], sky[3*i+1]);
      if (k < 0) {
        inside = false;
        continue;
      }
      const std::size_t j = (std::size_t)k;
      double ra = sky[3*i];
      while (ra < c.cellRaMin[j]) ra += kTwoPi;
      while (ra >= c.cellRaMax[j]) ra -= kTwoPi;
      held = held && d >= c.cellDistanceMin[j] && d < c.cellDistanceMax[j] && sky[3*i+1] >= c.cellDecMin[j] &&
             sky[3*i+1] < c.cellDecMax[j] && ra >= c.cellRaMin[j];
    }
    check(inside, what + ": every object lies in a cell");
    check(held, what + ": the cell's limits hold the object");
  }

  void check_volumes (const internal::LightconeGrid& g, const std::string& what)
  {
    const DivergenceGrid& c = g.cells;
    bool ok = true;
    for (std::size_t i = 0; i + 1 < g.shellEdge.size(); ++i) {
      const double r0 = g.shellEdge[i], r1 = g.shellEdge[i + 1];
      for (std::int64_t r = g.shellFirstRing[i]; r < g.shellFirstRing[i + 1]; ++r) {
        const std::size_t first = (std::size_t)g.ringFirstCell[(std::size_t)r];
        const std::size_t end = (std::size_t)g.ringFirstCell[(std::size_t)r + 1];
        double sum = 0.;
        for (std::size_t k = first; k < end; ++k) sum += c.volume[k];
        const double span = c.cellRaMax[end - 1] - c.cellRaMin[first];
        const double expected = 4. * kPi / 3. * (r1*r1*r1 - r0*r0*r0) *
          (std::sin(g.ringDecMax[(std::size_t)r]) - std::sin(g.ringDecMin[(std::size_t)r])) / 2. * span / kTwoPi;
        ok = ok && std::fabs(sum - expected) <= 1e-9 * expected;
      }
    }
    check(ok, what + ": the cell volumes of each ring add up to its share of the shell");
  }

  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  template <typename T>
  void hash_vector (std::uint64_t& h, const std::vector<T>& v)
  {
    hash_bytes(h, v.data(), v.size() * sizeof(T));
  }

  std::uint64_t grid_hash ()
  {
    const std::vector<double> tracers = objects(3000, 5.5, 8., -0.4, 0.5, 0.2, 0.45, 3);
    const std::vector<double> randoms = objects(6000, 5.5, 8., -0.4, 0.5, 0.2, 0.45, 4);
    std::vector<double> pixels(12, 1.);
    pixels[4] = 0.;
    const otswap::Mask mask(pixels, otswap::PixelOrdering::Ring);
    const Built b = build(tracers, randoms, 14., mask);
    const DivergenceGrid& c = b.grid.cells;
    std::vector<std::int64_t> cells;
    for (std::size_t i = 0; i < tracers.size() / 3; ++i)
      cells.push_back(internal::lightcone_cell(b.grid, table().distanceAt(tracers[3*i+2]), tracers[3*i],
                                               tracers[3*i+1]));
    std::uint64_t h = 0xcbf29ce484222325ull;
    hash_bytes(h, &b.extent, sizeof b.extent);
    hash_vector(h, b.grid.shellEdge);
    hash_vector(h, b.grid.ringDecMin);
    hash_vector(h, b.grid.ringDecMax);
    hash_vector(h, b.grid.ringFirstCell);
    hash_vector(h, c.cellRaMin);
    hash_vector(h, c.cellRaMax);
    hash_vector(h, c.centre);
    hash_vector(h, c.centreSky);
    hash_vector(h, c.volume);
    hash_vector(h, c.masked);
    hash_vector(h, cells);
    return h;
  }

}

int main ()
{
  // The extent covers the objects read, with the right ascensions cut at their largest gap on the circle.
  group("the extent is that of the objects read, cut at the largest gap of right ascension");
  {
    const std::vector<double> sky {6.2, 0.1, 0.3,  0.1, -0.2, 0.35,  0.3, 0.05, 0.4,  3., 0.4, 0.2};
    const std::vector<std::uint8_t> keep {1, 1, 1, 0};
    const internal::LightconeExtent e = internal::lightcone_extent(sky, keep, {}, {}, table());
    check(e.raMin == 6.2 && e.raMax == 0.3 + kTwoPi, "a footprint across RA = 0 starts after the gap");
    check(e.decMin == -0.2 && e.decMax == 0.1 && e.redshiftMin == 0.3 && e.redshiftMax == 0.4,
          "a tracer not read changes nothing");
    check(e.distanceMin == table().distanceAt(0.3) && e.distanceMax == table().distanceAt(0.4),
          "distances are those of the extreme redshifts");
    const std::vector<double> randoms {3.5, 0., 0.3};
    const internal::LightconeExtent u = internal::lightcone_extent(sky, keep, {}, randoms, table());
    check(u.raMin == 3.5 && u.raMax == 0.3 + kTwoPi, "a random in the gap moves the cut");
    const std::vector<double> one {1., 0.2, 0.3};
    const internal::LightconeExtent s = internal::lightcone_extent(one, {}, {}, {}, table());
    check(s.raMin == 1. && s.raMax == 1., "a single object");
  }

  // Shells start a tenth of a cell inside the nearest object and end at least a tenth of a cell beyond the farthest.
  group("the shells: first edge, widths from mps(z) half a step ahead, last edge past the margin");
  {
    const Built b = build(objects(100, 1., 2., -0.2, 0.2, 0.25, 0.4, 5), objects(300, 1., 2., -0.2, 0.2, 0.25, 0.4, 6),
                          9., full_mask());
    const internal::LightconeGrid& g = b.grid;
    const otswap::MpsProfile p = flat_profile();
    const double s = 9.;
    const double zMin = table().redshiftAt(b.extent.distanceMin), zMax = table().redshiftAt(b.extent.distanceMax);
    check(!g.centralSphere && g.shellEdge[0] == b.extent.distanceMin - 0.1 * s * p.at(zMin), "the first edge");
    bool widths = true;
    double width = s * p.at(table().redshiftAt(g.shellEdge[0]));
    for (std::size_t i = 0; i + 1 < g.shellEdge.size(); ++i) {
      const double next = g.shellEdge[i] + s * p.at(table().redshiftAt(g.shellEdge[i] + width / 2.));
      widths = widths && g.shellEdge[i + 1] == next && g.shellCentre[i] == (g.shellEdge[i] + next) / 2. &&
               g.shellRedshift[i] == table().redshiftAt(g.shellCentre[i]);
      width = next - g.shellEdge[i];
    }
    check(widths, "each width is s mps(z) at the previous edge plus half the previous width");
    const double outer = b.extent.distanceMax + 0.1 * s * p.at(zMax);
    const std::size_t n = g.shellEdge.size();
    check(g.shellEdge[n - 1] >= outer && g.shellEdge[n - 2] < outer, "the last edge is the first past the margin");
    check(g.cells.nCells == (std::size_t)g.ringFirstCell.back() && g.shellFirstRing.size() == n,
          "the cell count and the shell layout");
  }

  // Near the observer, a central sphere covers the full sky out to (3 / (4 pi))^(1/3) cells.
  group("a central sphere near the observer, one cell over the full sky");
  {
    const Built b = build(objects(50, 0., 6.2, -1., 1., 0.0005, 0.02, 7), objects(200, 0., 6.2, -1., 1., 0.0005, 0.02, 8),
                          4., full_mask());
    const internal::LightconeGrid& g = b.grid;
    const double zMin = table().redshiftAt(b.extent.distanceMin);
    check(g.centralSphere && g.shellEdge[0] == 0. &&
          g.shellEdge[1] == otswap::maths::cbrt(3. / (4. * kPi)) * 4. * flat_profile().at(zMin),
          "the sphere's radius");
    check(g.shellFirstRing[1] == 1 && g.ringFirstCell[1] == 1 && g.ringFull[0] == 1, "one ring of one cell");
    check(g.cells.cellDecMin[0] == -kPi / 2 && g.cells.cellDecMax[0] == kPi / 2 && g.cells.cellRaMin[0] == 0. &&
          g.cells.cellRaMax[0] == kTwoPi && g.cells.centre[0] == 0. && g.cells.centre[2] == 0.,
          "full-sky limits, centre at the observer");
    check(std::fabs(g.cells.volume[0] - 4. * kPi / 3. * std::pow(g.shellEdge[1], 3)) <= 1e-12 * g.cells.volume[0],
          "the sphere's volume");
    check(g.redshiftFirstEdge == 0., "the grid's range starts at z = 0");
    check_objects(g, objects(50, 0., 6.2, -1., 1., 0.0005, 0.02, 7), "central sphere");
  }

  // Rings cover the declination range padded by a tenth of a cell; near a pole the last ring is a cap on the pole.
  group("the rings: padded declination range, polar caps, one sector on the pole");
  {
    const Built b = build(objects(100, 1., 2., -0.2, 0.3, 0.25, 0.3, 9), objects(300, 1., 2., -0.2, 0.3, 0.25, 0.3, 10),
                          6., full_mask());
    const internal::LightconeGrid& g = b.grid;
    bool padded = true, contiguous = true;
    for (std::size_t i = 0; i + 1 < g.shellEdge.size(); ++i) {
      const double cell = 6. * g.shellMps[i], R = g.shellCentre[i];
      const std::size_t r0 = (std::size_t)g.shellFirstRing[i], r1 = (std::size_t)g.shellFirstRing[i + 1];
      padded = padded && g.ringDecMin[r0] == b.extent.decMin - 0.1 * cell / R &&
               g.ringDecMax[r1 - 1] == b.extent.decMax + 0.1 * cell / R;
      for (std::size_t r = r0; r + 1 < r1; ++r) contiguous = contiguous && g.ringDecMax[r] == g.ringDecMin[r + 1];
    }
    check(padded, "the first and last ring edges");
    check(contiguous, "consecutive rings share their edge");

    std::vector<double> near = objects(400, 0., 6.28, 0.5, 1.5705, 0.25, 0.3, 12);
    near.insert(near.end(), {1., 1.5705, 0.27});
    const Built c = build(objects(100, 0., 6.28, 0.5, 1.5705, 0.25, 0.3, 11), near, 10., full_mask());
    const internal::LightconeGrid& h = c.grid;
    bool caps = true;
    for (std::size_t i = 0; i + 1 < h.shellEdge.size(); ++i) {
      const std::size_t last = (std::size_t)h.shellFirstRing[i + 1] - 1;
      const std::size_t cell = (std::size_t)h.ringFirstCell[last];
      caps = caps && h.ringDecMax[last] == kPi / 2 && h.ringPole[last] == 1 &&
             h.ringFirstCell[last + 1] == h.ringFirstCell[last] + 1 && h.cells.centreSky[3*cell+1] == kPi / 2 &&
             h.cells.cellRaMin[cell] == 0. && h.cells.cellRaMax[cell] == kTwoPi;
    }
    check(caps, "each shell ends with a polar cap of one sector centred on the pole");
    check_objects(h, objects(100, 0., 6.28, 0.5, 1.5705, 0.25, 0.3, 11), "polar cap");
    const std::size_t cap = (std::size_t)h.ringFirstCell[(std::size_t)h.shellFirstRing[1] - 1];
    check(internal::lightcone_cell(h, h.shellCentre[0], 1., kPi / 2) == (std::int64_t)cap,
          "a point on the north pole is in the cap");
  }

  // A footprint narrow in right ascension gets sectors centred on its padded range; a wide one gets full rings from RA = 0.
  group("the sectors: partial rings on the padded RA range, full rings from 0 to 2 pi");
  {
    const Built b = build(objects(100, 6., 6.6, -0.2, 0.2, 0.25, 0.3, 13), objects(300, 6., 6.6, -0.2, 0.2, 0.25, 0.3, 14),
                          4., full_mask());
    const internal::LightconeGrid& g = b.grid;
    bool partial = true;
    for (std::size_t i = 0; i + 1 < g.shellEdge.size(); ++i)
      for (std::int64_t r = g.shellFirstRing[i]; r < g.shellFirstRing[i + 1]; ++r) {
        const std::size_t first = (std::size_t)g.ringFirstCell[(std::size_t)r];
        const std::size_t end = (std::size_t)g.ringFirstCell[(std::size_t)r + 1];
        const double cell = 4. * g.shellMps[i];
        const double RcosDec = g.shellCentre[i] * std::cos(g.cells.centreSky[3*first+1]);
        const double pad = 0.1 * cell / RcosDec;
        const double lo = g.cells.cellRaMin[first], hi = g.cells.cellRaMax[end - 1];
        partial = partial && !g.ringFull[(std::size_t)r] &&
                  std::fabs(lo - (b.extent.raMin - pad)) < 1e-12 && std::fabs(hi - (b.extent.raMax + pad)) < 1e-12 &&
                  (double)(end - first) == std::ceil(RcosDec * (hi - lo) / cell - 1e-9);
      }
    check(b.extent.raMin > 5.9 && b.extent.raMax > kTwoPi, "the footprint crosses RA = 0");
    check(partial, "every ring spans the padded range in ceil(arc / cell) sectors");
    check_objects(g, objects(100, 6., 6.6, -0.2, 0.2, 0.25, 0.3, 13), "partial rings");
    check_volumes(g, "partial rings");

    const Built w = build(objects(300, 0., 6.283, -0.2, 0.2, 0.25, 0.3, 15),
                          objects(900, 0., 6.283, -0.2, 0.2, 0.25, 0.3, 16), 10., full_mask());
    const internal::LightconeGrid& f = w.grid;
    bool full = true;
    for (std::size_t r = 0; r < f.ringFull.size(); ++r) {
      const std::size_t first = (std::size_t)f.ringFirstCell[r], end = (std::size_t)f.ringFirstCell[r + 1];
      full = full && f.ringFull[r] == 1 && f.cells.cellRaMin[first] == 0. && f.cells.cellRaMax[end - 1] == kTwoPi;
    }
    check(full, "every ring is full, from 0 to 2 pi");
    const std::size_t r = (std::size_t)f.shellFirstRing[1];
    const std::size_t first = (std::size_t)f.ringFirstCell[r];
    const double dec = (f.ringDecMin[r] + f.ringDecMax[r]) / 2.;
    check(internal::lightcone_cell(f, f.shellCentre[1], 0., dec) == (std::int64_t)first &&
          internal::lightcone_cell(f, f.shellCentre[1], kTwoPi, dec) == (std::int64_t)first &&
          internal::lightcone_cell(f, f.shellCentre[1], -1e-300, dec) == (std::int64_t)first,
          "RA = 0 and 2 pi fall in the first sector");
    check_volumes(f, "full rings");
  }

  // Over the whole sky, the cells fill the spherical shell between the first and the last edge.
  group("over the whole sky the volumes add up to the spherical shell");
  {
    std::vector<double> sky = objects(1500, 0., 6.283, -1.5706, 1.5706, 0.25, 0.3, 18);
    for (int k = 0; k < 2000; ++k) sky.insert(sky.end(), {kTwoPi * k / 2000., 0., 0.27});
    sky.insert(sky.end(), {1., 1.5706, 0.27,  2., -1.5706, 0.27});
    const Built b = build(objects(500, 0., 6.283, -1.5706, 1.5706, 0.25, 0.3, 17), sky, 14., full_mask());
    const internal::LightconeGrid& g = b.grid;
    double sum = 0.;
    for (const double v : g.cells.volume) sum += v;
    const double r0 = g.shellEdge.front(), r1 = g.shellEdge.back();
    const double expected = 4. * kPi / 3. * (r1*r1*r1 - r0*r0*r0);
    check(std::fabs(sum - expected) <= 1e-9 * expected, "sum of volumes");
    check_volumes(g, "whole sky");
  }

  // A point on an edge belongs to the cell above it; the last distance edge and anything beyond the grid are outside.
  group("the cell of a point: edges belong to the cell above, the grid's last edges are outside");
  {
    const Built b = build(objects(100, 1., 2., -0.2, 0.2, 0.25, 0.4, 19), objects(300, 1., 2., -0.2, 0.2, 0.25, 0.4, 20),
                          8., full_mask());
    const internal::LightconeGrid& g = b.grid;
    const DivergenceGrid& c = g.cells;
    const std::size_t r = (std::size_t)g.shellFirstRing[1] + 1;
    const std::size_t k = (std::size_t)g.ringFirstCell[r] + 1;
    const double ra = (c.cellRaMin[k] + c.cellRaMax[k]) / 2., dec = (c.cellDecMin[k] + c.cellDecMax[k]) / 2.;
    check(internal::lightcone_cell(g, g.shellEdge[1], ra, dec) == (std::int64_t)k, "a shell edge: the shell above");
    check(internal::lightcone_cell(g, g.shellCentre[1], ra, c.cellDecMin[k]) == (std::int64_t)k,
          "a ring edge: the ring above");
    check(internal::lightcone_cell(g, g.shellCentre[1], c.cellRaMin[k], dec) == (std::int64_t)k,
          "a sector edge: the sector above");
    check(internal::lightcone_cell(g, g.shellEdge.front(), ra, dec) >= 0, "the first edge is inside");
    check(internal::lightcone_cell(g, g.shellEdge.back(), ra, dec) == -1, "the last edge is outside");
    check(internal::lightcone_cell(g, g.shellCentre[1], ra, 1.) == -1 &&
          internal::lightcone_cell(g, g.shellCentre[1], 4., dec) == -1 &&
          internal::lightcone_cell(g, g.shellEdge.front() * 0.99, ra, dec) == -1, "outside in each coordinate");
    check(internal::lightcone_cell(g, g.shellCentre[1], ra + kTwoPi, dec) == (std::int64_t)k &&
          internal::lightcone_cell(g, g.shellCentre[1], ra - 3 * kTwoPi, dec) == (std::int64_t)k,
          "any turn of right ascension");
    check_objects(g, objects(100, 1., 2., -0.2, 0.2, 0.25, 0.4, 19), "tracers");
    check_objects(g, objects(300, 1., 2., -0.2, 0.2, 0.25, 0.4, 20), "randoms");
    bool centres = true;
    for (std::size_t j = 0; j < c.nCells; ++j)
      centres = centres && c.centreSky[3*j] >= 0. && c.centreSky[3*j] < kTwoPi &&
                internal::lightcone_cell(g, g.shellCentre[0], c.centreSky[3*j], c.centreSky[3*j+1]) >= -1;
    check(centres, "centres hold a folded right ascension");
  }

  // Declination edges that round to within 1e-12 of 0 are exactly 0.
  group("declination edges near 0 are exactly 0");
  {
    DivergenceGrid record;
    record.geometry = otswap::Geometry::Lightcone;
    record.distanceMin = 800.;
    record.distanceMax = 900.;
    record.raMin = 1.;
    record.raMax = 1.5;
    record.decMin = -0.3;
    record.decMax = 0.3;
    record.cellSizeMps = 4.;
    record.mpsProfile = flat_profile();
    const internal::LightconeGrid g = internal::lightcone_grid(record, table(), full_mask());
    bool snapped = true, zero = false;
    for (std::size_t r = 0; r < g.ringDecMin.size(); ++r) {
      snapped = snapped && (g.ringDecMin[r] == 0. || std::fabs(g.ringDecMin[r]) >= 1e-12);
      zero = zero || g.ringDecMin[r] == 0.;
    }
    check(snapped, "no edge within 1e-12 of 0 but 0 itself");
    check(zero, "a symmetric range with an even ring count has an edge at 0");
  }

  // A cell is masked when its centre falls on an unobserved pixel of the mask.
  group("masked cells: the mask at the cell centre");
  {
    std::vector<double> pixels(12, 1.);
    pixels[4] = 0.;
    pixels[5] = 0.;
    const otswap::Mask mask(pixels, otswap::PixelOrdering::Ring);
    const Built b = build(objects(300, 0., 6.283, -0.5, 0.5, 0.25, 0.3, 21),
                          objects(900, 0., 6.283, -0.5, 0.5, 0.25, 0.3, 22), 12., mask);
    const DivergenceGrid& c = b.grid.cells;
    bool same = true;
    std::size_t masked = 0;
    for (std::size_t j = 0; j < c.nCells; ++j) {
      same = same && c.masked[j] == (mask.allows(c.centreSky[3*j], c.centreSky[3*j+1]) ? 0 : 1);
      masked += c.masked[j];
    }
    check(same, "masked follows the mask at the centre");
    check(masked > 0 && masked < c.nCells, "some cells masked, not all");
  }

  // The grid and the cells of the objects do not depend on the number of threads, and match the recorded hash.
  group("the same grid with 1 and 4 threads, and the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = grid_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = grid_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "grid hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  // A profile that is not positive over the redshifts of the objects read is refused before the grid is built.
  group("a profile not positive over the objects' redshifts is refused");
  {
    Case c;
    c.profile.redshift = {0.2, 0.3};
    c.profile.mps = {6., 1.};
    fails(c, "the mps profile is not positive at redshift 0.3", "profile");
    fails(c, "(tracers and randoms span [", "profile range");
  }

  // A profile positive over the objects but not over the padded grid is refused, naming the grid's range.
  group("a profile not positive inside the grid's range is refused");
  {
    Case c;
    for (std::size_t i = 0; i < c.tracersSky.size() / 3; ++i) c.tracersSky[3*i+2] = 0.29 + 0.1 * (i % 2);
    for (std::size_t i = 0; i < c.randomsSky.size() / 3; ++i) c.randomsSky[3*i+2] = 0.29 + 0.1 * (i % 2);
    c.profile.redshift = {0.3, 0.4};
    c.profile.mps = {1., 6.};
    c.config.divergenceCellSize = 1000.;
    fails(c, "the mps profile is not positive at redshift 0.2", "grid range");
    fails(c, ", inside the grid's range [", "grid range text");
  }

  // A grid that reaches beyond the distance table is refused, advising a wider table.
  group("a grid reaching outside the distance table is refused");
  {
    Case c;
    c.distances = otswap::DistanceTable(0.3, 0.7, -1., 0., 0., 0.4, 1000);
    c.tracersSky[2] = 0.3995;
    fails(c, "outside the distance table [0, ", "outer edge");
    fails(c, "; extend the table", "advice");
    Case s;
    s.distances = otswap::DistanceTable(0.3, 0.7, -1., 0., 0.001, 0.1, 1000);
    s.tracersSky = objects(40, 0.5, 2., -0.3, 0.3, 0.002, 0.01, 23);
    s.randomsSky = objects(200, 0.5, 2., -0.3, 0.3, 0.002, 0.01, 24);
    fails(s, "the grid reaches distance 0, outside the distance table [", "central sphere");
  }

  // A footprint too thin to hold one declination ring is refused, naming the shell.
  group("a shell with no declination ring is refused");
  {
    Case c;
    for (std::size_t i = 0; i < c.tracersSky.size() / 3; ++i) c.tracersSky[3*i+1] = 0.2;
    for (std::size_t i = 0; i < c.randomsSky.size() / 3; ++i) c.randomsSky[3*i+1] = 0.2;
    fails(c, "shell 0 (from 0) holds no declination ring; use a smaller cell or a central sphere", "no ring");
  }

  // A cell size so small that a level of the cell index would exceed 2097151 is refused while building.
  group("too many shells, rings or sectors are refused");
  {
    Case c;
    c.config.divergenceCellSize = 1e-6;
    fails(c, "the grid would need more than 2097151 shells; increase the divergence cell size", "shells");
    Case r;
    for (std::size_t i = 0; i < r.tracersSky.size() / 3; ++i) r.tracersSky[3*i+2] = 0.3;
    for (std::size_t i = 0; i < r.randomsSky.size() / 3; ++i) r.randomsSky[3*i+2] = 0.3;
    r.config.divergenceCellSize = 1e-5;
    fails(r, "the grid would need more than 2097151 rings in shell 0 (from 0); increase the divergence cell size",
          "rings");
    Case s = r;
    for (std::size_t i = 0; i < s.tracersSky.size() / 3; ++i) s.tracersSky[3*i+1] = 0.2 + 5e-7 * (i % 2);
    for (std::size_t i = 0; i < s.randomsSky.size() / 3; ++i) s.randomsSky[3*i+1] = 0.2 + 5e-7 * (i % 2);
    s.config.divergenceCellSize = 1e-4;
    fails(s, "the grid would need more than 2097151 sectors in ring 0 (from 0) of shell 0 (from 0)", "sectors");
  }

  // At normal verbosity the grid line gives the cells, shells (a central sphere counted among them), cell size and extent; detailed adds masked cells and the redshift range.
  group("the grid's report lines");
  {
    Case c;
    std::ostringstream log;
    c.config.verbosity = otswap::Verbosity::Detailed;
    c.config.log = &log;
    passes_grid([&] { c.run(); }, "report case");
    check(log.str().find(" shells, cell size 12 mps (given), distances [") != std::string::npos &&
          log.str().find("] deg (from the tracers and randoms) in ") != std::string::npos, "grid line: " + log.str());
    check(log.str().find("bitvf: grid: 0 of ") != std::string::npos &&
          log.str().find(" cells masked, redshifts [") != std::string::npos, "detailed grid line: " + log.str());
    Case d = small_footprint();
    d.config.divergenceCellSize.reset();
    std::ostringstream normal;
    d.config.verbosity = otswap::Verbosity::Normal;
    d.config.log = &normal;
    passes_grid([&] { d.run(); }, "default cell size");
    check(normal.str().find("cell size 1 mps (default)") != std::string::npos, "default cell size: " + normal.str());
    check(normal.str().find("central sphere") == std::string::npos, "no central sphere: " + normal.str());
    Case s;
    s.tracersSky = objects(40, 0.5, 2., -0.3, 0.3, 0.002, 0.01, 25);
    s.randomsSky = objects(200, 0.5, 2., -0.3, 0.3, 0.002, 0.01, 26);
    std::ostringstream sphere;
    s.config.verbosity = otswap::Verbosity::Normal;
    s.config.log = &sphere;
    passes_grid([&] { s.run(); }, "central sphere case");
    check(sphere.str().find(" shells, central sphere included, cell size ") != std::string::npos,
          "central sphere: " + sphere.str());
  }

  // The grid builder gives the grid a search builds from the same data, with only the fields that define its cells.
  group("the grid builder gives the search's grid, defining fields only");
  {
    const Case c;
    const DivergenceGrid g = c.built();
    const Built b = build(c.tracersSky, c.randomsSky, kCaseCellSize, full_mask());
    check(same_record(g, b.grid.cells), "the array overload builds the grid of the data");
    check(g.cellRaMin.empty() && g.centre.empty() && g.volume.empty() && g.masked.empty() && std::isnan(g.cellSize) &&
          std::isnan(g.mps) && g.nCellsAxis[0] == 0, "no per-cell or box field is set");
    check(run_log(c).find("grid of " + std::to_string(g.nCells) + " cells in ") != std::string::npos,
          "the search's grid has the same cells");
    const otswap::Result rec = reconstruction(c);
    const DivergenceGrid r = lightconeGrid(rec, c.randomsSky, full_mask(), kCaseCellSize);
    check(same_record(r, g) && r.mpsProfile.count == rec.mpsProfile.count, "the result overload builds the same grid");
    Case d = small_footprint();
    d.config.divergenceCellSize = 1.;
    check(same_record(lightconeGrid(reconstruction(d), d.randomsSky, full_mask()), d.built()),
          "without a cell size, the default for the reconstruction's realizations");
    otswap::Result bad = rec;
    bad.displacement[4] = std::numeric_limits<double>::quiet_NaN();
    throws_naming([&] { lightconeGrid(bad, c.randomsSky, full_mask()); }, "the displacements hold a non-finite value",
                  "the displacements are checked although the grid does not read them");
    bad = rec;
    bad.geometry = otswap::Geometry::Box;
    throws_naming([&] { lightconeGrid(bad, c.randomsSky, full_mask()); },
                  "lightconeGrid: the reconstruction is not a lightcone result; use boxGrid", "geometry");
    bad = rec;
    bad.mpsProfile = {};
    throws_naming([&] { lightconeGrid(bad, c.randomsSky, full_mask()); },
                  "lightconeGrid: the reconstruction has no mpsProfile, which a lightcone search needs", "no profile");
    throws_naming([&] { lightconeGrid({}, c.randomsSky, full_mask(), table(), flat_profile(), 4.); },
                  "lightconeGrid: tracersSky is empty", "array overload prefix");
    throws_naming([&] { lightconeGrid(c.tracersSky, c.randomsSky, full_mask(), table(), flat_profile(), 0.); },
                  "the divergence cell size is 0", "cell size");
  }

  // A grid passed for reuse gives the search the same cells and the same selection as the grid it would build, with any number of threads.
  group("a reused grid gives the same grid and entry as one built in the call, with 1 and 4 threads");
  {
    std::vector<double> pixels(12, 1.);
    pixels[4] = 0.;
    const otswap::Mask mask(pixels, otswap::PixelOrdering::Ring);
    Case c;
    c.mask = &mask;
    c.tracersSky = objects(400, 0.2, 2., -0.3, 0.3, 0.25, 0.4, 27);
    c.displacement.assign(1200, 0.1);
    c.randomsSky = objects(2000, 0.2, 2., -0.3, 0.3, 0.25, 0.4, 28);
    c.profile.redshift = {0.1, 0.3, 0.5};
    c.profile.mps = {5., 5.5, 6.};
    std::string built[2], reused[2], fromResult[2];
    DivergenceGrid grids[2];
    for (int t = 0; t < 2; ++t) {
#ifdef _OPENMP
      omp_set_num_threads(t == 0 ? 1 : 4);
#endif
      built[t] = grid_lines(run_log(c), true);
      grids[t] = c.built();
      Case r = c;
      r.config.divergenceCellSize.reset();
      r.grid = &grids[t];
      reused[t] = grid_lines(run_log(r), false);
      const otswap::Result rec = reconstruction(c);
      const DivergenceGrid g = lightconeGrid(rec, c.randomsSky, mask, kCaseCellSize);
      otswap::Result noProfile = rec;
      noProfile.mpsProfile = {};
      std::ostringstream log;
      Config config;
      config.verbosity = otswap::Verbosity::Detailed;
      config.log = &log;
      passes_grid([&] { findVoidsLightcone(noProfile, c.randomsSky, mask, config, &g); },
                  "result overload with a grid");
      fromResult[t] = grid_lines(log.str(), false);
    }
    check(built[0].find("bitvf: kept ") != std::string::npos && built[0].find("bitvf: grid of ") != std::string::npos &&
          built[0].find("bitvf: grid: ") != std::string::npos && built[0].find(" 0 on unobserved") == std::string::npos,
          "the lines compared hold the selection and the grid: " + built[0]);
    check(reused[0] == built[0], "reused and built agree: " + reused[0] + " against " + built[0]);
    check(fromResult[0] == built[0], "the result overload reuses it without a profile: " + fromResult[0]);
    check(built[1] == built[0] && reused[1] == reused[0] && fromResult[1] == fromResult[0] &&
          same_record(grids[0], grids[1]), "1 and 4 threads agree");
  }

  // With a grid, a reconstruction's own profile is ignored, even one that differs from the grid's or would be refused.
  group("with a grid, the reconstruction's profile is ignored, whatever it holds");
  {
    const Case c;
    const otswap::Result rec = reconstruction(c);
    const DivergenceGrid g = lightconeGrid(rec, c.randomsSky, full_mask(), kCaseCellSize);
    auto lines = [&] (const otswap::Result& r) {
      std::ostringstream log;
      Config config;
      config.verbosity = otswap::Verbosity::Detailed;
      config.log = &log;
      passes_grid([&] { findVoidsLightcone(r, c.randomsSky, full_mask(), config, &g); },
                  "result overload with a grid and its own profile");
      return grid_lines(log.str(), false);
    };
    const std::string same = lines(rec);
    otswap::Result other = rec;
    other.mpsProfile.redshift = {0.2, 0.6};
    other.mpsProfile.mps = {7., 9.};
    check(lines(other) == same, "a profile that differs from the grid's changes nothing");
    otswap::Result invalid = rec;
    invalid.mpsProfile.mps = {0., -1.};
    check(lines(invalid) == same, "a profile that would be refused is not read");
    check(same.find("bitvf: grid of ") != std::string::npos, "the lines compared hold the grid: " + same);
  }

  // Every inconsistent field of a grid to reuse is refused before any row of data is read, naming the field.
  group("a grid to reuse with an inconsistent record is refused");
  {
    const Case base;
    const DivergenceGrid good = base.built();
    const std::string bad = "findVoidsLightcone: the grid to reuse ";
    auto refused = [&] (DivergenceGrid g, const std::string& text, const std::string& what) {
      Case c;
      c.config.divergenceCellSize.reset();
      c.tracersSky[1] = std::numeric_limits<double>::quiet_NaN();
      c.grid = &g;
      fails(c, text, what);
    };
    DivergenceGrid g = good;
    g.geometry = otswap::Geometry::Box;
    refused(g, bad + "is not a lightcone grid", "geometry");
    g = good;
    g.distanceMin = -1.;
    refused(g, bad + "has a distance range [-1, ", "negative distance");
    g = good;
    g.distanceMax = std::numeric_limits<double>::infinity();
    refused(g, "] in coordinate units; both bounds must be finite, with 0 <= min <= max", "infinite distance");
    g = good;
    g.raMax = g.raMin + 7.;
    refused(g, bad + "has an RA range [", "RA wider than 2 pi");
    g = good;
    g.raMax = g.raMin - 0.1;
    refused(g, "] radians; both bounds must be finite, with min <= max and a width of at most 2 pi", "RA reversed");
    g = good;
    g.decMax = 2.;
    refused(g, bad + "has a declination range [", "declination beyond the pole");
    g = good;
    g.decMin = g.decMax + 0.1;
    refused(g, "] radians; it must lie in [-pi/2, pi/2], with min <= max", "declination reversed");
    g = good;
    g.cellSizeMps = 0.;
    refused(g, bad + "has a cell size of 0 mps (cellSizeMps); it must be finite and positive", "cell size");
    g = good;
    g.mpsProfile = {};
    refused(g, "findVoidsLightcone: the mps profile of the grid to reuse is empty", "empty profile");
    g = good;
    g.mpsProfile.mps.pop_back();
    refused(g, "findVoidsLightcone: the redshift and mps lists of the profile of the grid to reuse hold 2 and 1 entries",
            "profile sizes");
    g = good;
    g.mpsProfile.mps[1] = -1.;
    refused(g, "findVoidsLightcone: the mps profile of the grid to reuse holds a mean particle separation of -1 at "
            "entry 1 (from 0); it must be positive", "profile value");
    const char* const fields[] = {"cellSize", "lower", "upper", "mps", "nCellsAxis"};
    for (const char* field : fields) {
      g = good;
      if (std::string(field) == "cellSize") g.cellSize = 1.;
      if (std::string(field) == "lower") g.lower[2] = 0.;
      if (std::string(field) == "upper") g.upper[0] = 0.;
      if (std::string(field) == "mps") g.mps = 1.;
      if (std::string(field) == "nCellsAxis") g.nCellsAxis[1] = 1;
      refused(g, bad + "is a lightcone grid but holds the box field " + field, field);
    }
  }

  // In the array overloads the profile given must be the grid's, bit for bit, in redshift and mps.
  group("a profile that differs from the reused grid's is refused");
  {
    Case c;
    const DivergenceGrid g = c.built();
    c.config.divergenceCellSize.reset();
    c.grid = &g;
    c.profile.mps[1] = 6.5;
    fails(c, "findVoidsLightcone: the mps profile given differs from the grid's at entry 1 (from 0): mps 6.5 against 6; "
          "pass the grid's profile", "mps");
    c.profile = flat_profile();
    c.profile.redshift[0] = std::nextafter(0.1, 1.);
    fails(c, "at entry 0 (from 0): redshift 0.10000000000000002 against 0.10000000000000001", "redshift");
    c.profile.redshift = {0.1, 0.5, 0.6};
    c.profile.mps = {5., 6., 6.};
    fails(c, "findVoidsLightcone: the mps profile given holds 3 entries and the grid's 2; pass the grid's profile",
          "sizes");
    c.profile = flat_profile();
    c.profile.count = {1, 2};
    c.profile.representative = 5.5;
    passes_grid([&] { c.run(); }, "the other fields are not compared");
    c.config.divergenceCellSize = kCaseCellSize;
    fails(c, "the divergence cell size is set together with a grid to reuse", "cell size and grid");
  }

  // A grid to reuse is rebuilt with the call's table and must give its recorded cell count; a rebuild that fails names the grid to reuse.
  group("the rebuild of a grid to reuse: its cell count, and its failures naming it");
  {
    const Case base;
    const DivergenceGrid good = base.built();
    auto refused = [&] (const DivergenceGrid& g, const std::string& text, const std::string& what,
                        const otswap::DistanceTable* distances = nullptr) {
      Case c;
      c.config.divergenceCellSize.reset();
      c.profile = g.mpsProfile;
      if (distances) c.distances = *distances;
      c.grid = &g;
      fails(c, text, what);
    };
    DivergenceGrid g = good;
    g.nCells += 1;
    refused(g, "findVoidsLightcone: the grid to reuse holds " + std::to_string(good.nCells + 1) +
            " cells but the rebuild gives " + std::to_string(good.nCells) +
            ": the record was modified or the distance table differs", "cell count");
    const otswap::DistanceTable shortTable(0.3, 0.7, -1., 0., 0., 0.4, 1000);
    refused(good, "the grid to reuse reaches distance ", "distance table", &shortTable);
    g = good;
    g.mpsProfile.redshift = {0.3, 0.4};
    g.mpsProfile.mps = {1., 6.};
    refused(g, "the mps profile of the grid to reuse is not positive at redshift ", "profile range");
    g = good;
    g.decMin = 0.2;
    g.decMax = 0.2;
    refused(g, "shell 0 (from 0) of the grid to reuse holds no declination ring", "no ring");
    g = good;
    g.cellSizeMps = 1e-6;
    refused(g, "the grid to reuse would need more than 2097151 shells", "cell limit");
    Case far;
    far.config.divergenceCellSize.reset();
    far.randomsSky = objects(200, 0.5, 2., -0.3, 0.3, 0.6, 0.7, 29);
    far.grid = &good;
    fails(far, "no random is left: of 200 randoms, 0 lie on unobserved pixels and 0 outside the redshift cut, 200 "
          "outside the reused grid", "randoms outside the reused grid");
  }

  return report("test_lightcone_grid");
}
