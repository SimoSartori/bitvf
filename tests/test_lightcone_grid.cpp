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

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;

  constexpr std::uint64_t kRecordedHash = 0x49f33c5c5cd90e3aull;

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
    Config config;
    Case ()
    {
      config.verbosity = otswap::Verbosity::Silent;
      config.divergenceCellSize = 4.;
    }
    void run () const
    {
      findVoidsLightcone(tracersSky, displacement, {}, randomsSky, full_mask(), distances, profile, config, {}, {});
    }
  };

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
    const Built b = build(tracers, randoms, 2., mask);
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
                          3., full_mask());
    const internal::LightconeGrid& g = b.grid;
    const otswap::MpsProfile p = flat_profile();
    const double s = 3.;
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
                          4., full_mask());
    const internal::LightconeGrid& g = b.grid;
    bool padded = true, contiguous = true;
    for (std::size_t i = 0; i + 1 < g.shellEdge.size(); ++i) {
      const double cell = 4. * g.shellMps[i], R = g.shellCentre[i];
      const std::size_t r0 = (std::size_t)g.shellFirstRing[i], r1 = (std::size_t)g.shellFirstRing[i + 1];
      padded = padded && g.ringDecMin[r0] == b.extent.decMin - 0.1 * cell / R &&
               g.ringDecMax[r1 - 1] == b.extent.decMax + 0.1 * cell / R;
      for (std::size_t r = r0; r + 1 < r1; ++r) contiguous = contiguous && g.ringDecMax[r] == g.ringDecMin[r + 1];
    }
    check(padded, "the first and last ring edges");
    check(contiguous, "consecutive rings share their edge");

    std::vector<double> near = objects(400, 0., 6.28, 0.5, 1.5705, 0.25, 0.3, 12);
    near.insert(near.end(), {1., 1.5705, 0.27});
    const Built c = build(objects(100, 0., 6.28, 0.5, 1.5705, 0.25, 0.3, 11), near, 4., full_mask());
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
                          objects(900, 0., 6.283, -0.2, 0.2, 0.25, 0.3, 16), 4., full_mask());
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
    const Built b = build(objects(500, 0., 6.283, -1.5706, 1.5706, 0.25, 0.3, 17), sky, 4., full_mask());
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
                          4., full_mask());
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
                          objects(900, 0., 6.283, -0.5, 0.5, 0.25, 0.3, 22), 4., mask);
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
    fails(c, "not implemented yet", "report case");
    check(log.str().find(" shells, cell size 4 mps (given), distances [") != std::string::npos &&
          log.str().find("] deg (from the tracers and randoms) in ") != std::string::npos, "grid line: " + log.str());
    check(log.str().find("bitvf: grid: 0 of ") != std::string::npos &&
          log.str().find(" cells masked, redshifts [") != std::string::npos, "detailed grid line: " + log.str());
    Case d;
    d.config.divergenceCellSize.reset();
    std::ostringstream normal;
    d.config.verbosity = otswap::Verbosity::Normal;
    d.config.log = &normal;
    fails(d, "not implemented yet", "default cell size");
    check(normal.str().find("cell size 1 mps (default)") != std::string::npos, "default cell size: " + normal.str());
    check(normal.str().find("central sphere") == std::string::npos, "no central sphere: " + normal.str());
    Case s;
    s.tracersSky = objects(40, 0.5, 2., -0.3, 0.3, 0.002, 0.01, 25);
    s.randomsSky = objects(200, 0.5, 2., -0.3, 0.3, 0.002, 0.01, 26);
    std::ostringstream sphere;
    s.config.verbosity = otswap::Verbosity::Normal;
    s.config.log = &sphere;
    fails(s, "not implemented yet", "central sphere case");
    check(sphere.str().find(" shells, central sphere included, cell size ") != std::string::npos,
          "central sphere: " + sphere.str());
  }

  return report("test_lightcone_grid");
}
