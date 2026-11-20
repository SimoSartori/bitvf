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
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "bitvf/VoidFinder.h"
#include "otswap/Maths.h"

#include "check.h"
#include "internal.h"

using namespace bitvf;
using internal::LightconeAdjacency;
using internal::LightconeBasins;
using internal::LightconeGrid;

namespace {

  constexpr double kPi = 3.14159265358979323846;
  constexpr double kTwoPi = 2. * kPi;
  constexpr double kHalfPi = kPi / 2.;

  constexpr std::uint64_t kRecordedHash = 0x3847df7da8b6be15ull;

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  /// Cells given only by their limits in distance, declination and right ascension.
  struct Boxes {

    DivergenceGrid grid;

    std::size_t add (const double rLo, const double rHi, const double dLo, const double dHi, const double aLo,
                     const double aHi)
    {
      grid.cellDistanceMin.push_back(rLo);
      grid.cellDistanceMax.push_back(rHi);
      grid.cellDecMin.push_back(dLo);
      grid.cellDecMax.push_back(dHi);
      grid.cellRaMin.push_back(aLo);
      grid.cellRaMax.push_back(aHi);
      return grid.nCells++;
    }

    double area (const std::size_t c, const std::size_t n) const
    {
      return internal::lightcone_shared_area(grid, c, n);
    }
  };

  bool near (const double a, const double b)
  {
    return std::abs(a - b) <= 1e-13 * std::abs(b);
  }

  LightconeGrid built_grid (const bool fullSky)
  {
    std::mt19937 rng(1);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < 4000; ++i) {
      if (fullSky) {
        sky.push_back(uniform(0., kTwoPi));
        sky.push_back(otswap::maths::asin(uniform(-1., 1.)));
        sky.push_back(uniform(0.001, 0.05));
      }
      else {
        sky.push_back(uniform(0.5, 2.));
        sky.push_back(otswap::maths::asin(uniform(otswap::maths::sin(-0.3), otswap::maths::sin(0.3))));
        sky.push_back(uniform(0.25, 0.4));
      }
    }
    const internal::LightconeExtent e = internal::lightcone_extent(sky, {}, {}, {}, table());
    DivergenceGrid record;
    record.geometry = otswap::Geometry::Lightcone;
    record.distanceMin = e.distanceMin;
    record.distanceMax = e.distanceMax;
    record.raMin = e.raMin;
    record.raMax = e.raMax;
    record.decMin = e.decMin;
    record.decMax = e.decMax;
    record.cellSizeMps = fullSky ? 8. : 9.;
    record.mpsProfile.redshift = {0., 0.5};
    record.mpsProfile.mps = {5., 6.};
    const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    return internal::lightcone_grid(record, table(), mask);
  }

  /// Whether a cell is on the border of the grid, from its shell, ring and sector.
  bool on_border (const LightconeGrid& g, const std::int64_t c)
  {
    const std::int64_t ring = std::upper_bound(g.ringFirstCell.begin(), g.ringFirstCell.end(), c) - g.ringFirstCell.begin() - 1;
    const std::int64_t shell = std::upper_bound(g.shellFirstRing.begin(), g.shellFirstRing.end(), ring) -
                               g.shellFirstRing.begin() - 1;
    const std::int64_t nShells = (std::int64_t)g.shellCentre.size();
    if (shell == 0 && !g.centralSphere) return true;
    if (shell == nShells - 1) return true;
    if (ring == g.shellFirstRing[(std::size_t)shell] && g.ringDecMin[(std::size_t)ring] != -kHalfPi) return true;
    if (ring == g.shellFirstRing[(std::size_t)shell + 1] - 1 && g.ringDecMax[(std::size_t)ring] != kHalfPi) return true;
    if (g.ringFull[(std::size_t)ring]) return false;
    return c == g.ringFirstCell[(std::size_t)ring] || c == g.ringFirstCell[(std::size_t)ring + 1] - 1;
  }

  /// Every unmasked cell its own void, numbered in ascending cell index, with the seeds of those voids.
  LightconeBasins one_void_per_cell (DivergenceGrid& cells)
  {
    LightconeBasins b;
    cells.divergence.assign(cells.nCells, -1.);
    cells.divergenceSmoothed.clear();
    cells.voidIndex.assign(cells.nCells, -1);
    cells.overdenseIndex.assign(cells.nCells, -1);
    for (std::size_t c = 0; c < cells.nCells; ++c)
      if (!cells.masked[c]) {
        cells.voidIndex[c] = (std::int64_t)b.voids.seed.size();
        b.voids.seed.push_back((std::int64_t)c);
      }
    return b;
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

  std::uint64_t basins_hash ()
  {
    LightconeGrid g = built_grid(false);
    DivergenceGrid& cells = g.cells;
    for (std::size_t c = 0; c < cells.nCells; c += 23) cells.masked[c] = 1;
    cells.divergence.assign(cells.nCells, 0.);
    cells.nTracers.assign(cells.nCells, 0);
    for (std::size_t c = 0; c < cells.nCells; ++c) {
      const double* x = &cells.centre[3*c];
      cells.divergence[c] = otswap::maths::sin(0.021 * x[0]) * otswap::maths::cos(0.017 * x[1]) +
                            otswap::maths::sin(0.013 * x[2] + 0.5) +
                            0.8 * otswap::maths::sin(0.11 * x[0]) * otswap::maths::sin(0.13 * x[1]) + 0.3;
      cells.nTracers[c] = (std::uint64_t)(c % 7);
    }
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    LightconeBasins b = internal::lightcone_watershed(a, cells);
    BasinCatalog voids, overdense;
    internal::lightcone_basins(g, a, b, voids, overdense);
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const BasinCatalog* cat : {&voids, &overdense}) {
      hash_vector(h, cat->radius);
      hash_vector(h, cat->radiusExt);
      hash_vector(h, cat->edgeFlag);
      hash_vector(h, cat->contamination);
      hash_vector(h, cat->nTracers);
    }
    hash_vector(h, b.voidCoreVolume);
    hash_vector(h, b.voidVolume);
    hash_vector(h, b.overdenseCoreVolume);
    hash_vector(h, b.overdenseVolume);
    return h;
  }

}

int main ()
{
  // Two cells share a face when they touch over a sphere, a cone or a half-plane of constant right ascension, with the area of that face over the overlapping ranges; an edge, a corner or a gap shares none.
  group("the shared area of two cells: spheres, cones, right-ascension planes, edges and gaps");
  {
    Boxes b;
    const std::size_t c = b.add(100., 110., 0., 0.2, 1., 1.2);
    const std::size_t out = b.add(110., 120., 0., 0.2, 1., 1.2);
    const std::size_t up = b.add(100., 110., 0.2, 0.4, 1.1, 1.3);
    const std::size_t side = b.add(100., 110., 0., 0.2, 1.2, 1.4);
    const std::size_t diagonal = b.add(110., 120., 0.2, 0.4, 1., 1.2);
    const std::size_t apart = b.add(100., 110., 0., 0.2, 1.2 + 1e-6, 1.4);
    const double sphere = 110. * 110. * (otswap::maths::sin(0.2) - otswap::maths::sin(0.)) * 0.2;
    check(near(b.area(c, out), sphere) && near(b.area(out, c), sphere), "a sphere, the same from both sides");
    check(near(b.area(c, up), 0.5 * (110. * 110. - 100. * 100.) * otswap::maths::cos(0.2) * (1.2 - 1.1)),
          "a cone over the overlapping right ascensions");
    check(near(b.area(c, side), 0.5 * (110. * 110. - 100. * 100.) * 0.2), "a right-ascension half-plane");
    check(b.area(c, diagonal) == 0., "cells touching along an edge share no face");
    check(b.area(c, apart) == 0., "cells with a gap share no face");
    check(b.area(c, c) == 0., "a cell shares no face with itself");
  }

  // On a full circle the right ascensions are compared across RA = 0, so the last sector shares a face with the first.
  group("the shared area across right ascension 0");
  {
    Boxes b;
    const std::size_t first = b.add(100., 110., 0., 0.2, 0., 0.5);
    const std::size_t last = b.add(100., 110., 0., 0.2, kTwoPi - 0.5, kTwoPi);
    const std::size_t cap = b.add(100., 110., 0.2, kHalfPi, 0., kTwoPi);
    const double plane = 0.5 * (110. * 110. - 100. * 100.) * 0.2;
    check(near(b.area(first, last), plane) && near(b.area(last, first), plane), "the two sectors share a plane");
    check(near(b.area(last, cap), 0.5 * (110. * 110. - 100. * 100.) * otswap::maths::cos(0.2) * 0.5),
          "a sector and the cap above it share a cone");
  }

  // The tolerance on each axis is relative to the cell's own extent: a gap far below the extent is a contact, a gap above it is not, whatever the units.
  group("the contact tolerance scales with the size of the cell");
  {
    Boxes b;
    const std::size_t large = b.add(1e6, 2e6, 0., 0.2, 1., 1.2);
    const std::size_t nearGap = b.add(2e6 + 1e-4, 3e6, 0., 0.2, 1., 1.2);
    const std::size_t small = b.add(100., 110., 0., 0.2, 1., 1. + 1e-3);
    const std::size_t smallGap = b.add(100., 110., 0., 0.2, 1. + 1e-3 + 1e-10, 1.1);
    check(b.area(large, nearGap) > 0., "a gap of 1e-10 of the extent is a contact");
    check(b.area(small, smallGap) == 0., "a gap of 1e-7 of the extent is not");
  }

  // A cell is on the border on shell 0 (unless the central sphere), on the last shell, on the first and last ring of a shell unless its edge is on a pole, and on the first and last sector of a ring that is not a full circle; with every cell its own void, the core touches the border exactly there.
  group("the border of the grid: partial rings, caps and the central sphere");
  {
    for (const bool fullSky : {false, true}) {
      LightconeGrid g = built_grid(fullSky);
      const LightconeAdjacency a = internal::lightcone_adjacency(g);
      LightconeBasins b = one_void_per_cell(g.cells);
      BasinCatalog voids, overdense;
      internal::lightcone_basins(g, a, b, voids, overdense);
      bool match = true, someInside = false, someBorder = false;
      for (std::size_t c = 0; c < g.cells.nCells; ++c) {
        const bool border = on_border(g, (std::int64_t)c);
        match = match && voids.edgeFlag[c] == (border ? 2 : 0);
        someInside = someInside || !border;
        someBorder = someBorder || border;
      }
      const std::string name = fullSky ? "full sky" : "partial sky";
      if (fullSky) check(g.centralSphere && g.ringDecMax.back() == kHalfPi, name + ": central sphere and caps");
      else check(!g.centralSphere && !g.ringFull[0], name + ": partial rings");
      check(someInside && someBorder && g.cells.nCells <= 3000, name + ": inner and border cells, small grid");
      check(match, name + ": edge flag 2 exactly on the border");
    }
  }

  // The faces with area of a cell inside the grid are all shared with its contact list, so a basin with no border and no masked contact has contamination 0, and a basin filling the grid has 1.
  group("contamination: 0 for an inner cell, 1 for a basin filling the grid");
  {
    for (const bool fullSky : {false, true}) {
      LightconeGrid g = built_grid(fullSky);
      const LightconeAdjacency a = internal::lightcone_adjacency(g);
      LightconeBasins b = one_void_per_cell(g.cells);
      BasinCatalog voids, overdense;
      internal::lightcone_basins(g, a, b, voids, overdense);
      bool inner = true, border = true;
      for (std::size_t c = 0; c < g.cells.nCells; ++c) {
        if (on_border(g, (std::int64_t)c)) border = border && voids.contamination[c] > 0. && voids.contamination[c] <= 1.;
        else inner = inner && std::abs(voids.contamination[c]) < 1e-12;
      }
      const std::string name = fullSky ? "full sky" : "partial sky";
      check(inner, name + ": inner cells have contamination 0");
      check(border, name + ": border cells have contamination in (0, 1]");

      LightconeBasins whole;
      whole.voids.seed = {0};
      g.cells.voidIndex.assign(g.cells.nCells, 0);
      internal::lightcone_basins(g, a, whole, voids, overdense);
      check(voids.nBasins == 1 && voids.contamination[0] == 1. && voids.edgeFlag[0] == 2,
            name + ": one void over the grid, contamination 1, edge flag 2");
    }
  }

  // A masked cell marks the basins whose cells share a face with it, and counts as boundary in their contamination; one touching only along an edge or a corner does not.
  group("a masked cell marks its face neighbours only");
  {
    LightconeGrid g = built_grid(false);
    const LightconeAdjacency all = internal::lightcone_adjacency(g);
    std::int64_t m = -1;
    for (std::size_t c = 0; c < g.cells.nCells && m < 0; ++c) {
      bool deep = !on_border(g, (std::int64_t)c);
      for (std::int64_t k = all.contactFirst[c]; k < all.contactFirst[c + 1]; ++k)
        deep = deep && !on_border(g, all.contactNeighbour[(std::size_t)k]);
      if (deep) m = (std::int64_t)c;
    }
    check(m >= 0, "a cell whose contacts are all inside the grid");
    g.cells.masked[(std::size_t)m] = 1;
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    LightconeBasins b = one_void_per_cell(g.cells);
    BasinCatalog voids, overdense;
    internal::lightcone_basins(g, a, b, voids, overdense);
    bool faces = true, edges = true;
    std::size_t nFace = 0, nEdge = 0;
    for (std::int64_t k = all.contactFirst[(std::size_t)m]; k < all.contactFirst[(std::size_t)m + 1]; ++k) {
      const std::int64_t n = all.contactNeighbour[(std::size_t)k];
      const std::int64_t v = g.cells.voidIndex[(std::size_t)n];
      if (internal::lightcone_shared_area(g.cells, (std::size_t)n, (std::size_t)m) > 0.) {
        ++nFace;
        faces = faces && voids.edgeFlag[(std::size_t)v] == 2 && voids.contamination[(std::size_t)v] > 0.;
      }
      else {
        ++nEdge;
        edges = edges && voids.edgeFlag[(std::size_t)v] == 0 && std::abs(voids.contamination[(std::size_t)v]) < 1e-12;
      }
    }
    check(nFace > 0 && nEdge > 0, "the masked cell has face and edge contacts");
    check(faces, "face neighbours: edge flag 2, contamination above 0");
    check(edges, "edge and corner neighbours: edge flag 0, contamination 0");
    check(voids.nBasins == g.cells.nCells - 1, "the masked cell is in no basin");
  }

  // The edge flag is 2 when a core cell touches the boundary and 1 when only cells outside the core do.
  group("edge flag 1 or 2 by whether the core touches the boundary");
  {
    LightconeGrid g = built_grid(false);
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    DivergenceGrid& cells = g.cells;
    std::int64_t border = -1, inner = -1;
    for (std::size_t c = 0; c < cells.nCells && inner < 0; ++c) {
      if (!on_border(g, (std::int64_t)c)) continue;
      for (std::int64_t k = a.contactFirst[c]; k < a.contactFirst[c + 1]; ++k) {
        const std::int64_t n = a.contactNeighbour[(std::size_t)k];
        if (!on_border(g, n) && internal::lightcone_shared_area(cells, c, (std::size_t)n) > 0.) {
          border = (std::int64_t)c;
          inner = n;
          break;
        }
      }
    }
    check(inner >= 0, "a border cell with an inner face neighbour");
    for (const double borderValue : {0.5, -0.5}) {
      LightconeBasins b;
      cells.divergence.assign(cells.nCells, 1.);
      cells.divergence[(std::size_t)inner] = -1.;
      cells.divergence[(std::size_t)border] = borderValue;
      cells.voidIndex.assign(cells.nCells, -1);
      cells.overdenseIndex.assign(cells.nCells, -1);
      cells.voidIndex[(std::size_t)inner] = 0;
      cells.voidIndex[(std::size_t)border] = 0;
      b.voids.seed = {inner};
      const std::int64_t rest = std::min<std::int64_t>(border, inner) == 0 ? std::max<std::int64_t>(border, inner) + 1 : 0;
      for (std::size_t c = 0; c < cells.nCells; ++c)
        if (cells.voidIndex[c] < 0) cells.overdenseIndex[c] = 0;
      b.overdense.seed = {rest};
      BasinCatalog voids, overdense;
      internal::lightcone_basins(g, a, b, voids, overdense);
      check(voids.edgeFlag[0] == (borderValue < 0. ? 2 : 1),
            borderValue < 0. ? "border cell in the core: 2" : "border cell outside the core: 1");
      check(b.voids.nCells[0] == 2 && b.voids.nCoreCells[0] == (borderValue < 0. ? 2u : 1u), "cells and core cells");
    }
  }

  // Counts, volumes and radii are the sums over each basin's cells in ascending cell index; the centres and densities are left for later steps.
  group("per-basin counts, volumes and radii on a built grid");
  {
    LightconeGrid g = built_grid(false);
    DivergenceGrid& cells = g.cells;
    for (std::size_t c = 0; c < cells.nCells; c += 17) cells.masked[c] = 1;
    cells.divergence.assign(cells.nCells, 0.);
    cells.nTracers.assign(cells.nCells, 0);
    for (std::size_t c = 0; c < cells.nCells; ++c) {
      const double* x = &cells.centre[3*c];
      cells.divergence[c] = otswap::maths::sin(0.021 * x[0]) * otswap::maths::cos(0.017 * x[1]) +
                            otswap::maths::sin(0.013 * x[2] + 0.5) +
                            0.8 * otswap::maths::sin(0.11 * x[0]) * otswap::maths::sin(0.13 * x[1]) + 1.0;
      cells.nTracers[c] = (std::uint64_t)(c % 5);
    }
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    LightconeBasins b = internal::lightcone_watershed(a, cells);
    BasinCatalog voids, overdense;
    internal::lightcone_basins(g, a, b, voids, overdense);
    bool sums = true;
    for (int kind = 0; kind < 2; ++kind) {
      const bool isVoid = kind == 0;
      const BasinCatalog& cat = isVoid ? voids : overdense;
      const internal::BasinSet& set = isVoid ? b.voids : b.overdense;
      const std::vector<double>& core = isVoid ? b.voidCoreVolume : b.overdenseCoreVolume;
      const std::vector<double>& ext = isVoid ? b.voidVolume : b.overdenseVolume;
      for (std::size_t k = 0; k < cat.nBasins; ++k) {
        double vCore = 0., vExt = 0.;
        std::uint64_t n = 0, nCore = 0, t = 0;
        for (std::size_t c = 0; c < cells.nCells; ++c) {
          if ((isVoid ? cells.voidIndex[c] : cells.overdenseIndex[c]) != (std::int64_t)k) continue;
          const bool inCore = isVoid ? cells.divergence[c] < 0. : cells.divergence[c] >= 0.;
          ++n;
          vExt += cells.volume[c];
          if (inCore) {
            ++nCore;
            vCore += cells.volume[c];
            t += cells.nTracers[c];
          }
        }
        sums = sums && set.nCells[k] == n && set.nCoreCells[k] == nCore && set.nCoreTracers[k] == t &&
               cat.nTracers[k] == t && core[k] == vCore && ext[k] == vExt &&
               cat.radius[k] == otswap::maths::cbrt(3. * vCore / (4. * kPi)) &&
               cat.radiusExt[k] == otswap::maths::cbrt(3. * vExt / (4. * kPi));
      }
    }
    check(voids.nBasins > 1 && overdense.nBasins > 1, "several basins of each kind");
    check(sums, "every count, volume and radius equals the serial sum");
    check(voids.centre.empty() && voids.centreSky.empty() && voids.deltaAv.empty() && voids.delta2Mps.empty(),
          "centres and densities left empty");
    cells.nTracers.clear();
    internal::lightcone_basins(g, a, b, voids, overdense);
    check(voids.nTracers.empty() && b.voids.nCoreTracers.empty(), "no tracer counts without the grid's tracers");
  }

  // A search whose divergence is 0 or positive in every unmasked cell finds no void and raises an error.
  group("no void: an error");
  {
    const std::vector<double> tracersSky {1., 0.2, 0.3,  2., -0.2, 0.35,  3., 0.1, 0.4,  4., 0.3, 0.25};
    const std::vector<double> displacement(12, 0.);
    const std::vector<double> randomsSky {1.5, 0.1, 0.3,  2.5, 0., 0.3,  3.5, 0.2, 0.35,  4.5, -0.1, 0.4,
                                          5., 0.2, 0.3,  6., 0.1, 0.25,  0.5, 0.3, 0.3,  0.2, -0.3, 0.35};
    const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    otswap::MpsProfile profile;
    profile.redshift = {0.1, 0.5};
    profile.mps = {5., 6.};
    Config config;
    config.divergenceCellSize = 15.;
    throws_naming([&] { findVoidsLightcone(tracersSky, displacement, {}, randomsSky, mask, table(), profile, config); },
                  "no void found: the divergence is 0 or positive in every unmasked cell", "zero displacements");
  }

  // Counts, volumes, edge flags and contamination of a fixed grid are identical bit for bit with any number of threads and match the recorded hash.
  group("the basins are the same with 1 and 4 threads, and match the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = basins_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = basins_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "basins hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_basins");
}
