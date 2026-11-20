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
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <utility>
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
using internal::LightconeFill;
using internal::LightconeGrid;

namespace {

  constexpr std::uint64_t kRecordedHash = 0xdf9a98030dc2551aull;
  constexpr double kTwoPi = 2. * 3.14159265358979323846;

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  struct Ring {
    double decMin, decMax;
    bool full;
    std::vector<double> raEdges;
  };

  /// A grid laid out by hand: shells of rings, rings of sectors between consecutive RA edges.
  LightconeGrid layout (const std::vector<std::vector<Ring>>& shells)
  {
    LightconeGrid g;
    DivergenceGrid& cells = g.cells;
    cells.geometry = otswap::Geometry::Lightcone;
    g.ringFirstCell.push_back(0);
    for (const std::vector<Ring>& shell : shells) {
      g.shellFirstRing.push_back((std::int64_t)g.ringDecMin.size());
      for (const Ring& ring : shell) {
        g.ringDecMin.push_back(ring.decMin);
        g.ringDecMax.push_back(ring.decMax);
        g.ringFull.push_back(ring.full ? 1 : 0);
        g.ringPole.push_back(0);
        for (std::size_t k = 0; k + 1 < ring.raEdges.size(); ++k) {
          cells.cellRaMin.push_back(ring.raEdges[k]);
          cells.cellRaMax.push_back(ring.raEdges[k + 1]);
          cells.cellDecMin.push_back(ring.decMin);
          cells.cellDecMax.push_back(ring.decMax);
        }
        g.ringFirstCell.push_back((std::int64_t)cells.cellRaMin.size());
      }
    }
    g.shellFirstRing.push_back((std::int64_t)g.ringDecMin.size());
    cells.nCells = cells.cellRaMin.size();
    cells.masked.assign(cells.nCells, 0);
    cells.volume.assign(cells.nCells, 1.);
    for (std::size_t c = 0; c < cells.nCells; ++c) cells.centre.insert(cells.centre.end(), {100. * (double)c, 0., 0.});
    return g;
  }

  /// Put the cells of a laid-out grid at the given centres, with the given volumes.
  void place (LightconeGrid& g, const std::vector<std::array<double, 3>>& centres, const std::vector<double>& volumes)
  {
    for (std::size_t c = 0; c < centres.size(); ++c) {
      for (std::size_t a = 0; a < 3; ++a) g.cells.centre[3*c+a] = centres[c][a];
      g.cells.volume[c] = volumes[c];
    }
  }

  bool lists (const LightconeAdjacency& a, const std::int64_t c, const std::int64_t n)
  {
    return std::binary_search(a.contactNeighbour.begin() + a.contactFirst[(std::size_t)c],
                              a.contactNeighbour.begin() + a.contactFirst[(std::size_t)c + 1], n);
  }

  std::size_t size (const LightconeAdjacency& a, const std::int64_t c)
  {
    return (std::size_t)(a.contactFirst[(std::size_t)c + 1] - a.contactFirst[(std::size_t)c]);
  }

  bool radius_lists (const LightconeAdjacency& a, const std::int64_t c, const std::int64_t n)
  {
    return std::binary_search(a.radiusNeighbour.begin() + a.radiusFirst[(std::size_t)c],
                              a.radiusNeighbour.begin() + a.radiusFirst[(std::size_t)c + 1], n);
  }

  /// Every list is ascending, without repeats or the cell itself; masked cells have none; between unmasked
  /// cells the relation is symmetric.
  bool well_formed (const LightconeAdjacency& a, const DivergenceGrid& cells)
  {
    if (a.contactFirst.size() != cells.nCells + 1 || a.contactFirst[0] != 0 ||
        a.contactFirst.back() != (std::int64_t)a.contactNeighbour.size())
      return false;
    for (std::int64_t c = 0; c < (std::int64_t)cells.nCells; ++c) {
      if (cells.masked[(std::size_t)c] && size(a, c) != 0) return false;
      for (std::int64_t k = a.contactFirst[(std::size_t)c]; k < a.contactFirst[(std::size_t)c + 1]; ++k) {
        const std::int64_t n = a.contactNeighbour[(std::size_t)k];
        if (n == c || n < 0 || n >= (std::int64_t)cells.nCells) return false;
        if (k > a.contactFirst[(std::size_t)c] && !(a.contactNeighbour[(std::size_t)k - 1] < n)) return false;
        if (!cells.masked[(std::size_t)n] && !lists(a, n, c)) return false;
      }
    }
    return true;
  }

  /// Every radius list is ascending, without repeats or the cell itself, and holds unmasked cells only; masked
  /// cells have none; the relation is symmetric.
  bool radius_well_formed (const LightconeAdjacency& a, const DivergenceGrid& cells)
  {
    if (a.radiusFirst.size() != cells.nCells + 1 || a.radiusFirst[0] != 0 ||
        a.radiusFirst.back() != (std::int64_t)a.radiusNeighbour.size())
      return false;
    for (std::int64_t c = 0; c < (std::int64_t)cells.nCells; ++c) {
      if (cells.masked[(std::size_t)c] && a.radiusFirst[(std::size_t)c + 1] != a.radiusFirst[(std::size_t)c])
        return false;
      for (std::int64_t k = a.radiusFirst[(std::size_t)c]; k < a.radiusFirst[(std::size_t)c + 1]; ++k) {
        const std::int64_t n = a.radiusNeighbour[(std::size_t)k];
        if (n == c || n < 0 || n >= (std::int64_t)cells.nCells || cells.masked[(std::size_t)n]) return false;
        if (k > a.radiusFirst[(std::size_t)c] && !(a.radiusNeighbour[(std::size_t)k - 1] < n)) return false;
        if (!radius_lists(a, n, c)) return false;
      }
    }
    return true;
  }

  /// The radius lists by testing every pair of unmasked cells.
  bool radius_brute_force (const LightconeAdjacency& a, const DivergenceGrid& cells)
  {
    const double rNeigh = 0.5 * (std::sqrt(3.) + 2.);
    for (std::size_t c = 0; c < cells.nCells; ++c)
      for (std::size_t n = 0; n < cells.nCells; ++n) {
        if (c == n || cells.masked[c] || cells.masked[n]) continue;
        const double dx = cells.centre[3*c] - cells.centre[3*n], dy = cells.centre[3*c+1] - cells.centre[3*n+1],
                     dz = cells.centre[3*c+2] - cells.centre[3*n+2];
        const double limit = rNeigh * 0.5 * (otswap::maths::cbrt(cells.volume[c]) + otswap::maths::cbrt(cells.volume[n]));
        if ((dx*dx + dy*dy + dz*dz < limit * limit) != radius_lists(a, (std::int64_t)c, (std::int64_t)n)) return false;
      }
    return true;
  }

  /// A grid of cells given by their centres, with a symmetric adjacency built from pairs.
  struct Cells {

    DivergenceGrid grid;
    std::vector<std::pair<std::int64_t, std::int64_t>> pairs;

    std::int64_t add (const double x, const double y, const double z, const double divergence, const bool empty,
                      const bool masked = false)
    {
      grid.centre.insert(grid.centre.end(), {x, y, z});
      grid.divergence.push_back(divergence);
      grid.nDisplacements.push_back(empty ? 0 : 2);
      grid.weight.push_back(empty ? 0. : 2.);
      grid.masked.push_back(masked ? 1 : 0);
      return (std::int64_t)grid.nCells++;
    }

    void link (const std::int64_t a, const std::int64_t b)
    {
      pairs.push_back({a, b});
    }

    LightconeFill fill ()
    {
      LightconeAdjacency a;
      a.radiusFirst.assign(grid.nCells + 1, 0);
      std::vector<std::vector<std::int64_t>> list(grid.nCells);
      for (const auto& p : pairs) {
        if (!grid.masked[(std::size_t)p.first]) list[(std::size_t)p.first].push_back(p.second);
        if (!grid.masked[(std::size_t)p.second]) list[(std::size_t)p.second].push_back(p.first);
      }
      for (std::size_t c = 0; c < grid.nCells; ++c) {
        std::sort(list[c].begin(), list[c].end());
        a.radiusNeighbour.insert(a.radiusNeighbour.end(), list[c].begin(), list[c].end());
        a.radiusFirst[c + 1] = (std::int64_t)a.radiusNeighbour.size();
      }
      return internal::lightcone_fill(a, grid);
    }
  };

  LightconeGrid built_grid ()
  {
    std::mt19937 rng(1);
    auto uniform = [&rng] (const double lo, const double hi) { return lo + (hi - lo) * (rng() / 4294967296.); };
    std::vector<double> sky;
    for (std::size_t i = 0; i < 4000; ++i) {
      sky.push_back(uniform(0.5, 2.));
      sky.push_back(otswap::maths::asin(uniform(otswap::maths::sin(-0.3), otswap::maths::sin(0.3))));
      sky.push_back(uniform(0.25, 0.4));
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
    record.cellSizeMps = 8.;
    record.mpsProfile.redshift = {0.1, 0.5};
    record.mpsProfile.mps = {5., 6.};
    const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    LightconeGrid g = internal::lightcone_grid(record, table(), mask);
    for (std::size_t c = 0; c < g.cells.nCells; c += 23) g.cells.masked[c] = 1;
    return g;
  }

  void hash_bytes (std::uint64_t& h, const void* data, const std::size_t size)
  {
    const unsigned char* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      h ^= p[i];
      h *= 0x100000001b3ull;
    }
  }

  std::uint64_t fill_hash ()
  {
    LightconeGrid g = built_grid();
    DivergenceGrid& cells = g.cells;
    std::mt19937 rng(7);
    cells.nDisplacements.assign(cells.nCells, 0);
    cells.divergence.assign(cells.nCells, 0.);
    cells.weight.assign(cells.nCells, 0.);
    for (std::size_t c = 0; c < cells.nCells; ++c)
      if (!cells.masked[c] && rng() % 5 < 2) {
        cells.nDisplacements[c] = 1 + rng() % 4;
        cells.divergence[c] = -1. + 2. * (rng() / 4294967296.);
        cells.weight[c] = (double)cells.nDisplacements[c];
      }
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    const LightconeFill f = internal::lightcone_fill(a, cells);
    std::uint64_t h = 0xcbf29ce484222325ull;
    hash_bytes(h, a.contactFirst.data(), a.contactFirst.size() * sizeof(std::int64_t));
    hash_bytes(h, a.contactNeighbour.data(), a.contactNeighbour.size() * sizeof(std::int64_t));
    hash_bytes(h, a.radiusFirst.data(), a.radiusFirst.size() * sizeof(std::int64_t));
    hash_bytes(h, a.radiusNeighbour.data(), a.radiusNeighbour.size() * sizeof(std::int64_t));
    hash_bytes(h, cells.divergence.data(), cells.divergence.size() * sizeof(double));
    hash_bytes(h, cells.weight.data(), cells.weight.size() * sizeof(double));
    hash_bytes(h, &f.filled, sizeof f.filled);
    hash_bytes(h, &f.leftEmpty, sizeof f.leftEmpty);
    hash_bytes(h, &f.passes, sizeof f.passes);
    return h;
  }

  std::string run_log (const otswap::Verbosity verbosity)
  {
    const std::vector<double> tracersSky {1., 0.2, 0.3,  2., -0.2, 0.35,  3., 0.1, 0.4,  4., 0.3, 0.25};
    const std::vector<double> displacement(12, 0.1);
    const std::vector<double> randomsSky {1.5, 0.1, 0.3,  2.5, 0., 0.3,  3.5, 0.2, 0.35,  4.5, -0.1, 0.4,
                                          5., 0.2, 0.3,  6., 0.1, 0.25,  0.5, 0.3, 0.3,  0.2, -0.3, 0.35};
    const otswap::Mask mask(std::vector<double>(12, 1.), otswap::PixelOrdering::Ring);
    otswap::MpsProfile profile;
    profile.redshift = {0.1, 0.5};
    profile.mps = {5., 6.};
    Config config;
    std::ostringstream log;
    config.divergenceCellSize = 15.;
    config.verbosity = verbosity;
    config.log = &log;
    throws_naming([&] { findVoidsLightcone(tracersSky, displacement, {}, randomsSky, mask, table(), profile, config); },
                  "not implemented yet", "run with a log");
    return log.str();
  }

}

int main ()
{
  // A masked cell whose window reaches an unmasked cell is in that cell's list even when the unmasked cell's own window does not reach it; the masked cell keeps no list.
  group("a masked cell reached only by its own window is listed by its unmasked neighbour");
  {
    LightconeGrid g = layout({{{0., 0.1, false, {0., 1.}}},
                              {{0., 0.1, false, {0., 0.5, 1.005, 1.5}}}});
    g.cells.masked[0] = 1;
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    check(lists(a, 3, 0), "the unmasked cell lists the masked one");
    check(size(a, 0) == 0, "the masked cell has no list");
    check(well_formed(a, g.cells), "lists sorted, symmetric among unmasked cells");
  }

  // A cell of a full ring at RA = 0 and a cell of a partial ring lying just past 2 pi are neighbours both ways, though only the partial cell's window wraps.
  group("a full ring at RA = 0 and a partial ring past 2 pi list each other");
  {
    const double q = kTwoPi / 4.;
    LightconeGrid g = layout({{{0., 0.1, false, {5., 5.5, 6., 6.3}}},
                              {{0., 0.1, true, {0., q, 2. * q, 3. * q, kTwoPi}}}});
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    check(lists(a, 2, 3), "the partial cell lists the full ring's first cell");
    check(lists(a, 3, 2), "the full ring's first cell lists the partial cell");
    check(lists(a, 2, 6) && lists(a, 6, 2), "the partial cell and the full ring's last cell");
    check(!lists(a, 4, 2) && !lists(a, 0, 3), "cells far apart in RA are not neighbours");
    check(well_formed(a, g.cells), "lists sorted, symmetric among unmasked cells");
  }

  // The 1 % tolerance of a wide cell reaches a narrow cell in the next shell whose own tolerance does not reach back; the two are still neighbours both ways.
  group("a wide cell's tolerance makes it a neighbour of a narrow cell both ways");
  {
    LightconeGrid g = layout({{{0., 0.1, false, {0., 1.}}},
                              {{0., 0.1, false, {0., 0.5, 1.005, 1.5}}},
                              {{0.1015, 0.5, false, {0., 1.}}}});
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    check(lists(a, 0, 3) && lists(a, 3, 0), "across RA: wide and narrow list each other");
    check(lists(a, 1, 4) && lists(a, 4, 1), "across declination: wide and narrow list each other");
    check(!lists(a, 0, 4), "cells two shells apart are not neighbours");
    check(lists(a, 0, 1) && lists(a, 0, 2), "cells sharing a face are neighbours");
    check(well_formed(a, g.cells), "lists sorted, symmetric among unmasked cells");
  }

  // On a grid as built, every list is sorted, holds no cell twice nor the cell itself, and the relation is symmetric among unmasked cells; masked cells keep no list but appear in their neighbours'.
  group("the adjacency of a built grid is well formed and symmetric");
  {
    const LightconeGrid g = built_grid();
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    check(well_formed(a, g.cells), "lists sorted, symmetric among unmasked cells");
    bool maskedListed = false;
    for (const std::int64_t n : a.contactNeighbour) maskedListed = maskedListed || g.cells.masked[(std::size_t)n];
    check(maskedListed, "masked cells appear in the lists");
  }

  // Two unmasked cells are radius neighbours when their centres are closer than (sqrt(3) + 2) / 2 times the mean of their sides; the bound is strict, scales with the sides, and masked cells are in no list.
  group("radius neighbours follow the distance bound, the cells' sides and the mask");
  {
    const double rNeigh = 0.5 * (std::sqrt(3.) + 2.);
    LightconeGrid g = layout({{{0., 0.1, false, {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}}}});
    place(g, {{0., 0., 0.}, {1.8, 0., 0.}, {0., 1.9, 0.}, {0., 0., 2.7}, {0.5, 0., 0.}, {0., 20., 0.}, {rNeigh, 20., 0.}},
          {1., 1., 1., 8., 1., 1., 1.});
    g.cells.masked[4] = 1;
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    check(radius_lists(a, 0, 1) && radius_lists(a, 1, 0), "closer than the bound: neighbours both ways");
    check(!radius_lists(a, 0, 2), "farther than the bound: not neighbours");
    check(radius_lists(a, 0, 3) && radius_lists(a, 3, 0), "a cell of side 2 reaches farther");
    check(!radius_lists(a, 2, 3), "two cells beyond the bound of their sides");
    check(a.radiusFirst[4] == a.radiusFirst[5], "a masked cell has no list");
    check(!radius_lists(a, 0, 4) && !radius_lists(a, 1, 4), "a masked cell is in no list");
    check(!radius_lists(a, 5, 6), "a distance equal to the bound is not a neighbour");
    check(radius_well_formed(a, g.cells), "lists sorted, symmetric, unmasked cells only");
  }

  // On a grid as built, with masked cells, the radius lists are exactly the pairs closer than the bound.
  group("the radius lists of a built grid are those of a search over every pair");
  {
    const LightconeGrid g = built_grid();
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    check(radius_well_formed(a, g.cells), "lists sorted, symmetric, unmasked cells only");
    check(radius_brute_force(a, g.cells), "lists equal the brute-force pairs");
    check(a.radiusNeighbour.size() > 4 * g.cells.nCells, "cells have several neighbours");
  }

  // A grid whose cells are all masked has empty lists in both definitions, and no error.
  group("a grid with every cell masked has no neighbours");
  {
    LightconeGrid g = built_grid();
    std::fill(g.cells.masked.begin(), g.cells.masked.end(), 1);
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    check(a.radiusNeighbour.empty() && a.contactNeighbour.empty(), "no list holds a cell");
  }

  // An empty cell with at least three neighbours holding a value takes their mean weighted by 1 / r^2, and keeps a weight of 0; with two it stays empty and ends at 0.
  group("three sources fill an empty cell with 1 / r^2 weights; two do not");
  {
    Cells c;
    const std::int64_t t = c.add(0., 0., 0., 9., true);
    c.link(t, c.add(1., 0., 0., 1., false));
    c.link(t, c.add(0., 1., 0., 2., false));
    c.link(t, c.add(0., 0., 2., 3., false));
    const std::int64_t u = c.add(10., 0., 0., 9., true);
    c.link(u, c.add(11., 0., 0., 1., false));
    c.link(u, c.add(10., 1., 0., 2., false));
    const LightconeFill f = c.fill();
    const double num = 1. * 1. + 2. * 1. + 3. * 0.25, den = 1. + 1. + 0.25;
    check(c.grid.divergence[(std::size_t)t] == num / den, "the 1 / r^2 mean");
    check(c.grid.weight[(std::size_t)t] == 0., "a filled cell keeps a weight of 0");
    check(c.grid.divergence[(std::size_t)u] == 0. && c.grid.weight[(std::size_t)u] == 0., "left empty: 0 and 0");
    check(c.grid.divergence[1] == 1. && c.grid.weight[1] == 2., "sources unchanged");
    check(f.filled == 1 && f.passes == 1 && f.leftEmpty == 1, "one filled in one pass, one left empty");
  }

  // Each pass reads only the values of the previous pass: a cell filled in a pass is a source from the next pass on.
  group("a cell filled in one pass is a source only in the next");
  {
    Cells c;
    const std::int64_t a = c.add(0., 0., 0., 0., true);
    c.link(a, c.add(1., 0., 0., 1., false));
    c.link(a, c.add(-1., 0., 0., 1., false));
    c.link(a, c.add(0., 1., 0., 1., false));
    const std::int64_t b = c.add(0., -1., 0., 0., true);
    c.link(b, a);
    c.link(b, c.add(1., -1., 0., 4., false));
    c.link(b, c.add(-1., -1., 0., 4., false));
    const LightconeFill f = c.fill();
    check(c.grid.divergence[(std::size_t)a] == 1., "the first cell in the first pass");
    const double expected = (1. * 1. + 4. * 1. + 4. * 1.) / 3.;
    check(c.grid.divergence[(std::size_t)b] == expected, "the second cell from the first cell's filled value");
    check(f.filled == 2 && f.passes == 2 && f.leftEmpty == 0, "two cells in two passes");
  }

  // Masked cells are never sources, whatever they hold, and are never filled.
  group("masked cells are neither sources nor filled");
  {
    Cells c;
    const std::int64_t t = c.add(0., 0., 0., 0., true);
    c.link(t, c.add(1., 0., 0., 1., false));
    c.link(t, c.add(-1., 0., 0., 1., false));
    c.link(t, c.add(0., 1., 0., 5., false, true));
    const std::int64_t m = c.add(0., -1., 0., 0., true, true);
    c.link(t, m);
    c.link(m, c.add(0., -2., 0., 1., false));
    const LightconeFill f = c.fill();
    check(c.grid.divergence[(std::size_t)t] == 0., "two unmasked sources do not fill");
    check(c.grid.divergence[(std::size_t)m] == 0., "a masked empty cell is not filled");
    check(f.filled == 0 && f.passes == 0 && f.leftEmpty == 1, "nothing filled, one left empty");
  }

  // With no empty cell the fill runs no pass and changes nothing.
  group("a grid with no empty cell is unchanged");
  {
    Cells c;
    const std::int64_t a = c.add(0., 0., 0., 0.5, false);
    c.link(a, c.add(1., 0., 0., -0.5, false));
    const LightconeFill f = c.fill();
    check(c.grid.divergence[0] == 0.5 && c.grid.divergence[1] == -0.5, "values unchanged");
    check(f.filled == 0 && f.passes == 0 && f.leftEmpty == 0, "no pass");
  }

  // The report gives the fill line at Detailed, and at Normal only when cells are left empty.
  group("the fill line is written at Detailed, and at Normal when cells are left empty");
  {
    const std::string detailed = run_log(otswap::Verbosity::Detailed);
    const std::size_t at = detailed.find("bitvf: filled ");
    check(at != std::string::npos, "Detailed: " + detailed);
    if (at != std::string::npos) {
      const std::string line = detailed.substr(at, detailed.find('\n', at) - at);
      const std::size_t comma = line.find(", ");
      const std::size_t leftEmpty = std::stoul(line.substr(comma + 2));
      check(line.find(" left empty in ") != std::string::npos, "the line: " + line);
      const std::string normal = run_log(otswap::Verbosity::Normal);
      check((normal.find("bitvf: filled ") != std::string::npos) == (leftEmpty > 0),
            "Normal writes the line exactly when cells are left empty: " + normal);
    }
  }

  // The lists and the fill of a fixed grid are identical bit for bit with any number of threads and match the recorded hash.
  group("the lists and the fill are the same with 1 and 4 threads, and match the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = fill_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = fill_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "fill hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_fill");
}
