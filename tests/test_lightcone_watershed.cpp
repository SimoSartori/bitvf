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
#include <sstream>
#include <random>
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
using internal::LightconeBasins;
using internal::LightconeGrid;

namespace {

  constexpr std::uint64_t kRecordedHash = 0x85ebc3ed24cc8904ull;

  const otswap::DistanceTable& table ()
  {
    static const otswap::DistanceTable t(0.3, 0.7, -1., 0., 0., 2., 1000);
    return t;
  }

  /// Cells with a value, a mask flag and radius neighbours given by pairs; pairs touching a masked cell are dropped.
  struct Field {

    DivergenceGrid grid;
    std::vector<std::pair<std::int64_t, std::int64_t>> pairs;

    std::int64_t add (const double value, const bool masked = false)
    {
      grid.divergence.push_back(value);
      grid.masked.push_back(masked ? 1 : 0);
      return (std::int64_t)grid.nCells++;
    }

    void link (const std::int64_t a, const std::int64_t b)
    {
      pairs.push_back({a, b});
    }

    LightconeBasins run ()
    {
      LightconeAdjacency a;
      std::vector<std::vector<std::int64_t>> list(grid.nCells);
      for (const auto& p : pairs)
        if (!grid.masked[(std::size_t)p.first] && !grid.masked[(std::size_t)p.second]) {
          list[(std::size_t)p.first].push_back(p.second);
          list[(std::size_t)p.second].push_back(p.first);
        }
      a.radiusFirst.assign(grid.nCells + 1, 0);
      for (std::size_t c = 0; c < grid.nCells; ++c) {
        std::sort(list[c].begin(), list[c].end());
        a.radiusNeighbour.insert(a.radiusNeighbour.end(), list[c].begin(), list[c].end());
        a.radiusFirst[c + 1] = (std::int64_t)a.radiusNeighbour.size();
      }
      return internal::lightcone_watershed(a, grid);
    }
  };

  bool labels (const DivergenceGrid& g, const std::vector<std::int64_t>& voids, const std::vector<std::int64_t>& overdense)
  {
    return g.voidIndex == voids && g.overdenseIndex == overdense;
  }

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
    DivergenceGrid& cells = g.cells;
    cells.divergence.assign(cells.nCells, 0.);
    for (std::size_t c = 0; c < cells.nCells; ++c) {
      const double* x = &cells.centre[3*c];
      cells.divergence[c] = otswap::maths::sin(0.021 * x[0]) * otswap::maths::cos(0.017 * x[1]) +
                            otswap::maths::sin(0.013 * x[2] + 0.5) +
                            0.8 * otswap::maths::sin(0.11 * x[0]) * otswap::maths::sin(0.13 * x[1]) + 0.3;
    }
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

  std::uint64_t watershed_hash ()
  {
    LightconeGrid g = built_grid();
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    const LightconeBasins b = internal::lightcone_watershed(a, g.cells);
    std::uint64_t h = 0xcbf29ce484222325ull;
    hash_bytes(h, g.cells.voidIndex.data(), g.cells.voidIndex.size() * sizeof(std::int64_t));
    hash_bytes(h, g.cells.overdenseIndex.data(), g.cells.overdenseIndex.size() * sizeof(std::int64_t));
    hash_bytes(h, b.voids.seed.data(), b.voids.seed.size() * sizeof(std::int64_t));
    hash_bytes(h, b.overdense.seed.data(), b.overdense.seed.size() * sizeof(std::int64_t));
    return h;
  }

}

int main ()
{
  // A cell flows to its first neighbour below it and below every earlier neighbour; a cell with none is a seed, a void below 0 and an overdense basin otherwise, each numbered by ascending seed cell.
  group("strict descent on a line of cells: seeds, kinds and labels");
  {
    Field f;
    for (const double v : {2., 1., 3., -1., 0.5}) f.add(v);
    for (std::int64_t c = 0; c < 4; ++c) f.link(c, c + 1);
    const LightconeBasins b = f.run();
    check(b.voids.seed == std::vector<std::int64_t>{3}, "one void, seeded at its minimum");
    check(b.overdense.seed == std::vector<std::int64_t>{1}, "one overdense basin, seeded at its minimum");
    check(labels(f.grid, {-1, -1, 0, 0, 0}, {0, 0, -1, -1, -1}), "every cell in the basin its flow reaches");
    check(b.voids.nCells.empty() && b.overdense.nCoreCells.empty(), "the counts are left empty");
  }

  // Of two neighbours of equal smallest value the first in ascending index wins, so a tie goes to the lower cell.
  group("a tie between two smaller neighbours goes to the lower index");
  {
    Field f;
    for (const double v : {-1., 5., 0., 5., -1.}) f.add(v);
    f.link(0, 1);
    f.link(1, 2);
    f.link(2, 3);
    f.link(3, 4);
    f.link(0, 2);
    f.link(2, 4);
    const LightconeBasins b = f.run();
    check(b.voids.seed == (std::vector<std::int64_t>{0, 4}), "two voids, numbered by seed cell");
    check(labels(f.grid, {0, 0, 0, 1, 1}, {-1, -1, -1, -1, -1}), "the middle cell flows to the lower of the equal minima");
  }

  // Two neighbouring cells of equal value, with nothing smaller, are both seeds.
  group("a plateau gives one seed per cell");
  {
    Field f;
    f.add(1.);
    f.add(1.);
    f.link(0, 1);
    const LightconeBasins b = f.run();
    check(b.overdense.seed == (std::vector<std::int64_t>{0, 1}) && b.voids.seed.empty(), "two overdense seeds, no void");
    check(labels(f.grid, {-1, -1}, {0, 1}), "each cell its own basin");
  }

  // A masked cell is never a seed, a source or a target, and takes -1 in both labels; a cell whose only lower neighbour is masked is a seed.
  group("masked cells belong to no basin and never attract a flow");
  {
    Field f;
    f.add(1.);
    f.add(-10., true);
    f.add(2.);
    f.add(-1.);
    f.link(0, 1);
    f.link(1, 2);
    f.link(0, 2);
    f.link(2, 3);
    const LightconeBasins b = f.run();
    check(b.voids.seed == std::vector<std::int64_t>{3}, "the masked minimum seeds nothing");
    check(b.overdense.seed == std::vector<std::int64_t>{0}, "a cell with only a masked lower neighbour is a seed");
    check(labels(f.grid, {-1, -1, 0, 0}, {0, -1, -1, -1}), "masked cell labelled -1 in both");
  }

  // The smoothed field is the one split when it is not empty, otherwise the raw divergence.
  group("the smoothed divergence is split when present");
  {
    Field f;
    for (const double v : {1., 2.}) f.add(v);
    f.link(0, 1);
    f.run();
    check(f.grid.voidIndex.size() == 2 && f.grid.overdenseIndex[1] == 0, "raw field: cell 1 flows to cell 0");
    f.grid.divergenceSmoothed = {-1., -2.};
    f.pairs.clear();
    f.link(0, 1);
    const LightconeBasins b = f.run();
    check(b.voids.seed == std::vector<std::int64_t>{1} && b.overdense.seed.empty(), "smoothed field: one void at cell 1");
    check(labels(f.grid, {0, 0}, {-1, -1}), "both cells in it");
  }

  // With no cell below 0 there is no void and no error.
  group("a field with no negative cell has no void");
  {
    Field f;
    for (const double v : {0., 1., 2.}) f.add(v);
    f.link(0, 1);
    f.link(1, 2);
    const LightconeBasins b = f.run();
    check(b.voids.seed.empty() && b.overdense.seed == std::vector<std::int64_t>{0}, "one overdense basin");
    check(labels(f.grid, {-1, -1, -1}, {0, 0, 0}), "all cells in it");
  }

  // On a built grid with masked cells and a smooth field, every cell reaches the seed that serial descent on the radius lists reaches, seeds ascend, and the kind follows the sign.
  group("a built grid: labels equal serial descent, seeds ascend, masked cells unlabelled");
  {
    LightconeGrid g = built_grid();
    const LightconeAdjacency a = internal::lightcone_adjacency(g);
    const LightconeBasins b = internal::lightcone_watershed(a, g.cells);
    const DivergenceGrid& cells = g.cells;
    bool ascending = std::is_sorted(b.voids.seed.begin(), b.voids.seed.end()) &&
                     std::is_sorted(b.overdense.seed.begin(), b.overdense.seed.end());
    bool kinds = true, serial = true, masked = true, exclusive = true;
    for (const std::int64_t s : b.voids.seed) kinds = kinds && cells.divergence[(std::size_t)s] < 0.;
    for (const std::int64_t s : b.overdense.seed) kinds = kinds && cells.divergence[(std::size_t)s] >= 0.;
    for (std::size_t c = 0; c < cells.nCells; ++c) {
      if (cells.masked[c]) {
        masked = masked && cells.voidIndex[c] == -1 && cells.overdenseIndex[c] == -1;
        continue;
      }
      std::int64_t cur = (std::int64_t)c;
      while (true) {
        std::int64_t next = cur;
        double minVal = cells.divergence[(std::size_t)cur];
        for (std::int64_t k = a.radiusFirst[(std::size_t)cur]; k < a.radiusFirst[(std::size_t)cur + 1]; ++k) {
          const std::int64_t n = a.radiusNeighbour[(std::size_t)k];
          if (cells.divergence[(std::size_t)n] < minVal) {
            minVal = cells.divergence[(std::size_t)n];
            next = n;
          }
        }
        if (next == cur) break;
        cur = next;
      }
      const bool isVoid = cells.divergence[(std::size_t)cur] < 0.;
      const std::vector<std::int64_t>& seeds = isVoid ? b.voids.seed : b.overdense.seed;
      const std::int64_t index = std::lower_bound(seeds.begin(), seeds.end(), cur) - seeds.begin();
      serial = serial && (isVoid ? cells.voidIndex[c] : cells.overdenseIndex[c]) == index &&
               (isVoid ? cells.overdenseIndex[c] : cells.voidIndex[c]) == -1;
      exclusive = exclusive && (cells.voidIndex[c] >= 0) != (cells.overdenseIndex[c] >= 0);
    }
    check(!b.voids.seed.empty() && !b.overdense.seed.empty(), "both kinds of basin exist");
    check(ascending && kinds, "seeds ascend and a seed's sign gives its kind");
    check(serial, "every label is the basin of the serial descent");
    check(exclusive, "each unmasked cell is in exactly one basin");
    check(masked, "masked cells are unlabelled");
  }

  // The search reports its watershed at Detailed, and nothing at Normal.
  group("the watershed line is written at Detailed only");
  {
    auto run_log = [] (const otswap::Verbosity verbosity) {
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
    };
    const std::string detailed = run_log(otswap::Verbosity::Detailed);
    check(detailed.find("bitvf: watershed into ") != std::string::npos && detailed.find(" overdense basin") != std::string::npos,
          "Detailed: " + detailed);
    check(run_log(otswap::Verbosity::Normal).find("watershed") == std::string::npos, "Normal writes no line");
  }

  // The labels and seeds of a fixed grid are identical bit for bit with any number of threads and match the recorded hash.
  group("the watershed is the same with 1 and 4 threads, and matches the recorded hash");
  {
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    const std::uint64_t one = watershed_hash();
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    const std::uint64_t four = watershed_hash();
    check(one == four, "1 and 4 threads give the same hash");
    if (one != kRecordedHash) std::cerr << "watershed hash: 0x" << std::hex << one << std::dec << std::endl;
    check(one == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_lightcone_watershed");
}
