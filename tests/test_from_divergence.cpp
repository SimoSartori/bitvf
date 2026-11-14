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
 *  @file tests/test_from_divergence.cpp
 *
 *  @brief findVoidsBoxFromDivergence: the search it repeats, the records
 *  it copies, the tracers it checks, and its errors.
 */

#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include "bitvf/io.h"
#include "check.h"
#include "fixtures.h"

using namespace bitvf;
using namespace fixtures;

namespace {

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "bitvf_test_from_divergence_" + name;
  }

  /// The assignment of two searches equal, bit for bit, but used.
  bool same_assignment (const TracerAssignment& a, const TracerAssignment& b)
  {
    return a.cellIndex == b.cellIndex && a.voidIndex == b.voidIndex && a.overdenseIndex == b.overdenseIndex &&
           same_bits(a.distance, b.distance) && a.inMask == b.inMask && a.inBox == b.inBox;
  }

  /// The search from a divergence equal to the search that produced it.
  bool repeats (const Result& f, const Result& r)
  {
    const DivergenceGrid& g = *f.divergence;
    const DivergenceGrid& d = *r.divergence;
    return same_catalog(f.voids, r.voids) && same_catalog(f.overdense, r.overdense) &&
           same_bits(g.divergence, d.divergence) && same_bits(g.divergenceSmoothed, d.divergenceSmoothed) &&
           same_bits(g.weight, d.weight) && same_bits(g.weightSmoothed, d.weightSmoothed) && g.voidIndex == d.voidIndex && g.overdenseIndex == d.overdenseIndex &&
           g.nTracers == d.nTracers && g.nDisplacements == d.nDisplacements && f.assignment && r.assignment &&
           same_assignment(*f.assignment, *r.assignment) && f.ids == r.ids;
  }

}

int main ()
{
  set_threads(4);
  const Translated in;
  const Result r = in.run();
  const DivergenceGrid& d = *r.divergence;

  group("the fixture finds voids and overdense basins");
  check(r.voids.nBasins > 0 && r.overdense.nBasins > 0, "both kinds of basin");

  group("from the Result, and from its divergence file in FITS and ASCII: the catalogues, the grid's per-cell "
        "fields and the assignment but used equal the search's, bit for bit, with 1 and 4 threads");
  for (const char* source : {"memory", ".fits", ".txt"}) {
    DivergenceGrid grid = d;
    const std::string where(source);
    if (where != "memory") {
      const std::string file = temporary("divergence" + where);
      io::writeDivergence(file, r);
      grid = io::readDivergence(file);
      std::remove(file.c_str());
    }
    for (const int t : {1, 4}) {
      set_threads(t);
      const Result f = findVoidsBoxFromDivergence(grid, in.tracers, in.config, in.ids);
      check(repeats(f, r), where + ", " + std::to_string(t) + " thread(s)");
    }
    set_threads(4);
  }

  group("records: NREC, NRECMEAN, nDisplacements, limits, mps and cell size copied from the grid; used empty, "
        "no displacement counted, the tracers outside the limits and the grid counted as by the search");
  {
    const Result f = findVoidsBoxFromDivergence(d, in.tracers, in.config, in.ids);
    check(f.nRealizations == 2 && d.nRealizations == 2 && f.divergence->nRealizations == 2,
          "NREC: the grid's 2 realizations");
    check(f.config.nRealizationsMean == 3 && f.divergence->nRealizationsMean == 3, "NRECMEAN: the grid's 3");
    check(f.divergence->nDisplacements == d.nDisplacements, "nDisplacements copied");
    check(f.limits && d.limits && f.limits->lower == d.limits->lower && f.limits->upper == d.limits->upper,
          "the grid's limits");
    check(f.mps == Translated::kMps && f.config.divergenceCellSize && *f.config.divergenceCellSize == d.cellSizeMps,
          "the grid's mps and cell size");
    check(f.assignment->used.empty(), "used empty");
    check(f.counts.displacements == 0 && f.counts.tracersUsed == 0, "no displacement, no tracer used");
    check(f.counts.tracers == r.counts.tracers && f.counts.tracersOutsideBox == Translated::kOutside &&
          f.counts.tracersOutsideBox == r.counts.tracersOutsideBox &&
          f.counts.tracersOutsideGrid == r.counts.tracersOutsideGrid, "tracers given, outside the limits, the grid");
    check(f.nTracers == r.nTracers && f.geometry == otswap::Geometry::Box, "nTracers and geometry");
    check(f.outsideGrid == std::vector<std::uint8_t>(r.nTracers, 0) && r.outsideGrid.empty(),
          "outsideGrid set, as with a reused grid, none outside");
  }

  group("another sigma: equal to a search with that sigma, sigma 0 included; the divergence file of a search "
        "at sigma 0, searched with the first sigma, equal to the first search");
  for (const double sigma : {0., 1.}) {
    Config c = in.config;
    c.smoothing.sigma = sigma;
    const Result s = in.run(c);
    const Result f = findVoidsBoxFromDivergence(d, in.tracers, c, in.ids);
    check(repeats(f, s) && s.overdense.nBasins > 0, "sigma " + std::to_string(sigma) + " mps");
  }
  {
    Config c = in.config;
    c.smoothing.sigma = 0.;
    const std::string file = temporary("sigma0.fits");
    io::writeDivergence(file, in.run(c));
    const DivergenceGrid zero = io::readDivergence(file);
    std::remove(file.c_str());
    check(zero.divergenceSmoothed.empty() && zero.weightSmoothed.empty() &&
          repeats(findVoidsBoxFromDivergence(zero, in.tracers, in.config, in.ids), r), "from sigma 0 to the first");
  }

  group("the tracers against the grid's nTracers: a tracer moved to another cell or dropped is refused; "
        "without nTracers the same tracers are read");
  {
    std::vector<double> moved = in.tracers;
    moved[0] += 2. * d.cellSize;
    throws_naming([&] { findVoidsBoxFromDivergence(d, moved, in.config); },
                  "the tracers do not match the divergence: cell ", "a tracer moved");
    const std::vector<double> dropped(in.tracers.begin() + 3, in.tracers.end());
    throws_naming([&] { findVoidsBoxFromDivergence(d, dropped, in.config); },
                  "pass the tracers of the search that produced the field", "a tracer dropped");
    DivergenceGrid bare = d;
    bare.nTracers.clear();
    const Result f = findVoidsBoxFromDivergence(bare, moved, in.config);
    check(f.voids.nBasins > 0 && f.divergence->nTracers.size() == d.nCells, "no nTracers: accepted, then counted");
    check(same_catalog(findVoidsBoxFromDivergence(bare, in.tracers, in.config, in.ids).voids, r.voids),
          "no nTracers: the same catalogues with the same tracers");
  }

  group("errors of a search from a divergence grid");
  {
    const Config& c = in.config;
    throws_naming([&] { findVoidsBoxFromDivergence(d, {}, c); }, "tracers is empty", "no tracer");
    throws_naming([&] { findVoidsBoxFromDivergence(d, {1., 2.}, c); }, "tracers", "size");
    throws_naming([&] { findVoidsBoxFromDivergence(d, in.tracers, c, {1, 2}); }, "ids", "ids");
    throws_naming([&] { Config x = c; x.smoothing.sigma = -1.; findVoidsBoxFromDivergence(d, in.tracers, x); },
                  "sigma", "config");
    throws_naming([&] { Config x = c; x.divergenceCellSize = 1.; findVoidsBoxFromDivergence(d, in.tracers, x); },
                  "the divergence cell size is set for a search from a divergence", "cell size");
    throws_naming([&] { DivergenceGrid x = d; x.nCellsAxis[2] -= 1; findVoidsBoxFromDivergence(x, in.tracers, c); },
                  "findVoidsBoxFromDivergence: the divergence grid has ", "non-cubic");
    throws_naming([&] { DivergenceGrid x = d; x.geometry = otswap::Geometry::Lightcone;
                        findVoidsBoxFromDivergence(x, in.tracers, c); },
                  "the divergence grid is not a box grid", "lightcone");
    throws_naming([&] { DivergenceGrid x = d; x.limits->upper[1] = x.limits->lower[1];
                        findVoidsBoxFromDivergence(x, in.tracers, c); },
                  "the grid's box limits are ", "the grid's box limits are ");
    throws_naming([&] { DivergenceGrid x = d; x.divergence.clear();
                        findVoidsBoxFromDivergence(x, in.tracers, c); },
                  "the grid's divergence holds 0 entries", "no divergence");
    throws_naming([&] { DivergenceGrid x = d; x.weight[7] = std::numeric_limits<double>::quiet_NaN();
                        findVoidsBoxFromDivergence(x, in.tracers, c); },
                  "the grid's weight holds a non-finite value at cell 7", "weight");
    throws_naming([&] { DivergenceGrid x = d; x.nTracers.pop_back(); findVoidsBoxFromDivergence(x, in.tracers, c); },
                  "the grid's nTracers holds ", "the grid's nTracers holds ");
    throws_naming([&] { Config x = c; x.smoothing.sigma = 1e-200; findVoidsBoxFromDivergence(d, in.tracers, x); },
                  "too small to build the kernel", "too small to build the kernel");
    throws_naming([&] { findVoidsBoxFromDivergence(d, {200., 0., 0., 0., 200., 0.}, c); },
                  "no tracer is left to read: of 2 tracers, 2 lie outside the grid's box limits", "no tracer is left to read: of 2 tracers, 2 lie outside the grid's box limits");
    throws_naming([&] { DivergenceGrid x = d; for (double& v : x.divergence) v = 1.;
                        findVoidsBoxFromDivergence(x, in.tracers, c); },
                  "no void found: the divergence is 0 or positive in every cell", "no void found: the divergence is 0 or positive in every cell");
  }

  return report("test_from_divergence");
}
