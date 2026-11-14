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
 *  @file tests/test_from_field.cpp
 *
 *  @brief findVoidsBoxFromField: the lattice it infers, the search on a
 *  bitvf divergence file, fieldCentre, the field result, its errors and
 *  its determinism.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "bitvf/io.h"
#include "internal.h"
#include "check.h"
#include "fixtures.h"

using namespace bitvf;
using namespace fixtures;

namespace {

  // The hash of the search on the default LatticeField, with its weights,
  // sigma 1.5 in coordinate units, recorded once. A different value on
  // some platform means results are not portable there.
  constexpr std::uint64_t kRecordedHash = 0x06eb8de77478e0b9ull;

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "bitvf_test_from_field_" + name;
  }

  Config coordinate (const double sigma)
  {
    Config c;
    c.log = nullptr;
    c.keepDivergence = true;
    c.smoothing.sigma = sigma;
    c.smoothing.unit = LengthUnit::Coordinate;
    return c;
  }

  /// The catalogues and the grid of a field result.
  std::uint64_t hash (const Result& r)
  {
    std::uint64_t h = 0xcbf29ce484222325ull;
    for (const BasinCatalog* c : {&r.voids, &r.overdense}) {
      fnv(h, &c->nBasins, sizeof c->nBasins);
      fnv_vector(h, c->centre);
      fnv_vector(h, c->radius);
      fnv_vector(h, c->radiusExt);
      fnv_vector(h, c->edgeFlag);
      fnv_vector(h, c->contamination);
      fnv_vector(h, c->fieldCentre);
    }
    const DivergenceGrid& g = *r.divergence;
    fnv(h, &g.cellSize, sizeof g.cellSize);
    fnv(h, g.lower.data(), sizeof g.lower);
    fnv(h, g.upper.data(), sizeof g.upper);
    fnv(h, g.nCellsAxis.data(), sizeof g.nCellsAxis);
    fnv_vector(h, g.divergenceSmoothed);
    fnv_vector(h, g.weight);
    fnv_vector(h, g.voidIndex);
    fnv_vector(h, g.overdenseIndex);
    return h;
  }

  /// The catalogues and the grid of two field results equal, but weight
  /// and weightSmoothed.
  bool same_search (const Result& a, const Result& b)
  {
    const DivergenceGrid& g = *a.divergence;
    const DivergenceGrid& d = *b.divergence;
    return same_catalog(a.voids, b.voids) && same_catalog(a.overdense, b.overdense) &&
           same_bits(g.divergence, d.divergence) && same_bits(g.divergenceSmoothed, d.divergenceSmoothed) &&
           g.voidIndex == d.voidIndex && g.overdenseIndex == d.overdenseIndex && g.cellSize == d.cellSize &&
           g.lower == d.lower && g.upper == d.upper && g.nCellsAxis == d.nCellsAxis;
  }

  bool same_result (const Result& a, const Result& b)
  {
    return same_search(a, b) && same_bits(a.divergence->weight, b.divergence->weight) &&
           same_bits(a.divergence->weightSmoothed, b.divergence->weightSmoothed);
  }

  /// The rows of a field in the order of perm.
  void permute (const std::vector<std::size_t>& perm, std::vector<double>& centres, std::vector<double>& values,
                std::vector<double>& weights)
  {
    const std::vector<double> c = centres, v = values, w = weights;
    for (std::size_t i = 0; i < perm.size(); ++i) {
      for (std::size_t a = 0; a < 3; ++a) centres[3*i+a] = c[3*perm[i]+a];
      values[i] = v[perm[i]];
      if (!w.empty()) weights[i] = w[perm[i]];
    }
  }

  /// A grid of side 1 and lower (0, 0, 0) for box_field_at.
  DivergenceGrid unit_grid (const std::array<std::size_t, 3>& n)
  {
    DivergenceGrid g;
    g.cellSize = 1.;
    g.lower = {0., 0., 0.};
    g.nCellsAxis = n;
    g.nCells = n[0] * n[1] * n[2];
    return g;
  }

}

int main ()
{
  set_threads(4);
  std::mt19937 rng(11);

  group("the fixtures find voids and overdense basins");
  const Translated in;
  const Result r = in.run();
  const DivergenceGrid& d = *r.divergence;
  const LatticeField lf;
  const Result lr = findVoidsBoxFromField(lf.centres, lf.values, coordinate(1.5), lf.weights);
  check(r.voids.nBasins > 0 && r.overdense.nBasins > 0, "the search on Translated");
  check(lr.voids.nBasins > 0 && lr.overdense.nBasins > 0, "the search on LatticeField");

  group("a bitvf divergence file through readField, FITS and ASCII, with the search's sigma and mps: the "
        "search's lattice and partition; centres and radii within the rounding of the inferred grid");
  for (const char* format : {".fits", ".txt"}) {
    const std::string file = temporary(std::string("divergence") + format);
    io::writeDivergence(file, r);
    const io::Field f = io::readField(file, {"cellX", "cellY", "cellZ"}, "divergence", "weight");
    std::remove(file.c_str());
    Config c = in.config;
    c.keepAssignment = false;
    c.nRealizationsMean = 1;
    const Result ff = findVoidsBoxFromField(f.centres, f.values, c, f.weights, Translated::kMps);
    const DivergenceGrid& g = *ff.divergence;
    const std::string what = std::string(format) + ": ";
    check(g.nCellsAxis == d.nCellsAxis && g.nCells == d.nCells, what + "the cell counts");
    bool grid = std::abs(g.cellSize - d.cellSize) <= 1e-12 * d.cellSize;
    for (std::size_t a = 0; a < 3; ++a)
      grid = grid && std::abs(g.lower[a] - d.lower[a]) <= 1e-12 * (std::abs(d.lower[a]) + d.cellSize);
    check(grid, what + "cell side and lower faces within rounding");
    check(g.voidIndex == d.voidIndex && g.overdenseIndex == d.overdenseIndex &&
          ff.voids.nBasins == r.voids.nBasins && ff.overdense.nBasins == r.overdense.nBasins,
          what + "the partition");
    check(ff.voids.edgeFlag == r.voids.edgeFlag && ff.overdense.edgeFlag == r.overdense.edgeFlag,
          what + "edgeFlag");
    bool close = true;
    double largest = 0.;
    const BasinCatalog* found[2] = {&ff.voids, &ff.overdense};
    const BasinCatalog* searched[2] = {&r.voids, &r.overdense};
    for (std::size_t kind = 0; kind < 2; ++kind) {
      const BasinCatalog& a = *found[kind];
      const BasinCatalog& b = *searched[kind];
      for (std::size_t i = 0; i < b.centre.size(); ++i) {
        const double e = std::abs(a.centre[i] - b.centre[i]);
        largest = std::max(largest, e);
        close = close && e <= 1e-12 * (std::abs(b.centre[i]) + d.cellSize);
      }
      for (std::size_t i = 0; i < b.nBasins; ++i)
        close = close && std::abs(a.radius[i] - b.radius[i]) <= 1e-12 * b.radius[i] &&
                std::abs(a.radiusExt[i] - b.radiusExt[i]) <= 1e-12 * b.radiusExt[i] &&
                std::abs(a.contamination[i] - b.contamination[i]) <= 1e-12;
    }
    std::cout << "    " << format << ": largest centre difference " << largest << std::endl;
    check(close, what + "centres, radii and contamination within 1e-12 relative");
    check(ff.voids.deltaAv.empty() && ff.voids.delta2Mps.empty() && ff.voids.nTracers.empty() &&
          ff.overdense.deltaAv.empty() && ff.overdense.nTracers.empty(), what + "no tracer columns");
    check(ff.voids.fieldCentre.size() == ff.voids.nBasins && ff.overdense.fieldCentre.size() == ff.overdense.nBasins,
          what + "fieldCentre, one per basin");
  }

  group("the lattice: counts, faces at the extreme centres -/+ s / 2, rows in lattice order, cell side the "
        "spacing of the axis with the most nodes, x before y before z on a tie");
  {
    const DivergenceGrid& g = *lr.divergence;
    check(g.nCellsAxis == lf.n && g.nCells == lf.values.size(), "(12, 17, 9) cells");
    bool faces = true;
    for (std::size_t a = 0; a < 3; ++a) {
      const double max = lf.origin[a] + lf.s * (double)(lf.n[a] - 1);
      faces = faces && g.lower[a] == lf.origin[a] - lf.s / 2 && g.upper[a] == max + lf.s / 2;
    }
    check(faces, "lower = min - s / 2, upper = max + s / 2");
    check(same_bits(g.divergence, lf.values), "rows in lattice order: the input's when it is in order");

    for (const std::array<std::size_t, 3> n : {std::array<std::size_t, 3>{6, 9, 4}, {9, 6, 4}, {4, 6, 9},
                                               {7, 7, 3}, {3, 7, 7}, {7, 3, 7}}) {
      const double spacing[3] = {1.5 * (1. + 3e-7), 1.5, 1.5 * (1. - 4e-7)};
      std::vector<double> centres, values;
      for (std::size_t i = 0; i < n[0]; ++i)
        for (std::size_t j = 0; j < n[1]; ++j)
          for (std::size_t k = 0; k < n[2]; ++k) {
            centres.insert(centres.end(), {spacing[0] * (double)i, 2. + spacing[1] * (double)j,
                                           -1. + spacing[2] * (double)k});
            values.push_back(i + j + k == 0 ? -1. : 1.);
          }
      std::size_t axis = 0;
      for (std::size_t a = 1; a < 3; ++a)
        if (n[a] > n[axis]) axis = a;
      const double base[3] = {0., 2., -1.};
      const double expected = (base[axis] + spacing[axis] * (double)(n[axis] - 1) - base[axis]) /
                              (double)(n[axis] - 1);
      const Result f = findVoidsBoxFromField(centres, values, coordinate(0.), {}, 3.);
      const DivergenceGrid& g = *f.divergence;
      check(g.cellSize == expected && g.cellSizeMps == expected / 3. && *f.config.divergenceCellSize == expected / 3.,
            "(" + std::to_string(n[0]) + ", " + std::to_string(n[1]) + ", " + std::to_string(n[2]) +
            "): the spacing of " + "xyz"[axis]);
    }
  }

  group("fieldCentre: the field the watershed split, interpolated at each centre; empty in the finder's "
        "catalogues");
  {
    const DivergenceGrid& g = *lr.divergence;
    bool ok = true;
    for (const BasinCatalog* c : {&lr.voids, &lr.overdense})
      for (std::size_t b = 0; b < c->nBasins; ++b)
        ok = ok && c->fieldCentre[b] == internal::box_field_at(g, g.divergenceSmoothed, &c->centre[3*b]);
    check(ok, "every basin, smoothed field");
    check(r.voids.fieldCentre.empty() && r.overdense.fieldCentre.empty(), "findVoidsBox: empty");
    const Result rd = findVoidsBoxFromDivergence(d, in.tracers, in.config);
    check(rd.voids.fieldCentre.empty() && rd.overdense.fieldCentre.empty(), "findVoidsBoxFromDivergence: empty");
  }

  group("fieldCentre interpolation: a cell centre gives its value; a linear field is reproduced up to "
        "rounding; a point beyond the extreme centres takes the clamped value; an axis of one cell");
  {
    DivergenceGrid g = unit_grid({4, 5, 3});
    std::vector<double> field(g.nCells), linear(g.nCells);
    for (std::size_t i = 0; i < 4; ++i)
      for (std::size_t j = 0; j < 5; ++j)
        for (std::size_t k = 0; k < 3; ++k) {
          const std::size_t c = (i * 5 + j) * 3 + k;
          field[c] = (double)(rng() % 1000) / 7. - 50.;
          linear[c] = 0.3 + 1.25 * (i + 0.5) - 2.5 * (j + 0.5) + 0.75 * (k + 0.5);
        }
    bool centres = true;
    for (std::size_t i = 0; i < 4; ++i)
      for (std::size_t j = 0; j < 5; ++j)
        for (std::size_t k = 0; k < 3; ++k) {
          const double x[3] = {i + 0.5, j + 0.5, k + 0.5};
          centres = centres && internal::box_field_at(g, field, x) == field[(i * 5 + j) * 3 + k];
        }
    check(centres, "at every cell centre, the cell's value exactly");

    bool exact = true;
    for (int t = 0; t < 200; ++t) {
      const double x[3] = {0.5 + 3. * (rng() / 4294967296.), 0.5 + 4. * (rng() / 4294967296.),
                           0.5 + 2. * (rng() / 4294967296.)};
      const double want = 0.3 + 1.25 * x[0] - 2.5 * x[1] + 0.75 * x[2];
      exact = exact && std::abs(internal::box_field_at(g, linear, x) - want) <= 1e-13;
    }
    check(exact, "a linear field, within 1e-13 at 200 points");

    bool clamped = true;
    for (int t = 0; t < 50; ++t) {
      const double x[3] = {-3. + 10. * (rng() / 4294967296.), -3. + 11. * (rng() / 4294967296.),
                           -3. + 9. * (rng() / 4294967296.)};
      const double y[3] = {std::min(3.5, std::max(0.5, x[0])), std::min(4.5, std::max(0.5, x[1])),
                           std::min(2.5, std::max(0.5, x[2]))};
      clamped = clamped && internal::box_field_at(g, field, x) == internal::box_field_at(g, field, y);
    }
    const double corner[3] = {-10., 100., 2.5};
    check(clamped && internal::box_field_at(g, field, corner) == field[(0 * 5 + 4) * 3 + 2],
          "beyond the extreme centres: the value at the clamped point");

    DivergenceGrid flat = unit_grid({4, 1, 3});
    std::vector<double> plane(flat.nCells);
    for (std::size_t c = 0; c < flat.nCells; ++c) plane[c] = (double)(rng() % 100);
    bool one = true;
    for (const double y : {-5., 0.2, 0.5, 0.9, 7.}) {
      const double x[3] = {1.25, y, 2.};
      const double ref[3] = {1.25, 0.5, 2.};
      one = one && internal::box_field_at(flat, plane, x) == internal::box_field_at(flat, plane, ref);
    }
    const double at[3] = {2.5, 3., 0.5};
    check(one && internal::box_field_at(flat, plane, at) == plane[2 * 3 + 0],
          "an axis of one cell: the coordinate on it is not read");
  }

  group("a field result: no tracer, no realization, no assignment; records with and without mps");
  {
    check(lr.nTracers == 0 && lr.nRealizations == 0 && !lr.assignment && lr.ids.empty() &&
          lr.counts.tracers == 0 && lr.counts.displacements == 0, "tracer fields empty, counts 0");
    const DivergenceGrid& g = *lr.divergence;
    check(g.nTracers.empty() && g.nDisplacements.empty() && g.nRealizations == 0 && g.nRealizationsMean == 0,
          "the grid: no nTracers, nDisplacements, realizations");
    check(std::isnan(lr.mps) && std::isnan(g.mps) && std::isnan(g.cellSizeMps) && !lr.config.divergenceCellSize,
          "without mps: mps and cell size in mps NaN, divergenceCellSize empty");
    check(g.volume.size() == g.nCells && g.volume[0] == g.cellSize * g.cellSize * g.cellSize &&
          g.centre.size() == 3 * g.nCells && same_bits(g.centre, lf.centres), "volumes and centres");
    Config c = coordinate(1.5);
    c.smoothing.unit = LengthUnit::Mps;
    c.smoothing.sigma = 0.5;
    const Result m = findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights, 3.);
    check(same_result(m, lr) && m.mps == 3. && m.divergence->mps == 3. &&
          m.divergence->cellSizeMps == lf.s / 3. && *m.config.divergenceCellSize == lf.s / 3.,
          "with mps 3: sigma 0.5 mps = 1.5, the same result; mps and cell size in mps recorded");
    Config k = coordinate(1.5);
    k.keepDivergence = false;
    const Result nk = findVoidsBoxFromField(lf.centres, lf.values, k, lf.weights);
    check(!nk.divergence && same_catalog(nk.voids, lr.voids), "without keepDivergence: no grid, same catalogues");
    check(same_bits(g.divergence, lf.values) && g.divergenceSmoothed.size() == g.nCells &&
          g.weightSmoothed.size() == g.nCells, "sigma > 0: the field as given, its smoothed field and weight");
    const Result zero = findVoidsBoxFromField(lf.centres, lf.values, coordinate(0.), lf.weights);
    const DivergenceGrid& z = *zero.divergence;
    check(same_bits(z.divergence, lf.values) && z.divergenceSmoothed.empty() && z.weightSmoothed.empty(),
          "sigma 0: the field as given, no smoothed field or weight");
  }

  group("without mps: a sigma in mps refused only when positive");
  {
    Config c = coordinate(0.);
    c.smoothing.unit = LengthUnit::Mps;
    const Result z = findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights);
    check(same_result(z, findVoidsBoxFromField(lf.centres, lf.values, coordinate(0.), lf.weights)),
          "sigma 0 in mps: runs, as sigma 0 in coordinate units");
    c.smoothing.sigma = 1.;
    throws_naming([&] { findVoidsBoxFromField(lf.centres, lf.values, c, lf.weights); },
                  "the smoothing sigma is in mps, but no mps is given", "sigma 1 in mps: refused");
  }

  group("weights: used as given, a scale by a power of 2 changes neither the search nor the threshold; none is "
        "weight 1 everywhere; all 0 accepted at sigma 0");
  {
    std::vector<double> scaled = lf.weights;
    for (double& w : scaled) w *= 1024.;
    const Result sr = findVoidsBoxFromField(lf.centres, lf.values, coordinate(1.5), scaled);
    bool times = sr.divergence->weightSmoothed.size() == lr.divergence->weightSmoothed.size();
    for (std::size_t c = 0; times && c < lr.divergence->nCells; ++c)
      times = sr.divergence->weight[c] == 1024. * lr.divergence->weight[c] &&
              sr.divergence->weightSmoothed[c] == 1024. * lr.divergence->weightSmoothed[c];
    check(same_search(sr, lr) && times, "weights times 1024: the same search, weight and weightSmoothed times 1024");
    const Result none = findVoidsBoxFromField(lf.centres, lf.values, coordinate(1.5));
    const std::vector<double> fours(lf.values.size(), 4.);
    const Result fr = findVoidsBoxFromField(lf.centres, lf.values, coordinate(1.5), fours);
    check(same_search(fr, none) && none.divergence->weight == std::vector<double>(lf.values.size(), 1.) &&
          fr.divergence->weight == fours,
          "no weights: the grid's weight 1; weights all 4: kept, the same search");
    const std::vector<double> zeros(lf.values.size(), 0.);
    check(findVoidsBoxFromField(lf.centres, lf.values, coordinate(0.), zeros).voids.nBasins > 0,
          "all 0 at sigma 0: runs");
  }

  group("the order of the rows: shuffled rows give the same result, bit for bit");
  {
    std::vector<std::size_t> perm(lf.values.size());
    for (std::size_t i = 0; i < perm.size(); ++i) perm[i] = i;
    for (std::size_t i = perm.size() - 1; i > 0; --i) std::swap(perm[i], perm[rng() % (i + 1)]);
    std::vector<double> centres = lf.centres, values = lf.values, weights = lf.weights;
    permute(perm, centres, values, weights);
    check(same_result(findVoidsBoxFromField(centres, values, coordinate(1.5), weights), lr), "LatticeField");

    const io::Field f = {d.centre, d.divergence, d.weight};
    std::vector<std::size_t> rows(d.nCells);
    for (std::size_t i = 0; i < rows.size(); ++i) rows[i] = i;
    for (std::size_t i = rows.size() - 1; i > 0; --i) std::swap(rows[i], rows[rng() % (i + 1)]);
    std::vector<double> dc = f.centres, dv = f.values, dw = f.weights;
    permute(rows, dc, dv, dw);
    Config c = in.config;
    c.keepAssignment = false;
    c.nRealizationsMean = 1;
    check(same_result(findVoidsBoxFromField(dc, dv, c, dw, Translated::kMps),
                      findVoidsBoxFromField(f.centres, f.values, c, f.weights, Translated::kMps)),
          "the Translated grid");
  }

  group("errors of a search on a field");
  {
    const std::vector<double>& x = lf.centres;
    const std::vector<double>& v = lf.values;
    const Config c = coordinate(1.5);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    throws_naming([&] { findVoidsBoxFromField({}, {}, c); }, "findVoidsBoxFromField: centres is empty", "empty");
    throws_naming([&] { findVoidsBoxFromField({1., 2.}, {1.}, c); }, "centres", "size");
    throws_naming([&] { auto y = x; y[5] = nan; findVoidsBoxFromField(y, v, c); }, "centres", "non-finite");
    throws_naming([&] { auto w = v; w.pop_back(); findVoidsBoxFromField(x, w, c); },
                  "findVoidsBoxFromField: values holds " + std::to_string(v.size() - 1) + " entries for", "size");
    throws_naming([&] { auto w = v; w[3] = nan; findVoidsBoxFromField(x, w, c); },
                  "the field values hold a non-finite value at cell 3 (from 0)", "non-finite");
    throws_naming([&] { findVoidsBoxFromField(x, v, c, {1., 2.}); }, "or none", "size");
    throws_naming([&] { auto w = lf.weights; w[4] = -1.; findVoidsBoxFromField(x, v, c, w); },
                  "at cell 4 (from 0); every weight must be finite and non-negative", "negative");
    throws_naming([&] { findVoidsBoxFromField(x, v, c, {}, 0.); }, "the mean particle separation is ", "the mean particle separation is ");
    throws_naming([&] { Config y = c; y.divergenceCellSize = 1.; findVoidsBoxFromField(x, v, y); },
                  "the divergence cell size is set for a search on a field", "cell size");
    throws_naming([&] { Config y = c; y.keepAssignment = true; findVoidsBoxFromField(x, v, y); },
                  "a search on a field has no tracers to assign to cells; do not ask for the assignment", "keepAssignment");
    throws_naming([&] { Config y = c; y.nRealizationsMean = 2; findVoidsBoxFromField(x, v, y); },
                  "for a search on a field, which has no realizations", "nRealizationsMean");
    throws_naming([&] { Config y = c; y.smoothing.unit = LengthUnit::Mps; findVoidsBoxFromField(x, v, y); },
                  "no mps is given", "no mps is given");
    throws_naming([&] { findVoidsBoxFromField(x, v, coordinate(1e-200)); }, "too small to build the kernel", "too small to build the kernel");
    throws_naming([&] { findVoidsBoxFromField(x, v, c, std::vector<double>(v.size(), 0.)); },
                  "every weight is 0", "every weight is 0");
    throws_naming([&] { findVoidsBoxFromField({-1e308, 0., 0., 0., 0., 0., 1e308, 0., 0.}, {-1., 1., 1.}, c); },
                  "too wide a range for the cell side to be finite", "too wide a range for the cell side to be finite");
    throws_naming([&] { findVoidsBoxFromField({1., 2., 3.}, {-1.}, c); },
                  "the centres span one node on every axis", "the centres span one node on every axis");
    throws_naming([&] { auto y = x; for (std::size_t i = 0; i < v.size(); ++i) y[3*i+2] *= 1.01;
                        findVoidsBoxFromField(y, v, c); },
                  "the cells are not cubic", "the cells are not cubic");
    throws_naming([&] { auto y = x; auto w = v; y.erase(y.begin() + 30, y.begin() + 33); w.erase(w.begin() + 10);
                        findVoidsBoxFromField(y, w, c); },
                  "the centres form a lattice of (12, 17, 9) nodes but are only", "hole");
    throws_naming([&] { auto y = x; y[3*7] += 0.3 * lf.s; findVoidsBoxFromField(y, v, c); },
                  "the centre of cell 7 (from 0) lies ", "off node");
    throws_naming([&] { auto y = x; for (std::size_t a = 0; a < 3; ++a) y[3*7+a] = y[3*8+a];
                        findVoidsBoxFromField(y, v, c); },
                  "the centres of cells 7 and 8 (from 0) lie at the same node", "duplicate");
    throws_naming([&] { findVoidsBoxFromField(x, std::vector<double>(v.size(), 0.5), c); },
                  "no void found: the field is 0 or positive in every cell", "no void found: the field is 0 or positive in every cell");
  }

  group("non-cubic lattices with an axis of one cell or two, the field lowered by 2");
  for (const std::array<std::size_t, 3> n : {std::array<std::size_t, 3>{1, 17, 9}, {12, 1, 9}, {12, 17, 1},
                                             {2, 17, 9}}) {
    const LatticeField f(n);
    std::vector<double> lowered = f.values;
    for (double& value : lowered) value -= 2.;
    const Result one = findVoidsBoxFromField(f.centres, lowered, coordinate(1.5), f.weights);
    check(one.divergence->nCellsAxis == n && one.voids.nBasins > 0,
          "(" + std::to_string(n[0]) + ", " + std::to_string(n[1]) + ", " + std::to_string(n[2]) + ")");
  }

  group("the same result with 1 and 4 threads, and the recorded hash");
  {
    set_threads(1);
    const Result one = findVoidsBoxFromField(lf.centres, lf.values, coordinate(1.5), lf.weights);
    set_threads(4);
    const Result four = findVoidsBoxFromField(lf.centres, lf.values, coordinate(1.5), lf.weights);
    check(same_result(one, four) && same_result(one, lr), "1 and 4 threads");
    const std::uint64_t h = hash(one);
    if (h != kRecordedHash) std::cerr << "field search hash: 0x" << std::hex << h << std::dec << std::endl;
    check(h == kRecordedHash, "the hash equals the recorded one");
  }

  return report("test_from_field");
}
