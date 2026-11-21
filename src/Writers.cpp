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
 *  @file src/Writers.cpp
 *
 *  @brief The writers of bitvf::io.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "otswap/io.h"

#include "bitvf/io.h"
#include "internal.h"

namespace {

  using bitvf::Error;
  using otswap::io::Column;
  using otswap::io::Keyword;

  constexpr int kPrecision = 9;
  constexpr int kLossless = 17;

  struct Source {
    Column column;
    std::function<double(std::size_t)> value;
    std::function<std::int64_t(std::size_t)> integer;
  };

  Source real (const std::string& name, const std::string& description, std::function<double(std::size_t)> value)
  {
    return {Column{name, 'D', description, "", {}, {}}, std::move(value), {}};
  }

  Source flag (const std::string& name, const std::string& description, const std::vector<std::uint8_t>& values)
  {
    return {Column{name, 'J', description, "", {}, {}}, [&values] (std::size_t i) { return (double)values[i]; }, {}};
  }

  Source integer (const std::string& name, const std::string& description,
                  std::function<std::int64_t(std::size_t)> value)
  {
    return {Column{name, 'K', description, "", {}, {}}, {}, std::move(value)};
  }

  std::string number (const double value)
  {
    for (int digits = 15; ; ++digits) {
      std::ostringstream out;
      out.imbue(std::locale::classic());
      out << std::setprecision(digits) << value;
      if (digits == kLossless || std::strtod(out.str().c_str(), nullptr) == value) return out.str();
    }
  }

  std::string exact (const double value)
  {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(kLossless) << value;
    return out.str();
  }

  void check_size (const std::size_t size, const std::size_t expected, const std::string& field,
                   const std::string& writer)
  {
    if (size != expected)
      throw Error(writer + ": cannot write " + field + ": it holds " + std::to_string(size) + " entries where " +
                  std::to_string(expected) + " are expected");
  }

  template <typename T>
  void check_int64 (const std::vector<T>& values, const std::string& field, const std::string& index,
                    const std::string& writer)
  {
    for (std::size_t i = 0; i < values.size(); ++i)
      if (values[i] > (T)std::numeric_limits<std::int64_t>::max())
        throw Error(writer + ": cannot write " + field + ": " + index + " " + std::to_string(i) + " (from 0) is " +
                    std::to_string(values[i]) + ", above the largest 64-bit integer");
  }

  /// Whether a result is a field result (findVoidsBoxFromField).
  bool is_field (const bitvf::Result& result)
  {
    return result.geometry == otswap::Geometry::Box && result.nTracers == 0;
  }

  void check_result (const bitvf::Result& result, const std::string& writer)
  {
    if (!result.config.divergenceCellSize && !(is_field(result) && std::isnan(result.mps)))
      throw Error(writer + ": the result is malformed: it records no divergence cell size");
    if (result.config.nRealizationsMean == 0)
      throw Error(writer + ": the result is malformed: it records 0 realizations averaged");
  }

  std::vector<Keyword> keywords (const std::string& product, const bitvf::Result& result,
                                 const char* countName, const std::size_t count, const char* countComment)
  {
    const bitvf::Config& config = result.config;
    const bool inMps = config.smoothing.unit == bitvf::LengthUnit::Mps;
    const bool field = is_field(result);
    const bool lightcone = result.geometry == otswap::Geometry::Lightcone;
    std::vector<Keyword> k = {
      {"PRODUCT", product, "what the file holds"},
      {"BITVFV", bitvf::internal::version, "bitvf version"},
      {"GEOMETRY", lightcone ? "lightcone" : "box", "geometry of the search"},
      {countName, std::to_string(count), countComment}};
    if (!field) {
      k.push_back({"NREC", std::to_string(result.nRealizations), "realizations"});
      k.push_back({"NRECMEAN", std::to_string(config.nRealizationsMean),
                   "realizations averaged into each displacement"});
    }
    k.push_back({"SIGMA", number(config.smoothing.sigma), "smoothing sigma, in SIGUNIT"});
    k.push_back({"SIGUNIT", inMps ? "mps" : "coordinate", "unit of SIGMA"});
    if (!field || !std::isnan(result.mps)) {
      k.push_back({"CELLSIZE", number(*config.divergenceCellSize), "divergence cell size, in units of the mps"});
      if (!lightcone)
        k.push_back({"MPS", number(result.mps), "mean particle separation, in coordinate units"});
    }
    if (result.randomSeed)
      k.push_back({"RANDSEED", std::to_string(*result.randomSeed), "seed of the randoms bitvf drew"});
    if (lightcone) {
      if (std::isfinite(result.redshiftCut.min))
        k.push_back({"ZCUTMIN", number(result.redshiftCut.min), "lower bound of the redshift cut"});
      if (std::isfinite(result.redshiftCut.max))
        k.push_back({"ZCUTMAX", number(result.redshiftCut.max), "upper bound of the redshift cut"});
      k.push_back({"NRANDUSE", std::to_string(result.counts.randoms), "randoms used"});
      k.push_back({"MASKNSID", std::to_string(result.maskNside), "NSIDE of the mask"});
    }
    if (result.limits) {
      k.push_back({"BOXLIM", "1", "the search had box limits"});
      static const char* const axis[3] = {"X", "Y", "Z"};
      for (int c = 0; c < 3; ++c) {
        k.push_back({std::string("BOX") + axis[c] + "LO", number(result.limits->lower[c]),
                     std::string("box limit, lower ") + (char)('x' + c) + ", coordinate units"});
        k.push_back({std::string("BOX") + axis[c] + "HI", number(result.limits->upper[c]),
                     std::string("box limit, upper ") + (char)('x' + c) + ", coordinate units"});
      }
    }
    return k;
  }

  void write_table (const std::string& file, const std::vector<Source>& sources, const std::size_t nRows,
                    const int precision, std::vector<Keyword> keywords)
  {
    std::vector<Column> columns;
    columns.reserve(sources.size());
    for (const Source& s : sources) columns.push_back(s.column);

    otswap::io::WriteOptions options;
    options.precision = precision;
    options.keywords = std::move(keywords);

    const auto fill = [&sources] (const std::size_t row, std::vector<double>& values,
                                  std::vector<std::int64_t>& integers) {
      for (std::size_t c = 0; c < sources.size(); ++c) {
        if (sources[c].column.type == 'K') integers[c] = sources[c].integer(row);
        else values[c] = sources[c].value(row);
      }
    };

    try {
      otswap::io::write(file, columns, nRows, fill, options);
    }
    catch (const otswap::Error& e) {
      std::throw_with_nested(Error(std::string("otswap: ") + e.what()));
    }
  }

  void write_catalog (const std::string& file, const bitvf::Result& result, const bitvf::BasinCatalog& cat,
                      const std::string& name, const std::string& product,
                      const std::string& writer)
  {
    check_result(result, writer);
    const bool field = is_field(result);
    const std::size_t n = cat.nBasins;
    const std::size_t nt = field ? 0 : n;
    check_size(cat.centre.size(), 3 * n, name + ".centre", writer);
    check_size(cat.centreSky.size(), result.geometry == otswap::Geometry::Lightcone ? 3 * n : 0,
               name + ".centreSky", writer);
    check_size(cat.radius.size(), n, name + ".radius", writer);
    check_size(cat.radiusExt.size(), n, name + ".radiusExt", writer);
    check_size(cat.deltaAv.size(), nt, name + ".deltaAv", writer);
    check_size(cat.delta2Mps.size(), nt, name + ".delta2Mps", writer);
    check_size(cat.edgeFlag.size(), n, name + ".edgeFlag", writer);
    check_size(cat.contamination.size(), n, name + ".contamination", writer);
    check_size(cat.nTracers.size(), nt, name + ".nTracers", writer);
    check_size(cat.fieldCentre.size(), field ? n : 0, name + ".fieldCentre", writer);
    check_int64(cat.nTracers, name + ".nTracers", "basin", writer);

    std::vector<Source> sources = {
      integer("index", "row of the basin", [] (std::size_t b) { return (std::int64_t)b; })};
    if (result.geometry == otswap::Geometry::Lightcone) {
      constexpr double kDegrees = 57.295779513082320876798154814105;
      sources.push_back(real("centreRA", "right ascension of the centre, in degrees",
                             [&cat] (std::size_t b) { return cat.centreSky[3*b] * kDegrees; }));
      sources.push_back(real("centreDec", "declination of the centre, in degrees",
                             [&cat] (std::size_t b) { return cat.centreSky[3*b+1] * kDegrees; }));
      sources.push_back(real("centreRed", "redshift of the centre",
                             [&cat] (std::size_t b) { return cat.centreSky[3*b+2]; }));
    }
    sources.insert(sources.end(), {
      real("centreX", "x of the centre, in coordinate units", [&cat] (std::size_t b) { return cat.centre[3*b]; }),
      real("centreY", "y of the centre, in coordinate units", [&cat] (std::size_t b) { return cat.centre[3*b+1]; }),
      real("centreZ", "z of the centre, in coordinate units", [&cat] (std::size_t b) { return cat.centre[3*b+2]; }),
      real("radius", "equal-volume radius of the core, in coordinate units",
           [&cat] (std::size_t b) { return cat.radius[b]; }),
      real("radiusExt", "equal-volume radius of the extended basin, in coordinate units",
           [&cat] (std::size_t b) { return cat.radiusExt[b]; })});
    if (!field) {
      sources.push_back(real("deltaAv", "mean density contrast of the tracers in the core",
                             [&cat] (std::size_t b) { return cat.deltaAv[b]; }));
      sources.push_back(real("delta2Mps", "density contrast within 2 mps of the centre",
                             [&cat] (std::size_t b) { return cat.delta2Mps[b]; }));
    }
    sources.push_back(flag("edgeFlag", "contact with the boundary: 0 none, 1 extended basin only, 2 core",
                           cat.edgeFlag));
    sources.push_back(real("contamination", "fraction of the surface on the boundary",
                           [&cat] (std::size_t b) { return cat.contamination[b]; }));
    if (field)
      sources.push_back(real("fieldCentre", "field at the centre, interpolated between the cell centres",
                             [&cat] (std::size_t b) { return cat.fieldCentre[b]; }));
    else
      sources.push_back(integer("nTracers", "tracers starting in the core",
                                [&cat] (std::size_t b) { return (std::int64_t)cat.nTracers[b]; }));

    write_table(file, sources, n, kPrecision, keywords(product, result, "NBASINS", n, "number of basins"));
  }

}


// ============================================================================


std::string bitvf::io::fileName (const Product product, const std::string& directory, const std::string& name,
                                 const bool zspace)
{
  if (name.empty()) throw Error("the file name is empty");

  const char* stem = product == Product::Voids ? "voids_"
                   : product == Product::Overdense ? "overdense_basins_"
                   : product == Product::Divergence ? "divergence_" : "galaxies_";

  std::string path = directory;
  if (!path.empty() && path.back() != '/') path += '/';
  return path + stem + (zspace ? "zspace_" : "") + name;
}


// ============================================================================


void bitvf::io::writeVoids (const std::string& file, const Result& result)
{
  write_catalog(file, result, result.voids, "voids", "voids", "writeVoids");
}


// ============================================================================


void bitvf::io::writeOverdense (const std::string& file, const Result& result)
{
  write_catalog(file, result, result.overdense, "overdense", "overdense_basins", "writeOverdense");
}


// ============================================================================


void bitvf::io::writeDivergence (const std::string& file, const Result& result)
{
  if (result.geometry != otswap::Geometry::Box)
    throw Error("writeDivergence: writing a lightcone result is not implemented");
  if (!result.divergence)
    throw Error("writeDivergence: the result holds no divergence grid; run the search with the divergence grid kept");
  check_result(result, "writeDivergence");

  const DivergenceGrid& grid = *result.divergence;
  const bool field = is_field(result);
  const std::string cannot = "writeDivergence: cannot write the divergence: ";
  if (field) {
    if (!std::isfinite(grid.cellSize) || grid.cellSize <= 0.)
      throw Error(cannot + "the grid has a cell size of " + std::to_string(grid.cellSize) +
                  " in coordinate units; it must be finite and positive");
    for (std::size_t k = 0; k < 3; ++k)
      if (!std::isfinite(grid.lower[k]) || !std::isfinite(grid.upper[k]))
        throw Error(cannot + "the grid has a bound that is not finite on axis " + "xyz"[k]);
    std::size_t cells = 1;
    bool fits = true;
    for (const std::size_t count : grid.nCellsAxis) {
      if (count == 0 || cells > std::numeric_limits<std::size_t>::max() / count) fits = false;
      else cells *= count;
    }
    if (!fits || cells != grid.nCells)
      throw Error(cannot + "the grid has " + std::to_string(grid.nCells) + " cells for (" +
                  std::to_string(grid.nCellsAxis[0]) + ", " + std::to_string(grid.nCellsAxis[1]) + ", " +
                  std::to_string(grid.nCellsAxis[2]) + ") along the axes; it must have their product");
  }
  else
    try {
      internal::check_box_grid(grid, result.mps, false, "");
    }
    catch (const Error& e) {
      throw Error(cannot + e.what());
    }

  const std::size_t n = grid.nCells;
  const std::size_t nt = field ? 0 : n;
  check_size(grid.centre.size(), 3 * n, "divergence.centre", "writeDivergence");
  check_size(grid.divergence.size(), n, "divergence.divergence", "writeDivergence");
  if (!grid.divergenceSmoothed.empty())
    check_size(grid.divergenceSmoothed.size(), n, "divergence.divergenceSmoothed", "writeDivergence");
  check_size(grid.volume.size(), n, "divergence.volume", "writeDivergence");
  check_size(grid.weight.size(), n, "divergence.weight", "writeDivergence");
  if (!grid.weightSmoothed.empty())
    check_size(grid.weightSmoothed.size(), n, "divergence.weightSmoothed", "writeDivergence");
  check_size(grid.nDisplacements.size(), nt, "divergence.nDisplacements", "writeDivergence");
  check_size(grid.nTracers.size(), nt, "divergence.nTracers", "writeDivergence");
  check_size(grid.voidIndex.size(), n, "divergence.voidIndex", "writeDivergence");
  check_size(grid.overdenseIndex.size(), n, "divergence.overdenseIndex", "writeDivergence");
  check_int64(grid.nDisplacements, "divergence.nDisplacements", "cell", "writeDivergence");
  check_int64(grid.nTracers, "divergence.nTracers", "cell", "writeDivergence");

  const std::string malformed = "writeDivergence: the result is malformed: ";
  if (!field) {
    if (grid.mps != result.mps)
      throw Error(malformed + "its divergence grid has an mps of " + exact(grid.mps) +
                  " in coordinate units, but the search's is " + exact(result.mps) + "; they must be equal");
    if (grid.cellSizeMps != *result.config.divergenceCellSize)
      throw Error(malformed + "its divergence grid has a cell size of " + exact(grid.cellSizeMps) +
                  " mps, but its configuration records " + exact(*result.config.divergenceCellSize) +
                  "; they must be equal");
    if (grid.limits.has_value() != result.limits.has_value() ||
        (grid.limits && (grid.limits->lower != result.limits->lower || grid.limits->upper != result.limits->upper)))
      throw Error(malformed + "the box limits of its divergence grid differ from those of the search; "
                  "they must be equal");
    if (grid.nRealizations != result.nRealizations)
      throw Error(malformed + "its divergence grid records " + std::to_string(grid.nRealizations) +
                  " realizations, but the search's are " + std::to_string(result.nRealizations) +
                  "; they must be equal");
    if (grid.nRealizationsMean != result.config.nRealizationsMean)
      throw Error(malformed + "its divergence grid records " + std::to_string(grid.nRealizationsMean) +
                  " realizations averaged, but its configuration records " +
                  std::to_string(result.config.nRealizationsMean) + "; they must be equal");
  }

  std::vector<Source> sources = {
    integer("index", "row of the cell", [] (std::size_t c) { return (std::int64_t)c; }),
    real("cellX", "x of the cell centre, in coordinate units", [&grid] (std::size_t c) { return grid.centre[3*c]; }),
    real("cellY", "y of the cell centre, in coordinate units",
         [&grid] (std::size_t c) { return grid.centre[3*c+1]; }),
    real("cellZ", "z of the cell centre, in coordinate units",
         [&grid] (std::size_t c) { return grid.centre[3*c+2]; }),
    real("divergence", field ? "field, before smoothing" : "divergence of the displacement field, before smoothing",
         [&grid] (std::size_t c) { return grid.divergence[c]; })};
  if (!grid.divergenceSmoothed.empty())
    sources.push_back(real("divergenceSmoothed",
                           field ? "field, after smoothing" : "divergence of the displacement field, after smoothing",
                           [&grid] (std::size_t c) { return grid.divergenceSmoothed[c]; }));
  sources.push_back(real("volume", "volume of the cell, in coordinate units cubed",
                         [&grid] (std::size_t c) { return grid.volume[c]; }));
  sources.push_back(real("weight", field ? "weight of the cell in the smoothing" : "tracer weight of the cell",
                         [&grid] (std::size_t c) { return grid.weight[c]; }));
  if (!grid.weightSmoothed.empty())
    sources.push_back(real("weightSmoothed",
                           field ? "weight of the cell, after smoothing" : "tracer weight of the cell, after smoothing",
                           [&grid] (std::size_t c) { return grid.weightSmoothed[c]; }));
  if (!field) {
    sources.push_back(integer("nDispl", "displacement ends counted in the cell",
                              [&grid] (std::size_t c) { return (std::int64_t)grid.nDisplacements[c]; }));
    sources.push_back(integer("nTracers", "tracers starting in the cell",
                              [&grid] (std::size_t c) { return (std::int64_t)grid.nTracers[c]; }));
  }
  sources.push_back(integer("voidIndex", "row of the void of the cell, -1 none",
                            [&grid] (std::size_t c) { return grid.voidIndex[c]; }));
  sources.push_back(integer("overdenseIndex", "row of the overdense basin of the cell, -1 none",
                            [&grid] (std::size_t c) { return grid.overdenseIndex[c]; }));

  std::vector<Keyword> k = keywords("divergence", result, "NCELLS", n, "number of cells");
  static const char* const axis[3] = {"X", "Y", "Z"};
  for (int c = 0; c < 3; ++c)
    k.push_back({std::string("LOWER") + axis[c], number(grid.lower[c]),
                 std::string(field ? "grid" : "padded box") + ", lower " + (char)('x' + c) + ", coordinate units"});
  for (int c = 0; c < 3; ++c)
    k.push_back({std::string("UPPER") + axis[c], number(grid.upper[c]),
                 std::string(field ? "grid" : "padded box") + ", upper " + (char)('x' + c) + ", coordinate units"});
  k.push_back({"CELLSIDE", number(grid.cellSize), "cell side, in coordinate units"});
  for (int c = 0; c < 3; ++c)
    k.push_back({std::string("NCELL") + axis[c], std::to_string(grid.nCellsAxis[c]),
                 std::string("cells along ") + (char)('x' + c)});

  write_table(file, sources, n, kLossless, std::move(k));
}


// ============================================================================


void bitvf::io::writeGalaxies (const std::string& file, const Result& result)
{
  if (is_field(result))
    throw Error("writeGalaxies: the result is a search on a field (findVoidsBoxFromField), which has no tracers: "
                "it has no galaxies file");
  if (!result.assignment)
    throw Error("writeGalaxies: the result holds no tracer assignment; run the search with the tracer assignment kept");
  check_result(result, "writeGalaxies");

  const TracerAssignment& a = *result.assignment;
  const std::size_t n = result.nTracers;
  if (!result.ids.empty()) check_size(result.ids.size(), n, "ids", "writeGalaxies");
  if (!result.uncorrected.empty()) check_size(result.uncorrected.size(), n, "uncorrected", "writeGalaxies");
  if (!result.outsideGrid.empty()) check_size(result.outsideGrid.size(), n, "outsideGrid", "writeGalaxies");
  if (!result.outsideCut.empty()) check_size(result.outsideCut.size(), n, "outsideCut", "writeGalaxies");
  check_size(a.cellIndex.size(), n, "assignment.cellIndex", "writeGalaxies");
  check_size(a.voidIndex.size(), n, "assignment.voidIndex", "writeGalaxies");
  check_size(a.overdenseIndex.size(), n, "assignment.overdenseIndex", "writeGalaxies");
  check_size(a.distance.size(), n, "assignment.distance", "writeGalaxies");
  const bool lightcone = result.geometry == otswap::Geometry::Lightcone;
  const bool boundedCut = std::isfinite(result.redshiftCut.min) || std::isfinite(result.redshiftCut.max);
  check_size(a.inMask.size(), lightcone ? n : 0, "assignment.inMask", "writeGalaxies");
  if (lightcone && (boundedCut || !a.inCut.empty())) check_size(a.inCut.size(), n, "assignment.inCut", "writeGalaxies");
  if (!lightcone) check_size(a.inCut.size(), 0, "assignment.inCut", "writeGalaxies");
  check_size(a.inBox.size(), lightcone ? 0 : n, "assignment.inBox", "writeGalaxies");
  if (!a.used.empty()) check_size(a.used.size(), n, "assignment.used", "writeGalaxies");

  const std::vector<std::int64_t>& ids = result.ids;
  std::vector<Source> sources = {
    integer("id", "ID of the tracer, or its input row",
            [&ids] (std::size_t i) { return ids.empty() ? (std::int64_t)i : ids[i]; }),
    integer("cellIndex", "row of the cell the tracer starts in, -1 none",
            [&a] (std::size_t i) { return a.cellIndex[i]; }),
    integer("voidIndex", "row of the void of the tracer, -1 none", [&a] (std::size_t i) { return a.voidIndex[i]; }),
    integer("overdenseIndex", "row of the overdense basin of the tracer, -1 none",
            [&a] (std::size_t i) { return a.overdenseIndex[i]; }),
    real("distance", "distance from the centre of its basin, in coordinate units",
         [&a] (std::size_t i) { return a.distance[i]; })};
  if (lightcone) sources.push_back(flag("inMask", "0 on an unobserved pixel of the mask", a.inMask));
  if (lightcone && boundedCut) sources.push_back(flag("inCut", "0 outside the redshift cut, skipped", a.inCut));
  if (result.limits) sources.push_back(flag("inBox", "0 outside the box limits, skipped", a.inBox));
  if (!a.used.empty())
    sources.push_back(flag("used", "1 when a displacement of the tracer entered the divergence", a.used));
  if (!result.uncorrected.empty())
    sources.push_back(flag("uncorrected", "1 when left unmoved by the redshift-space correction, excluded",
                           result.uncorrected));
  if (!result.outsideGrid.empty())
    sources.push_back(flag("outsideGrid", "1 when outside the reused grid, excluded", result.outsideGrid));
  if (!result.outsideCut.empty())
    sources.push_back(flag("outsideCut", "1 when moved outside the redshift cut, excluded", result.outsideCut));

  write_table(file, sources, n, kPrecision, keywords("galaxies", result, "NOBJECTS", n, "number of input tracers"));
}
