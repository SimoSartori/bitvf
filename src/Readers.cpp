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
 *  @file src/Readers.cpp
 *
 *  @brief The readers of bitvf::io.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <exception>
#include <fstream>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "otswap/io.h"

#include "bitvf/io.h"
#include "internal.h"

namespace {

  using bitvf::Error;
  using otswap::io::Keyword;

  constexpr std::uint64_t kUnsigned = std::numeric_limits<unsigned>::max();
  constexpr std::uint64_t kSize = std::numeric_limits<std::size_t>::max();

  const char* const kBoxKeywords[6] = {"BOXXLO", "BOXXHI", "BOXYLO", "BOXYHI", "BOXZLO", "BOXZHI"};
  const char* const kAxis[3] = {"X", "Y", "Z"};

  std::string lower_case (std::string s)
  {
    std::transform(s.begin(), s.end(), s.begin(), [] (unsigned char c) { return (char)std::tolower(c); });
    return s;
  }

  /**
   *  @brief The keywords of a table, looked up by name, with errors that
   *  name the file and the keyword.
   */
  struct Header {

    std::string file;
    std::vector<Keyword> keywords;

    const std::string* find (const std::string& name) const
    {
      for (const Keyword& k : keywords)
        if (k.name == name) return &k.value;
      return nullptr;
    }

    std::string text (const std::string& name) const
    {
      const std::string* value = find(name);
      if (value == nullptr) throw Error(file + ": keyword " + name + " is missing");
      return *value;
    }

    double real (const std::string& name) const
    {
      const std::string v = text(name);
      std::istringstream in(v);
      in.imbue(std::locale::classic());
      double x = 0.;
      if (v.empty() || !(in >> x) || in.get() != std::char_traits<char>::eof())
        throw Error(file + ": keyword " + name + " is '" + v + "', not a number");
      return x;
    }

    std::uint64_t whole (const std::string& name, const std::uint64_t max) const
    {
      const std::string v = text(name);
      std::uint64_t n = 0;
      const auto conversion = std::from_chars(v.data(), v.data() + v.size(), n);
      if (conversion.ec != std::errc() || conversion.ptr != v.data() + v.size() || n > max)
        throw Error(file + ": keyword " + name + " is '" + v + "', not a whole number in [0, " +
                    std::to_string(max) + "]");
      return n;
    }

    void product (const std::string& expected, const std::string& what, const std::string& source) const
    {
      const std::string* product = find("PRODUCT");
      if (product == nullptr) throw Error(file + " has no PRODUCT keyword: not " + source + " file");
      const std::string found = *product;
      if (found != expected)
        throw Error(file + " is not a " + what + " file: its PRODUCT is '" + found + "', where '" + expected +
                    "' is expected" +
                    (found == "complete displacements"
                     ? "; read it with otswap::io::readDisplacementComplete and pass the otswap::Result to "
                       "findVoidsBox"
                     : ""));
    }

    void box_geometry (const std::string& what) const
    {
      const std::string geometry = text("GEOMETRY");
      if (geometry == "lightcone")
        throw Error("reading a lightcone " + what + " file is not implemented");
      if (geometry != "box")
        throw Error(file + ": GEOMETRY is '" + geometry + "'; it must be box or lightcone");
    }

    /**
     *  @brief The box limits: the six box keywords when BOXLIM = 1, else
     *  none. A box keyword without BOXLIM is an error.
     */
    std::optional<otswap::BoxLimits> limits () const
    {
      const std::string* boxlim = find("BOXLIM");
      if (boxlim != nullptr && *boxlim != "1")
        throw Error(file + ": keyword BOXLIM is '" + *boxlim + "', where only 1 is written");
      if (boxlim == nullptr) {
        for (const char* k : kBoxKeywords)
          if (find(k) != nullptr)
            throw Error(file + ": keyword " + std::string(k) + " is present without BOXLIM = 1, "
                        "which bitvf writes with the box limits");
        return std::nullopt;
      }
      for (const char* k : kBoxKeywords)
        if (find(k) == nullptr)
          throw Error(file + ": BOXLIM = 1 records box limits, but keyword " + std::string(k) + " is missing");
      otswap::BoxLimits box;
      for (int c = 0; c < 3; ++c) {
        box.lower[c] = real(kBoxKeywords[2*c]);
        box.upper[c] = real(kBoxKeywords[2*c+1]);
      }
      return box;
    }

    /**
     *  @brief Whether @p names, the columns of the table, holds that
     *  name, ignoring case.
     */
    static bool has_column (const std::vector<std::string>& names, const std::string& name)
    {
      const std::string wanted = lower_case(name);
      for (const std::string& n : names)
        if (lower_case(n) == wanted) return true;
      return false;
    }
  };

}


// ============================================================================


bitvf::DivergenceGrid bitvf::io::readDivergence (const std::string& file)
{
  const Header h {file, internal::call_otswap([&file] { return otswap::io::readKeywords(file); })};
  h.product("divergence", "divergence", "a bitvf divergence");
  h.box_geometry("divergence");

  DivergenceGrid grid;
  grid.geometry = otswap::Geometry::Box;
  grid.nCells = (std::size_t)h.whole("NCELLS", kSize);
  grid.nRealizations = (unsigned)h.whole("NREC", kUnsigned);
  if (grid.nRealizations == 0) throw Error(file + ": NREC is 0; it must be at least 1");
  grid.nRealizationsMean = (unsigned)h.whole("NRECMEAN", kUnsigned);
  if (grid.nRealizationsMean == 0) throw Error(file + ": NRECMEAN is 0; it must be at least 1");
  grid.mps = h.real("MPS");
  grid.cellSizeMps = h.real("CELLSIZE");
  grid.cellSize = h.real("CELLSIDE");
  for (int c = 0; c < 3; ++c) {
    grid.lower[c] = h.real(std::string("LOWER") + kAxis[c]);
    grid.upper[c] = h.real(std::string("UPPER") + kAxis[c]);
    grid.nCellsAxis[c] = (std::size_t)h.whole(std::string("NCELL") + kAxis[c], kSize);
  }
  grid.limits = h.limits();

  try {
    internal::check_box_grid(grid, grid.mps, true, "");
  }
  catch (const Error& e) {
    throw Error(file + ": " + e.what());
  }

  const std::vector<std::string> names = internal::call_otswap([&file] { return otswap::io::columnNames(file); });
  const bool smoothed = Header::has_column(names, "divergenceSmoothed");
  const bool weightSmoothed = Header::has_column(names, "weightSmoothed");
  std::vector<std::string> columns = {"cellX", "cellY", "cellZ", "divergence"};
  if (smoothed) columns.push_back("divergenceSmoothed");
  columns.push_back("volume");
  columns.push_back("weight");
  if (weightSmoothed) columns.push_back("weightSmoothed");
  const otswap::io::Table t = internal::call_otswap([&] {
    return otswap::io::read(file, columns, {"index", "nDispl", "nTracers", "voidIndex", "overdenseIndex"});
  });

  const std::size_t n = grid.nCells;
  if (t.nRows != n)
    throw Error(file + ": it holds " + std::to_string(t.nRows) + " rows, but NCELLS is " + std::to_string(n));

  const std::size_t nc = columns.size();
  const std::size_t s = smoothed ? 1 : 0;
  grid.centre.resize(3 * n);
  grid.divergence.resize(n);
  if (smoothed) grid.divergenceSmoothed.resize(n);
  grid.volume.resize(n);
  grid.weight.resize(n);
  if (weightSmoothed) grid.weightSmoothed.resize(n);
  grid.nDisplacements.resize(n);
  grid.nTracers.resize(n);
  grid.voidIndex.resize(n);
  grid.overdenseIndex.resize(n);
  for (std::size_t c = 0; c < n; ++c) {
    const double* v = &t.values[nc*c];
    const std::int64_t* k = &t.integers[5*c];
    if (k[0] != (std::int64_t)c)
      throw Error(file + ": data row " + std::to_string(c+1) + " (from 1) has index " + std::to_string(k[0]) +
                  "; the rows must hold the cells 0 to NCELLS - 1 in order");
    if (k[1] < 0 || k[2] < 0)
      throw Error(file + ": data row " + std::to_string(c+1) + " (from 1) has a negative count (nDispl " +
                  std::to_string(k[1]) + ", nTracers " + std::to_string(k[2]) + ")");
    if (k[3] < -1 || k[4] < -1)
      throw Error(file + ": data row " + std::to_string(c+1) + " (from 1) has a basin index below -1 (voidIndex " +
                  std::to_string(k[3]) + ", overdenseIndex " + std::to_string(k[4]) + ")");
    for (std::size_t a = 0; a < 3; ++a) grid.centre[3*c+a] = v[a];
    grid.divergence[c] = v[3];
    if (smoothed) grid.divergenceSmoothed[c] = v[4];
    grid.volume[c] = v[4+s];
    grid.weight[c] = v[5+s];
    if (weightSmoothed) grid.weightSmoothed[c] = v[6+s];
    grid.nDisplacements[c] = (std::uint64_t)k[1];
    grid.nTracers[c] = (std::uint64_t)k[2];
    grid.voidIndex[c] = k[3];
    grid.overdenseIndex[c] = k[4];
  }

  return grid;
}


// ============================================================================


bitvf::io::Field bitvf::io::readField (const std::string& file, const std::array<std::string, 3>& centreColumns,
                                       const std::string& valueColumn, const std::string& weightColumn,
                                       const char delimiter, const char comment)
{
  std::vector<std::string> columns = {centreColumns[0], centreColumns[1], centreColumns[2], valueColumn};
  const bool weighted = !weightColumn.empty();
  if (weighted) columns.push_back(weightColumn);
  const otswap::io::Table t = internal::call_otswap([&] { return otswap::io::read(file, columns, delimiter, comment); });

  const std::size_t n = t.nRows;
  const std::size_t nc = columns.size();
  for (std::size_t i = 0; i < nc * n; ++i)
    if (!std::isfinite(t.values[i])) {
      const double v = t.values[i];
      throw Error(file + ": data row " + std::to_string(i / nc + 1) + " (from 1) has " +
                  (v != v ? "nan" : v > 0. ? "inf" : "-inf") + " in column " + columns[i % nc] + "; every value must be finite");
    }
  Field f;
  f.centres.resize(3 * n);
  f.values.resize(n);
  if (weighted) f.weights.resize(n);
  for (std::size_t c = 0; c < n; ++c) {
    const double* v = &t.values[nc*c];
    for (std::size_t a = 0; a < 3; ++a) f.centres[3*c+a] = v[a];
    f.values[c] = v[3];
    if (weighted) f.weights[c] = v[4];
  }

  return f;
}
