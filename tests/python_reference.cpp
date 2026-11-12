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
 *  @file tests/python_reference.cpp
 *
 *  @brief Runs the C++ library on inputs written to raw files, and writes
 *  its outputs to raw files and its report to a text file, so that the
 *  Python tests can check that the bindings return and print exactly what
 *  the library returns and prints.
 *
 *  Not a test by itself. Usage, with key=value arguments:
 *
 *    python_reference box tracers=F displacement=F n_realizations=R [valid=F] mps=X
 *                         [box=F] [ids=F] [reuse=1] [write=DIR] CONFIG out=P
 *    python_reference reconstruction complete=FITS CONFIG out=P
 *    python_reference divergence grid=FILE tracers=F [ids=F] CONFIG out=P
 *    python_reference field centres=F values=F [weights=F] [mps=X] CONFIG out=P
 *    python_reference mean file=FILE CONFIG out=P
 *    python_reference read_field file=FILE x=NAME y=NAME z=NAME value=NAME [weight=NAME] CONFIG out=P
 *    python_reference watershed field=F nx=N ny=N nz=N verbosity=V out=P
 *    python_reference grid_box box=F mps=X cell=X out=P
 *    python_reference grid_reconstruction complete=FITS [cell=X] out=P
 *    python_reference run tracers=F [randoms=F] [ids=F] n_realizations=R seed=N [box=F]
 *                         [run_mps=X] [keep_reconstructions=1] [axis=A redshift=X bias=X
 *                         rsd_sigma=X zspace=0|1 omega_m=X h=X z_max=X] CONFIG out=P
 *
 *  CONFIG is sigma=X unit=mps|coordinate [cell=X] keep=0|1 nrecmean=N
 *  verbosity=silent|normal|detailed; keep sets keepDivergence and
 *  keepAssignment (keepAssignment not for field and read_field). box
 *  reuse=1 searches twice, the second time on the grid of the first,
 *  which keeps its divergence and prints nothing;
 *  write=DIR writes the four products of the result, FITS and ASCII, into
 *  DIR. A run with axis corrects redshift-space distortions with the
 *  table DistanceTable(omega_m, h, -1, 0, 0, z_max, 20000).
 *
 *  Input files hold native float64 values, flat and row-major, ids native
 *  int64. A search writes P.report, the text of its report; P.scalars,
 *  one "name value" line per scalar, doubles in hexadecimal; and P.<name>
 *  for each array of the Result, raw, in its C++ type: voids.<column>,
 *  overdense.<column>, grid.<field> with Config::keepDivergence,
 *  assignment.<field> with Config::keepAssignment, ids, uncorrected,
 *  outside_grid, outside_cut. A run writes the same for the voids under
 *  P.voids, for the redshift-space voids under P.zspace, and the
 *  displacement of the kept reconstruction as P.reconstruction.displacement.
 *  watershed writes P.partition.<field> and P.report; the grids write
 *  P.scalars and P.grid.<field>.
 */

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "otswap/OT.h"
#include "otswap/io.h"

#include "bitvf/BitVF.h"
#include "bitvf/VoidFinder.h"
#include "bitvf/io.h"

namespace {

  std::map<std::string, std::string> g_args;

  bool has (const std::string& key)
  {
    return g_args.count(key) != 0;
  }

  std::string arg (const std::string& key)
  {
    if (!has(key)) throw std::runtime_error("missing argument " + key);
    return g_args.at(key);
  }

  double number (const std::string& key)
  {
    return std::stod(arg(key));
  }

  template <typename T>
  std::vector<T> load (const std::string& key)
  {
    std::ifstream in(arg(key), std::ios::binary);
    std::vector<T> v;
    T x;
    while (in.read(reinterpret_cast<char*>(&x), sizeof(T))) v.push_back(x);
    return v;
  }

  std::string g_out;

  template <typename T>
  void dump (const std::string& name, const std::vector<T>& v)
  {
    std::ofstream out(g_out + "." + name, std::ios::binary);
    out.write(reinterpret_cast<const char*>(v.data()), (std::streamsize)(v.size() * sizeof(T)));
  }

  std::ofstream g_scalars;

  void scalar (const std::string& name, const double v)
  {
    char text[64];
    std::snprintf(text, sizeof(text), "%a", v);
    g_scalars << name << " " << text << "\n";
  }

  void scalar (const std::string& name, const std::uint64_t v)
  {
    g_scalars << name << " " << v << "\n";
  }

  otswap::BoxLimits box_of (const std::string& key)
  {
    const std::vector<double> b = load<double>(key);
    otswap::BoxLimits box;
    for (std::size_t k = 0; k < 3; ++k) {
      box.lower[k] = b[k];
      box.upper[k] = b[3 + k];
    }
    return box;
  }

  otswap::Verbosity verbosity ()
  {
    const std::string v = arg("verbosity");
    if (v == "silent") return otswap::Verbosity::Silent;
    if (v == "detailed") return otswap::Verbosity::Detailed;
    return otswap::Verbosity::Normal;
  }

  bitvf::Config config (std::ostream* log, const bool assignment = true)
  {
    bitvf::Config c;
    c.smoothing.sigma = number("sigma");
    c.smoothing.unit = arg("unit") == "mps" ? bitvf::LengthUnit::Mps : bitvf::LengthUnit::Coordinate;
    if (has("cell")) c.divergenceCellSize = number("cell");
    const bool keep = arg("keep") == "1";
    c.keepDivergence = keep;
    c.keepAssignment = keep && assignment;
    c.nRealizationsMean = (unsigned)number("nrecmean");
    c.verbosity = verbosity();
    c.log = log;
    return c;
  }

  void dump_catalog (const std::string& prefix, const bitvf::BasinCatalog& c)
  {
    scalar(prefix + ".n_basins", (std::uint64_t)c.nBasins);
    dump(prefix + ".centre", c.centre);
    dump(prefix + ".centre_sky", c.centreSky);
    dump(prefix + ".radius", c.radius);
    dump(prefix + ".radius_ext", c.radiusExt);
    dump(prefix + ".delta_av", c.deltaAv);
    dump(prefix + ".delta_2mps", c.delta2Mps);
    dump(prefix + ".edge_flag", c.edgeFlag);
    dump(prefix + ".contamination", c.contamination);
    dump(prefix + ".n_tracers", c.nTracers);
    dump(prefix + ".field_centre", c.fieldCentre);
  }

  void dump_grid (const std::string& prefix, const bitvf::DivergenceGrid& g)
  {
    scalar(prefix + ".n_cells", (std::uint64_t)g.nCells);
    scalar(prefix + ".cell_size", g.cellSize);
    for (std::size_t k = 0; k < 3; ++k) {
      scalar(prefix + ".lower" + std::to_string(k), g.lower[k]);
      scalar(prefix + ".upper" + std::to_string(k), g.upper[k]);
      scalar(prefix + ".n_cells_axis" + std::to_string(k), (std::uint64_t)g.nCellsAxis[k]);
    }
    scalar(prefix + ".n_realizations", (std::uint64_t)g.nRealizations);
    scalar(prefix + ".n_realizations_mean", (std::uint64_t)g.nRealizationsMean);
    scalar(prefix + ".cell_size_mps", g.cellSizeMps);
    scalar(prefix + ".mps", g.mps);
    scalar(prefix + ".has_limits", (std::uint64_t)g.limits.has_value());
    dump(prefix + ".centre", g.centre);
    dump(prefix + ".divergence", g.divergence);
    dump(prefix + ".divergence_smoothed", g.divergenceSmoothed);
    dump(prefix + ".volume", g.volume);
    dump(prefix + ".weight", g.weight);
    dump(prefix + ".weight_smoothed", g.weightSmoothed);
    dump(prefix + ".n_displacements", g.nDisplacements);
    dump(prefix + ".n_tracers", g.nTracers);
    dump(prefix + ".void_index", g.voidIndex);
    dump(prefix + ".overdense_index", g.overdenseIndex);
  }

  void dump_result (const std::string& prefix, const bitvf::Result& r)
  {
    const std::string p = prefix.empty() ? "" : prefix + ".";
    scalar(p + "n_tracers", (std::uint64_t)r.nTracers);
    scalar(p + "n_realizations", (std::uint64_t)r.nRealizations);
    scalar(p + "mps", r.mps);
    scalar(p + "random_seed", (std::uint64_t)(r.randomSeed ? *r.randomSeed + 1ull : 0ull));
    scalar(p + "cell_size", r.config.divergenceCellSize.value_or(std::numeric_limits<double>::quiet_NaN()));
    scalar(p + "has_limits", (std::uint64_t)r.limits.has_value());
    const bitvf::Counts& c = r.counts;
    scalar(p + "counts.tracers", (std::uint64_t)c.tracers);
    scalar(p + "counts.tracers_used", (std::uint64_t)c.tracersUsed);
    scalar(p + "counts.tracers_uncorrected", (std::uint64_t)c.tracersUncorrected);
    scalar(p + "counts.tracers_outside_box", (std::uint64_t)c.tracersOutsideBox);
    scalar(p + "counts.tracers_outside_grid", (std::uint64_t)c.tracersOutsideGrid);
    scalar(p + "counts.displacements", (std::uint64_t)c.displacements);
    dump_catalog(p + "voids", r.voids);
    dump_catalog(p + "overdense", r.overdense);
    if (r.divergence) dump_grid(p + "grid", *r.divergence);
    if (r.assignment) {
      dump(p + "assignment.cell_index", r.assignment->cellIndex);
      dump(p + "assignment.void_index", r.assignment->voidIndex);
      dump(p + "assignment.overdense_index", r.assignment->overdenseIndex);
      dump(p + "assignment.distance", r.assignment->distance);
      dump(p + "assignment.in_mask", r.assignment->inMask);
      dump(p + "assignment.in_box", r.assignment->inBox);
      dump(p + "assignment.used", r.assignment->used);
    }
    dump(p + "ids", r.ids);
    dump(p + "uncorrected", r.uncorrected);
    dump(p + "outside_grid", r.outsideGrid);
    dump(p + "outside_cut", r.outsideCut);
  }

  void write_products (const std::string& directory, const bitvf::Result& r)
  {
    for (const std::string extension : {".fits", ".dat"}) {
      bitvf::io::writeVoids(directory + "/voids" + extension, r);
      bitvf::io::writeOverdense(directory + "/overdense" + extension, r);
      bitvf::io::writeDivergence(directory + "/divergence" + extension, r);
      bitvf::io::writeGalaxies(directory + "/galaxies" + extension, r);
    }
  }

  void run (const std::string& mode, std::ostream& log)
  {
    if (mode == "box") {
      const std::vector<double> tracers = load<double>("tracers");
      const std::vector<double> displacement = load<double>("displacement");
      const std::vector<std::uint8_t> valid = has("valid") ? load<std::uint8_t>("valid") : std::vector<std::uint8_t>();
      const std::optional<otswap::BoxLimits> limits = has("box") ? std::optional<otswap::BoxLimits>(box_of("box"))
                                                                 : std::nullopt;
      const std::vector<std::int64_t> ids = has("ids") ? load<std::int64_t>("ids") : std::vector<std::int64_t>();
      const bitvf::Config c = config(&log);
      std::optional<bitvf::Result> first;
      if (has("reuse") && arg("reuse") == "1") {
        bitvf::Config quiet = c;
        quiet.keepDivergence = true;
        quiet.verbosity = otswap::Verbosity::Silent;
        first = bitvf::findVoidsBox(tracers, displacement, valid, number("mps"), quiet, limits, ids);
      }
      const bitvf::Result r = bitvf::findVoidsBox(tracers, displacement, valid, number("mps"), c, limits, ids,
                                                  first ? &*first->divergence : nullptr);
      dump_result("", r);
      if (has("write")) write_products(arg("write"), r);
    }
    else if (mode == "reconstruction") {
      const otswap::Result rec = otswap::io::readDisplacementComplete(arg("complete"));
      dump_result("", bitvf::findVoidsBox(rec, config(&log)));
    }
    else if (mode == "divergence") {
      const bitvf::DivergenceGrid grid = bitvf::io::readDivergence(arg("grid"));
      const std::vector<std::int64_t> ids = has("ids") ? load<std::int64_t>("ids") : std::vector<std::int64_t>();
      dump_result("", bitvf::findVoidsBoxFromDivergence(grid, load<double>("tracers"), config(&log), ids));
    }
    else if (mode == "field") {
      const std::vector<double> weights = has("weights") ? load<double>("weights") : std::vector<double>();
      const std::optional<double> mps = has("mps") ? std::optional<double>(number("mps")) : std::nullopt;
      dump_result("", bitvf::findVoidsBoxFromField(load<double>("centres"), load<double>("values"),
                                                    config(&log, false), weights, mps));
    }
    else if (mode == "mean") {
      const otswap::Result m = otswap::io::readDisplacementMean(arg("file"));
      bitvf::Config c = config(&log);
      c.nRealizationsMean = m.nRealizationsMean;
      dump_result("", bitvf::findVoidsBox(m, c));
    }
    else if (mode == "read_field") {
      const bitvf::io::Field f = bitvf::io::readField(arg("file"), {arg("x"), arg("y"), arg("z")}, arg("value"),
                                                       has("weight") ? arg("weight") : "");
      dump_result("", bitvf::findVoidsBoxFromField(f.centres, f.values, config(&log, false), f.weights));
    }
    else if (mode == "watershed") {
      const bitvf::Partition p = bitvf::watershedBox(load<double>("field"),
                                                     {(std::size_t)number("nx"), (std::size_t)number("ny"),
                                                      (std::size_t)number("nz")}, verbosity(), &log);
      scalar("partition.n_voids", (std::uint64_t)p.nVoids);
      scalar("partition.n_overdense", (std::uint64_t)p.nOverdense);
      dump("partition.void_index", p.voidIndex);
      dump("partition.overdense_index", p.overdenseIndex);
      dump("partition.void_seed", p.voidSeed);
      dump("partition.overdense_seed", p.overdenseSeed);
      dump("partition.void_cells", p.voidCells);
      dump("partition.void_core_cells", p.voidCoreCells);
      dump("partition.overdense_cells", p.overdenseCells);
      dump("partition.overdense_core_cells", p.overdenseCoreCells);
    }
    else if (mode == "grid_box") {
      dump_grid("grid", bitvf::boxGrid(box_of("box"), number("mps"), number("cell")));
    }
    else if (mode == "grid_reconstruction") {
      const otswap::Result rec = otswap::io::readDisplacementComplete(arg("complete"));
      dump_grid("grid", bitvf::boxGrid(rec, has("cell") ? std::optional<double>(number("cell")) : std::nullopt));
    }
    else if (mode == "run") {
      bitvf::BoxConfig c;
      c.ot.nRealizations = (unsigned)number("n_realizations");
      c.ot.seed = (unsigned)number("seed");
      c.ot.verbosity = verbosity();
      c.ot.log = &log;
      if (has("box")) c.limits = box_of("box");
      if (has("run_mps")) c.mps = number("run_mps");
      c.voids = config(&log);
      c.keepReconstructions = has("keep_reconstructions");
      const std::vector<double> randoms = has("randoms") ? load<double>("randoms") : std::vector<double>();
      const std::vector<std::int64_t> ids = has("ids") ? load<std::int64_t>("ids") : std::vector<std::int64_t>();
      bitvf::BitVFResult r;
      if (has("axis")) {
        bitvf::BoxRsd rsd;
        rsd.axis = (unsigned)number("axis");
        rsd.redshift = number("redshift");
        rsd.bias = number("bias");
        rsd.sigma = number("rsd_sigma");
        rsd.zspaceVoids = arg("zspace") == "1";
        rsd.correction.verbosity = verbosity();
        rsd.correction.log = &log;
        c.rsd = rsd;
        const otswap::DistanceTable table(number("omega_m"), number("h"), -1., 0., 0., number("z_max"), 20000);
        r = bitvf::runBox(load<double>("tracers"), randoms, table, c, ids);
      }
      else
        r = bitvf::runBox(load<double>("tracers"), randoms, c, ids);
      dump_result("voids", r.voids);
      if (r.zspaceVoids) dump_result("zspace", *r.zspaceVoids);
      if (r.reconstruction) dump("reconstruction.displacement", r.reconstruction->displacement);
    }
    else
      throw std::runtime_error("unknown mode " + mode);
  }

}

int main (int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "usage: python_reference <mode> key=value ... out=P\n";
    return 2;
  }
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    const std::size_t eq = a.find('=');
    if (eq == std::string::npos) {
      std::cerr << "python_reference: not key=value: " << a << "\n";
      return 2;
    }
    g_args[a.substr(0, eq)] = a.substr(eq + 1);
  }
  try {
    g_out = arg("out");
    g_scalars.open(g_out + ".scalars");
    std::ostringstream log;
    run(argv[1], log);
    std::ofstream(g_out + ".report") << log.str();
  }
  catch (const std::exception& e) {
    std::cerr << "python_reference: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
