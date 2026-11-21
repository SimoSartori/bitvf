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
 *  @file python/bitvf_python.cpp
 *
 *  @brief nanobind module bitvf._bitvf: the bindings of
 *  include/bitvf/VoidFinder.h and include/bitvf/BitVF.h, box geometry, as
 *  specified by python/bitvf/__init__.pyi, and of include/bitvf/io.h, as
 *  the submodule bitvf.io specified by python/bitvf/io.pyi.
 *
 *  Every argument is taken as a Python object and converted here, so that
 *  any invalid argument raises bitvf.Error; every named argument accepts
 *  None, so that None reaches the conversions and raises bitvf.Error like
 *  any other invalid value. Arrays are converted with numpy.asarray to
 *  C-contiguous float64 (int64 for ids) and their shapes checked; the C++
 *  library then works on its own copy. otswap's types (Config,
 *  CorrectionConfig, BoxLimits, DistanceTable, MpsProfile, Result,
 *  RealSpaceCatalog) are otswap's Python types themselves, bound by otswap
 *  and imported first.
 *
 *  Smoothing, Config, BoxRsd and BoxConfig are the C++ structs; their
 *  setters convert each value to its C++ type, and the library checks it
 *  at the call. A nested configuration is read as a reference into the
 *  object that holds it, so that setting one of its fields changes the
 *  holder; assigning one copies it in. The results are the C++ structs,
 *  read-only and without a constructor: their arrays are read-only views
 *  on the C++ object, kept alive by the Python object that holds it; ids is
 *  a new copy on each access. Every per-cell array of a box grid has shape
 *  (nx, ny, nz), C order, so that a[i, j, k] is row (i ny + j) nz + k; 0/1
 *  flags are bool.
 *
 *  Each call converts its arguments, releases the GIL and runs the
 *  library with its report written to a stream whose buffer, on each
 *  complete line, takes the GIL, checks for signals and prints the line
 *  with builtins.print, so the report appears on sys.stdout as the call
 *  runs. A signal handler that raises, as Ctrl+C's does, stops the call at
 *  its next line and the call raises that exception; any other exception
 *  raised by print stops the printing, and the call raises it once the
 *  library returns, its result discarded.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <numeric>
#include <optional>
#include <ostream>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include "bitvf/BitVF.h"
#include "bitvf/VoidFinder.h"
#include "bitvf/io.h"

#include "internal.h"

namespace nb = nanobind;
using namespace nb::literals;

static_assert(sizeof(bool) == 1 && sizeof(std::uint8_t) == 1, "flags are exposed as bool");
static_assert(sizeof(std::size_t) == 8, "counts are exposed as uint64");

namespace {

  using CArray = nb::ndarray<const double, nb::c_contig, nb::device::cpu>;

  std::string repr_of (nb::handle value)
  {
    return nb::repr(value).c_str();
  }

  std::string shape_of (const CArray& a)
  {
    std::string s = "(";
    for (std::size_t i = 0; i < a.ndim(); ++i)
      s += (i ? ", " : "") + std::to_string(a.shape(i));
    return s + (a.ndim() == 1 ? ",)" : ")");
  }

  // ==========================================================================
  // Argument conversion
  // ==========================================================================

  /// The argument as a C-contiguous float64 numpy array of the same shape.
  CArray as_float64 (nb::handle value, const std::string& name)
  {
    const nb::module_ np = nb::module_::import_("numpy");
    nb::object a;
    try {
      a = np.attr("asarray")(value, "dtype"_a = np.attr("float64"), "order"_a = "C");
    }
    catch (nb::python_error& e) {
      throw bitvf::Error(name + " cannot be converted to a float64 array: " + nb::str(e.value()).c_str());
    }
    return nb::cast<CArray>(a);
  }

  /// An (N, 3) array, copied into the flat layout the library takes.
  std::vector<double> rows3 (nb::handle value, const std::string& name)
  {
    const CArray a = as_float64(value, name);
    if (a.ndim() != 2 || a.shape(1) != 3)
      throw bitvf::Error(name + " must have shape (N, 3); it has shape " + shape_of(a));
    return std::vector<double>(a.data(), a.data() + a.size());
  }

  /// A one-dimensional array.
  std::vector<double> column (nb::handle value, const std::string& name)
  {
    const CArray a = as_float64(value, name);
    if (a.ndim() != 1)
      throw bitvf::Error(name + " must be one-dimensional; it has shape " + shape_of(a));
    return std::vector<double>(a.data(), a.data() + a.size());
  }

  double to_double (nb::handle value, const std::string& name)
  {
    try {
      return nb::cast<double>(value);
    }
    catch (nb::cast_error&) {
      throw bitvf::Error(name + " must be a real number; got " + repr_of(value));
    }
  }

  std::optional<double> to_optional_double (nb::handle value, const std::string& name)
  {
    if (value.is_none()) return std::nullopt;
    return to_double(value, name);
  }

  unsigned to_unsigned (nb::handle value, const std::string& name)
  {
    long long v;
    try {
      v = nb::cast<long long>(value);
    }
    catch (nb::cast_error&) {
      throw bitvf::Error(name + " must be an integer; got " + repr_of(value));
    }
    if (v < 0)
      throw bitvf::Error(name + " is " + std::to_string(v) + "; it must be a non-negative integer");
    if (v > 0xffffffffLL)
      throw bitvf::Error(name + " is " + std::to_string(v) + "; it must lie in [0, 4294967295]");
    return (unsigned)v;
  }

  bool to_bool (nb::handle value, const std::string& name)
  {
    if (!PyBool_Check(value.ptr()))
      throw bitvf::Error(name + " must be True or False; got " + repr_of(value));
    return value.ptr() == Py_True;
  }

  std::string to_string (nb::handle value, const std::string& name)
  {
    if (!nb::isinstance<nb::str>(value))
      throw bitvf::Error(name + " must be a str; got " + repr_of(value));
    return nb::str(value).c_str();
  }

  /// A path: a str or an os.PathLike.
  std::string to_path (nb::handle value, const std::string& name)
  {
    nb::object path;
    try {
      path = nb::module_::import_("os").attr("fspath")(value);
    }
    catch (nb::python_error&) {
      throw bitvf::Error(name + " must be a str or an os.PathLike; got " + repr_of(value));
    }
    if (!nb::isinstance<nb::str>(path))
      throw bitvf::Error(name + " must be a str or an os.PathLike naming a str; got " + repr_of(value));
    return nb::str(path).c_str();
  }

  char to_char (nb::handle value, const std::string& name)
  {
    const std::string s = to_string(value, name);
    if (s.size() != 1)
      throw bitvf::Error(name + " must be a single character; got " + repr_of(value));
    return s[0];
  }

  /// An instance of a bound class, named with its module in type; the
  /// reference stays valid while the caller holds the argument.
  template <typename T>
  T& instance (nb::handle value, const std::string& name, const std::string& type)
  {
    if (!nb::isinstance<T>(value))
      throw bitvf::Error(name + " must be " + type + "; got " + repr_of(value));
    return nb::cast<T&>(value);
  }

  bool in_degrees (nb::handle unit)
  {
    if (nb::isinstance<nb::str>(unit)) {
      const std::string u = nb::str(unit).c_str();
      if (u == "deg") return true;
      if (u == "rad") return false;
    }
    throw bitvf::Error("angle_unit must be \"deg\" or \"rad\"; got " + repr_of(unit));
  }

  otswap::Verbosity to_verbosity (nb::handle value)
  {
    if (nb::isinstance<nb::str>(value)) {
      const std::string v = nb::str(value).c_str();
      if (v == "silent") return otswap::Verbosity::Silent;
      if (v == "normal") return otswap::Verbosity::Normal;
      if (v == "detailed") return otswap::Verbosity::Detailed;
    }
    throw bitvf::Error("verbosity must be \"silent\", \"normal\" or \"detailed\"; got " + repr_of(value));
  }

  nb::str verbosity_name (const otswap::Verbosity v)
  {
    if (v == otswap::Verbosity::Silent) return nb::str("silent");
    if (v == otswap::Verbosity::Detailed) return nb::str("detailed");
    return nb::str("normal");
  }

  /// A NotImplemented for the comparison operators.
  nb::object not_implemented ()
  {
    return nb::borrow(nb::handle(Py_NotImplemented));
  }

  /// A one-dimensional array of 64-bit integers: integer or boolean values
  /// taken as they are, floating ones only when each is an integer in range.
  std::vector<std::int64_t> int64_column (nb::handle value, const std::string& name)
  {
    const nb::module_ np = nb::module_::import_("numpy");
    nb::object raw;
    try {
      raw = np.attr("asarray")(value);
    }
    catch (nb::python_error& e) {
      throw bitvf::Error(name + " cannot be converted to an array: " + nb::str(e.value()).c_str());
    }
    const std::string kind = nb::str(raw.attr("dtype").attr("kind")).c_str();
    if (nb::cast<std::size_t>(raw.attr("ndim")) != 1)
      throw bitvf::Error(name + " must be one-dimensional; it has shape " +
                         repr_of(nb::object(raw.attr("shape"))));
    if (kind == "f") {
      const std::vector<double> c = column(raw, name);
      std::vector<std::int64_t> out(c.size());
      for (std::size_t i = 0; i < c.size(); ++i) {
        if (!(c[i] >= -9223372036854775808. && c[i] < 9223372036854775808.) || c[i] != std::floor(c[i]))
          throw bitvf::Error(name + " holds " + std::to_string(c[i]) + " for object " + std::to_string(i) +
                             " (from 0); every entry must be a 64-bit integer");
        out[i] = (std::int64_t)c[i];
      }
      return out;
    }
    if (kind != "i" && kind != "u" && kind != "b")
      throw bitvf::Error(name + " must hold integers; its dtype is " + repr_of(raw.attr("dtype")));
    if (kind == "u" && nb::cast<std::size_t>(raw.attr("size")) > 0 &&
        nb::cast<unsigned long long>(raw.attr("max")()) > 9223372036854775807ULL)
      throw bitvf::Error(name + " holds a value above 9223372036854775807, the largest 64-bit integer");
    using IArray = nb::ndarray<const std::int64_t, nb::c_contig, nb::device::cpu>;
    const IArray a = nb::cast<IArray>(np.attr("ascontiguousarray")(raw, "dtype"_a = np.attr("int64")));
    return std::vector<std::int64_t>(a.data(), a.data() + a.size());
  }

  /// The IDs to pass for n objects: none for None, an empty array or
  /// exactly 0..n-1, else the values, which the library checks.
  std::vector<std::int64_t> ids_value (nb::handle value, const std::size_t n)
  {
    if (value.is_none()) return {};
    std::vector<std::int64_t> ids = int64_column(value, "ids");
    if (ids.size() != n) return ids;
    for (std::size_t i = 0; i < ids.size(); ++i)
      if (ids[i] != (std::int64_t)i) return ids;
    return {};
  }

  std::string dims (const std::vector<std::size_t>& shape)
  {
    std::string s = "(";
    for (std::size_t i = 0; i < shape.size(); ++i) s += (i ? ", " : "") + std::to_string(shape[i]);
    return s + (shape.size() == 1 ? ",)" : ")");
  }

  /// 0 or 1 flags of the given shape, from bool or numbers already
  /// converted; the shape error names what the shape follows.
  std::vector<std::uint8_t> flags (const CArray& a, const std::string& name,
                                   const std::vector<std::size_t>& shape, const std::string& as)
  {
    bool same = a.ndim() == shape.size();
    for (std::size_t i = 0; same && i < shape.size(); ++i) same = a.shape(i) == shape[i];
    if (!same)
      throw bitvf::Error(name + " must have shape " + dims(shape) + ", " + as + "; it has shape " + shape_of(a));
    std::vector<std::uint8_t> out(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      const double v = a.data()[i];
      if (v != 0. && v != 1.)
        throw bitvf::Error("the validity flags hold " + std::to_string(v) + " for object " + std::to_string(i) +
                           " (from 0); each must be 0 or 1");
      out[i] = (std::uint8_t)v;
    }
    return out;
  }

  /// The box limits from None, an otswap.BoxLimits, or a pair of corners
  /// ((x0, y0, z0), (x1, y1, z1)). The values are checked by the library.
  std::optional<otswap::BoxLimits> to_box (nb::handle value, const std::string& name)
  {
    if (value.is_none()) return std::nullopt;
    if (nb::isinstance<otswap::BoxLimits>(value)) return nb::cast<otswap::BoxLimits>(value);
    const CArray a = as_float64(value, name);
    if (a.ndim() != 2 || a.shape(0) != 2 || a.shape(1) != 3)
      throw bitvf::Error(name + " must be an otswap.BoxLimits or a pair of corners ((x0, y0, z0), "
                         "(x1, y1, z1)); it has shape " + shape_of(a));
    otswap::BoxLimits box;
    for (std::size_t k = 0; k < 3; ++k) {
      box.lower[k] = a.data()[k];
      box.upper[k] = a.data()[3 + k];
    }
    return box;
  }

  template <typename T>
  nb::object triple (const std::array<T, 3>& c)
  {
    return nb::make_tuple(c[0], c[1], c[2]);
  }

  nb::object box_or_none (const std::optional<otswap::BoxLimits>& box)
  {
    if (!box) return nb::none();
    return nb::cast(*box);
  }

  bool same_box (const std::optional<otswap::BoxLimits>& a, const std::optional<otswap::BoxLimits>& b)
  {
    if (!a || !b) return !a && !b;
    return a->lower == b->lower && a->upper == b->upper;
  }

  nb::str optional_repr (const std::optional<double>& v)
  {
    if (!v) return nb::str("None");
    return nb::repr(nb::float_(*v));
  }

  // ==========================================================================
  // Returned arrays
  // ==========================================================================

  template <typename T>
  using View = nb::ndarray<nb::numpy, const T>;

  /// A read-only view of the given shape on the elements of a vector of
  /// owner, reinterpreted as T, which keeps owner alive.
  template <typename T, typename V>
  View<T> view (nb::handle owner, const std::vector<V>& values, const std::vector<std::size_t>& shape)
  {
    static_assert(sizeof(T) == sizeof(V), "a view reinterprets the vector's elements");
    return View<T>(reinterpret_cast<const T*>(values.data()), shape.size(), shape.data(), owner);
  }

  /// A column, shape (n,).
  template <typename T, typename V>
  View<T> rows (nb::handle owner, const std::vector<V>& values)
  {
    return view<T>(owner, values, {values.size()});
  }

  /// Triplets, shape (n, 3).
  View<double> triplets (nb::handle owner, const std::vector<double>& values)
  {
    return view<double>(owner, values, {values.size() / 3, 3});
  }

  /// A per-cell array of a grid of n cells along each axis, width values
  /// per cell: shape (nx, ny, nz) or (nx, ny, nz, width); (0,) or
  /// (0, width) when empty; flat, (cells,) or (cells, width), when the
  /// counts do not describe it.
  template <typename T, typename V>
  View<T> cells (nb::handle owner, const std::vector<V>& values, const std::array<std::size_t, 3>& n,
                 const std::size_t width = 1)
  {
    std::vector<std::size_t> shape;
    if (!values.empty() && n[0] * n[1] * n[2] * width == values.size()) shape = {n[0], n[1], n[2]};
    else shape = {values.size() / width};
    if (width != 1) shape.push_back(width);
    return view<T>(owner, values, shape);
  }

  /// A new read-only array of the IDs of n objects: a copy of ids, or 0..n-1
  /// when it is empty.
  View<std::int64_t> ids_array (const std::vector<std::int64_t>& ids, const std::size_t n)
  {
    auto* held = new std::vector<std::int64_t>(ids);
    if (held->empty()) {
      held->resize(n);
      std::iota(held->begin(), held->end(), std::int64_t(0));
    }
    nb::capsule owner(held, [] (void* p) noexcept { delete static_cast<std::vector<std::int64_t>*>(p); });
    return View<std::int64_t>(held->data(), {held->size()}, owner);
  }

  // ==========================================================================
  // Live report
  // ==========================================================================

  /// Thrown by the report's buffer to stop the call when an exception that
  /// stops Python, a raising signal handler or a KeyboardInterrupt from
  /// print, has been stored. Not a std::exception, so that no handler of
  /// the library or otswap catches it.
  struct Interrupted {};

  /**
   *  @brief The buffer of a call's report: on each complete line it takes
   *  the GIL, runs the pending signal handlers, then prints the line,
   *  without its '\n', with builtins.print.
   *
   *  A handler that raises stores its exception and throws Interrupted,
   *  as does a KeyboardInterrupt raised by print. Any other exception
   *  raised by print is stored, and nothing more is printed. Lines come
   *  from the calling thread only, which holds no GIL while the library
   *  runs.
   */
  class LineBuffer : public std::streambuf {
  public:

    /// The first exception stored, or null.
    nb::object error;

    /// Print what is left after the last '\n'. Needs the GIL.
    void finish ()
    {
      if (!error.is_valid() && !line.empty()) emit();
      line.clear();
    }

  protected:

    int_type overflow (const int_type c) override
    {
      if (!traits_type::eq_int_type(c, traits_type::eof())) put(traits_type::to_char_type(c));
      return traits_type::not_eof(c);
    }

    std::streamsize xsputn (const char* s, const std::streamsize n) override
    {
      for (std::streamsize i = 0; i < n; ++i) put(s[i]);
      return n;
    }

  private:

    std::string line;

    void put (const char c)
    {
      if (c != '\n') {
        line += c;
        return;
      }
      emit();
      line.clear();
    }

    void emit ()
    {
      if (error.is_valid()) return;
      bool interrupted = false;
      {
        nb::gil_scoped_acquire held;
        if (PyErr_CheckSignals() != 0) {
          nb::python_error raised;
          error = nb::borrow(raised.value());
          interrupted = true;
        }
        else {
          try {
            nb::module_::import_("builtins").attr("print")(nb::str(line.data(), line.size()));
          }
          catch (nb::python_error& raised) {
            error = nb::borrow(raised.value());
            interrupted = raised.matches(PyExc_KeyboardInterrupt);
          }
        }
      }
      if (interrupted) throw Interrupted();
    }
  };

  /// The report of one call: a stream over a LineBuffer that passes on
  /// the exceptions its buffer throws.
  class LiveReport {
  public:

    LiveReport () : stream(&buffer)
    {
      stream.exceptions(std::ios::badbit);
    }

    LiveReport (const LiveReport&) = delete;
    LiveReport& operator= (const LiveReport&) = delete;

    std::ostream* log () { return &stream; }

    /// Run f with the GIL released, and return its value. An exception
    /// stored by the buffer is raised after f returns or throws, in place
    /// of its value or its exception.
    template <typename F>
    auto run (F&& f) -> decltype(f())
    {
      std::optional<decltype(f())> value;
      std::exception_ptr failure;
      {
        nb::gil_scoped_release released;
        try {
          value.emplace(f());
        }
        catch (const Interrupted&) {
        }
        catch (...) {
          failure = std::current_exception();
        }
      }
      buffer.finish();
      if (buffer.error.is_valid()) {
        PyErr_SetObject((PyObject*)Py_TYPE(buffer.error.ptr()), buffer.error.ptr());
        throw nb::python_error();
      }
      if (failure) std::rethrow_exception(failure);
      return std::move(*value);
    }

  private:

    LineBuffer buffer;
    std::ostream stream;
  };

  // ==========================================================================
  // Results
  // ==========================================================================

  /// A Partition with the shape of the field it splits.
  struct PartitionObject {
    bitvf::Partition partition;
    std::array<std::size_t, 3> n;
  };

  /// A Result as its Python object, config.log reset to &std::clog, so
  /// that the stream of the call never outlives the call.
  nb::object result_object (bitvf::Result&& r)
  {
    r.config.log = &std::clog;
    return nb::cast(std::move(r), nb::rv_policy::move);
  }

  nb::object pipeline_object (bitvf::BitVFResult&& r)
  {
    r.voids.config.log = &std::clog;
    if (r.zspaceVoids) r.zspaceVoids->config.log = &std::clog;
    if (r.reconstruction) r.reconstruction->config.log = &std::clog;
    if (r.zspaceReconstruction) r.zspaceReconstruction->config.log = &std::clog;
    return nb::cast(std::move(r), nb::rv_policy::move);
  }

  using GridHandle = nb::handle_t<bitvf::DivergenceGrid>;
  using CatalogHandle = nb::handle_t<bitvf::BasinCatalog>;
  using AssignmentHandle = nb::handle_t<bitvf::TracerAssignment>;
  using ResultHandle = nb::handle_t<bitvf::Result>;
  using PartitionHandle = nb::handle_t<PartitionObject>;
  using FieldHandle = nb::handle_t<bitvf::io::Field>;

  const bitvf::DivergenceGrid& grid_of (GridHandle self) { return nb::cast<const bitvf::DivergenceGrid&>(self); }
  const bitvf::BasinCatalog& catalog_of (CatalogHandle self) { return nb::cast<const bitvf::BasinCatalog&>(self); }
  const bitvf::TracerAssignment& assignment_of (AssignmentHandle self)
  {
    return nb::cast<const bitvf::TracerAssignment&>(self);
  }
  const bitvf::Result& result_of (ResultHandle self) { return nb::cast<const bitvf::Result&>(self); }
  const PartitionObject& partition_of (PartitionHandle self) { return nb::cast<const PartitionObject&>(self); }

  nb::str geometry_name (const otswap::Geometry g)
  {
    return nb::str(g == otswap::Geometry::Box ? "box" : "lightcone");
  }

  // ==========================================================================
  // Configuration text
  // ==========================================================================

  std::string unit_name (const bitvf::LengthUnit u)
  {
    return u == bitvf::LengthUnit::Mps ? "LengthUnit.MPS" : "LengthUnit.COORDINATE";
  }

  std::string smoothing_repr (const bitvf::Smoothing& s)
  {
    return "Smoothing(sigma=" + repr_of(nb::float_(s.sigma)) + ", unit=" + unit_name(s.unit) + ")";
  }

  std::string config_repr (const bitvf::Config& c)
  {
    return "Config(smoothing=" + smoothing_repr(c.smoothing) +
      ", divergence_cell_size=" + std::string(optional_repr(c.divergenceCellSize).c_str()) +
      ", keep_divergence=" + (c.keepDivergence ? "True" : "False") +
      ", keep_assignment=" + (c.keepAssignment ? "True" : "False") +
      ", n_realizations_mean=" + std::to_string(c.nRealizationsMean) +
      ", verbosity=" + repr_of(verbosity_name(c.verbosity)) + ")";
  }

  std::string rsd_repr (const bitvf::BoxRsd& r)
  {
    return "BoxRsd(axis=" + (r.axis ? std::to_string(*r.axis) : std::string("None")) +
      ", redshift=" + repr_of(nb::float_(r.redshift)) + ", bias=" + repr_of(nb::float_(r.bias)) +
      ", sigma=" + repr_of(nb::float_(r.sigma)) + ", correction=" + repr_of(nb::cast(r.correction)) +
      ", zspace_voids=" + (r.zspaceVoids ? "True" : "False") + ")";
  }

  bool same_smoothing (const bitvf::Smoothing& a, const bitvf::Smoothing& b)
  {
    return a.sigma == b.sigma && a.unit == b.unit;
  }

  bool same_config (const bitvf::Config& a, const bitvf::Config& b)
  {
    return same_smoothing(a.smoothing, b.smoothing) && a.divergenceCellSize == b.divergenceCellSize &&
      a.keepDivergence == b.keepDivergence && a.keepAssignment == b.keepAssignment &&
      a.nRealizationsMean == b.nRealizationsMean && a.verbosity == b.verbosity;
  }

  /// Two otswap configurations equal by otswap's own __eq__.
  template <typename T>
  bool same_otswap (const T& a, const T& b)
  {
    return nb::cast(a).equal(nb::cast(b));
  }

  bool same_rsd (const bitvf::BoxRsd& a, const bitvf::BoxRsd& b)
  {
    return a.axis == b.axis && a.redshift == b.redshift && a.bias == b.bias && a.sigma == b.sigma &&
      same_otswap(a.correction, b.correction) && a.zspaceVoids == b.zspaceVoids;
  }

  // ==========================================================================
  // Calls
  // ==========================================================================

  /// The names given, as "a", "a and b" or "a, b and c".
  std::string listed (const std::vector<std::string>& names)
  {
    std::string s;
    for (std::size_t i = 0; i < names.size(); ++i)
      s += (i == 0 ? "" : i + 1 == names.size() ? " and " : ", ") + names[i];
    return s;
  }

  const bitvf::DivergenceGrid* grid_or_null (nb::handle grid)
  {
    if (grid.is_none()) return nullptr;
    return &instance<bitvf::DivergenceGrid>(grid, "grid", "a bitvf.DivergenceGrid");
  }

}

NB_MODULE(_bitvf, m)
{
  m.doc() = "Compiled part of the bitvf package; import bitvf instead.";

  nb::module_::import_("otswap");

  m.attr("__version__") = bitvf::internal::version;

  nb::exception<bitvf::Error> error(m, "Error", PyExc_Exception);
  error.attr("__doc__") = "Raised by every bitvf function on invalid input or failure.";

  // ==========================================================================
  // Configuration
  // ==========================================================================

  nb::enum_<bitvf::LengthUnit>(m, "LengthUnit", "Unit of a length the user chooses.")
    .value("MPS", bitvf::LengthUnit::Mps, "Multiples of the mean particle separation.")
    .value("COORDINATE", bitvf::LengthUnit::Coordinate, "The unit of the tracer coordinates.");

  auto to_unit = [] (nb::handle value) {
    if (!nb::isinstance<bitvf::LengthUnit>(value))
      throw bitvf::Error("unit must be a bitvf.LengthUnit; got " + repr_of(value));
    return nb::cast<bitvf::LengthUnit>(value);
  };

  const bitvf::Smoothing smoothingDefaults;

  nb::class_<bitvf::Smoothing>(m, "Smoothing",
      "Gaussian smoothing of the divergence field before the watershed, as bitvf::Smoothing.")
    .def("__init__", [to_unit] (bitvf::Smoothing* self, nb::handle sigma, nb::handle unit) {
          bitvf::Smoothing s;
          s.sigma = to_double(sigma, "sigma");
          s.unit = to_unit(unit);
          new (self) bitvf::Smoothing(s);
        }, nb::kw_only(), "sigma"_a.none() = smoothingDefaults.sigma, "unit"_a.none() = smoothingDefaults.unit,
        "Every field a keyword, with the defaults of bitvf::Smoothing.")
    .def_prop_rw("sigma",
        [] (const bitvf::Smoothing& s) { return s.sigma; },
        [] (bitvf::Smoothing& s, nb::handle v) { s.sigma = to_double(v, "sigma"); },
        nb::for_setter(nb::arg("value").none()),
        "Width of the Gaussian, in unit; 0 disables the smoothing. 1 mps is a reasonable starting value.")
    .def_prop_rw("unit",
        [] (const bitvf::Smoothing& s) { return s.unit; },
        [to_unit] (bitvf::Smoothing& s, nb::handle v) { s.unit = to_unit(v); },
        nb::for_setter(nb::arg("value").none()), "Unit of sigma.")
    .def("__repr__", [] (const bitvf::Smoothing& s) { return nb::str(smoothing_repr(s).c_str()); })
    .def("__eq__", [] (const bitvf::Smoothing& a, nb::handle other) -> nb::object {
        if (!nb::isinstance<bitvf::Smoothing>(other)) return not_implemented();
        return nb::bool_(same_smoothing(a, nb::cast<const bitvf::Smoothing&>(other)));
      }, "other"_a.none())
    .attr("__hash__") = nb::none();

  const bitvf::Config configDefaults;

  nb::class_<bitvf::Config>(m, "Config",
      "Free parameters of the void finder, as bitvf::Config. Its fields are checked by the call that "
      "reads them, not here.")
    .def("__init__",
        [] (bitvf::Config* self, nb::handle smoothing, nb::handle divergenceCellSize, nb::handle keepDivergence,
            nb::handle keepAssignment, nb::handle nRealizationsMean, nb::handle verbosity) {
          bitvf::Config c;
          c.smoothing = instance<bitvf::Smoothing>(smoothing, "smoothing", "a bitvf.Smoothing");
          c.divergenceCellSize = to_optional_double(divergenceCellSize, "divergence_cell_size");
          c.keepDivergence = to_bool(keepDivergence, "keep_divergence");
          c.keepAssignment = to_bool(keepAssignment, "keep_assignment");
          c.nRealizationsMean = to_unsigned(nRealizationsMean, "n_realizations_mean");
          c.verbosity = to_verbosity(verbosity);
          new (self) bitvf::Config(c);
        },
        nb::kw_only(), "smoothing"_a.none() = configDefaults.smoothing,
        "divergence_cell_size"_a.none() = nb::none(),
        "keep_divergence"_a.none() = configDefaults.keepDivergence,
        "keep_assignment"_a.none() = configDefaults.keepAssignment,
        "n_realizations_mean"_a.none() = configDefaults.nRealizationsMean,
        "verbosity"_a.none() = verbosity_name(configDefaults.verbosity),
        "Every field a keyword, with the defaults of bitvf::Config.")
    .def_prop_rw("smoothing",
        [] (bitvf::Config& c) -> bitvf::Smoothing& { return c.smoothing; },
        [] (bitvf::Config& c, nb::handle v) {
          c.smoothing = instance<bitvf::Smoothing>(v, "smoothing", "a bitvf.Smoothing");
        },
        nb::for_getter(nb::rv_policy::reference_internal), nb::for_setter(nb::arg("value").none()),
        "Smoothing of the divergence field: a reference into this Config, so that setting its fields "
        "changes it; assigning one copies it in.")
    .def_prop_rw("divergence_cell_size",
        [] (const bitvf::Config& c) -> nb::object {
          if (!c.divergenceCellSize) return nb::none();
          return nb::float_(*c.divergenceCellSize);
        },
        [] (bitvf::Config& c, nb::handle v) {
          c.divergenceCellSize = to_optional_double(v, "divergence_cell_size");
        },
        nb::for_setter(nb::arg("value").none()),
        "Side of the divergence-grid cell, in units of the mean particle separation; None for "
        "min(1, (20 / n_realizations)^(1/3)).")
    .def_prop_rw("keep_divergence",
        [] (const bitvf::Config& c) { return c.keepDivergence; },
        [] (bitvf::Config& c, nb::handle v) { c.keepDivergence = to_bool(v, "keep_divergence"); },
        nb::for_setter(nb::arg("value").none()), "Fill Result.divergence with the grid and the field.")
    .def_prop_rw("keep_assignment",
        [] (const bitvf::Config& c) { return c.keepAssignment; },
        [] (bitvf::Config& c, nb::handle v) { c.keepAssignment = to_bool(v, "keep_assignment"); },
        nb::for_setter(nb::arg("value").none()), "Fill Result.assignment, where each tracer ended up.")
    .def_prop_rw("n_realizations_mean",
        [] (const bitvf::Config& c) { return c.nRealizationsMean; },
        [] (bitvf::Config& c, nb::handle v) { c.nRealizationsMean = to_unsigned(v, "n_realizations_mean"); },
        nb::for_setter(nb::arg("value").none()),
        "Realizations averaged into each displacement given; a record written as NRECMEAN.")
    .def_prop_rw("verbosity",
        [] (const bitvf::Config& c) { return verbosity_name(c.verbosity); },
        [] (bitvf::Config& c, nb::handle v) { c.verbosity = to_verbosity(v); },
        nb::for_setter(nb::arg("value").none()), "What the call prints.")
    .def("__repr__", [] (const bitvf::Config& c) { return nb::str(config_repr(c).c_str()); })
    .def("__eq__", [] (const bitvf::Config& a, nb::handle other) -> nb::object {
        if (!nb::isinstance<bitvf::Config>(other)) return not_implemented();
        return nb::bool_(same_config(a, nb::cast<const bitvf::Config&>(other)));
      }, "other"_a.none())
    .attr("__hash__") = nb::none();

  nb::class_<bitvf::BoxRsd>(m, "BoxRsd",
      "Correction of redshift-space distortions in a box, as bitvf::BoxRsd. axis, redshift, bias "
      "and sigma have no default.")
    .def("__init__",
        [] (bitvf::BoxRsd* self, nb::handle axis, nb::handle redshift, nb::handle bias, nb::handle sigma,
            nb::handle correction, nb::handle zspaceVoids) {
          bitvf::BoxRsd r;
          if (!axis.is_none()) r.axis = to_unsigned(axis, "axis");
          r.redshift = to_double(redshift, "redshift");
          r.bias = to_double(bias, "bias");
          r.sigma = to_double(sigma, "sigma");
          r.correction = instance<otswap::CorrectionConfig>(correction, "correction",
                                                            "an otswap.CorrectionConfig");
          r.zspaceVoids = to_bool(zspaceVoids, "zspace_voids");
          new (self) bitvf::BoxRsd(r);
        },
        nb::kw_only(), "axis"_a.none(), "redshift"_a.none(), "bias"_a.none(), "sigma"_a.none(),
        "correction"_a.none() = otswap::CorrectionConfig(), "zspace_voids"_a.none() = false,
        "axis, redshift, bias and sigma required; axis None is refused at the call.")
    .def_prop_rw("axis",
        [] (const bitvf::BoxRsd& r) -> nb::object {
          if (!r.axis) return nb::none();
          return nb::int_(*r.axis);
        },
        [] (bitvf::BoxRsd& r, nb::handle v) {
          if (v.is_none()) r.axis.reset();
          else r.axis = to_unsigned(v, "axis");
        },
        nb::for_setter(nb::arg("value").none()), "The line of sight: 0, 1 or 2, for x, y or z.")
    .def_prop_rw("redshift",
        [] (const bitvf::BoxRsd& r) { return r.redshift; },
        [] (bitvf::BoxRsd& r, nb::handle v) { r.redshift = to_double(v, "redshift"); },
        nb::for_setter(nb::arg("value").none()), "The redshift of the box, where f is read.")
    .def_prop_rw("bias",
        [] (const bitvf::BoxRsd& r) { return r.bias; },
        [] (bitvf::BoxRsd& r, nb::handle v) { r.bias = to_double(v, "bias"); },
        nb::for_setter(nb::arg("value").none()), "The linear bias of the tracers.")
    .def_prop_rw("sigma",
        [] (const bitvf::BoxRsd& r) { return r.sigma; },
        [] (bitvf::BoxRsd& r, nb::handle v) { r.sigma = to_double(v, "sigma"); },
        nb::for_setter(nb::arg("value").none()),
        "Width of the neighbour average, in coordinate units; 0 for none. For coordinates in Mpc/h, "
        "10 is a reasonable starting value, to check in each analysis.")
    .def_prop_rw("correction",
        [] (bitvf::BoxRsd& r) -> otswap::CorrectionConfig& { return r.correction; },
        [] (bitvf::BoxRsd& r, nb::handle v) {
          r.correction = instance<otswap::CorrectionConfig>(v, "correction", "an otswap.CorrectionConfig");
        },
        nb::for_getter(nb::rv_policy::reference_internal), nb::for_setter(nb::arg("value").none()),
        "The options of otswap's correction: a reference into this BoxRsd; assigning one copies it in.")
    .def_prop_rw("zspace_voids",
        [] (const bitvf::BoxRsd& r) { return r.zspaceVoids; },
        [] (bitvf::BoxRsd& r, nb::handle v) { r.zspaceVoids = to_bool(v, "zspace_voids"); },
        nb::for_setter(nb::arg("value").none()), "Also search voids on the redshift-space reconstruction.")
    .def("__repr__", [] (const bitvf::BoxRsd& r) { return nb::str(rsd_repr(r).c_str()); })
    .def("__eq__", [] (const bitvf::BoxRsd& a, nb::handle other) -> nb::object {
        if (!nb::isinstance<bitvf::BoxRsd>(other)) return not_implemented();
        return nb::bool_(same_rsd(a, nb::cast<const bitvf::BoxRsd&>(other)));
      }, "other"_a.none())
    .attr("__hash__") = nb::none();

  auto to_rsd = [] (nb::handle value) -> std::optional<bitvf::BoxRsd> {
    if (value.is_none()) return std::nullopt;
    return instance<bitvf::BoxRsd>(value, "rsd", "a bitvf.BoxRsd or None");
  };

  nb::class_<bitvf::BoxConfig>(m, "BoxConfig", "The parameters of run_box, as bitvf::BoxConfig.")
    .def("__init__",
        [to_rsd] (bitvf::BoxConfig* self, nb::handle ot, nb::handle limits, nb::handle mps, nb::handle rsd,
                  nb::handle voids, nb::handle keepReconstructions) {
          bitvf::BoxConfig c;
          c.ot = instance<otswap::Config>(ot, "ot", "an otswap.Config");
          c.limits = to_box(limits, "limits");
          c.mps = to_optional_double(mps, "mps");
          c.rsd = to_rsd(rsd);
          c.voids = instance<bitvf::Config>(voids, "voids", "a bitvf.Config");
          c.keepReconstructions = to_bool(keepReconstructions, "keep_reconstructions");
          new (self) bitvf::BoxConfig(c);
        },
        nb::kw_only(), "ot"_a.none() = otswap::Config(), "limits"_a.none() = nb::none(),
        "mps"_a.none() = nb::none(), "rsd"_a.none() = nb::none(), "voids"_a.none() = bitvf::Config(),
        "keep_reconstructions"_a.none() = false,
        "Every field a keyword, with the defaults of bitvf::BoxConfig.")
    .def_prop_rw("ot",
        [] (bitvf::BoxConfig& c) -> otswap::Config& { return c.ot; },
        [] (bitvf::BoxConfig& c, nb::handle v) { c.ot = instance<otswap::Config>(v, "ot", "an otswap.Config"); },
        nb::for_getter(nb::rv_policy::reference_internal), nb::for_setter(nb::arg("value").none()),
        "The reconstruction, for both reconstructions of a run with correction: a reference into this "
        "BoxConfig; assigning one copies it in.")
    .def_prop_rw("limits",
        [] (const bitvf::BoxConfig& c) { return box_or_none(c.limits); },
        [] (bitvf::BoxConfig& c, nb::handle v) { c.limits = to_box(v, "limits"); },
        nb::for_setter(nb::arg("value").none()),
        "The box of the run, an otswap.BoxLimits (a pair of corners is accepted), or None for the "
        "bounding box of the tracers, and of the randoms when given.")
    .def_prop_rw("mps",
        [] (const bitvf::BoxConfig& c) -> nb::object {
          if (!c.mps) return nb::none();
          return nb::float_(*c.mps);
        },
        [] (bitvf::BoxConfig& c, nb::handle v) { c.mps = to_optional_double(v, "mps"); },
        nb::for_setter(nb::arg("value").none()),
        "The mean particle separation of both reconstructions, or None for otswap's.")
    .def_prop_rw("rsd",
        [] (bitvf::BoxConfig& c) -> bitvf::BoxRsd* { return c.rsd ? &*c.rsd : nullptr; },
        [to_rsd] (bitvf::BoxConfig& c, nb::handle v) { c.rsd = to_rsd(v); },
        nb::for_getter(nb::rv_policy::reference_internal), nb::for_setter(nb::arg("value").none()),
        "The correction of redshift-space distortions, or None: a reference into this BoxConfig; "
        "assigning one copies it in.")
    .def_prop_rw("voids",
        [] (bitvf::BoxConfig& c) -> bitvf::Config& { return c.voids; },
        [] (bitvf::BoxConfig& c, nb::handle v) { c.voids = instance<bitvf::Config>(v, "voids", "a bitvf.Config"); },
        nb::for_getter(nb::rv_policy::reference_internal), nb::for_setter(nb::arg("value").none()),
        "The void search, on every pass: a reference into this BoxConfig; assigning one copies it in.")
    .def_prop_rw("keep_reconstructions",
        [] (const bitvf::BoxConfig& c) { return c.keepReconstructions; },
        [] (bitvf::BoxConfig& c, nb::handle v) { c.keepReconstructions = to_bool(v, "keep_reconstructions"); },
        nb::for_setter(nb::arg("value").none()),
        "Keep the otswap results and the corrected catalogue in BitVFResult.")
    .def("__repr__", [] (const bitvf::BoxConfig& c) {
        const std::string text = "BoxConfig(ot=" + repr_of(nb::cast(c.ot)) + ", limits=" +
          repr_of(box_or_none(c.limits)) + ", mps=" + optional_repr(c.mps).c_str() + ", rsd=" +
          (c.rsd ? rsd_repr(*c.rsd) : std::string("None")) + ", voids=" + config_repr(c.voids) +
          ", keep_reconstructions=" + (c.keepReconstructions ? "True" : "False") + ")";
        return nb::str(text.c_str());
      })
    .def("__eq__", [] (const bitvf::BoxConfig& a, nb::handle other) -> nb::object {
        if (!nb::isinstance<bitvf::BoxConfig>(other)) return not_implemented();
        const bitvf::BoxConfig& b = nb::cast<const bitvf::BoxConfig&>(other);
        const bool rsd = (!a.rsd && !b.rsd) || (a.rsd && b.rsd && same_rsd(*a.rsd, *b.rsd));
        return nb::bool_(same_otswap(a.ot, b.ot) && same_box(a.limits, b.limits) && a.mps == b.mps && rsd &&
                         same_config(a.voids, b.voids) && a.keepReconstructions == b.keepReconstructions);
      }, "other"_a.none())
    .attr("__hash__") = nb::none();

  // ==========================================================================
  // Results
  // ==========================================================================

  nb::class_<bitvf::BasinCatalog>(m, "BasinCatalog",
      "One catalogue of watershed basins, the voids or the overdense basins; row b is basin b, in the "
      "order of its seed cell. Read-only.")
    .def_prop_ro("n_basins", [] (const bitvf::BasinCatalog& c) { return c.nBasins; })
    .def_prop_ro("centre", [] (CatalogHandle self) { return triplets(self, catalog_of(self).centre); },
                 "Cartesian centre, in coordinate units, shape (n_basins, 3).")
    .def_prop_ro("centre_sky", [] (CatalogHandle self) { return triplets(self, catalog_of(self).centreSky); },
                 "Sky coordinates of the centre; shape (0, 3) in a box.")
    .def_prop_ro("radius", [] (CatalogHandle self) { return rows<double>(self, catalog_of(self).radius); },
                 "Radius of the sphere with the volume of the core, in coordinate units.")
    .def_prop_ro("radius_ext", [] (CatalogHandle self) { return rows<double>(self, catalog_of(self).radiusExt); },
                 "Radius of the sphere with the volume of the extended basin, in coordinate units.")
    .def_prop_ro("delta_av", [] (CatalogHandle self) { return rows<double>(self, catalog_of(self).deltaAv); },
                 "Mean density contrast of the tracers in the core; (0,) after a field search.")
    .def_prop_ro("delta_2mps", [] (CatalogHandle self) { return rows<double>(self, catalog_of(self).delta2Mps); },
                 "Density contrast within 2 mps of the centre; (0,) after a field search.")
    .def_prop_ro("edge_flag", [] (CatalogHandle self) { return rows<std::uint8_t>(self, catalog_of(self).edgeFlag); },
                 "Contact with the boundary: 0 none, 1 only the extended basin, 2 the core.")
    .def_prop_ro("contamination",
                 [] (CatalogHandle self) { return rows<double>(self, catalog_of(self).contamination); },
                 "Fraction of the basin's surface on the boundary.")
    .def_prop_ro("n_tracers", [] (CatalogHandle self) { return rows<std::uint64_t>(self, catalog_of(self).nTracers); },
                 "Tracers starting in the cells of the core; (0,) after a field search.")
    .def_prop_ro("field_centre",
                 [] (CatalogHandle self) { return rows<double>(self, catalog_of(self).fieldCentre); },
                 "The field at the centre, after a field search; (0,) otherwise.");

  nb::class_<bitvf::DivergenceGrid>(m, "DivergenceGrid",
      "The grid of the void finder and the divergence field on it. Read-only: it comes from a Result "
      "or from io.read_divergence, and can be passed back as grid. Every per-cell array of a box grid "
      "has shape (nx, ny, nz), n_cells_axis, C order: a[i, j, k] is the cell (i, j, k); (0,) when empty.")
    .def_prop_ro("geometry", [] (const bitvf::DivergenceGrid& g) { return geometry_name(g.geometry); },
                 "\"box\" or \"lightcone\".")
    .def_prop_ro("angle_unit", [] (const bitvf::DivergenceGrid& g) -> nb::object {
        if (g.geometry == otswap::Geometry::Box) return nb::none();
        return nb::str("rad");
      }, "The unit of the angle fields; None for a box.")
    .def_prop_ro("n_cells", [] (const bitvf::DivergenceGrid& g) { return g.nCells; })
    .def_prop_ro("cell_size", [] (const bitvf::DivergenceGrid& g) { return g.cellSize; },
                 "Cell side, in coordinate units; NaN in a lightcone.")
    .def_prop_ro("lower", [] (const bitvf::DivergenceGrid& g) { return triple(g.lower); },
                 "x, y, z of the lower grid face, in coordinate units.")
    .def_prop_ro("upper", [] (const bitvf::DivergenceGrid& g) { return triple(g.upper); },
                 "x, y, z of the upper face of the padded box as placed, in coordinate units.")
    .def_prop_ro("n_cells_axis", [] (const bitvf::DivergenceGrid& g) { return triple(g.nCellsAxis); },
                 "Cells along x, y, z.")
    .def_prop_ro("limits", [] (const bitvf::DivergenceGrid& g) { return box_or_none(g.limits); },
                 "The box limits of the search that filled the grid, or None.")
    .def_prop_ro("n_realizations", [] (const bitvf::DivergenceGrid& g) { return g.nRealizations; })
    .def_prop_ro("n_realizations_mean", [] (const bitvf::DivergenceGrid& g) { return g.nRealizationsMean; })
    .def_prop_ro("distance_min", [] (const bitvf::DivergenceGrid& g) { return g.distanceMin; })
    .def_prop_ro("distance_max", [] (const bitvf::DivergenceGrid& g) { return g.distanceMax; })
    .def_prop_ro("ra_min", [] (const bitvf::DivergenceGrid& g) { return g.raMin; })
    .def_prop_ro("ra_max", [] (const bitvf::DivergenceGrid& g) { return g.raMax; })
    .def_prop_ro("dec_min", [] (const bitvf::DivergenceGrid& g) { return g.decMin; })
    .def_prop_ro("dec_max", [] (const bitvf::DivergenceGrid& g) { return g.decMax; })
    .def_prop_ro("cell_size_mps", [] (const bitvf::DivergenceGrid& g) { return g.cellSizeMps; },
                 "The divergence cell size used, in units of the mps.")
    .def_prop_ro("mps", [] (const bitvf::DivergenceGrid& g) { return g.mps; },
                 "Box: the mean particle separation the grid was built with, in coordinate units.")
    .def_prop_ro("mps_profile", [] (const bitvf::DivergenceGrid& g) -> const otswap::MpsProfile* {
        return g.mpsProfile.redshift.empty() ? nullptr : &g.mpsProfile;
      }, nb::rv_policy::reference_internal, "Lightcone: the mps(z) profile; None for a box.")
    .def_prop_ro("cell_ra_min", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.cellRaMin, g.nCellsAxis); })
    .def_prop_ro("cell_ra_max", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.cellRaMax, g.nCellsAxis); })
    .def_prop_ro("cell_dec_min", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.cellDecMin, g.nCellsAxis); })
    .def_prop_ro("cell_dec_max", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.cellDecMax, g.nCellsAxis); })
    .def_prop_ro("cell_distance_min", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.cellDistanceMin, g.nCellsAxis); })
    .def_prop_ro("cell_distance_max", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.cellDistanceMax, g.nCellsAxis); })
    .def_prop_ro("centre", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.centre, g.nCellsAxis, 3);
      }, "Cartesian centre of each cell, in coordinate units, shape (nx, ny, nz, 3).")
    .def_prop_ro("centre_sky", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.centreSky, g.nCellsAxis, 3);
      }, "Sky coordinates of each cell centre; shape (0, 3) in a box.")
    .def_prop_ro("divergence", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.divergence, g.nCellsAxis);
      }, "The divergence of each cell, before smoothing.")
    .def_prop_ro("divergence_smoothed", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.divergenceSmoothed, g.nCellsAxis);
      }, "The divergence after smoothing, the field the watershed splits, with sigma > 0; (0,) otherwise.")
    .def_prop_ro("volume", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.volume, g.nCellsAxis);
      }, "Volume of each cell, in coordinate units cubed.")
    .def_prop_ro("weight", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.weight, g.nCellsAxis);
      }, "Tracer weight of each cell, the weight of the smoothing.")
    .def_prop_ro("weight_smoothed", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<double>(self, g.weightSmoothed, g.nCellsAxis);
      }, "Smoothed tracer weight of each cell, with sigma > 0; (0,) otherwise.")
    .def_prop_ro("n_displacements", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<std::uint64_t>(self, g.nDisplacements, g.nCellsAxis);
      }, "Displacement ends counted in each cell, before the fill.")
    .def_prop_ro("n_tracers", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<std::uint64_t>(self, g.nTracers, g.nCellsAxis);
      }, "Tracers starting in each cell.")
    .def_prop_ro("void_index", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<std::int64_t>(self, g.voidIndex, g.nCellsAxis);
      }, "Row of the cell's void in Result.voids; -1 none.")
    .def_prop_ro("overdense_index", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<std::int64_t>(self, g.overdenseIndex, g.nCellsAxis);
      }, "Row of the cell's overdense basin in Result.overdense; -1 none.")
    .def_prop_ro("masked", [] (GridHandle self) {
        const auto& g = grid_of(self); return cells<bool>(self, g.masked, g.nCellsAxis);
      }, "Lightcone: whether the cell centre is unobserved; (0,) in a box.");

  nb::class_<bitvf::TracerAssignment>(m, "TracerAssignment",
      "Where each input tracer ended up, one row per input row. Read-only.")
    .def_prop_ro("cell_index", [] (AssignmentHandle self) {
        return rows<std::int64_t>(self, assignment_of(self).cellIndex);
      }, "Row of the cell the tracer starts in, as a flat row of the grid (numpy.unravel_index with "
         "n_cells_axis gives (i, j, k)); -1 when outside the grid or left out.")
    .def_prop_ro("void_index", [] (AssignmentHandle self) {
        return rows<std::int64_t>(self, assignment_of(self).voidIndex);
      }, "Row of the tracer's void; -1 none.")
    .def_prop_ro("overdense_index", [] (AssignmentHandle self) {
        return rows<std::int64_t>(self, assignment_of(self).overdenseIndex);
      }, "Row of the tracer's overdense basin; -1 none.")
    .def_prop_ro("distance", [] (AssignmentHandle self) {
        return rows<double>(self, assignment_of(self).distance);
      }, "Distance from the centre of the tracer's basin, in coordinate units; NaN when in none.")
    .def_prop_ro("in_mask", [] (AssignmentHandle self) {
        return rows<bool>(self, assignment_of(self).inMask);
      }, "False when the tracer falls on an unobserved pixel of the mask; (0,) in a box.")
    .def_prop_ro("in_cut", [] (AssignmentHandle self) {
        return rows<bool>(self, assignment_of(self).inCut);
      }, "False when the tracer lies outside the redshift cut, and so was skipped; (0,) in a box.")
    .def_prop_ro("in_box", [] (AssignmentHandle self) {
        return rows<bool>(self, assignment_of(self).inBox);
      }, "False when the tracer lies outside the box limits, and so was skipped; (0,) in a lightcone.")
    .def_prop_ro("used", [] (AssignmentHandle self) {
        return rows<bool>(self, assignment_of(self).used);
      }, "True when a valid displacement of the tracer entered the divergence; (0,) after "
         "find_voids_box_from_divergence.");

  nb::class_<bitvf::Counts>(m, "Counts", "What the void finder counted. Read-only.")
    .def_ro("tracers", &bitvf::Counts::tracers)
    .def_ro("tracers_used", &bitvf::Counts::tracersUsed)
    .def_ro("tracers_uncorrected", &bitvf::Counts::tracersUncorrected)
    .def_ro("tracers_outside_box", &bitvf::Counts::tracersOutsideBox)
    .def_ro("tracers_outside_grid", &bitvf::Counts::tracersOutsideGrid)
    .def_ro("tracers_outside_cut", &bitvf::Counts::tracersOutsideCut)
    .def_ro("randoms", &bitvf::Counts::randoms)
    .def_ro("tracers_in_masked_cells", &bitvf::Counts::tracersInMaskedCells,
            "Lightcone: tracers on an observed pixel whose cell is masked; they keep their cell but are in no "
            "basin. 0 in a box.")
    .def_ro("randoms_outside_grid", &bitvf::Counts::randomsOutsideGrid)
    .def_ro("displacements", &bitvf::Counts::displacements)
    .def_ro("displacements_rejected", &bitvf::Counts::displacementsRejected);

  nb::class_<bitvf::Result>(m, "Result", "Output of one void search. Read-only.")
    .def_prop_ro("geometry", [] (const bitvf::Result& r) { return geometry_name(r.geometry); },
                 "\"box\" or \"lightcone\".")
    .def_prop_ro("angle_unit", [] (const bitvf::Result& r) -> nb::object {
        if (r.geometry == otswap::Geometry::Box) return nb::none();
        return nb::str("rad");
      }, "The unit of the sky arrays; None for a box.")
    .def_prop_ro("config", [] (const bitvf::Result& r) {
        bitvf::Config c = r.config;
        c.log = &std::clog;
        return c;
      }, "The configuration the search ran with, divergence_cell_size the value used: a new copy on "
         "each access, so changing it changes nothing in the result.")
    .def_prop_ro("n_tracers", [] (const bitvf::Result& r) { return r.nTracers; },
                 "Tracers given; 0 for a field result.")
    .def_prop_ro("n_realizations", [] (const bitvf::Result& r) { return r.nRealizations; })
    .def_prop_ro("ids", [] (const bitvf::Result& r) { return ids_array(r.ids, r.nTracers); },
                 nb::rv_policy::reference,
                 "The ID of each tracer, shape (n_tracers,), a new read-only copy on each access; "
                 "0..n_tracers-1 when none were set.")
    .def_prop_ro("random_seed", [] (const bitvf::Result& r) -> nb::object {
        if (!r.randomSeed) return nb::none();
        return nb::int_(*r.randomSeed);
      }, "The seed of the randoms bitvf drew in a box pipeline; None when it drew none.")
    .def_prop_ro("voids", [] (bitvf::Result& r) -> bitvf::BasinCatalog& { return r.voids; },
                 nb::rv_policy::reference_internal, "The voids.")
    .def_prop_ro("overdense", [] (bitvf::Result& r) -> bitvf::BasinCatalog& { return r.overdense; },
                 nb::rv_policy::reference_internal, "The overdense basins.")
    .def_prop_ro("divergence", [] (bitvf::Result& r) -> bitvf::DivergenceGrid* {
        return r.divergence ? &*r.divergence : nullptr;
      }, nb::rv_policy::reference_internal,
      "The grid and the divergence field, with Config.keep_divergence; None otherwise.")
    .def_prop_ro("assignment", [] (bitvf::Result& r) -> bitvf::TracerAssignment* {
        return r.assignment ? &*r.assignment : nullptr;
      }, nb::rv_policy::reference_internal,
      "Where each tracer ended up, with Config.keep_assignment; None otherwise.")
    .def_prop_ro("uncorrected", [] (ResultHandle self) { return rows<bool>(self, result_of(self).uncorrected); },
                 "Real-space pass of a pipeline with correction: True for a tracer the correction left "
                 "unmoved, excluded; (0,) elsewhere, which reads as all False.")
    .def_prop_ro("outside_grid", [] (ResultHandle self) { return rows<bool>(self, result_of(self).outsideGrid); },
                 "True for a tracer outside the reused grid, excluded; (0,) when no grid was reused.")
    .def_prop_ro("outside_cut", [] (ResultHandle self) { return rows<bool>(self, result_of(self).outsideCut); },
                 "Lightcone pipeline with correction: True for a tracer moved outside the redshift cut; "
                 "(0,) elsewhere.")
    .def_prop_ro("limits", [] (const bitvf::Result& r) { return box_or_none(r.limits); },
                 "Box: the box limits of the search, or None.")
    .def_prop_ro("mps", [] (const bitvf::Result& r) { return r.mps; },
                 "Box: the mean particle separation the search used, in coordinate units; NaN in a "
                 "lightcone and in a field result given none.")
    .def_prop_ro("mps_profile", [] (const bitvf::Result& r) -> const otswap::MpsProfile* {
        return r.mpsProfile.redshift.empty() ? nullptr : &r.mpsProfile;
      }, nb::rv_policy::reference_internal, "Lightcone: the mps(z) profile; None for a box.")
    .def_prop_ro("mask_nside", [] (const bitvf::Result& r) { return r.maskNside; },
                 "Lightcone: the NSIDE of the mask the search was given; 0 in a box.")
    .def_prop_ro("redshift_cut", [] (const bitvf::Result& r) -> nb::object {
        if (r.geometry == otswap::Geometry::Box) return nb::none();
        return nb::cast(r.redshiftCut);
      }, "Lightcone: the closed range of observed redshift the search kept, an otswap.RedshiftCut; None for "
         "a box.")
    .def_prop_ro("distances", [] (const bitvf::Result& r) -> nb::object {
        if (!r.distances) return nb::none();
        return nb::cast(otswap::DistanceTable(*r.distances));
      }, "Lightcone: the distance table of the search, a copy; None for a box.")
    .def_prop_ro("counts", [] (bitvf::Result& r) -> bitvf::Counts& { return r.counts; },
                 nb::rv_policy::reference_internal, "What the search counted.")
    .def_prop_ro("elapsed_seconds", [] (const bitvf::Result& r) { return r.elapsedSeconds; },
                 "Wall time of the call, in seconds.");

  nb::class_<PartitionObject>(m, "Partition",
      "The basins of a field, as watershed_box finds them, with no catalogue. Read-only.")
    .def_prop_ro("n_voids", [] (const PartitionObject& p) { return p.partition.nVoids; })
    .def_prop_ro("n_overdense", [] (const PartitionObject& p) { return p.partition.nOverdense; })
    .def_prop_ro("void_index", [] (PartitionHandle self) {
        const auto& p = partition_of(self); return cells<std::int64_t>(self, p.partition.voidIndex, p.n);
      }, "Row of each cell's void, -1 none; the shape of the field.")
    .def_prop_ro("overdense_index", [] (PartitionHandle self) {
        const auto& p = partition_of(self); return cells<std::int64_t>(self, p.partition.overdenseIndex, p.n);
      }, "Row of each cell's overdense basin, -1 none; the shape of the field.")
    .def_prop_ro("void_seed", [] (PartitionHandle self) {
        return rows<std::int64_t>(self, partition_of(self).partition.voidSeed);
      }, "Flat row of each void's seed cell, ascending (numpy.unravel_index gives (i, j, k)).")
    .def_prop_ro("overdense_seed", [] (PartitionHandle self) {
        return rows<std::int64_t>(self, partition_of(self).partition.overdenseSeed);
      }, "Flat row of each overdense basin's seed cell, ascending.")
    .def_prop_ro("void_cells", [] (PartitionHandle self) {
        return rows<std::uint64_t>(self, partition_of(self).partition.voidCells);
      }, "Cells of each extended void.")
    .def_prop_ro("void_core_cells", [] (PartitionHandle self) {
        return rows<std::uint64_t>(self, partition_of(self).partition.voidCoreCells);
      }, "Cells of each void's core, where the field is below 0.")
    .def_prop_ro("overdense_cells", [] (PartitionHandle self) {
        return rows<std::uint64_t>(self, partition_of(self).partition.overdenseCells);
      }, "Cells of each overdense basin.")
    .def_prop_ro("overdense_core_cells", [] (PartitionHandle self) {
        return rows<std::uint64_t>(self, partition_of(self).partition.overdenseCoreCells);
      }, "Cells of each overdense basin's core, where the field is at or above 0.");

  nb::class_<bitvf::BitVFResult>(m, "BitVFResult",
      "Output of run_box. Read-only. Its otswap results are otswap's own types.")
    .def_prop_ro("voids", [] (bitvf::BitVFResult& r) -> bitvf::Result& { return r.voids; },
                 nb::rv_policy::reference_internal,
                 "The voids of the run, one row per input tracer, in input order.")
    .def_prop_ro("zspace_voids", [] (bitvf::BitVFResult& r) -> bitvf::Result* {
        return r.zspaceVoids ? &*r.zspaceVoids : nullptr;
      }, nb::rv_policy::reference_internal,
      "The voids of the redshift-space pass, with BoxRsd.zspace_voids; None otherwise.")
    .def_prop_ro("reconstruction", [] (bitvf::BitVFResult& r) -> otswap::Result* {
        return r.reconstruction ? &*r.reconstruction : nullptr;
      }, nb::rv_policy::reference_internal,
      "With keep_reconstructions: the reconstruction the voids were found on; None otherwise.")
    .def_prop_ro("zspace_reconstruction", [] (bitvf::BitVFResult& r) -> otswap::Result* {
        return r.zspaceReconstruction ? &*r.zspaceReconstruction : nullptr;
      }, nb::rv_policy::reference_internal,
      "With keep_reconstructions and a correction: the redshift-space reconstruction; None otherwise.")
    .def_prop_ro("real_space", [] (bitvf::BitVFResult& r) -> otswap::RealSpaceCatalog* {
        return r.realSpace ? &*r.realSpace : nullptr;
      }, nb::rv_policy::reference_internal,
      "With keep_reconstructions and a correction: the corrected catalogue; None otherwise.");

  // ==========================================================================
  // Void finder
  // ==========================================================================

  m.def("find_voids_box",
      [] (nb::handle first, nb::handle displacement, nb::handle valid, nb::handle mps, nb::handle box,
          nb::handle configArg, nb::handle ids, nb::handle grid) {
        bitvf::Config config = instance<bitvf::Config>(configArg, "config", "a bitvf.Config");
        const bitvf::DivergenceGrid* reused = grid_or_null(grid);
        LiveReport report;
        config.log = report.log();
        if (nb::isinstance<otswap::Result>(first)) {
          std::vector<std::string> given;
          if (!displacement.is_none()) given.push_back("displacement");
          if (!valid.is_none()) given.push_back("valid");
          if (!mps.is_none()) given.push_back("mps");
          if (!box.is_none()) given.push_back("box");
          if (!ids.is_none()) given.push_back("ids");
          if (!given.empty())
            throw bitvf::Error("find_voids_box: " + listed(given) + (given.size() == 1 ? " does" : " do") +
                               " not apply to a reconstruction (otswap.Result)");
          const otswap::Result& reconstruction = nb::cast<const otswap::Result&>(first);
          return result_object(report.run([&] { return bitvf::findVoidsBox(reconstruction, config, reused); }));
        }
        if (displacement.is_none())
          throw bitvf::Error("find_voids_box: displacement is required with arrays");
        if (mps.is_none())
          throw bitvf::Error("find_voids_box: mps is required with arrays");
        const std::vector<double> tracers = rows3(first, "tracers");
        const std::size_t n = tracers.size() / 3;
        const CArray d = as_float64(displacement, "displacement");
        if (d.ndim() != 3 || d.shape(2) != 3 || d.shape(0) == 0 || d.shape(1) == 0)
          throw bitvf::Error("displacement must have shape (n_realizations, n_tracers, 3), neither 0; it has "
                             "shape " + shape_of(d));
        if (d.shape(1) != n)
          throw bitvf::Error("displacement must have shape (n_realizations, " + std::to_string(n) + ", 3), " +
                             std::to_string(n) + " the rows of tracers; it has shape " + shape_of(d));
        const std::vector<double> displacements(d.data(), d.data() + d.size());
        std::vector<std::uint8_t> validity;
        if (!valid.is_none()) {
          const CArray v = as_float64(valid, "valid");
          if (v.size() == 0)
            throw bitvf::Error("valid is empty; omit it to take every displacement as valid");
          validity = flags(v, "valid", {d.shape(0), n}, "as displacement");
        }
        const double separation = to_double(mps, "mps");
        const std::optional<otswap::BoxLimits> limits = to_box(box, "box");
        const std::vector<std::int64_t> identity = ids_value(ids, n);
        return result_object(report.run([&] {
          return bitvf::findVoidsBox(tracers, displacements, validity, separation, config, limits, identity, reused);
        }));
      },
      "tracers"_a.none(), "displacement"_a.none() = nb::none(), "valid"_a.none() = nb::none(), nb::kw_only(),
      "mps"_a.none() = nb::none(), "box"_a.none() = nb::none(), "config"_a.none() = bitvf::Config(),
      "ids"_a.none() = nb::none(), "grid"_a.none() = nb::none(),
      "Find voids and overdense basins in a box, from arrays or from an otswap.Result.");

  m.def("find_voids_box_from_divergence",
      [] (nb::handle grid, nb::handle tracers, nb::handle configArg, nb::handle ids) {
        const bitvf::DivergenceGrid& g = instance<bitvf::DivergenceGrid>(grid, "grid", "a bitvf.DivergenceGrid");
        const std::vector<double> t = rows3(tracers, "tracers");
        bitvf::Config config = instance<bitvf::Config>(configArg, "config", "a bitvf.Config");
        const std::vector<std::int64_t> identity = ids_value(ids, t.size() / 3);
        LiveReport report;
        config.log = report.log();
        return result_object(report.run([&] {
          return bitvf::findVoidsBoxFromDivergence(g, t, config, identity);
        }));
      },
      "grid"_a.none(), "tracers"_a.none(), nb::kw_only(), "config"_a.none() = bitvf::Config(),
      "ids"_a.none() = nb::none(),
      "Find voids and overdense basins in a box from a divergence computed earlier, with the tracers "
      "given.");

  m.def("find_voids_box_from_field",
      [] (nb::handle centres, nb::handle values, nb::handle configArg, nb::handle weights, nb::handle mps) {
        const std::vector<double> c = rows3(centres, "centres");
        const std::vector<double> v = column(values, "values");
        bitvf::Config config = instance<bitvf::Config>(configArg, "config", "a bitvf.Config");
        std::vector<double> w;
        if (!weights.is_none()) {
          w = column(weights, "weights");
          if (w.empty()) throw bitvf::Error("weights is empty; omit it for unit weights");
        }
        const std::optional<double> separation = to_optional_double(mps, "mps");
        LiveReport report;
        config.log = report.log();
        return result_object(report.run([&] { return bitvf::findVoidsBoxFromField(c, v, config, w, separation); }));
      },
      "centres"_a.none(), "values"_a.none(), nb::kw_only(), "config"_a.none() = bitvf::Config(),
      "weights"_a.none() = nb::none(), "mps"_a.none() = nb::none(),
      "Find voids and overdense basins in a field given at the centres of the cells of a box lattice.");

  m.def("box_grid",
      [] (nb::handle reconstruction, nb::handle box, nb::handle mps, nb::handle cellSize) {
        if (!reconstruction.is_none()) {
          const otswap::Result& r = instance<otswap::Result>(reconstruction, "reconstruction", "an otswap.Result");
          if (!box.is_none() || !mps.is_none())
            throw bitvf::Error("box_grid: box and mps do not apply to a reconstruction");
          const std::optional<double> size = to_optional_double(cellSize, "cell_size");
          nb::gil_scoped_release released;
          return bitvf::boxGrid(r, size);
        }
        if (box.is_none() || mps.is_none() || cellSize.is_none())
          throw bitvf::Error("box_grid: give a reconstruction, or box, mps and cell_size");
        const otswap::BoxLimits limits = *to_box(box, "box");
        const double separation = to_double(mps, "mps");
        const double size = to_double(cellSize, "cell_size");
        nb::gil_scoped_release released;
        return bitvf::boxGrid(limits, separation, size);
      },
      "reconstruction"_a.none() = nb::none(), nb::kw_only(), "box"_a.none() = nb::none(),
      "mps"_a.none() = nb::none(), "cell_size"_a.none() = nb::none(),
      "The grid find_voids_box builds on box limits, or from a reconstruction, bit for bit.");

  m.def("watershed_box",
      [] (nb::handle field, nb::handle verbosity) {
        const CArray a = as_float64(field, "field");
        if (a.ndim() != 3)
          throw bitvf::Error("field must have shape (nx, ny, nz); it has shape " + shape_of(a));
        const std::vector<double> values(a.data(), a.data() + a.size());
        const std::array<std::size_t, 3> n = {a.shape(0), a.shape(1), a.shape(2)};
        const otswap::Verbosity level = to_verbosity(verbosity);
        LiveReport report;
        PartitionObject p;
        p.partition = report.run([&] { return bitvf::watershedBox(values, n, level, report.log()); });
        p.n = n;
        return p;
      },
      "field"_a.none(), nb::kw_only(), "verbosity"_a.none() = "normal",
      "Split a field on a box grid, shape (nx, ny, nz), into basins by steepest descent.");

  // ==========================================================================
  // Pipeline
  // ==========================================================================

  m.def("run_box",
      [] (nb::handle tracers, nb::handle randoms, nb::handle configArg, nb::handle distances, nb::handle ids) {
        const std::vector<double> t = rows3(tracers, "tracers");
        std::vector<double> r;
        if (!randoms.is_none()) {
          r = rows3(randoms, "randoms");
          if (r.empty()) throw bitvf::Error("randoms is empty; omit it to have bitvf draw the randoms in the box");
        }
        bitvf::BoxConfig config = instance<bitvf::BoxConfig>(configArg, "config", "a bitvf.BoxConfig");
        const otswap::DistanceTable* table = distances.is_none() ? nullptr
          : &instance<otswap::DistanceTable>(distances, "distances", "an otswap.DistanceTable");
        const std::vector<std::int64_t> identity = ids_value(ids, t.size() / 3);
        LiveReport report;
        config.ot.log = report.log();
        config.voids.log = report.log();
        if (config.rsd) config.rsd->correction.log = report.log();
        return pipeline_object(report.run([&] {
          return table ? bitvf::runBox(t, r, *table, config, identity) : bitvf::runBox(t, r, config, identity);
        }));
      },
      "tracers"_a.none(), "randoms"_a.none() = nb::none(), nb::kw_only(),
      "config"_a.none() = bitvf::BoxConfig(), "distances"_a.none() = nb::none(), "ids"_a.none() = nb::none(),
      "Reconstruct a box and find its voids, with the correction of redshift-space distortions when "
      "config.rsd is set.");

  // ==========================================================================
  // bitvf.io
  // ==========================================================================

  nb::module_ io = m.def_submodule("io", "Writing the products of a search, and reading files back.");

  nb::enum_<bitvf::io::Product>(io, "Product", "The products, each with its file-name stem.")
    .value("VOIDS", bitvf::io::Product::Voids)
    .value("OVERDENSE", bitvf::io::Product::Overdense)
    .value("DIVERGENCE", bitvf::io::Product::Divergence)
    .value("GALAXIES", bitvf::io::Product::Galaxies);

  io.def("file_name",
      [] (nb::handle product, nb::handle directory, nb::handle name, nb::handle zspace) {
        if (!nb::isinstance<bitvf::io::Product>(product))
          throw bitvf::Error("product must be a bitvf.io.Product; got " + repr_of(product));
        const bitvf::io::Product p = nb::cast<bitvf::io::Product>(product);
        const std::string d = to_path(directory, "directory");
        const std::string n = to_string(name, "name");
        const bool z = to_bool(zspace, "zspace");
        const std::string file = bitvf::io::fileName(p, d, n, z);
        return nb::str(file.data(), file.size());
      },
      "product"_a.none(), "directory"_a.none(), "name"_a.none(), nb::kw_only(), "zspace"_a.none(),
      "The conventional name of the file of a product: <directory>/<stem>_[zspace_]<name>.");

  auto writer = [&io] (const char* name, void (*write) (const std::string&, const bitvf::Result&),
                       const char* doc) {
    io.def(name,
        [write] (nb::handle file, nb::handle result) {
          const std::string path = to_path(file, "file");
          const bitvf::Result& r = instance<bitvf::Result>(result, "result", "a bitvf.Result");
          nb::gil_scoped_release released;
          write(path, r);
        },
        "file"_a.none(), "result"_a.none(), doc);
  };
  writer("write_voids", &bitvf::io::writeVoids, "Write the voids, one row per void.");
  writer("write_overdense", &bitvf::io::writeOverdense, "Write the overdense basins, one row per basin.");
  writer("write_divergence", &bitvf::io::writeDivergence, "Write the divergence grid, one row per cell.");
  writer("write_galaxies", &bitvf::io::writeGalaxies, "Write where each tracer ended up, one row per tracer.");

  nb::class_<bitvf::io::Field>(io, "Field", "A field at the centres of cells, for find_voids_box_from_field. "
                               "Read-only.")
    .def_prop_ro("centres", [] (FieldHandle self) {
        return triplets(self, nb::cast<const bitvf::io::Field&>(self).centres);
      }, "The centre of each cell, in coordinate units, shape (n, 3).")
    .def_prop_ro("values", [] (FieldHandle self) {
        return rows<double>(self, nb::cast<const bitvf::io::Field&>(self).values);
      }, "The field at each centre, shape (n,).")
    .def_prop_ro("weights", [] (FieldHandle self) -> nb::object {
        const auto& f = nb::cast<const bitvf::io::Field&>(self);
        if (f.weights.empty()) return nb::none();
        return nb::cast(rows<double>(self, f.weights));
      }, "The weight of each cell, shape (n,); None when no weight column was named.");

  io.def("read_divergence",
      [] (nb::handle file, nb::handle angleUnit) {
        const std::string path = to_path(file, "file");
        const bool given = !angleUnit.is_none();
        if (given) in_degrees(angleUnit);
        bitvf::DivergenceGrid grid;
        {
          nb::gil_scoped_release released;
          grid = bitvf::io::readDivergence(path);
        }
        if (grid.geometry == otswap::Geometry::Box && given)
          throw bitvf::Error("angle_unit must not be given for a box file; a box file holds no angles");
        if (grid.geometry != otswap::Geometry::Box && !given)
          throw bitvf::Error("angle_unit is required for a lightcone file");
        return grid;
      },
      "file"_a.none(), nb::kw_only(), "angle_unit"_a.none() = nb::none(),
      "Read a divergence file written by write_divergence; angle_unit (\"deg\" or \"rad\") is that of a "
      "lightcone grid's angles, required for a lightcone file and an error for a box file.");

  io.def("read_field",
      [] (nb::handle file, nb::handle centres, nb::handle value, nb::handle weight, nb::handle delimiter,
          nb::handle comment) {
        const std::string path = to_path(file, "file");
        std::array<std::string, 3> names;
        std::size_t count = 0;
        bool named = !nb::isinstance<nb::str>(centres);
        if (named) {
          try {
            for (nb::handle item : nb::iter(centres)) {
              if (count == 3 || !nb::isinstance<nb::str>(item)) { named = false; break; }
              names[count++] = nb::str(item).c_str();
            }
          }
          catch (nb::python_error&) {
            named = false;
          }
        }
        if (!named || count != 3)
          throw bitvf::Error("centres must be a sequence of three column names (str); got " + repr_of(centres));
        const std::string v = to_string(value, "value");
        const std::string w = weight.is_none() ? std::string() : to_string(weight, "weight");
        const char d = to_char(delimiter, "delimiter"), k = to_char(comment, "comment");
        nb::gil_scoped_release released;
        return bitvf::io::readField(path, names, v, w, d, k);
      },
      "file"_a.none(), nb::kw_only(), "centres"_a.none(), "value"_a.none(), "weight"_a.none() = nb::none(),
      "delimiter"_a.none() = " ", "comment"_a.none() = "#",
      "Read a field at cell centres from a table, for find_voids_box_from_field.");
}
