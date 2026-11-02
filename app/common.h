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
 *  @file app/common.h
 *
 *  @brief What the two programs, bitvf and bitvf-field, share: the
 *  command line, the typed reads of the parameter file with their
 *  messages, the scan of the keys that do not apply, and the output
 *  directory.
 *
 *  Compiled into both programs; not installed. Uses the public headers
 *  only, as an outside program would.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef BITVF_COMMON_H
#define BITVF_COMMON_H

#include <algorithm>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <functional>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <vector>

#include "otswap/OT.h"
#include "otswap/ParameterFile.h"

namespace app {

  // ==========================================================================
  // Errors
  // ==========================================================================

  /// An error of the program itself: printed after "<program>: error: ".
  class Error : public std::runtime_error {
  public:
    explicit Error (const std::string& what) : std::runtime_error(what) {}
  };

  /// Calls f, an otswap function the program calls itself; an otswap::Error
  /// is rethrown with its message prefixed "otswap: ", as the library does.
  template <typename F>
  auto call_otswap (F&& f)
  {
    try {
      return f();
    }
    catch (const otswap::Error& e) {
      throw Error(std::string("otswap: ") + e.what());
    }
  }

  // ==========================================================================
  // The keys that do not apply
  // ==========================================================================

  /// Every key no read has marked, in file order, each with the message of
  /// rule: empty when the key applies to the run, "unknown key" or
  /// "does not apply to ..." when it does not.
  inline std::vector<std::string> irrelevant_keys (const otswap::io::ParameterFile& parameters,
                                                   const std::string& file,
                                                   const std::function<std::string (const std::string&)>& rule)
  {
    std::vector<std::string> messages;
    for (const otswap::io::ParameterFile::Entry& entry : parameters.unread()) {
      const std::string message = rule(entry.key);
      if (!message.empty())
        messages.push_back(file + ":" + std::to_string(entry.line) + ": key '" + entry.key + "': " + message);
    }
    return messages;
  }

  /// Prints each message as an error of the program; true if there was one.
  inline bool report_errors (const std::string& program, const std::vector<std::string>& messages)
  {
    for (const std::string& message : messages) std::cerr << program << ": error: " << message << '\n';
    return !messages.empty();
  }

  /// The keys left unread after every typed read: a key the program knows
  /// but did not read is a bug of the program, reported as unknown.
  inline bool report_unread (const std::string& program, const otswap::io::ParameterFile& parameters,
                             const std::string& file)
  {
    return report_errors(program, irrelevant_keys(parameters, file, [] (const std::string&) {
      return std::string("unknown key");
    }));
  }

  // ==========================================================================
  // Typed reads, with the programs' own checks
  // ==========================================================================

  /// The parameter file with the helpers of the programs' messages, which
  /// name the file, the line and the key, and quote a value as written in
  /// the file.
  class Parameters {
  public:

    Parameters (const otswap::io::ParameterFile& parameters, const std::string& file)
      : m_parameters(parameters), m_file(file) {}

    const otswap::io::ParameterFile& file () const { return m_parameters; }

    const std::string& name () const { return m_file; }

    std::string at (const std::string& key) const
    {
      return m_file + ":" + std::to_string(m_parameters.line(key)) + ": key '" + key + "': ";
    }

    std::string line (const std::string& key) const { return std::to_string(m_parameters.line(key)); }

    [[noreturn]] void fail (const std::string& key, const std::string& message) const
    {
      throw Error(at(key) + message);
    }

    /// The value of a scalar key as written, read after its typed read.
    std::string written (const std::string& key) const { return m_parameters.get<std::string>(key); }

    /// Element i, from 0, of a list key as written.
    std::string written (const std::string& key, const std::size_t i) const
    {
      return m_parameters.get<std::vector<std::string>>(key)[i];
    }

    /// A value outside the range the library accepts.
    void require (const std::string& key, const bool inRange, const std::string& rule) const
    {
      if (!inRange) fail(key, written(key) + " is out of range; it must be " + rule);
    }

    /// A list key of three elements.
    template <typename T>
    std::vector<T> three (const std::string& key) const
    {
      const std::vector<T> values = m_parameters.get<std::vector<T>>(key);
      if (values.size() != 3)
        fail(key, "expected 3 values, found " + std::to_string(values.size()));
      return values;
    }

    /// A key that may be left without a value, "key =", which means none,
    /// as leaving it out: empty when absent or empty.
    template <typename T>
    std::optional<T> optional (const std::string& key) const
    {
      if (!m_parameters.has(key)) return std::nullopt;
      const std::vector<T> values = m_parameters.get<std::vector<T>>(key);
      if (values.size() > 1)
        fail(key, "expected 1 value or none, found " + std::to_string(values.size()));
      if (values.empty()) return std::nullopt;
      return values[0];
    }

    /// A list key of three elements, or none: empty when absent or empty.
    std::vector<double> optional_three (const std::string& key) const
    {
      if (!m_parameters.has(key)) return {};
      const std::vector<double> values = m_parameters.get<std::vector<double>>(key);
      if (!values.empty() && values.size() != 3)
        fail(key, "expected 3 values or none, found " + std::to_string(values.size()));
      return values;
    }

    /// A key whose value is one of a few words; without a default the key
    /// is required.
    std::string word (const std::string& key, const std::vector<std::string>& words,
                      const char* defaultWord = nullptr) const
    {
      const std::string value = defaultWord == nullptr ? m_parameters.get<std::string>(key)
                                                       : m_parameters.get<std::string>(key, defaultWord);
      if (std::find(words.begin(), words.end(), value) != words.end()) return value;
      std::string list;
      for (std::size_t i = 0; i < words.size(); ++i)
        list += (i == 0 ? "" : i+1 == words.size() ? " or " : ", ") + words[i];
      fail(key, value + " is not " + list);
    }

    /// The verbosity key, normal when absent.
    otswap::Verbosity verbosity () const
    {
      const std::string v = word("verbosity", {"silent", "normal", "detailed"}, "normal");
      return v == "silent" ? otswap::Verbosity::Silent
           : v == "normal" ? otswap::Verbosity::Normal : otswap::Verbosity::Detailed;
    }

    /// The output name: non-empty and without '/'.
    std::string output_file () const
    {
      const std::string name = m_parameters.get<std::string>("output_file");
      if (name.empty() || name.find('/') != std::string::npos)
        fail("output_file", "\"" + name + "\" is not a file name; it must be non-empty and hold no '/'");
      return name;
    }

  private:

    const otswap::io::ParameterFile& m_parameters;
    const std::string& m_file;
  };

  /// True if otswap::io::read reads the file as FITS: the rule of its
  /// extension.
  inline bool is_fits (const std::string& file)
  {
    auto ends_with = [&file] (const std::string& tail) {
      if (file.size() < tail.size()) return false;
      for (std::size_t i = 0; i < tail.size(); ++i) {
        char c = file[file.size()-tail.size()+i];
        if (c >= 'A' && c <= 'Z') c = char(c-'A'+'a');
        if (c != tail[i]) return false;
      }
      return true;
    };
    return ends_with(".fits") || ends_with(".fit") || ends_with(".fits.gz");
  }

  /// True if two column selectors of otswap::io::read name the same column,
  /// as far as that can be told without reading the file. In ASCII two
  /// indices are equal when their values are, two names when they are
  /// identical; an index and a name cannot be compared here. In FITS every
  /// selector is a name, compared ignoring case.
  inline bool same_column (const std::string& a, const std::string& b, const bool fits)
  {
    auto lower = [] (std::string s) {
      for (char& c : s) if (c >= 'A' && c <= 'Z') c = char(c-'A'+'a');
      return s;
    };
    auto all_digits = [] (const std::string& s) {
      return !s.empty() && std::all_of(s.begin(), s.end(), [] (const char c) { return c >= '0' && c <= '9'; });
    };
    if (fits) return lower(a) == lower(b);
    if (all_digits(a) && all_digits(b)) {
      auto strip = [] (const std::string& s) {
        const std::size_t first = s.find_first_not_of('0');
        return first == std::string::npos ? std::string("0") : s.substr(first);
      };
      return strip(a) == strip(b);
    }
    return a == b;
  }

  // ==========================================================================
  // The run
  // ==========================================================================

  /// The first ID met twice, in row order, with both data rows counted
  /// from 1.
  inline void check_ids (const std::vector<std::int64_t>& ids, const std::string& file, const std::string& column)
  {
    std::unordered_map<std::int64_t, std::size_t> first;
    for (std::size_t row = 0; row < ids.size(); ++row) {
      const auto found = first.emplace(ids[row], row);
      if (!found.second)
        throw Error(file + ": the column " + column + " holds the ID " + std::to_string(ids[row]) +
                    " twice, at data rows " + std::to_string(found.first->second+1) + " and " +
                    std::to_string(row+1));
    }
  }

  /// Creates the output directory and its parents, if missing.
  inline void create_directory (const Parameters& p, const std::string& directory)
  {
    std::error_code error;
    const std::filesystem::file_status status = std::filesystem::status(directory, error);
    if (std::filesystem::exists(status) && !std::filesystem::is_directory(status))
      p.fail("output_dir", "cannot create the directory " + directory + ": not a directory");
    std::filesystem::create_directories(directory, error);
    if (error)
      p.fail("output_dir", "cannot create the directory " + directory + ": " + error.message());
  }

  /// The help of a program, on stdout.
  inline std::string help (const std::string& program)
  {
    return "usage: " + program + " <parameter-file>\n"
           "  --help     print this help and exit\n"
           "  --version  print the version and exit\n"
           "The keys of the parameter file are listed at https://simosartori.github.io/bitvf/executable.html\n";
  }

  /// The command line, then run_file on the parameter file. Exit status 0
  /// when run_file returns 0, 1 on any error, each printed as one line,
  /// "<program>: error: <message>".
  inline int main (const int argc, char** argv, const std::string& program, const std::string& version,
                   const std::function<int (const std::string&)>& run_file)
  {
    const std::string usage = "usage: " + program + " <parameter-file>";
    if (argc != 2) {
      std::cerr << usage << '\n';
      return 1;
    }
    const std::string argument = argv[1];
    if (argument == "--help") {
      std::cout << help(program);
      return 0;
    }
    if (argument == "--version") {
      std::cout << program << ' ' << version << '\n';
      return 0;
    }
    if (argument.compare(0, 2, "--") == 0) {
      std::cerr << usage << '\n';
      return 1;
    }

    try {
      return run_file(argument);
    }
    catch (const std::exception& e) {
      std::cerr << program << ": error: " << e.what() << '\n';
    }
    catch (...) {
      std::cerr << program << ": error: unknown exception\n";
    }
    return 1;
  }

}

#endif
