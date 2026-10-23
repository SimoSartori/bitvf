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
 *  @file tests/check.h
 *
 *  @brief The checks shared by the C++ tests.
 */

#ifndef BITVF_TESTS_CHECK_H
#define BITVF_TESTS_CHECK_H

#include <exception>
#include <iostream>
#include <string>

#include "bitvf/VoidFinder.h"

namespace {

  int failures = 0;

  void check (const bool condition, const std::string& what)
  {
    if (condition) return;
    std::cerr << "FAILED: " << what << std::endl;
    ++failures;
  }

  /// @p f must raise a bitvf::Error whose message contains @p text.
  template <typename F>
  void throws_naming (F&& f, const std::string& text, const std::string& what)
  {
    try {
      f();
    }
    catch (const bitvf::Error& e) {
      if (std::string(e.what()).find(text) != std::string::npos) return;
      std::cerr << "FAILED: " << what << " (the message is \"" << e.what() << "\", without \""
                << text << "\")" << std::endl;
      ++failures;
      return;
    }
    catch (const std::exception& e) {
      std::cerr << "FAILED: " << what << " (threw the wrong type: " << e.what() << ")"
                << std::endl;
      ++failures;
      return;
    }
    std::cerr << "FAILED: " << what << " (did not throw)" << std::endl;
    ++failures;
  }

  void group (const std::string& what)
  {
    std::cout << "  " << what << std::endl;
  }

  int report (const std::string& name)
  {
    if (failures == 0) {
      std::cout << name << ": all groups passed" << std::endl;
      return 0;
    }
    std::cerr << name << ": " << failures << " failure(s)" << std::endl;
    return 1;
  }

}

#endif
