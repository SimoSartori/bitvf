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
 *  @file src/Support.cpp
 *
 *  @brief Definitions shared by the bitvf translation units.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <iomanip>
#include <locale>

#include "internal.h"

const char* const bitvf::internal::version = BITVF_VERSION;


// ============================================================================


bitvf::internal::Report::Report (const otswap::Verbosity verbosity, std::ostream* log)
  : m_verbosity(verbosity), m_log(log), m_start(std::chrono::steady_clock::now()), m_last(m_start)
{}


// ============================================================================


void bitvf::internal::Report::write (const std::string& text) const
{
  if (!normal()) return;
  *m_log << ("bitvf: " + text + "\n") << std::flush;
}


// ============================================================================


double bitvf::internal::Report::lap ()
{
  const auto now = std::chrono::steady_clock::now();
  const double seconds = std::chrono::duration<double>(now - m_last).count();
  m_last = now;
  return seconds;
}


// ============================================================================


double bitvf::internal::Report::elapsed () const
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - m_start).count();
}


// ============================================================================


void bitvf::internal::Report::step (std::ostringstream& out)
{
  const double seconds = lap();
  out << ' ';
  time(out, seconds);
  write(out.str());
}


// ============================================================================


std::ostringstream bitvf::internal::Report::line ()
{
  std::ostringstream out;
  out.imbue(std::locale::classic());
  return out;
}


// ============================================================================


void bitvf::internal::Report::time (std::ostream& out, const double t)
{
  const std::ios_base::fmtflags flags = out.flags();
  const std::streamsize precision = out.precision();
  out << "in " << std::fixed << std::setprecision(2) << t << " s";
  out.flags(flags);
  out.precision(precision);
}


// ============================================================================


void bitvf::internal::Report::box (std::ostream& out, const std::array<double, 3>& lower,
                                   const std::array<double, 3>& upper)
{
  for (std::size_t k = 0; k < 3; ++k)
    out << (k ? " x [" : "[") << lower[k] << ", " << upper[k] << "]";
}


// ============================================================================


void bitvf::internal::Report::shape (std::ostream& out, const std::array<std::size_t, 3>& n)
{
  out << n[0] << " x " << n[1] << " x " << n[2] << " cells";
}


// ============================================================================


void bitvf::internal::Report::basins (std::ostream& out, const std::size_t nVoids, const std::size_t nOverdense)
{
  out << nVoids << (nVoids == 1 ? " void and " : " voids and ") << nOverdense
      << (nOverdense == 1 ? " overdense basin" : " overdense basins");
}
