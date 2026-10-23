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
 *  @brief nanobind module bitvf._bitvf, the compiled part of the bitvf
 *  package. It exposes the library version as __version__.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include <nanobind/nanobind.h>

#include "internal.h"

NB_MODULE(_bitvf, m)
{
  m.doc() = "Compiled part of the bitvf package; import bitvf instead.";

  m.attr("__version__") = bitvf::internal::version;
}
