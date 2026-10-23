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
 *  @file src/VoidFinderBox.cpp
 *
 *  @brief The box void finder, findVoidsBox.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#include "otswap/Maths.h"

#include "internal.h"

namespace {

  /**
   *  @brief The search shared by both overloads, once the arguments are
   *  checked: the grid, then the steps of the search.
   *  @param grid the grid to reuse, or null to build one
   */
  bitvf::Result search_box (const bitvf::internal::BoxSearch& search, const bitvf::Config& config,
                            const bitvf::DivergenceGrid* grid)
  {
    if (grid == nullptr) {
      const double cellSizeMps = config.divergenceCellSize
        ? *config.divergenceCellSize
        : otswap::maths::cbrt(20. / search.nRealizations);
      bitvf::internal::box_grid(search, cellSizeMps);
    }

    throw bitvf::Error("findVoidsBox is not implemented yet: the arguments are checked and the "
                       "grid is built, the search is not");
  }

}


// ============================================================================


bitvf::Result bitvf::findVoidsBox (const std::vector<double>& tracers,
                                   const std::vector<double>& displacement,
                                   const std::vector<std::uint8_t>& valid,
                                   const std::vector<double>& randoms,
                                   const double mps,
                                   const Config& config,
                                   const std::optional<otswap::BoxLimits>& limits,
                                   const std::vector<std::int64_t>& ids,
                                   const DivergenceGrid* grid)
{
  const internal::BoxSearch search =
    internal::check_box(tracers, displacement, valid, randoms, mps, config, limits, ids, grid);
  return search_box(search, config, grid);
}


// ============================================================================


bitvf::Result bitvf::findVoidsBox (const otswap::Result& reconstruction,
                                   const Config& config,
                                   const std::vector<double>& randoms,
                                   const DivergenceGrid* grid)
{
  const internal::BoxSearch search = internal::check_box(reconstruction, config, randoms, grid);
  return search_box(search, config, grid);
}
