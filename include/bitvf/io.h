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
 *  @file include/bitvf/io.h
 *
 *  @brief The writers of the results of the void finder, and the readers
 *  of a divergence file and of otswap's mean displacement file.
 *
 *  Each writer writes one product of a Result, a table with one row per
 *  void, overdense basin, cell or input tracer. The format follows the
 *  extension of the file name, as in otswap::io: FITS for .fits, .fit and
 *  .fits.gz, ASCII for anything else, with the layouts otswap::io::write
 *  produces. An existing file is overwritten. Every check runs before the
 *  file is opened, so a call that throws leaves no file.
 *
 *  Lengths are in coordinate units, which bitvf cannot know, so their
 *  columns carry no unit and their descriptions say "in coordinate units".
 *  Index columns hold rows of the products: index is the row of the file
 *  itself, voidIndex, overdenseIndex and cellIndex the rows of the voids,
 *  overdense basins and divergence files of the same result, -1 for none.
 *  Counts and indices are 64-bit integers (FITS K), flags 32-bit integers
 *  (FITS J) holding 0 or 1, and every other column a double (FITS D),
 *  written in ASCII with 9 significant digits, or 17 in the divergence
 *  file, which then reads back as the same double. FITS stores the doubles
 *  themselves. NaN is written as NaN.
 *
 *  Every file starts with these keywords, in this order:
 *  - PRODUCT, what the file holds: voids, overdense_basins, divergence or
 *    galaxies;
 *  - BITVFV, the bitvf version that wrote it;
 *  - GEOMETRY, box;
 *  - the row count: NBASINS (voids, overdense basins), NCELLS (divergence)
 *    or NOBJECTS (galaxies, the input tracers);
 *  - NREC, Result::nRealizations;
 *  - NRECMEAN, Config::nRealizationsMean of the search;
 *  - SIGMA and SIGUNIT, the smoothing sigma and its unit, mps or
 *    coordinate;
 *  - CELLSIZE, the divergence cell size used, in units of the mps;
 *  - MPS, the mean particle separation, in coordinate units;
 *  - RANDSEED, Result::randomSeed, when set;
 *  - BOXLIM = 1, then BOXXLO, BOXXHI, BOXYLO, BOXYHI, BOXZLO, BOXZHI, the
 *    box limits of the search in coordinate units, when Result::limits is
 *    set.
 *
 *  A number is written as the shortest decimal of at most 17 significant
 *  digits that reads back as the same double, an integer as an integer;
 *  otswap::io::readKeywords returns them as written, in both formats.
 *
 *  Only box results are written for now.
 *
 *  @author Simone Sartori <simone.sartori@inaf.it>
 */

#ifndef BITVF_IO_H
#define BITVF_IO_H

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "bitvf/VoidFinder.h"

namespace bitvf::io {

  /// The products, each with its file-name stem: voids,
  /// overdense_basins, divergence, galaxies.
  enum class Product { Voids, Overdense, Divergence, Galaxies };

  /**
   *  @brief The conventional name of the file of a product.
   *
   *  <directory>/<stem>_[zspace_]<name>, with the stem of the product. A
   *  '/' is added only when directory is not empty and does not end in
   *  one; an empty directory gives <stem>_[zspace_]<name>. The extension,
   *  which sets the format, is part of name.
   *
   *  @param zspace true only for the redshift-space pass of a run with
   *  redshift-space correction, whose files then do not collide with those
   *  of the real-space pass.
   *
   *  @exception Error if name is empty.
   */
  std::string fileName (Product product, const std::string& directory, const std::string& name, bool zspace);

  /**
   *  @brief Write the voids, one row per void (Result::voids).
   *
   *  Columns: index (K); centreX, centreY, centreZ, in coordinate units;
   *  radius, radiusExt, in coordinate units; deltaAv; delta2Mps; edgeFlag
   *  (J); contamination; nTracers (K). Keywords as listed in the file
   *  description, the row count NBASINS. ASCII: 9 significant digits.
   *
   *  @exception Error
   *  - if the result is not a box result;
   *  - if it is malformed: Config::divergenceCellSize empty,
   *    Config::nRealizationsMean 0, a column of the catalogue not of
   *    nBasins entries (centre 3 * nBasins), centreSky not empty, a count
   *    above the largest 64-bit integer;
   *  - if the file cannot be written: the message starts with "otswap: ",
   *    the otswap::Error nested.
   */
  void writeVoids (const std::string& file, const Result& result);

  /**
   *  @brief Write the overdense basins, one row per basin
   *  (Result::overdense), with the columns, keywords and errors of
   *  writeVoids.
   */
  void writeOverdense (const std::string& file, const Result& result);

  /**
   *  @brief Write the divergence grid, one row per cell
   *  (Result::divergence).
   *
   *  Columns: index (K); cellX, cellY, cellZ, the cell centre, in
   *  coordinate units; divergence, after smoothing; divergenceUnsmoothed,
   *  only when DivergenceGrid::divergenceUnsmoothed is not empty; volume, in
   *  coordinate units cubed; weight; nDispl (K), the displacement ends
   *  counted; nTracers (K); voidIndex (K); overdenseIndex (K). ASCII: 17
   *  significant digits.
   *
   *  Keywords as listed in the file description, the row count NCELLS,
   *  followed by the record of the grid: LOWERX, LOWERY, LOWERZ, UPPERX,
   *  UPPERY, UPPERZ, CELLSIDE, in coordinate units, and NCELLX, NCELLY,
   *  NCELLZ. With MPS and CELLSIZE, they hold every field that defines the
   *  cells at 17 significant digits, so the grid read back is the same bit
   *  for bit.
   *
   *  @exception Error
   *  - if the result is not a box result, or holds no divergence grid
   *    (Config::keepDivergence);
   *  - if the grid fails the checks of a grid to reuse (see findVoidsBox),
   *    the message starting with "cannot write the divergence: ";
   *  - if it is malformed: as writeVoids; a column of the grid not of
   *    nCells entries (centre 3 * nCells; divergenceUnsmoothed nCells or
   *    none); the mps, cellSizeMps or limits of the grid different from
   *    Result::mps, Config::divergenceCellSize or Result::limits, bit for
   *    bit; its nRealizations or nRealizationsMean different from
   *    Result::nRealizations or Config::nRealizationsMean;
   *  - if the file cannot be written, as writeVoids.
   */
  void writeDivergence (const std::string& file, const Result& result);

  /**
   *  @brief Write where each tracer ended up, one row per input tracer,
   *  in input order (Result::assignment).
   *
   *  Columns: id (K), Result::ids, or the row from 0 when it is empty;
   *  cellIndex (K); voidIndex (K); overdenseIndex (K); distance, from the
   *  centre of the tracer's basin, in coordinate units, NaN for a tracer
   *  not in a basin; inMask (J); inBox (J), only when Result::limits is
   *  set; used (J); uncorrected (J), only when Result::uncorrected is not
   *  empty; outsideGrid (J), only when Result::outsideGrid is not empty.
   *  Keywords as listed in the file description, the row count NOBJECTS.
   *  ASCII: 9 significant digits.
   *
   *  @exception Error
   *  - if the result is not a box result, or holds no tracer assignment
   *    (Config::keepAssignment);
   *  - if it is malformed: Config::divergenceCellSize empty,
   *    Config::nRealizationsMean 0; ids, uncorrected or outsideGrid
   *    neither empty nor of nTracers entries; a column of the assignment
   *    not of nTracers entries;
   *  - if the file cannot be written, as writeVoids.
   */
  void writeGalaxies (const std::string& file, const Result& result);

  /**
   *  @brief otswap's mean displacement file
   *  (otswap::io::writeDisplacementMean), ready for the array overload of
   *  findVoidsBox as one realization.
   */
  struct DisplacementMean {

    /// Tracer positions, in coordinate units: columns tracX, tracY, tracZ.
    /// Flat, [tracer][xyz].
    std::vector<double> tracers;

    /// Mean displacement of each tracer, in coordinate units: columns
    /// displX, displY, displZ, as stored, NaN where no realization was
    /// valid. Flat, [tracer][xyz].
    std::vector<double> displacement;

    /// Validity of each displacement: 1 where column nValidRec is
    /// positive, 0 otherwise. [tracer].
    std::vector<std::uint8_t> valid;

    /// Mean particle separation of the reconstruction, in coordinate units:
    /// keyword MPS.
    double mps = std::numeric_limits<double>::quiet_NaN();

    /// The box limits of the reconstruction, keywords BOXXLO to BOXZHI,
    /// when the file records them as limits (BOXLIM = 1); empty otherwise.
    std::optional<otswap::BoxLimits> limits;

    /// Column id when the file has one; empty otherwise.
    std::vector<std::int64_t> ids;

    /// Realizations averaged into each displacement: keyword NREC.
    unsigned nRealizationsAveraged = 0;
  };

  /**
   *  @brief Read otswap's mean displacement file of a box reconstruction.
   *
   *  The fields are filled for the array overload of findVoidsBox, which
   *  searches them as one realization:
   *
   *      const bitvf::io::DisplacementMean m = bitvf::io::readDisplacementMean("displacement_mean.fits");
   *      config.nRealizationsMean = m.nRealizationsAveraged;
   *      const bitvf::Result r = bitvf::findVoidsBox(m.tracers, m.displacement, m.valid, m.mps, config,
   *                                                  m.limits, m.ids);
   *
   *  FITS stores the doubles themselves; ASCII holds tracers and
   *  displacements with 9 significant digits, so they are rounded, and a
   *  tracer within rounding of a box limit can fall on the other side of
   *  it. Keywords are exact in both formats. The box keywords are limits
   *  only with BOXLIM = 1; otherwise the search builds its grid on the
   *  tracers and the endpoints of their displacements. The default
   *  Config::divergenceCellSize of such a search is 1, whatever the
   *  realizations averaged.
   *
   *  The divergence is not linear in the displacement, so a search on the
   *  mean of R realizations does not give the result of the search on the
   *  R realizations. Set Config::nRealizationsMean to nRealizationsAveraged,
   *  so that bitvf's files record it.
   *
   *  Reading a lightcone mean displacement file is not implemented.
   *
   *  @exception Error
   *  - if the file cannot be read, or lacks a column, or an integer column
   *    holds anything but an integer: the message starts with "otswap: ",
   *    the otswap::Error nested;
   *  - if it has no PRODUCT keyword, its PRODUCT is not displacements, or
   *    GEOMETRY is not box;
   *  - if NREC or MPS is missing or not a number, NREC is 0, BOXLIM holds
   *    anything but 1, or BOXLIM = 1 comes without the six box keywords;
   *  - if nValidRec lies outside 0 to NREC in a row.
   */
  DisplacementMean readDisplacementMean (const std::string& file);

  /**
   *  @brief Read a divergence file written by writeDivergence.
   *
   *  Every field the writer wrote is read back bit for bit, in both
   *  formats: the cells (lower, upper, cellSize, nCellsAxis, nCells, mps,
   *  cellSizeMps), the box limits, nRealizations, nRealizationsMean and
   *  every per-cell column, divergenceUnsmoothed when the file has it.
   *  The grid passes the checks of a grid to reuse, so it can be passed to
   *  findVoidsBox as grid; the per-search fields are records, which a
   *  search reusing the grid does not read.
   *
   *  Reading a lightcone divergence file is not implemented.
   *
   *  @exception Error
   *  - if the file cannot be read, or lacks a column, as
   *    readDisplacementMean;
   *  - if it has no PRODUCT keyword, its PRODUCT is not divergence, or
   *    GEOMETRY is not box;
   *  - if a keyword of the grid record is missing, not a number, or not a
   *    whole number in range; if NREC or NRECMEAN is 0; if BOXLIM holds
   *    anything but 1, or BOXLIM and the six box keywords do not come
   *    together;
   *  - if the rows are not NCELLS, index is not 0 to NCELLS - 1 in order,
   *    a count is negative, or a basin index is below -1;
   *  - if the grid fails the checks of a grid to reuse.
   */
  DivergenceGrid readDivergence (const std::string& file);

}

#endif
