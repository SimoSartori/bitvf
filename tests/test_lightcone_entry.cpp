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
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "bitvf/VoidFinder.h"
#include "otswap/io.h"

#include "check.h"
#include "internal.h"

using namespace bitvf;

namespace {

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
  constexpr double kInf = std::numeric_limits<double>::infinity();

  struct Case {

    std::vector<double> tracersSky {1., 0.2, 0.3,  2., -0.2, 0.35,  3., 0.1, 0.4,  4., 0.3, 0.25};
    std::vector<double> displacement = std::vector<double>(12, 0.1);
    std::vector<std::uint8_t> valid;
    std::vector<double> randomsSky {1.5, 0.1, 0.3,  2.5, 0., 0.3,  3.5, 0.2, 0.35,  4.5, -0.1, 0.4,
                                    5., 0.2, 0.3,  6., 0.1, 0.25,  0.5, 0.3, 0.3,  0.2, -0.3, 0.35};
    std::vector<double> tracers;
    std::vector<double> randoms;
    otswap::Mask mask {std::vector<double>(12, 1.), otswap::PixelOrdering::Ring};
    otswap::DistanceTable distances {0.3, 0.7, -1., 0., 0., 2., 1000};
    otswap::MpsProfile profile;
    Config config;
    otswap::RedshiftCut cut;
    std::vector<std::int64_t> ids;

    Case ()
    {
      profile.redshift = {0.1, 0.5};
      profile.mps = {5., 6.};
      config.verbosity = otswap::Verbosity::Silent;
      tracers = otswap::toCartesian(tracersSky, distances);
      randoms = otswap::toCartesian(randomsSky, distances);
    }

    void sky () const
    {
      findVoidsLightcone(tracersSky, displacement, valid, randomsSky, mask, distances, profile, config, cut, ids);
    }

    void cartesian () const
    {
      findVoidsLightcone(tracers, tracersSky, displacement, valid, randoms, randomsSky, mask, distances, profile,
                         config, cut, ids);
    }
  };

  void fails (const Case& c, const std::string& text, const std::string& what)
  {
    throws_naming([&] { c.sky(); }, text, what + " (sky)");
    throws_naming([&] { c.cartesian(); }, text, what + " (Cartesian)");
  }

  void fails_sky (const Case& c, const std::string& text, const std::string& what)
  {
    throws_naming([&] { c.sky(); }, text, what);
  }

  struct FromResult {

    otswap::Result rec;
    std::vector<double> randomsSky;
    otswap::Mask mask {std::vector<double>(12, 1.), otswap::PixelOrdering::Ring};
    Config config;

    FromResult ()
    {
      const Case c;
      rec.geometry = otswap::Geometry::Lightcone;
      rec.nObjects = 4;
      rec.nRealizations = 1;
      rec.config.nRealizations = 1;
      rec.tracersSky = c.tracersSky;
      rec.tracers = c.tracers;
      rec.displacement = c.displacement;
      rec.valid.assign(4, 1);
      rec.outsideMask.assign(4, 0);
      rec.outsideRedshiftCut.assign(4, 0);
      rec.mpsProfile = c.profile;
      rec.mpsProfile.count = {2, 2};
      rec.distances = std::make_shared<const otswap::DistanceTable>(c.distances);
      otswap::recomputeMeans(rec);
      randomsSky = c.randomsSky;
      config.verbosity = otswap::Verbosity::Silent;
    }

    void run () const
    {
      findVoidsLightcone(rec, randomsSky, mask, config);
    }
  };

  void fails (const FromResult& r, const std::string& text, const std::string& what)
  {
    throws_naming([&] { r.run(); }, text, what);
  }

  std::string temporary (const std::string& name)
  {
    const char* dir = std::getenv("TMPDIR");
    std::string base = dir != nullptr ? dir : "/tmp";
    if (!base.empty() && base.back() != '/') base += '/';
    return base + "bitvf_test_lightcone_entry_" + name;
  }

  FromResult mean_read_back ()
  {
    FromResult r;
    r.rec.nRealizations = 3;
    r.rec.config.nRealizations = 3;
    r.rec.displacement = std::vector<double>(36, 0.1);
    r.rec.valid.assign(12, 1);
    r.rec.valid[1] = 0;
    otswap::recomputeMeans(r.rec);
    const std::string file = temporary("mean.fits");
    std::remove(file.c_str());
    otswap::io::writeDisplacementMean(file, r.rec);
    r.rec = otswap::io::readDisplacementMean(file);
    std::remove(file.c_str());
    return r;
  }

}

int main ()
{
  // A complete, consistent entry passes every check and reaches the search.
  group("a valid entry passes every check and stops at the search");
  {
    Case c;
    fails(c, "not implemented yet", "valid entry");
  }

  // Empty position arrays are refused for tracers and for randoms.
  group("an empty array of sky positions is refused");
  {
    Case c;
    c.tracersSky.clear();
    fails_sky(c, "findVoidsLightcone: tracersSky is empty", "empty tracersSky");
    Case r;
    r.randomsSky.clear();
    fails_sky(r, "findVoidsLightcone: randomsSky is empty", "empty randomsSky");
  }

  // Position arrays whose length is not a multiple of 3 are refused.
  group("a size that is not a multiple of 3 is refused");
  {
    Case c;
    c.tracersSky.pop_back();
    fails_sky(c, "findVoidsLightcone: tracersSky holds 11 entries; it must hold 3 per object", "tracersSky size");
    Case r;
    r.randomsSky.pop_back();
    fails_sky(r, "findVoidsLightcone: randomsSky holds 23 entries; it must hold 3 per object", "randomsSky size");
  }

  // Non-finite angles or redshifts and declinations outside [-pi/2, pi/2] are refused on every row, even rows the cut would drop.
  group("a non-finite component or a declination out of range is refused on any row");
  {
    Case c;
    c.tracersSky[7] = kNaN;
    fails_sky(c, "the sky positions of the tracers hold a non-finite value at entry 7 (object 2, component 1)", "NaN");
    Case a;
    a.tracersSky[3 * 3 + 2] = kInf;
    fails_sky(a, "the sky positions of the tracers hold a non-finite value at entry 11 (object 3, component 2)", "inf");
    Case d;
    d.tracersSky[4] = 1.75;
    fails_sky(d, "the declination of tracer 1 is 1.75 radians; it must lie in [-pi/2, pi/2]", "declination");
    Case r;
    r.randomsSky[1] = -2.;
    fails_sky(r, "the declination of random 0 is -2 radians; it must lie in [-pi/2, pi/2]", "random declination");
    Case m;
    m.cut.min = 0.35;
    m.tracersSky[0] = kNaN;
    fails_sky(m, "the sky positions of the tracers hold a non-finite value at entry 0", "row outside the cut");
  }

  // Cartesian arrays must hold one position per object of the sky arrays.
  group("the Cartesian arrays must match the sky arrays in size");
  {
    Case c;
    c.tracers.pop_back();
    throws_naming([&] { c.cartesian(); }, "findVoidsLightcone: tracers holds 11 entries for 4 objects in tracersSky",
                  "tracers size");
    Case r;
    r.randoms.push_back(0.);
    throws_naming([&] { r.cartesian(); }, "findVoidsLightcone: randoms holds 25 entries for 8 objects in randomsSky",
                  "randoms size");
  }

  // Displacements must be a whole number of realizations of the tracers, and validity flags must be 0 or 1 and fit them.
  group("the displacements and their validity flags must fit the tracers");
  {
    Case c;
    c.displacement.clear();
    fails(c, "findVoidsLightcone: displacement holds 0 entries for 4 tracers", "empty displacement");
    Case d;
    d.displacement.push_back(0.);
    fails(d, "findVoidsLightcone: displacement holds 13 entries for 4 tracers", "displacement size");
    Case v;
    v.valid.assign(3, 1);
    fails(v, "findVoidsLightcone: valid holds 3 entries for 1 realizations of 4 tracers", "valid size");
    Case w;
    w.valid.assign(4, 1);
    w.valid[2] = 2;
    fails(w, "the validity flags hold 2 at entry 2; each must be 0 or 1", "valid value");
  }

  // A redshift cut with a NaN bound or reversed bounds is refused.
  group("a redshift cut with a NaN bound or a lower bound above the upper one is refused");
  {
    Case c;
    c.cut.min = kNaN;
    fails(c, "neither bound may be NaN and the lower bound must not exceed the upper one", "NaN cut");
    Case d;
    d.cut.min = 0.5;
    d.cut.max = 0.25;
    fails(d, "the redshift cut is [0.5, 0.25]; neither bound", "reversed cut");
  }

  // An empty, inconsistent, non-finite, non-increasing or non-positive mps profile is refused.
  group("the mps profile must be non-empty, consistent, finite, increasing and positive");
  {
    Case e;
    e.profile = {};
    fails(e, "the mps profile is empty; it must hold at least one entry", "empty profile");
    Case s;
    s.profile.mps.pop_back();
    fails(s, "findVoidsLightcone: the redshift and mps lists of the profile hold 2 and 1 entries", "sizes");
    Case z;
    z.profile.redshift[1] = kNaN;
    fails(z, "the mps profile holds a redshift of nan at entry 1; it must be finite", "redshift not finite");
    Case m;
    m.profile.mps[0] = kInf;
    fails(m, "the mps profile holds a mean particle separation of inf at entry 0; it must be finite", "mps not finite");
    Case i;
    i.profile.redshift = {0.5, 0.5};
    fails(i, "the redshifts must increase strictly", "not increasing");
    Case p;
    p.profile.mps[1] = 0.;
    fails(p, "the mps profile holds a mean particle separation of 0 at entry 1; it must be positive", "mps not positive");
    Case b;
    b.profile.mps[1] = -1.;
    fails(b, "it must be positive", "negative mps");
  }

  // Identities must number one per tracer and be distinct.
  group("the identities must be one per tracer and distinct");
  {
    Case c;
    c.ids = {1, 2, 3};
    fails(c, "findVoidsLightcone: ids holds 3 entries for 4 tracers; give one per tracer, or none", "ids size");
    Case d;
    d.ids = {1, 2, 3, 2};
    fails(d, "the ids hold 2 twice, at entries 1 and 3", "ids repeated");
  }

  // The configuration is validated by the same rules as in a box.
  group("the configuration is checked as in a box");
  {
    Case s;
    s.config.smoothing.sigma = kNaN;
    fails(s, "the smoothing sigma is", "sigma");
    Case c;
    c.config.divergenceCellSize = -1.;
    fails(c, "the divergence cell size is", "cell size");
    Case n;
    n.config.nRealizationsMean = 0;
    fails(n, "the number of realizations averaged is 0", "realizations averaged");
  }

  // Rows used by the search must lie inside the distance table; rows masked or cut away are never read.
  group("a kept row with a redshift outside the distance table is refused; rows left out are not read");
  {
    Case c;
    c.tracersSky[2] = 3.;
    fails(c, "the redshift of tracer 0 is 3, outside the distance table [0, 2]", "tracer redshift");
    Case r;
    r.randomsSky[8] = -0.5;
    fails(r, "the redshift of random 2 is -0.5, outside the distance table [0, 2]", "random redshift");
    Case m;
    m.tracersSky[2] = 3.;
    m.cut.max = 1.;
    fails(m, "not implemented yet", "tracer outside the cut is not read");
  }

  // Rows used by the search must have finite Cartesian coordinates; masked or cut rows may hold NaN.
  group("a kept Cartesian row that is not finite is refused; rows left out may hold NaN");
  {
    Case c;
    c.tracers[4] = kNaN;
    throws_naming([&] { c.cartesian(); }, "the tracers hold a non-finite value at entry 4 (object 1, component 1)",
                  "tracers NaN");
    Case r;
    r.randoms[17] = kInf;
    throws_naming([&] { r.cartesian(); }, "the randoms hold a non-finite value at entry 17 (object 5, component 2)",
                  "randoms inf");
    Case m;
    m.cut.min = 0.28;
    m.tracers[9] = kNaN;
    throws_naming([&] { m.cartesian(); }, "not implemented yet", "tracer outside the cut may hold NaN");
  }

  // Valid displacements of used tracers must be finite; invalid ones and those of dropped tracers are never read.
  group("a non-finite valid displacement of a kept tracer is refused; the others are not read");
  {
    Case c;
    c.displacement[5] = kNaN;
    fails(c, "the displacements hold a non-finite value at entry 5 (realization 0, tracer 1, component 2)",
          "valid displacement");
    Case i;
    i.displacement[5] = kNaN;
    i.valid = {1, 0, 1, 1};
    fails(i, "not implemented yet", "invalid displacement is not read");
    Case m;
    m.displacement[5] = kNaN;
    m.cut.max = 0.3;
    m.cut.min = 0.28;
    fails(m, "not implemented yet", "displacement of a tracer outside the cut is not read");
  }

  // A call that leaves no tracer, or no valid displacement, to read is refused with the counts of rows dropped.
  group("no tracer kept, or no valid displacement left to read, is refused");
  {
    Case c;
    c.mask = otswap::Mask(std::vector<double>(12, 0.), otswap::PixelOrdering::Ring);
    fails(c, "no tracer is left to read: of 4 tracers, 4 lie on unobserved pixels and 0 outside the redshift cut",
          "all masked");
    Case d;
    d.cut.min = 0.6;
    fails(d, "no tracer is left to read: of 4 tracers, 0 lie on unobserved pixels and 4 outside the redshift cut",
          "all cut");
    Case v;
    v.valid = {0, 0, 0, 1};
    v.cut.min = 0.28;
    fails(v, "no valid displacement is left to read: of 1 valid displacements", "no displacement read");
  }

  // A call that leaves no random is refused with the counts of rows dropped.
  group("no random kept is refused, naming the mask and the cut");
  {
    Case c;
    c.randomsSky = {1.5, 0.1, 0.9,  2.5, 0., 1.};
    c.randoms = otswap::toCartesian(c.randomsSky, c.distances);
    c.cut.max = 0.6;
    fails(c, "no random is left: of 2 randoms, 0 lie on unobserved pixels and 2 outside the redshift cut",
          "randoms cut");
  }

  // At normal verbosity the report counts, for tracers and randoms, the rows kept, masked and cut, masked rows first.
  group("the Normal report gives the selection lines, counting masked rows before cut rows");
  {
    Case c;
    c.cut.min = 0.28;
    std::ostringstream log;
    c.config.verbosity = otswap::Verbosity::Normal;
    c.config.log = &log;
    fails_sky(c, "not implemented yet", "report case");
    check(log.str().find("bitvf: kept 3 of 4 tracers: 0 on unobserved pixels, 1 outside the redshift cut\n") !=
          std::string::npos, "tracer selection line: " + log.str());
    check(log.str().find("bitvf: kept 7 of 8 randoms: 0 on unobserved pixels, 1 outside the redshift cut\n") !=
          std::string::npos, "random selection line: " + log.str());
  }

  // A consistent lightcone reconstruction, or a mean read back from its file, passes every check and reaches the search.
  group("a valid reconstruction or mean passes every check and stops at the search");
  {
    FromResult r;
    fails(r, "not implemented yet", "valid reconstruction");
    const FromResult m = mean_read_back();
    check(m.rec.nRealizations == 1 && m.rec.nRealizationsMean == 3, "mean read back as one realization of 3");
    fails(m, "not implemented yet", "valid mean");
  }

  // A box reconstruction is refused by the lightcone void finder.
  group("a reconstruction that is not a lightcone is refused");
  {
    FromResult r;
    r.rec.geometry = otswap::Geometry::Box;
    fails(r, "findVoidsLightcone: the reconstruction is not a lightcone result; use the box void finder", "box result");
  }

  // A reconstruction without sky positions, distance table or, without a grid, mps profile is refused, naming what is missing.
  group("a reconstruction lacking the sky positions, the distance table or the profile is refused");
  {
    FromResult r;
    r.rec.tracersSky.clear();
    r.rec.distances.reset();
    r.rec.mpsProfile = {};
    fails(r, "findVoidsLightcone: the reconstruction has no tracersSky, distances and mpsProfile, which a lightcone "
             "search needs", "all missing");
    FromResult d;
    d.rec.distances.reset();
    d.rec.nObjects = 0;
    fails(d, "findVoidsLightcone: the reconstruction has no distances, which a lightcone search needs",
          "missing before sizes");
    FromResult p;
    p.rec.mpsProfile = {};
    fails(p, "findVoidsLightcone: the reconstruction has no mpsProfile, which a lightcone search needs; give a grid to "
             "search without mpsProfile", "profile only");
  }

  // Every array of a reconstruction must match its numbers of objects and realizations.
  group("a reconstruction whose sizes disagree is refused as malformed");
  {
    const std::string malformed = "findVoidsLightcone: the reconstruction is malformed: ";
    FromResult n;
    n.rec.nObjects = 0;
    fails(n, malformed + "it holds no objects", "no objects");
    FromResult z;
    z.rec.nRealizations = 0;
    fails(z, malformed + "it holds no realizations", "no realizations");
    FromResult t;
    t.rec.tracers.pop_back();
    fails(t, malformed + "tracers holds 11 entries for 1 realizations of 4 objects", "tracers");
    FromResult s;
    s.rec.tracersSky.resize(15, 0.);
    fails(s, malformed + "tracersSky holds 15 entries", "tracersSky");
    FromResult d;
    d.rec.displacement.push_back(0.);
    fails(d, malformed + "displacement holds 13 entries", "displacement");
    FromResult v;
    v.rec.valid.pop_back();
    fails(v, malformed + "valid holds 3 entries", "valid");
    FromResult m;
    m.rec.outsideMask.assign(3, 0);
    fails(m, malformed + "outsideMask holds 3 entries", "outsideMask");
    FromResult c;
    c.rec.outsideRedshiftCut.assign(5, 0);
    fails(c, malformed + "outsideRedshiftCut holds 5 entries", "outsideRedshiftCut");
    FromResult i;
    i.rec.ids = {1, 2};
    fails(i, malformed + "ids holds 2 entries", "ids");
  }

  // A mean holds one realization, averages at least one, and counts at most that many valid realizations per tracer.
  group("a mean with inconsistent counts of realizations is refused as malformed");
  {
    const std::string malformed = "findVoidsLightcone: the reconstruction is malformed: ";
    FromResult z = mean_read_back();
    z.rec.nRealizationsMean = 0;
    fails(z, malformed + "nRealizationsMean is 0; it must be at least 1", "nothing averaged");
    FromResult r;
    r.rec.nRealizations = 2;
    r.rec.nRealizationsMean = 3;
    fails(r, malformed + "it holds 2 realizations, each the mean of 3; a mean holds one realization", "two means");
    FromResult s = mean_read_back();
    s.rec.validRealizations.pop_back();
    fails(s, malformed + "validRealizations holds 3 entries for 4 objects", "count size");
    FromResult a = mean_read_back();
    a.rec.validRealizations[2] = 4;
    fails(a, malformed + "validRealizations holds 4 for object 2, above nRealizations * nRealizationsMean = 3",
          "count above the mean");
  }

  // Validity flags of a reconstruction must be 0 or 1 and its identities distinct.
  group("a reconstruction with a validity flag other than 0 and 1 or a repeated identity is refused");
  {
    FromResult v;
    v.rec.valid[3] = 2;
    fails(v, "findVoidsLightcone: the reconstruction's validity flags hold 2 at entry 3; each must be 0 or 1",
          "valid value");
    FromResult i;
    i.rec.ids = {5, 6, 5, 7};
    fails(i, "findVoidsLightcone: the reconstruction's ids hold 5 twice, at entries 0 and 2", "ids repeated");
  }

  // A reconstruction filtered against a mask must be searched with a mask of the same resolution.
  group("a reconstruction filtered at another mask resolution is refused");
  {
    FromResult r;
    r.rec.filteredNside = 2;
    fails(r, "findVoidsLightcone: the reconstruction was filtered against a mask of NSIDE 2 "
             "(otswap::Result::filteredNside), and the mask given has NSIDE 1", "NSIDE");
  }

  // Sky positions of a reconstruction, its profile when read, and the randoms are checked as in the array overloads.
  group("the sky positions, the profile and the randoms are checked on every row");
  {
    FromResult n;
    n.rec.tracersSky[5] = kNaN;
    n.rec.outsideRedshiftCut[1] = 1;
    fails(n, "findVoidsLightcone: the reconstruction's tracersSky holds a non-finite value for object 1 (from 0), "
             "component 2", "NaN on a cut row");
    FromResult d;
    d.rec.tracersSky[7] = -1.75;
    fails(d, "findVoidsLightcone: the reconstruction's tracersSky holds a declination of -1.75 radians for object 2 "
             "(from 0); it must lie in [-pi/2, pi/2]", "declination");
    FromResult p;
    p.rec.mpsProfile.mps[0] = 0.;
    fails(p, "the mps profile holds a mean particle separation of 0 at entry 0; it must be positive", "profile");
    FromResult e;
    e.randomsSky.clear();
    fails(e, "findVoidsLightcone: randomsSky is empty", "no randoms");
    FromResult s;
    s.randomsSky.pop_back();
    fails(s, "findVoidsLightcone: randomsSky holds 23 entries; it must hold 3 per object", "randoms size");
    FromResult r;
    r.randomsSky[4] = 2.;
    fails(r, "the declination of random 1 is 2 radians; it must lie in [-pi/2, pi/2]", "random declination");
  }

  // The redshift cut a reconstruction records must have no NaN bound and its lower bound not above the upper one.
  group("a reconstruction whose redshift cut is NaN or reversed is refused as malformed");
  {
    FromResult n;
    n.rec.selection.redshiftCut = otswap::RedshiftCut{kNaN, 1.};
    fails(n, "findVoidsLightcone: the reconstruction is malformed: its redshift cut (otswap::Result::selection) "
             "is [nan, 1]", "NaN bound");
    FromResult r;
    r.rec.selection.redshiftCut = otswap::RedshiftCut{0.5, 0.25};
    fails(r, "its redshift cut (otswap::Result::selection) is [0.5, 0.25]",
          "reversed");
  }

  // The mask given must leave out every tracer the reconstruction left out as unobserved.
  group("a mask that allows a tracer the reconstruction masked is refused");
  {
    FromResult r;
    r.rec.outsideMask[1] = 1;
    fails(r, "findVoidsLightcone: the mask given allows tracer 1, which the reconstruction excluded as outside its "
             "mask (Result::outsideMask); the mask must not be more permissive than the reconstruction's", "mask");
  }

  // The configuration and the rows the search reads are checked as in the array overloads, with the reconstruction's table.
  group("the configuration and the kept rows are checked as in the array overloads");
  {
    FromResult c;
    c.config.divergenceCellSize = -1.;
    fails(c, "the divergence cell size is", "cell size");
    FromResult z;
    z.rec.tracersSky[5] = 2.5;
    fails(z, "the redshift of tracer 1 is 2.5, outside the distance table [0, 2]", "tracer redshift");
    FromResult t;
    t.rec.tracers[6] = kNaN;
    fails(t, "findVoidsLightcone: the reconstruction's tracers hold a non-finite value for object 2 (from 0), "
             "component 0", "tracers NaN");
    FromResult d;
    d.rec.displacement[10] = kInf;
    fails(d, "the displacements hold a non-finite value at entry 10 (realization 0, tracer 3, component 1)",
          "displacement");
    FromResult k;
    k.rec.outsideRedshiftCut = {1, 1, 0, 0};
    k.rec.outsideMask = {0, 0, 1, 1};
    k.mask = otswap::Mask(std::vector<double>(12, 0.), otswap::PixelOrdering::Ring);
    fails(k, "no tracer is left to read: of 4 tracers, 4 lie on unobserved pixels and 0 outside the redshift cut",
          "no tracer kept");
    FromResult v;
    v.rec.outsideRedshiftCut = {1, 1, 1, 0};
    v.rec.valid[3] = 0;
    fails(v, "no valid displacement is left to read: of 3 valid displacements", "no displacement read");
    FromResult r;
    r.rec.selection.redshiftCut = otswap::RedshiftCut{0.5, 1.};
    fails(r, "no random is left: of 8 randoms, 0 lie on unobserved pixels and 8 outside the redshift cut",
          "no random kept");
  }

  return report("test_lightcone_entry");
}
