# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Tests of the Python bindings of the lightcone entry points: the angle
unit of the arrays and of the results, the forms of each function, the
configuration classes, the pipeline and the errors the bindings raise.

The bit-for-bit tests compare the arrays, scalars and report lines of the
bindings with those of the C++ library on the same inputs, written by
tests/python_reference.cpp, as the box tests do; they are skipped without
BITVF_PYTHON_REFERENCE."""

from pathlib import Path

import numpy as np
import pytest

import otswap
import bitvf
from test_bitvf import (assert_same_arrays, assert_same_grid, assert_same_report, assert_same_result,
                        config_params, kept as box_kept, needs_cpp, run_cpp, same_double, scalars)

DATA = Path(__file__).parent / "executable" / "data"


def quiet(**fields):
    return bitvf.Config(verbosity="silent", **fields)


def kept(**fields):
    return quiet(keep_divergence=True, keep_assignment=True, **fields)


# --------------------------------------------------------------------------
# inputs

@pytest.fixture(scope="module")
def sky():
    tracers = np.loadtxt(DATA / "lightcone_tracers.dat")[:, :3]
    randoms = np.loadtxt(DATA / "lightcone_randoms.dat")[:, :3]
    return tracers, randoms


@pytest.fixture(scope="module")
def mask():
    return otswap.Mask(str(DATA / "lightcone_mask.fits"))


@pytest.fixture(scope="module")
def distances():
    return otswap.DistanceTable.flat(0.3, 0.7)


@pytest.fixture(scope="module")
def reconstruction(sky, mask, distances):
    tracers, randoms = sky
    return otswap.reconstruct_lightcone(
        tracers, randoms, mask=mask, n_bins=1, distances=distances, angle_unit="deg",
        config=otswap.Config(n_realizations=2, seed=1, verbosity="silent"))


@pytest.fixture(scope="module")
def result(reconstruction, sky, mask):
    return bitvf.find_voids_lightcone(reconstruction, sky[1], mask=mask, angle_unit="deg",
                                      config=kept(divergence_cell_size=0.4))


def rad(a):
    out = np.array(a, dtype=float)
    out[:, :2] = np.deg2rad(out[:, :2])
    return out


def error(call, text):
    with pytest.raises(bitvf.Error) as caught:
        call()
    assert text in str(caught.value)


# --------------------------------------------------------------------------
# the angle unit of the results follows the call that made them

def test_a_result_reports_the_unit_of_its_call(reconstruction, sky, mask, result):
    config = kept(divergence_cell_size=0.4)
    in_rad = bitvf.find_voids_lightcone(reconstruction, rad(sky[1]), mask=mask, angle_unit="rad", config=config)
    assert result.geometry == "lightcone" and result.angle_unit == "deg"
    assert in_rad.angle_unit == "rad"
    assert result.divergence.angle_unit == "deg" and in_rad.divergence.angle_unit == "rad"
    assert np.array_equal(result.voids.centre, in_rad.voids.centre)
    assert np.array_equal(result.voids.centre_sky[:, 2], in_rad.voids.centre_sky[:, 2])
    assert np.array_equal(result.voids.centre_sky[:, :2], np.rad2deg(in_rad.voids.centre_sky[:, :2]))


def test_the_grid_angles_follow_the_unit(result, reconstruction, sky, mask):
    in_rad = bitvf.find_voids_lightcone(reconstruction, rad(sky[1]), mask=mask, angle_unit="rad",
                                        config=kept(divergence_cell_size=0.4))
    g, h = result.divergence, in_rad.divergence
    for name in ("ra_min", "ra_max", "dec_min", "dec_max"):
        assert getattr(g, name) == getattr(h, name) * (180.0 / np.pi)
    for name in ("cell_ra_min", "cell_ra_max", "cell_dec_min", "cell_dec_max"):
        assert np.array_equal(getattr(g, name), getattr(h, name) * (180.0 / np.pi))
    sky_g, sky_h = g.centre_sky, h.centre_sky
    assert np.array_equal(sky_g[:, 2], sky_h[:, 2])
    assert np.array_equal(sky_g[:, :2], np.rad2deg(sky_h[:, :2]))
    assert np.all(sky_g[:, 0] >= 0.0) and np.all(sky_g[:, 0] < 360.0)
    assert g.centre_sky is sky_g
    for a in (sky_g, g.cell_ra_min, h.centre_sky):
        assert not a.flags.writeable


def test_lightcone_grid_arrays_are_flat(result):
    g = result.divergence
    n = g.n_cells
    assert g.geometry == "lightcone" and n > 1
    for name in ("divergence", "volume", "weight", "cell_ra_min", "cell_distance_max", "void_index",
                 "overdense_index", "n_tracers", "n_displacements", "masked"):
        assert getattr(g, name).shape == (n,), name
    assert g.centre.shape == (n, 3) and g.centre_sky.shape == (n, 3)
    assert g.masked.dtype == bool
    assert result.voids.centre_sky.shape == (result.voids.n_basins, 3)


def test_the_parts_of_a_pipeline_result_follow_the_unit(sky, mask, distances):
    for unit, factor in (("deg", 1.0), ("rad", np.pi / 180.0)):
        tracers, randoms = (sky[0].copy(), sky[1].copy())
        if unit == "rad":
            tracers, randoms = rad(tracers), rad(randoms)
        run = bitvf.run_lightcone(tracers, randoms, mask=mask, distances=distances, angle_unit=unit,
                                  config=bitvf.LightconeConfig(
                                      ot=otswap.Config(n_realizations=2, seed=1, verbosity="silent"), n_bins=1,
                                      voids=kept(divergence_cell_size=0.4)))
        assert run.voids.angle_unit == unit and run.voids.divergence.angle_unit == unit
        assert abs(run.voids.divergence.ra_min) < 360.0 * factor


def test_a_box_result_has_no_angle_unit():
    rng = np.random.default_rng(3)
    tracers = rng.uniform(0.0, 50.0, (300, 3))
    d = rng.uniform(-1.0, 1.0, (1, 300, 3))
    r = bitvf.find_voids_box(tracers, d, mps=3.0, config=kept())
    assert r.angle_unit is None and r.divergence.angle_unit is None
    assert r.counts.tracers_outside_mask == 0 and r.counts.randoms_outside_redshift_cut == 0


# --------------------------------------------------------------------------
# the forms of find_voids_lightcone and lightcone_grid

def test_arrays_give_the_search_of_the_reconstruction(reconstruction, sky, mask, distances, result):
    tracers, randoms = sky
    config = kept(divergence_cell_size=0.4)
    a = bitvf.find_voids_lightcone(tracers, reconstruction.displacement, randoms, reconstruction.valid, mask=mask,
                                   distances=distances, mps_profile=reconstruction.mps_profile, angle_unit="deg",
                                   config=config)
    assert a.n_tracers == result.n_tracers and a.n_realizations == result.n_realizations
    assert np.array_equal(a.voids.centre, result.voids.centre)
    assert np.array_equal(a.divergence.divergence, result.divergence.divergence)
    assert a.counts.randoms == result.counts.randoms
    assert a.redshift_cut.min == -np.inf and a.mask_nside == mask.nside
    cut = (0.26, 0.4)
    c = bitvf.find_voids_lightcone(tracers, reconstruction.displacement, randoms, mask=mask, distances=distances,
                                   mps_profile=reconstruction.mps_profile, angle_unit="deg", config=config,
                                   redshift_cut=cut, ids=np.arange(500)[::-1])
    assert (c.redshift_cut.min, c.redshift_cut.max) == cut
    assert c.counts.tracers_outside_redshift_cut > 0
    assert not c.assignment.in_cut.all() and c.assignment.in_mask.dtype == bool
    assert np.array_equal(c.ids, np.arange(500)[::-1])


def test_cartesian_positions(reconstruction, sky, mask, distances, result):
    tracers = otswap.to_cartesian(sky[0], distances, angle_unit="deg")
    randoms = otswap.to_cartesian(sky[1], distances, angle_unit="deg")
    c = bitvf.find_voids_lightcone_from_cartesian(
        tracers, reconstruction.displacement, randoms, mask=mask, distances=distances,
        mps_profile=reconstruction.mps_profile, angle_unit="deg", config=kept(divergence_cell_size=0.4))
    assert c.geometry == "lightcone" and c.angle_unit == "deg" and c.n_tracers == 500
    assert c.divergence.n_cells == result.divergence.n_cells


def test_a_reused_grid_may_come_back_in_another_unit(reconstruction, sky, mask, result):
    again = bitvf.find_voids_lightcone(reconstruction, rad(sky[1]), mask=mask, angle_unit="rad",
                                       config=kept(), grid=result.divergence)
    assert again.angle_unit == "rad" and again.divergence.angle_unit == "rad"
    assert again.divergence.n_cells == result.divergence.n_cells
    assert np.array_equal(again.voids.centre, result.voids.centre)


def test_lightcone_grid_builds_the_grid_of_the_search(reconstruction, sky, mask, distances, result):
    g = bitvf.lightcone_grid(reconstruction, sky[1], mask=mask, angle_unit="deg", cell_size=0.4)
    h = bitvf.lightcone_grid(sky[0], sky[1], mask=mask, distances=distances, mps_profile=reconstruction.mps_profile,
                             cell_size=0.4, angle_unit="deg")
    r = result.divergence
    for grid in (g, h):
        assert grid.geometry == "lightcone" and grid.angle_unit == "deg"
        assert grid.n_cells == r.n_cells and grid.cell_size_mps == r.cell_size_mps
        assert (grid.distance_min, grid.distance_max, grid.ra_min, grid.ra_max, grid.dec_min, grid.dec_max) == \
            (r.distance_min, r.distance_max, r.ra_min, r.ra_max, r.dec_min, r.dec_max)
        assert grid.divergence.shape == (0,) and grid.centre_sky.shape == (0, 3)
    default = bitvf.lightcone_grid(reconstruction, rad(sky[1]), mask=mask, angle_unit="rad")
    assert default.angle_unit == "rad" and default.cell_size_mps == min(1.0, (20.0 / 2) ** (1.0 / 3.0))
    shared = bitvf.find_voids_lightcone(reconstruction, sky[1], mask=mask, angle_unit="deg", config=quiet(), grid=g)
    assert shared.voids.n_basins == result.voids.n_basins


# --------------------------------------------------------------------------
# the configuration classes

def test_lightcone_rsd_defaults_repr_and_equality():
    bias = otswap.BiasTable([0.0, 2.0], [1.0, 1.5])
    rsd = bitvf.LightconeRsd(bias=bias, sigma=10.0)
    assert rsd.second_pass_randoms == bitvf.SecondPassRandoms.RESAMPLE and rsd.zspace_voids is False
    assert rsd.sigma == 10.0 and np.array_equal(rsd.bias.bias, [1.0, 1.5])
    assert rsd == bitvf.LightconeRsd(bias=otswap.BiasTable([0.0, 2.0], [1.0, 1.5]), sigma=10.0)
    assert rsd != bitvf.LightconeRsd(bias=bias, sigma=11.0)
    assert "LightconeRsd(bias=" in repr(rsd) and "SecondPassRandoms.RESAMPLE" in repr(rsd)
    with pytest.raises(TypeError):
        bitvf.LightconeRsd(sigma=1.0)
    with pytest.raises(TypeError):
        bitvf.LightconeRsd(bias=bias)
    rsd.second_pass_randoms = bitvf.SecondPassRandoms.EXTERNAL
    rsd.zspace_voids = True
    assert rsd.second_pass_randoms == bitvf.SecondPassRandoms.EXTERNAL and rsd.zspace_voids is True
    error(lambda: setattr(rsd, "second_pass_randoms", "none"),
          "second_pass_randoms must be a bitvf.SecondPassRandoms; got 'none'")
    error(lambda: setattr(rsd, "sigma", "wide"), "sigma must be a real number; got 'wide'")
    error(lambda: bitvf.LightconeRsd(bias=1.5, sigma=1.0), "bias must be an otswap.BiasTable; got 1.5")


def test_lightcone_config_defaults_references_and_equality():
    cfg = bitvf.LightconeConfig(n_bins=5)
    assert cfg.n_bins == 5 and cfg.rsd is None and cfg.keep_reconstructions is False
    assert cfg.cut.min == -np.inf and cfg.cut.max == np.inf
    assert cfg == bitvf.LightconeConfig(n_bins=5) and cfg != bitvf.LightconeConfig(n_bins=6)
    with pytest.raises(TypeError):
        bitvf.LightconeConfig()
    cfg.voids.keep_divergence = True
    assert cfg.voids.keep_divergence is True
    cfg.ot.seed = 99
    assert cfg.ot.seed == 99
    cfg.cut = (0.2, 0.5)
    assert (cfg.cut.min, cfg.cut.max) == (0.2, 0.5)
    cfg.cut.max = 0.6
    assert cfg.cut.max == 0.6
    cfg.cut = None
    assert cfg.cut.max == np.inf
    cfg.rsd = bitvf.LightconeRsd(bias=otswap.BiasTable([0.0, 2.0], [1.0, 1.5]), sigma=5.0)
    cfg.rsd.sigma = 7.0
    assert cfg.rsd.sigma == 7.0
    assert cfg != bitvf.LightconeConfig(n_bins=5)
    assert repr(cfg).startswith("LightconeConfig(ot=") and "n_bins=5" in repr(cfg)
    error(lambda: setattr(cfg, "n_bins", -1), "n_bins is -1; it must be a non-negative integer")
    error(lambda: setattr(cfg, "rsd", "yes"), "rsd must be a bitvf.LightconeRsd or None; got 'yes'")
    error(lambda: setattr(cfg, "cut", (0.1, 0.2, 0.3)), "cut must be an otswap.RedshiftCut or a (min, max) pair")


# --------------------------------------------------------------------------
# run_lightcone

def lightcone_config(**fields):
    return bitvf.LightconeConfig(ot=otswap.Config(n_realizations=2, seed=1, verbosity="silent"), n_bins=1,
                                 voids=kept(divergence_cell_size=0.4), **fields)


def test_run_lightcone_without_correction(sky, mask, distances):
    run = bitvf.run_lightcone(sky[0], sky[1], mask=mask, distances=distances, angle_unit="deg",
                              config=lightcone_config(keep_reconstructions=True), ids=np.arange(500) + 7)
    assert run.voids.geometry == "lightcone" and run.zspace_voids is None
    assert run.reconstruction.geometry == "lightcone" and run.real_space is None
    assert np.array_equal(run.voids.ids, np.arange(500) + 7)
    assert run.voids.divergence.n_cells > 1
    lean = bitvf.run_lightcone(sky[0], sky[1], mask=mask, distances=distances, angle_unit="deg",
                               config=lightcone_config())
    assert lean.reconstruction is None and np.array_equal(lean.voids.voids.centre, run.voids.voids.centre)


@pytest.mark.parametrize("mode", [bitvf.SecondPassRandoms.RESAMPLE, bitvf.SecondPassRandoms.NONE])
def test_run_lightcone_with_correction(sky, mask, distances, mode):
    rsd = bitvf.LightconeRsd(bias=otswap.BiasTable([0.0, 1.5], [1.0, 1.5]), sigma=10.0, second_pass_randoms=mode,
                             zspace_voids=True,
                             correction=otswap.CorrectionConfig(verbosity="silent"))
    run = bitvf.run_lightcone(sky[0], sky[1], mask=mask, distances=distances, angle_unit="deg",
                              config=lightcone_config(rsd=rsd, keep_reconstructions=True))
    assert run.zspace_voids is not None and run.zspace_voids.angle_unit == "deg"
    assert run.voids.angle_unit == "deg" and run.voids.n_tracers == 500
    assert run.real_space.geometry == "lightcone" and run.zspace_reconstruction.geometry == "lightcone"
    assert run.voids.uncorrected.dtype == bool


def test_run_lightcone_with_external_randoms(sky, mask, distances):
    rsd = bitvf.LightconeRsd(bias=otswap.BiasTable([0.0, 1.5], [1.0, 1.5]), sigma=10.0,
                             second_pass_randoms=bitvf.SecondPassRandoms.EXTERNAL,
                             correction=otswap.CorrectionConfig(verbosity="silent"))
    config = lightcone_config(rsd=rsd)
    error(lambda: bitvf.run_lightcone(sky[0], sky[1], mask=mask, distances=distances, angle_unit="deg",
                                      config=config), "second")
    run = bitvf.run_lightcone(sky[0], sky[1], sky[1], mask=mask, distances=distances, angle_unit="deg", config=config)
    assert run.voids.n_tracers == 500
    error(lambda: bitvf.run_lightcone(sky[0], sky[1], np.empty((0, 3)), mask=mask, distances=distances,
                                      angle_unit="deg", config=config),
          "randoms_sky_second_pass is empty; omit it unless the second pass takes external randoms")


# --------------------------------------------------------------------------
# files

def test_divergence_file_round_trip_in_both_units(result, tmp_path):
    file = tmp_path / "divergence.fits"
    bitvf.io.write_divergence(file, result)
    g = bitvf.io.read_divergence(file, angle_unit="deg")
    h = bitvf.io.read_divergence(file, angle_unit="rad")
    assert g.geometry == "lightcone" and g.angle_unit == "deg" and h.angle_unit == "rad"
    assert g.n_cells == result.divergence.n_cells
    assert np.allclose(g.ra_min, result.divergence.ra_min, rtol=1e-12)
    assert np.allclose(h.ra_min, np.deg2rad(result.divergence.ra_min), rtol=1e-12)
    assert g.centre_sky.shape == (g.n_cells, 3)
    error(lambda: bitvf.io.read_divergence(file), "angle_unit is required for a lightcone file")
    error(lambda: bitvf.io.read_divergence(file, angle_unit="turns"),
          "angle_unit must be \"deg\" or \"rad\"; got 'turns'")
    bitvf.io.write_voids(tmp_path / "voids.dat", result)
    bitvf.io.write_galaxies(tmp_path / "galaxies.dat", result)


# --------------------------------------------------------------------------
# errors

def test_argument_errors(reconstruction, sky, mask, distances):
    tracers, randoms = sky
    profile = reconstruction.mps_profile
    d, v = reconstruction.displacement, reconstruction.valid

    def search(**over):
        args = dict(mask=mask, distances=distances, mps_profile=profile, angle_unit="deg", config=quiet())
        args.update(over)
        return bitvf.find_voids_lightcone(tracers, d, randoms, v, **args)

    error(lambda: search(angle_unit="degrees"), "angle_unit must be \"deg\" or \"rad\"; got 'degrees'")
    error(lambda: search(angle_unit=None), "angle_unit must be \"deg\" or \"rad\"; got None")
    with pytest.raises(TypeError):
        bitvf.find_voids_lightcone(tracers, d, randoms, v, mask=mask, distances=distances, mps_profile=profile)
    error(lambda: search(mask=None), "mask must be an otswap.Mask; got None")
    error(lambda: search(distances=None), "find_voids_lightcone: distances is required with arrays")
    error(lambda: search(mps_profile=None), "find_voids_lightcone: mps_profile is required with arrays")
    error(lambda: search(config="x"), "config must be a bitvf.Config; got 'x'")
    error(lambda: search(redshift_cut=(0.1, 0.2, 0.3)), "redshift_cut must be an otswap.RedshiftCut or a (min, max) pair")
    error(lambda: search(grid=3), "grid must be a bitvf.DivergenceGrid; got 3")
    bad = tracers.copy()
    bad[7, 1] = 120.0
    error(lambda: bitvf.find_voids_lightcone(bad, d, randoms, v, mask=mask, distances=distances, mps_profile=profile,
                                             angle_unit="deg"),
          "tracers_sky holds a declination of 120.000000 degrees at object 7 (from 0), outside [-90, 90]")
    error(lambda: bitvf.find_voids_lightcone(tracers, d[:, :10], randoms, mask=mask, distances=distances,
                                             mps_profile=profile, angle_unit="deg"),
          "displacement must have shape (n_realizations, 500, 3), 500 the rows of the tracers")
    error(lambda: bitvf.find_voids_lightcone(tracers, d, None, mask=mask, distances=distances, mps_profile=profile,
                                             angle_unit="deg"), "find_voids_lightcone: randoms_sky is required")
    error(lambda: bitvf.find_voids_lightcone(tracers, d, randoms, v, mask=mask, distances=distances,
                                             mps_profile=profile, angle_unit="deg", ids=[1, 2]), "ids")


def test_form_errors(reconstruction, sky, mask, distances):
    randoms = sky[1]
    error(lambda: bitvf.find_voids_lightcone(reconstruction, randoms, mask=mask, angle_unit="deg", distances=distances,
                                             ids=[0]),
          "find_voids_lightcone: distances and ids do not apply to a reconstruction (otswap.Result)")
    error(lambda: bitvf.find_voids_lightcone(reconstruction, mask=mask, angle_unit="deg"),
          "find_voids_lightcone: randoms_sky is required")
    error(lambda: bitvf.find_voids_lightcone(reconstruction, reconstruction.displacement, randoms, mask=mask,
                                             angle_unit="deg"),
          "displacement does not apply to a reconstruction")
    error(lambda: bitvf.lightcone_grid(reconstruction, randoms, mask=mask, angle_unit="deg", distances=distances),
          "lightcone_grid: distances does not apply to a reconstruction (otswap.Result)")
    error(lambda: bitvf.lightcone_grid(sky[0], randoms, mask=mask, angle_unit="deg", distances=distances,
                                       mps_profile=reconstruction.mps_profile),
          "lightcone_grid: cell_size is required with arrays")
    error(lambda: bitvf.lightcone_grid(sky[0], randoms, mask=mask, angle_unit="deg", cell_size=0.4),
          "lightcone_grid: distances is required with arrays")
    error(lambda: bitvf.lightcone_grid(reconstruction, randoms, mask=mask, angle_unit="deg", cell_size="big"),
          "cell_size must be a real number; got 'big'")


def test_pipeline_errors(sky, mask, distances):
    config = lightcone_config()
    error(lambda: bitvf.run_lightcone(sky[0], sky[1], mask=mask, distances=distances, angle_unit="deg", config=None),
          "config must be a bitvf.LightconeConfig; got None")
    error(lambda: bitvf.run_lightcone(sky[0], sky[1], mask=mask, distances=None, angle_unit="deg", config=config),
          "distances must be an otswap.DistanceTable; got None")
    error(lambda: bitvf.run_lightcone(sky[0], sky[1], mask=mask, distances=distances, angle_unit="deg",
                                      config=bitvf.LightconeConfig(n_bins=0)), "bins")
    with pytest.raises(TypeError):
        bitvf.run_lightcone(sky[0], sky[1], mask=mask, distances=distances, angle_unit="deg")


# --------------------------------------------------------------------------
# bit identity with the C++ library, and the same report

FITS = DATA / "lightcone_mask.fits"
UNITS = ("deg", "rad")
CUT = (0.26, 0.4)


def lc_config(verbosity="normal"):
    return box_kept(verbosity, divergence_cell_size=0.4)


def in_unit(a, unit):
    return rad(a) if unit == "rad" else np.array(a, dtype=float)


def assert_same_lightcone_grid(g, out, prefix, values):
    assert g.geometry == "lightcone"
    for name in ("distance_min", "distance_max", "ra_min", "ra_max", "dec_min", "dec_max"):
        assert same_double(getattr(g, name), values[f"{prefix}.{name}"]), name
    assert_same_arrays({name: getattr(g, name) for name in
                        ("cell_ra_min", "cell_ra_max", "cell_dec_min", "cell_dec_max", "cell_distance_min",
                         "cell_distance_max", "centre_sky", "masked")}, out, prefix + ".")
    profile = g.mps_profile
    assert_same_arrays({"redshift": profile.redshift, "mps": profile.mps, "count": profile.count}, out,
                       prefix + ".mps_profile.")


def assert_same_lightcone(r, out, prefix=""):
    assert_same_result(r, out, prefix)
    p = prefix + "." if prefix else ""
    values = scalars(out)
    for name in ("tracers_outside_cut", "tracers_outside_mask", "tracers_outside_redshift_cut", "randoms",
                 "randoms_outside_mask", "randoms_outside_redshift_cut", "tracers_in_masked_cells",
                 "randoms_outside_grid", "displacements_rejected"):
        assert getattr(r.counts, name) == values[f"{p}counts.{name}"], name
    assert r.mask_nside == values[p + "mask_nside"]
    assert same_double(r.redshift_cut.min, values[p + "redshift_cut_min"])
    assert same_double(r.redshift_cut.max, values[p + "redshift_cut_max"])
    assert_same_arrays({"redshift": r.mps_profile.redshift, "mps": r.mps_profile.mps,
                        "count": r.mps_profile.count}, out, p + "mps_profile.")
    assert_same_arrays({"in_cut": r.assignment.in_cut}, out, p + "assignment.")
    assert_same_lightcone_grid(r.divergence, out, p + "grid", values)


def lc_inputs(tmp_path, sky, unit, **extra):
    tracers, randoms = (in_unit(sky[0], unit), in_unit(sky[1], unit))
    return {"tracers": tracers, "randoms": randoms, **extra}


def complete_file(tmp_path, reconstruction):
    complete = tmp_path / "complete.fits"
    otswap.io.write_displacement_complete(complete, reconstruction)
    return complete


@needs_cpp
@pytest.mark.parametrize("unit", UNITS)
@pytest.mark.parametrize("verbosity", ("normal", "detailed"))
def test_find_voids_lightcone_arrays(tmp_path, capsys, reconstruction, sky, mask, distances, unit, verbosity):
    config = lc_config(verbosity)
    ids = np.arange(500, dtype=np.int64)[::-1] + 10
    tracers, randoms = (in_unit(sky[0], unit), in_unit(sky[1], unit))
    r = bitvf.find_voids_lightcone(tracers, reconstruction.displacement, randoms, reconstruction.valid, mask=mask,
                                   distances=distances, mps_profile=reconstruction.mps_profile, angle_unit=unit,
                                   config=config, redshift_cut=CUT, ids=ids)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "lightcone", {"tracers": tracers, "randoms": randoms}, {"ids": ids}, form="arrays",
                  complete=complete_file(tmp_path, reconstruction), mask=FITS, omega_m=0.3, h=0.7, angles=unit,
                  cut_min=CUT[0], cut_max=CUT[1], **config_params(config))
    assert_same_lightcone(r, out)
    assert_same_report(printed, out)
    assert r.counts.tracers_outside_redshift_cut > 0


@needs_cpp
@pytest.mark.parametrize("unit", UNITS)
def test_find_voids_lightcone_cartesian(tmp_path, capsys, reconstruction, sky, mask, distances, unit):
    config = lc_config()
    tracers = otswap.to_cartesian(in_unit(sky[0], unit), distances, angle_unit=unit)
    randoms = otswap.to_cartesian(in_unit(sky[1], unit), distances, angle_unit=unit)
    r = bitvf.find_voids_lightcone_from_cartesian(
        tracers, reconstruction.displacement, randoms, reconstruction.valid, mask=mask, distances=distances,
        mps_profile=reconstruction.mps_profile, angle_unit=unit, config=config, redshift_cut=CUT)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "lightcone", {"tracers": tracers, "randoms": randoms}, form="cartesian",
                  complete=complete_file(tmp_path, reconstruction), mask=FITS, omega_m=0.3, h=0.7, angles=unit,
                  cut_min=CUT[0], cut_max=CUT[1], **config_params(config))
    assert_same_lightcone(r, out)
    assert_same_report(printed, out)


@needs_cpp
@pytest.mark.parametrize("unit", UNITS)
def test_find_voids_lightcone_on_a_reconstruction(tmp_path, capsys, reconstruction, sky, mask, unit):
    complete = complete_file(tmp_path, reconstruction)
    rec = otswap.io.read_displacement_complete(complete, angle_unit="rad")
    config = lc_config()
    randoms = in_unit(sky[1], unit)
    r = bitvf.find_voids_lightcone(rec, randoms, mask=mask, angle_unit=unit, config=config)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "lightcone", {"randoms": randoms}, form="reconstruction", complete=complete,
                  mask=FITS, omega_m=0.3, h=0.7, angles=unit, **config_params(config))
    assert_same_lightcone(r, out)
    assert_same_report(printed, out)


@needs_cpp
@pytest.mark.parametrize("unit", UNITS)
def test_lightcone_grid_is_the_library_grid(tmp_path, reconstruction, sky, mask, distances, unit):
    complete = complete_file(tmp_path, reconstruction)
    rec = otswap.io.read_displacement_complete(complete, angle_unit="rad")
    tracers, randoms = (in_unit(sky[0], unit), in_unit(sky[1], unit))
    common = dict(complete=complete, mask=FITS, omega_m=0.3, h=0.7, angles=unit)
    for cell in (None, 0.4):
        g = bitvf.lightcone_grid(rec, randoms, mask=mask, angle_unit=unit, cell_size=cell)
        out = run_cpp(tmp_path, "lightcone_grid", {"randoms": randoms}, form="reconstruction", name=f"rec{cell}",
                      **common, **({} if cell is None else {"cell": cell}))
        assert_same_grid(g, out, "grid", scalars(out))
        assert_same_lightcone_grid(g, out, "grid", scalars(out))
        assert g.divergence.shape == (0,) and g.centre_sky.shape == (0, 3)
    g = bitvf.lightcone_grid(tracers, randoms, mask=mask, distances=distances,
                             mps_profile=reconstruction.mps_profile, cell_size=0.4, angle_unit=unit,
                             redshift_cut=CUT)
    out = run_cpp(tmp_path, "lightcone_grid", {"tracers": tracers, "randoms": randoms}, form="arrays",
                  cut_min=CUT[0], cut_max=CUT[1], cell=0.4, **common)
    assert_same_grid(g, out, "grid", scalars(out))
    assert_same_lightcone_grid(g, out, "grid", scalars(out))


def rsd_for(mode, verbosity, zspace):
    return bitvf.LightconeRsd(bias=otswap.BiasTable([0.0, 1.5], [1.0, 1.5]), sigma=10.0, second_pass_randoms=mode,
                              zspace_voids=zspace, correction=otswap.CorrectionConfig(verbosity=verbosity))


@needs_cpp
@pytest.mark.parametrize("unit", UNITS)
@pytest.mark.parametrize("verbosity", ("normal", "detailed"))
def test_run_lightcone_without_correction_is_the_library_run(tmp_path, capsys, sky, mask, distances, unit,
                                                             verbosity):
    config = bitvf.LightconeConfig(ot=otswap.Config(n_realizations=2, seed=7, verbosity=verbosity), n_bins=1,
                                   cut=CUT, voids=lc_config(verbosity), keep_reconstructions=True)
    ids = np.arange(500, dtype=np.int64) + 3
    tracers, randoms = (in_unit(sky[0], unit), in_unit(sky[1], unit))
    r = bitvf.run_lightcone(tracers, randoms, mask=mask, distances=distances, angle_unit=unit, config=config,
                            ids=ids)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "run_lightcone", {"tracers": tracers, "randoms": randoms}, {"ids": ids}, mask=FITS,
                  omega_m=0.3, h=0.7, angles=unit, n_realizations=2, seed=7, n_bins=1, cut_min=CUT[0],
                  cut_max=CUT[1], keep_reconstructions=1, **config_params(config.voids))
    assert_same_lightcone(r.voids, out, "voids")
    assert_same_report(printed, out)
    assert r.zspace_voids is None
    assert np.asarray(r.reconstruction.displacement).tobytes() == \
        Path(f"{out}.reconstruction.displacement").read_bytes()


@needs_cpp
@pytest.mark.parametrize("mode", ["none", "resample", "external"])
@pytest.mark.parametrize("unit", UNITS)
def test_run_lightcone_with_correction_is_the_library_run(tmp_path, capsys, sky, mask, distances, mode, unit):
    verbosity = "normal"
    modes = {"none": bitvf.SecondPassRandoms.NONE, "resample": bitvf.SecondPassRandoms.RESAMPLE,
             "external": bitvf.SecondPassRandoms.EXTERNAL}
    config = bitvf.LightconeConfig(ot=otswap.Config(n_realizations=2, seed=11, verbosity=verbosity), n_bins=1,
                                   cut=CUT, rsd=rsd_for(modes[mode], verbosity, True), voids=lc_config(verbosity),
                                   keep_reconstructions=True)
    tracers, randoms = (in_unit(sky[0], unit), in_unit(sky[1], unit))
    second = np.ascontiguousarray(randoms[::-1])
    arrays = {"tracers": tracers, "randoms": randoms}
    if mode == "external":
        arrays["randoms2"] = second
    r = bitvf.run_lightcone(tracers, randoms, second if mode == "external" else None, mask=mask,
                            distances=distances, angle_unit=unit, config=config)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "run_lightcone", arrays, mask=FITS, omega_m=0.3, h=0.7, angles=unit,
                  n_realizations=2, seed=11, n_bins=1, cut_min=CUT[0], cut_max=CUT[1],
                  bias=DATA / "lightcone_bias.dat", rsd_sigma=10.0, second=mode, zspace=1, keep_reconstructions=1,
                  **config_params(config.voids))
    assert_same_lightcone(r.voids, out, "voids")
    assert_same_lightcone(r.zspace_voids, out, "zspace")
    assert_same_report(printed, out)
    assert r.voids.uncorrected.dtype == bool and r.voids.outside_cut.dtype == bool
    assert np.asarray(r.reconstruction.displacement).tobytes() == \
        Path(f"{out}.reconstruction.displacement").read_bytes()
