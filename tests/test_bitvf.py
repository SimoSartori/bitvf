# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Tests of the Python bindings of the box entry points.

The bit-for-bit tests compare every array, scalar and report line of the
bindings with those of the C++ library on the same inputs, written by
tests/python_reference.cpp. They need the path of that executable in the
environment variable BITVF_PYTHON_REFERENCE, and are skipped without it.
The other tests check what the bindings guarantee on their own: shapes,
dtypes, ownership, the configuration classes, the errors they raise and the
live report.
"""

import _thread
import gc
import io as pyio
import os
import re
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pytest

import otswap
import bitvf

REFERENCE = os.environ.get("BITVF_PYTHON_REFERENCE")

needs_cpp = pytest.mark.skipif(
    not REFERENCE or not os.access(REFERENCE, os.X_OK),
    reason="BITVF_PYTHON_REFERENCE does not name the python_reference executable")

VERBOSITIES = ("normal", "detailed")


# --------------------------------------------------------------------------
# inputs

def translated():
    """3000 tracers in [0, 100)^3 and four beyond the box, over two
    realizations: those with x below 50 move by about 30 along x, the others
    by at most 2. Some displacements invalid; IDs negative and positive."""
    rng = np.random.default_rng(2026)
    n = 3004
    tracers = rng.uniform(0.0, 100.0, (n, 3))
    tracers[3000:] = 150.0 + np.arange(12.0).reshape(4, 3)
    moved = tracers[:, 0] < 50.0
    displacement = np.empty((2, n, 3))
    for r in range(2):
        d = rng.uniform(-2.0, 2.0, (n, 3))
        d[moved] = rng.uniform(-3.0, 3.0, (moved.sum(), 3))
        d[moved, 0] += 30.0
        displacement[r] = d
    valid = np.ones((2, n), dtype=bool)
    valid[0, ::97] = False
    ids = np.where(np.arange(n) % 2 == 0, 500 + np.arange(n), -(500 + np.arange(n))).astype(np.int64)
    box = ((-5.0, -5.0, -5.0), (135.0, 105.0, 105.0))
    return tracers, displacement, valid, ids, box


def lattice_field(n=(12, 17, 9), s=1.5, lowest=(-3.25, 10.5, 0.75)):
    """A field at the centres of an (nx, ny, nz) lattice, both kinds of
    basin; rows shuffled."""
    rng = np.random.default_rng(7)
    i, j, k = np.meshgrid(*(np.arange(m) for m in n), indexing="ij")
    centres = np.stack([lowest[0] + s * i, lowest[1] + s * j, lowest[2] + s * k], axis=-1).reshape(-1, 3)
    kk = 2.0 * np.pi / (6.0 * s)
    values = (np.cos(kk * i * s) + np.cos(kk * j * s) + np.cos(kk * k * s) + 6.0 / (n[0] * s) * i * s).ravel()
    values = values + rng.uniform(-0.01, 0.01, values.size)
    weights = rng.uniform(0.5, 2.0, values.size)
    order = rng.permutation(values.size)
    return centres[order], values[order], weights[order]


def quiet(**fields):
    return bitvf.Config(verbosity="silent", **fields)


def kept(verbosity="normal", **fields):
    return bitvf.Config(smoothing=bitvf.Smoothing(sigma=0.5), keep_divergence=True, keep_assignment=True,
                        n_realizations_mean=3, verbosity=verbosity, **fields)


@pytest.fixture(scope="module")
def reconstruction():
    rng = np.random.default_rng(11)
    tracers = rng.uniform(0.0, 80.0, (1500, 3))
    return otswap.reconstruct_box(tracers, box=((0.0, 0.0, 0.0), (80.0, 80.0, 80.0)),
                                  config=otswap.Config(n_realizations=3, seed=4321, verbosity="silent"))


# --------------------------------------------------------------------------
# the C++ reference

def run_cpp(tmp_path, mode, arrays=(), ints=(), name="out", **params):
    args = [REFERENCE, mode]
    for key, array in dict(arrays).items():
        path = tmp_path / f"{name}.{key}.in"
        np.ascontiguousarray(array, dtype=np.float64).tofile(path)
        args.append(f"{key}={path}")
    for key, array in dict(ints).items():
        path = tmp_path / f"{name}.{key}.in"
        np.ascontiguousarray(array, dtype=np.int64).tofile(path)
        args.append(f"{key}={path}")
    args += [f"{key}={float(value)!r}" if isinstance(value, float) else f"{key}={value}"
             for key, value in params.items()]
    out = tmp_path / name
    args.append(f"out={out}")
    run = subprocess.run(args, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    return out


def config_params(config):
    return dict(sigma=config.smoothing.sigma,
                unit="mps" if config.smoothing.unit == bitvf.LengthUnit.MPS else "coordinate",
                keep=1 if config.keep_divergence else 0, nrecmean=config.n_realizations_mean,
                verbosity=config.verbosity,
                **({} if config.divergence_cell_size is None else {"cell": config.divergence_cell_size}))


def scalars(out):
    values = {}
    for line in Path(f"{out}.scalars").read_text().splitlines():
        key, value = line.split(" ")
        values[key] = float.fromhex(value) if ("0x" in value or "inf" in value or "nan" in value) else int(value)
    return values


def same_double(a, b):
    return np.float64(a).tobytes() == np.float64(b).tobytes()


def catalog_arrays(c):
    return {"centre": c.centre, "centre_sky": c.centre_sky, "radius": c.radius, "radius_ext": c.radius_ext,
            "delta_av": c.delta_av, "delta_2mps": c.delta_2mps, "edge_flag": c.edge_flag,
            "contamination": c.contamination, "n_tracers": c.n_tracers, "field_centre": c.field_centre}


def grid_arrays(g):
    return {"centre": g.centre, "divergence": g.divergence, "divergence_smoothed": g.divergence_smoothed,
            "volume": g.volume, "weight": g.weight, "weight_smoothed": g.weight_smoothed,
            "n_displacements": g.n_displacements,
            "n_tracers": g.n_tracers, "void_index": g.void_index, "overdense_index": g.overdense_index}


def assert_same_arrays(arrays, out, prefix):
    for name, a in arrays.items():
        reference = Path(f"{out}.{prefix}{name}").read_bytes()
        assert np.ascontiguousarray(a).tobytes() == reference, f"{prefix}{name} differs from the library's"


def assert_same_grid(g, out, prefix, values):
    assert_same_arrays(grid_arrays(g), out, prefix + ".")
    assert g.n_cells == values[prefix + ".n_cells"]
    assert same_double(g.cell_size, values[prefix + ".cell_size"])
    for k in range(3):
        assert same_double(g.lower[k], values[f"{prefix}.lower{k}"])
        assert same_double(g.upper[k], values[f"{prefix}.upper{k}"])
        assert g.n_cells_axis[k] == values[f"{prefix}.n_cells_axis{k}"]
    assert g.n_realizations == values[prefix + ".n_realizations"]
    assert g.n_realizations_mean == values[prefix + ".n_realizations_mean"]
    assert same_double(g.cell_size_mps, values[prefix + ".cell_size_mps"])
    assert same_double(g.mps, values[prefix + ".mps"])
    assert (g.limits is not None) == bool(values[prefix + ".has_limits"])


def assert_same_result(r, out, prefix=""):
    p = prefix + "." if prefix else ""
    values = scalars(out)
    assert r.n_tracers == values[p + "n_tracers"]
    assert r.n_realizations == values[p + "n_realizations"]
    assert same_double(r.mps, values[p + "mps"])
    assert (0 if r.random_seed is None else r.random_seed + 1) == values[p + "random_seed"]
    cell = r.config.divergence_cell_size
    assert same_double(cell, values[p + "cell_size"]) if cell is not None else np.isnan(values[p + "cell_size"])
    assert (r.limits is not None) == bool(values[p + "has_limits"])
    for name in ("tracers", "tracers_used", "tracers_uncorrected", "tracers_outside_box",
                 "tracers_outside_grid", "displacements"):
        assert getattr(r.counts, name) == values[f"{p}counts.{name}"], name
    for kind in ("voids", "overdense"):
        catalog = getattr(r, kind)
        assert catalog.n_basins == values[f"{p}{kind}.n_basins"]
        assert_same_arrays(catalog_arrays(catalog), out, f"{p}{kind}.")
    if r.divergence is not None:
        assert_same_grid(r.divergence, out, p + "grid", values)
    else:
        assert not Path(f"{out}.{p}grid.divergence").exists()
    if r.assignment is not None:
        a = r.assignment
        assert_same_arrays({"cell_index": a.cell_index, "void_index": a.void_index,
                            "overdense_index": a.overdense_index, "distance": a.distance, "in_mask": a.in_mask,
                            "in_box": a.in_box, "used": a.used}, out, p + "assignment.")
    ids = np.fromfile(f"{out}.{p}ids", dtype=np.int64)
    assert np.array_equal(r.ids, ids if ids.size else np.arange(r.n_tracers))
    assert_same_arrays({"uncorrected": r.uncorrected, "outside_grid": r.outside_grid,
                        "outside_cut": r.outside_cut}, out, p)


def masked(text):
    return re.sub(r"in [0-9]+\.[0-9]{2} s", "in T s", text)


def assert_same_report(captured, out):
    assert masked(captured) == masked(Path(f"{out}.report").read_text())


# --------------------------------------------------------------------------
# bit identity with the C++ library, and the same report

@needs_cpp
@pytest.mark.parametrize("verbosity", VERBOSITIES)
def test_find_voids_box_arrays(tmp_path, capsys, verbosity):
    tracers, displacement, valid, ids, box = translated()
    config = kept(verbosity)
    r = bitvf.find_voids_box(tracers, displacement, valid, mps=7.0, box=box, config=config, ids=ids)
    printed = capsys.readouterr().out
    written = tmp_path / "products"
    written.mkdir()
    valid.astype(np.uint8).tofile(f"{tmp_path}/out.valid.in")
    out = run_cpp(tmp_path, "box", {"tracers": tracers, "displacement": displacement, "box": np.array(box)},
                  {"ids": ids}, n_realizations=2, mps=7.0, write=written,
                  valid=f"{tmp_path}/out.valid.in", **config_params(config))
    assert_same_result(r, out)
    assert_same_report(printed, out)
    # The writers of bitvf.io write the library's files, byte for byte.
    for extension in (".fits", ".dat"):
        for product, write in (("voids", bitvf.io.write_voids), ("overdense", bitvf.io.write_overdense),
                               ("divergence", bitvf.io.write_divergence),
                               ("galaxies", bitvf.io.write_galaxies)):
            mine = tmp_path / f"{product}{extension}"
            write(mine, r)
            assert mine.read_bytes() == (written / f"{product}{extension}").read_bytes(), product + extension


@needs_cpp
@pytest.mark.parametrize("verbosity", VERBOSITIES)
def test_find_voids_box_on_a_reused_grid(tmp_path, capsys, verbosity):
    tracers, displacement, _, ids, box = translated()
    first = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box, config=kept("silent"), ids=ids)
    config = kept(verbosity)
    r = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box, config=config, ids=ids,
                             grid=first.divergence)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "box", {"tracers": tracers, "displacement": displacement, "box": np.array(box)},
                  {"ids": ids}, n_realizations=2, mps=7.0, reuse=1, **config_params(config))
    assert_same_result(r, out)
    assert_same_report(printed, out)


@needs_cpp
@pytest.mark.parametrize("verbosity", VERBOSITIES)
def test_find_voids_box_on_a_reconstruction(tmp_path, capsys, reconstruction, verbosity):
    complete = tmp_path / "complete.fits"
    otswap.io.write_displacement_complete(complete, reconstruction)
    rec = otswap.io.read_displacement_complete(complete, angle_unit="rad")
    config = kept(verbosity)
    r = bitvf.find_voids_box(rec, config=config)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "reconstruction", complete=complete, **config_params(config))
    assert_same_result(r, out)
    assert_same_report(printed, out)
    # The in-memory reconstruction gives the same result as the one read back.
    again = bitvf.find_voids_box(reconstruction, config=kept("silent"))
    assert np.array_equal(again.voids.centre, r.voids.centre)


@needs_cpp
@pytest.mark.parametrize("verbosity", VERBOSITIES)
def test_find_voids_box_from_divergence(tmp_path, capsys, verbosity):
    tracers, displacement, _, ids, box = translated()
    source = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box, config=kept("silent"), ids=ids)
    file = tmp_path / "divergence.fits"
    bitvf.io.write_divergence(file, source)
    grid = bitvf.io.read_divergence(file)
    config = kept(verbosity)
    config.smoothing.sigma = 1.5
    r = bitvf.find_voids_box_from_divergence(grid, tracers, config=config, ids=ids)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "divergence", {"tracers": tracers}, {"ids": ids}, grid=file, **config_params(config))
    assert_same_result(r, out)
    assert_same_report(printed, out)
    assert r.assignment.used.shape == (0,) and r.assignment.used.dtype == np.bool_


@needs_cpp
@pytest.mark.parametrize("verbosity", VERBOSITIES)
def test_find_voids_box_from_field(tmp_path, capsys, verbosity):
    centres, values, weights = lattice_field()
    config = bitvf.Config(smoothing=bitvf.Smoothing(sigma=0.8), keep_divergence=True, verbosity=verbosity)
    r = bitvf.find_voids_box_from_field(centres, values, config=config, weights=weights, mps=1.2)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "field", {"centres": centres, "values": values, "weights": weights}, mps=1.2,
                  **config_params(config))
    assert_same_result(r, out)
    assert_same_report(printed, out)
    assert r.divergence.divergence.shape == (12, 17, 9)
    assert r.divergence.n_displacements.shape == (0,)


@needs_cpp
def test_read_field_and_search_on_a_mean_file(tmp_path, capsys, reconstruction):
    tracers, displacement, _, ids, box = translated()
    source = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box, config=kept("silent"), ids=ids)
    file = tmp_path / "divergence.dat"
    bitvf.io.write_divergence(file, source)
    f = bitvf.io.read_field(file, centres=("cellX", "cellY", "cellZ"), value="divergence", weight="weight")
    config = bitvf.Config(smoothing=bitvf.Smoothing(sigma=3.5, unit=bitvf.LengthUnit.COORDINATE),
                          keep_divergence=True)
    r = bitvf.find_voids_box_from_field(f.centres, f.values, config=config, weights=f.weights)
    out = run_cpp(tmp_path, "read_field", file=file, x="cellX", y="cellY", z="cellZ",
                  value="divergence", weight="weight", name="field", **config_params(config))
    assert_same_result(r, out)
    assert bitvf.io.read_field(file, centres=["cellX", "cellY", "cellZ"], value="divergence").weights is None

    capsys.readouterr()
    mean = tmp_path / "mean.fits"
    otswap.io.write_displacement_mean(mean, reconstruction)
    m = otswap.io.read_displacement_mean(mean, angle_unit="deg")
    assert m.n_realizations == 1 and m.n_realizations_mean == reconstruction.n_realizations
    config = kept()
    config.n_realizations_mean = m.n_realizations_mean
    r = bitvf.find_voids_box(m, config=config)
    assert r.n_realizations == 1 and r.config.n_realizations_mean == m.n_realizations_mean
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "mean", file=mean, name="mean", **config_params(config))
    assert_same_result(r, out)
    assert_same_report(printed, out)


@needs_cpp
@pytest.mark.parametrize("verbosity", VERBOSITIES)
def test_watershed_box(tmp_path, capsys, verbosity):
    tracers, displacement, _, ids, box = translated()
    source = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box, config=kept("silent"))
    field = np.asarray(source.divergence.divergence)[:, :10, 3:]
    p = bitvf.watershed_box(field, verbosity=verbosity)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "watershed", {"field": np.ascontiguousarray(field)}, nx=field.shape[0],
                  ny=field.shape[1], nz=field.shape[2], verbosity=verbosity)
    values = scalars(out)
    assert p.n_voids == values["partition.n_voids"] and p.n_overdense == values["partition.n_overdense"]
    assert p.void_index.shape == field.shape and p.overdense_index.shape == field.shape
    assert_same_arrays({name: getattr(p, name) for name in (
        "void_index", "overdense_index", "void_seed", "overdense_seed", "void_cells", "void_core_cells",
        "overdense_cells", "overdense_core_cells")}, out, "partition.")
    assert_same_report(printed, out)
    # A field in Fortran order is the same field.
    q = bitvf.watershed_box(np.asfortranarray(field), verbosity="silent")
    assert np.array_equal(q.void_index, p.void_index)


@needs_cpp
def test_box_grid(tmp_path, reconstruction):
    g = bitvf.box_grid(box=((0.0, -1.0, 2.0), (50.0, 30.0, 40.0)), mps=3.5, cell_size=0.9)
    out = run_cpp(tmp_path, "grid_box", {"box": np.array([0.0, -1.0, 2.0, 50.0, 30.0, 40.0])}, mps=3.5,
                  cell=0.9)
    assert_same_grid(g, out, "grid", scalars(out))
    complete = tmp_path / "complete.fits"
    otswap.io.write_displacement_complete(complete, reconstruction)
    rec = otswap.io.read_displacement_complete(complete, angle_unit="rad")
    for cell in (None, 0.75):
        g = bitvf.box_grid(rec, cell_size=cell)
        out = run_cpp(tmp_path, "grid_reconstruction", complete=complete, name=f"rec{cell}",
                      **({} if cell is None else {"cell": cell}))
        assert_same_grid(g, out, "grid", scalars(out))
        assert g.divergence.shape == (0,) and g.centre.shape == (0, 3)


@needs_cpp
@pytest.mark.parametrize("verbosity", VERBOSITIES)
def test_run_box(tmp_path, capsys, verbosity):
    tracers, _, _, ids, _ = translated()
    tracers = tracers[:1200]
    ids = ids[:1200]
    config = bitvf.BoxConfig(ot=otswap.Config(n_realizations=2, seed=99, verbosity=verbosity),
                             voids=kept(verbosity), keep_reconstructions=True)
    r = bitvf.run_box(tracers, config=config, ids=ids)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "run", {"tracers": tracers}, {"ids": ids}, n_realizations=2, seed=99,
                  keep_reconstructions=1, **config_params(config.voids))
    assert_same_result(r.voids, out, "voids")
    assert_same_report(printed, out)
    assert r.zspace_voids is None and r.zspace_reconstruction is None and r.real_space is None
    assert isinstance(r.reconstruction, otswap.Result)
    assert np.asarray(r.reconstruction.displacement).tobytes() == \
        Path(f"{out}.reconstruction.displacement").read_bytes()
    assert r.voids.random_seed == 99 + 1299709


@needs_cpp
@pytest.mark.parametrize("verbosity", VERBOSITIES)
def test_run_box_with_correction(tmp_path, capsys, verbosity):
    tracers, _, _, _, _ = translated()
    tracers = tracers[:1200]
    randoms = np.random.default_rng(5).uniform(0.0, 100.0, (2600, 3))
    table = otswap.DistanceTable.flat(0.31, 0.68, z_max=2.0, n_samples=20000)
    rsd = bitvf.BoxRsd(axis=2, redshift=0.5, bias=1.5, sigma=8.0,
                       correction=otswap.CorrectionConfig(verbosity=verbosity), zspace_voids=True)
    config = bitvf.BoxConfig(ot=otswap.Config(n_realizations=2, seed=17, verbosity=verbosity),
                             limits=((0.0, 0.0, 0.0), (100.0, 100.0, 100.0)), rsd=rsd, voids=kept(verbosity),
                             keep_reconstructions=True)
    r = bitvf.run_box(tracers, randoms, config=config, distances=table)
    printed = capsys.readouterr().out
    out = run_cpp(tmp_path, "run", {"tracers": tracers, "randoms": randoms,
                                    "box": np.array([0.0, 0.0, 0.0, 100.0, 100.0, 100.0])},
                  n_realizations=2, seed=17, keep_reconstructions=1, axis=2, redshift=0.5, bias=1.5,
                  rsd_sigma=8.0, zspace=1, omega_m=0.31, h=0.68, z_max=2.0, **config_params(config.voids))
    assert_same_result(r.voids, out, "voids")
    assert_same_result(r.zspace_voids, out, "zspace")
    assert_same_report(printed, out)
    assert isinstance(r.real_space, otswap.RealSpaceCatalog)
    assert isinstance(r.zspace_reconstruction, otswap.Result)
    assert r.voids.uncorrected.dtype == np.bool_ and r.voids.uncorrected.shape == (1200,)
    assert r.voids.random_seed is None


# --------------------------------------------------------------------------
# the live report

def box_search(verbosity="normal", **fields):
    tracers, displacement, _, _, box = translated()
    return lambda: bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box,
                                        config=kept(verbosity, **fields))


def test_silent_prints_nothing(capsys):
    box_search("silent")()
    bitvf.watershed_box(np.ones((3, 3, 3)) - 2.0, verbosity="silent")
    assert capsys.readouterr().out == ""


def test_reports_go_to_sys_stdout_line_by_line(capsys):
    box_search("detailed")()
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) > 3 and all(line.startswith("bitvf: ") for line in lines)
    assert lines[-1].startswith("bitvf: findVoidsBox: ")


def test_a_failing_call_keeps_its_lines_and_has_no_call_line(capsys):
    tracers, _, _, _, box = translated()
    still = np.zeros((1, len(tracers), 3))
    with pytest.raises(bitvf.Error, match="no void"):
        bitvf.find_voids_box(tracers, still, mps=7.0, box=box, config=kept("detailed"))
    lines = capsys.readouterr().out.splitlines()
    assert lines and lines[0].startswith("bitvf: ")
    assert any(line.startswith("bitvf: divergence") for line in lines)
    assert not any("findVoidsBox:" in line for line in lines)


class Stdout:
    """A sys.stdout that records what it is given and runs action on the
    write of the first newline."""

    def __init__(self, action):
        self.action = action
        self.writes = []

    def write(self, text):
        self.writes.append(text)
        if text == "\n" and len(self.writes) == 2:
            self.action()
        return len(text)

    def flush(self):
        pass


def with_stdout(monkeypatch, action, call):
    out = Stdout(action)
    monkeypatch.setattr(sys, "stdout", out)
    try:
        call()
    finally:
        monkeypatch.undo()
    return out


# An exception raised by print stops the printing; the call runs to its end
# and raises it.
def test_an_exception_from_print_is_raised_after_the_call(monkeypatch):
    def fail():
        raise ValueError("stdout is closed")

    with pytest.raises(ValueError, match="stdout is closed"):
        with_stdout(monkeypatch, fail, box_search("detailed"))


# A KeyboardInterrupt raised by print, and a pending signal whose handler
# raises it, stop the call at the line, and the call raises it.
@pytest.mark.parametrize("raised_by_print", (True, False))
def test_an_interrupt_stops_the_call(monkeypatch, raised_by_print):
    def interrupt():
        if raised_by_print:
            raise KeyboardInterrupt
        _thread.interrupt_main()

    out = Stdout(interrupt)
    monkeypatch.setattr(sys, "stdout", out)
    try:
        with pytest.raises(KeyboardInterrupt):
            box_search("detailed")()
    finally:
        monkeypatch.undo()
    assert len(out.writes) == 2, "nothing is printed after the interrupt"


def test_searches_release_the_gil():
    rng = np.random.default_rng(3)
    tracers = rng.uniform(0.0, 100.0, (40000, 3))
    displacement = rng.uniform(-2.0, 2.0, (4, 40000, 3))
    config = bitvf.Config(divergence_cell_size=0.2, verbosity="silent")
    elapsed = []

    def run():
        t0 = time.perf_counter()
        bitvf.find_voids_box(tracers, displacement, mps=2.0, config=config)
        elapsed.append(time.perf_counter() - t0)

    worker = threading.Thread(target=run)
    last = time.perf_counter()
    largest = 0.0
    worker.start()
    while worker.is_alive():
        now = time.perf_counter()
        largest = max(largest, now - last)
        last = now
    worker.join()
    assert elapsed[0] > 0.02, "the call is too short to tell"
    assert largest < 0.5 * elapsed[0]


# --------------------------------------------------------------------------
# results

@pytest.fixture(scope="module")
def result():
    tracers, displacement, valid, ids, box = translated()
    return bitvf.find_voids_box(tracers, displacement, valid, mps=7.0, box=box, config=kept("silent"), ids=ids)


# Every array is a read-only view, every field read-only, and no result
# class can be made from Python.
def test_results_are_read_only(result):
    for a in (result.voids.centre, result.divergence.divergence, result.assignment.in_box, result.ids):
        assert not a.flags.writeable
    with pytest.raises(AttributeError):
        result.mps = 1.0
    with pytest.raises(AttributeError):
        result.voids.n_basins = 3
    for cls in (bitvf.Result, bitvf.BasinCatalog, bitvf.DivergenceGrid, bitvf.TracerAssignment, bitvf.Counts,
                bitvf.Partition, bitvf.BitVFResult, bitvf.io.Field):
        with pytest.raises(TypeError):
            cls()


# A view keeps the result alive.
def test_views_keep_their_owner_alive():
    tracers, displacement, _, _, box = translated()
    r = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box, config=kept("silent"))
    centre, copy = r.voids.centre, np.array(r.voids.centre)
    divergence, cells = r.divergence.divergence, np.array(r.divergence.divergence)
    del r
    gc.collect()
    assert np.array_equal(centre, copy) and np.array_equal(divergence, cells)


# Per-cell arrays are (nx, ny, nz), C order: cell (i, j, k) of the C++ row
# (i ny + j) nz + k, at the centre lower + (i, j, k + 1/2) cell_size.
def test_per_cell_arrays_are_indexed_by_cell(result):
    g = result.divergence
    n = g.n_cells_axis
    assert g.divergence.shape == n and g.centre.shape == n + (3,) and g.void_index.shape == n
    i, j, k = 3, 5, 7
    expected = np.array(g.lower) + (np.array([i, j, k]) + 0.5) * g.cell_size
    assert np.allclose(g.centre[i, j, k], expected)
    flat = np.asarray(g.divergence).ravel()
    assert flat[(i * n[1] + j) * n[2] + k] == g.divergence[i, j, k]
    # cell_index is a flat row: unravel_index gives the cell, which holds the tracer.
    tracers = translated()[0]
    a = result.assignment
    for row in (0, 10, 2999):
        cell = np.unravel_index(a.cell_index[row], n)
        assert np.all(np.abs(tracers[row] - g.centre[cell]) <= 0.5 * g.cell_size + 1e-9)
        assert a.void_index[row] == g.void_index[cell]


def test_dtypes_and_shapes(result):
    v, g, a = result.voids, result.divergence, result.assignment
    n = result.n_tracers
    assert v.centre.shape == (v.n_basins, 3) and v.centre.dtype == np.float64
    assert v.centre_sky.shape == (0, 3)
    assert v.edge_flag.dtype == np.uint8 and v.n_tracers.dtype == np.uint64
    assert v.field_centre.shape == (0,)
    for flag in (a.in_box, a.used):
        assert flag.dtype == np.bool_ and flag.shape == (n,)
    assert a.in_mask.dtype == np.bool_ and a.in_mask.shape == (0,)
    assert result.counts.tracers_in_masked_cells == 0
    assert a.in_box.sum() == n - 4
    for flag in (result.uncorrected, result.outside_grid, result.outside_cut, g.masked):
        assert flag.dtype == np.bool_ and flag.shape == (0,)
    assert g.n_displacements.dtype == np.uint64 and g.void_index.dtype == np.int64
    assert g.geometry == "box" and g.angle_unit is None and g.mps_profile is None
    assert g.cell_ra_min.shape == (0,) and g.centre_sky.shape == (0, 3)
    assert result.geometry == "box" and result.angle_unit is None
    assert result.mps_profile is None and result.distances is None and result.random_seed is None
    assert isinstance(result.limits, otswap.BoxLimits)
    assert result.limits == otswap.BoxLimits((-5.0, -5.0, -5.0), (135.0, 105.0, 105.0))


def test_ids(result):
    assert result.ids.dtype == np.int64 and result.ids[1] == -501
    assert result.ids is not result.ids
    tracers, displacement, _, _, box = translated()
    plain = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box, config=quiet())
    assert np.array_equal(plain.ids, np.arange(len(tracers)))
    same = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=box, config=quiet(),
                                ids=np.arange(len(tracers), dtype=float))
    assert np.array_equal(same.ids, plain.ids)


# Result.config is a copy, divergence_cell_size the value used.
def test_result_config_is_a_copy(result):
    c = result.config
    assert c.divergence_cell_size == 1.0
    c.keep_divergence = False
    assert result.config.keep_divergence


# --------------------------------------------------------------------------
# configuration

def test_config_defaults_repr_and_equality():
    c = bitvf.Config()
    assert c.smoothing == bitvf.Smoothing(sigma=1.0, unit=bitvf.LengthUnit.MPS)
    assert c.divergence_cell_size is None and c.n_realizations_mean == 1 and c.verbosity == "normal"
    assert eval(repr(c), vars(bitvf)) == c
    d = bitvf.Config(divergence_cell_size=0.5, verbosity="detailed")
    assert d != c and eval(repr(d), vars(bitvf)) == d
    assert bitvf.BoxConfig() == bitvf.BoxConfig()
    rsd = bitvf.BoxRsd(axis=None, redshift=0.5, bias=1.2, sigma=10.0)
    assert rsd.axis is None and rsd.correction == otswap.CorrectionConfig()
    assert "BoxRsd(axis=None" in repr(rsd)


# Nested configurations are references into their holder; assigning one
# copies it in.
def test_nested_configurations_are_references():
    c = bitvf.BoxConfig()
    c.voids.keep_divergence = True
    c.voids.smoothing.sigma = 2.5
    c.ot.seed = 7
    assert c.voids.keep_divergence and c.voids.smoothing.sigma == 2.5 and c.ot.seed == 7
    assert c.rsd is None
    c.rsd = bitvf.BoxRsd(axis=0, redshift=0.5, bias=1.2, sigma=10.0)
    c.rsd.sigma = 3.0
    c.rsd.correction.weight_by_realizations = True
    assert c.rsd.sigma == 3.0 and c.rsd.correction.weight_by_realizations
    s = bitvf.Smoothing(sigma=4.0)
    c.voids.smoothing = s
    s.sigma = 5.0
    assert c.voids.smoothing.sigma == 4.0
    c.limits = ((0, 0, 0), (1, 2, 3))
    assert c.limits == otswap.BoxLimits((0, 0, 0), (1, 2, 3))
    c.rsd = None
    assert c.rsd is None


# --------------------------------------------------------------------------
# errors raised by the bindings

def test_error_class():
    assert issubclass(bitvf.Error, Exception) and not issubclass(bitvf.Error, otswap.Error)
    assert not issubclass(bitvf.Error, RuntimeError)


def error(call, message):
    with pytest.raises(bitvf.Error) as e:
        call()
    assert str(e.value) == message, str(e.value)


def test_argument_errors(tmp_path):
    tracers, displacement, valid, ids, box = translated()
    n = len(tracers)
    f = bitvf.find_voids_box
    with pytest.raises(bitvf.Error, match="^tracers cannot be converted to a float64 array: "):
        f("abc", displacement, mps=7.0)
    error(lambda: f(tracers[:, :2], displacement, mps=7.0),
          f"tracers must have shape (N, 3); it has shape ({n}, 2)")
    error(lambda: f(tracers, displacement[0], mps=7.0),
          "displacement must have shape (n_realizations, n_tracers, 3), neither 0; it has shape "
          f"({n}, 3)")
    error(lambda: f(tracers, displacement[:, :10], mps=7.0),
          f"displacement must have shape (n_realizations, {n}, 3), {n} the rows of tracers; it has shape "
          "(2, 10, 3)")
    error(lambda: f(tracers, displacement, valid[:1], mps=7.0),
          f"valid must have shape (2, {n}), as displacement; it has shape (1, {n})")
    error(lambda: bitvf.find_voids_box_from_field(np.zeros((4, 3)), np.zeros((2, 2))),
          "values must be one-dimensional; it has shape (2, 2)")
    error(lambda: bitvf.watershed_box(np.zeros(8)), "field must have shape (nx, ny, nz); it has shape (8,)")
    error(lambda: f(tracers, displacement, np.zeros((0,)), mps=7.0),
          "valid is empty; omit it to take every displacement as valid")
    error(lambda: bitvf.run_box(tracers, np.zeros((0, 3))),
          "randoms is empty; omit it to have bitvf draw the randoms in the box")
    error(lambda: bitvf.find_voids_box_from_field(np.zeros((4, 3)), np.zeros(4), weights=[]),
          "weights is empty; omit it for unit weights")
    error(lambda: f(tracers, displacement, valid * 2, mps=7.0),
          "the validity flags hold 2.000000 for object 1 (from 0); each must be 0 or 1")
    error(lambda: f(tracers, displacement, mps=7.0, ids=np.full(n, 0.5)),
          "ids holds 0.500000 for object 0 (from 0); every entry must be a 64-bit integer")
    error(lambda: f(tracers, displacement, mps=7.0, ids=["a"] * n),
          "ids must hold integers; its dtype is dtype('<U1')")
    error(lambda: bitvf.Config(keep_divergence=1), "keep_divergence must be True or False; got 1")
    error(lambda: bitvf.Config(n_realizations_mean=-1),
          "n_realizations_mean is -1; it must be a non-negative integer")
    error(lambda: bitvf.Smoothing(sigma="x"), "sigma must be a real number; got 'x'")
    error(lambda: bitvf.Config(verbosity="loud"),
          "verbosity must be \"silent\", \"normal\" or \"detailed\"; got 'loud'")
    error(lambda: bitvf.io.read_divergence(tmp_path / "x.fits", angle_unit="degrees"),
          "angle_unit must be \"deg\" or \"rad\"; got 'degrees'")
    divergence_file = tmp_path / "py14.fits"
    bitvf.io.write_divergence(divergence_file, bitvf.find_voids_box(tracers, displacement, mps=7.0, config=kept("silent")))
    error(lambda: bitvf.io.read_divergence(divergence_file, angle_unit="deg"),
          "angle_unit must not be given for a box file; a box file holds no angles")
    error(lambda: bitvf.io.read_field("x.dat", centres=("a", "b", "c"), value="v", delimiter=",,"),
          "delimiter must be a single character; got ',,'")
    error(lambda: bitvf.io.read_field("x.dat", centres="abc", value="v"),
          "centres must be a sequence of three column names (str); got 'abc'")
    error(lambda: bitvf.io.write_voids(3, None), "file must be a str or an os.PathLike; got 3")
    error(lambda: f(tracers, displacement, mps=7.0, config=otswap.Config()),
          "config must be a bitvf.Config; got " + repr(otswap.Config()))
    error(lambda: f(tracers, displacement, mps=7.0, grid=3), "grid must be a bitvf.DivergenceGrid; got 3")
    error(lambda: f(tracers, displacement, mps=7.0, box=(1, 2, 3)),
          "box must be an otswap.BoxLimits or a pair of corners ((x0, y0, z0), (x1, y1, z1)); it has shape (3,)")
    error(lambda: bitvf.Smoothing(unit="mps"), "unit must be a bitvf.LengthUnit; got 'mps'")
    error(lambda: bitvf.BoxConfig(rsd=1), "rsd must be a bitvf.BoxRsd or None; got 1")
    error(lambda: bitvf.run_box(tracers, distances=1), "distances must be an otswap.DistanceTable; got 1")
    error(lambda: bitvf.io.file_name("voids", "d", "n.fits", zspace=False),
          "product must be a bitvf.io.Product; got 'voids'")
    error(lambda: f(tracers, mps=7.0), "find_voids_box: displacement is required with arrays")
    error(lambda: f(tracers, displacement), "find_voids_box: mps is required with arrays")
    error(lambda: bitvf.box_grid(box=box, mps=1.0), "box_grid: give a reconstruction, or box, mps and cell_size")


def test_form_errors(reconstruction):
    error(lambda: bitvf.find_voids_box(reconstruction, mps=3.0),
          "find_voids_box: mps does not apply to a reconstruction (otswap.Result)")
    error(lambda: bitvf.find_voids_box(reconstruction, np.zeros((1, 1, 3)), mps=3.0, ids=[1]),
          "find_voids_box: displacement, mps and ids do not apply to a reconstruction (otswap.Result)")
    error(lambda: bitvf.box_grid(reconstruction, mps=1.0), "box_grid: box and mps do not apply to a reconstruction")
    error(lambda: bitvf.box_grid(3), "reconstruction must be an otswap.Result; got 3")


# The library's errors reach Python unchanged as bitvf.Error, otswap's
# prefixed.
def test_library_errors():
    tracers, displacement, _, _, _ = translated()
    with pytest.raises(bitvf.Error, match="mean particle separation"):
        bitvf.find_voids_box(tracers, displacement, mps=-1.0, config=quiet())
    error(lambda: bitvf.run_box(tracers, config=bitvf.BoxConfig(ot=otswap.Config(convergence=2.0), voids=quiet())),
          "otswap: the convergence threshold is 2.000000; it must lie strictly between 0 and 1")
    error(lambda: bitvf.io.file_name(bitvf.io.Product.VOIDS, "d", "", zspace=False), "the file name is empty")


def test_file_name():
    assert bitvf.io.file_name(bitvf.io.Product.VOIDS, "out", "box.fits", zspace=True) == "out/voids_zspace_box.fits"
    assert bitvf.io.file_name(bitvf.io.Product.GALAXIES, Path("out"), "b.dat", zspace=False) == "out/galaxies_b.dat"
