# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Bit-for-bit reproducibility across platforms and thread counts.

Every case below hashes the bytes of one output array, and the hashes must
equal the ones recorded in tests/determinism_hashes.json: the same inputs
give the same bits on every platform the wheels are built for, with one
thread and with four.

The inputs are identical in bytes everywhere: they come from a fixed 64-bit
linear congruential generator in pure Python, mapped to doubles by exact
division by 2^53 and then by single IEEE operations, with no libm call and
no NumPy random generator.

Each thread count runs in its own interpreter, since OpenMP reads
OMP_NUM_THREADS once. Run as a script, this file prints the hashes as JSON,
the content of the "hashes" entry of the reference file.
"""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "determinism_hashes.json"
THREADS = (1, 4)


# --------------------------------------------------------------------------
# inputs

class Lcg:
    """Knuth's MMIX generator, x <- 6364136223846793005 x + 1442695040888963407
    mod 2^64; a double from the top 53 bits."""

    def __init__(self, seed):
        self.state = seed & 0xFFFFFFFFFFFFFFFF

    def next(self):
        self.state = (6364136223846793005 * self.state + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        return self.state

    def unit(self):
        return (self.next() >> 11) / 9007199254740992.0

    def uniform(self, lo, hi):
        return lo + self.unit() * (hi - lo)


def box_points(n, size, seed):
    g = Lcg(seed)
    return np.array([[g.uniform(0.0, size) for _ in range(3)] for _ in range(n)])


def translated_displacement(tracers, n_realizations, seed):
    """Tracers with x below half the box move by about 30 along x, the
    others by at most 2."""
    g = Lcg(seed)
    out = np.empty((n_realizations, len(tracers), 3))
    for r in range(n_realizations):
        for i, t in enumerate(tracers):
            moved = t[0] < 50.0
            out[r, i] = [g.uniform(27.0, 33.0) if moved else g.uniform(-2.0, 2.0),
                         g.uniform(-3.0, 3.0) if moved else g.uniform(-2.0, 2.0),
                         g.uniform(-3.0, 3.0) if moved else g.uniform(-2.0, 2.0)]
    return out


# --------------------------------------------------------------------------
# cases

def digest(a):
    a = np.asarray(a)
    if a.dtype == np.bool_:
        a = a.astype("|u1")
    elif a.dtype.kind == "f":
        a = a.astype("<f8")
    elif a.dtype.kind in "iu":
        a = a.astype("<i8")
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def search_digests(out, key, r):
    for kind in ("voids", "overdense"):
        c = getattr(r, kind)
        out[f"{key}.{kind}.centre"] = digest(c.centre)
        out[f"{key}.{kind}.radius"] = digest(c.radius)
        out[f"{key}.{kind}.edge_flag"] = digest(c.edge_flag)
        out[f"{key}.{kind}.contamination"] = digest(c.contamination)
        if c.delta_av.size:
            out[f"{key}.{kind}.delta_av"] = digest(c.delta_av)
            out[f"{key}.{kind}.delta_2mps"] = digest(c.delta_2mps)
    if r.divergence is not None:
        out[f"{key}.divergence"] = digest(r.divergence.divergence)
        if r.divergence.divergence_smoothed.size:
            out[f"{key}.divergence_smoothed"] = digest(r.divergence.divergence_smoothed)
            out[f"{key}.weight_smoothed"] = digest(r.divergence.weight_smoothed)
        out[f"{key}.void_index"] = digest(r.divergence.void_index)
    if r.assignment is not None:
        out[f"{key}.assignment.distance"] = digest(r.assignment.distance)


def compute():
    """The hash of every case, by name."""
    import otswap
    import bitvf

    out = {}
    keep = dict(keep_divergence=True, keep_assignment=True, verbosity="silent")

    tracers = box_points(2000, 100.0, 31)
    displacement = translated_displacement(tracers, 2, 32)
    r = bitvf.find_voids_box(tracers, displacement, mps=7.0, box=((-5.0, -5.0, -5.0), (135.0, 105.0, 105.0)),
                             config=bitvf.Config(smoothing=bitvf.Smoothing(sigma=0.5), **keep))
    search_digests(out, "find_voids_box", r)
    again = bitvf.find_voids_box_from_divergence(
        bitvf.find_voids_box(tracers, displacement, mps=7.0,
                             config=bitvf.Config(**keep)).divergence,
        tracers, config=bitvf.Config(smoothing=bitvf.Smoothing(sigma=1.5), **keep))
    search_digests(out, "find_voids_box_from_divergence", again)

    g = Lcg(41)
    centres = np.array([[1.5 * i, 1.5 * j, 1.5 * k] for i in range(12) for j in range(15) for k in range(9)])
    values = np.array([g.uniform(-1.0, 1.0) + 0.1 * c[0] - 1.0 for c in centres])
    weights = np.array([g.uniform(0.5, 2.0) for _ in centres])
    f = bitvf.find_voids_box_from_field(centres, values, weights=weights, mps=1.5,
                                        config=bitvf.Config(smoothing=bitvf.Smoothing(sigma=1.2),
                                                            keep_divergence=True, verbosity="silent"))
    search_digests(out, "find_voids_box_from_field", f)
    p = bitvf.watershed_box(values.reshape(12, 15, 9), verbosity="silent")
    out["watershed_box.void_index"] = digest(p.void_index)
    out["watershed_box.overdense_index"] = digest(p.overdense_index)

    table = otswap.DistanceTable.flat(0.31, 0.68, z_max=2.0, n_samples=20000)
    config = bitvf.BoxConfig(ot=otswap.Config(n_realizations=2, seed=777, verbosity="silent"),
                             rsd=bitvf.BoxRsd(axis=2, redshift=0.5, bias=1.5, sigma=8.0,
                                              correction=otswap.CorrectionConfig(verbosity="silent"),
                                              zspace_voids=True),
                             voids=bitvf.Config(**keep))
    run = bitvf.run_box(box_points(800, 100.0, 51), config=config, distances=table)
    search_digests(out, "run_box", run.voids)
    search_digests(out, "run_box.zspace", run.zspace_voids)
    out["run_box.uncorrected"] = digest(run.voids.uncorrected)
    return out


def compute_with_threads(threads):
    """compute() in a fresh interpreter with OMP_NUM_THREADS = threads."""
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    run = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve())], env=env,
                         capture_output=True, text=True, check=False)
    if run.returncode != 0:
        raise RuntimeError(f"the cases failed with {threads} threads:\n{run.stderr}")
    return json.loads(run.stdout)


# --------------------------------------------------------------------------
# tests

def test_hashes_do_not_depend_on_the_thread_count():
    one, four = (compute_with_threads(t) for t in THREADS)
    assert one == four


def test_hashes_match_the_reference():
    expected = json.loads(REFERENCE.read_text())["hashes"]
    for threads in THREADS:
        got = compute_with_threads(threads)
        differ = sorted(k for k in expected if k in got and got[k] != expected[k])
        assert not differ, f"with {threads} threads these outputs changed: {differ}"
        assert set(got) == set(expected), (
            f"the cases differ from the recorded ones: new {sorted(set(got) - set(expected))}, "
            f"gone {sorted(set(expected) - set(got))}")


if __name__ == "__main__":
    print(json.dumps(compute(), indent=1, sort_keys=True))
