# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Reports the OpenMP runtimes loaded by one process that imports otswap and
bitvf and runs a parallel call of each, for tools/wheels/test_wheel.sh.

The process must exit 0. On macOS, dyld must also list a single libomp.dylib
image, the one shipped in the otswap package; when dyld prints nothing, the
check was not made and fails. On Linux each wheel carries its
own libgomp, so the number of libgomp images mapped is reported, without a
verdict.

Usage: python check_openmp.py
"""

import os
import subprocess
import sys

import otswap

CHILD = """
import sys
import numpy as np
import otswap
import bitvf

tracers = np.random.default_rng(1).uniform(0.0, 60.0, (3000, 3))
ot = otswap.Config(n_realizations=2, seed=3, verbosity="silent")
rec = otswap.reconstruct_box(tracers, config=ot)
voids = bitvf.find_voids_box(rec, config=bitvf.Config(verbosity="silent"))
again = otswap.reconstruct_box(tracers, config=ot)
assert np.array_equal(again.displacement, rec.displacement)
assert voids.voids.n_basins > 0
print("otswap and bitvf parallel calls done;", voids.voids.n_basins, "basins")
if sys.platform.startswith("linux"):
    maps = {line.split()[-1] for line in open("/proc/self/maps") if "libgomp" in line}
    print("LIBGOMP", *sorted(maps))
"""


def main():
    env = dict(os.environ, OMP_NUM_THREADS="4")
    if sys.platform == "darwin":
        env["DYLD_PRINT_LIBRARIES"] = "1"
    proc = subprocess.run([sys.executable, "-c", CHILD], env=env, capture_output=True, text=True)
    errors = [line for line in proc.stderr.splitlines() if not line.startswith("dyld[")]
    sys.stdout.write("".join(line + "\n" for line in proc.stdout.splitlines() if not line.startswith("LIBGOMP")))
    if errors:
        print("stderr, without dyld's lines:\n    " + "\n    ".join(errors))
    if proc.returncode != 0 or "OMP: Error #15" in proc.stderr:
        sys.exit(f"check_openmp: failed: exit {proc.returncode}")

    if sys.platform == "darwin":
        images = [line.split(": ", 1)[1].split("> ", 1)[-1].strip() for line in proc.stderr.splitlines()
                  if line.startswith("dyld[") and ": " in line]
        if not images:
            sys.exit("check_openmp: FAILED: NOT MADE: DYLD_PRINT_LIBRARIES printed nothing, so the runtimes loaded are unknown")
        loaded = [p for p in images if os.path.basename(p) == "libomp.dylib"]
        shipped = os.path.join(otswap.openmp_root() or "", "lib", "libomp.dylib")
        ok = len(loaded) == 1 and os.path.realpath(loaded[0]) == os.path.realpath(shipped)
        print(f"check_openmp: {'ok' if ok else 'FAILED'}: {len(images)} images loaded; libomp.dylib: {loaded}")
        if not ok:
            sys.exit(1)
    else:
        gomp = [w for w in proc.stdout.splitlines() if w.startswith("LIBGOMP")]
        count = len(gomp[0].split()) - 1 if gomp else 0
        print(f"check_openmp: {count} libgomp image(s) mapped in the process; exit 0")


if __name__ == "__main__":
    main()
