# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""otswap and bitvf in one process.

bitvf takes and returns otswap's Python types, and both packages run
OpenMP: imported together, the two must share one OpenMP runtime, or the
first parallel call of the second one aborts the process. Each test runs
in its own interpreter, so that an abort fails the test instead of the
session.
"""

import os
import subprocess
import sys
import textwrap

SCRIPT = textwrap.dedent("""
    import numpy as np
    import otswap
    import bitvf

    rng = np.random.default_rng(8)
    tracers = rng.uniform(0.0, 60.0, (3000, 3))
    rec = otswap.reconstruct_box(tracers, config=otswap.Config(n_realizations=2, seed=3, verbosity="silent"))
    voids = bitvf.find_voids_box(rec, config=bitvf.Config(verbosity="silent"))
    again = otswap.reconstruct_box(tracers, config=otswap.Config(n_realizations=2, seed=3, verbosity="silent"))
    assert np.array_equal(again.displacement, rec.displacement)
    assert voids.voids.n_basins > 0 and voids.n_tracers == 3000
    run = bitvf.run_box(tracers, config=bitvf.BoxConfig(
        ot=otswap.Config(n_realizations=2, seed=3, verbosity="silent"),
        voids=bitvf.Config(verbosity="silent"), keep_reconstructions=True))
    assert isinstance(run.reconstruction, otswap.Result)
    print(voids.voids.n_basins, run.voids.voids.n_basins)
""")


def run(threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    return subprocess.run([sys.executable, "-B", "-c", SCRIPT], env=env, capture_output=True, text=True)


# A parallel call of otswap, then of bitvf, then of otswap again, on several
# threads, in one process; the counts do not depend on the thread count.
def test_parallel_calls_of_both_packages_in_one_process():
    four = run(4)
    assert four.returncode == 0, four.stderr
    one = run(1)
    assert one.returncode == 0, one.stderr
    assert four.stdout == one.stdout
