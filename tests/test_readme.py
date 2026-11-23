# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Tests of the code shown in README.md."""

import contextlib
import io
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


# The Python code of the README's "A first run" runs as shown from the root
# of a clone, prints the number of voids and writes the voids file.
def test_readme_first_run(tmp_path):
    readme = (ROOT / "README.md").read_text()
    section = readme.split("## A first run", 1)[1].split("\n## ", 1)[0]
    blocks = re.findall(r"```python\n(.*?)```", section, re.S)
    assert len(blocks) == 1

    (tmp_path / "examples").mkdir()
    (tmp_path / "examples" / "data").symlink_to(ROOT / "examples" / "data")
    printed = io.StringIO()
    here = os.getcwd()
    os.chdir(tmp_path)
    try:
        with contextlib.redirect_stdout(printed):
            exec(compile(blocks[0], "README.md", "exec"), {"__name__": "readme"})
    finally:
        os.chdir(here)

    assert re.search(r"^[1-9][0-9]* voids$", printed.getvalue(), re.M)
    assert (tmp_path / "voids.dat").stat().st_size > 0
