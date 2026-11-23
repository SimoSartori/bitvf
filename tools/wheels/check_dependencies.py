# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
"""Fails unless every repaired wheel in a directory depends on nothing
outside itself but the allowed libraries, and prints what each binary in it
links.

macOS: the wheel ships no library of its own: no .dylib, and no .dylibs/
directory, where delocate would put a copy of libomp or cfitsio. Every load
command of the module names a library under /usr/lib or /System/Library, or
@rpath/libomp.dylib, which otswap's wheel ships; the module's only rpath is
@loader_path/../otswap/openmp/lib, where that libomp is.
Linux: auditwheel show must find the wheel consistent with a manylinux tag
among its own, which admits no external library beyond that policy's;
libgomp is the copy auditwheel put in the wheel.

Usage: python check_dependencies.py WHEEL_DIRECTORY
"""

import json
import os
import subprocess
import sys
import tempfile
import zipfile

SYSTEM_PREFIXES = ("/usr/lib/", "/System/Library/")
OTSWAP_OMP_RPATH = "@loader_path/../otswap/openmp/lib"


def check_macos(wheel):
    problems = []
    with tempfile.TemporaryDirectory() as root:
        with zipfile.ZipFile(wheel) as z:
            z.extractall(root)
        modules = []
        for folder, _, files in os.walk(root):
            for name in files:
                path = os.path.join(folder, name)
                rel = os.path.relpath(path, root)
                if name.endswith(".dylib") or ".dylibs" in rel.split(os.sep):
                    problems.append(f"{rel}: a library bundled in the wheel")
                elif name.endswith(".so"):
                    modules.append(path)
        if not modules:
            problems.append("no extension module in the wheel")
        for module in modules:
            rel = os.path.relpath(module, root)
            print(rel)
            out = subprocess.run(["otool", "-L", module], check=True,
                                 capture_output=True, text=True).stdout.splitlines()
            for dep in (line.split(" (compatibility")[0].strip() for line in out[1:]):
                print("    " + dep)
                if not dep.startswith(SYSTEM_PREFIXES) and dep != "@rpath/libomp.dylib":
                    problems.append(f"{rel} links {dep}")
            load = subprocess.run(["otool", "-l", module], check=True,
                                  capture_output=True, text=True).stdout.splitlines()
            rpaths = [load[i + 2].split()[1] for i, line in enumerate(load)
                      if line.strip() == "cmd LC_RPATH"]
            print("    rpath " + ", ".join(rpaths))
            if rpaths != [OTSWAP_OMP_RPATH]:
                problems.append(f"{rel} has rpaths {rpaths}; the only one allowed is {OTSWAP_OMP_RPATH}")
    return problems


def check_linux(wheel):
    shown = subprocess.run(["auditwheel", "show", wheel], check=True,
                           capture_output=True, text=True).stdout
    print(shown)
    report = json.loads(subprocess.run(["auditwheel", "show", "--json", wheel], check=True,
                                       capture_output=True, text=True).stdout)
    tags = os.path.basename(wheel)[:-len(".whl")].split("-")[-1].split(".")
    overall = report["overall_tag"]
    if not overall.startswith("manylinux_") or overall not in tags:
        return [f"auditwheel finds {os.path.basename(wheel)} consistent with {overall}, "
                "which is not a manylinux tag of the wheel"]
    return []


def main():
    directory = sys.argv[1]
    wheels = [os.path.join(directory, f) for f in sorted(os.listdir(directory)) if f.endswith(".whl")]
    if not wheels:
        sys.exit(f"check_dependencies: no wheel in {directory}")
    problems = []
    for wheel in wheels:
        print(f"== {os.path.basename(wheel)}")
        problems += check_macos(wheel) if sys.platform == "darwin" else check_linux(wheel)
    if problems:
        sys.exit("check_dependencies: " + "; ".join(problems))


if __name__ == "__main__":
    main()
