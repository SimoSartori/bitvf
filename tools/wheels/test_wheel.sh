#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
#
# Runs the Python tests against an installed wheel: the comparison tests
# against python_reference, which is built for the wheel's architecture from
# the same sources and against the same dependencies, and whose absence is an
# error rather than a skip; the reproducibility test, which compares the
# wheel's outputs with the hashes recorded in tests/determinism_hashes.json;
# the README's code; and the tests that run otswap and bitvf in one process.
# Last, tools/wheels/check_openmp.py reports which OpenMP runtimes that process
# loads.
#
# Usage: test_wheel.sh PROJECT WHEEL
#
# python_reference is built once per architecture, in
# PROJECT/build/wheel-deps/python-reference-ARCH.

set -euo pipefail

PROJECT=$1
WHEEL=$2
DEPS="$PROJECT/build/wheel-deps"

case "$WHEEL" in
  *_x86_64.whl)  arch=x86_64 ;;
  *_arm64.whl)   arch=arm64 ;;
  *_aarch64.whl) arch=aarch64 ;;
  *) echo "test_wheel: no architecture in $WHEEL" >&2; exit 1 ;;
esac

build="$DEPS/python-reference-$arch"
reference="$build/python_reference"

if [ ! -x "$reference" ]; then
  args=(-DCMAKE_BUILD_TYPE=Release -DBITVF_BUILD_TESTS=ON -DBITVF_BUILD_EXECUTABLE=OFF
        -DCFITSIO_ROOT="$DEPS")
  if [ "$(uname)" = Darwin ]; then
    # the wheel's own deployment target, from its tag: macosx_10_13_x86_64; and
    # the libomp of otswap's wheel of this architecture, which the wheel loads
    target=$(basename "$WHEEL" | sed -E 's/.*-macosx_([0-9]+)_([0-9]+)_.*/\1.\2/')
    args+=(-DCMAKE_OSX_ARCHITECTURES="$arch" -DCMAKE_OSX_DEPLOYMENT_TARGET="$target"
           -DOpenMP_ROOT="$DEPS/otswap-omp-$arch")
  fi
  ZLIB_ROOT="$DEPS" cmake -S "$PROJECT" -B "$build" "${args[@]}"
  cmake --build "$build" --target python_reference --parallel
fi

tests="$PROJECT/tests"
BITVF_PYTHON_REFERENCE="$reference" python -m pytest \
  "$tests/test_bitvf.py" "$tests/test_lightcone_python.py" "$tests/test_crossmodule.py" \
  "$tests/test_import.py" "$tests/test_determinism.py" "$tests/test_readme.py" \
  -p no:cacheprovider -rs

python "$PROJECT/tools/wheels/check_openmp.py"
