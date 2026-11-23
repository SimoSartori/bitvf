#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
#
# macOS, before each wheel build: unpacks the libomp that otswap's wheel of
# the architecture named in ARCHFLAGS ships (otswap/openmp/), and points
# PREFIX/omp, the OpenMP_ROOT of the build, at it. bitvf then links the
# runtime that otswap loads, whichever architecture the build host has.
#
# The wheel is the newest otswap 0.3.x for that architecture, for the
# CPython of the build.
#
# Usage: select_otswap_omp.sh PREFIX

set -euo pipefail

PREFIX=$1

case "${ARCHFLAGS:-}" in
  "-arch arm64")  arch=arm64;  platforms=(macosx_11_0_arm64) ;;
  "-arch x86_64") arch=x86_64; platforms=(macosx_10_13_x86_64 macosx_11_0_x86_64) ;;
  *) echo "select_otswap_omp: ARCHFLAGS is '${ARCHFLAGS:-}'; one of arm64 or x86_64 is required" >&2
     exit 1 ;;
esac

abi=$(python -c 'import sys; print(f"cp{sys.version_info[0]}{sys.version_info[1]}")')
dir="$PREFIX/otswap-omp-$arch"
work=$(mktemp -d "${TMPDIR:-/tmp}/otswap-wheel.XXXXXX")
trap 'rm -rf "$work"' EXIT

args=()
for p in "${platforms[@]}"; do args+=(--platform "$p"); done
python -m pip download --quiet --no-deps --only-binary :all: "${args[@]}" \
  --implementation cp --python-version "${abi#cp}" --abi "$abi" \
  --dest "$work" "otswap>=0.3,<0.4"

wheel=$(ls "$work"/otswap-*.whl)
rm -rf "$dir"
mkdir -p "$dir"
unzip -q "$wheel" 'otswap/openmp/*' -d "$work/unpacked"
cp -R "$work/unpacked/otswap/openmp/." "$dir/"

if [ ! -f "$dir/lib/libomp.dylib" ] || [ ! -f "$dir/include/omp.h" ]; then
  echo "select_otswap_omp: $(basename "$wheel") holds no otswap/openmp/{lib/libomp.dylib,include/omp.h}" >&2
  exit 1
fi

ln -sfn "otswap-omp-$arch" "$PREFIX/omp"
echo "select_otswap_omp: $PREFIX/omp -> otswap-omp-$arch, from $(basename "$wheel")"
