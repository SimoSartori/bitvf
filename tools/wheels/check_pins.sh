#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2026 by Simone Sartori, simone.sartori@inaf.it
#
# Fails unless the zlib and cfitsio versions and checksums of build_deps.sh
# are those of tools/wheels/build_deps.sh in the otswap tag
# that bitvf is built against (BITVF_OTSWAP_TAG in CMakeLists.txt), so that
# both packages' wheels link the same sources. Where the files are
# downloaded from is not compared.
#
# Usage: check_pins.sh PROJECT

set -euo pipefail

PROJECT=$1

tag=$(sed -n 's/^set(BITVF_OTSWAP_TAG "\(v[0-9.]*\)".*/\1/p' "$PROJECT/CMakeLists.txt")
repository=$(sed -n 's/^set(BITVF_OTSWAP_REPOSITORY "https:\/\/github.com\/\(.*\)\.git".*/\1/p' "$PROJECT/CMakeLists.txt")
if [ -z "$tag" ] || [ -z "$repository" ]; then
  echo "check_pins: no otswap tag or repository in $PROJECT/CMakeLists.txt" >&2
  exit 1
fi

theirs=$(mktemp "${TMPDIR:-/tmp}/otswap-build-deps.XXXXXX")
trap 'rm -f "$theirs"' EXIT
curl --fail --silent --show-error --location --retry 3 --output "$theirs" \
  "https://raw.githubusercontent.com/$repository/$tag/tools/wheels/build_deps.sh"

# the pins: the versions and checksums
pins () {
  grep -E '^(ZLIB_(VERSION|SHA256)|CFITSIO_(VERSION|SHA256|SHA512))=' "$1"
}

if [ -z "$(pins "$theirs")" ]; then
  echo "check_pins: no pins found in otswap $tag's build_deps.sh" >&2
  exit 1
fi
if ! diff <(pins "$theirs") <(pins "$PROJECT/tools/wheels/build_deps.sh"); then
  echo "check_pins: build_deps.sh differs from otswap $tag's (diff above: otswap's <, bitvf's >)" >&2
  exit 1
fi
echo "check_pins: zlib and cfitsio versions and checksums are those of otswap $tag"
