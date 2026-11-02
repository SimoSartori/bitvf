# bitvf

<!-- docs:about-start -->
bitvf is the Back-in-time Void Finder (BitVF): it finds cosmic voids in the
displacement field that optimal transport reconstructs from a catalogue of
large-scale structure tracers. The divergence of the displacement field is
computed on a grid, smoothed, and split by a watershed into basins around its
local minima: voids, where the reconstructed flow diverges, and overdense
basins. Each basin comes with a centre, radii and densities. The method is
described in Sartori et al. 2026,
[arXiv:2601.15378](https://arxiv.org/abs/2601.15378).

The reconstruction is done by [otswap](https://github.com/SimoSartori/otswap);
bitvf runs it, optionally corrects the redshift-space distortions and
reconstructs again in real space, and searches the voids. Results are
identical bit for bit on every supported platform and compiler, whatever the
number of threads. C++17, with a Python package on the same code.
<!-- docs:about-end -->

It is the last of three packages: it depends on otswap and on
[meshsearch](https://github.com/SimoSartori/meshsearch). A version of BitVF is
included in CosmoBolognaLib.

Documentation, for both APIs: **<https://simosartori.github.io/bitvf/>**

## Install

```
pip install bitvf
```

From a clone, `pip install .` builds the package from source, with the
requirements below.

## Build from source

Needs CMake 3.16 or later, GCC or Clang (C, and C++17) and cfitsio; OpenMP
when found, from the same runtime family otswap was built with. otswap 0.2.1
or a later 0.2.x is used when installed, and 0.2.1 is fetched from GitHub
otherwise, with meshsearch.

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /path/to/prefix
```

This installs the library, its headers and CMake package, and the `bitvf` and
`bitvf-field` executables.

```cmake
find_package(bitvf 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE bitvf::bitvf)
```

The options: `BITVF_BUILD_TESTS`, `BITVF_BUILD_PYTHON`, `BITVF_INSTALL`,
`BITVF_BUILD_EXECUTABLE`, and `BITVF_OTSWAP_REPOSITORY`, `BITVF_OTSWAP_TAG`
for an otswap fetched from elsewhere.

## Tests

```
ctest --test-dir build
BITVF_PYTHON_REFERENCE=build/python_reference pytest tests/
```

## Citing

Please cite Sartori et al. 2026,
[arXiv:2601.15378](https://arxiv.org/abs/2601.15378); `CITATION.cff` has the
details.

## Licence

GPL-2.0-or-later; see `LICENSE`. The example data in `examples/data` are not
covered by this licence: they may be redistributed unchanged with the examples
and used only to run them (see `examples/data/README.md`).
