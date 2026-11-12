# Quickstart

## Install

```
pip install bitvf
```

installs the Python package, with otswap. On macOS the package uses the OpenMP
runtime shipped inside otswap's wheel, found at `../otswap/openmp/lib` from
the bitvf package: install both in the same `site-packages`. For the C++ library and the
`bitvf` and `bitvf-field` executables, build and install with CMake; it needs
GCC or Clang (C, and C++17) and cfitsio, uses OpenMP when found, and fetches
otswap and meshsearch when they are not installed:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /path/to/prefix
```

and use it from your own CMake project:

```cmake
find_package(bitvf 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE bitvf::bitvf)
```

## Voids in a few lines

The voids of tracers in a box, given as an array: otswap reconstructs their
displacements, with randoms drawn in the box, and bitvf searches the voids
and writes them to a file.

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```python
import bitvf
import otswap

config = bitvf.BoxConfig(ot=otswap.Config(n_realizations=8, seed=12345))

run = bitvf.run_box(tracers, config=config)
print(run.voids.voids.n_basins, "voids")
bitvf.io.write_voids("voids.dat", run.voids)
```
:::

:::{tab-item} C++
:sync: cpp

```cpp
#include "bitvf/BitVF.h"
#include "bitvf/io.h"

bitvf::BoxConfig config;
config.ot.nRealizations = 8;
config.ot.seed = 12345;

const bitvf::BitVFResult run = bitvf::runBox(tracers, {}, config);
std::cout << run.voids.voids.nBasins << " voids\n";
bitvf::io::writeVoids("voids.dat", run.voids);
```
:::

::::

`tracers` is an array of shape (N, 3) in Python, a `std::vector<double>` of
`3 * N` entries in C++, in coordinate units. The [box guide](box/guide.md)
explains the parameters, and the [examples](box/examples.md) read real
catalogues and show every step.
