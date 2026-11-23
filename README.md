# bitvf

<!-- docs:about-start -->
BitVF, the Back-in-time Void Finder (Sartori et al. 2026), is a dynamical
void finder that identifies cosmic voids as non-overlapping regions of
coherent mass outflow. It reconstructs by optimal transport the displacement
field that moves each tracer of a catalogue back in time, from its observed
position to its initial position in a nearly uniform distribution. The
divergence of this field is negative where matter has flowed out: BitVF
computes it on a grid and divides the grid into basins around its local
minima, so that the basins whose minimum is negative are the voids and the
others are overdense regions. Each basin comes with a centre, radii and
densities. The method is described in
[arXiv:2601.15378](https://arxiv.org/abs/2601.15378).

BitVF works in two geometries. The **lightcone geometry** is that of a
galaxy survey, or of a mock catalogue built like one: tracers in right
ascension, declination and redshift, a random catalogue that samples the
survey volume and its systematics, and a HEALPix mask of the observed sky.
The **box geometry** is that of a simulation snapshot: tracers in Cartesian
coordinates filling a rectangular volume.

This package, bitvf, implements BitVF. The reconstruction is done by
[otswap](https://github.com/SimoSartori/otswap); bitvf runs it, can correct
the redshift-space distortions and reconstruct the corrected tracers again,
and searches the voids. Results are identical bit for bit on Linux (x86_64,
aarch64) and macOS (arm64, x86_64), whatever the number of threads. bitvf is a C++17 library, with a Python package that runs
the same code.
<!-- docs:about-end -->

It is the last of three packages: it depends on otswap and on
[meshsearch](https://github.com/SimoSartori/meshsearch). A version of BitVF is
included in CosmoBolognaLib.

Documentation, for both APIs: **<https://simosartori.github.io/bitvf/>**

## Install

```
pip install bitvf
```

installs the Python package, with otswap. From a clone, `pip install .`
builds the package from source, with the requirements below.

## A first run

The voids of a survey-like catalogue shipped with the repository, in
Python, from the root of a clone: otswap reconstructs the displacements of
the tracers, and bitvf searches the voids of the reconstruction and writes
them to a file.

```python
import bitvf
import otswap

# Right ascension and declination in degrees, and redshift: arrays of shape (N, 3)
tracers = otswap.io.read("examples/data/lightcone_tracers.dat", [0, 1, 2]).values
randoms = otswap.io.read("examples/data/lightcone_randoms.dat", [0, 1, 2]).values
mask = otswap.Mask("examples/data/lightcone_mask.fits.gz")
distances = otswap.DistanceTable.flat(0.286, 0.7)

reconstruction = otswap.reconstruct_lightcone(tracers, randoms, mask=mask, n_bins=30,
                                              distances=distances, angle_unit="deg",
                                              config=otswap.Config(n_realizations=8, seed=12345))
result = bitvf.find_voids_lightcone(reconstruction, randoms, mask=mask, angle_unit="deg")

print(result.voids.n_basins, "voids")
bitvf.io.write_voids("voids.dat", result)
```

The randoms must follow the footprint, the redshift distribution and the
systematics of the survey, and number at least the realizations times the
tracers; this catalogue has about eight per tracer, enough for 8
realizations. An analysis uses otswap's default of 30 realizations, and so
30 randoms per tracer as the starting point. The
[quickstart](https://simosartori.github.io/bitvf/quickstart.html) runs the
same search in C++ and Python, and the
[guides](https://simosartori.github.io/bitvf/lightcone/guide.html) explain
the parameters, the box geometry and the executable.

## Build from source

Needs CMake 3.16 or later, GCC or Clang (C, and C++17) and cfitsio; OpenMP
when found, from the same runtime family otswap was built with. otswap 0.3.x
is used when installed, and 0.3.0 is fetched from GitHub
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
