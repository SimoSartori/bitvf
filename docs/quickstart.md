# Quickstart

## Install

```
pip install bitvf
```

installs the Python package, with otswap, which bitvf uses for the
reconstruction. bitvf and otswap must be installed in the same environment,
as pip does.

The C++ library and the `bitvf` and `bitvf-field` executables are built and
installed with CMake. The build needs GCC or Clang (C, and C++17) and
cfitsio, uses OpenMP when found, and fetches otswap and meshsearch when they
are not installed:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /path/to/prefix
```

A CMake project of your own then uses the library with

```cmake
find_package(bitvf 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE bitvf::bitvf)
```

## A first lightcone run

The lightcone geometry is that of a galaxy survey: tracers observed on a part
of the sky, each with a right ascension, a declination and a redshift. This
first run finds the voids of a small survey-like catalogue shipped with the
repository, in {repo-tree}`examples/data`: about 30,000 tracers between
redshifts 0.885 and 1.10, a *random catalogue*, a set of points with the
same footprint and redshift distribution as the tracers but no clustering,
and the footprint itself as a HEALPix *mask*, a map of the pixels of the sky
that the survey observed. The data come with a clone of the repository, not
with the Python package:

```
git clone https://github.com/SimoSartori/bitvf.git
cd bitvf
```

The run has two steps. otswap first reconstructs the *displacements*, the
vectors that move each tracer back in time to a nearly uniform
distribution; bitvf then computes their divergence on a grid and finds the
voids. The code below is the core of the [lightcone
example](lightcone/examples.md), where `data` is the data folder and
`output` the folder the files are written to:

::::{tab-set}
:sync-group: lang

:::{tab-item} Python
:sync: python

```{literalinclude} ../examples/lightcone/lightcone.py
:language: python
:start-at: "# Right ascension and declination in degrees"
:end-at: bitvf.io.write_galaxies(
:dedent: 4
```
:::

:::{tab-item} C++
:sync: cpp

```{literalinclude} ../examples/lightcone/lightcone.cpp
:language: cpp
:start-at: "// Right ascension, declination and redshift for each object"
:end-at: bitvf::io::writeGalaxies(galaxiesFile, result);
:dedent: 4
```
:::

::::

A few points of the code deserve a word.

- **Angles.** The catalogues give the angles in degrees. Python takes them
  as they are, with `angle_unit="deg"`, and returns angles in the same unit;
  the C++ library works in radians, so the program converts them with
  `otswap::skyToRadians`.
- **Randoms.** They describe the survey, so they must follow its footprint,
  its redshift distribution and its systematics. Each realization of the
  reconstruction uses as many randoms as there are tracers, so at least the
  realizations times the tracers are needed; this catalogue has about eight
  per tracer, enough for the 8 realizations of the run. For an analysis, 30
  randoms per tracer, the minimum with otswap's default of 30 realizations,
  is a reasonable start.
- **Cosmology.** The distance table converts redshifts into comoving
  distances, here in Mpc/h, the unit of every length in the result.
- **Mean particle separation.** The *mps*, the cube root of the volume per
  tracer, changes with redshift in a survey; the reconstruction measures it
  in 30 redshift bins, and it sets the size of the grid cells and the width
  of the smoothing.
- **Search.** The search takes everything else from the reconstruction; the
  randoms and the mask are given again, because the grid is built on them.

To run the whole example in Python:

```
python examples/lightcone/lightcone.py
```

or in C++, with bitvf installed under `<prefix>`:

```
cmake -S examples/lightcone -B build/examples/lightcone -DCMAKE_PREFIX_PATH=<prefix>
cmake --build build/examples/lightcone
build/examples/lightcone/lightcone
```

Both print the steps as they run, then the number of voids and of
overdense basins, and write three files to `examples/lightcone/output/`:
`voids_lightcone.dat`, `overdense_basins_lightcone.dat` and
`galaxies_lightcone.dat`, the last with the basin of every tracer. The two
programs write the same files, byte for byte, with any number of threads.

## Where to go next

- [The method](method.md) explains what each step computes.
- The [lightcone guide](lightcone/guide.md) explains the inputs and the
  parameters, the correction of redshift-space distortions, and how to read
  and check the results; the [full reconstruction
  example](lightcone/examples.md#full-reconstruction) runs the whole
  pipeline, correction included, in one call.
- The [lightcone executable](lightcone/executable.md) runs the same search
  from a parameter file, with no code to write.

## The box geometry

For a simulation snapshot, with tracers in Cartesian coordinates filling a
rectangular volume, bitvf searches in the box geometry. There the randoms
are optional, since otswap can draw them uniformly in the box, and no mask
is needed: {cpp:func}`bitvf::runBox` ({py:func}`bitvf.run_box`) goes from
the tracers to the voids in one call. The [box guide](box/guide.md) and the
[box examples](box/examples.md) show how.
