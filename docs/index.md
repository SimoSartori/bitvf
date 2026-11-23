# bitvf

```{include} ../README.md
:start-after: <!-- docs:about-start -->
:end-before: <!-- docs:about-end -->
```

The two give the same results, bit for bit, and [names side by
side](naming.md) maps the C++ API to the Python one.

## Install

```
pip install bitvf
```

installs the Python package, with otswap. The C++ library, with its
executables, is built and installed with CMake, as the
[quickstart](quickstart.md) shows.

## Where to go next

- The [quickstart](quickstart.md): install, and a first search of the voids
  of a survey-like catalogue shipped with the repository.
- [The method](method.md): what each step computes, in both geometries.
- **Lightcone geometry**, for surveys and their mocks: the
  [guide](lightcone/guide.md) explains the inputs, the parameters and the
  results; the [examples](lightcone/examples.md) run a search and the whole
  pipeline; the [executable](lightcone/executable.md) runs them from a
  parameter file; the [file formats](lightcone/formats.md) describe what is
  written.
- **Box geometry**, for simulation snapshots: the [guide](box/guide.md), the
  [examples](box/examples.md), the [executables](box/executable.md) and the
  [file formats](box/formats.md).
- Reference: every function and structure of the [C++](cpp.rst) and
  [Python](python.rst) libraries.
- [Citing](citing.md): the paper to cite when using bitvf.

```{toctree}
:maxdepth: 2
:hidden:

quickstart
method
```

```{toctree}
:caption: Lightcone
:maxdepth: 2
:hidden:

lightcone/guide
lightcone/formats
lightcone/examples
lightcone/executable
```

```{toctree}
:caption: Box
:maxdepth: 2
:hidden:

box/guide
box/formats
box/examples
box/executable
```

```{toctree}
:caption: Reference
:maxdepth: 2
:hidden:

cpp
python
naming
```

```{toctree}
:caption: About
:maxdepth: 1
:hidden:

changelog
citing
```
