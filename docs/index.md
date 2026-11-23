# bitvf

```{include} ../README.md
:start-after: <!-- docs:about-start -->
:end-before: <!-- docs:about-end -->
```

The method is described in the paper, Sartori et al. 2026,
[arXiv:2601.15378](https://arxiv.org/abs/2601.15378). The C++ library and the
Python package run the same code and give the same results, bit for bit;
[names side by side](naming.md) maps one API to the other.

## Install

```
pip install bitvf
```

The C++ library, with its executables, is built and installed with CMake, as
the [quickstart](quickstart.md) shows.

## Where to go next

- [Quickstart](quickstart.md): install, and voids in a few lines.
- [The method](method.md), for both geometries.
- Lightcone: the [guide](lightcone/guide.md), the [file
  formats](lightcone/formats.md), the [examples](lightcone/examples.md) and
  the [executable](lightcone/executable.md).
- Box: the [guide](box/guide.md), the [file formats](box/formats.md), the
  [examples](box/examples.md) and the [executables](box/executable.md).
- Reference: [C++](cpp.rst), [Python](python.rst).

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
