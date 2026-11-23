# Shared parts of the example pages

Not a page of its own: the box and lightcone example pages include the parts
between the markers.

<!-- intro-start -->
Each example is a C++ program and a Python script that do the same thing:
the same inputs, parameters and seed, and the same output files, byte for
byte. Each folder in {repo-tree}`examples` is a CMake project of its own,
built against an installed bitvf.
<!-- intro-end -->

<!-- data-start -->
The examples read the files of {repo-tree}`examples/data`, which come with a
clone of the repository; they are not in the Python package or the source
distribution. They are provided for running the examples, may be
redistributed unchanged with them, and may not be used for any other purpose:
{repo-file}`examples/data/README.md` gives the terms and the content of each
file.
<!-- data-end -->

<!-- building-start -->
For an example `<name>`, with bitvf installed under `<prefix>`:

```
cmake -S examples/<name> -B build/examples/<name> -DCMAKE_PREFIX_PATH=<prefix>
cmake --build build/examples/<name>
build/examples/<name>/<name> [data_dir]
python examples/<name>/<name>.py [data_dir]
```

Both read the examples' data folder unless `data_dir` is given, and write to
`examples/<name>/output/`, created if missing. Where cfitsio is not in a
system path, give its prefix with `-DCFITSIO_ROOT=...` (and `-DZLIB_ROOT=...`
for a static cfitsio).
<!-- building-end -->
