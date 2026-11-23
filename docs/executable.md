---
orphan: true
---

# The executables

`bitvf` runs a void search from a parameter file, with no code to write;
`bitvf-field` searches a field given on a box lattice. Both are built and
installed with the C++ library. The keys depend on the geometry, set by the
key `geometry`; each page lists every key of its geometry, with the output
files and a commented parameter file:

- [Lightcone executable](lightcone/executable.md): `geometry = lightcone`.
- [Box executables](box/executable.md): `geometry = box`, and
  [`bitvf-field`](box/executable.md#bitvf-field).
