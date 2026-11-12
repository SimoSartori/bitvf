# Names side by side

The two APIs are one library. A C++ name in camelCase is the Python name in
snake_case; a C++ configuration struct is a Python class of the same name,
passed as `config=`; a flat `std::vector<double>` of `3 * N` entries becomes
a NumPy array of shape (N, 3). The types of otswap pass between the two
packages as they are: an `otswap.Result` from `otswap.reconstruct_box` goes
straight into `bitvf.find_voids_box`. The rows below differ from that rule,
or are listed for completeness.

## Errors

| C++ | Python | Notes |
|---|---|---|
| `bitvf::Error` | `bitvf.Error` | C++ records the throw location: `file()`, `function()`, `line()` |
| `otswap::Error` inside a bitvf call | `bitvf.Error` | the message starts `otswap: `; in C++ the original is nested |
| — | `KeyboardInterrupt` | Ctrl+C, at the next report line of a search |

## Configuration

| C++ | Python | Notes |
|---|---|---|
| `LengthUnit::Mps`, `LengthUnit::Coordinate` | `LengthUnit.MPS`, `LengthUnit.COORDINATE` | |
| `Smoothing`, `sigma`, `unit` | `Smoothing`, `sigma`, `unit` | |
| `Config`, its fields | `Config`, `smoothing`, `divergence_cell_size`, `keep_divergence`, `keep_assignment`, `n_realizations_mean`, `verbosity` | `divergence_cell_size` `None` for the default |
| `Config::verbosity` (`otswap::Verbosity`) | `"silent"`, `"normal"`, `"detailed"` | |
| `Config::log` | — | Python prints to `sys.stdout`, each line as its step finishes |
| `BoxConfig`, `ot`, `limits`, `mps`, `rsd`, `voids`, `keepReconstructions` | `BoxConfig`, `ot`, `limits`, `mps`, `rsd`, `voids`, `keep_reconstructions` | the nested configurations read as references into the `BoxConfig` |
| `BoxRsd`, `axis`, `redshift`, `bias`, `sigma`, `correction`, `zspaceVoids` | `BoxRsd`, `axis`, `redshift`, `bias`, `sigma`, `correction`, `zspace_voids` | `axis`, `redshift`, `bias`, `sigma` have no default |
| `otswap::BoxLimits` | `otswap.BoxLimits`; keyword `box`, field `limits` | also a pair `((x0, y0, z0), (x1, y1, z1))` |

## Functions

| C++ | Python | Notes |
|---|---|---|
| `findVoidsBox(tracers, displacement, valid, mps, config, limits, ids, grid)` | `find_voids_box(tracers, displacement, valid, *, mps, box, config, ids, grid)` | `displacement` (R, N, 3), `valid` (R, N) or `None` |
| `findVoidsBox(reconstruction, config, grid)` | `find_voids_box(reconstruction, *, config, grid)` | an `otswap.Result` as the first argument |
| `findVoidsBoxFromDivergence(divergence, tracers, config, ids)` | `find_voids_box_from_divergence(grid, tracers, *, config, ids)` | |
| `findVoidsBoxFromField(centres, values, config, weights, mps)` | `find_voids_box_from_field(centres, values, *, config, weights, mps)` | |
| `boxGrid(box, mps, divergenceCellSize)` | `box_grid(*, box, mps, cell_size)` | |
| `boxGrid(reconstruction, divergenceCellSize)` | `box_grid(reconstruction, *, cell_size)` | |
| `watershedBox(field, nCellsAxis, verbosity, log)` | `watershed_box(field, *, verbosity)` | `field` of shape (nx, ny, nz) carries `nCellsAxis` |
| `runBox(tracers, randoms, config, ids)` | `run_box(tracers, randoms, *, config, ids)` | `randoms` `None` to have them drawn |
| `runBox(tracers, randoms, distances, config, ids)` | `run_box(tracers, randoms, *, config, distances, ids)` | |

## Results

| C++ | Python | Notes |
|---|---|---|
| `Result`, `BasinCatalog`, `DivergenceGrid`, `TracerAssignment`, `Counts`, `Partition`, `BitVFResult` | the same names | read-only; no constructor |
| `std::optional` fields | `None` when empty | |
| per-cell arrays of a box grid, `[cell]` | shape (nx, ny, nz), C order | `a[i, j, k]` is cell (i, j, k); `centre` (nx, ny, nz, 3) |
| per-cell arrays of `Partition` | the shape of the field | |
| a cell row stored as a value (`cellIndex`, `voidSeed`, `overdenseSeed`) | a flat row | `numpy.unravel_index(row, n_cells_axis)` gives the cell |
| flags, `uint8` 0 or 1 | `bool` | inputs accept `bool` or 0/1 numbers |
| `BasinCatalog::edgeFlag` | `edge_flag`, `uint8` | 0, 1 or 2 |
| `Result::ids` | `ids` | a new copy on each access |
| `Result::randomSeed` | `random_seed` | `None` when no randoms were drawn |
| `Result::geometry` (`otswap::Geometry`) | `"box"`, `"lightcone"` | |
| `Result::elapsedSeconds` | `elapsed_seconds` | |

## Files

| C++ | Python | Notes |
|---|---|---|
| `io::Product` | `io.Product` | `VOIDS`, `OVERDENSE`, `DIVERGENCE`, `GALAXIES` |
| `io::fileName(product, directory, name, zspace)` | `io.file_name(product, directory, name, *, zspace)` | |
| `io::writeVoids`, `writeOverdense`, `writeDivergence`, `writeGalaxies` | `io.write_voids`, `write_overdense`, `write_divergence`, `write_galaxies` | the same files, byte for byte |
| `io::readDivergence(file)` | `io.read_divergence(file, *, angle_unit=None)` | `angle_unit` required for a lightcone file, an error for a box file |
| `io::readField(file, centreColumns, valueColumn, weightColumn, delimiter, comment)`, `io::Field` | `io.read_field(file, *, centres, value, weight, delimiter, comment)`, `io.Field` | `weights` `None` when no column is named |
| a path, `std::string` | `str` or `os.PathLike` | |
