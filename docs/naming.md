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
| `LightconeConfig`, `ot`, `nBins`, `cut`, `rsd`, `voids`, `keepReconstructions` | `LightconeConfig`, `ot`, `n_bins`, `cut`, `rsd`, `voids`, `keep_reconstructions` | `n_bins` has no default; the nested configurations read as references |
| `LightconeRsd`, `bias`, `sigma`, `correction`, `secondPassRandoms`, `zspaceVoids` | `LightconeRsd`, `bias`, `sigma`, `correction`, `second_pass_randoms`, `zspace_voids` | `bias`, `sigma` have no default |
| `SecondPassRandoms::None`, `Resample`, `External` | `SecondPassRandoms.NONE`, `RESAMPLE`, `EXTERNAL` | |
| `otswap::RedshiftCut` | `otswap.RedshiftCut`; keyword `redshift_cut`, field `cut` | also a pair `(min, max)` |

## Lightcone inputs

| C++ | Python | Notes |
|---|---|---|
| sky arrays in radians, flat `[object][ra, dec, z]` | arrays of shape (N, 3) in `angle_unit`, `"deg"` or `"rad"` | `angle_unit` is required in every lightcone call |
| `otswap::skyToRadians(skyDegrees)` | — | Python converts at the call, from `angle_unit` |
| `tracersSky`, `randomsSky`, `randomsSkySecondPass` | `tracers_sky`, `randoms_sky`, `randoms_sky_second_pass` | |
| `otswap::Mask` | `otswap.Mask`; keyword `mask` | |
| `otswap::DistanceTable` | `otswap.DistanceTable`; keyword `distances` | |
| `otswap::MpsProfile`, argument `mps` | `otswap.MpsProfile`; keyword `mps_profile` | |
| `otswap::BiasTable`, `otswap::io::readBiasTable` | `otswap.BiasTable`, `otswap.io.read_bias_table` | the bias of `LightconeRsd` |

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
| `findVoidsLightcone(tracersSky, displacement, valid, randomsSky, mask, distances, mps, config, redshiftCut, ids, grid)` | `find_voids_lightcone(tracers_sky, displacement, randoms_sky, valid, *, mask, distances, mps_profile, angle_unit, config, redshift_cut, ids, grid)` | sky arrays in `angle_unit`; `randoms_sky` before `valid` |
| `findVoidsLightcone(reconstruction, randomsSky, mask, config, grid)` | `find_voids_lightcone(reconstruction, randoms_sky, *, mask, angle_unit, config, grid)` | an `otswap.Result` as the first argument |
| `findVoidsLightconeFromCartesian(tracers, displacement, valid, randoms, mask, distances, mps, config, redshiftCut, ids, grid)` | `find_voids_lightcone_from_cartesian(tracers, displacement, randoms, valid, *, mask, distances, mps_profile, angle_unit, config, redshift_cut, ids, grid)` | `angle_unit` is that of the sky arrays of the result |
| `lightconeGrid(reconstruction, randomsSky, mask, divergenceCellSize)` | `lightcone_grid(reconstruction, randoms_sky, *, mask, angle_unit, cell_size)` | |
| `lightconeGrid(tracersSky, randomsSky, mask, distances, mps, divergenceCellSize, redshiftCut)` | `lightcone_grid(tracers_sky, randoms_sky, *, mask, distances, mps_profile, cell_size, angle_unit, redshift_cut)` | |
| `runLightcone(tracersSky, randomsSky, mask, distances, config, randomsSkySecondPass, ids)` | `run_lightcone(tracers_sky, randoms_sky, randoms_sky_second_pass, *, mask, distances, angle_unit, config, ids)` | `config` required |

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
| `BasinCatalog::delta2Mps` | `delta_2mps` | |
| `BasinCatalog::centreSky`, `DivergenceGrid::centreSky`, flat `[row][ra, dec, z]` | `centre_sky`, shape (n, 3) | angles in `angle_unit`; shape (0, 3) in a box |
| `Result::ids` | `ids` | a new copy on each access |
| `Result::randomSeed` | `random_seed` | `None` when no randoms were drawn |
| `Result::geometry` (`otswap::Geometry`) | `"box"`, `"lightcone"` | |
| angles in a lightcone result or grid (radians) | the `angle_unit` of the call that made it, `"deg"` or `"rad"` | `angle_unit` is `None` for a box |
| per-cell arrays of a lightcone grid | flat, `(n_cells,)` and `(n_cells, 3)` | |
| `DivergenceGrid::distanceMin`, `distanceMax`, `raMin`, `raMax`, `decMin`, `decMax` | `distance_min`, `distance_max`, `ra_min`, `ra_max`, `dec_min`, `dec_max` | NaN in a box |
| `DivergenceGrid::cellRaMin`, `cellRaMax`, `cellDecMin`, `cellDecMax`, `cellDistanceMin`, `cellDistanceMax` | `cell_ra_min`, `cell_ra_max`, `cell_dec_min`, `cell_dec_max`, `cell_distance_min`, `cell_distance_max` | shape (0,) in a box and in a grid read from a file |
| `DivergenceGrid::masked` | `masked`, `bool` | |
| `TracerAssignment::inMask`, `inCut` | `in_mask`, `in_cut`, `bool` | shape (0,) in a box |
| `Result::mpsProfile`, `DivergenceGrid::mpsProfile` | `mps_profile` | empty in C++, `None` in Python, for a box |
| `Result::distances` (`std::shared_ptr`) | `distances` | null in C++, `None` in Python, for a box |
| `Result::redshiftCut` | `redshift_cut` | no cut in C++, `None` in Python, for a box |
| `Result::maskNside` | `mask_nside` | 0 in a box |
| `Counts::tracersOutsideMask`, `tracersOutsideRedshiftCut`, `tracersInMaskedCells`, `randomsOutsideMask`, `randomsOutsideRedshiftCut`, `randomsOutsideGrid` | `tracers_outside_mask`, `tracers_outside_redshift_cut`, `tracers_in_masked_cells`, `randoms_outside_mask`, `randoms_outside_redshift_cut`, `randoms_outside_grid` | |
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
