# Builds the example against an installed bitvf, from this folder alone.
# Configure with -DCMAKE_PREFIX_PATH=<the prefix bitvf was installed to>,
# and CFITSIO_ROOT (and ZLIB_ROOT for a static cfitsio) where cfitsio is not
# in a system path:
#
#   cmake -S examples/full_recons_lightcone -B build/examples/full_recons_lightcone -DCMAKE_PREFIX_PATH=<prefix>
#   cmake --build build/examples/full_recons_lightcone
#   build/examples/full_recons_lightcone/full_recons_lightcone [data_dir]
#
# The program reads the examples' data folder, unless data_dir is given,
# and writes to the folder output/ next to this file.

cmake_minimum_required(VERSION 3.16)

project(bitvf_example_full_recons_lightcone LANGUAGES CXX)

find_package(bitvf 0.1 CONFIG REQUIRED)

add_executable(full_recons_lightcone full_recons_lightcone.cpp)
target_link_libraries(full_recons_lightcone PRIVATE bitvf::bitvf)
target_compile_definitions(full_recons_lightcone PRIVATE
  BITVF_EXAMPLES_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../data"
  BITVF_EXAMPLE_OUTPUT="${CMAKE_CURRENT_SOURCE_DIR}/output")
