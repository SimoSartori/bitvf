# Builds the example against an installed bitvf, from this folder alone.
# Configure with -DCMAKE_PREFIX_PATH=<the prefix bitvf was installed to>,
# and CFITSIO_ROOT (and ZLIB_ROOT for a static cfitsio) where cfitsio is not
# in a system path:
#
#   cmake -S examples/field -B build/examples/field -DCMAKE_PREFIX_PATH=<prefix>
#   cmake --build build/examples/field
#   build/examples/field/field [field_file]
#
# The program reads the divergence file the full_recons_box example writes,
# so run that example first, unless field_file is given; it writes to the
# folder output/ next to this file.

cmake_minimum_required(VERSION 3.16)

project(bitvf_example_field LANGUAGES CXX)

find_package(bitvf 0.1 CONFIG REQUIRED)

add_executable(field field.cpp)
target_link_libraries(field PRIVATE bitvf::bitvf)
target_compile_definitions(field PRIVATE
  BITVF_EXAMPLE_FIELD="${CMAKE_CURRENT_SOURCE_DIR}/../full_recons_box/output/divergence_box.dat"
  BITVF_EXAMPLE_OUTPUT="${CMAKE_CURRENT_SOURCE_DIR}/output")
