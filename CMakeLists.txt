cmake_minimum_required(VERSION 3.16)

project(bitvf VERSION 0.1.0 LANGUAGES C CXX)

# PROJECT_IS_TOP_LEVEL needs CMake 3.21; this is the 3.16 equivalent, true
# for a standalone build and false under add_subdirectory or FetchContent.
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(BITVF_IS_TOP_LEVEL ON)
else()
  set(BITVF_IS_TOP_LEVEL OFF)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Only a top-level build may choose these for the whole build tree; as a
# subproject, the parent's choices stand.
if (BITVF_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  endif()
endif()

# scikit-build-core defines SKBUILD. The wheel carries the Python module only,
# so the C++ install rules are off there.
if (DEFINED SKBUILD)
  set(BITVF_INSTALL_DEFAULT OFF)
else()
  set(BITVF_INSTALL_DEFAULT ON)
endif()

option(BITVF_BUILD_PYTHON "Build the Python bindings" OFF)
option(BITVF_BUILD_TESTS "Build the C++ test suite" ${BITVF_IS_TOP_LEVEL})
option(BITVF_INSTALL "Generate the C++ install rules" ${BITVF_INSTALL_DEFAULT})

# The executables bitvf and bitvf-field: built by a standalone build, not by
# the wheel nor as a subproject.
if (BITVF_IS_TOP_LEVEL AND NOT DEFINED SKBUILD)
  set(BITVF_BUILD_EXECUTABLE_DEFAULT ON)
else()
  set(BITVF_BUILD_EXECUTABLE_DEFAULT OFF)
endif()
option(BITVF_BUILD_EXECUTABLE "Build the bitvf and bitvf-field executables" ${BITVF_BUILD_EXECUTABLE_DEFAULT})

include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# Floating-point flags
#
# No a*b+c is fused into one FMA anywhere: compilers differ in what they
# contract by default, and so do the instruction sets, so contraction would
# make the results depend on the platform. otswap and meshsearch set it on
# their own targets; this covers bitvf's and a fetched otswap.
#
# Maths functions do not set errno, so that the compiler may expand the ones
# that map to an instruction (llrint, sqrt) inline instead of calling the C
# library; GCC sets errno by default, clang on macOS does not. The results
# are the same either way. This covers bitvf's targets and a fetched otswap
# and meshsearch.
#
# Neither is exported: all of bitvf's arithmetic is in its own archive.
# ---------------------------------------------------------------------------

if (NOT MSVC)
  add_compile_options(-ffp-contract=off -fno-math-errno)
endif()

# ---------------------------------------------------------------------------
# OpenMP
#
# AppleClang needs to be pointed at a libomp. In a Python build (scikit-build-core),
# the otswap package it builds against is imported with the interpreter that
# configures the build, and if it bundles a libomp (its macOS wheels do,
# otswap.openmp_root()) that one is used, so that bitvf and otswap load a
# single runtime. Other builds never look at a Python package: they use
# Homebrew's libomp, as the otswap they build or find does. OpenMP_ROOT, as a
# CMake variable or in the environment, takes precedence over both.
# bitvfConfig.cmake repeats the Homebrew lookup for consumers.
# ---------------------------------------------------------------------------

if (SKBUILD AND Python_EXECUTABLE AND
    NOT DEFINED OpenMP_ROOT AND NOT DEFINED ENV{OpenMP_ROOT})
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c
            "import otswap; print(otswap.openmp_root() or '')"
    OUTPUT_VARIABLE bitvf_omp_root OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE bitvf_omp_status ERROR_QUIET)
  if (bitvf_omp_status EQUAL 0 AND bitvf_omp_root AND
      EXISTS "${bitvf_omp_root}/include/omp.h")
    set(OpenMP_ROOT "${bitvf_omp_root}")
    message(STATUS "bitvf: OpenMP from the otswap package, ${OpenMP_ROOT}")
  endif()
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND
    NOT DEFINED OpenMP_ROOT AND NOT DEFINED ENV{OpenMP_ROOT})
  foreach (candidate /opt/homebrew/opt/libomp /usr/local/opt/libomp)
    if (EXISTS ${candidate}/include/omp.h)
      set(OpenMP_ROOT ${candidate})
      break()
    endif()
  endforeach()
endif()

find_package(OpenMP COMPONENTS CXX)

# ---------------------------------------------------------------------------
# otswap
#
# An installed otswap 0.3.x is used when found, through
# CMAKE_PREFIX_PATH or otswap_DIR; otherwise the pinned tag is fetched and
# built from source, and brings meshsearch the same way. A parent project that already defines
# otswap::otswap provides its own.
# ---------------------------------------------------------------------------

set(BITVF_OTSWAP_VERSION 0.3)
set(BITVF_OTSWAP_TAG "v0.3.0" CACHE STRING
    "otswap tag to fetch when no installed otswap is found")
set(BITVF_OTSWAP_REPOSITORY "https://github.com/SimoSartori/otswap.git"
    CACHE STRING "otswap repository to fetch from; may be a local path")

if (NOT TARGET otswap::otswap)
  find_package(otswap ${BITVF_OTSWAP_VERSION} CONFIG QUIET)
endif()

if (TARGET otswap::otswap)
  if (otswap_FOUND)
    message(STATUS "bitvf: using otswap ${otswap_VERSION} from ${otswap_DIR}")
  else()
    message(STATUS "bitvf: using the parent project's otswap::otswap")
  endif()
else()
  message(STATUS "bitvf: no installed otswap ${BITVF_OTSWAP_VERSION}; "
                 "fetching ${BITVF_OTSWAP_TAG}")

  # otswap's install rules follow bitvf's, so that an installed bitvf comes
  # with the otswap it was built against, and the wheel installs neither. A
  # normal variable overrides otswap's option (CMP0077).
  set(OTSWAP_INSTALL ${BITVF_INSTALL})

  include(FetchContent)
  FetchContent_Declare(otswap
    GIT_REPOSITORY ${BITVF_OTSWAP_REPOSITORY}
    GIT_TAG        ${BITVF_OTSWAP_TAG}
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(otswap)
endif()

# ---------------------------------------------------------------------------
# The library
# ---------------------------------------------------------------------------

add_library(bitvf STATIC
  src/Basins.cpp
  src/BitVFBox.cpp
  src/BitVFLightcone.cpp
  src/BoxGrid.cpp
  src/Checks.cpp
  src/Divergence.cpp
  src/LightconeBasins.cpp
  src/LightconeCentres.cpp
  src/LightconeDensities.cpp
  src/LightconeDivergence.cpp
  src/LightconeGrid.cpp
  src/LightconeSmoothing.cpp
  src/LightconeWatershed.cpp
  src/Readers.cpp
  src/Smoothing.cpp
  src/Support.cpp
  src/VoidFinderBox.cpp
  src/VoidFinderField.cpp
  src/VoidFinderLightcone.cpp
  src/Watershed.cpp
  src/Writers.cpp)
add_library(bitvf::bitvf ALIAS bitvf)

set_target_properties(bitvf PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(bitvf
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_features(bitvf PUBLIC cxx_std_17)

# The version the writers record in the files they write.
target_compile_definitions(bitvf PRIVATE BITVF_VERSION="${PROJECT_VERSION}")

# otswap is public: bitvf's interface takes and returns otswap's types.
# OpenMP is private: a consumer of the static archive gets its link
# requirement only.
target_link_libraries(bitvf PUBLIC otswap::otswap)

set(BITVF_WITH_OPENMP OFF)
if (OpenMP_CXX_FOUND)
  target_link_libraries(bitvf PRIVATE OpenMP::OpenMP_CXX)
  set(BITVF_WITH_OPENMP ON)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(bitvf PRIVATE -Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# The executables
#
# app/bitvf.cpp and app/bitvf_field.cpp, on the public headers only, with
# app/common.h. The targets are bitvf_app and bitvf_field_app, since bitvf
# is the library's; the files they build are bitvf and bitvf-field.
# ---------------------------------------------------------------------------

if (BITVF_BUILD_EXECUTABLE)
  foreach (app bitvf bitvf_field)
    add_executable(${app}_app app/${app}.cpp)
    target_compile_definitions(${app}_app PRIVATE BITVF_VERSION="${PROJECT_VERSION}")
    target_link_libraries(${app}_app PRIVATE bitvf)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
      target_compile_options(${app}_app PRIVATE -Wall -Wextra)
    endif()
  endforeach()
  set_target_properties(bitvf_app PROPERTIES OUTPUT_NAME bitvf)
  set_target_properties(bitvf_field_app PROPERTIES OUTPUT_NAME bitvf-field)
endif()

# ---------------------------------------------------------------------------
# Tests
#
# One executable per file of tests/, each a CTest test. They include src/
# to reach the internal functions they check.
# ---------------------------------------------------------------------------

if (BITVF_BUILD_TESTS)
  include(CTest)
  enable_testing()

  foreach (test_name test_checks test_box_grid test_divergence test_smoothing test_watershed test_basins test_writers test_readers
                    test_from_divergence test_from_field test_report test_runbox test_lightcone_entry
                    test_lightcone_grid test_lightcone_faces test_lightcone_divergence test_lightcone_fill
                    test_lightcone_smoothing test_lightcone_watershed test_lightcone_basins test_lightcone_centres
                    test_lightcone_densities test_lightcone_pipeline)
    add_executable(${test_name} tests/${test_name}.cpp)
    target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(${test_name} PRIVATE bitvf)
    if (OpenMP_CXX_FOUND)
      target_link_libraries(${test_name} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
      target_compile_options(${test_name} PRIVATE -Wall -Wextra)
    endif()
    add_test(NAME ${test_name} COMMAND ${test_name})
  endforeach()

  # Runs the library on raw input files for the Python tests, which compare
  # its outputs and its report with those of the bindings; not a test by
  # itself. tests/test_bitvf.py finds it through BITVF_PYTHON_REFERENCE.
  add_executable(python_reference tests/python_reference.cpp)
  target_link_libraries(python_reference PRIVATE bitvf)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(python_reference PRIVATE -Wall -Wextra)
  endif()

  # The executables: one test per parameter file of tests/executable/cases,
  # run by run_case.cmake in a fresh directory of its own, which checks
  # what the file's "#!" lines expect. executable_reference first writes,
  # with the library calls of each kind of run, the files the runs must
  # reproduce byte for byte. The two shipped parameter files run from the
  # source directory, as their comments say.
  if (BITVF_BUILD_EXECUTABLE)
    set(BITVF_EXECUTABLE_REFERENCE ${CMAKE_CURRENT_BINARY_DIR}/executable/reference)
    add_executable(executable_reference tests/executable/reference.cpp)
    target_link_libraries(executable_reference PRIVATE bitvf)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
      target_compile_options(executable_reference PRIVATE -Wall -Wextra)
    endif()
    add_test(NAME executable_reference
      COMMAND executable_reference ${CMAKE_CURRENT_SOURCE_DIR}/tests/executable/data ${BITVF_EXECUTABLE_REFERENCE})
    set_tests_properties(executable_reference PROPERTIES FIXTURES_SETUP executable_reference)

    file(GLOB bitvf_executable_cases CONFIGURE_DEPENDS
         ${CMAKE_CURRENT_SOURCE_DIR}/tests/executable/cases/*.ini)
    foreach (case_file ${bitvf_executable_cases})
      get_filename_component(case_name ${case_file} NAME_WE)
      add_test(NAME executable_${case_name}
        COMMAND ${CMAKE_COMMAND}
          -DBITVF=$<TARGET_FILE:bitvf_app>
          -DBITVF_FIELD=$<TARGET_FILE:bitvf_field_app>
          -DCASE=${case_file}
          -DWORK=${CMAKE_CURRENT_BINARY_DIR}/executable/${case_name}
          -DTEST_DATA=${CMAKE_CURRENT_SOURCE_DIR}/tests/executable/data
          -DREFERENCE=${BITVF_EXECUTABLE_REFERENCE}
          -DPROJECT_VERSION=${PROJECT_VERSION}
          -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/executable/run_case.cmake)
      set_tests_properties(executable_${case_name} PROPERTIES FIXTURES_REQUIRED executable_reference)
    endforeach()

    # The commented parameter files of the documentation, run as shipped:
    # field.ini reads the divergence file box.ini writes.
    add_test(NAME executable_example_box COMMAND bitvf_app examples/executable/box.ini
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(executable_example_box PROPERTIES FIXTURES_SETUP executable_example_box)
    add_test(NAME executable_example_field COMMAND bitvf_field_app examples/executable/field.ini
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    set_tests_properties(executable_example_field PROPERTIES FIXTURES_REQUIRED executable_example_box)
  endif()
endif()

# ---------------------------------------------------------------------------
# Python bindings
#
# The extension module bitvf._bitvf, installed into the bitvf package of
# python/bitvf. Built by pip through scikit-build-core (see pyproject.toml).
# It is built as otswap's own module is, static nanobind and no NB_DOMAIN,
# so that both share nanobind's type registry and bitvf can take and return
# otswap's Python types. src/ is on its include path for the version.
# ---------------------------------------------------------------------------

if (BITVF_BUILD_PYTHON)
  find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
  find_package(nanobind CONFIG REQUIRED)

  nanobind_add_module(bitvf_python NB_STATIC python/bitvf_python.cpp)
  set_target_properties(bitvf_python PROPERTIES OUTPUT_NAME _bitvf)
  target_include_directories(bitvf_python PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(bitvf_python PRIVATE bitvf)
  if (APPLE)
    set_target_properties(bitvf_python PROPERTIES
      INSTALL_RPATH "@loader_path/../otswap/openmp/lib")
  endif()
  install(TARGETS bitvf_python LIBRARY DESTINATION bitvf)
endif()

# ---------------------------------------------------------------------------
# Install and export
#
# find_package(bitvf 0.1 CONFIG) then bitvf::bitvf. Its dependencies are
# found again on the consumer's machine by bitvfConfig.cmake: otswap 0.3.x,
# which brings meshsearch and cfitsio, and OpenMP when bitvf
# was built with it.
# ---------------------------------------------------------------------------

if (BITVF_INSTALL)
  include(CMakePackageConfigHelpers)

  set(BITVF_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/bitvf)

  install(TARGETS bitvf EXPORT bitvfTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY include/bitvf DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  if (BITVF_BUILD_EXECUTABLE)
    install(TARGETS bitvf_app bitvf_field_app RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()

  install(EXPORT bitvfTargets
          NAMESPACE bitvf::
          DESTINATION ${BITVF_CMAKE_DIR})

  configure_package_config_file(cmake/bitvfConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/bitvfConfig.cmake
    INSTALL_DESTINATION ${BITVF_CMAKE_DIR})

  # Before 1.0 a minor release may change the interface, so a request for 0.1
  # accepts any 0.1.x and nothing else.
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/bitvfConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/bitvfConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/bitvfConfigVersion.cmake
          DESTINATION ${BITVF_CMAKE_DIR})
endif()
