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

include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# Floating-point contraction
#
# No a*b+c is fused into one FMA anywhere: compilers differ in what they
# contract by default, and so do the instruction sets, so contraction would
# make the results depend on the platform. otswap and meshsearch set it on
# their own targets; this covers bitvf's and a fetched otswap. It is not
# exported: all of bitvf's arithmetic is in its own archive.
# ---------------------------------------------------------------------------

if (NOT MSVC)
  add_compile_options(-ffp-contract=off)
endif()

# ---------------------------------------------------------------------------
# OpenMP
#
# AppleClang needs to be pointed at a libomp; Homebrew's is the usual one,
# used unless OpenMP_ROOT is set, as a CMake variable or in the environment.
# otswap looks for it the same way, so both libraries use the same runtime.
# bitvfConfig.cmake repeats this for consumers.
# ---------------------------------------------------------------------------

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
# An installed otswap 0.2.x is used when found, through CMAKE_PREFIX_PATH or
# otswap_DIR; otherwise the pinned tag is fetched and built from source, and
# brings meshsearch the same way. A parent project that already defines
# otswap::otswap provides its own.
# ---------------------------------------------------------------------------

set(BITVF_OTSWAP_VERSION 0.2)
set(BITVF_OTSWAP_TAG "v0.2.0" CACHE STRING
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
  src/BoxGrid.cpp
  src/Checks.cpp
  src/Support.cpp
  src/VoidFinderBox.cpp)
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
# Tests
#
# One executable per file of tests/, each a CTest test. They include src/
# to reach the internal functions they check.
# ---------------------------------------------------------------------------

if (BITVF_BUILD_TESTS)
  include(CTest)
  enable_testing()

  foreach (test_name test_checks test_box_grid)
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
  install(TARGETS bitvf_python LIBRARY DESTINATION bitvf)
endif()

# ---------------------------------------------------------------------------
# Install and export
#
# find_package(bitvf 0.1 CONFIG) then bitvf::bitvf. Its dependencies are
# found again on the consumer's machine by bitvfConfig.cmake: otswap 0.2.x,
# which brings meshsearch and cfitsio, and OpenMP when bitvf was built with
# it.
# ---------------------------------------------------------------------------

if (BITVF_INSTALL)
  include(CMakePackageConfigHelpers)

  set(BITVF_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/bitvf)

  install(TARGETS bitvf EXPORT bitvfTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY include/bitvf DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
