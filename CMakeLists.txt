#
# CMake build for CPMD
#
# Quick start (native configuration, MPI + LAPACK found by CMake):
#
#   cmake -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j
#
# Reusing a machine configuration from the configure/ directory:
#
#   cmake -B build -DCPMD_CONFIG=LINUX-X86_64-GFORTRAN-MPI -DCPMD_ENABLE_OPENMP=ON
#   CPMD_CONFIG=LINUX-X86_64-GFORTRAN-MPI cmake -B build      # same, via environment
#
# See doc/cmake.md for all options.
#
cmake_minimum_required(VERSION 3.18)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" CPMD_VERSION LIMIT_COUNT 1)

# ---------------------------------------------------------------------------
# Options (defined before project() as a configuration may pick compilers)
# ---------------------------------------------------------------------------
set(CPMD_CONFIG "$ENV{CPMD_CONFIG}" CACHE STRING
    "Machine configuration from configure/ (name or path) providing compilers, flags and libraries")

option(CPMD_ENABLE_MPI        "Build the parallel (MPI) version, -D__PARALLEL" ON)
option(CPMD_ENABLE_OPENMP     "Enable OpenMP (configure.sh -omp)" OFF)
option(CPMD_ENABLE_CUDA       "Enable the CUDA code, -D__HAS_CUDA" OFF)
option(CPMD_ENABLE_LIBXC      "Use libxc for exchange-correlation functionals, -D__HAS_LIBXC" OFF)
option(CPMD_ENABLE_QMMM       "Build the QM/MM interface (configure.sh -qmmm); needs the Gromos sources in modules/Gromos" OFF)
option(CPMD_ENABLE_IPHIGENIE  "Build libiffiinter for IPHIGENIE (configure.sh -iphigenie)" OFF)
option(CPMD_ENABLE_SCRATCHLIB "Use the scratch library in modules/scratchmodule_lib (configure.sh -scr)" OFF)
option(CPMD_ENABLE_VDWLIB     "Use the Grimme vdW library in modules/vdw_lib (configure.sh -vdw)" OFF)
option(CPMD_ENABLE_COVERAGE   "Instrument for coverage analysis (configure.sh -coverage)" OFF)
set(CPMD_FFT "DEFAULT" CACHE STRING "FFT library without CPMD_CONFIG: DEFAULT (built-in) or FFTW3")
set_property(CACHE CPMD_FFT PROPERTY STRINGS DEFAULT FFTW3)
set(CPMD_MINPACK_LIBRARY "" CACHE FILEPATH "MINPACK library to link (configure.sh -minpack=), -D__MINPACK")
option(CPMD_BUILD_TESTS       "Add the smoke tests in tests/ (run with ctest)" ON)

if(CPMD_CONFIG)
  include(cmake/CPMDImportConfig.cmake)
elseif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  # A configuration file brings its own optimisation flags; without one,
  # build optimised code unless asked otherwise.
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

project(CPMD
  VERSION ${CPMD_VERSION}
  DESCRIPTION "Car-Parrinello Molecular Dynamics"
  HOMEPAGE_URL "https://github.com/OpenCPMD/CPMD"
  LANGUAGES Fortran C)

if(CPMD_ENABLE_CUDA OR CPMD_CFG_NVCC)
  enable_language(CUDA)
endif()

include(GNUInstallDirs)

# Let the compiler run the preprocessor on all Fortran sources (*.F, *.F90)
set(CMAKE_Fortran_PREPROCESS ON)

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
endif()

include(cmake/CPMDFlags.cmake)

# Version string printed in the output header (scripts/getversion.sh)
execute_process(
  COMMAND "${PROJECT_SOURCE_DIR}/scripts/getversion.sh" "${PROJECT_SOURCE_DIR}"
  OUTPUT_VARIABLE CPMD_GIT_REV
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT CPMD_GIT_REV)
  set(CPMD_GIT_REV "${CPMD_VERSION}")
endif()

# gromos.h is included by the force-matching code; with QM/MM it pulls in
# the Gromos topology headers (see configure.sh).
if(CPMD_ENABLE_QMMM)
  set(_gromos_comment " ")
else()
  set(_gromos_comment "!")
endif()
file(CONFIGURE OUTPUT "${PROJECT_BINARY_DIR}/include/gromos.h" CONTENT
"${_gromos_comment}     include 'Gromos/toposz.h'
${_gromos_comment}     include 'Gromos/topoar.h'
${_gromos_comment}     include 'Gromos/box.h'
")

add_subdirectory(modules)
add_subdirectory(src)

if(CPMD_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "CPMD ${CPMD_GIT_REV}")
if(CPMD_CONFIG)
  message(STATUS "  Configuration     : ${CPMD_CONFIG}")
endif()
message(STATUS "  Fortran compiler  : ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION})")
message(STATUS "  C compiler        : ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION})")
message(STATUS "  Build type        : ${CMAKE_BUILD_TYPE}")
message(STATUS "  MPI / OpenMP      : ${CPMD_ENABLE_MPI} / ${CPMD_ENABLE_OPENMP}")
message(STATUS "  CUDA              : ${CPMD_ENABLE_CUDA}")
message(STATUS "  QM/MM / IPHIGENIE : ${CPMD_ENABLE_QMMM} / ${CPMD_ENABLE_IPHIGENIE}")
message(STATUS "  Definitions       : ${CPMD_DEFINITIONS}")
message(STATUS "")
