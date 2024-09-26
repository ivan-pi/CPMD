cmake_minimum_required(VERSION 3.1)

project(my_project VERSION 0.1
    DESCRIPTION "CPMD"
    LANGUAGES Fortran C)

#
# CONFIGURATION PART
#

# Modules
option(QMMM "QMMM module" OFF)
option(SCR "Use scratch library" ON)
option(VDW "Van-der-Waals library" ON)
option(IPHIGENIE "Support for external interface to iphigenie" OFF)

# Threading
option(OMP "OpenMP support" ON)
option(DISABLE_OMP3 "Disable OMP3 instructions" OFF)

# Optional librares
option(MINPACK "Use MINPACK library" OFF)

#--------------------------------------------------------------

#
# BUILD PART
#
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_definitions(-D__GIT_REV="cmake")

add_definitions(-D__HAS_FFT_FFTW3)
add_definitions(-D__ALPHALINUX)
add_definitions(-D__HPC)
add_definitions(-D__PARALLEL)
add_definitions(-D__HAS_SIZEOF)
add_definitions(-D__HAS_BF_STREAM_IO)
add_definitions(-D__INTEL_MKL)
add_definitions(-D__USE_SCRATCH_LIBRARY)
add_definitions(-D__HAS_LIBGRIMMEVDW)

if(SCR)
    add_definitions(-D__USE_SCRATCHLIBRARY)
endif()

add_subdirectory(modules)

# libcpmd.a and cpmd.x are defined in here
add_subdirectory(src)