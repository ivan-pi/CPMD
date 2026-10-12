# Compile and link settings shared by all CPMD targets.
#
# They are collected on two INTERFACE targets:
#   cpmd_options      compiler flags and link libraries (FFLAGS, CFLAGS, LFLAGS)
#   cpmd_definitions  preprocessor definitions and include paths (CPPFLAGS)
# The Gromos library only uses cpmd_options, everything else uses both.
#
# The settings come either from a machine configuration (CPMD_CONFIG, see
# cmake/CPMDImportConfig.cmake) or, without one, from CMake's own compiler
# detection and find_package().

add_library(cpmd_options INTERFACE)
add_library(cpmd_definitions INTERFACE)

# Definitions selected by the CMake options, used in both modes
set(CPMD_DEFINITIONS)
if(CPMD_ENABLE_QMMM)
  list(APPEND CPMD_DEFINITIONS __GROMOS)
endif()
if(CPMD_ENABLE_IPHIGENIE)
  list(APPEND CPMD_DEFINITIONS __IFFIINTER)
endif()
if(CPMD_ENABLE_SCRATCHLIB)
  list(APPEND CPMD_DEFINITIONS __USE_SCRATCHLIBRARY)
endif()
if(CPMD_ENABLE_VDWLIB)
  list(APPEND CPMD_DEFINITIONS __HAS_LIBGRIMMEVDW)
endif()
if(CPMD_MINPACK_LIBRARY)
  list(APPEND CPMD_DEFINITIONS __MINPACK)
  target_link_libraries(cpmd_options INTERFACE "${CPMD_MINPACK_LIBRARY}")
endif()

if(CPMD_CONFIG)
  # -------------------------------------------------------------------------
  # Settings from configure/<CPMD_CONFIG>
  # -------------------------------------------------------------------------
  # CPPFLAGS: -D -> definitions, -I -> include directories, rest -> options
  separate_arguments(_cppflags UNIX_COMMAND "${CPMD_CFG_CPPFLAGS}")
  set(_defs ${_cppflags})
  list(FILTER _defs INCLUDE REGEX "^-D")
  list(TRANSFORM _defs REPLACE "^-D" "")
  set(_incs ${_cppflags})
  list(FILTER _incs INCLUDE REGEX "^-I")
  list(TRANSFORM _incs REPLACE "^-I" "")
  list(FILTER _cppflags EXCLUDE REGEX "^-[DI]")
  list(APPEND CPMD_DEFINITIONS ${_defs})
  target_include_directories(cpmd_definitions INTERFACE ${_incs})
  target_compile_options(cpmd_definitions INTERFACE ${_cppflags})

  separate_arguments(_fflags UNIX_COMMAND "${CPMD_CFG_FFLAGS}")
  separate_arguments(_cflags UNIX_COMMAND "${CPMD_CFG_CFLAGS}")
  separate_arguments(_lflags UNIX_COMMAND "${CPMD_CFG_LFLAGS}")
  list(PREPEND _fflags ${CPMD_CFG_FC_EXTRA_FLAGS})
  list(PREPEND _cflags ${CPMD_CFG_CC_EXTRA_FLAGS})
  target_compile_options(cpmd_options INTERFACE
    "$<$<COMPILE_LANGUAGE:Fortran>:${_fflags}>"
    "$<$<COMPILE_LANGUAGE:C>:${_cflags}>")
  # The generated Makefile links with "$(LD) $(FFLAGS) ... $(LFLAGS)"
  target_link_options(cpmd_options INTERFACE ${_fflags})
  target_link_libraries(cpmd_options INTERFACE ${_lflags})

  if(CPMD_CFG_NVCC)
    separate_arguments(_nvccflags UNIX_COMMAND "${CPMD_CFG_NVCCFLAGS}")
    list(PREPEND _nvccflags ${CPMD_CFG_NVCC_EXTRA_FLAGS})
    target_compile_options(cpmd_options INTERFACE
      "$<$<COMPILE_LANGUAGE:CUDA>:${_nvccflags}>")
  endif()

  if(CPMD_ENABLE_OPENMP AND CPMD_CFG_OMP3_DISABLED STREQUAL "true")
    list(APPEND CPMD_DEFINITIONS __HASNT_OMP_COLLAPSE)
  endif()
else()
  # -------------------------------------------------------------------------
  # Native CMake configuration
  # -------------------------------------------------------------------------
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    list(APPEND CPMD_DEFINITIONS __Linux)
  elseif(APPLE)
    list(APPEND CPMD_DEFINITIONS __OSX)
  endif()

  set(_fflags)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    list(APPEND CPMD_DEFINITIONS __GNU)
    list(APPEND _fflags -ffree-line-length-none -falign-commons)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      list(APPEND _fflags -fallow-argument-mismatch)
    endif()
  endif()
  target_compile_options(cpmd_options INTERFACE
    "$<$<COMPILE_LANGUAGE:Fortran>:${_fflags}>")

  if(CPMD_ENABLE_COVERAGE)
    if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      message(FATAL_ERROR "CPMD_ENABLE_COVERAGE requires the GNU compilers")
    endif()
    target_compile_options(cpmd_options INTERFACE --coverage)
    target_link_options(cpmd_options INTERFACE --coverage)
  endif()

  if(CPMD_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    if(NOT MPI_Fortran_HAVE_F08_MODULE)
      message(FATAL_ERROR "CPMD requires an MPI library with the mpi_f08 module")
    endif()
    list(APPEND CPMD_DEFINITIONS __PARALLEL)
    target_link_libraries(cpmd_options INTERFACE MPI::MPI_Fortran)
  endif()

  if(CPMD_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran C)
    target_link_libraries(cpmd_options INTERFACE
      OpenMP::OpenMP_Fortran OpenMP::OpenMP_C)
  endif()

  if(CPMD_FFT STREQUAL "FFTW3")
    find_library(CPMD_FFTW3_LIBRARY NAMES fftw3 REQUIRED
                 HINTS ENV FFTW_ROOT ENV FFTW_HOME PATH_SUFFIXES lib lib64)
    if(CPMD_ENABLE_OPENMP)
      # dfftw_plan_with_nthreads() is called in OpenMP builds
      find_library(CPMD_FFTW3_OMP_LIBRARY NAMES fftw3_omp REQUIRED
                   HINTS ENV FFTW_ROOT ENV FFTW_HOME PATH_SUFFIXES lib lib64)
    endif()
    list(APPEND CPMD_DEFINITIONS __HAS_FFT_FFTW3)
    target_link_libraries(cpmd_options INTERFACE
      ${CPMD_FFTW3_OMP_LIBRARY} ${CPMD_FFTW3_LIBRARY})
  else()
    list(APPEND CPMD_DEFINITIONS __HAS_FFT_DEFAULT)
  endif()

  if(CPMD_ENABLE_LIBXC)
    find_package(Libxc CONFIG QUIET)
    if(TARGET Libxc::xcf03)
      target_link_libraries(cpmd_options INTERFACE Libxc::xcf03)
    else()
      find_package(PkgConfig REQUIRED)
      pkg_check_modules(LIBXCF03 REQUIRED IMPORTED_TARGET libxcf03)
      target_link_libraries(cpmd_options INTERFACE PkgConfig::LIBXCF03)
    endif()
    list(APPEND CPMD_DEFINITIONS __HAS_LIBXC)
  endif()

  if(CPMD_ENABLE_CUDA)
    find_package(CUDAToolkit REQUIRED)
    list(APPEND CPMD_DEFINITIONS __HAS_CUDA)
    target_link_libraries(cpmd_options INTERFACE
      CUDA::cudart CUDA::cufft CUDA::cublas CUDA::cusolver)
    if(TARGET CUDA::nvToolsExt)
      list(APPEND CPMD_DEFINITIONS __HAS_NVTX_TIMER)
      target_link_libraries(cpmd_options INTERFACE CUDA::nvToolsExt)
    endif()
  endif()

  find_package(LAPACK REQUIRED)
  target_link_libraries(cpmd_options INTERFACE LAPACK::LAPACK)
endif()

list(REMOVE_DUPLICATES CPMD_DEFINITIONS)
target_compile_definitions(cpmd_definitions INTERFACE ${CPMD_DEFINITIONS})
target_include_directories(cpmd_definitions INTERFACE
  "${PROJECT_SOURCE_DIR}/src"
  "${PROJECT_BINARY_DIR}/include")
