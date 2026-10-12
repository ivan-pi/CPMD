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

# Split a configure-style CPPFLAGS string into definitions, include
# directories and other options.
function(cpmd_split_cppflags flags out_defs out_incs out_opts)
  separate_arguments(_words UNIX_COMMAND "${flags}")
  set(_defs)
  set(_incs)
  set(_opts)
  foreach(_w IN LISTS _words)
    if(_w MATCHES "^-D(.+)$")
      list(APPEND _defs "${CMAKE_MATCH_1}")
    elseif(_w MATCHES "^-I(.+)$")
      list(APPEND _incs "${CMAKE_MATCH_1}")
    else()
      list(APPEND _opts "${_w}")
    endif()
  endforeach()
  set(${out_defs} "${_defs}" PARENT_SCOPE)
  set(${out_incs} "${_incs}" PARENT_SCOPE)
  set(${out_opts} "${_opts}" PARENT_SCOPE)
endfunction()

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
  cpmd_split_cppflags("${CPMD_CFG_CPPFLAGS}" _defs _incs _opts)
  list(APPEND CPMD_DEFINITIONS ${_defs})
  target_include_directories(cpmd_definitions INTERFACE ${_incs})
  target_compile_options(cpmd_definitions INTERFACE ${_opts})

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

  if("__PARALLEL" IN_LIST CPMD_DEFINITIONS)
    set(CPMD_ENABLE_MPI ON CACHE BOOL "Build with MPI" FORCE)
  else()
    set(CPMD_ENABLE_MPI OFF CACHE BOOL "Build with MPI" FORCE)
  endif()
  if(CPMD_CFG_NVCC AND "__HAS_CUDA" IN_LIST CPMD_DEFINITIONS)
    set(CPMD_ENABLE_CUDA ON CACHE BOOL "Build with CUDA" FORCE)
  else()
    set(CPMD_ENABLE_CUDA OFF CACHE BOOL "Build with CUDA" FORCE)
  endif()

  # Sources the Makefile compiles with NOOPT_FLAG instead of FFLAGS. Here
  # NOOPT_FLAG is appended, so it has to override the optimisation level.
  set(CPMD_NOOPT_SOURCES)
  separate_arguments(_objs UNIX_COMMAND "${CPMD_CFG_NOOPT_OBJS}")
  foreach(_o IN LISTS _objs)
    string(REGEX REPLACE "\\.o$" ".F90" _src "${_o}")
    list(APPEND CPMD_NOOPT_SOURCES "${_src}")
  endforeach()
  separate_arguments(CPMD_NOOPT_FLAGS UNIX_COMMAND "${CPMD_CFG_NOOPT_FLAG}")

  separate_arguments(CPMD_GROMOS_FFLAGS UNIX_COMMAND "${CPMD_CFG_FFLAGS_GROMOS}")
  separate_arguments(CPMD_GROMOS_MODULES_FFLAGS UNIX_COMMAND
                     "${CPMD_CFG_FFLAGS_GROMOS_MODULES}")
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
