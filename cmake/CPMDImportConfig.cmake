# Import a machine configuration from configure/<NAME> (the files used by
# configure.sh) into the CMake build.
#
# This file is included *before* project(), because a configuration selects
# the compilers. The configuration file is evaluated by
# scripts/cmake/import_config.sh with the same switches as configure.sh
# (-omp, -debug, -coverage, -qmmm), and the variables it defines are made
# available as CPMD_CFG_<VAR>. They are turned into compile options,
# definitions and link libraries in cmake/CPMDFlags.cmake.
#
# Mapping of configuration variables:
#   FC, CC, NVCC     -> CMAKE_Fortran_COMPILER, CMAKE_C_COMPILER,
#                       CMAKE_CUDA_COMPILER (only on the first configure and
#                       only if not given explicitly)
#   FFLAGS           -> Fortran compile options (also used when linking,
#                       as the generated Makefile does)
#   CFLAGS           -> C compile options
#   NVCCFLAGS        -> CUDA compile options
#   CPPFLAGS         -> preprocessor definitions / include directories
#   LFLAGS           -> link libraries and flags of cpmd.x
#   FFLAGS_GROMOS*   -> extra Fortran options for the QM/MM modules
#   NOOPT_OBJS/FLAG  -> extra options for the listed sources
#   OMP3_DISABLED    -> -D__HASNT_OMP_COLLAPSE with OpenMP enabled
#   CPP, LD, AR, RANLIB, IRAT are not imported: CMake preprocesses with the
#   compiler and links with the Fortran compiler.

set(_cpmd_cfg_args)
if(CPMD_ENABLE_OPENMP)
  list(APPEND _cpmd_cfg_args -omp)
endif()
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  list(APPEND _cpmd_cfg_args -debug)
endif()
if(CPMD_ENABLE_COVERAGE)
  list(APPEND _cpmd_cfg_args -coverage)
endif()
if(CPMD_ENABLE_QMMM)
  list(APPEND _cpmd_cfg_args -qmmm)
endif()

find_program(CPMD_BASH bash REQUIRED)
mark_as_advanced(CPMD_BASH)
execute_process(
  COMMAND "${CPMD_BASH}" "${CMAKE_CURRENT_LIST_DIR}/../scripts/cmake/import_config.sh"
          ${_cpmd_cfg_args} "${CPMD_CONFIG}"
  OUTPUT_VARIABLE _cpmd_cfg_out
  ERROR_VARIABLE _cpmd_cfg_err
  RESULT_VARIABLE _cpmd_cfg_res)
if(NOT _cpmd_cfg_res EQUAL 0)
  message(FATAL_ERROR
    "Could not import configuration '${CPMD_CONFIG}' (${_cpmd_cfg_args}):\n"
    "${_cpmd_cfg_err}")
endif()
if(_cpmd_cfg_err)
  message(STATUS "Configuration ${CPMD_CONFIG} says:\n${_cpmd_cfg_err}")
endif()

cmake_language(EVAL CODE "${_cpmd_cfg_out}")
message(STATUS "Using CPMD configuration ${CPMD_CFG_FILE}")

# Features fixed by the configuration; these normal variables take
# precedence over the cache options of the same name.
if(" ${CPMD_CFG_CPPFLAGS} " MATCHES " -D__PARALLEL ")
  set(CPMD_ENABLE_MPI ON)
else()
  set(CPMD_ENABLE_MPI OFF)
endif()
if(CPMD_CFG_NVCC AND " ${CPMD_CFG_CPPFLAGS} " MATCHES " -D__HAS_CUDA ")
  set(CPMD_ENABLE_CUDA ON)
  # architectures are chosen with -gencode in NVCCFLAGS
  set(CMAKE_CUDA_ARCHITECTURES OFF)
else()
  set(CPMD_ENABLE_CUDA OFF)
endif()

# Compilers: the first word is the compiler, the rest are extra flags.
# They can only be set before the first project() call of a build tree.
function(_cpmd_cfg_compiler lang var)
  if(NOT CPMD_CFG_${var})
    return()
  endif()
  separate_arguments(_words UNIX_COMMAND "${CPMD_CFG_${var}}")
  list(POP_FRONT _words _prog)
  set(CPMD_CFG_${var}_EXTRA_FLAGS ${_words} PARENT_SCOPE)
  if(NOT CMAKE_${lang}_COMPILER)
    set(CMAKE_${lang}_COMPILER "${_prog}" PARENT_SCOPE)
  elseif(NOT CMAKE_${lang}_COMPILER MATCHES "(^|/)${_prog}$")
    message(STATUS "${lang} compiler is ${CMAKE_${lang}_COMPILER}; "
                   "${var}=${_prog} from ${CPMD_CONFIG} is not used "
                   "(compilers can only be chosen in a fresh build tree)")
  endif()
endfunction()
_cpmd_cfg_compiler(Fortran FC)
_cpmd_cfg_compiler(C CC)
_cpmd_cfg_compiler(CUDA NVCC)
