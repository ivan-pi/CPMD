# Tests for the CMake build, run with ctest.
#
# Each test runs cpmd.x on one or more small inputs in its own directory
# below the build tree and compares the last total energy with a reference
# value. The smoke test is in this directory; the regression tests run by
# the CI workflows are in .github/tests (see .github/tests/tests.cmake).

if(CPMD_ENABLE_MPI)
  # Already done without a configuration; with one it provides MPIEXEC_*
  find_package(MPI QUIET)
  set(CPMD_TEST_NPROCS 2 CACHE STRING "Number of MPI processes for the tests")
  set(_launcher ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${CPMD_TEST_NPROCS}
                ${MPIEXEC_PREFLAGS})
endif()
set(CPMD_TEST_NTHREADS 1 CACHE STRING "Number of OpenMP threads per process for the tests")

# cpmd_add_test(<name> <input>... ENERGY <reference> TOLERANCE <tolerance>
#               [PP_DIR <directory>])
#
# The inputs (relative to the calling directory) are run one after the
# other in the same work directory, so later ones can restart from the
# files written by earlier ones. The last total energy printed by the last
# run is compared with the reference. Pseudopotentials are read from
# PP_DIR, by default the directory of the first input.
function(cpmd_add_test name)
  cmake_parse_arguments(PARSE_ARGV 1 _arg "" "ENERGY;TOLERANCE;PP_DIR" "")
  set(_inputs ${_arg_UNPARSED_ARGUMENTS})
  if(NOT _inputs OR NOT DEFINED _arg_ENERGY OR NOT DEFINED _arg_TOLERANCE)
    message(FATAL_ERROR "cpmd_add_test(${name}): inputs, ENERGY and TOLERANCE are required")
  endif()
  list(TRANSFORM _inputs PREPEND "${CMAKE_CURRENT_LIST_DIR}/")
  if(DEFINED _arg_PP_DIR)
    get_filename_component(_pp_dir "${_arg_PP_DIR}" ABSOLUTE
                           BASE_DIR "${CMAKE_CURRENT_LIST_DIR}")
  else()
    list(GET _inputs 0 _first)
    get_filename_component(_pp_dir "${_first}" DIRECTORY)
  endif()
  # Lists cannot be passed through -D unescaped
  string(REPLACE ";" "|" _inputs "${_inputs}")
  add_test(NAME ${name}
    COMMAND ${CMAKE_COMMAND}
      "-DLAUNCHER=${_launcher}"
      "-DCPMD=$<TARGET_FILE:cpmd>"
      "-DPOSTFLAGS=${MPIEXEC_POSTFLAGS}"
      "-DINPUTS=${_inputs}"
      "-DPP_DIR=${_pp_dir}"
      "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${name}"
      "-DENERGY=${_arg_ENERGY}"
      "-DTOLERANCE=${_arg_TOLERANCE}"
      -P "${PROJECT_SOURCE_DIR}/tests/run_cpmd.cmake")
  set_tests_properties(${name} PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=${CPMD_TEST_NTHREADS}"
    TIMEOUT 600)
endfunction()

# H2 molecule, LDA, wavefunction optimisation
cpmd_add_test(h2_wfopt h2/h2.inp ENERGY -1.11731768 TOLERANCE 0.000001)

include("${PROJECT_SOURCE_DIR}/.github/tests/tests.cmake")
