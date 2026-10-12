# Smoke tests for the CMake build, run with ctest.
#
# Each test runs cpmd.x on a small input in its own directory below the
# build tree and compares the final total energy with a reference value.

if(CPMD_ENABLE_MPI)
  # Already done without a configuration; with one it provides MPIEXEC_*
  find_package(MPI QUIET)
  set(CPMD_TEST_NPROCS 2 CACHE STRING "Number of MPI processes for the tests")
  set(_launcher ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${CPMD_TEST_NPROCS}
                ${MPIEXEC_PREFLAGS})
endif()

# cpmd_add_test(<name> <input> <reference energy> <tolerance>)
function(cpmd_add_test name input energy tolerance)
  get_filename_component(_dir "${CMAKE_CURRENT_SOURCE_DIR}/${input}" DIRECTORY)
  add_test(NAME ${name}
    COMMAND ${CMAKE_COMMAND}
      "-DLAUNCHER=${_launcher}"
      "-DCPMD=$<TARGET_FILE:cpmd>"
      "-DPOSTFLAGS=${MPIEXEC_POSTFLAGS}"
      "-DINPUT=${CMAKE_CURRENT_SOURCE_DIR}/${input}"
      "-DPP_DIR=${_dir}"
      "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${name}"
      "-DENERGY=${energy}"
      "-DTOLERANCE=${tolerance}"
      -P "${CMAKE_CURRENT_SOURCE_DIR}/run_cpmd.cmake")
  set_tests_properties(${name} PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1"
    TIMEOUT 300)
endfunction()

# H2 molecule, LDA, wavefunction optimisation
cpmd_add_test(h2_wfopt h2/h2.inp -1.11731768 0.000001)
