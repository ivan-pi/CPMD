# Run cpmd.x and compare the last total energy with a reference value.
#
# Variables: LAUNCHER (MPI launcher, may be empty), CPMD, INPUT, PP_DIR,
#            WORK_DIR, ENERGY, TOLERANCE

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

execute_process(
  COMMAND ${LAUNCHER} "${CPMD}" "${INPUT}" "${PP_DIR}"
  WORKING_DIRECTORY "${WORK_DIR}"
  OUTPUT_FILE "${WORK_DIR}/output"
  ERROR_FILE "${WORK_DIR}/error"
  RESULT_VARIABLE result)
file(READ "${WORK_DIR}/output" output)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "cpmd.x failed (${result}), see ${WORK_DIR}:\n${output}")
endif()
if(NOT output MATCHES "PROGRAM CPMD ENDED")
  message(FATAL_ERROR "cpmd.x did not finish, see ${WORK_DIR}/output")
endif()

string(REGEX MATCHALL "TOTAL ENERGY = +[-0-9.]+" lines "${output}")
if(NOT lines)
  message(FATAL_ERROR "No total energy found in ${WORK_DIR}/output")
endif()
list(POP_BACK lines last)
string(REGEX REPLACE ".* " "" energy "${last}")

# CMake has no floating point arithmetic: compare in units of 1e-8 Hartree
function(to_units value out)
  if(NOT value MATCHES "^(-?)([0-9]*)\\.?([0-9]*)$")
    message(FATAL_ERROR "Not a number: ${value}")
  endif()
  set(_sign "${CMAKE_MATCH_1}")
  set(_int "${CMAKE_MATCH_2}")
  string(SUBSTRING "${CMAKE_MATCH_3}00000000" 0 8 _frac)
  math(EXPR _units "${_sign}(${_int}0 / 10 * 100000000 + 1${_frac} - 100000000)")
  set(${out} ${_units} PARENT_SCOPE)
endfunction()
to_units("${energy}" energy_u)
to_units("${ENERGY}" reference_u)
to_units("${TOLERANCE}" tolerance_u)
math(EXPR diff "${energy_u} - ${reference_u}")
if(diff LESS 0)
  math(EXPR diff "-(${diff})")
endif()

message(STATUS "Total energy ${energy}, reference ${ENERGY}")
if(diff GREATER tolerance_u)
  message(FATAL_ERROR "Total energy ${energy} differs from reference ${ENERGY} "
                      "by more than ${TOLERANCE}")
endif()
