# Run cpmd.x on one or more inputs and compare the last total energy with a
# reference value.
#
# Variables: LAUNCHER (MPI launcher, may be empty), CPMD, POSTFLAGS,
#            INPUTS ("|"-separated, run in order), PP_DIR, WORK_DIR, ENERGY,
#            TOLERANCE

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")

string(REPLACE "|" ";" inputs "${INPUTS}")
list(LENGTH inputs ninputs)
set(step 0)
foreach(input IN LISTS inputs)
  math(EXPR step "${step} + 1")
  if(ninputs GREATER 1)
    set(out "output.${step}")
    set(err "error.${step}")
  else()
    set(out "output")
    set(err "error")
  endif()
  execute_process(
    COMMAND ${LAUNCHER} "${CPMD}" ${POSTFLAGS} "${input}" "${PP_DIR}"
    WORKING_DIRECTORY "${WORK_DIR}"
    OUTPUT_FILE "${WORK_DIR}/${out}"
    ERROR_FILE "${WORK_DIR}/${err}"
    RESULT_VARIABLE result)
  file(READ "${WORK_DIR}/${out}" output)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "cpmd.x failed (${result}) on ${input}, see ${WORK_DIR}:\n${output}")
  endif()
  if(NOT output MATCHES "PROGRAM CPMD ENDED")
    message(FATAL_ERROR "cpmd.x did not finish on ${input}, see ${WORK_DIR}/${out}")
  endif()
endforeach()

string(REGEX MATCHALL "TOTAL ENERGY = +[-0-9.]+" lines "${output}")
if(NOT lines)
  message(FATAL_ERROR "No total energy found in ${WORK_DIR}/${out}")
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
