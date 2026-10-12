# Optional external modules of CPMD.
#
# Each block mirrors one switch of configure.sh:
#   -qmmm       Gromos + MM_Interface   (libgromos.a, libinterface.a)
#   -iphigenie  IPhigenie_Interface     (libiffiinter.a)
#   -scr        scratchmodule_lib       (libscratch_module.a)
#   -vdw        vdw_lib                 (libgrimme.a)

# cpmd_module_library(<target> <output name> <description> <sources>...)
#
# Static library from optional module sources: stops with a clear message
# if sources are missing, keeps its Fortran modules in mod/<target> and
# exports that directory to its users.
function(cpmd_module_library target output description)
  set(_missing)
  foreach(_f IN LISTS ARGN)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_f}")
      list(APPEND _missing "${_f}")
    endif()
  endforeach()
  if(_missing)
    list(JOIN _missing "\n  " _missing)
    message(FATAL_ERROR "${description}: missing source files in ${CMAKE_CURRENT_SOURCE_DIR}:\n  ${_missing}")
  endif()

  set(_moddir "${CMAKE_CURRENT_BINARY_DIR}/mod/${target}")
  add_library(${target} STATIC ${ARGN})
  set_target_properties(${target} PROPERTIES
    OUTPUT_NAME ${output}
    Fortran_MODULE_DIRECTORY "${_moddir}")
  target_include_directories(${target} PUBLIC "${_moddir}")
endfunction()

# ---------------------------------------------------------------------------
# QM/MM: Gromos and MM_Interface (modules/QMMM_SOURCES)
# ---------------------------------------------------------------------------
if(CPMD_ENABLE_QMMM)
  # GROMOS_MODULES
  set(CPMD_GROMOS_MODULES
    Gromos/coordsz.mod.F90
    Gromos/coordar.mod.F90
  )
  # GROMOS_SRC
  set(CPMD_GROMOS_SRC
    Gromos/blockio.F
    Gromos/c06xxx.F
    Gromos/cobond.F
    Gromos/conat.F
    Gromos/fileio.F
    Gromos/gauss.F
    Gromos/latsum.F
    Gromos/mdutils.F
    Gromos/mm_concg.F
    Gromos/mm_for.F
    Gromos/mm_force.F
    Gromos/mm_setup.F
    Gromos/mm_setup_dr.F
    Gromos/mm_shake.F
    Gromos/nbpml.F
    Gromos/nonbml.F
    Gromos/posio.F
    Gromos/random.F
    Gromos/rdmd.F
    Gromos/mm_restx.F
    Gromos/rdtopo.F
    Gromos/schifezza.F
    Gromos/string.F
    Gromos/traco.F
    Gromos/wrtopo.F
    Gromos/allocate_gromos.F
  )
  # INTERFACE_SRC
  set(CPMD_INTERFACE_SRC
    MM_Interface/mm_add_hydrogen.F
    MM_Interface/mm_center.F
    MM_Interface/mm_reflex.F
    MM_Interface/mm_charges.F
    MM_Interface/mm_detit.F
    MM_Interface/mm_elstat.F
    MM_Interface/mm_elstat_atdens.F
    MM_Interface/mm_elstat_sr.F
    MM_Interface/mm_excl.F
    MM_Interface/mm_excl_mech.F
    MM_Interface/mm_flex_solv.F
    MM_Interface/mm_gcp.F
    MM_Interface/mm_get_NSX.F
    MM_Interface/mm_invert.F
    MM_Interface/mm_long_range_classic.F
    MM_Interface/mm_long_range_potential.F
    MM_Interface/mm_min_im.F
    MM_Interface/mm_multipole.F
    MM_Interface/mm_nlist.F
    MM_Interface/mm_pardef.F
    MM_Interface/mm_quantum_topo.F
    MM_Interface/mm_readgromos.F
    MM_Interface/mm_short_range_classic.F
    MM_Interface/mm_solv_constr.F
    MM_Interface/mm_transl_c0.F
    MM_Interface/mm_write_gromos_coord.F
    MM_Interface/mm_read_qmmm_input.F
    MM_Interface/mm_write_potential.F
    MM_Interface/mm_restart_traj.F
    MM_Interface/mm_add_dummy.F
    MM_Interface/mm_cap_H.F
    MM_Interface/mm_charge_rest.F
  )
  # configure.sh: CPPFLAGS_GROMOS
  set(_gromos_definitions EWALD EWATCUT HAT_SHAPE UNPACKED_GRID)
  if("__HASNT_OMP_COLLAPSE" IN_LIST CPMD_DEFINITIONS)
    list(APPEND _gromos_definitions __HASNT_OMP_COLLAPSE)
  endif()
  set(_gromos_includes
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/Gromos"
    "${CMAKE_CURRENT_SOURCE_DIR}/MM_Interface")
  # FFLAGS_GROMOS[_MODULES] of the configuration, on top of FFLAGS
  separate_arguments(_fflags_gromos UNIX_COMMAND "${CPMD_CFG_FFLAGS_GROMOS}")
  separate_arguments(_fflags_gromos_modules UNIX_COMMAND "${CPMD_CFG_FFLAGS_GROMOS_MODULES}")
  set_source_files_properties(${CPMD_GROMOS_SRC} ${CPMD_INTERFACE_SRC}
    PROPERTIES Fortran_FORMAT FIXED COMPILE_OPTIONS "${_fflags_gromos}")
  set_source_files_properties(${CPMD_GROMOS_MODULES}
    PROPERTIES COMPILE_OPTIONS "${_fflags_gromos_modules}")

  set(_missing_gromos "QM/MM (the Gromos sources are not distributed with CPMD, see Gromos/README.md)")

  # Gromos is compiled with the Fortran flags, but not with the CPMD
  # preprocessor definitions (as in the generated Makefile)
  cpmd_module_library(cpmd_gromos gromos "${_missing_gromos}"
    ${CPMD_GROMOS_MODULES} ${CPMD_GROMOS_SRC})
  target_compile_definitions(cpmd_gromos PRIVATE ${_gromos_definitions})
  target_include_directories(cpmd_gromos PRIVATE ${_gromos_includes})
  target_link_libraries(cpmd_gromos PUBLIC cpmd_options)

  # The MM interface uses the CPMD modules; libcpmd calls back into it.
  cpmd_module_library(cpmd_mm_interface interface "${_missing_gromos}"
    ${CPMD_INTERFACE_SRC})
  target_compile_definitions(cpmd_mm_interface PRIVATE ${_gromos_definitions})
  target_include_directories(cpmd_mm_interface PRIVATE ${_gromos_includes})
  target_link_libraries(cpmd_mm_interface PRIVATE cpmd_lib cpmd_gromos)

  # Gromos headers included through gromos.h by the CPMD sources
  target_include_directories(cpmd_definitions INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# ---------------------------------------------------------------------------
# IPHIGENIE interface (IPhigenie_Interface/IFFIINTER_SOURCES)
# ---------------------------------------------------------------------------
if(CPMD_ENABLE_IPHIGENIE)
  set(CPMD_IFFIINTER_SRC
    IPhigenie_Interface/iffi_comm.mod.F90
    IPhigenie_Interface/iffi_types.mod.F90
    IPhigenie_Interface/iffi_mgmt.mod.F90
    IPhigenie_Interface/iffi_elstat.mod.F90
    IPhigenie_Interface/iffi_inter.mod.F90
    IPhigenie_Interface/iffi_wfct.mod.F90
  )
  # libiffiinter.a is linked into IPHIGENIE together with libcpmd.a
  cpmd_module_library(cpmd_iffiinter iffiinter "IPHIGENIE interface"
    ${CPMD_IFFIINTER_SRC})
  target_sources(cpmd_iffiinter PRIVATE "${PROJECT_SOURCE_DIR}/src/timetag.F90")
  target_include_directories(cpmd_iffiinter PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/IPhigenie_Interface")
  target_link_libraries(cpmd_iffiinter PUBLIC cpmd_lib)
  install(TARGETS cpmd_iffiinter ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

# ---------------------------------------------------------------------------
# Scratch library (scratchmodule_lib/SCRATCHLIB_SOURCES)
# ---------------------------------------------------------------------------
if(CPMD_ENABLE_SCRATCHLIB)
  set(CPMD_SCRATCHMODULE_SRC
    scratchmodule_lib/data_managment_utils.F90
    scratchmodule_lib/segment_managment_utils.F90
    scratchmodule_lib/pool_managment_utils.F90
    scratchmodule_lib/scratch_interface.F90
  )
  cpmd_module_library(cpmd_scratch scratch_module "Scratch library"
    ${CPMD_SCRATCHMODULE_SRC})
  target_include_directories(cpmd_scratch PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/scratchmodule_lib")
  target_link_libraries(cpmd_scratch PUBLIC cpmd_options cpmd_definitions)
endif()

# ---------------------------------------------------------------------------
# Grimme vdW library (vdw_lib/VDWLIB_SOURCES)
# ---------------------------------------------------------------------------
if(CPMD_ENABLE_VDWLIB)
  set(CPMD_GRIMME_SRC
    vdw_lib/vdw_calculator.F90
    vdw_lib/vdw_interface.F90
    vdw_lib/vdw_param.F90
  )
  cpmd_module_library(cpmd_grimme grimme "Grimme vdW library"
    ${CPMD_GRIMME_SRC})
  target_include_directories(cpmd_grimme PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/vdw_lib")
  target_link_libraries(cpmd_grimme PUBLIC cpmd_options cpmd_definitions)
endif()
