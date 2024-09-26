if(VDW)
    message(STATUS "VDW enabled")
    add_definitions(-D__HAS_LIBGRIMMEVDW)
    add_library(vdw
        vdw_lib/vdw_calculator.F90
        vdw_lib/vdw_interface.F90
        vdw_lib/vdw_param.F90)
endif()

if(SCR)
    message(STATUS "Scratch library enabled")
    add_definitions(-D__USE_SCRATCH_LIBRARY)
    add_library(scratchmodule
        scratchmodule_lib/data_managment_utils.F90
        scratchmodule_lib/segment_managment_utils.F90
        scratchmodule_lib/pool_managment_utils.F90
        scratchmodule_lib/scratch_interface.F90)
endif()


add_library(mimic
        MiMiC/mimic.F90)


if (IPHIGENIE)
    message(STATUS "Building interface to iphigenie")
    add_library(iffiinter
        IPhigenie_Interface/iffi_comm.mod.F90
        IPhigenie_Interface/iffi_elstat.mod.F90
        IPhigenie_Interface/iffi_inter.mod.F90
        IPhigenie_Interface/iffi_mgmt.mod.F90
        IPhigenie_Interface/iffi_types.mod.F90
        IPhigenie_Interface/iffi_wfct.mod.F90)
    target_include_directories(iffiinter PRIVATE src)
    target_link_libraries(iffiinter PRIVATE cpmd)
endif()


if(QMMM)
    message("Building QMMM interface")
    # TODO: Just GLOB these
    add_library(mm_interface
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
endif()
