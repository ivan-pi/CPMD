if(VDW)
    add_library(vdw
        modules/vdw_lib/vdw_calculator.F90
        modules/vdw_lib/vdw_grimme.F90
        modules/vdw_lib/vdw_interface.F90
        modules/vdw_lib/vdw_param-D3mod.inc
        modules/vdw_lib/vdw_param-D3org.inc
        modules/vdw_lib/vdw_param-R0ab.inc
        modules/vdw_lib/vdw_param.F90)
endif()


if(SCR)
    add_library(scratchmodule
        modules/scratchmodule_lib/data_managment_utils.F90
        modules/scratchmodule_lib/segment_managment_utils.F90
        modules/scratchmodule_lib/pool_managment_utils.F90
        modules/scratchmodule_lib/scratch_interface.F90)
endif()


add_library(mimic
        modules/MiMiC/mimic.F90)


if (IPHIGENIE)
    message("Building interface to iphigenie")
    add_library(iffiinter
        modules/IPhigenie_Interface/iffi_comm.mod.F90
        modules/IPhigenie_Interface/iffi_elstat.mod.F90
        modules/IPhigenie_Interface/iffi_inter.mod.F90
        modules/IPhigenie_Interface/iffi_mgmt.mod.F90
        modules/IPhigenie_Interface/iffi_types.mod.F90
        modules/IPhigenie_Interface/iffi_wfct.mod.F90)
    target_include_directories(iffiinter PRIVATE src)
    target_link_libraries(iffiinter PRIVATE cpmd)
endif()


if(QMMM)
    message("Building QMMM interface")
    # TODO: Just GLOB these
    add_library(mm_interface
        modules/MM_Interface/mm_add_hydrogen.F
        modules/MM_Interface/mm_center.F
        modules/MM_Interface/mm_reflex.F
        modules/MM_Interface/mm_charges.F
        modules/MM_Interface/mm_detit.F
        modules/MM_Interface/mm_elstat.F
        modules/MM_Interface/mm_elstat_atdens.F
        modules/MM_Interface/mm_elstat_sr.F
        modules/MM_Interface/mm_excl.F
        modules/MM_Interface/mm_excl_mech.F
        modules/MM_Interface/mm_flex_solv.F
        modules/MM_Interface/mm_gcp.F
        modules/MM_Interface/mm_get_NSX.F
        modules/MM_Interface/mm_invert.F
        modules/MM_Interface/mm_long_range_classic.F
        modules/MM_Interface/mm_long_range_potential.F
        modules/MM_Interface/mm_min_im.F
        modules/MM_Interface/mm_multipole.F
        modules/MM_Interface/mm_nlist.F
        modules/MM_Interface/mm_pardef.F
        modules/MM_Interface/mm_quantum_topo.F
        modules/MM_Interface/mm_readgromos.F
        modules/MM_Interface/mm_short_range_classic.F
        modules/MM_Interface/mm_solv_constr.F
        modules/MM_Interface/mm_transl_c0.F
        modules/MM_Interface/mm_write_gromos_coord.F
        modules/MM_Interface/mm_read_qmmm_input.F
        modules/MM_Interface/mm_write_potential.F 
        modules/MM_Interface/mm_restart_traj.F
        modules/MM_Interface/mm_add_dummy.F
        modules/MM_Interface/mm_cap_H.F
        modules/MM_Interface/mm_charge_rest.F
        )
endif()
