cmake_minimum_required(VERSION 3.1)

project(my_project VERSION 0.1
    DESCRIPTION "CPMD"
    LANGUAGES Fortran C)

#option(BOOTSTRAP "Build bootstrap version of Fpm" OFF)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

add_library(iphigenie
    modules/IPhigenie_Interface/iffi_comm.mod.F90
    modules/IPhigenie_Interface/iffi_elstat.mod.F90
    modules/IPhigenie_Interface/iffi_inter.mod.F90
    modules/IPhigenie_Interface/iffi_mgmt.mod.F90
    modules/IPhigenie_Interface/iffi_types.mod.F90
    modules/IPhigenie_Interface/iffi_wfct.mod.F90)

add_library(mimic
        modules/MiMiC/mimic.F90)

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

add_subdirectory(src)


