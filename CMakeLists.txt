cmake_minimum_required(VERSION 3.1)

project(my_project VERSION 0.1
    DESCRIPTION "CPMD"
    LANGUAGES Fortran C)


#  -D__HAS_FFT_DEFAULT Use default FFT routine in the code (from Goedecker)#
#  -D__HAS_FFT_ESSL    Use FFT ESSL (need ESSL library)                    #
#  -D__HAS_FFT_HP      Use FFT HP (Hewlet-Packard)                         #
#  -D__HAS_FFT_FFTW3   Use FFTW FFT routines                               #
#  -D__PARALLEL        For PARALLEL computers, you need to have a          #
#                      MPI parallel library installed                      #
#  -D__VECTOR          For vectorial computers                             #
#  -D__NOINT8          If the compiler does not use INTEGER*8 or           #
#                      if the function BTEST and IBTEST do not like it     #
#  -D__NO_MEMSIZE      to omit the detection of the current program size  

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
add_definitions(-DSVN_REV="cmake")
add_definitions(-D__GIT_REV="cmake")

add_subdirectory(src)

add_library(iffiinter
    modules/IPhigenie_Interface/iffi_comm.mod.F90
    modules/IPhigenie_Interface/iffi_elstat.mod.F90
    modules/IPhigenie_Interface/iffi_inter.mod.F90
    modules/IPhigenie_Interface/iffi_mgmt.mod.F90
    modules/IPhigenie_Interface/iffi_types.mod.F90
    modules/IPhigenie_Interface/iffi_wfct.mod.F90)
target_include_directories(iffiinter PRIVATE src)
target_link_libraries(iffiinter PRIVATE cpmd)

add_library(mimic
        modules/MiMiC/mimic.F90)

# add_library(mm_interface
#     modules/MM_Interface/mm_add_hydrogen.F
#     modules/MM_Interface/mm_center.F
#     modules/MM_Interface/mm_reflex.F
#     modules/MM_Interface/mm_charges.F
#     modules/MM_Interface/mm_detit.F
#     modules/MM_Interface/mm_elstat.F
#     modules/MM_Interface/mm_elstat_atdens.F
#     modules/MM_Interface/mm_elstat_sr.F
#     modules/MM_Interface/mm_excl.F
#     modules/MM_Interface/mm_excl_mech.F
#     modules/MM_Interface/mm_flex_solv.F
#     modules/MM_Interface/mm_gcp.F
#     modules/MM_Interface/mm_get_NSX.F
#     modules/MM_Interface/mm_invert.F
#     modules/MM_Interface/mm_long_range_classic.F
#     modules/MM_Interface/mm_long_range_potential.F
#     modules/MM_Interface/mm_min_im.F
#     modules/MM_Interface/mm_multipole.F
#     modules/MM_Interface/mm_nlist.F
#     modules/MM_Interface/mm_pardef.F
#     modules/MM_Interface/mm_quantum_topo.F
#     modules/MM_Interface/mm_readgromos.F
#     modules/MM_Interface/mm_short_range_classic.F
#     modules/MM_Interface/mm_solv_constr.F
#     modules/MM_Interface/mm_transl_c0.F
#     modules/MM_Interface/mm_write_gromos_coord.F
#     modules/MM_Interface/mm_read_qmmm_input.F
#     modules/MM_Interface/mm_write_potential.F 
#     modules/MM_Interface/mm_restart_traj.F
#     modules/MM_Interface/mm_add_dummy.F
#     modules/MM_Interface/mm_cap_H.F
#     modules/MM_Interface/mm_charge_rest.F
#     )

add_library(vdw
    modules/vdw_lib/vdw_calculator.F90
    modules/vdw_lib/vdw_grimme.F90
    modules/vdw_lib/vdw_interface.F90
    modules/vdw_lib/vdw_param-D3mod.inc
    modules/vdw_lib/vdw_param-D3org.inc
    modules/vdw_lib/vdw_param-R0ab.inc
    modules/vdw_lib/vdw_param.F90)

add_library(scratchmodule
    modules/scratchmodule_lib/data_managment_utils.F90
    modules/scratchmodule_lib/segment_managment_utils.F90
    modules/scratchmodule_lib/pool_managment_utils.F90
    modules/scratchmodule_lib/scratch_interface.F90)