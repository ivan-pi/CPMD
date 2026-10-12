# Regression tests run by the CI workflows (included from tests/CMakeLists.txt).
#
# Each directory holds the input(s) of one test; the pseudopotentials are
# in pseudo/. The reference energies were computed with gfortran 13 and
# OpenMPI on 2 MPI processes and agree between the built-in FFT, FFTW3,
# OpenMP and the scratch library to well within the tolerances. They test
# that the code paths keep working and give the same results; the small
# cell and cutoff make them unsuitable as physical reference values.
#
# Add a test with
#   cpmd_add_test(<name> <dir>/<input>... ENERGY <hartree> TOLERANCE <hartree>
#                 PP_DIR pseudo)

set(_pp PP_DIR pseudo)

# Wavefunction optimisation
cpmd_add_test(h2o_wfopt_odiis h2o_wfopt_odiis/input.inp
  ENERGY -16.68988805 TOLERANCE 0.000001 ${_pp})
cpmd_add_test(h2o_wfopt_pcg h2o_wfopt_pcg/input.inp
  ENERGY -16.68988805 TOLERANCE 0.000001 ${_pp})
cpmd_add_test(h2o_wfopt_lsd h2o_wfopt_lsd/input.inp
  ENERGY -16.32965709 TOLERANCE 0.000001 ${_pp})
cpmd_add_test(h2o_wfopt_pbe h2o_wfopt_pbe/input.inp
  ENERGY -16.78776004 TOLERANCE 0.000001 ${_pp})
cpmd_add_test(h2_kpoints h2_kpoints/input.inp
  ENERGY -1.12385552 TOLERANCE 0.000001 ${_pp})
# CP groups split the MPI processes into 2 groups
if(CPMD_ENABLE_MPI)
  math(EXPR _odd "${CPMD_TEST_NPROCS} % 2")
endif()
if(CPMD_ENABLE_MPI AND NOT _odd)
  cpmd_add_test(h2o_wfopt_cp_groups h2o_wfopt_cp_groups/input.inp
    ENERGY -16.68988805 TOLERANCE 0.000001 ${_pp})
endif()

# Restart from the files of a previous run
cpmd_add_test(h2o_restart h2o_restart/step1.inp h2o_restart/step2.inp
  ENERGY -16.68988805 TOLERANCE 0.000001 ${_pp})

# Geometry optimisation and molecular dynamics
cpmd_add_test(h2o_geoopt h2o_geoopt/input.inp
  ENERGY -16.69188215 TOLERANCE 0.000001 ${_pp})
cpmd_add_test(h2o_bomd h2o_bomd/input.inp
  ENERGY -16.69174211 TOLERANCE 0.00001 ${_pp})
cpmd_add_test(h2o_cpmd h2o_cpmd/input.inp
  ENERGY -16.69160812 TOLERANCE 0.00001 ${_pp})

# Optional libraries
if(CPMD_ENABLE_VDWLIB)
  cpmd_add_test(h2o_dimer_vdw_d3 h2o_dimer_vdw_d3/input.inp
    ENERGY -33.13303085 TOLERANCE 0.000001 ${_pp})
endif()
