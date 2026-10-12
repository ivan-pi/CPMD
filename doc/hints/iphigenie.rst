QM/(P)MM Interface to IPHIGENIE
===============================

The interface of CPMD to the PMM-MD program
IPHIGENIE\ [Schwoerer2013]_ [Schwoerer2015]_
(https://sourceforge.net/projects/iphigenie) supersedes the CPMD interface of
Eichinger et al. [egoqmmm]_ to the MD program
EGO.\ [ego1]_ It provides a Hamiltonian coupling to polarizable MM
force fields (PMM). Here, similar to the Eichinger implementation and the
derived CPMD/Gromacs coupling\ [gmxqmmm]_, CPMD is used only to
compute energies and forces but the propagation of the equations of motion is
done externally by IPHIGENIE.

The implementation is based on a single executable, which links CPMD and the
interface routines as libraries to the IPHIGENIE MD executable ’iffi’. For
CPMD-3.17 a patch is provided on `www.cpmd.org <www.cpmd.org>`__ with example
configurations for the compilation. For CPMD-4.0 and later the interface to
IPHIGENIE is provided as a module. Here, using the configure.sh script with the
option ``-iphigenie`` generates a makefile which builds the interface library
along with a regular cpmd.x executable.

For further instructions on compiling and example inputs for QM/PMM runs please
consult the main iphigenie page https://sourceforge.net/projects/iphigenie and
the associated Wiki pages.
