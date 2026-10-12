.. _sec:gmxqmmm:

Gromacs/CPMD QM/MM Calculations
===============================

As of version 3.3 the Gromacs\ [gmx3]_ classical MD code contains
a generic API for QM/MM calculations. So far this API has been used to QM/MM
interfaces for GAMESS, Gaussian, and …CPMD. Unlike in the Gromos/CPMD QM/MM
interface code (see above) the main MD driver is in the classical code. The
Gromacs/CPMD interface code is based on the EGO/CPMD interface and was
developed by Pradip Kumar Biswas in the group of Valentin Gogonea at Cleveland
State University. For additional information and downloads see
http://comppsi.csuohio.edu/groups/qmmm.html, and the respective
publication [gmxqmmm]_.

Technical Introduction
----------------------

The whole interface code is divided into two parts: one part is embedded in the
Gromacs code and the other in CPMD. Both, the modified Gromacs and the CPMD
codes are compiled independently and communicate via files in the current
working directory.

Since the Gromacs code acts as the driver, you first have to set up a regular
Gromacs classical MD simulation in the usual way by building/providing a
.pdf/.gro file and a .top file. Before running grompp, you also need to create
an index file (usually named index.ndx) that lists the atoms of the QM
subsystem and provide further parameters for the CPMD calculation like the size
of QM-simulation box, plane-wave cutoff for CPMD, Coulomb cutoff, if any, etc
(for details, see the rgmx script in the QM/MM examples).

During mdrun, the interface is controlled by two function calls:

a) init_cpmd() prepares the ground for the QM/MM interface. It sets the flags
for the QM and MM atoms finds LINK atoms from the topology and prepares
temporary structures to process the QM/MM data etc.

b) call_cpmd() first creates the CPMD input file "CPMD_inp.run" using a
template "CPMD_inp.tmpl" and then kickstarts the CPMD code via a "fork/exec" or
"system" call. The interface is set to use "fork". If system call is preferred,
you need to set the defined variable NOFORK to 1. call_cpmd() gets forces and
energy from CPMD and appends them to Gromacs structures. Gromacs then moves the
atoms and while evaluating the forces, calls CPMD again (this is the QM/MM
loop). Thus this interface essentially performs a QM/MM Born-Oppenheimer MD
simulation.

Compilation of Gromacs
----------------------

In its present state, you need to use CFLAGS = -DGMX_QMMM_CPMD in configure to
include the CPMD interface code into Gromacs. In the adapted Gromacs package, a
script "build" is provided in the gromacs folder that takes care of the QM/MM
configuration and compilation. Please adapt as needed.

Execution of QM/MM runs
-----------------------

It is like running Gromacs with the additional needs are given by:

a) having an index.ndx file specifying the QM atoms (see example index.ndx
file). You can create a usual Gromacs index.ndx file and then append to it the
QM group.

b) specifying other QM information like planewave cutoff, qmbox size etc in
the grompp setup (for the mdp file).

c) having a CPMD input file template "CPMD_inp.tmpl" where essential keywords
for CPMD run need to be mentioned. "INTERFACE GMX" is essential for QMMM; it
ensures a single-point calculation inside CPMD each time it is invoked. Inside
the interface, all the QM & MM atoms are translated in such a way that the QM
system be at the center of the QM box. Thus the keyword "MOLECULE CENTER OFF"
is required to avoid any further movements of the QM atoms.

Right now the :kw:`ODIIS` minimizer and :kw:`PCG` minimizer
(including PCG MINIMIZE) are allowed to be used inside CPMD. There also is a
hybrid scheme where for the MD first step it will use the "PCG MINIMIZE" but
for all subsequent steps it will use the faster ODIIS minimizer.). Other
sections of CPMD input structures need to be kept as usual though the final
values for the CELL size and CUTOFF will be those provided by you in the mdp
file.

QM/MM Examples
--------------

Example inputs for a H\ :math:`{}_2`\ O-dimer and an ethane molecule are
bundled with the modified gromacs distribution from.

http://comppsi.csuohio.edu/groups/qmmm.html,
