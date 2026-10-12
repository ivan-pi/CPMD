.. _intro:

History
=======

CPMD Version 1
--------------

In summer 1993 a project was started to combine the two different ab initio
molecular dynamics codes [Allen87]_ that were used in the group
for computational physics of the IBM Research Laboratory in Rüschlikon. There
was the IBM-AIX version (ported by J. Kohanoff and F. Buda) of the IBM-VM
version (by W. Andreoni and P. Ballone) of the original
Car-Parrinello [CP85]_ code and a version of the code by
K. Laasonen and F. Buda that could handle ultra-soft
pseudopotentials [Vanderbilt]_. Further goals were to provide a
common platform for future developments, as new integration techniques or
parallelization. The original Car–Parrinello code [CP85]_ [Galli91]_
was about 8000 lines of Fortran. A first parallel version using the IBM MPL
library was finished in 1993. Many people contributed to this effort in
different ways: M. Parrinello, J. Hutter, W. Andreoni, A. Curioni,
P. Giannozzi, E. Fois, D. Marx and M. Tuckerman.

CPMD Version 2
--------------

Version 2.0
~~~~~~~~~~~

The first major update of the code was finished in summer 1993. New features of
the code included a keyword driven input, an initial guess from atomic
pseudo-wavefunctions, a module for geometry optimization, several new types of
molecular dynamics, Nosé–Hoover [Nose84]_ [Hoover85]_ thermostats and
a diagonalization routine to get Kohn-Sham energies [KS]_. This
code had 17’000 lines.

Version 2.5
~~~~~~~~~~~

In 1994 many additions were made to the code. The communication was improved
and a library interface for MPI was introduced. The code reached its most
stable version at the end of the year with version number 2.5. At this stage a
working version of *ab initio* path integrals [Marx94]_ [Marx96]_
based on a one level parallelization was implemented in a separate branch of
the code by Dominik Marx.

CPMD Version 3
--------------

Version 3.0
~~~~~~~~~~~

This major update included changes to improve the portability of the code to
other platforms. Most notable was the shmem interface for optimal parallel
performance on Cray computers. New features of this version were constant
pressure molecular dynamics using the Parrinello-Rahman
Lagrangian [parrah]_ [parrah2]_, the possibility for symmetry
constraints and Stefan Goedecker’s dual space
pseudopotentials [goe_pp]_. The library concept for the
pseudopotentials had been changed. The code had grown to 55’000 lines.

Version 3.1
~~~~~~~~~~~

Only minor updates were made for this version. However, it served as a starting
point for two major new developments. The free energy
functional [Alavi94]_ code with **k** points was developed by Ali
Alavi and Thierry Deutsch in Belfast. An efficient path integral version using
two level parallelism was put together by Mark Tuckerman, Dominik Marx, and
Jürg Hutter [Tuckerman96]_ .

Version 3.2
~~~~~~~~~~~

This version included several new algorithms. Some of these were lost in the
transfer to the next version.

Version 3.3
~~~~~~~~~~~

This version was developed using the free energy functional version (based on
3.1) as a basis. The path integral version was fully included but only part of
the changes from the "main" version 3.2 were taken over. The QM/MM interface to
the EGO code was included [egoqmmm]_. Development of the linear
response [LinRes]_ parts of the code started. Maximally localized
Wannier functions [Marzari97]_ were implemented. This version was
finished in 1998, the code was about 115’000 lines long.

Version 3.4
~~~~~~~~~~~

The most notable change to this version was the inclusion of the QM/MM
interface developed by Alessandro Laio, Joost VandeVondele and Ursula
Röthlisberger [qmmm02]_ [qmmm03]_ [qmmm04]_. Besides that only minor
changes to the functionality of the code were done. This version included
mostly bug fixes and was finished in 2000.

Version 3.5
~~~~~~~~~~~

This was the first version made generally available at www.cpmd.org in early
2002. Many bugs were fixed, most notably the code for the ultra-soft
pseudopotentials was working again. The new size of the code was 136’000 lines.

Version 3.6
~~~~~~~~~~~

This developers version included the final versions of the linear response
methods for the calculation of the polarizability and the chemical NMR shifts
developed by Anna Putrino and Daniel
Sebastiani [apdsmp]_ [apmp]_ [dsmp]_. Marcella Iannuzzi contributed a
:math:`{\bf k} \cdot {\bf p}` module [mimp]_. Time-dependent
density functional response theory was implemented and forces for excited state
energies programmed. Salomon Billeter, Alessandro Curioni and Wanda Andreoni
implemented new linear scaling geometry optimizers that allow to locate
geometrical transition states in a clean way [LSCAL]_. Fine
grained parallelism with OpenMP was added (by Alessandro Curioni and Jürg
Hutter) and can be used together with the distributed memory MPI version.

Version 3.7
~~~~~~~~~~~

The stable version of the developers code was made publicly available in early
2003. The code has 150’000 lines.

Version 3.8
~~~~~~~~~~~

Developer’s version.

Version 3.9
~~~~~~~~~~~

Many new developments, improvements, cleanups, and bug fixes have been added
since the last public version of the code. Most notably, the methodology for
reactive Car-Parrinello metadynamics [alaio]_ [iannuzzi]_ is made
available in this version.

Other new functionality includes G-space localization of wavefunctions,
Hockney-type Poisson Solver [Hockney70]_ for slabs with influence
function in G-Space, code to determine molecular KS states from Wannier
functions, code for trajectory analysis, calculation of dipole moments using
the Berry phase and in real space, transition matrix elements between orbitals,
growth function for constraints and restraints, new code for applying static
electrical fields, periodic or final diagonalization of WF, van der Waals force
field according to Elstner’s formula [Elstner]_ and dumping files
for PDOS.

Improvements of the code include performance and OpenMP improvements, improved
code for keeping wavefunction in real space, updated TDDFT, SAOP TDDFT
functional, a much improved configure script, bug fixes for HF exchange,
screened exchange, cleanup of memory management, more checks on unsupported
options, fixed constraints in geometry optimization.

Modified ROKS [Frank98]_, Ports to MacOS-X/PPC, Cray X1, and Intel
EM64T, k-points with swapfiles are working again on many platforms, detection
of incompatible Vanderbilt pseudopotentials.

Version 3.10
~~~~~~~~~~~~

Developer’s version.

Version 3.11
~~~~~~~~~~~~

Many improvements, cleanups, bug fixes and some new features have been added
since the last public version of the code. New functionalities include
calculation of the electric field gradient tensor along MD trajectory, EPR
calculations, efficient wavefunction extrapolation for BOMD, distance screening
for HFX calculation and hybrid functional with PBC, interaction perturbation
method, molecular states in TDDFT calculations, analytic second derivatives of
gradient corrected functionals [xcder]_, Born charge tensor during
finite difference vibrational analysis, Gromacs QM/MM
interface [gmxqmmm]_, and distributed linear algebra support.

New supported platforms include, IBM Blue Gene/L [bgl]_, Cray XT3,
NEC-SX6 Earth Simulator (Vector-Parallel) and Windows NT/XT using GNU Gfortran.
Performance tunings for existing platforms include FFTW interface, 16 Byte
memory, alignment for Blue Gene/L, extension of the taskgroup implementation to
cartesian taskgroups (Blue Gene/L), parallel distributed linear algebra,
alltoall communication in either single (to reduce communication bandwidth) or
double precision, special parallel OPEIGR, improved OpenMP support
[mixed]_, and improved metadynamics.

Version 3.12
~~~~~~~~~~~~

Developer’s version.

Version 3.13
~~~~~~~~~~~~

Several improvements, cleanups, bug fixes and a few new features have been
added since the last public version of the code. New functionalities include
additional distributed linear algebra code for initialization, final
wavefunction projection and Friesner diagonalization, mean free energy path
search method, multiscale shock method [shock]_, Langevin
integrator for metadynamics with extended Lagrangian, calculation of
non-adiabatic couplings, Landau-Zener Surface hopping, ROKS-based Slater
transition-state density, linear-response DFPT with a ROKS-based reference
state [GrimmJCP2003]_, simplified van der Waals correction
according to Grimme [Grimme06]_, simplified ROKS input options
with hard-wired variants of modified Goedecker algorithms for ROKS, PBEsol
functional, ports to IBM Blue Gene/P, MacOS-X/x86 and PACS-CS / T2K, support
for fftw-3, improved ultrasoft pseudopotential parallel code (VDB) (MPI and
OpenMP), optimizations for scalar CPUs, new collective variables for
metadynamics, variable cell support in DCD output, isotropic and zflexible cell
for Parrinello-Rahman dynamics, damped dynamics and Berendsen thermostats for
electrons, ions and cell, path-integral support for BO-MD, support for
completely reproducible outputs for CPMD TestSuite, consistent and updated unit
conversions throughout the code, spin-density Mulliken analysis, aClimax format
output of vibrational frequencies, optimization scheme for Goedecker
pseudopotential parameters for use as link atoms in QM/MM applications, support
for QUENCH BO with PCG MINIMIZE when using VDB potentials, corrections for a
number of serious bugs in the Gromos QM/MM code, use of PDB format coordinate
files for Amber2Gromos, Taskgroup support for Gromos QM/MM with SPLIT option,
BO-MD with EXTRAPOLATE WFN fully restartable, access to QM and MM energy in
QM/MM calculations, and improvements of the manual.

Version 3.14
~~~~~~~~~~~~

Developer’s version.

Version 3.15
~~~~~~~~~~~~

New features, performance improvement and bug fixes have been added since the
latest version of the code. These include a new efficient orthogonalization
scheme [be_cur]_, Constrained Density Functional
Theory [ober]_, force matching in QM/MM runs, the new generalized
Langevin thermostat  [ceriotti]_, the screened hybrid functional
HSE06 from the group of Scuseria  [HSE06a]_ [HSE06b]_, multiple
walkers in metadynamics, surface hopping dynamics with non-adiabatic coupling
vectors in TDDFT [taver]_, extensions of Grimme vdW corrections
and initial support for the Ehrenfest Molecular
Dynamics [taver1]_, and Kerker mixing for handling metallic
systems [Kerker]_. The new version includes also ports to IBM
POWER7, Fujitsu-Primergy BX900, several Linux updates and an updated
pseudopotential library.

Version 3.17
~~~~~~~~~~~~

New features, scalability improvements and bug fixes have been added since the
latest version of the code. Among the several novel implementations we number a
new highly parallel scheme for the evaluation of the Hartree–Fock exact
exchange with an efficient thresholding, support for parallel I/O in writing
and reading the RESTART file, the introduction of a second parallelization
layer over the molecular states which replaces the original TASKGROUP
parallelization, increasing further the scalability together with the already
available parallelization over plane-waves. Finally, improved OpenMP 3.0
support has been introduced both in the main code and in the QM/MM interface.
This allows for a substantial speedup of the code in a hybrid OMP/MPI
distribution (see note below). Between the novel methodologies we mention the
availability of the *ab-initio* vdW corrections, of a fully functional
Ehrenfest Dynamics module and of several improvements to the FO-DFT scheme.
Additional developments include novel schemes for treating QM/MM link-atoms and
the possibility of post-processing Nosé-Hoover trajectories and energies. The
new version includes also porting and optimization to the IBM BlueGene/Q and
different updates to the Linux architecture files.

*Note on OMP 3.0*: The *collapse(n)* clause, although powerful in well nested
loops, may not work in old compilers and is known to have troubles e.g. in old
Intel Fortran Compilers (versions 11 and former ones). Please, check carefully
your OS and your compiler before compiling blindly and refer to the discussion
in http://openmp.org/forum/viewtopic.php?f=3&t=1023&start=10.

CPMD Version 4
--------------

Version 4.0
~~~~~~~~~~~

A new development version started in 2012 (and developed concurrently with
version 3.17.1) with the purpose to refactor CPMD.

Version 4.1
~~~~~~~~~~~

Release 4.1 is the new version made available after almost two years from the
last public release. It sets a turning point to CPMD and its new position as a
modern software (object oriented) while still retaining all the features and
functionalities from the previous versions.

Among the major technical breakthroughs: all memory allocations are converted
to dynamical allocations using standard fortran; extensive usage of modules so
to promote a code re-usage policy; fixed several instabilities due to
arithmetical exceptions; all floats operations are now converted to arbitrary
precision; support of few more compilers; all COMMONs have been converted to
TYPES; extensive usage of IMPLICIT NONE; removed a massive amounts of unused
variables and procedures. Finally, we have build-up a regression tester that
covers more than half of the functionalities (700 tests and keeps increasing
day by day), providing a quality control for developers and users.

Beyond cleanup and bug fixes improving the stability of all functionalities, we
also report the implementation of new features such as: Spin-Orbit-Coupling,
local control, inter-system-crossing SH dynamics together, Adiabatic Bohmian
dynamics; coupling with IPHIGENIE to perform QM/MM calculations with
polarizable schemes and a wider usage of CP_GROUPS to improve scalability of
different other functionalities. The new version includes also porting,
optimisation and different updates for various architectures. The work
initiated with version 4.0 will continue tirelessly on the forthcoming version.

Version 4.3
~~~~~~~~~~~

With the version 4.3, we have ported some of the core-procedures to gpu. The XC
procedures have been refactored and it was added the possibility to link to
libxc. Coulomb-attenuated functionals (CAM) and the Minnesota functionals (M05,
M06, M08 and M11 families) are now available as internal xc procedures. A new
development for Ehrenfest dynamics with ionized states is deployed. Few novel
implementations are available: a multiple time step scheme for MD, accelerated
exact exchange for isolated systems, the stress tensor, LSD and UltraSoft
Vanderbilt pseudopotentials and NPT. Also a new bicanonical ensemble method is
made available with version 4.3. Extensive cleanup and bug fixes and
strengthening of the regression test suite.

Version 4.5
~~~~~~~~~~~

The Minnesota non-separable exchange functional families (MN12, MN15) as well
as the revM06-L functional are now available.
