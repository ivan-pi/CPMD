.. _faq:

Questions and Answers
=====================

The following section is a slightly edited collection of questions and answers
from the cpmd mailing list, cpmd-list@cpmd.org.

How to Report Problems
----------------------

Up front a few remarks on how to report problems (and how to respond), so that
the chances to solve the problem (permanently) are as high as possible.

If you have compilation problems, please always state what version of CPMD you
are trying to compile and what kind of machine you are using, i.e. what
operating system, what compiler (particularly important on linux machines),
which compilation flags, and what libraries you are using. Best you include the
first part of your makefile (up to ‘End of Personal Configuration’, please
**don’t** post the whole makefile) as this contains most of the required
information. Also include the **relevant** part of the make output (again, the
full output usually is very long and rarely needed).

If you have problems with a specific calculation, please include your input and
the output of the run, so that others can try to reproduce the error. Again,
please state the version of CPMD you are using and the platform you are running
on.

A useful article about how to ask questions the smart way is at:

http://www.catb.org/~esr/faqs/smart-questions.html

.. _sec:expl-warn-error:

Explanation of Warnings and Error Messages
------------------------------------------

| If execution of CPMD aborts with an error message, please do not only consult
  the standard output, but also the **LocalError** files which are produced in
  the run directory. They may contain more valuable information of the problem
  at hand.
|

FAQ on warnings and errors:

-  :ref:`The code crashes without leaving any error message. <faqcrash>`

-  :ref:`WANNIER CODE WARNING: GRADIENT FOR RESTA FUNCTIONAL <faqwannier>`
   (during CP molecular dynamics)

-  :ref:`WARNING! XC FUNCTIONALS INCONSISTENT <faqxcpp>`

-  :ref:`Warning! Spline region smaller than maximum G-value <faqspline>`
   (during geometry optimization with variable cell)

-  :ref:`GORDER\| PROGRAMING ERROR. INFORM THE PROGRAMER <faqgorder>` (during
   startup of WAVEFUNCTION or GEOMETRY optimization)

.. _faqcrash:

**Q:** My CPMD run suddenly aborted without any error message. What happened?

**A:** CPMD will only print error messages to the output file in some very
specific cases. Usually, if an error is encountered, every task writes to a
file called *LocalError-X-X-X.log* (where :math:`x` denotes a task ID), where
you will find an *error message*, the *procedure* and *file name* where the
error occured, as well as a *call stack* (although the latter may sometimes be
incomplete). In case the error message is empty, you may have to consult the
source code directly.

If no *LocalError* files are produced, the processes may have exited due to an
insufficient amount of RAM available for a single task.

.. _faqwannier:

**Q:** Could anybody explain to me the follwing error in the cpmd output during
CP Molecular Dynamics runs with the flag :kw:`WANNIER WFNOUT`
LIST DENSITY?

::

   WANNIER CODE WARNING: GRADIENT FOR RESTA FUNCTIONAL IS GMAX=0.118E-02

Does it mean any serious error in the calculation?

**A:** The default spreadfunctional used in CPMD is the Vanderbilt type. At the
end of the calculation the convergence with respect to the Resta type
functional is also checked. For large cells both should be converged at the
same time. However, for typical application this is not the case and you get
the warning. This is not serious and you can ignore it.

.. _faqxcpp:

**Q:** A warning message appeared in the output file:

::

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !WARNING! XC FUNCTIONALS INCONSISTENT FOR h.pp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Does it mean that my pseudopotential file is wrong? I used fhi98PP to create
this pp file, and just recast it. So if it is wrong, what more should I do to
the output file of fhi98PP?

**A:** It means that the XC functional used to generate the pseudo potential
and the functional which you want to use in CPMD are not the same. This could
be more or less serious depending on the specific situation: Some of the XC
functionals in CPMD are just minor variations of each other, e.g. the LDA part
may evaluated using Perdew-Wang, Perdew-Zunger, VWN or the Pade-interpolation
formula. If this is the case, the resulting error is usually small. However, if
your pseudo potential was generated *e.g.* with PBE and you try to use it with
the BLYP functional, the error might become non-negligible. You may compare the
two functionals (pp & CPMD-input), and judge for yourself if the difference is
significant or not.

For **hybrid** and **meta** functionals, there are usually no pseudo potentials
available. In these cases, it is common practice to use the pseudo potentials
of related GGA xc functionals (*e.g.* PBE for TPSS, BLYP for CAM-B3LYP, etc.),
or, if there is no ‘close relative’, LDA (this is possible *e.g.* in
combination with the Minnesota xc functionals).

.. _faqspline:

**Q:** I am trying to optimize a crystal structure (both ion positions and cell
volume) using CPMD and get the warning message:

Warning! Spline region smaller than maximum G-value.

The optimization seems to converge nicely but what does this warning
mean/imply?

**A:** This in fact touches several points. The following applies also to other
variable cell calculations with CPMD (e.g. constant pressure simulations).

-  Pseudopotential functions in CPMD are calculated on a radial grid in G-space
   and then used in spline interpolations. This speeds up variable cell
   calculations considerably. The maximum grid point is given by the cutoff.
   With the keyword

   ::

         SPLINE RANGE
           x.xx

   you can enlarge the grid to x.xx times the cutoff. Also, you should make
   sure, that the number of spline points is large enough. Older version of
   CPMD defaulted to as little as 501. This is at the lower limit for accuracy
   with a fixed cell. Especially if you have high cutoffs it is better to
   increase this value, e.g.

   ::

         SPLINE POINTS
           2500

   or larger. The current default (5000) should be large enough.

-  In a variable cell calculation, CPMD uses a constant number of plane waves.
   Therefore if your cell contracts the cutoff increases, if the cell gets
   larger the cutoff decreases. So if you have the first case the spline
   interpolation needs points above the original cutoff and you get a warning.
   Depending on the amount of change in the cell you expect a value for the
   :kw:`SPLINE` **RANGE** of 2–5 is needed.

-  Coming back to the constant number of plane waves. If your cell gets larger
   the effective cutoff decreases. This may have very undesirable effects and
   it is better to define the plane waves on a box larger than any box you
   anticipate the simulation will reach. In order to not have to start with an
   unreasonable cell you can define the plane waves not with the actual box but
   with a reference cell, use the keyword

   ::

         REFERENCE CELL
           a b c xa xb xc

-  However, you really want to do a constant cutoff calculation, not a constant
   number of plane waves. For technical reasons this is not possible and in
   principle you should do the calculation at a high enough cutoff in order
   that the calculation is converged all along the simulation path (with the
   slight changes in cutoff).

   To avoid these very high cutoffs the group in Trieste came up with a method
   that allows to perform pseudo constant cutoff calculations. This method is
   implemented in CPMD (keyword :kw:`CONSTANT CUTOFF`) and
   explained in the paper [bernasconi95]_

.. _faqgorder:

**Q:** I was optimizing wavefunctions for my system. After a successful run I
modified CELL VECTORS. This latter calculation crashed with the following error
either when I restarted (only wavefunctions) or if I started from scratch.

::

    GORDER| PROGRAMING ERROR. INFORM THE PROGRAMER

Is there something wrong with the code?

**A:** This is a known problem in CPMD. The message comes from a test, which
probes whether the G vectors are in a “safe” order, namely such that after a
restart with a different number of processes or on a different machine the
results agree. Usually this error only occurs in large systems with a high
cut-off energy and/or large unit cells, i. e. where one gets lots of
close-lying G vectors.

There are two possible workarounds:

#. Slightly change your computational box in one dimension e. g. from
   :math:`10.000000` to :math:`10.000001`. This helps some times.

#. The check is not 100 % accurate. This means by just ignoring the message you
   will most likely get correct results. The error would only appear in
   restarts where you could see a small inconsistency in energy in the first
   step. Final results should not be affected (except for MD if you do
   restarts).

   To avoid the stop, comment out the two lines at the end of file
   loadpa_utils.mod.F90 of the form

   ::

              CALL STOPGM('GORDER','ERROR IN G-VEC ORDERING (NHG)')

   However, be sure to check that the results are reasonable.

.. _sec:FAQPP:

.. _pseudopotentials-1:

Pseudopotentials
----------------

FAQ on pseudopotentials:

-  :ref:`Which type of pseudopotentials to use? <faqwhichPP>`

-  :ref:`Which is the correct value of LMAX? <faqwhichLMAX>`

.. _faqwhichPP:

**Q:** I’m confused about how to select a pseudopotential type
(Troullier-Martins, Goedecker, etc.). What makes one choose say a Goedecker
potential instead of a Vanderbilt potential?

**A:** The choice of a pseudopotential in CPMD calculations depends on needs,
available resources and taste.

Troullier-Martins norm-conserving pseudopotentials are probably the
most-commonly used type of pseudopotentials in CPMD calculations. They work
very well for non-problematic elements and they are quite easy to create (note,
that it is also easy to create a bad pseudopotential). When using the
Kleinman-Bylander separation, one also has to be careful to avoid so called
ghost states (e.g. many transition metals need LOC=\ :math:`l` with :math:`l`
being an angular momentum smaller than default value which is highest).

Goedecker pseudopotentials are stored in an analytical form, that ensures the
separability, but they usually need a higher (sometimes much higher) plane wave
cutoff than the corresponding Troullier-Martins counterparts. Also the creation
procedure is more complicated, but there is a very large library of tested
pseudopotentials (mostly LDA but also some GGA pseudopotentials).

Vanderbilt pseudopotentials have the advantage of needing a much reduced plane
wave cutoff. The drawback is, that only a limited subset of the functionality
in CPMD is actually implemented with uspps (MD, wavefunction/geometry
optimization and related stuff and only at the gamma point and you have to make
sure, that your real space grid is tight enough). Also due to sacrificing
norm-conservation for softer pseudopotentials, your wavefunction has very
limited meaning, so that not all features available for norm-conserving
pseudopotentials can actually be easily implemented or implemented at all.

For some elements it can be rather difficult to generate good (i.e.
transferable) pseudopotentials, so you should always check out the available
literature.

.. _faqwhichLMAX:

**Q:** How do I choose the correct value of LMAX?

**A:** If you use a Vanderbilt or Goedecker type potential only the format of
the LMAX-line has to be valid. The actual value is read from the
pseudopotential file and the value in the input file will be ignored. It is
highly recommended to still use values that make sense, in case you want to do
a quick test with a numerical (Troullier-Martins) pseudopotential.

Generally, the highest possible value of LMAX depends on the highest angular
momentum for which a “channel” was created in the pseudopotential. In the
pseudopotential file you can see this from the number of coloumns in the
&POTENTIAL section. The first is the radius, the next ones are the orbital
angular momenta (s, p, d, f,…). As an example you can determine a potential for
carbon using f-electrons and set LMAX=F. Since the f-state is not occupied in
this case there is very little advantage but it costs calculation time. In
short, you can use values as high as there is data in the pseudopotential file,
but you don’t have to if it is not needed by the physics of the problem.

A fact that causes confusion is that Hamann’s code for pseudopotential
generation always produces output for the d-channel, even if you only request
channels s and p. You should be cautious if :math:`r_{c}` and the energy
eigenvalue of the p- and d-channels are equal. In most of these cases LMAX=P
should be used.

.. _sec:file-formats:

File Formats and Interpretation of Data
---------------------------------------

FAQ on interpretation of data:

-  :ref:`Why is my total energy so much different from a Gaussian calculation? <faqetot>`

-  What is the meaning of the :ref:`energies <faqenergies>` reported in a
   molecular dynamics simulation?

-  What do :ref:`GNMAX, GNORM, and CNSTR <faqgnmax>` mean in a geometry
   optimization?

-  All :ref:`IR intensities are zero <faqir-int>` in VIB.LOG when I try to
   calculate the IR of NH\ :math:`_4^+` ion by CPMD

-  Why can some atoms move far away from the central box in a long MD
   simulation of :ref:`bulk liquid using periodic boundary conditions <faqliquidpbcs>`?

-  Why did the :ref:`electron energy continuously increase <faqeneele>` during
   the MD simulation of my metal?

-  What is the meaning of the files :ref:`APT, POLARIZATION and POLARIZABILITY <faqraman>` resulting from my linear response
   calculation of Raman spectra?

-  What is the meaning of :ref:`columns 2 to 7 in the DIPOLE file <faqdipole>`?

-  How do I use the :ref:`binary RESTART file across different platforms <faqrestart>` with different encoding?

.. _faqetot:

**Q:** Why is my total energy so much different from a Gaussian calculation?

**A:** With CPMD you are using **pseudopotentials** to describe the atoms.
Since the total energy describes only the interactions between the pseudocores
and the valence electrons (and **some** core electrons in the case of so-called
semi-core pseudopotentials), you are missing the contribution of the core
electrons and the full core charges of a regular all-electron calculation.
**Energy differences** between two configurations, on the other hand, should be
comparable, provided you use the same number of atoms, the same plane wave
cutoff, the same pseudopotentials, and the same supercell geometry in the CPMD
calculation.

.. _faqenergies:

**Q:** In a molecular dynamics simulation, CPMD prints out a list of energies
for each integration step. Does anyone know the meaning of the individual
values.

**A:** Some explanations to the energy terms:

EKINC
   fictitious kinetic energy of the electrons in a.u. this quantity should
   oscillate but not increase during a simulation.

TEMPP
   Temperature of the ions, calculated from the kinetic energy of the ions
   (EKIONS).

EKS
   Kohn-Sham energy (the equivalent of the potential energy in classical MD).

ECLASSIC
   = EKS + EKIONS

EHAM
   = ECLASSIC + EKINC. Hamiltonian energy, this is the conserved quantity,
   depending on the time step and the electron mass, this might oscillate but
   should not drift.

DIS
   mean square displacement of the ions with respect to the initial positions.
   Gives some information on the diffusion.

You can modify the list of individual energies to be displayed with the :kw:`PRINT ENERGY` keyword.

.. _faqgnmax:

**Q:** What do GNMAX, GNORM and CNSTR mean in a geometry optimization?

**A:** These are abbreviations for the following quantities:

GNMAX
   max\ :math:`{}_{I,a}` (:math:`|F_{Ia}|`) = largest absolute component
   (:math:`a=x,y,z`) of the force on any atom :math:`I`.

GNORM
   :math:`\left< F_{I}^2\right>_{I}` = average force on the atoms :math:`I`

CNSTR
   max\ :math:`{}_{I,a} { F^{constr}_{Ia} }` = largest absolute component
   (:math:`a=x,y,z`) of force due to constraints on any atom :math:`I`.

.. _faqir-int:

**Q:** I found all the IR intensities in VIB.log file are zero when I try to
calculate the IR of NH\ :math:`{}_4^+` ion by CPMD.

::

    Harmonic frequencies (cm**-1), IR intensities (KM/Mole),
    Raman scattering activities (A**4/AMU), Raman depolarization ratios,
    reduced masses (AMU), force constants (mDyne/A) and normal coordinates:
                        1                      2                      3
                       ?A                     ?A                     ?A
    Frequencies --   142.9800               188.9340               237.2614
    Red. masses --     0.0000                 0.0000                 0.0000
    Frc consts  --     0.0000                 0.0000                 0.0000
    IR Inten    --     0.0000                 0.0000                 0.0000
    Raman Activ --     0.0000                 0.0000                 0.0000
    Depolar     --     0.0000                 0.0000                 0.0000
    Atom AN      X      Y      Z        X      Y      Z        X      Y      Z
      1   7     0.00   0.00   0.00     0.00   0.00   0.00     0.00   0.00   0.00
      2   1     0.00  -0.35  -0.50    -0.35   0.00   0.00    -0.50   0.00   0.00
      3   1     0.00  -0.35   0.50    -0.35   0.00   0.00     0.50   0.00   0.00
      4   1     0.00   0.35   0.00     0.35   0.00   0.50     0.00  -0.50   0.00
      5   1     0.00   0.35   0.00     0.35   0.00  -0.50     0.00   0.50   0.00
                        4                      5                      6

**A:** That’s not a problem of your calculation. The keyword **VIBRATIONAL
ANALYSIS** does not calculate intensities. The calculation of intensities is
currently not possible in CPMD. The intensities in the ’VIBx.log’ files are
arbitrarily set to zero. The entries have to be there so that visualisation
programs, that are able to read output of the Gaussian program, can be also
used to visualize the CPMD results.

.. _faqliquidpbcs:

**Q:** I am trying to simulate a bulk liquid in CPMD and supposing that
periodic boundary conditions are built into the program. But after several
thousand MD steps, I found some particles are far away from the central
simulation box.

Why it is so if periodic boundary conditions (PBC) on particle coordinates are
imposed in all three directions?

**A:** If you are not using the

::

   SYMMETRY
     0

options your calculations are actually using periodic boundary conditions
(PBC). PBC are imposed within CPMD for all calculations. However, the particle
positions are not folded back to the original computational box. The reason for
this is that most people prefer to have “smooth” trajectories without jumps of
particles. This allows for easier tracking of special particles and nicer
graphics. In addition it is easy (with a little script) to apply PBC afterwards
yourself, if needed.

.. _faqeneele:

**Q:** I am trying to simulate a bulk sodium and I found electron energy is
increasing continuously and it is in the range of 0.07 a.u. at the end of 20000
steps.

**A:** Sodium is a metal, and therefore missing an important feature that
allows for stable CP dynamics: the band gap. Using Nosé thermostats (on
electrons and ions) it might still be possible to perform meaningful CP
simulations [Blochl92]_.

The choice of parameters for the thermostats, however, will be nontrivial,
highly system dependent and require extensive testing. Without thermostats you
will have strong coupling between electronic degrees of freedom and ionic
degrees of freedom. Adiabaticity is not maintained and a steady increase of the
fictitious kinetic energy will occur.

.. _faqraman:

**Q:** I have computed RAMAN by LINEAR RESPONSE, and get three files: APT,
POLARIZATION and POLARIZABILITY with lots of data in these files. I want to
know the meaning of the data, please give me some answer in detail.

**A:** The POLARIZABILITY file simply contains the polarizability tensor of the
whole system in atomic units. The POLARIZATION file contains the total dipole
moment (electronic + ionic) of the whole system in atomic units. As for the
file APT, it contains the atomic polar tensors for each atom in the system. The
atomic polar tensor is the derivative of the forces on the atoms with respect
to an applied external electric field. Equivalently it is, from a Maxwell
relation, the derivative of the total dipole of the system with respect to the
nuclei positions. It is thus an important ingredient of the calculation of
infrared spectra intensities, for example used in an harmonic approximation.
The trace of this tensor is the so-called Born charge of the considered atom.
The data is arranged in the following order (still in a.u.): the APT tensor is
:math:`\frac{\mathrm{d}F_{I,i}}{\mathrm{d}E_{j}}` where :math:`F_{I,i}` is the
force on atom I along :math:`i=x,y,z` and :math:`E_j` is the electric field
along :math:`j=x,y,z`. :math:`(I,i)` are the indices of the :math:`3N` atoms
lines in the APT file, one atom after the other, and :math:`j` is the column
index in the APT file.

.. _faqdipole:

**Q:** I was wondering what columns 2 to 7 in the DIPOLE file correspond to?
When I run CPMD v.3.5.1, columns 2 to 4 come out identical to columns 5 to 7
respectively. When I run with CPMD v.3.4.1, the columns come out different. Is
there an explanation for this?

**A:** Columns 2 to 4 in the DIPOLE file are the electronic contribution to the
dipole moment, columns 5 to 7 are the total (electronic + ionic) dipole moment.
All dipole moments are divided by the volume of the box.

In CPMD version 3.5.1 we have changed the reference point of the calculation.
Now the reference point is chosen such that the ionic contribution is zero and
the electronic contribution minimal (=total dipole). This avoids a problem that
occasionally was seen in older versions. The electronic dipole is calculated
modulo(2\ :math:`\pi`/L). Now if the electronic dipole became too large,
because the ionic contribution was large (bad choice of reference point) the
total dipole made jumps of 2\ :math:`\pi`.

.. _faqrestart:

**Q:** As you know, the cpmd RESTART file is saved as binary. But I want to
change it to ASCII and vice versa, because I use several machines of different
architecture, for example COMPAQ, IBM, and LINUX machine. Please help me with
any comments.

**A:** The code to read and write the RESTART file is in the files
rv30_utils.mod.F90 and wv30_utils.mod.F90. Feel free to implement an ASCII
version of the restart, but be aware that the file will be **huge**.

But you may not need to do that. Let’s say you decide to use big-endian binary
encoding (this is what e.g. IBM, Sun and SGI machines do natively).

With Compaq machines there is a compiler flag, -convert, which you could set to
big_endian (we only have here linuxalpha, but the compaq compiler should be
essentially the same).

On a Linux PC you can use the use the -Mbyteswapio or the -byteswapio flag, if
you have the PGI compiler.

| For the Intel compiler (ifc/ifort/efc) you simply set the environment
  variable F_UFMTENDIAN to big (i.e.
| ’export F_UFMTENDIAN=big’
| if you are in a bourne/korn shell and
| ’setenv F_UFMTENDIAN big’
| if you are in a (t)csh).

Now even your cpmd executables will read and write big-endian restart files.

Check your compiler documentation for more details (search for endian).

.. _sec:input-param-valu:

Input Parameter Values
----------------------

FAQ on input parameters:

-  If I set the keyword :ref:`RESTART WAVEFUNCTION COORDINATES <faqsystematomrestart>`, would I have to write the &SYSTEM
   and &ATOM section again?

-  Could anybody tell me :ref:`how to choose the energy cutoff in &SYSTEM section <faqcutoff>`?

-  I have a problem with :ref:`visualising unoccupied orbitals <faqunocc>`.
   Either I get only occupied orbitals, or, if I add one empty state when
   optimizing wavefunction the program never reaches convergence.

-  Is there any way to force :ref:`CPMD to dump DENSITY files every N steps of a molecular dynamics run <faqdumpdensity>`?

-  :ref:`How do I calculate a Band structure with CPMD <faqbandstruct>`?

-  I want to collide fast atoms against surfaces. Why are the :ref:`initial velocities given using the VELOCITIES keyword <faqvelocities>` not
   considered in the calculation?

-  I want to to run CPMD with :ref:`basis sets equivalent to Gaussian <faqgaussianbasis>` 6-31+G(d) and 6-311+G(2d,p). How do i set
   up the &BASIS section?

-  :ref:`How do I add support for a new functional? <faqnewdft>`

.. _faqsystematomrestart:

**Q:** If I set the keyword RESTART WAVEFUNCTION COORDINATES, would I have to
write the &SYSTEM and &ATOM section again?

**A:** Yes, you have to include the &SYSTEM and &ATOM sections even if you are
restarting. If you write RESTART COORDINATES, the coordinates in the RESTART
file override the ones in the input. RESTART WAVEFUNCTION alone does not select
the coordinates in the RESTART file, but does use those in the &ATOMS section.

.. _faqcutoff:

**Q:** Could anybody tell me how to choose the energy cutoff in &SYSTEM
section?

**A:** The best way to choose the cutoff for CPMD calculations is by running
first a series of tests. Select a test system and a representative quantity
(bond length, reaction energy, etc.), perform a series of calculations with
increasing cutoff, pick the lowest cutoff with satisfactory results.

It’s always a good idea to make checks at some critical points of the
calculations by increasing the cutoff. See also section
:ref:`hints:cutoff`.

.. _faqunocc:

**Q:** I have a problem with visualising unoccupied orbitals. When I use
:kw:`RHOOUT` BANDS or :kw:`CUBEFILE` ORBITALS after the
wavefunction optimization I get only occupied orbitals. If I add one empty
state when optimizing wavefunction the program never reaches convergence.

**A:** The most efficient way to calculate unoccupied orbitals is to first
optimize the occupied orbitals and then restart the calculation using the run
option

::

    KOHN-SHAM ENERGIES
     n

where n ist the number of unoccupied orbitals. This will diagonalize the
Kohn-Sham Potential (defined by the occupied orbitals alone).

To test if everything goes fine, you can check the total energy printed at the
beginning of this job, it should be exactly the one at the end of the
optimization. In addition, if you don’t change the default convergence
criteria, the number of converged Kohn-Sham states should be equal to the
number of occupied states in the first step.

.. _faqdumpdensity:

**Q:** Is there any way to force CPMD to dump DENSITY files every N steps of
molecular dynamics run instead (or except) of the end of the job?

**A:** Short of modifying the source code, you could set the parameter
:kw:`RESTFILE` to a large number and than have CPMD write a restart
file every N steps via the :kw:`STORE` keyword. Now you rename each
restart in turn from RESTART.# to RESTART and do a single step calculation
using the :kw:`RESTART` keyword without the LATEST modifier which
will write the DENSITY file (or run a :kw:`PROPERTIES` job using
:kw:`CUBEFILE` DENSITY to get the cube file directly).

.. _faqbandstruct:

**Q:** How do I calculate a Band structure with CPMD? To calculate a band
structure with CPMD, You first calculate the correct density for your system
with a Monkhorst-Pack Mesh.

**A:** Then you use: OPTIMIZE WAVEFUNCTIONS with MAXSTEP 1 (no
self-consistency) and RESTART DENSITY.

In the section KPOINTS, you should use for instance a bcc:

::

      KPOINTS BANDS
           51   0   0   0      0   0   1             Gamma to H
           51   0   0   1      0  .5  .5             H to N
           51   0  .5  .5     .5  .5  .5             N to P
           51  .5  .5  .5      0   0   0             P to Gamma
           51   0   0   0     .5  .5   0             Gamma to N
           51   0   0   1     .5  .5  .5             H to P
           0    0 0 0  0 0 0

You say that you want 51 points from (0,0,0) and (0,0,1) and so on. The last
line with many zeros is to stop.

If the memory of your computer is not enough, you can add in the line KPOINTS
the option BLOCK=50 that means you want to have only 50 kpoints in memory. This
options worked some time ago.

.. _faqvelocities:

**Q:** I’ve been recently trying to use the VELOCITIES keyword in a molecular
dynamics run. I want to collide fast atoms against surfaces. Despite the code
seems to read the input velocities properly, when the run starts the initial
velocities are always the same (apparently coming from a thermal distribution),
no matter what is the velocity you specify for the incoming atom. I’m not using
QUENCH IONS, so I don’t understand why the input initial velocities are not
considered in the calculation.

**A:** There is no straightforward way in CPMD to achieve what you want. I
suggest to follow this procedure:

1) Run a single step of MD with the following set up

::

   MOLECULAR DYNAMICS
   MAXSTEP
    1
   RESTART WAVEFUNCTION COORDINATES
   TEMPERATURE
     300   <- or whatever your surface should be

This generates RESTART and GEOMETRY files. Now edit the GEOMETRY file to change
the velocities of the particles according to your experiment. Now restart the
MD with the options

::

   MOLECULAR DYNAMICS
   MAXSTEP
    1000
   RESTART WAVEFUNCTION COORDINATES VELOCITIES GEOFILE
   QUENCH ELECTRONS

The effect of this is: IONIC coordinates and velocities are read from GEOMETRY,
ELECTRON wavefunctions and velocities are read from RESTART, ELECTRON
velocities are set to zero.

.. _faqgaussianbasis:

**Q:** I want to to run CPMD with basis sets equivalent to Gaussian 6-31+G(d)
and 6-311+G(2d,p). How do I set up the &BASIS section?

**A:** You should be able to construct inputs from the description in this
manual (see section :ref:`input:basis`).

Please note, that the basis set generated from the &BASIS section is used in
CPMD for two purposes:

#. Analyzing orbitals:

   We usually use the atomic pseudo-wavefunctions to analyze the orbitals from
   CPMD. The 6-31G type Gaussian basis sets are for all electron calculations.
   Don’t expect very good results when analyzing wavefunctions from a
   pseudopotential calculation.

#. Generating orbitals for an initial guess:

   By default we use a Slater minimal basis. In most cases the effort to
   produce a better initial guess using “better” wavefunctions does not pay
   off.

.. _faqnewdft:

**Q:** How do I add support for a new functional?

**A:** The exchange-correlation functionals and scratch variables are set and
handled in the module cp_xc_utils.mod.F90. The individual functionals are
grouped in modules (cp_gga_exchange_utils.mod.F90 etc.) which are in turn used
in cp_xc_utils.mod.F90. If you wish to implement higher-order analytical
derivatives for LR-TDDFT, cp_dxc_utils.mod.F90, and the cp_dgga…routines may be
of interest.

The inclusion of hybrids will require additional modifications to
dftin_utils.mod.F90 and cpfunc_func_init (in cp_xc_utils.mod.F90). The
variables for range separation are currently still found in func.mod.F90 and
the range separated exact exchange is handled in initclust_utils.mod.F90.

.. _faqminnesotaconv:

**Q:** I am trying to use the Minnesota functionals (M06, M11, …), but the
results I get are most erratic, or they require very high values for the
:kw:`CUTOFF`. What is wrong?

**A:** The Minnesota functionals require unusually dense integration grids.
Simply increasing the cutoff will introduce numerical noise and make the
wavefunction optimisation unstable. This issue can be overcome by increasing
the :kw:`DUAL` up to 8-12. A detailled assessment of these effects
along with recommendations for appropriate values can be found in the following
reference: M.P. Bircher, P. Lopez-Tarifa and U. Rothlisberger: *J. Chem. Theory
Comput.*, **15** (1), 557 (2019). DOI: 10.1021/acs.jctc.8b00897

.. _faqminnesotacomp:

**Q:** Results obtained from the Minnesota functionals (M05, M06, M08, M11)
seem well converged, but exhibit significant differences with respect to
Gaussian bases. Is there something wrong?

**A:** Due to their highly flexible functional form, the Minnesota functionals
are extremely sensitive to the choice and flexibility of the underlying basis
set. Since Gaussian functions always impose a certain rigidity, results
obtained in Gaussian basis sets may occasionally differ by more than just
chemical accuracy from the plane wave results. This effect is not limited to
plane waves, it can also be observed - to a lesser extent - when comparing
plane waves and Slater functions, and can be rather pronounced when comparing
Gaussian and Slater functions. A detailled study of basis set effects is
provided in the following article: M.P. Bircher, P. Lopez-Tarifa and U.
Rothlisberger: *J. Chem. Theory Comput.*, **15** (1), 557 (2019). DOI:
10.1021/acs.jctc.8b00897

.. _faqhfxslow:

**Q:** I am running a calculation with exact exchange, but it is exceedingly
slow compared to a GGA or MGGA run. Can this issue be resolved?

**A:** The implementation of exact (or Hartree-Fock) exchange in CPMD is highly
efficient and parallel. If the number of available nodes is limited, increasing
the value in :kw:`HFX_BLOCK_SIZE <HFX BLOCK SIZE>` can substantially speed up
the calculation. In the presence of ample computational resources, the
calculation can be sped up by using :kw:`CP_GROUPS <CP GROUPS>`.

In isolated systems, the calculation can be sped up by up to an additional
order of magnitude by using the coordinate-scaling scheme by Bircher and
Rothlisberger. For details, see :kw:`SCEX`.
