.. _sec:qmmm:

CPMD/Gromos QM/MM Calculations
==============================

.. _sec:qmmm-overview:

General Overview
----------------

An additional interface code (**MM_Interface** folder) and an adapted classical
force field code [gromos96]_ (**Gromos** folder) are needed to run
CPMD in fully Hamiltonian hybrid QM/MM mode\ [qmmm02]_. To use
this code a *Gromos license* is required and therefore it is **not** included
in the standard CPMD code. The interface code and the adapted classical force
field code can be obtained by directly contacting the CPMD developers.

To create a makefile for compilation of a QM/MM enabled CPMD binary, you have
to copy the two above folders (or create symbolic links to them) in the CPMD
source directory and then to add the ``-qmmm`` flag when executing the
``config.sh`` script (see section :ref:`installation`). The resulting
binary can be used for normal CPMD runs as well as for QM/MM simulations.

.. _sec:qmmm-input:

Input files for QM/MM CPMD
--------------------------

A QM/MM run requires a modified CPMD input file, some additional input files,
and creates the normal CPMD output file and some new ones. The input file
consists of a standard CPMD input with the :kw:`QMMM` keyword in
the &CPMD section, a modified &ATOMS section and a mandatory &QMMM section.
Furthermore three files for the classical code are needed (coordinates,
topology and input file). These can be taken from previous fully classical
simulations and have to be in Gromos format. Topologies and coordinates files
created with the Amber\ [amber7]_ package are also supported. A
converter to Gromos format\ [gromos96]_ is available.

.. _sec:qmmm-start:

Starting a QM/MM run
--------------------

To start a QM/MM simulation, you first do a simulation of your system with a
regular classical MD-code to get an equilibrated configuration. The tricky part
in this is usually the treatment of (metal-)ion or special molecules, that are
not parameterized within a given force field but are in the active center of
your molecule (one of the prominent reasons why you want to do a QM/MM run in
the first place). It is usually easiest to keep that part rigid throughout the
equilibration, until after you have defined the QM-subsystem.

Starting from the classically equilibrated structure, you have to create a
topology, a coordinate and an input file in Gromos format (either by using the
Gromos tools or a converter). Now you need to define your QM system by
assigning pseudopotentials to selected atoms in your CPMD input file (see
:ref:`sec:defining-qm-system`).

You can now start to continue the classical equilibration with CPMD using
:kw:`MOLECULAR DYNAMICS` CLASSICAL. Please note, that
there are several special constraints available to ease the transition in case
of strong interactions within the QM part or between the QM and the MM part.
Finally, a wavefunction optimization (either directly or via
:kw:`QUENCH` BO) and a normal :kw:`MOLECULAR DYNAMICS` CP or BO can be performed.

.. _sec:defining-arraysizes:

Defining internal Gromos array dimensions
-----------------------------------------

One rather new feature of this QM/MM interface code is the :kw:`ARRAYSIZES ... END ARRAYSIZES` block in the &QMMM section
which allows to change the internal array dimensions of the Gromos part
dynamically. Previously one had to change some include files and recompile
everything to adapt the code for a different system.

These settings have to be consistent during a series of calculations, or else
you may not be able to read your restart files correctly.

.. _sec:defining-qm-system:

Defining the QM system
----------------------

For a QM/MM calculation a subset of atoms are selected from the classical
restart and then for this QM part an isolated system (SYMMETRY 0) calculation
is performed. The supercell size has to follow the requirements of the various
Poisson solvers, as listed in the hints section (:ref:`hints:symm0`).

If not otherwise specified, the QM system (atoms and wavefunction) is always
re-centered in the given supercell (the current offset of the QM cell is
recorded in the file MM_CELL_TRANS).

The quantum atoms are specified in the &ATOMS section similar to normal CPMD
calculations. Instead of explicit coordinates one has to provide the atom index
as given in the Gromos topology and coordinates files.

The keywords of the ``&QMMM`` input section are described in
:doc:`/input/qmmm`.

.. _sec:qmmm-gromos-inp:

Keywords in the Gromos Input and Topology files
-----------------------------------------------

For a detailed description of the Gromos file formats please have a look at the
Gromos documentation\ [gromos96]_. Note, that not all keyword are
actually active in QM/MM simulations, but the files still have to be
syntactically correct. Both, the input and the topology file are structured in
sections starting with a keyword in the first column and ending with the
keyword END. Lines starting with a pound sign ’#’ may contain comments and are
ignored. Both files are required to have a TITLE section as the first section.
The rest can be in almost any order. Here is a short list of some important
flags and their meaning.

.. container:: description

   TITLE

      Text that identifies this input file. Will be copied into the CPMD
      output.

   SYSTEM

      This section contains two integer numbers. The first is the number of
      (identical) solute molecules (NPM) and the second the number of
      (identical) solvent molecules (NSM).

   BOUNDARY

      This section defines the classical simulation cell. It contains 6
      numbers. The first (NTB) defines the type of boundary conditions ( NTB
      :math:`<0` means truncated octahedron boundary conditions, NTB=0 vacuum,
      and NTB :math:`>0` rectangular boundary conditions).

      The next three numbers (BOX(1..3)) define the size of the classical cell.
      The fifth number (BETA) is the angle between the x- and z-axes and the
      last number usually determines whether the cell dimensions are taken from
      the input file (NRDBOX=0) or from the BOX section of the
      :kw:`COORDINATES` file (NRDBOX=1), but is ignored for
      QM/MM simulations.

      Note: that even for vacuum simulations valid simulation cell sizes must
      be provided.

   PRINT

      This section determines how often some properties are monitored. Here
      only the first number (NTPR) matters, as it determines the number of MD
      steps between printing the various energies to the CPMD output.

      Note many old Gromos input files created by the amber2gromos program
      default to NTPR=1, which makes the CPMD output huge.

   SUBMOLECULES

      Defines number of submolecules in the solute. The first number is the
      number of submolecules followed by the index number of the last atom of
      each submolecule. The last number must be identical to the number of
      atoms in the solute.

   FORCE

      Contains two groups of numbers, that controls the various force component
      and the partitioning of the resulting energies. The first group of 1/0
      flags turn the various force components on or off. The second group
      defines energy groups (the first number is the number of groups followed
      by the index number of the last atom in each group). The last number must
      be identical to the number of all atoms.

   TITLE

      Text that identifies this topology file. Will be copied into the CPMD
      output.

   ATOMTYPENAME

      This section contains the number of classical atom types (NRATT) followed
      by the respective labels, one per line. Note that the :kw:`ARRAYSIZES ... END ARRAYSIZES` MAXATT must be large
      enough to accommodate all defined atom types.

   RESNAME

      This section contains the number of residues in the solute (NRAA2)
      followed by the respective residue names.

   SOLUTEATOM

      This section defines the number (NRP) and sequence of atoms in the
      solute, their names, residue numbers, non-bonded interaction codes,
      masses, charges, charge groups and their full + scaled 1-4 exclusions.

   BONDTYPE

      This section contains the list of parameters for bonded interactions. You
      have to pick the two matching entries from this list for the O-H and H-H
      potential, when using the :kw:`FLEXIBLE WATER` keyword
      to convert solvent water back into solute (e.g. to included them into the
      QM part).

   SOLVENTATOM

      This section defines the number of atoms (NRAM) in the solvent and their
      respective names, non-bonded interactions types, masses, and charges.

   SOLVENCONTSTR

      This section defines the number (NCONS) and parameters for the distance
      constraints, that are used to keep the solvent rigid.

.. _sec:qmmm-files:

Files generated by the interface code
-------------------------------------

-  **QMMM_ORDER**

   The first line specifies the total number of atoms (NAT) and the number of
   quantum atoms (NATQ). The subsequent NAT lines contain, for every atom, the
   gromos atom number, the internal CPMD atom number, the CP species number isp
   and the number in the list of atoms for this species NA(isp). The quantum
   atoms are specified in the first NATQ lines.

-  **CRD_INI.grm**

   Contains the positions of all atoms in the first frame of the simulation in
   Gromos extended format.

-  **CRD_FIN.grm**

   Contains the positions of all atoms in the last frame of the simulation in
   Gromos extended format.

-  **INTERACTING.pdb**

   Contains (in pdb format) all the QM atoms and all the MM atoms in the
   electrostatic coupling NN list. The 5-th column in this file specifies the
   gromos atom number as defined in the topology file and in the coordinates
   file. The 10-th column specifies the CPMD atom number as in the TRAJECTORY
   file. The quantum atoms are labelled by the residue name QUA.

-  **INTERACTING_NEW.pdb**

   The same as INTERACTING.pdb, but it is created if the file INTERACTING.pdb
   is detected in the current working directory of the CPMD run.

-  **TRAJECTORY_INTERACTING**

   Contains the coordinates and the velocities (in TRAJECTORY format) of the
   atoms listed in INTERACTING.pdb. The format is the same as in the files
   TRAJECTORY and MOVIE, hence frames belonging to different runs are separated
   by the line <<<<< NEW DATA >>>>>.

   The atoms in this file do not necessarily coincide with the NN atoms, that
   are written at every update of the pair list in the file
   INTERACTING_NEW.pdb.

-  **MOVIE_INTERACTING**

   The MOVIE-like file corresponding to TRAJECTORY_INTERACTING.

-  **ESP**

   Contains the ESP charges of the QM atoms in CPMD order (the corresponding
   Gromos numbers can be found in QMMM_ORDER). The first column is the frame
   number.

-  **EL_ENERGY**

   Contains the electrostatic interaction energy. First column: frame number.
   Second column: total electrostatic interaction energy. Other columns:
   interaction energy of the NN atoms with the QM system; interaction energy
   with the ESP coupled atoms; multipolar interaction energy; electrostatic
   interaction energy evaluated using the classical force field charges for the
   QM atoms.

-  **MULTIPOLE**

   Contains, for every frame (specified in the first column), the three
   components of the dipole D(ix) and the five independent components of the
   quadrupole Q(ix,jx) of the quantum system in a.u. The order is:
   D(1),D(2),D(3),Q(1,1),Q(2,2),Q(1,2),Q(1,3),Q(2,3).

-  **MM_CELL_TRANS**

   Contains, the Trajectory of the re-centering offset for the QM-box. The
   first column is the frame number (NFI) followed by the x-, y-, and
   z-component of the cell-shift vector.

.. _sec:qmmm-cut-bonds:

Hydrogen Capping vs. Link Atoms
-------------------------------

Whenever the QM/MM-boundary cuts through an existing bond, special care has to
be taken to make sure that the electronic structure of the QM-subsystem is a
good representation of an all-QM calculation and also the structure in the
boundary region is preserved. So far, two methods methods are available to do
this: using special link-atom pseudopotentials and hydrogen capping.

Link Atom

The simplest way is to use a link-atom pseudopotential. In the simplest case,
this would be a scaled down pseudopotential with the required valence change
(e.g. ZV=1 when cutting through a single carbon-carbon bond). However in this
case it is required to constrain the distance between the link atom and the
(full-QM) neighbor atom to the (full-QM) equilibrium distance, to preserve the
electronic structure in the center of the QM subsystem. You should be aware of
the fact, that this is a rather crude approximation and that the constraint
will create a small imbalance in the forces between the QM and MM subsystems,
that can result in a drift in the total energy, if the length of the constraint
is badly chosen.

A more rigorous approach would be to use an optimized pseudopotential
constructed with the method described in ref. [opt-ecp04]_, that
should take care of the need for the constraint.

Hydrogen Capping

An alternative way would be to use the :kw:`CAPPING` flag in order
to introduce additional (dummy) hydrogen atoms to saturate the dangling bonds.
These capping hydrogen atoms have to be hidden from the MM hamiltonian so the
Gromos :kw:`INPUT` and :kw:`TOPOLOGY` files have to be
modified for subsequent runs, in addition to adding an explicit
:kw:`EXCLUSION` LIST to the cpmd input. The whole procedure is a
bit complicated so here is a short protocol of the required steps.

-  Set up a normal QM/MM run as for using link atoms with the additional
   keyword :kw:`CAPPING` and instead of the link-atom potential use
   a hydrogen potential with the additional flag ``ADD_H``. Note that you have
   to provide the correct number of hydrogens, but no atom index number.

-  Run a short MD (a couple of steps) and use the resulting ``CRD_FIN.grm``
   file (under a different name) in the :kw:`COORDINATES`
   section.

-  Modify the Gromos input file to match the new coordinate file.

   -  Increase the number of atoms per solute molecule in the SUBMOLECULES
      section.

   -  increase the total number of atoms in the FORCE section.

-  Modify the Gromos topology file to match the new coordinate file.

   -  Add a DUM atom type at the end of the ATOMTYPENAME section if not already
      present (and increase NRATT accordingly).

   -  if you have added a new atom type, you have to add the corresponding
      entries in the LJPARAMETERS section as well. Since the capping hydrogen
      atoms should be invisible from the MM Hamiltonian all Lennard-Jones
      parameters are set to 0.0 for those new entries. This section is a
      triangular matrix, so you have to add NRATT lines (and increase NRATT2
      accordingly).

   -  Add new residues named DUM (one for each capping hydrogen) to the RESNAME
      section (and increase NRAA2).

   -  In the SOLUTEATOM section you have to increase NRP and add the dummy
      hydrogens at the end of the solute. The structure of the entry is:
      ``<atom nr> <residue nr> <atom type name> <vdw type index> <mass> <charge> 1 0``
      and a single ’0’ on the next line. Use a mass of 1.008 and a charge of
      0.000.

-  Modify the CPMD input file.

   -  Make sure that the :kw:`TOPOLOGY`, :kw:`INPUT`, and
      :kw:`COORDINATES` keywords in the &QMMM section match the
      newly created or modified files.

   -  Add the capping hydrogens to the &ATOMS section as normal QM atoms, but
      add the DUMMY flag to the pseudopotential line.

   -  Build an :kw:`EXCLUSION` LIST entry that lists for each
      capping atom the respective QM/MM atoms pairs that should be excluded
      from the electrostatic coupling (all other MM interactions are set to
      zero already in the topology file). For consistency only full charge
      groups should be excluded. In the supplementary material should be a
      script genexcl.tcl which can help you in building that list (it needs the
      modified Gromos coordinate and topology file as well as the QMMM_ORDER
      file as input).

   -  Update the :kw:`ARRAYSIZES ... END ARRAYSIZES` entry to match the new
      topology.

With the three modified files you should be able to run a regular QM/MM run.
Note, that you may have to update the exclusion list occasionally, depending on
your system and that you should pick the bond(s) to cut very carefully.

.. _sec:what-type-qmmm:

What type of QM/MM calculations are available?
----------------------------------------------

The QM/MM interface only supports a subset of the functionality of CPMD. Please
note, that although there are some tests and warnings included into the source
code, not every job that runs without a warning will be automatically correct.
So far, the interface code requires the use of norm-conserving
pseudopotentials. Tested and supported job types are: :kw:`MOLECULAR DYNAMICS` (CLASSICAL, CP and BO), :kw:`OPTIMIZE WAVEFUNCTION`, :kw:`KOHN-SHAM ENERGIES`, and :kw:`ELECTRONIC SPECTRA`. Supported are closed shell systems as well
as :kw:`LSD` and :kw:`LOW SPIN EXCITATION`
calculations.

:kw:`OPTIMIZE GEOMETRY` is experimental and currently
supports optimization of the QM atom positions only. Use of the linear scaling
geometry optimizer (:kw:`LBFGS`) is highly recommended and the
currently also the default.

:kw:`PROPERTIES` calculations with QM/MM are experimental. Most
properties (WF projection, population analysis, localization) that only need
the plain QM wavefunction work.

:kw:`LINEAR RESPONSE` calculations are currently at an
twofold experimental status. Both, the isolated system setup
(:kw:`SYMMETRY` 0) and the QM/MM coupling of the response
calculations itself are not yet fully tested.

Options that are known to be incompatible with :kw:`QMMM` are
:kw:`VIBRATIONAL ANALYSIS`, :kw:`PATH INTEGRAL`, and all calculations that require a
wavefunction optimization via a diagonalization method at some point.

.. _sec:forcematch-desc:

QM/MM Force Matching
--------------------

This tool allows the automated (re)parametrization of classical force fields
from QM/MM reference calculations via a force matching protocol as published in
[FM-maurer]_. Thereby only MM parameters among the atoms comprised
in the QM subsystem are reparametrized. In this first release VdW parameters
are excluded from the optimization and kept constant. Fitting of these
parameters will be a feature of a future release.

The jobs requires a QM/MM reference TRAJECTORY file and the corresponding
gromos topology, input and coordinate files and the cpmd input that were used
to generate the reference TRAJECTORY file. The initial topology can be a
reasonable guess and will be refined during the actual force matching
procedure. Currently the QM/MM reference trajectory has to be generated prior
to the force matching job.

The actual forcematching job is invoked by the :kw:`FORCEMATCH`
in &CPMD and the :kw:`FORCEMATCH ... END FORCEMATCH` block in &QMMM. Besides, the
&CPMD and &QMMM sections should contain sensible keywords and parameters for
high quality reference forces (e.g. convergence orbitals 1.0d-7).

| The parametrization protocol consists of a three-step process: First, the
  reference trajectory is read with a given stride. On each of the selected
  frames QM/MM reference forces (BO) are calculated. The forces on the atoms of
  the QM subsystem are stored (FM_REF_FORCES) along with the Hirshfeld charges
  (FM_REF_CHJ) as well as the electrostatic potential and field on the nearby
  MM atoms (FM_REF_PIP). Second, a set of atomic point charges that reproduce
  the electrostatic potential and forces that the QM system exerts on the
  surrounding classical atoms is derived. Third, the nonbonded contributions,
  computed with the charges obtained in the second step and given Lennard-Jones
  parameters, are subtracted from the total reference forces on the QM atoms.
  The remaining forces are assumed to be derived from bonded interactions. The
  parameters for bonded interactions (torsions, bending and bonds) are thus
  adjusted in order to reproduce the remaining forces. See reference
  [FM-maurer]_ for details. An updated topology file FMATCH.top is
  written at the end of the run. To check the quality of the fitting procedure
  a section with the absolute and relative force RMS per atom is printed at the
  end to standard output.
| Files generated by the force matching code (Some of the following default
  filenames can be changed via the respective keywords in the FORCEMATCH
  block):

-  **FMATCH.top**

   Updated topology file at the end of the job.

-  **FM_REF_CHJ**

   Hirshfeld charges on the QM atoms from the reference force calculations.
   Format: Two lines per frame. First line contains frame index from the
   original reference trajectory file. Second line gives the Hirshfeld charges
   on the QM atoms in cpmd ordering.

-  **FM_REF_PIP**

   Electrostatic potential and field on the NN atoms from the reference force
   calculations.

-  **FM_REF_FORCES**

   | For each frame extracted from the reference TRAJECTORY file and for which
     QM/MM forces were calculated, the QM/MM forces on the QM atoms are dumped
     into this file. One line per frame with the original frame index and the
     number of QM atoms. Then for each QM atom in cpmd ordering:
   | atom index,x,y,z,fx,fy,fz

All force matching related information written to standard output are labeled
with ’ fm ’. After the initialization the reference trajectory file is parsed
and the line ’fm extracting total number of frames for SPs’ is printed. For
each frame you should find the following lines: ’frame number’, ’computing
reference forces’ and ’Total nr. of iterations:’. The beginning of the second
part of the force matching protocol is marked with the line ’Reading values for
charge fitting from file’. At the end the RMS deviation of the charges,
electrostatic potential and field are printed. The third part starts with ’Will
now loop over reference frames again’ to compute the non-bonded interactions.
After ’Done with classical loop’ the covalent parameters are fitted and you can
monitor the change of the absolute and relative RMS deviation from the
reference covalent forces during the optimization. ’Optimization successful’
indicates the end of the fitting of the covalent parameters, FMATCH.top is
written and, finally, the total (non-bonded plus covalent) forces are
calculated with the updated topology to get the RMS deviation of the total
force. The force matching related output ends with a block containing
’computing RMS per atom’.
