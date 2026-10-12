.. _input-qmmm:

``&QMMM`` section
=================

Input for the Gromos QM/MM interface (see :ref:`sec:qmmm`). **Required** if the :kw:`QMMM` keyword is given in the ``&CPMD`` section.

The section starts with ``&QMMM`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &QMMM

Keywords
--------

.. keyword:: ADD_HYDROGEN
   :target: ADD-HYDROGEN
   :section: &QMMM
   :notes: (mandatory)

   This keyword is used to add hydrogens to the QM system if a united atom
   topology is used (like in Gromos). On the next line the number of atoms to be
   “hydrogenized” has to be given and in the line following that, the
   corresponding gromos atom numbers. A number of hydrogens consistent with the
   hybridization of the “hydrogenized” carbons are added.

.. keyword:: AMBER
   :section: &QMMM
   :notes: (mandatory)

   An Amber functional form for the classical force field is used. In this case
   coordinates and topology files as obtained by Amber have to be converted in
   Gromos format just for input/read consistency. This is done with the tool
   amber2gromos availabe with the CPMD/QMMM package.

   This keyword is mutually exclusive with the :kw:`GROMOS` keyword
   (which is used by default).

.. keyword:: ARRAYSIZES ... END ARRAYSIZES
   :section: &QMMM
   :notes: (mandatory)

   Parameters for the dimensions of various internal arrays can be given in this
   block. The syntax is one label and the according dimension per line. The
   suitable parameters can be estimated using the script ``estimate_gomos_size``
   bundled with the QM/MM-code distribution. Example:

   ::

       ARRAYSIZES
         MAXATT 20
         MAXAA2 17
         MXEX14 373
       END ARRAYSIZES

.. keyword:: BOX TOLERANCE
   :section: &QMMM
   :notes: (mandatory)

   The value for the box tolerance is read from the next line. In a QM/MM
   calculation the size of the QM-box is fixed and the QM-atoms must not come to
   close to the walls of this box. On top of always recentering the QM-box around
   the center of the distribution of the atoms, CPMD prints a warning message to
   the output when the distribution extends too much to fit into the QM-box
   properly anymore. This value may need to be adjusted to the requirements of the
   Poisson solver used (see section :ref:`hints:symm0`).

   **Default** value is 8 a.u.

.. keyword:: BOX WALLS
   :target: BOX WALLS QMMM
   :section: &QMMM
   :notes: (mandatory)

   The thickness parameter for soft, reflecting QM-box walls is read from the next
   line. This keyword allows to reverse the momentum of the particles
   (:math:`{\bf p}_I \rightarrow -{\bf p}_I`) when they reach the walls of the
   simulation supercell similar to the full quantum case, but acting along all the
   three directions :math:`x,y,z`. In the case this keyword is used in the &QMMM
   section,QM particles are reflected back in the QM box. Contrary to the normal
   procedure of re-centering the QM-box, a soft, reflecting confinement potential
   is applied if atoms come too close to the border of the QM
   box [box-walls]_. It is highly recommended to also use
   :kw:`SUBTRACT` COMVEL in combination with this feature. **NOTE:**
   to have your QM-box properly centered, it is best to run a short MD with this
   feature turned off and then start from the resulting restart with the soft
   walls turned on. Since the reflecting walls reverse the sign of the velocities,
   :math:`{\bf p}_I \to -{\bf p}_I` (:math:`I` = QM atom index), be aware that
   this options affects the momentum conservation in your QM subsystem.

   This feature is **disabled by default**

.. keyword:: CAPPING
   :section: &QMMM
   :notes: (mandatory)

   Add (dummy) hydrogen atoms to the QM-system to saturate dangling bonds when
   cutting between MM- and QM-system. This needs a special pseudopotential entry
   in the &ATOMS section (see section :ref:`sec:qmmm-cut-bonds` for more
   details).

.. keyword:: CAP_HYDROGEN
   :target: CAP-HYDROGEN
   :section: &QMMM
   :notes: (mandatory)

   same as :kw:`CAPPING`.

.. keyword:: COORDINATES
   :section: &QMMM
   :notes: (mandatory)

   On the next line the name of a Gromos96 format coordinate file has to be given.
   Note, that this file must match the corresponding input and topology files.
   Note, that in case of hydrogen capping, this file has to be modified to also
   contain the respective dummy hydrogen atoms.

.. keyword:: ELECTROSTATIC COUPLING
   :section: &QMMM
   :options: [LONG RANGE]
   :notes: (mandatory)

   The electrostatic interaction of the quantum system with the classical system
   is explicitly kept into account for all classical atoms at a distance
   :math:`r \leq` :kw:`RCUT_NN <RCUT-NN>` from any quantum atom and for all the
   MM atoms at a distance of
   :kw:`RCUT_NN <RCUT-NN>` :math:`< r \leq` :kw:`RCUT_MIX <RCUT-MIX>` and a
   charge larger than :math:`0.1 e_0` (NN atoms).

   MM-atoms with a charge smaller than :math:`0.1 e_0` and a distance of
   :kw:`RCUT_NN <RCUT-NN>` :math:`< r \leq` :kw:`RCUT_MIX <RCUT-MIX>` and all
   MM-atoms with
   :kw:`RCUT_MIX <RCUT-MIX>` :math:`< r \leq` :kw:`RCUT_ESP <RCUT-ESP>` are
   coupled to the QM system by a ESP coupling Hamiltonian (EC atoms).

   If the additional ``LONG RANGE`` keyword is specified, the interaction of the
   QM-system with the rest of the classical atoms is explicitly kept into account
   via interacting with a multipole expansion for the QM-system up to quadrupolar
   order. A file named ``MULTIPOLE`` is produced.

   If ``LONG RANGE`` is omitted the quantum system is coupled to the classical
   atoms not in the NN-area and in the EC-area list via the force-field charges.

   If the keyword ``ELECTROSTATIC COUPLING`` is omitted, all classical atoms are
   coupled to the quantum system by the force-field charges (mechanical coupling).

   The files INTERACTING.pdb, TRAJECTORY_INTERACTING, MOVIE_INTERACTING,
   TRAJ_INT.dcd, and ESP (or some of them) are created. The list of NN and EC
   atoms is updated every 100 MD steps. This can be changed using the keyword
   :kw:`UPDATE LIST`.

   The default values for the cut-offs are RCUT_NN=RCUT_MIX=RCUT_ESP=10 a.u..
   These values can be changed by the keywords :kw:`RCUT_NN <RCUT-NN>`,
   :kw:`RCUT_MIX <RCUT-MIX>`, and :kw:`RCUT_ESP <RCUT-ESP>` with
   :math:`r_{nn} \leq r_{mix} \leq r_{esp}`.

.. keyword:: ESPWEIGHT
   :section: &QMMM
   :notes: (mandatory)

   The ESP-charg fit weighting parameter is read from the next line.

   **Default** value is :math:`0.1 e_0`.

.. keyword:: EXCLUSION
   :section: &QMMM
   :options: {GROMOS,LIST{NORESP}}
   :notes: (mandatory)

   Specify charge interactions that should be excluded from the QM/MM Hamiltonian.
   With the additional flag GROMOS, the exclusions from the Gromos topology are
   used. With the additional flag LIST, an explicit list is read from following
   lines. The format of that list has the number of exclusions in the first line
   and then the exclusions listed in pairs of numbers of the QM atom and the MM
   atom in Gromos ordering; the optional flag NORESP in this case requests usage
   of MM point charges for the QM atoms instead of the D-RESP charges (default).

.. keyword:: FLEXIBLE WATER
   :section: &QMMM
   :options: [ALL,BONDTYPE]
   :notes: (mandatory)

   Convert some solven water molecules into solute molecules and thus using a
   flexible potential.

   With the BONDTYPE flag, the three bond potentials (OH1, OH2, and H1H2) can be
   given as index in the BONDTYPE section of the Gromos topology file. Note that
   the **non-bonded** parameters are taken from the SOLVENATOM section of the
   :kw:`TOPOLOGY` file. **Default** is to use the values: 35, 35, 41.

   With the additional flag ALL this applies to all solvent water molecules,
   otherwise on the next line the number of flexible water molecules has to be
   given with the Gromos index numbers of their respective Oxygen atoms on the
   following line(s).

   | On successful conversion a new, adapted topology file, MM_TOPOLOGY, is
     written that has to be used with the :kw:`TOPOLOGY` keyword for
     subsequent restarts. Also the :kw:`INPUT` file has to be adapted: in
     the SYSTEM section the number of solvent molecules has to be reduced by the
     number of converted molecules, and in the SUBMOLECULES section the new solute
     atoms have to be added accordingly.
   | Example:

   ::

           FLEXIBLE WATER BONDTYPE
            4 4 5
            26
              32   101   188   284   308   359   407   476   506   680
             764   779   926  1082  1175  1247  1337  1355  1607  1943
            1958  1985  2066  2111  2153  2273

.. keyword:: FORCEMATCH ... END FORCEMATCH
   :section: &QMMM
   :notes: (mandatory)

   Input block for the QM/MM forcematching. A general description is given in
   section :ref:`sec:forcematch-desc`.

   **READ REF FORCES [FILE,COVALENT]**

   | Flag to read the QM/MM reference forces directly from the file FM_REF_FORCES,
     i.e. no QM/MM SPs are computed. **Default**: false. An alternative file name
     can be specified on the next line with the option **FILE**. With the option
     **COVALENT** covalent forces are read from the file FM_REF_COVFORCES.
   | **READ REF TRAJ [FILE]**

   Read reference trajectory from file TRAJECTORY_REF (or set the **FILE** option
   to read a non-default file name from the next line) with a given stride and
   compute single points on the respective frames.

   **RESTART SP**

   | If in a previous force matching run not all of the SPs could be computed
     (e.g. limited wall time) this flag indicates cpmd to restart the SP
     calculations. The FM_REF\* files from the previous run have to be present and
     they will be appended. With this option make sure that the frames contained
     in the already existing FM_REF\* files are consistent. **Default**: false.
   | **READ REF STRIDE**

   Stride to apply when reading the TRAJECTORY_REF file is read from the next
   line. Default=1, i.e. every frame is used for the SP calculations.

   **TOPOL OUT**

   Filename for the final topology file. **Default**: FMATCH.top.

   **INITWF [OFF]**

   Generate an initial guess for the wfkt for the SP calculations based on AOs
   (default). With the **OFF** option the wfkt of the previous frame is used as an
   intial guess.

   **CHARGES [ONLY,NO],[FIX]**

   | Charge fitting is on by default and can be switched off with the **NO**
     option. In this case the charges from the initial topology will not be
     modified. **ONLY** will let the program stop after the charge fitting and the
     other parameters are not updated.
   | With the **FIX** option target values for the restraints in the charge
     fitting on specific atoms can be specified by the user. Usually the charges
     are restraint to the respective Hirshfeld values. On the next line the number
     of charges to be fixed has to be given and then the corresponding number of
     lines with: gromos index     charge.

   **MV**

   Weight on the potential in the charge fitting. **Default**\ =0.1.

   **WF**

   Weight on the field in the charge fitting. **Default**\ =0.0.

   **WQ INDIVIDUAL**

   Weights on the charge restraints can be given individually here. From the next
   line the total number of individual weights is read. Then the lines with:
   gromos index     weight.

   **WQ GENERAL**

   The weight for all the charge restraints that were not specified by individual
   weights can be given on the next line. **Default**\ =0.1.

   **WTOT**

   Weight of the total charge contribution in the charge fitting.
   **Default**\ =1.0E7.

   **EQUIV**

   Specify equivalent atoms. Syntax:

   ``EQUIV``

   | ``n_equiv``
   | ``atom1   atom4``
   | ``atom1   atom3``
   | ``...``
   | ``atom5   atom7``

   There are *``n_equiv``* equivalencies specified (``n_equiv`` lines are read
   from the input). For each pair of equivalencies the gromos indexes have to be
   specified on one separate line. The lower index has to be given first!

   If an atom is equivalent to more then one other atom. E.g. atom1,atom3 and
   atom4 are equivalent. Then this has to be encoded by:

   | ``atom1   atom3``
   | ``atom1   atom4``

   and not by:

   | ``atom1   atom4``
   | ``atom3   atom4``

   where atom1 has a lower gromos index then atom3 and atom3 has a lower one then
   atom4. Per default no equivalencies are assumed.

   **OPT FC ONLY**

   Serves as a flag to remove the equilibrium values of the bonded interactions
   from the list of fitted parameters. I. e. only force constants are fitted for
   the bonded interactions.

   **NO BONDS**

   Do not fit bonds. **Default**\ =false.

   **NO ANGLES**

   Do not fit angles. **Default**\ =false.

   **NO DIHEDRALS**

   Do not fit dihedrals. **Default**\ =false.

   **NO IMPROPERS**

   Do not fit improper dihedrals. **Default**\ =false.

   **MAXITER**

   Give on the next line the maximal number of iterations for the non-linear
   fitting procedure of the bonded interactions. **Default**\ =500.

   **COMPUTE RMS [NO]**

   Per default the RMS on the forces is computed after the fitting has been
   completed. Switch it off with the **NO** option. Example:

   ::

       FORCEMATCH
         READ REF TRAJ FILE
           TRAJECTORY_REF
         READ REF STRIDE
           10
         WV
         1.0
         WOT
         1000000.0
         WQ GENERAL
         0.1
       END FORCEMATCH

.. keyword:: GROMOS
   :section: &QMMM
   :notes: (mandatory)

   A Gromos functional form for the classical force field is used (this is the
   default).

   This keyword is mutually exclusive with the :kw:`AMBER` keyword.

.. keyword:: HIRSHFELD
   :section: &QMMM
   :options: [ON,OFF]
   :notes: (mandatory)

   With this option, restraints to Hirshfeld charges [Hirshfeld77]_
   can be turned on or off

   **Default** value is ON.

.. keyword:: INPUT
   :section: &QMMM
   :notes: (mandatory)

   On the next line the name of a Gromos input file has to be given. A short
   summary of the input file syntax and some keywords are in section
   :ref:`sec:qmmm-gromos-inp`. Note, that it has to be a correct input
   file, even though many options do not apply for QM/MM runs.

.. keyword:: MAXNN
   :section: &QMMM
   :notes: (mandatory)

   Then maximum number of NN atoms, i.e. the number of atoms coupled to the QM
   system via :kw:`ELECTROSTATIC COUPLING` is read from
   the next line. (Note: This keyword was renamed from MAXNAT in older versions of
   the QM/MM interface code to avoid confusion with the MAXNAT keyword in the
   :kw:`ARRAYSIZES ... END ARRAYSIZES` block.)

   **Default** value is 5000.

.. keyword:: NOSPLIT
   :section: &QMMM
   :notes: (mandatory)

   If the program is run on more than one node, the MM forces calculation is
   performed on all nodes. Since the MM part is not parallelized, this is mostly
   useful for systems with a small MM-part and for runs using only very few nodes.
   Usually the QM part of the calculation needs the bulk of the cpu-time in the
   QM/MM.

   This setting is the default. See also under :kw:`SPLIT`.

.. keyword:: RCUT_ESP
   :target: RCUT-ESP
   :section: &QMMM
   :notes: (mandatory)

   The cutoff distance for atoms in the ESP-area
   (:math:`r_{mix} < r \leq r_{esp}`) is read from the next line. (see
   :kw:`ELECTROSTATIC COUPLING` for more details).

   **Default** value is 10 a.u.

.. keyword:: RCUT_MIX
   :target: RCUT-MIX
   :section: &QMMM
   :notes: (mandatory)

   The cutoff distance for atoms in the intermediate region
   (:math:`r_{nn} < r \leq r_{mix}`) is read from the next line. (see
   :kw:`ELECTROSTATIC COUPLING` for more details).

   **Default** value is 10 a.u.

.. keyword:: RCUT_NN
   :target: RCUT-NN
   :section: &QMMM
   :notes: (mandatory)

   The cutoff distance for atoms in the nearest neighbor region from the QM-system
   (:math:`r \leq r_{nn}`) is read from the next line. (see :kw:`ELECTROSTATIC COUPLING` for more details).

   **Default** value is 10 a.u.

.. keyword:: RESTART TRAJECTORY
   :section: &QMMM
   :options: [FRAME {num},FILE '{fname}',REVERSE]
   :notes: (mandatory)

   Restart the MD with coordinates and velocities from a previous run. With the
   additional flag FRAME followed by the frame number the trajectory frame can be
   selected. With the flag FILE followed by the name of the trajectory file, the
   filename can be set (Default is TRAJECTORY). Finally the flag REVERSE will
   reverse the sign of the velocities, so the system will move backwards from the
   selected point in the trajectory.

.. keyword:: SAMPLE INTERACTING
   :section: &QMMM
   :options: [OFF,DCD]
   :notes: (mandatory)

   The sampling rate for writing a trajectory of the interacting subsystem is read
   from the next line. With the additional keyword OFF or a sampling rate of 0,
   those trajectories are not written. The coordinates of the atoms atoms
   contained in the file INTERACTING.pdb are written, in the same order, on the
   file TRAJECTORY_INTERACTING every. If the :kw:`MOVIE` output is turned
   on, a file MOVIE_INTERACTING is written as well. With the additional keyword
   DCD the file TRAJ_INT.dcd is also written to. if the sampling rate is negative,
   then **only** the TRAJ_INT.dcd is written.

   **Default** value is 5 for MD calculations and OFF for others.

.. keyword:: SPLIT
   :section: &QMMM
   :notes: (mandatory)

   If the program is run on more than one node, the MM forces calculation is
   performed on a separate node. This is mostly useful for systems with a large
   MM-part and runs with many nodes where the accumulated time used for the
   classical part has a larger impact on the performance than losing one node for
   the (in total) much more time consuming QM-part.

   **Default** is :kw:`NOSPLIT`.

.. keyword:: TIMINGS
   :section: &QMMM
   :notes: (mandatory)

   Display timing information about the various parts of the QM/MM interface code
   in the output file. Also a file ``TIMINGS`` with even more details is written.
   This option is off by **default**.

.. keyword:: TOPOLOGY
   :section: &QMMM
   :notes: (mandatory)

   On the next line the name of a Gromos topology file has to be given. Regardless
   of the force field, this topology file has to be in Gromos
   format\ [gromos96]_. Topologies created with Amber can be
   converted using the respective conversion tools shipped with the interface
   code. A short summary of the topology file syntax and some keywords are in
   section :ref:`sec:qmmm-gromos-inp`.

   **Other keywords:**

.. keyword:: UPDATE LIST
   :section: &QMMM
   :notes: (mandatory)

   On the next line the number of MD steps between updates of the various lists of
   atoms for :kw:`ELECTROSTATIC COUPLING` is given. At
   every list update a file INTERACTING_NEW.pdb is created (and overwritten).

   **Default** value is 100.

.. keyword:: VERBOSE
   :section: &QMMM
   :notes: (mandatory)

   The progress of the QM/MM simulation is reported more verbosely in the output.
   This option is off by **default**.

.. keyword:: WRITE LOCALTEMP
   :section: &QMMM
   :options: [STEP {nfi_lt}]
   :notes: (mandatory)

   The Temperatures of the QM subsystem, the MM solute (without the QM atoms) and
   the solvent (if present) are calculated separately and written to the standard
   output and a file ``QM_TEMP``. The file has 5 columns containing the QM
   temperature, the MM temperature, the solvent temperature (or 0.0 if the solvent
   is part of the solute), and the total temperature in that order. With the
   optional parameters STEP followed by an integer, this is done only every
   ``nfi_lt`` timesteps.
