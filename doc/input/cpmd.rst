.. _input-cpmd:

``&CPMD`` section
=================

General control parameters for the calculation (**required**).

The section starts with ``&CPMD`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &CPMD

Keywords
--------

.. keyword:: ALEXANDER MIXING
   :section: &CPMD

   Mixing used during optimization of geometry or molecular dynamics. Parameter
   read in the next line.

   **Default** value is **0.9**

.. keyword:: ALLTOALL
   :section: &CPMD
   :options: {SINGLE,DOUBLE}

   Perform the matrix transpose (AllToAll communication) in the 3D FFT using
   single/double precision numbers. Default is to use double precision numbers.

.. keyword:: ANDERSON MIXING
   :section: &CPMD
   :alternatives: N=n

   Anderson mixing for the electronic density during self-consistent iterations.
   In the next line the parameter (between 0 and 1) for the Anderson mixing is
   read.

   **Default** is **0.2**.

   With the additional option :math:`N=n` a mixing parameter can be specified for
   different threshold densities. :math:`n` different thresholds can be set. The
   program reads :math:`n` lines, each with a threshold density and an Anderson
   mixing parameter.

.. keyword:: ANNEALING
   :section: &CPMD
   :options: {IONS,ELECTRONS,CELL}

   Scale the ionic, electronic, or cell velocities every time step. The scaling
   factor is read from the next line.

.. keyword:: BENCHMARK
   :section: &CPMD

   This keyword is used to control some special features related to benchmarks. If
   you want to know more, have a look in the source code.

.. keyword:: BERENDSEN
   :section: &CPMD
   :options: {IONS,ELECTRONS,CELL}

   Use a simple Berendsen-type thermostat\ [Berendsen84]_ to control
   the respective temperature of ions, electrons, or cell. The target temperature
   and time constant :math:`\tau` (in a.u.) are read from the next line.

   These thermostats are a gentler alternative to the
   :kw:`TEMPCONTROL` mechanism to thermalize a system. For
   production runs, please use the corresponding :kw:`NOSE` or
   :kw:`LANGEVIN` thermostats, as the Berendsen scheme does not
   represent any defined statistical mechanical ensemble.

.. keyword:: BFGS
   :section: &CPMD

   Use a quasi-Newton method for optimization of the ionic positions. The
   approximated Hessian is updated using the Broyden-Fletcher-Goldfarb-Shano
   procedure [Fletcher80]_.

.. keyword:: BICANONICAL ENSEMBLE
   :section: &CPMD
   :options: {CHEMICALPOTENTIAL, XWEIGHT}
   :alternatives: [INFO]

   Use the bicanonical ensemble AIMD for simulation of open
   systems [Frenzel2017]_.

   Note, requesting this runtype cpmd must be executed on more then one MPI
   thread.

   Within the bicanonocal ensemble one out of two sampling modes must be specified
   in the input:

   |
   | **CHEMICALPOTENTIAL** requests to sample the bicanonical ensemble at constant
     chemical potential (:math:`\mu^\dagger_1`) or
   | **XWEIGHT** requests to sample the bicanonical ensemble at constant
     bicanonical weight (:math:`x_1`).
   | **Default** is to use no INFO printing at each AIMD step.
   | Two parameters read in from next line:
   | [**chemicalPotential**\ :math:`~|~`\ **xWeight**] [**temperature**]
   | **chemicalPotential** (:math:`\mu^\dagger_1`) if CHEMICALPOTENTIAL run type
     was requested.
   | A reasonable value can be obtained taking the ensemble average of
     :math:`\langle
         \mu^\dagger_1 \rangle_{x_1=0.5}` of a converged bicanonical AIMD
     simulation using the XWEIGHT run type (see Ref. [Frenzel2017]_).
   | **xWeight** (:math:`x_1`) if XWEIGHT run type was requested.
   | The value of **xWeight** must be within the interval :math:`[0,1]`. Setting
     the value of **xWeight** to 0 and 1 corresponds to the canonical limits of
     canonical system 2 and system 1, *i.e.* small and the large system,
     respectively.
   | **temperature** specifies the value of the temperature used to calculate the
     :math:`\gamma_1` prefactor (see Ref. [Frenzel2017]_).
   | Preparation of the INPUT (see also CPMD regtests BICAN):
   | The simulation setup requires two input files with the second one having the
     same name plus the extension "\_2". The two files are for the large and the
     small canonical system, system 1 and system 2, respectively. The two input
     files must differ in one atom. This atom must be declared as the \_last\_
     atomic species in the input of system 1. Thus the INPUT for system 2 is
     created from 1 by removing that species and then apply optional changes,
     *e.g.* modify the number of electrons in that system via the CHARGE keyword.
     Note, positions and velocities of the atomic nuclei must be the same in both
     systems. A restart from preoptimized electronic wavefunctions of the
     individual canonical systems is recommended (see Sec. :ref:`FILES` and
     CPMD regtests BICAN/mdX0_50, BICAN/c-cnf1, BICAN/c-cnf2). The output file of
     system 2 is written on OUTPUT_CNF2 while that of system 1 is stdout. The
     suffixes "\_CNF1" and "\_CNF2" are used to specify the canonical systems for
     the other cpmd files. For example the content of the file ENERGIES of a
     canoncal cpmd run is written to the files ENERGIES_CNF1 and ENERGIES_CNF2. In
     the bicanonical AIMD the ENERGIES then stores: nfi, chemical potential,
     xweight, :math:`\Delta U`, and :math:`\mu_{x_1=0.5}`. Where :math:`\Delta U`
     is the difference of the total energy of system 1 and system 2 at timestep
     nfi.

   Note: For the moment bicanonical ensemble is only implemented for MOLECULAR
   DYNAMICS. Other MD types BO, EH, PT, CLASSICAL or FILE as well as runtypes
   QMMM, path integral or surface hopping are not supported.

.. keyword:: BLOCKSIZE STATES
   :section: &CPMD

   Parameter read in from next line.

   *NSTBLK*

   Defines the minimal number of states used per processor in the distributed
   linear algebra calculations.

   **Default** is to equally distribute states over all processors.

.. keyword:: BOGOLIUBOV CORRECTION
   :section: &CPMD
   :alternatives: [OFF]

   Computes the Bogoliubov correction for the energy of the Trotter approximation
   or not.

   **Default** is **no Bogoliubov correction**.

   The keyword has to appear after :kw:`FREE ENERGY FUNCTIONAL`.

.. keyword:: BOX WALLS
   :section: &CPMD

   The thickness parameter for soft, reflecting QM-box walls is read from the next
   line. This keyword allows to reverse the momentum of the particles
   (:math:`{\bf p}_I \rightarrow -{\bf p}_I`) when they reach the walls of the
   simulation supercell in the case in which no periodic boundary conditions are
   applied. Specifically, in the unidimensional surface-like case, molecules
   departing from the surface are reflected back along the direction orthogonal to
   the surface, whereas in the bidimensional polymer-like case, they are reflected
   back in the two dimensions orthogonal to the "polymer" axis. Warning: This
   procedure, although keeping your particles inside the cell, affect the momentum
   conservation.

   This feature is **disabled by default**

.. keyword:: BROYDEN MIXING
   :section: &CPMD

   Parameters read in from next line.

   *BROYMIX, ECUTBROY, W02BROY, NFRBROY, IBRESET, KERMIX*

   These mean:

   *BROYMIX*:

   Initial mixing, e.g. :math:`0.1`; **default** value is **0.5**

   *ECUTBROY:*

   Cutoff for Broyden mixing. **DUAL*ECUT** is the best choice and the **default**

   *W02BROY:*

   :math:`w_0^2` parameter of Johnson [Johnson88]_. **Default**
   **0.01**

   *NFRBROY:*

   Number of Anderson mixing steps done before Broyden mixing. **Default** is
   **0**

   *IBRESET:*

   Number of Broyden vectors. :math:`5` is usually a good value and the default.

   *KERMIX:*

   Kerker mixing according to the original definition of
   Ref. [Kerker]_. By default the mixing parameter is set to 0.

   You can also specify some parameters with the following syntax:

   **[BROYMIX=**\ *BROYMIX*\ **]** **[ECUTBROY=**\ *ECUTBROY*\ **]**

   **[W02BROY=**\ *W02BROY*\ **]** **[NFRBROY=**\ *NFRBROY*\ **]**

   **[IBRESET=**\ *IBRESET*\ **]**

   **[KERMIX=**\ *KERMIX*\ **]**

   Finally, you can use the keyword **DEFAULT** to use the default values.

.. keyword:: CAYLEY
   :section: &CPMD

   Used to propagate the Kohn-Sham orbitals in :kw:`MOLECULAR DYNAMICS` EH and :kw:`PROPAGATION SPECTRA`. At present is the only propagation scheme
   available.

.. keyword:: CDFT
   :section: &CPMD
   :alternatives: [NEWTON, DEKKER],[SPIN, ALL, PCGFI, RESWF, NOCCOR, HDA[AUTO, PHIOUT, PROJECT]]

   The main switch for constrained DFT. Parameters :math:`N_\text{c}`,
   :math:`V_\text{init}`, and MAXSTEP are read from the next line.

   **NEWTON**, **DEKKER** (**off**) are switches to enable either the Newton or
   the Dekker optimisation scheme for the constraint. If neither of those are set
   a simple gradient scheme is used.

   **SPIN** (**off**) if activated the constraint will act on the spin density
   instead of the charge density. This may help against excessive spin
   contamination.

   **ALL** (**off**) activates dual spin and charge constraint, all inputs for
   :math:`N_\text{c}` and :math:`V_\text{init}` have to be given twice (first for
   charge then for spin)

   **PCGFI** (**off**) instructs CPMD to do PCG for the first V optimisation cycle
   regardless of the choice of optimiser.

   **RESWF** (**off**) if activated this switch re-initialises the wavefunction
   after each :math:`V` optimisation step. This is useful if the wavefunction
   convergence between the optimisation steps is slow. Usage in conjunction with
   :kw:`INITIALIZE WAVEFUNCTION` RANDOM may help.

   **NOCCOR** (**off**) if activated this switch turns off cutoff correction for
   the forces.

   **HDA** (**off**) if activated this switch turns on the calculation of the
   transition matrix element between the constrained states given by
   :math:`N_\text{c}` and :math:`\hat{N}_\text{c}` which is then read from the
   second line. For this keyword to take effect the :kw:`OPTIMIZE WAVEFUNCTION` option has to be activated.

   Sub-options of **HDA**

   **AUTO** (**off**) if activated this switch lets CPMD choose the constraint
   values for the transition matrix calculation. :math:`N_\text{c}` is chosen from
   the initial charge distribution and :math:`\hat{N}_\text{c}=-N_\text{c}`. It
   might be a good idea to use :kw:`INITIALIZE WAVEFUNCTION` ATOMS and :kw:`ATOMIC CHARGES` (&ATOM section) so that CPMD initialises the
   wavefunction with the desired pseudo wavefunction.

   **PHIOUT** (**off**) if activated this switch tells CPMD to write out the
   overlap matrices :math:`\Phi_\text{AA},\Phi_\text{BB},\Phi_\text{AB},` and
   :math:`\Phi_\text{BA}` to the file PHI_MAT.

   **PROJECT** (**off**) if activated this switch lets CPMD read in two reference
   states from RESTART.REF1 and RESTART.REF2 after the actual HDA calculation in
   order to project the two constrained states on them and thus calculate the
   diabatic transition matrix element in an orthogonalised “dressed” basis.

   If CDFT is activated the program writes the current :math:`V` value to
   CDFT_RESTART every time the RESTART file is written.

.. keyword:: CENTER MOLECULE
   :section: &CPMD
   :alternatives: [OFF]

   The center of mass is moved/not moved to the center of the computational box in
   a calculation with the cluster option. This is only done when the coordinates
   are read from the input file.

.. keyword:: CHECK MEMORY
   :section: &CPMD

   Check sanity of all dynamically allocated arrays whenever a change in the
   allocation is done. By default memory is checked only at break points.

.. keyword:: CLASSTRESS
   :section: &CPMD

   Not documented.

.. keyword:: CMASS
   :section: &CPMD

   The fictitious mass of the cell in atomic units is read from the next line.

   **Default** value is **200**

.. keyword:: COMBINE SYSTEMS
   :section: &CPMD
   :alternatives: [REFLECT,NONORTH,SAB]

   Read in two wavefunctions from RESTART.R1 and RESTART.R2 and combine them into
   RESTART.1 which can then be used in an FODFT calculations. The option NONORTH
   disables orthogonalisation of the combined WF’s. Parameters NTOT1, NTOT2,
   NSUP1, NSUP2 are read from the next line.

   NTOT1/NTOT2 total number of electrons in state 1/2 (mandatory).

   NSUP1/NSUP2 number of alpha electrons in state 1/2 (only LSD).

   If the option REFLECT is given a fifth parameter (CM_DIR) is read and the WF
   given in RESTART.R2 will be either mirrored through the centre of the box
   (CM_DIR=0), mirrored through the central yz-plane of the box (CM_DIR=1) or if
   CM_DIR=4 mirrored through the central yz-plane and translated in x direction by
   CM_DR (sixth parameter to be read).

   If the option SAB is set, write out the overlap matrix element between orbitals
   K and L. Parameters K and L are read from the next line.

   After combining the wavefunctions CPMD will exit. For this option to work the
   RESTART option and :kw:`OPTIMIZE WAVEFUNCTION` have to
   be activated.

.. keyword:: COMPRESS
   :section: &CPMD
   :alternatives: [WRITEnn]

   Write the wavefunctions with nn bytes precision to the restart file.

   Possible choices are ``WRITE32``, ``WRITE16``, ``WRITE8`` and ``WRITEAO``.

   ``WRITE32`` corresponds to the compress option in older versions. ``WRITEAO``
   stores the wavefunction as a projection on atomic basis sets. The atomic basis
   set can be specified in the section &BASIS … &END. If this input section is
   missing a default basis from Slater type orbitals is constructed. See
   section :ref:`input:basis` for more details.

.. keyword:: CONJUGATE GRADIENTS
   :section: &CPMD
   :alternatives: [ELECTRONS, IONS, NOPRECONDITIONING]

   For the electrons, the keyword is equivalent to :kw:`PCG`. The
   ``NOPRECONDITIONING`` parameter only applies for electrons. For the ions the
   conjugate gradients scheme is used to relax the atomic positions.

.. keyword:: CONVERGENCE
   :section: &CPMD
   :alternatives: [ADAPT, ENERGY, CALFOR, RELAX, INITIAL]

   The adaptive convergence criteria for the wavefunction during a geometry
   optimization are specified. For more informations, see [LSCAL]_.
   The ratio *TOLAD* between the smallest maximum component of the nuclear
   gradient reached so far and the maximum allowed component of the electronic
   **gradient** is specified with **CONVERGENCE ADAPT**. This criterion is
   switched off once the value *TOLOG* given with **CONVERGENCE ORBITALS** is
   reached. By default, the adaptive gradient criterion is not active. A
   reasonable value for the parameter *TOLAD* is 0.02.

   If the parameter *TOLENE* is given with **CONVERGENCE ENERGY**, in addition to
   the gradient criterion for the wavefunction, the energy change between two
   wavefunction optimization cycles must be smaller than the energy change of the
   last accepted geometry change multiplied by *TOLENE* for the wavefunction to be
   considered converged. By default, the adaptive energy criterion is not active.
   It is particularly useful for **transition state search** with P-RFO, where the
   trust radius is based on the quality of energy prediction. A reasonable value
   for *TOLENE* is 0.05.

   To save CPU time, the gradient on the ions is only calculated if the
   wavefunction is almost converged. The parameter *TOLFOR* given with
   **CONVERGENCE CALFOR** is the ratio between the convergence criteria for the
   wavefunction and the criteria whether the gradient on the ions is to be
   calculated. **Default** value for *TOLFOR* is **3.0**.

   If the wavefunction is very slowly converging during a geometry optimization, a
   small nuclear displacement can help. The parameter *NSTCNV* is given with
   **CONVERGENCE RELAX**. Every *NSTCNV* wavefunction optimization cycles, the
   convergence criteria for the wavefunction are relaxed by a factor of two. A
   geometry optimization step resets the criteria to the unrelaxed values. By
   default, the criteria for wavefunction convergence are never relaxed.

   When starting a geometry optimization from an unconverged wavefunction, the
   nuclear gradient and therefore the adaptive tolerance of the electronic
   gradient is not known. To avoid the **full convergence** criterion to be
   applied at the beginning, a convergence criterion for the wavefunction of the
   initial geometry can be supplied with **CONVERGENCE INITIAL**. By default, the
   initial convergence criterion is equal to the full convergence criterion.

.. keyword:: CONVERGENCE
   :target: CONVERGENCE 2
   :section: &CPMD
   :alternatives: [ORBITALS, GEOMETRY, CELL]

   The convergence criteria for optimization runs is specified.

   The maximum value for the biggest element of the gradient of the wavefunction
   (**ORBITALS**), of the ions (**GEOMETRY**), or the cell (**CELL**) is read from
   the next line.

   **Default** values are **10\ :math:`^{-5}`** for the wavefunction,
   **5\ :math:`\times`\ 10\ :math:`^{-4}`** for the ions and **1.0** for the cell.
   For diagonalization schemes the first value is the biggest variation of a
   density component. **Defaults** are **10\ :math:`^{-3}`** and
   **10\ :math:`^{-3}`**.

.. keyword:: CONVERGENCE
   :target: CONVERGENCE 3
   :section: &CPMD
   :options: CONSTRAINT

   Set constraint convergence parameters. Parameters VCCON and VCCONU are read
   from the next line:

   VCCON :math:`\in \mathbb{R}_+` is the maximally allowed total deviation of the
   constraint from the desired value :math:`\text{N}_\text{c}`.

   VCCONU :math:`\in \mathbb{R}_+` is the upper bound for the deviation in MD
   runs, excess of which triggers a new optimisation of :math:`V`.

   **Defaults** are **10\ :math:`^{-5}`** and **VCCON**.

.. keyword:: CP_GROUPS
   :target: CP GROUPS
   :section: &CPMD

   Set the number of groups to be used in the calculation. Default is 1 group. The
   number of groups is read from the next line and shall be a divisor of the
   number of nodes in a parallel run.

.. keyword:: CZONES
   :section: &CPMD
   :alternatives: [SET]

   Activates convergence zones for the wavefunction during the :kw:`CDFT`
   constraint minimisation. If SET is set the parameters CZONE1, CZONE2, and
   CZONE3 are read from the next line and CZLIMIT1 and CZLIMIT2 from the line
   after.

   CZONE1 **:math:`10^{-3}`**,CZONE2 **:math:`10^{-4}`**,CZONE3
   **:math:`10^{-5}`** :math:`\in \mathbb{R}_+` are the orbital convergences in
   zones 1-3, respectively.

   CZLIMIT1 **0.3**, CZLIMIT2 **0.1** :math:`\in \mathbb{R}_+` define the
   boundaries between zone 1-2 and 2-3, respectively.

.. keyword:: DAMPING
   :section: &CPMD
   :options: {IONS,ELECTRONS,CELL}

   Add a damping factor :math:`f_{damp}(x) = - \gamma \cdot v(x)` to the ionic,
   electronic, or cell forces in every time step. The scaling factor
   :math:`\gamma` is read from the next line. Useful values depend on the employed
   masses are generally in the range :math:`5.0 \to 50.0`.

   Damping can be used as a more efficient alternative to
   :kw:`ANNEALING` for wavefunction, geometry or cell optimization
   (and particularly combinations thereof) of systems where the faster methods
   (e.g. :kw:`ODIIS`, :kw:`PCG`, :kw:`LBFGS`,
   :kw:`GDIIS`) fail to converge or may converge to the wrong state.

.. keyword:: DAVIDSON DIAGONALIZATION
   :section: &CPMD

   Use Davidson diagonalization scheme.\ [davidson75]_

.. keyword:: DAVIDSON PARAMETER
   :section: &CPMD

   This keyword controls the Davidson diagonalization routine used to determine
   the Kohn-Sham energies.

   The maximum number of additional vectors to construct the Davidson matrix, the
   convergence criterion and the maximum number of steps are read from the next
   line.

   **Defaults** are **10\ :math:`^{-5}`** and the same number as states to be
   optimized. If the system has 20 occupied states and you ask for 5 unoccupied
   states, the default number of additional vectors is 25. By using less than 25
   some memory can be saved but convergence might be somewhat slower.

.. keyword:: DCACP
   :section: &CPMD

   Activates the DCACP scheme  [opt-ecp04]_ (dispersion-corrected
   atom-centred potentials) in order to account for van der Waals-interactions.
   Advanced settings can be specified in the &VDW section.

.. keyword:: DEBUG CODE
   :section: &CPMD

   Very verbose output concerning subroutine calls for debugging purpose.

.. keyword:: DEBUG FILEOPEN
   :section: &CPMD

   Very verbose output concerning opening files for debugging purpose.

.. keyword:: DEBUG FORCES
   :section: &CPMD

   Very verbose output concerning the calculation of each contribution to the
   forces for debugging purpose.

.. keyword:: DEBUG MEMORY
   :section: &CPMD

   Very verbose output concerning memory for debugging purpose.

.. keyword:: DEBUG NOACC
   :section: &CPMD

   Do not read/write accumulator information from/to the :kw:`RESTART`
   file. This avoids putting timing information to the restart and makes restart
   files identical for otherwise identical runs.

.. keyword:: DIIS MIXING
   :section: &CPMD

   Use the direct inversion iterative scheme to mix density.

   Read in the next line the number of previous densities (NRDIIS) for the mixing
   (however not useful).

.. keyword:: DIIS MIXING
   :target: DIIS MIXING N=n
   :section: &CPMD
   :alternatives: [N=n]

   Like DIIS MIXING, but number of previous densities for the mixing can be
   specified as a function of the density.

   :math:`n` different thresholds for the density can be set. The program reads
   :math:`n` lines with a threshold density and a NRDIIS number (number of
   previous densities for the mixing). Numbers NRDIIS have to increase. If the
   NRDIIS is equal to 0, Anderson mixing is used. Very efficient is to use
   Anderson mixing and afterwards DIIS mixing.

.. keyword:: DIPOLE DYNAMICS
   :section: &CPMD
   :options: {SAMPLE,WANNIER}

   Calculate the dipole moment [vdb_berry]_ [resta]_ every *NSTEP*
   iteration in MD.

   *NSTEP* is read from the next line if the keyword SAMPLE is present.

   **Default** is **every** time step.

   The keyword **Wannier** allows the calculation of optimally localized Wannier
   functions\ [Marzari97]_ [PSil99]_ [berghold]_. The procedure used is
   equivalent (for single k-point) to Boys localization.

   The produced output is IONS+CENTERS.xyz, IONS+CENTERS, DIPOLE, WANNIER_CENTER
   and WANNIER_DOS. The localization procedure is controlled by the following
   keywords.

.. keyword:: DISTRIBUTE FNL
   :section: &CPMD

   The array ``FNL`` is distributed in parallel runs.

.. keyword:: DISTRIBUTED LINALG
   :section: &CPMD
   :options: {ON, OFF}

   Perform linear algebra calculations using distributed memory algorithms.
   Setting this option ON will also enable (distributed) initialization from
   atomic wavefunctions using a parallel Lanczos algorithm
   [distrib.lanczos.07]_. Note that distributed initialization is not
   available with **KPOINTS** calculations. In this case, initialization from
   atomic wavefunctions will involve replicated calculations.

   When setting **LINALG ON** the keyword :kw:`BLOCKSIZE STATES` becomes relevant (see entry). The number of
   :kw:`BLOCKSIZE STATES` must be an **exact divisor** of the
   number of :kw:`STATES`.

.. keyword:: ELECTRONIC SPECTRA
   :section: &CPMD

   Perform a TDDFT calculation [tddft_all]_ [tddft_pw]_ to determine the
   electronic spectra. See below under :ref:`Electronic Spectra <sec:ElectronicSpectra>` and under the other keywords for the
   input sections :doc:`&LINRES <linres>` and
   :doc:`&TDDFT <tddft>` for further options.

.. keyword:: ELECTROSTATIC POTENTIAL
   :section: &CPMD
   :alternatives: [SAMPLE=nrhoout]

   Store the electrostatic potential on file. The resulting file is written in
   platform specific binary format. You can use the cpmd2cube tool to convert it
   into a Gaussian cube file for visualization. Note that this flag automatically
   activates the :kw:`RHOOUT` flag as well.

   With the optional keyword **SAMPLE** the file will be written every *nrhoout*
   steps during an MD trajectory. The corresponding time step number will be
   appended to the filename.

.. keyword:: ELF
   :section: &CPMD
   :alternatives: [PARAMETER]

   Store the total valence density and the valence electron localization function
   ELF [Becke90]_ [silsav]_ [marx-savin-97]_ [homeofelf]_ on files. The
   default smoothing parameters for ELF can be changed optionally when specifying
   in addition the PARAMETER keyword. Then the two parameters “elfcut” and
   “elfeps” are read from the next line. The particular form of ELF that is
   implemented is defined in the header of the subroutine elf_utils.mod.F90.

   Note 1: it is a *very* good idea to increase the planewave cutoff and then
   specify “elfcut” :math:`=0.0` and “elfeps” :math:`=0.0` if you want to obtain a
   smooth ELF for a given nuclear configuration. In the case of a spin–polarized
   (i.e. spin unrestricted) DFT calculation (see keyword :kw:`LSD`) in
   addition the spin–polarized average of ELF as well as the separate
   :math:`\alpha`– and :math:`\beta`–orbital parts are written to the files
   LSD_ELF, ELF_ALPHA and ELF_BETA, respectively; see Ref. [Kohut96]_
   for definitions and further infos.

   Note 2: ELF does not make much sense when using Vanderbilt’s ultra-soft
   pseudopotentials!

.. keyword:: EMASS
   :section: &CPMD

   The fictitious electron mass in atomic units is read from the next line.

   **Default** is **400 a.u.**.

.. keyword:: ENERGYBANDS
   :section: &CPMD

   Write the band energies (eigenvalues) for k points in the file ENERGYBANDS.

.. keyword:: EXTERNAL POTENTIAL
   :section: &CPMD
   :options: {ADD}

   Read an external potential from file. With ADD specified, its effects is added
   to the forces acting on the ions.

.. keyword:: EXTRAPOLATE CONSTRAINT
   :section: &CPMD

   In a CDFT MD run extrapolate the next value for :math:`V` using a Lagrange
   polynomial. The order :math:`k` of the polynomial is read from the next line.
   **Default** is **k=5**, but it pays off to use the orderfinder.py python script
   on the ENERGIES file of a former run to estimate the optimal extrapolation
   order :math:`k_\text{opt}`.

.. keyword:: EXTRAPOLATE WFN
   :section: &CPMD
   :options: {STORE}

   Read the number of wavefunctions to retain from the next line.

   These wavefunctions are used to extrapolate the initial guess wavefunction in
   Born-Oppenheimer MD. This can help to speed up BO-MD runs significantly by
   reducing the number of wavefunction optimization steps needed through two
   effects: the initial guess wavefunction is much improved and also you do not
   need to converge as tightly to conserve energy. BO-MD without needs CONVERGENCE
   ORBITALS to be set to 1.0e-7 or smaller to maintain good energy conservation.

   With the additional keyword **STORE** the wavefunction history is also written
   to restart files. See :kw:`RESTART` for how to read it back.

.. keyword:: FILE FUSION
   :section: &CPMD

   Reads in two separate :kw:`RESTART` files for ground state and
   :kw:`ROKS` excited state and writes them into a single restart file.

   Required to start :kw:`SURFACE HOPPING` calculations.

.. keyword:: FILEPATH
   :section: &CPMD

   The path to the files written by CPMD (RESTART.x, MOVIE, ENERGIES, DENSITY.x
   etc.) is read from the next line. This overwrites the value given in the
   environment variable **CPMD_FILEPATH**. **Default** is the **current
   directory**.

.. keyword:: FINITE DIFFERENCES
   :section: &CPMD

   The step length in a finite difference run for vibrational frequencies
   (VIBRATIONAL ANALYSIS keywords) is read from the next line.

   With the keywords **COORD=**\ *coord_fdiff(1..3)* and **RADIUS=**\ *radius* put
   in the same line as the step length, you can specify a sphere in order to
   calculate the finite differences only for the atoms inside it. The sphere is
   centered on the position *coord_fdiff(1..3)* with a radius *radius* (useful for
   a point defect).

   **NOTE:** The the step length for the finite difference is **always** in Bohr
   (default is 1.0d-2 a.u.). The (optional) coordinates of the center and the
   radius are read in either Angstrom or Bohr, depending on whether the
   :kw:`ANGSTROM` keyword is specified or not.

.. keyword:: FIXRHO UPWFN
   :section: &CPMD
   :alternatives: [VECT LOOP WFTOL]

   Wavefunctions optimization with the method of direct inversion of the iterative
   subspace (DIIS), performed without updating the charge density at each step.
   When the orbital energy gradients are below the given tolerance or when the
   maximum number of iterations is reached, the KS energies and the occupation
   numbers are calculated, the density is updated, and a new wavefunction
   optimization is started. The variations of the density coefficients are used as
   convergence criterion. The default electron temperature is 1000 K and 4
   unoccupied states are added. Implemented also for k-points. Only one sub-option
   is allowed per line and the respective parameter is read from the next line.
   The parameters mean:

   *VECT*:

   The number of DIIS vectors is read from the next line. (ODIIS with 4 vectors is
   the default).

   *LOOP:*

   the minimum and maximum number of DIIS iterations per each wfn optimization is
   read from the following line. Default values are 4 and 20.

   *WFTOL:*

   The convergence tolerance for the wfn optimization during the ODIIS is read
   from the following line. The default value is :math:`10^{-7}`. The program
   adjusts this criterion automatically, depending on the convergence status of
   the density. As the density improves (when the density updates become smaller),
   also the wavefunction convergence criterion is set to its final value.

.. keyword:: FORCEMATCH
   :section: &CPMD

   Activates the QM/MM force matching procedure. This keywords requires the
   presence of a &QMMM ... &END section with a corresponding :kw:`FORCEMATCH ... END FORCEMATCH` block. See
   sections :ref:`sec:qmmm` and :ref:`sec:forcematch-desc` for more
   details.

.. keyword:: FREE ENERGY FUNCTIONAL
   :section: &CPMD

   Calculates the electronic free energy using free energy density
   functional [Alavi94]_ [PSil]_ [mbaops]_ from DFT at finite temperature.

   This option needs additional keywords (free energy keywords).

   By **default** we use **Lanczos diagonalization** with **Trotter
   factorization** and **Bogoliubov correction**. If the number of states is not
   specified, use :math:`N_{electrons}/2+4`.

.. keyword:: GDIIS
   :section: &CPMD

   Use the method of direct inversion in the iterative subspace combined with a
   quasi-Newton method (using BFGS) for optimization of the ionic
   positions [Csaszar84]_. The number of DIIS vectors is read from
   the next line.

   GDIIS with **5 vectors** is the **default** method in optimization runs.

.. keyword:: GSHELL
   :section: &CPMD

   Write a file **GSHELL** with the information on the plane waves for further use
   in S(q) calculations.

.. keyword:: HAMILTONIAN CUTOFF
   :section: &CPMD

   The lower cutoff for the diagonal approximation to the Kohn-Sham
   matrix [Tuckerman94]_ is read from the next line.

   **Default** is **0.5** atomic units.

   For variable cell dynamics only the kinetic energy as calculated for the
   reference cell is used.

.. keyword:: HARMONIC REFERENCE SYSTEM
   :section: &CPMD
   :alternatives: [OFF]

   Switches harmonic reference system integration [Tuckerman94]_
   on/off.

   The number of shells included in the analytic integration is controlled with
   the keyword :kw:`HAMILTONIAN CUTOFF`.

   By **default** this option is switched **off**.

.. keyword:: HESSCORE
   :section: &CPMD

   Calculates the partial Hessian after relaxation of the environment, equivalent
   to *NSMAXP=0* (**PRFO NSMAXP**).

.. keyword:: HESSIAN
   :section: &CPMD
   :alternatives: [DISCO,SCHLEGEL,UNIT,PARTIAL]

   The initial approximate **Hessian** for a **geometry optimization** is
   constructed using empirical rules with the DISCO [Fischer92]_ or
   Schlegel’s [Schlegel84]_ parametrization or simply a unit matrix
   is used.

   If the option **PARTIAL** is used the initial approximate Hessian for a
   geometry optimization is constructed from a block matrix formed of the
   parametrized Hessian and the partial Hessian (of the reaction core). If the
   reaction core spans the entire system, its Hessian is simply copied. The
   keywords :kw:`RESTART` PHESS are required.

.. keyword:: IMPLICIT NEWTON RAPHSON
   :section: &CPMD
   :options: {PREC, CONTINUE, VERBOSE, ALTERNATIVE, STEP}
   :alternatives: [N=nreg]

   Not documented.

.. keyword:: INITIALIZE WAVEFUNCTION
   :section: &CPMD
   :alternatives: [RANDOM, ATOMS]
   :notes: [PRIMITIVE]

   The initial guess for wavefunction optimization are either random functions or
   functions derived from the atomic pseudo-wavefunctions.

   For INITIALIZE WAVEFUNCTION ATOMS PRIMITIVE, CPMD will use the occupation
   information given in the &BASIS section in order to construct a minimum spin
   multiplicity (i.e. doublet or singlet) initial wavefunction from the pseudo
   atomic orbitals. This option may be helpful to avoid excessive spin
   contamination in CDFT calculations (together with an already good initial guess
   for :math:`V`) as it allows a strict initial localisation of excess spins on
   any atom species.

   **Default** is to use the **atomic pseudo-wavefunctions**.

.. keyword:: INTERFACE
   :section: &CPMD
   :options: {EGO,GMX}
   :alternatives: {[MULLIKEN, LOWDIN, ESP, HIRSHFELD],PCGFIRST}

   Use CPMD together with a classical molecular dynamics code. CPMD and the
   classical MD code are run simultaneously and communicate via a file based
   protocol. See the file egointer_utils.mod.F90 for more details. This needs a
   specially adapted version of the respective classical MD code. So far, there is
   an interface\ [egoqmmm]_ [gmxqmmm]_ to the MD programs
   ego\ [ego1]_ [ego2]_ and Gromacs\ [gmx3]_.

   When using the suboption PCGFIRST the code will use :kw:`PCG` MINIMIZE
   on the very first wavefunction optimization and then switch back to DIIS.

.. keyword:: INTFILE
   :section: &CPMD
   :options: [READ,WRITE,FILENAME]

   This keyword means *Interface File* and allows to select a special file name in
   the reading and writing stages. The file name (max 40 characters) must be
   supplied in the next line.

.. keyword:: ISOLATED MOLECULE
   :section: &CPMD

   Calculate the ionic temperature assuming that the system consists of an
   isolated molecule or cluster.

   Note: This keyword affects exclusively the determination of the number of
   dynamical degrees of freedom. This keyword does **not** activate the ’cluster
   option’ :kw:`SYMMETRY` 0, but it is activated if SYMMETRY 0 is used
   **unless** the keyword :kw:`QMMM` is set as well. It allows studying an
   isolated molecule or cluster within periodic boundary conditions.

.. keyword:: KOHN-SHAM ENERGIES
   :section: &CPMD
   :alternatives: [OFF,NOWAVEFUNCTION]

   Calculation of the Kohn-Sham energies and the corresponding
   orbitals [KS]_.

   The number of empty states that have to be calculated in addition to the
   occupied states is read from the next line.

   The Kohn-Sham orbitals are stored on the file **RESTART.x** except if the
   keyword **NOWAVEFUNCTION** is used. In this case, the program does not allocate
   memory for wavefunctions for all k points. It computes eigenvalues k point per
   k point losing information about wavefunctions. This keyword is used for band
   structure calculation to compute the eigenvalues for many k points.

   **Default** is not to calculate Kohn-Sham energies (**OFF**).

   **Warning:** The usage of this keyword needs special care (especially
   restarts).

.. keyword:: KSHAM
   :section: &CPMD
   :alternatives: [MATRIX,ROUT,STATE]

   | Write out the Kohn-Sham Hamiltonian Matrix in the orbital basis given in the
     RESTART file to KS_HAM. For this option to work the :kw:`RESTART`
     option and :kw:`OPTIMIZE WAVEFUNCTION` have to be
     activated. This option is useful for fragment orbital DFT (FODFT)
     calculations. Orbitals for the output of the FO-DFT matrix element can be
     given with the option **STATE**, then indics of the two orbitals are read
     from the next line. **ROUT** controls printing of involved orbitals.
   | **MATRIX** instructs CPMD to read a transformation matrix from the file
     LOWDIN_A to transform the KS-Hamiltonian to the non-orthogonal orbital basis
     .

.. keyword:: LANCZOS DIAGONALIZATION
   :section: &CPMD
   :options: {ALL}

   Use **Lanczos diagonalization** scheme.

   **Default** with **free energy functional**.

.. keyword:: LANCZOS DIAGONALIZATION
   :target: LANCZOS DIAGONALIZATION OPT
   :section: &CPMD
   :options: {OPT,RESET=n}

   Use **Lanczos diagonalization** scheme after (OPT) or periodically during
   (RESET=n) direct wavefunction optimization using :kw:`ODIIS`. The
   number n specifies the number of DIIS resets (ODIIS NO_RESET=nreset) due to
   poor progress until the wavefunction is diagonalized. This can be helpful if
   the wavefunction is converging very slowly.

.. keyword:: LANCZOS PARAMETER
   :section: &CPMD
   :options: [N=n]
   :alternatives: [ALL]

   Give four parameters for Lanczos diagonalization in the next line:

   -  Maximal number of Lanczos iterations (50 is enough),

   -  Maximal number for the Krylov sub-space (8 best value),

   -  Blocking dimension ( :math:`\leq NSTATE`, best in range 20-100) If you put a
      negative or zero number, this parameter is fixed by the program in function
      of the number of states (:math:`(n+1)/(int(n/100+1))`).

   -  Tolerance for the accuracy of wavefunctions

      (:math:`10^{-8}` otherwise :math:`10^{-12}` with Trotter approximation)

   If n is specified, read :math:`n-1` lines after the first one, containing a
   threshold density and a tolerance. See the hints section
   :ref:`hints:lanczos` for more information.

.. keyword:: LANGEVIN
   :section: &CPMD
   :options: {WHITE, CPMD, OPTIMAL, SMART, CUSTOM, CENTROIDOFF}
   :alternatives: [MOVECM]
   :notes: [W0, NS]

   Use a (generalized) Langevin equation to thermostat the
   simulation\ [Ceriotti10]_. By default, the component of the noise
   parallel to the center of mass velocity is removed at each step of the
   thermostat. Removal can be disabled by the option *MOVECM*.

   *CUSTOM:*

   The **number of additional momenta** of the generalized Langevin equation *NS*
   is read from the next line. The drift matrix (dimension
   :math:`(NS+1)\times(NS+1)`) is read from the file ``GLE-A``, which must be in
   the same directory in which the program is run. Optionally, the static
   covariance for the GLE dynamics can be provided in the file ``GLE-C``, so as to
   generate **non-canonical sampling**. A library of GLE parameters can be
   downloaded from http://gle4md.berlios.de/

   |
   | A few **presets** are provided, and are activated by the keywords:

   *WHITE:*

   A simple **white-noise** Langevin dynamics is used. The optimally-sampled
   frequency *W0* (in cm\ :math:`^{-1}`) is read from the next line. Note that use
   of *LANGEVIN WHITE* in conjunction with *MOLECULAR DYNAMICS CPMD* will most
   likely cause a large drift of the electronic temperature.

   *OPTIMAL:*

   An **optimal-sampling** generalized Langevin dynamics is used. The frequencies
   in the range from :math:`10^{-4}`\ *W0* up to *W0* will be sampled efficiently.
   Note that use of *LANGEVIN OPTIMAL* in conjunction with *MOLECULAR DYNAMICS
   CPMD* will cause a large drift of the electronic temperature. This option is
   suggested for use in Born-Oppenheimer MD.

   *CPMD:*

   A generalized Langevin dynamics is used which is designed to work in
   conjunction with Car-Parrinello MD. The highest ionic frequency *W0* (in
   cm\ :math:`^{-1}`) is read from the next line. Ionic frequencies down to
   :math:`10^{-4}`\ *W0* will be sampled efficiently, but not as much as for the
   *OPTIMAL* keyword.

   *SMART:*

   A generalized Langevin dynamics that aims to be as efficient as possible on the
   slowest time scale accessible to a typical ab initio simulation. In practice,
   vibrations with a time scale which is about 10000 time steps will be sampled
   optimally, and faster modes will be sampled as efficiently as possible without
   disturbing slower modes. The highest ionic frequency *W0* (in
   cm\ :math:`^{-1}`) is read from the next line. Will be about 50% more efficient
   than *OPTIMAL* for slow modes, but less efficient for fast vibrations. Use only
   with Born-Oppenheimer dynamics.

   *CENTROIDOFF:*

   For centroid and ring-polymer dynamics the generalized Langevin thermostat for
   the centroids can be switched off by using the optional keyword CENTROIDOFF
   otherwise the thermostat with the frequency read from the next line is attached
   to the centroids. For each non-centroid mode the frequency of the thermostat is
   automatically determined depending on the frequency of its associated harmonic
   potential representing the quantum kinetic energy term.

.. keyword:: LBFGS
   :section: &CPMD
   :alternatives: [NREM, NTRUST, NRESTT, TRUSTR]

   Use the limited-memory BFGS method (L-BFGS) for linear scaling **optimization**
   of the **ionic positions**. For more informations, see [LSCAL]_.
   The information about the Hessian for the quasi-Newton method employed is
   derived from the history of the optimization [LSCAL]_ [Liu89]_.

   Only one sub-option is allowed per line and the respective parameter is read
   from the next line. The parameters mean:

   *NREM*:

   **Number** of **ionic gradients** and **displacements remembered** to
   approximate the Hessian. The default is either 40 or the number of ionic
   degrees of freedom, whichever is smaller. Values greater the number of degrees
   of freedom are not advisable.

   *NTRUST:*

   *NTRUST=1* switches from a trust radius algorithm to a **line search**
   algorithm. The default value of 0 (**trust radius**) is recommended.

   *NRESTT:*

   *NRESTT\ :math:`>`\ 0* demands a **periodic reset** of the optimizer every
   *NRESTT* steps. Default is 0 (no periodic reset). This option makes only sense
   if the ionic gradient is not accurate.

   *TRUSTR:*

   Maximum and initial **trust radius**. Default is 0.5 atomic units.

   It can be useful to combine these keywords with the keywords
   :kw:`PRFO`, :kw:`CONVERGENCE` ADAPT,
   :kw:`RESTART` LSSTAT, :kw:`PRINT` LSCAL ON and others.

.. keyword:: LINEAR RESPONSE
   :section: &CPMD

   A perturbation theory calculation is done, according to the (required) further
   input in the &RESP section. In the latter, one of the possible perturbation
   types (PHONONS, LANCZOS, RAMAN, FUKUI, KPERT, NMR, EPR, see section
   :ref:`sec:resp-section`) can be chosen, accompanied by further options.

.. keyword:: LOCAL SPIN DENSITY
   :section: &CPMD

   Alias for :kw:`LSD`.

.. keyword:: LSD
   :section: &CPMD

   Use the local spin density approximation.

   **Warning:** Not all functionals in **OLDCODE** or **NEWCODE** are implemented
   for this option. Only in the **XC_DRIVER**, all functionals are available both
   for spin-restricted and -unrestricted calculations.

.. keyword:: MAXITER
   :section: &CPMD

   The maximum number of iteration steps for the self-consistency of
   wavefunctions. Recommended use instead of :kw:`MAXSTEP` for pure
   wavefunction optimisation. The value is read from the next line.

   **Default** is **10000** steps.

.. keyword:: MAXRUNTIME
   :section: &CPMD

   The maximum RUN TIME (ELAPSED TIME) in seconds to be used is read from the next
   line. The calculation will stop after the given amount of time.

   **Default** is no limit.

.. keyword:: MAXSTEP
   :section: &CPMD

   The maximum number of steps for geometry optimization or molecular dynamics to
   be performed. In the case of pure wavefunction optimisation, this keyword may
   be used instead of :kw:`MAXITER`. The value is read from the next
   line.

   **Default** is **10000** steps.

.. keyword:: MEMORY
   :section: &CPMD
   :options: {SMALL, BIG}

   Using **BIG**, the structure factors for the density cutoff are only calculated
   once and stored for reuse.

   This option allows for considerable time savings in connection with Vanderbilt
   pseudopotentials.

   **Default** is (**SMALL**) to **recalculate** them whenever needed.

.. keyword:: MIRROR
   :section: &CPMD

   Write the input file to the output.

.. keyword:: MIXDIIS
   :section: &CPMD

   Not documented

.. keyword:: MIXSD
   :section: &CPMD

   Not documented

.. keyword:: MODIFIED GOEDECKER
   :section: &CPMD
   :alternatives: [PARAMETERS]

   To be used in combination with :kw:`LOW SPIN EXCITATION` **ROKS**.

   Calculation of the off-diagonal Kohn-Sham matrix elements :math:`F_{AB}` and
   :math:`F_{BA}` (with A, B: ROKS-SOMOs) is performed according to a modified
   Goedecker-Umrigar scheme (
   :math:`F_{AB} := (1-\lambda _{AB})F_{AB} + \lambda _{AB} F_{BA}` and
   :math:`F_{BA} := (1-\lambda _{BA})F_{BA} + \lambda _{BA} F_{AB}` ). Default
   values are :math:`\lambda _{AB}=-0.5` and :math:`\lambda _{BA}=0.5`. see
   Ref. [GrimmJCP2003]_.

   With the optional keyword **PARAMETERS**: :math:`\lambda _{AB}` and
   :math:`\lambda _{BA}` are read from the next line. Can be used to avoid
   unphysical rotation of the SOMOs. Always check the orbitals!

   See also :ref:`hints:roks`.

.. keyword:: MOLECULAR DYNAMICS
   :section: &CPMD
   :alternatives: [CP, BO, ET, PT, CLASSICAL, FILE [XYZ, NSKIP=N, NSAMPLE=M]]

   Perform a molecular dynamics (MD) run. **CP** stands for a Car-Parrinello type
   MD. With the option **BO** a Born-Oppenheimer MD is performed where the
   wavefunction is re-converged after each MD-step. **EH** specifies Ehrenfest
   type dynamics according to which the Kohn-Sham orbitals are propagated in time
   (real electronic dynamics coupled to the nuclear dynamics). In this case the
   time step has to be decreased accordingly due to the small mass of the
   electrons (typical values between 0.01 and 0.1 au). If you use EH dynamics an
   additional input section &PTDDFT has to be specified. You need to start the
   dynamics with well converged KS orbitals from the RESTART file (before starting
   the EH dynamics do an optimization of the wavefunction with a convergence of
   1.D-8 or 1.D-9, if possible. An additional file called "wavefunctions" is
   produced, which contains the complex KS orbitals needed for the restart of the
   EH dynamics (see restart options in &PTDDFT). Typical (minimal) input &CPMD and
   &PTDDFT sections to be used with EH dynmiacs

   &CPMD

   MOLECULAR DYNAMICS EH

   RESTART WAVEFUNCTION COORDINATES LATEST

   CAYLEY

   RUNGE-KUTTA

   TIMESTEP

   0.01

   MAXSTEP

   10000

   &END

   &PTDDFT

   ACCURACY

   1.0D-8

   RESTART

   2

   &END

   The keywords CAYLEY and RUNGE-KUTTA specifies the algorithms used for the
   propagation of the KS orbitals (are the default and recommended options).

   **CLASSICAL** means that a MD that includes classical atoms is performed.

   **BD** stands for Bohminan Dynamics. In the next line the number of
   trajectories (:math:`N_{traj}`) used to describe the nuclear wavepacket are
   specified. The initial spread of the Gaussians associated to each element is
   guessed from the De Broglie wavelength. A special tuning of this value for the
   different elements may be necessary (not yet automatized). During the dynamics
   the :math:`N_{traj}` trajectories (each corresponding to a configuration space
   point or fluid element) are interacting within each other through the action of
   the quantum potential. For the more details about the implementation see
   Ref. [ABDY]_.

   If **FILE** is set, then the trajectory is reread from a file instead of being
   calculated. This is useful for performing analysis on a previous trajectory.
   Can be used in conjunction with the standard MD options like DIPOLE DYNAMICS
   and WANNIER; some other features like LINEAR RESPONSE are also enabled. The
   trajectory is read from a file named TRAJSAVED (usually a copy of a previous
   TRAJECTORY file), or TRAJSAVED.xyz if **XYZ** is set. **NSKIP** and **NSAMPLE**
   control the selection of frames read: the frame read at step ISTEP is
   NSKIP+ISTEP*NSAMPLE.

   **Default** is **CP**.

.. keyword:: MOVERHO
   :section: &CPMD

   Mixing used during optimization of geometry or molecular dynamics. Use atomic
   or pseudowavefunctions to project wavefunctions in order to calculate the new
   ones with movement of atoms. Read in the next line the parameter (typically
   0.2).

.. keyword:: MOVIE
   :section: &CPMD
   :alternatives: [OFF, SAMPLE]

   Write the atomic coordinates without applying periodic boundary conditions in
   MOVIE format every *IMOVIE* time steps on file *MOVIE*. *IMOVIE* is read from
   the next line.

   **Default** is **not** to write a movie file.

.. keyword:: NOGEOCHECK
   :section: &CPMD

   Default is to check all atomic distances and stop the program if the smallest
   disctance is below 0.5 Bohr. This keyword requests not to perform the check.

.. keyword:: NONORTHOGONAL ORBITALS
   :section: &CPMD
   :alternatives: [OFF]

   Use the norm constraint method [HutterIP]_ for molecular dynamics
   or nonorthogonal orbitals in an optimization run.

   On the next line the limit of the off diagonal elements of the overlap matrix
   is defined. **Warning:** Adding or deleting this option during a MD run needs
   special care.

.. keyword:: NOSE
   :section: &CPMD
   :options: {IONS, ELECTRONS, CELL}
   :alternatives: [ULTRA,MASSIVE,CAFES,LOCAL,CENTROIDOFF]
   :notes: [T0]

   **Nosé-Hoover chains** [Nose84]_ [Hoover85]_ for the **ions**,
   **electrons**, or **cell parameters** are used.

   The **target temperature** in Kelvin and the **thermostat frequency** in
   :math:`cm^{-1}`, respectively the **fictitious kinetic energy** in atomic units
   and the **thermostat frequency** in :math:`cm^{-1}` are read from the next
   line. Two files NOSE_ENERGY and NOSE_TRAJEC are written at each step containing
   the Nosé-Hoover kinetic, potential and total energies along the dynamics
   (NOSE_ENERGY) and the Nosé-Hoover variables and their velocities (NOSE_TRAJEC);
   these are useful in a wealth of post-processing calculations such as, e. g.
   heat transfer problems\ [heat1]_ [heat2]_. For the ionic case the
   additional keyword **ULTRA** selects a thermostat for each species, the keyword
   **MASSIVE** selects a thermostat for each degree of freedom, and the keyword
   **CAFES** can be used to give different temperatures to different groups of
   atoms\ [cafes02]_. The syntax in the **CAFES** case is:

   ``NOSE IONS CAFES``

       *ncafesgrp*

     *cpnumber_a_1*  *cpnumber_a_2*  Temperature Frequency

   …

     *cpnumber_n_1*  *cpnumber_n_2*  Temperature Frequency

   There are *ncafesgrp* groups, specified by giving their first CPMD atom number
   (*cpnumber_X_1*) and last CPMD atom number (*cpnumber_X_2*). In the case of
   hybrid QM/MM simulations, you have to consult the QMMM_ORDER file to find those
   numbers. The temperature and frequency can be different for each group. All
   atoms of the system have to be in a CAFES group. A new file, ``CAFES`` is
   created containing the temperature of each group (cols. 2 …\ *ncafesgrp+1*) and
   the energy of the Nose-Hoover chains of that group (last columns).

   Using CAFES with different temperatures only makes sense if the different
   groups are decoupled from each other by increasing the masses of the involved
   atoms. The mass can be specified in the topology / or with the
   :kw:`ISOTOPE` keyword. However, you can only change the mass of a
   complete CPMD species at a time. Hence, the topology and/or the input should be
   such that atoms of different CAFES group are in different species.

   **NOTE:** CAFES is currently not restartable.

   The keyword **LOCAL** collects groups of atoms to separate thermostats, each
   having its own Nosé-Hoover chain. Specify the local thermostats as follows:

   =========================================== ======================= ==========
   ``NOSE IONS LOCAL``
   :math:`n_l` *(number of local thermostats)*
   *temperature 1*                             *frequency 1*
   ⋮
   *temperature :math:`n_l`*                   *frequency :math:`n_l`*
   :math:`n_r` *(number of atom ranges)*
   *thermostat number*                         *start atom*            *end atom*
   ⋮                                           *(:math:`n_r` entries)*
   =========================================== ======================= ==========

   The parser for the atom ranges uses either the CPMD ordering or the GROMOS
   ordering in case of classical or QM/MM runs. Multiple ranges may be specified
   for the same thermostat. Atoms belonging to the same CPMD constraint or the
   same solvent molecule in QM/MM runs must belong to the same local thermostat.

   If **T0** option is present, the initial temperature for the Nosé-Hoover chains
   are read soon after the thermostat frequencies in the same line (also for the
   LOCAL thermostat). By default it is same as the target temperature of the
   thermostat. Note: This is not implemented for the CAFES thermostat.

   The keyword **CENTROIDOFF** is intended to be used in connection with the
   centroid and ring-polymer dynamics. In this type of path integral MD, you can
   request not to attach the Nosé-Hoover chains for ions to the centroids by
   explicit use of the keyword CENTROIDOFF, otherwise the thermostat with the
   frequency read from the next line is attached to the centroids. For each
   non-centroid mode the frequency with the Nosé-Hoover chains is automatically
   determined depending on the frequency of its associated harmonic potential
   representing the *quantum* kinetic energy term.

.. keyword:: NOSE PARAMETERS
   :section: &CPMD

   The **parameters** controlling the **Nosé
   thermostats** [Nose84]_ [Hoover85]_ are read in the following order
   from the next line:

   The **length** of the Nosé-Hoover chain for the **ions**,

   the **length** of the Nosé-Hoover chain for the **electrons**,

   the **length** of the Nosé-Hoover chain for the **cell parameters**.

   (The respective **default** values are **4**.)

   The **multiplication factor** (NEDOF0, a real number) for the number of
   **electronic** degrees of freedom. The used degrees of freedom (NEDOF) are
   defined as :math:`NEDOF=NEDOF0*X` If NEDOF0 is a negative number X is the true
   number of DOFs, if it’s a positive number, X is the number of electronic states
   (**default** for NEDOF0 is **6**).

   The order of the **Suzuki/Yoshida integrator** (**default** is **7**, choices
   are 3, 5, 7, 9, 15, 25, 125 and 625),

   and the **decomposition ratio** of the time step (**default** is **1**).

   If this keyword is omitted, the defaults are used.

   **If the keyword is used all parameters have to be specified.**

.. keyword:: ODIIS
   :section: &CPMD
   :alternatives: [NOPRECONDITIONING,NO_RESET=nreset]

   Use the method of **direct inversion** in the iterative subspace for
   **optimization** of the **wavefunction** [Hutter94a]_.

   The number of DIIS vectors is read from the next line.

   (ODIIS with **10 vectors** is the **default** method in optimization runs.)

   The preconditioning is controlled by the keyword :kw:`HAMILTONIAN CUTOFF`. Optionally preconditioning can be disabled.

   By default, the number of wavefunction optimization cycles until DIIS is
   **reset** on poor progress, is the number of DIIS vectors. With **ODIIS
   NO_RESET**, this number can be changed, or DIIS resets can be **disabled**
   altogether with a value of -1.

.. keyword:: OPTIMIZE GEOMETRY
   :section: &CPMD
   :options: [XYZ, SAMPLE]

   This option causes the program to optimize the geometry of the system through a
   sequence of wavefunction optimizations and position updates. The additional
   keyword XYZ requests writing the “trajectory” of the geometry additionally in
   xmol/xyz-format in a file *GEO_OPT.xyz*. If the keyword SAMPLE is given,
   *NGXYZ* is read from the next line, and then only every *NGXTZ* step is written
   to the xmol/xyz file. The **default** is to write every step (*NGXYZ* =
   :math:`1`).

   By default the a BFGS/DIIS algorithm is used (see :kw:`GDIIS`) to
   updated the ionic positions. Other options are: :kw:`LBFGS`,
   :kw:`PRFO`, and :kw:`STEEPEST DESCENT` IONS. See
   :kw:`OPTIMIZE WAVEFUNCTION` for details on the
   corresponding options for wavefunction optimizations.

.. keyword:: OPTIMIZE WAVEFUNCTION
   :section: &CPMD

   Request a single point energy calculation through a wavefunction optimization.
   The resulting total energy is printed (for more output options see, e.g.,:
   :kw:`PRINT`, :kw:`RHOOUT`, :kw:`ELF`) and a
   :kw:`RESTART` file is written. This restart file is a prerequisite
   for many other subsequent calculation types in CPMD, e.g. :kw:`MOLECULAR DYNAMICS` CP or :kw:`PROPERTIES`. By
   default a DIIS optimizer is used (see :kw:`ODIIS`), but other options
   are: :kw:`PCG` (optionally with MINIMIZE), :kw:`LANCZOS DIAGONALIZATION`, :kw:`DAVIDSON DIAGONALIZATION`, and :kw:`STEEPEST DESCENT` ELECTRONS.

.. keyword:: ORBITAL HARDNESS
   :section: &CPMD
   :alternatives: [LR,FD]

   Perform an orbital hardness calculation. See section &Hardness for further
   input options.

.. keyword:: ORTHOGONALIZATION
   :section: &CPMD
   :options: {LOWDIN, GRAM-SCHMIDT}
   :alternatives: [MATRIX]

   Orthogonalization in optimization runs is done either by a Löwdin (symmetric)
   or Gram-Schmidt procedure.

   **Default** is Gram-Schmidt except for parallel runs where Löwdin
   orthogonalization is used with the conjugate-gradient scheme.

   With the additional keyword **MATRIX** the Löwdin transformation matrix is
   written to a file named LOWDIN_A.

.. keyword:: PARA_BUFF_SIZE
   :target: PARA BUFF SIZE
   :section: &CPMD

   Set the buffer size for parallel operation (sum, ...). Default is 2**16. The
   size is read from the next line.

.. keyword:: PARA_STACK_BUFF_SIZE
   :target: PARA STACK BUFF SIZE
   :section: &CPMD

   Set the stack buffer size for parallel operation (sum, ...). Default is 2**8.
   The size is read from the next line.

.. keyword:: PARA_USE_MPI_IN_PLACE
   :target: PARA USE MPI IN PLACE
   :section: &CPMD

   Use MPI_IN_PLACE for parallel operation (sum, ...). Default is FALSE.

.. keyword:: PARRINELLO-RAHMAN
   :section: &CPMD
   :options: {NPT,SHOCK}

   To be used together with :kw:`MOLECULAR DYNAMICS`.

   A **variable cell MD** with the **Parrinello-Rahman Lagrangian** is
   performed [parrah]_ [parrah2]_. With the additional keyword a
   **constant NPT MD** using the method of Martyna, Tobias, and
   Klein [Martyna94]_.

   If this keyword is used together with other run options like OPTIMIZE
   WAVEFUNCTIONS, calculations with different reference cells can be performed.

   With the additional keywork **SHOCK** a MD simulation using the multiscale
   shock method [shock]_ is performed.

.. keyword:: PATH INTEGRAL
   :section: &CPMD

   Perform a **path integral molecular dynamics**
   calculation [Marx94]_ [Marx96]_.

   This keyword requires further input in the section &PIMD ... &END.

.. keyword:: PATH MINIMIZATION
   :section: &CPMD

   Perform a **mean free energy path** search [Eijnden06]_.

   This keyword requires further input in the section &PATH ... &END.

.. keyword:: PATH SAMPLING
   :section: &CPMD

   Use CPMD together with a reaction path sampling [tps]_ program.
   This needs special software. Note: this keyword has *nothing* to do with path
   integral MD as activated by the keyword PATH INTEGRAL and as specified in the
   section &PIMD ... &END.

.. keyword:: PCG
   :section: &CPMD
   :alternatives: [MINIMIZE,NOPRECONDITIONING]

   Use the method of **preconditioned conjugate gradients** for **optimization**
   of the **wavefunction**.

   The fixed step length is controlled by the keywords :kw:`TIMESTEP ELECTRONS` and :kw:`EMASS`.

   If the additional option **MINIMIZE** is chosen, then additionally line
   searches are performed to improve the preconditioning.

   The preconditioning is controlled by the keyword :kw:`HAMILTONIAN CUTOFF`. Optionally preconditioning can be disabled.

.. keyword:: PRFO
   :section: &CPMD
   :alternatives: [MODE, MDLOCK, TRUSTP, OMIN, PRJHES, DISPLACEMENT, HESSTYPE]

   Use the partitioned rational function optimizer (P-RFO) with a quasi-Newton
   method for **optimization** of the **ionic positions**. For more informations,
   see [LSCAL]_. The approximated Hessian is updated using the Powell
   method [Powell71]_. This method is used to find **transition
   states** by **following eigenmodes** of the approximated
   Hessian [Banerjee85]_ [LSCAL]_.

   Only one suboption is allowed per line and the respective parameter is read
   from the next line. The suboption **PRJHES** does not take any parameter. If it
   is present, the translational and rotational modes are removed from the
   Hessian. This is only meaningful for conventional (not microiterative)
   transition state search. The parameters mean:

   *MODE*:

   Number of the initial Hessian **eigenmode** to be followed. Default is 1
   (lowest eigenvalue).

   *MDLOCK:*

   *MDLOCK=1* switches from a mode following algorithm to a **fixed eigenvector**
   to be maximized. The default value of 0 (**mode following**) is recommended.

   *TRUSTP:*

   Maximum and initial **trust radius**. Default is 0.2 atomic units.

   *OMIN:*

   This parameter is the minimum **overlap** between the maximized mode of the
   previous step and the most overlapping eigenvector of the current Hessian. The
   trust radius is reduced until this requirement is fulfilled. The default is
   0.5.

   *DISPLACEMENT:*

   Finite-difference **displacement** for initial partial Hessian. The default is
   0.02.

   *HESSTYPE:*

   **Type** of initial partial Hessian. 0: Finite-difference. 1: Taken from the
   full Hessian assuming a block-diagonal form. See keyword
   :kw:`HESSIAN`. The default is 0.

   It can be useful to combine these keywords with the keywords
   :kw:`CONVERGENCE` ENERGY, :kw:`RESTART` LSSTAT,
   :kw:`RESTART` PHESS, :kw:`PRFO` NSVIB, :kw:`PRINT`
   LSCAL ON and others.

.. keyword:: PRFO
   :target: PRFO NVAR
   :section: &CPMD
   :alternatives: [NVAR, CORE, TOLENV, NSMAXP]

   If any of these suboptions is present, the **microiterative transition state
   search** scheme for **optimization** of the **ionic positions** is used. For
   more informations, see [LSCAL]_. A combination of the **L-BFGS**
   and **P-RFO** methods is employed for linear scaling search for transition
   states [LSCAL]_ [Turner99]_. Before each P-RFO step in the reaction
   core towards the transition state, the **environment** is fully **relaxed**
   using L-BFGS.

   Only one suboption is allowed per line. The **reaction core** can be selected
   using the **NVAR** or **CORE=ncore** suboptions. The value in the line after
   **PRFO NVAR** sets the number of ionic **degrees of freedom** in the reaction
   core. The *ncore* values following the line **PRFO CORE=ncore** select the
   **member atoms** of the reaction core. If unspecified, the *NVAR/3* first atoms
   form the reaction core. The parameters read with the two remaining suboptions
   are:

   *TOLENV*:

   **Convergence criterion** for the maximum component of the gradient acting on
   the ions of the **environment** until a P-RFO step within the reaction core is
   performed. Default is one third of the convergence criterion for the gradient
   of the ions (**CONVERGENCE GEOMETRY**).

   *NSMAXP:*

   Maximum number of P-RFO **steps** to be performed in the reaction core. The
   keyword **HESSCORE** corresponds to **PRFO NSMAXP** with *NSMAXP=0*.

   It can be useful to combine these keywords with the keywords
   :kw:`LBFGS`, :kw:`CONVERGENCE` ADAPT,
   :kw:`CONVERGENCE` ENERGY, :kw:`RESTART` LSSTAT,
   :kw:`RESTART` PHESS, :kw:`PRFO NSVIB`,
   :kw:`PRINT` LSCAL ON, the other suboptions of PRFO, and others.

.. keyword:: PRFO NSVIB
   :section: &CPMD

   Perform a **vibrational analysis** every NSVIB P-RFO steps **on the fly**. This
   option only works with the P-RFO and microiterative transition state search
   algorithms. In case of microiterative TS search, only the reaction core is
   analyzed.

.. keyword:: PRINT
   :section: &CPMD
   :options: {ON,OFF}
   :alternatives: [INFO, EIGENVALUES, COORDINATES, LSCAL, FORCES, WANNIER]

   A **detailed output** is printed every *IPRINT* iterations. Either only
   different contribution to the energy or in addition the atomic coordinates and
   the forces are printed. *IPRINT* is read from the next line if the keywords
   **ON** or **OFF** are not specified.

   **Default** is **only energies** after the first step and at the end of the
   run. OFF switches the output off.

.. keyword:: PRINT ENERGY
   :section: &CPMD
   :options: {ON, OFF}
   :alternatives: [EKIN, ELECTROSTATIC, ESR, ESELF, EFREE, EBAND, ENTROPY, EPSEU, EHEP, EHEE, EHII, ENL, EXC, VXC, EGC, EBOGO]

   Display or not information about energies.

.. keyword:: PRNGSEED
   :section: &CPMD

   The seed for the random number generator is read as an integer number from the
   next line.

.. keyword:: PROJECT
   :section: &CPMD
   :alternatives: [NONE, DIAGONAL, FULL]

   This keyword is controlling the calculation of the constraint force in
   optimization runs.

.. keyword:: PROPAGATION SPECTRA
   :section: &CPMD

   Calculates the electronic absorption spectra using the TDDFT propagation of the
   Kohn-Sham orbitals. Use the section &PTDDFT to define the parameters. Use this
   principal keyword always with CAYLEY (in &CPMD). The program produces a file
   "dipole.dat" with the time series of the variation of the dipole in x, y, and z
   directions. After Fourier transform of this file one gets the desired
   absorption spectra.

   Typical (minimal) input file (for the sections &CPMD and &PTDDFT) &CPMD

   PROPAGATION SPECTRA RESTART WAVEFUNCTION COORDINATES LATEST CAYLEY &END

   &PTDDFT

   ACCURACY

   1.0D-8

   N_CYCLES

   100000

   PROP_TSTEP

   0.01

   EXT_PULSE

   1.D-5

   PERT_DIRECTION

   1

   RESTART

   2

   &END

   The time step is specified by setting :kw:`PROP-TSTEP`. The total
   number of iteration is controlled by :kw:`N-CYCLES`.

.. keyword:: PROPERTIES
   :section: &CPMD

   Calculate some properties.

   This keyword requires further input in the section &PROP …&END.

.. keyword:: QMMM
   :section: &CPMD
   :alternatives: [QMMMEASY]

   Activate the hybrid QM/MM code. This keyword requires further input in the
   section &QMMM …&END.

   The QM driver is the standard CPMD. An interface program (**MM_Interface**) and
   a classic force field
   (Gromos\ [gromos96]_/Amber\ [amber7]_-like) are
   needed to run the code in hybrid
   mode\ [qmmm02]_ [qmmm03]_ [qmmm04]_ [qmmm05]_ [qmmm06]_. This code requires a
   *special licence* and is **not** included in the standard CPMD code. (see
   section :ref:`sec:qmmm` for more information on the available options
   and the input format).

.. keyword:: QUENCH
   :section: &CPMD
   :alternatives: [IONS, ELECTRONS, CELL, BO]

   The **velocities** of the **ions**, **wavefunctions** or the **cell** are set
   to zero at the beginning of a run.

   With the option **BO** the wavefunctions are converged at the beginning of the
   MD run.

.. keyword:: RANDOMIZE
   :section: &CPMD
   :alternatives: [COORDINATES, WAVEFUNCTION, DENSITY, CELL]

   The **ionic positions** or the **wavefunction** or the **cell parameters** are
   **randomly displaced** at the beginning of a run.

   The maximal amplitude of the displacement is read from the next line.

.. keyword:: RATTLE
   :section: &CPMD

   This option can be used to set the maximum number of iterations and the
   tolerance for the **iterative orthogonalization**. These two numbers are read
   from the next line.

   **Defaults** are 30 and :math:`10^{-6}`.

.. keyword:: REAL SPACE WFN KEEP
   :section: &CPMD
   :alternatives: [SIZE]

   The real space wavefunctions are kept in memory for later reuse. This minimizes
   the number of Fourier transforms and can result in a significant speedup at the
   expense of a larger memory use. With the option **SIZE** the maximum available
   memory for the storage of wavefunctions is read from the next line (in MBytes).
   The program stores as many wavefunctions as possible within the given memory
   allocation.

.. keyword:: RESCALE OLD VELOCITIES
   :section: &CPMD

   Rescale **ionic** velocities after :kw:`RESTART` to the temperature
   specified by either :kw:`TEMPERATURE`,
   :kw:`TEMPCONTROL` **IONS**, or :kw:`NOSE` **IONS**.
   Useful if the type of ionic thermostatting is changed, (do not use RESTART
   NOSEP in this case).

   Note only for path integral runs: the scaling is only applied to the first
   (centroid) replica.

.. keyword:: RESTART
   :section: &CPMD
   :notes: [{OPTIONS}]

   This keyword controls what data is read (at the beginning) from the file
   RESTART.x.

   **Warning:** You can only read data that has been previously written into the
   RESTART-file.

   A list of different *OPTIONS* can be specified. List of valid options:

   WAVEFUNCTION
      Read old **wavefunction** from restart file.

   OCCUPATION
      Read old **occupation numbers** (useful for free energy functional.

   COORDINATES
      Read old **coordinates** from restart file.

   VELOCITIES
      Read old **ionic, wavefunction and (cell) velocities** from restart file.

   CELL
      Read old **cell parameters** from restart file.

   GEOFILE
      Read old **ionic positions and velocities** from file **GEOMETRY**. This
      file is updated every time step. It has higher priority than the COORDINATES
      option.

   ACCUMULATORS
      Read old **accumulator values**, for example the time step number, from
      restart file.

   HESSIAN
      Read old **approximate Hessian** from file *HESSIAN*.

   NOSEE
      Restart **Nosé thermostats** for **electrons** with values stored on restart
      file.

   NOSEP
      | Restart **Nosé thermostats** for **ions** with values stored on
      | restart file.

   NOSEC
      Restart **Nosé thermostats** for **cell** parameters with values stored on
      restart file.

   LATEST
      Restart from the **latest restart** file as indicated in file *LATEST*.

   PHESS
      Read partial Hessian (Hessian of the reaction core) for **transition state
      search** or **vibrational analysis** from restart file. Useful with the
      keywords :kw:`PRFO` or :kw:`HESSIAN`
      [DISCO,SCHLEGEL,UNIT] PARTIAL.

   LSSTAT
      Read all **status information** of the **linear scaling optimizers** (L-BFGS
      and P-RFO) including L-BFGS history but excluding partial Hessian for P-RFO
      from restart file. The **partial Hessian** is read separately using
      **RESTART PHESS**. Useful with the keywords :kw:`LBFGS` and/or
      :kw:`PRFO`.

   ADPTTL
      Read **wavefunction convergence criteria** at the current point of geometry
      optimization from restart file. Useful with the keywords
      :kw:`CONVERGENCE` [ADAPT, ENERGY, CALFOR].

   VIBANALYSIS
      Use the information on finite differences stored in the file *FINDIF*. This
      option requires a valid restart file for the wavefunctions, even when
      wavefunctions and coordinates are recalculated or read from the input file.

   POTENTIAL
      Read an old potential from the restart file. This applies to restarts for
      Kohn-Sham energy calculations.

   KPOINTS
      Restart with k points.

   DENSITY
      Restart with electronic density.

   GLE
      Restart the extended variables of the GLE dynamics

   PRNG
      Restart the internal state of the Marsaglia random number generator

   CONSTRAINTS
      Restart with old values for constraints. This option is mainly for
      restraints with GROWTH option.

   EXTRAP
      Restart from a previously saved wavefunction history. See :kw:`EXTRAPOLATE WFN` for details.

   CVDFT
      Read CDFT multiplier :math:`V` from ``CDFT_RESTART``, for HDA run read
      :math:`V` and wavefunction of the appropriate state.

   ALL
      Restart with all fields of RESTART file

.. keyword:: RESTFILE
   :section: &CPMD

   The number of distinct :kw:`RESTART` files generated during CPMD
   runs is read from the next line.

   The restart files are written in turn. **Default is 1**. If you specify e.g. 3,
   then the files RESTART.1, RESTART.2, RESTART.3 are used in rotation.

.. keyword:: REVERSE VELOCITIES
   :section: &CPMD

   Reverse the ionic and electronic (if applicable) velocities after the initial
   setup of an MD run. This way one can, e.g., go “backwards” from a given
   :kw:`RESTART` to improve sampling of a given MD “path”.

.. keyword:: RFO ORDER=nsorder
   :section: &CPMD

   Rational function approximation combined with a quasi-Newton method (using
   BFGS) for **optimization** of the **ionic positions** is
   used [Banerjee85]_. A saddle point of order nsorder is searched
   for.

.. keyword:: RHOOUT
   :section: &CPMD
   :alternatives: [BANDS,SAMPLE=nrhoout]

   **Store** the **density** at the end of the run on file *DENSITY*.

   If the keyword BANDS is defined then on the following lines the number of bands
   (or orbitals) to be plotted and their index (starting from 1) have to be given.
   If the position specification is a negative number, then the wavefunction
   instead of the density is written. Each band is stored on its own file
   *DENSITY.num*. For spin polarized calculations besides the total density also
   the spin density is stored on the file *SPINDEN*. The following example will
   request output of the orbitals or bands number 5, 7, and 8 as wavefunctions:

   ::

                  RHOOUT BANDS
                    3
                    -5 -7 -8

   With the optional keyword **SAMPLE** the requested file(s) will be written
   every *nrhoout* steps during an MD trajectory. The corresponding time step
   number will be appended to the filename.

.. keyword:: ROKS
   :section: &CPMD
   :alternatives: {SINGLET, TRIPLET},{DELOCALIZED, LOCALIZED, GOEDECKER}

   Calculates the first excited state using Restricted Open-shell Kohn-Sham
   theory [Frank98]_. By default, the singlet state is calculated
   using the delocalized variant of the modified Goedecker-Umrigar scheme, which
   is supposed to work in most cases. That is, for doing a ROKS simulation, it is
   usually sufficient to just include this keyword in the CPMD section (instead of
   using the :kw:`LSE <LOW SPIN EXCITATION>` input). See :ref:`hints:roks`
   for further information.

.. keyword:: RUNGE_KUTTA
   :target: RUNGE-KUTTA
   :section: &CPMD

   Defines the integration schemes used in the Ehrenfest :kw:`MOLECULAR DYNAMICS`. Always used this option.

.. keyword:: SCALED MASSES
   :section: &CPMD
   :alternatives: [OFF]

   Switches the usage of g-vector dependent masses on/off.

   The number of shells included in the analytic integration is controlled with
   the keyword **HAMILTONIAN CUTOFF**.

   By **default** this option is switched **off**.

.. keyword:: SHIFT POTENTIAL
   :section: &CPMD

   After this keyword, useful in Hamiltonian diagonalization, the shift value
   :math:`V_{\rm shift}` must be provided in the next line. This option is used in
   the Davidson diagonalization subroutine and shifts rigidly the total electronic
   potential as
   :math:`V_{\rm pot}({\bf r}) \to V_{\rm pot}({\bf r})+V_{\rm shift}` then it is
   subtracted again at the end of the main loop, restoring back the original
   :math:`V_{\rm pot}({\bf r})` that remains basically unaffected once that the
   calculation is completed.

.. keyword:: SOC
   :section: &CPMD

   Compute the spin-orbit coupling (SOC) elements between the two states specified
   in the next line. To be used together with the principal keyword SPECTRA. The
   calculation is based on the Auxiliary Many-Electron Wavefunction Approach for
   TDDFT developed in Ref. [taver]_. For more information see
   Ref. [SOC]_.

.. keyword:: SPLINE
   :section: &CPMD
   :alternatives: [POINTS, QFUNCTION, INIT, RANGE]

   This option controls the generation of the pseudopotential functions in
   g-space.

   All pseudopotential functions are first initialized on a evenly spaced grid in
   g-space and then calculated at the needed positions with a spline
   interpolation. The number of spline points is read from the next line when
   **POINTS** is specified.

   ( The **default** number is **5000**.) For calculations with the small cutoffs
   typically used together with Vanderbilt PP a much smaller value, like 1500 or
   2000, is sufficient.

   In addition it is possible to keep the Q-functions of the Vanderbilt
   pseudopotentials on the spline grid during the whole calculation and do the
   interpolation whenever needed. This option may be useful to save time during
   the initialization phase and memory in the case of Vanderbilt pseudopotentials
   when the number of shells is not much smaller than the total number of plane
   waves, i.e. for all cell symmetries except simple cubic and fcc.

.. keyword:: SSIC
   :section: &CPMD

   Apply an *ad hoc* Self Interaction Correction (SIC) to the ordinary DFT
   calculation expressed in terms of total energy as

   .. math:: E^{\rm tot}-a\cdot  E_H[m]- b\cdot E_{xc}[m, 0]

   where :math:`m({\bf x}) = \rho_\alpha({\bf x})-\rho_\beta({\bf x})`. The value
   of :math:`a` must be supplied in the next line, while in the present
   implementation :math:`b` is not required, being the optimal values
   :math:`a=0.2` and :math:`b=0.0` according to Ref. [SSIC]_. These
   are assumed as default values although it is not always the case
   [dna_sic]_. Note that if you select negative :math:`\{a, b \}`
   parameters, the signs in the equation above will be reversed. The Hartree
   electronic potential is changed accordingly as
   :math:`V_H[\rho] \to V_H[\rho] \pm a\cdot V_{\rm SIC}[m]`, being

   .. math:: V_{\rm SIC}[m]=\frac{\delta E_H[m]}{\delta m({\bf x})}

   where the sign is :math:`+` for :math:`\alpha` spin and :math:`-` for
   :math:`\beta` spin components, respectively. Be aware that this keyword should
   be used together with :math:`LSD` (set by default).

.. keyword:: STEEPEST DESCENT
   :section: &CPMD
   :alternatives: [ELECTRONS, IONS, CELL],[NOPRECONDITION\-ING],[LINE]

   NOPRECONDITIONING works only for electrons and LINE only for ions. Use the
   method of **steepest descent** for the **optimization** of wavefunction and/or
   atomic positions and/or cell.

   If both options are specified in a geometry optimization run, a simultaneous
   optimization is performed.

   Preconditioning of electron masses (scaled masses) is used by default. The
   preconditioning is controlled by the keyword **HAMILTONIAN CUTOFF**. Optionally
   preconditioning can be disabled.

   For ions optimization, the step length is controlled by the keywords
   **TIMESTEP** and **EMASS**.

.. keyword:: STORE
   :section: &CPMD
   :options: {OFF}
   :alternatives: [WAVEFUNCTIONS, DENSITY, POTENTIAL]

   The :kw:`RESTART` file is **updated** every *ISTORE* steps. *ISTORE*
   is read from the next line. **Default** is at the **end of the run**.

   Moreover, in the same line of the number ISTORE, you can specify the number of
   self-consistent iterations (with SC=number) between two updates of restart
   file. If OFF is specified , do not store wavefunctions and/or density (*ISTORE*
   is not necessary).

.. keyword:: STRESS TENSOR
   :target: STRESS TENSOR CPMD
   :section: &CPMD

   Calculate the **stress tensor** every *NSTEP* iteration in a constant volume
   MD.

   *NSTEP* is read from the next line. Works also for wavefunction or geometry
   optimisation. In this case NSTEP is meaningless.

.. keyword:: STRUCTURE
   :section: &CPMD
   :alternatives: [BONDS, ANGLES, DIHEDRALS, SELECT]

   Print **structure information** at the end of the run.

   Bonds, angles and dihedral angles can be printed. Dihedral angles are defined
   between 0 and 180 degrees. This might change in the future.

   If the option **SELECT** is used the output is restricted to a set of atoms.
   The number of atoms and a list of the selected atoms has to be given on the
   next lines.

.. keyword:: SUBTRACT
   :section: &CPMD
   :alternatives: [COMVEL, ROTVEL]

   If COMVEL is selected, the total momentum of the system is removed, if ROTVEL
   is selected the global angular momentum of the system is removed. Both options
   can be used separately and simultaneously. The subtraction is done each
   **ncomv** or **nrotv** steps, where the value is read in the next line.

   If this key is activated but no number provided, the **default** is
   :math:`10000` steps.

   **Note**: The use of these keywords is strongly recommended for long runs (e.g.
   :math:`t>10` ps) and/or low density systems (e.g. isolated molecules, gas phase
   & Co.). Otherwise the whole system will start to translate and/or rotate toward
   a (random) direction with increasing speed and spinning. The “relative”
   translation within the system slows down correspondingly and thus the system
   effectively cools down. As a consequence dynamic properties, like
   self-diffusion coefficients will be wrong.

   This option should not be used for systems, where some atoms are kept at fixed
   positions, e.g. slab configurations. Here the center of mass may (or should)
   move. Due to the interactions with the fixed atoms, a drift of the whole system
   should be much reduced, anyways.

   **Note**: since the subtracted kinetic energy is put back into the system by
   simple rescaling of the ionic velocities, these options is not fully compatible
   with :kw:`NOSE` thermostats.

.. keyword:: SURFACE HOPPING
   :section: &CPMD

   Nonadiabatic dynamics involving the ground state and a :kw:`ROKS`
   excited state\ [surfhop]_. Do NOT use this keyword together with
   :kw:`T-SHTDDFT`, which invokes the surface hopping MD scheme based
   on TDDFT [TDDFT-SH]_ (see :kw:`T-SHTDDFT`).

.. keyword:: TDDFT
   :section: &CPMD

   Calculate the energy according to TDDFT. This keyword can be used together with
   :kw:`OPTIMIZE GEOMETRY` or :kw:`MOLECULAR DYNAMICS` BO. Use the &TDDFT section to set
   parameters for the calculation. This keyword requires :kw:`RESTART`
   LINRES.

.. keyword:: TEMPCONTROL
   :section: &CPMD
   :alternatives: [IONS, ELECTRONS, CELL]

   The **temperature** of the **ions** in Kelvin or the **fictitious kinetic
   energy** of the **electrons** in atomic units or the **kinetic energy** of the
   **cell** in atomic units (?) is controlled by scaling.

   The **target** temperature and the **tolerance** for the ions or the target
   kinetic energy and the tolerance for the electrons or the cell are read from
   the next line.

   As a gentler alternative you may want to try the :kw:`BERENDSEN`
   scheme instead.

.. keyword:: TEMPERATURE
   :section: &CPMD
   :alternatives: [RAMP]

   The **initial temperature** in Kelvin of the **system** is read from the next
   line. With the additional keyword **RAMP** the temperature can be linearly
   ramped to a target value and two more numbers are read, the ramping target
   temperature in Kelvin and the ramping speed in Kelvin per atomic time unit (to
   get the change per timestep you have to multiply it with the value of
   :kw:`TIMESTEP`). Note that this ramping affects the target
   temperatures for :kw:`TEMPCONTROL`, :kw:`BERENDSEN`
   and the global :kw:`NOSE` thermostats.

.. keyword:: TEMPERATURE ELECTRON
   :section: &CPMD

   The **electronic temperature** is read from the next line.

   **Default** is :math:`1000`\ K.

.. keyword:: TIMESTEP
   :section: &CPMD

   The time step in atomic units is read from the next line.

   **Default** is a time step of **5 a.u.** (:math:`1\, a.u. = 0.0241888428` fs).

.. keyword:: TIMESTEP ELECTRONS
   :section: &CPMD

   The time step for electron dynamics in atomic units is read from the next line.
   This is can be used to tweak the convergence behavior of the wavefunction
   optimization in Born-Oppenheimer dynamics, where the default time step may be
   too large. see, e.g. :kw:`PCG`

.. keyword:: TIMESTEP IONS
   :section: &CPMD

   The time step in atomic units is read from the next line.

.. keyword:: TRACE
   :section: &CPMD
   :alternatives: [ALL,MASTER]

   Activate the tracing of the procedures. *ALL* specifies that all the mpi tasks
   are traced. *ALL* specifies that only the master is traced.

.. keyword:: TRACE_MAX_CALLS
   :target: TRACE MAX CALLS
   :section: &CPMD

   Set the maximal number of calls for tracing. The number is read from the next
   line.

.. keyword:: TRACE_MAX_DEPTH
   :target: TRACE MAX DEPTH
   :section: &CPMD

   Set the maximal depth for tracing. The depth is read from the next line.

.. keyword:: TRACE_PROCEDURE
   :target: TRACE PROCEDURE
   :section: &CPMD

   Select a procedure to be traced. The procedure is read from the next line.

.. keyword:: TRAJECTORY
   :section: &CPMD
   :alternatives: [OFF, XYZ, DCD, SAMPLE, BINARY, RANGE, FORCES]

   Store the atomic positions, velocities and optionally forces at every *NTRAJ*
   time step on file *TRAJECTORY*. This is the **default for MD runs**. With the
   additional keyword XYZ the trajectory is also writthen in xyz-format on the
   file *TRAJEC.xyz*, similarly with the additional keyword DCD a trajectory in
   dcd-format (binary and single precision, as used by CHARMM, X-PLOR and other
   programs) is written on the file TRAJEC.dcd. If the keyword SAMPLE is given
   *NTRAJ* is read from the next line, otherwise the default value for *NTRAJ* is
   :math:`1`. A negative value of *NTRAJ* will disable output of the *TRAJECTORY*
   file, but e.g. TRAJEC.xyz will still be written every *-NTRAJ* steps. A value
   of 0 for *NTRAJ* will disable writing of the trajectory files altogether.

   The TRAJECTORY file is written in binary format if the keyword BINARY is
   present. If FORCES is specified also the forces are written together with the
   positions and velocities into the file FTRAJECTORY. It is possible to store the
   data of a subset of atoms by specifying the suboption RANGE, the smallest and
   largest index of atoms is read from the next line. If both, SAMPLE and RANGE
   are given, the RANGE parameters have to come before the SAMPLE parameter.

.. keyword:: TROTTER FACTOR
   :section: &CPMD

   Solve :math:`e^{-H/k_BT}` directly using **Trotter approximation**

   :math:`\left( e^{-pH} \simeq  e^{-pK/2}e^{-pV}e^{-pK/2}\right)`.

   The Trotter approximation is twice as fast.

   The Trotter factor is read from the next line (typically 0.001 is very
   accurate).

.. keyword:: TROTTER FACTORIZATION OFF
   :section: &CPMD

   Do not use Trotter factorization to calculate free energy functional.

   Remark: Place this keywords only after FREE ENERGY FUNCTIONAL; before it has no
   effect. Note: this keyword has *nothing* to do with path integral MD as
   activated by the keyword PATH INTEGRAL and as specified in the section &PIMD
   ... &END.

.. keyword:: USE_IN_STREAM
   :target: USE IN STREAM
   :section: &CPMD

   Specify that the RESTART file shall be read in the stream mode. This allows for
   parallel restarting (to be used with *USE_MPI_IO*).

.. keyword:: USE_MPI_IO
   :target: USE MPI IO
   :section: &CPMD

   Specify that MPI shall be used for parallel reading/writing of the RESTART
   file. This shall be used with *USE_IN_STREAM* and /or *USE_OUT_STREAM*.

.. keyword:: USE_MTS
   :target: USE MTS
   :section: &CPMD

   Switch on the Multiple Time-Step scheme for molecular dynamics (see &MTS
   section).\ [Liberatore2018]_

.. keyword:: USE_OUT_STREAM
   :target: USE OUT STREAM
   :section: &CPMD

   Specify that the RESTART file shall be written in the stream mode. This allows
   for parallel restarting (to be used with *USE_MPI_IO*).

.. keyword:: VDW CORRECTION
   :section: &CPMD
   :options: [ON,OFF]

   An empirical van der Waals correction scheme is applied to pairs of atom types
   specified with this keyword. This activates reading the corresponding
   parameters from the &VDW ... & END in which you have to specify all the VDW
   parameters between the opening and closing section keywords EMPIRICAL
   CORRECTION and END EMPIRICAL CORRECTION. Note that the two possible vdW
   options, EMPIRICAL CORRECTION and WANNIER CORRECTION are mutually exclusive.
   See :kw:`VDW PARAMETERS` for more details.

.. keyword:: VDW WANNIER
   :section: &CPMD

   A first-principle van der Waals correction scheme
   [psil1]_ [psil2]_ [psil3]_ is applied to selected groups of atoms on
   which maximally localized Wannier functions (WF) and centers (WFC) have been
   previously computed. The file WANNIER-CENTER generated upon WFC calculation
   must be present. This activates the reading procedure of the corresponding
   parameters from the &VDW ... &END section.

.. keyword:: VGFACTOR
   :section: &CPMD

   For :kw:`CDFT` runs read the inverse of the gradient optimiser step
   size (:math:`1/dx`) from the next line. The standard value of **10.0** should
   be fine in most situations.

.. keyword:: VIBRATIONAL ANALYSIS
   :section: &CPMD
   :alternatives: [FD, LR, IN], [GAUSS, SAMPLE, ACLIMAX]

   Calculate harmonic frequencies by finite differences of first derivatives
   **(FD)** (see also keyword :kw:`FINITE DIFFERENCES`), by
   **linear response** to ionic displacements **(LR)** or from a
   **pre-calculated** Hessian **(IN)**. K-point sampling is currently possible
   using finite differences. If the option GAUSS is specified, additional output
   is written on the file *VIB1.log* which contains the modes in a style similar
   to GAUSSIAN 98 output. This file can be read in and visualized with programs
   like MOLDEN or MOLEKEL. The option SAMPLE reads an integer from the next line.
   If this number is 2 an additional file *VIB2.log* containing the lowest modes
   is written. The **default** value is 1. If the option ACLIMAX is specified,
   additional output is written on the file VIB.aclimax which contains the modes
   in a style readable by aClimax
   (http://www.isis.rl.ac.uk/molecularspectroscopy/aclimax/). If a section
   **&PROP** is present with the keyword :kw:`DIPOLE MOMENT`\ [BERRY] or :kw:`DIPOLE MOMENT`\ [RS], the Born charge tensor is calculated on the
   fly. See also the block &LINRES ... &END and the keywords
   :kw:`RESTART` PHESS and :kw:`HESSIAN`
   {DISCO,SCHLEGEL,UNIT} PARTIAL.

.. keyword:: VMIRROR
   :section: &CPMD

   For :kw:`CDFT` HDA runs initialise :math:`V` for the second state as
   the negative final :math:`V` value of the first state. Useful in symmetric
   systems.

.. keyword:: WANNIER DOS
   :section: &CPMD

   Outputs the projected density of states of the Wannier orbitals (file
   WANNIER_DOS) and the KS Hamiltonian in the Wannier states representation (file
   WANNIER_HAM).

   When running :kw:`MOLECULAR DYNAMICS` CP the files
   WANNIER_DOS and WANNIER_HAM solely written at the last step.

.. keyword:: WANNIER MOLECULAR
   :section: &CPMD

   Generates effective molecular orbitals from the Wannier representation. It
   first attributes Wannier orbitals to molecules and then diagonalizes by
   molecular blocks the KS Hamiltonian.

   Does not work with :kw:`MOLECULAR DYNAMICS` CP.

.. keyword:: WANNIER NPROC
   :section: &CPMD

   Set the number of mpi tasks to be used for localization. Default is to use all
   the tasks available. The number of tasks is read from the next line and shall
   be a divisor of the number of tasks in a parallel run.

.. keyword:: WANNIER OPTIMIZATION
   :section: &CPMD
   :options: {SD,JACOBI,SVD}

   Use steepest descent or Jacobi rotation method for the orbital localization.

   Default are Jacobi rotations.

.. keyword:: WANNIER PARAMETER
   :section: &CPMD

   *W_STEP, W_EPS, W_RAN, W_MAXS* are read from the next line. *W_STEP* is the
   step size of the steepest descent algorithm used in the optimization procedure
   (default value 0.1). *W_EPS* the convergence criteria for the gradient (default
   value :math:`1.e-7`). *W_RAN* the amplitude for the initial random rotation of
   the states (default value 0.0). *W_MAXS* is the maximum steps allowed in the
   optimization (default value 200).

.. keyword:: WANNIER REFERENCE
   :section: &CPMD

   The vector *W_REF* is read from the next line, which consists of 3 coordinates
   :math:`x, y, z`. These are assumed as the origin for the WFCs positions and
   related ionic coordinates (i.e. :math:`{\bf R}_I \to {\bf R}_I-(x, y, z)`). The
   default value is the center of the supercell, if :kw:`CENTER MOLECULE` keyword is active (Note, that this is
   implicitly turned on, for calculations with :kw:`SYMMETRY` 0).
   Otherwise it is set to (0,0,0), which is usually not the center of the box. In
   order to get the best results displaying the IONS+CENTERS.xyz file this
   parameter should be set explicitly.

.. keyword:: WANNIER RELOCALIZE_EVERY
   :target: WANNIER RELOCALIZE EVERY
   :section: &CPMD

   If present, relocalize the wavefunction at every SCF step with the Jacobi
   method. The stride is read from the next line.

.. keyword:: WANNIER RELOCALIZE_IN_SCF
   :target: WANNIER RELOCALIZE IN SCF
   :section: &CPMD

   If present, relocalize/project the wavefunction at every SCF step.

.. keyword:: WANNIER SERIAL
   :section: &CPMD

   Requests that the calculation of Wannier functions is performed using the
   serial code, even in parallel runs.

.. keyword:: WANNIER TYPE
   :section: &CPMD
   :options: {VANDERBILT,RESTA}

   Indicates the type of Wannier functions. Vanderbilt type is the default.

.. keyword:: WANNIER WFNOUT
   :section: &CPMD
   :alternatives: [ALL,PARTIAL,LIST,DENSITY]

   Controls the printing of Wannier functions. Either all or only some of the
   functions can be printed. This will be done at the end of each calculation of
   Wannier functions. For **PARTIAL** output you have to give the indices of the
   first and the last Wannier function to print; the *LIST* directive follows the
   syntax of :kw:`RHOOUT` *BANDS*.

   ::

                  WANNIER WFNOUT PARTIAL
                     5  8

.. keyword:: WOUT
   :section: &CPMD
   :alternatives: [FULL]

   Controls the printing of the CDFT weight(s). If the keyword FULL is set the
   full weight is written out in the form of a density to WEIGHT-(suff), where
   (suff) is defined by the kind of the CDFT job. (suff)=WFOPT for single point
   calculations, while for geometry optimisations and MD two weights are written,
   (suff)=INIT at the beginning and (suff)=FINAL for the last step. If FULL is not
   set write out a slice of the weight in gnuplot readable form to
   WEIGHT-(suff).dat. Parameters WSLICE and WSTEP are read from the next line.

   WSLICE **0.5** is if larger than zero the z coordinate of the x-y weight plane
   to write out divided by the total box height. If WSLICE\ :math:`<0` the weight
   at the z coordinate of the first acceptor atom will be used.

   WSTEP **1** is the grid point step size for the output.

.. keyword:: {[TSDE, TSDP, TSDC]}
   :target: TSDE, TSDP, TSDC
   :section: &CPMD
   :alternatives: [NOPRECONDITIONING]
   :notes: NOPRECONDITIONING only electrons

   Short forms for the different :kw:`STEEPEST DESCENT`
   options.
