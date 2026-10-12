.. _input-tddft:

``&TDDFT`` section
==================

Input for TDDFT calculations. See also :doc:`/hints/tddft`.

The section starts with ``&TDDFT`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &TDDFT

Keywords
--------

.. keyword:: DAVIDSON PARAMETER
   :target: DAVIDSON PARAMETER TDDFT
   :section: &TDDFT

   The maximum number of Davidson iterations, the convergence criteria for the
   eigenvectors and the maximal size of the Davidson subspace are set. The three
   parameters *ndavmax, epstdav, ndavspace* are read from the next line.

   **Default** values are **100**, **10\ :math:`^{-10}`** and **10**.

.. keyword:: DAVIDSON RDIIS
   :section: &TDDFT

   This keyword controls the residual DIIS method for TDDFT diagonalization. This
   method is used at the end of a DAVIDSON diagonalization for roots that are not
   yet converged. The first number gives the maximum iterations, the second the
   maximum allowed restarts, and the third the maximum residual allowed when the
   method is invoked.

   **Default** values are **20**, **3** and **:math:`10^{-3}`**.

.. keyword:: DIAGONALIZER
   :section: &TDDFT
   :options: {DAVIDSON,NONHERMIT,PCG}
   :alternatives: [MINIMIZE]

   Specify the iterative diagonalizer to be used.

   **Defaults** are *DAVIDSON* for the Tamm–Dancoff method, *NONHERMIT* (a
   non-hermitian Davidson method) for TDDFT LR and *PCG* (Conjugate gradients) for
   the optimized subspace method. The additional keyword *MINIMIZE* applies to the
   PCG method only. It forces a line minimization with quadratic search.

   **Default** is **not to use line minimization**.

.. keyword:: EXTPOT
   :section: &TDDFT

   Non adiabatic (nonadiabatic, non-adiabatic) Tully’s trajectory surface hopping
   dynamics using TDDFT energies and forces, coupled with an external
   field [tavernelli2010]_. To be used together with the keywords
   :kw:`MOLECULAR DYNAMICS` BO, :kw:`TDDFT` in the
   &CPMD section, and :kw:`T-SHTDDFT` in the &TDDFT section. Do NOT
   use the keyword :kw:`T-SHTDDFT` together with the keyword :kw:`SURFACE HOPPING` in &CPMD, which invokes the SH scheme based on
   :kw:`ROKS` [surfhop]_ (see :kw:`SURFACE HOPPING`).

   This keyword follow the same principle as described for the keyword
   :kw:`T-SHTDDFT`, except that, in the present dynamics, the
   trajectory starts on the ground state and is coupled with an external field
   through the equations of motion for the amplitudes of Tully’s trajectory
   surface hopping. According to the evolution of the amplitudes of the different
   excited states, the running trajectory can jump on an excited state. From
   there, deactivation through nonradiative processes is possible, within the
   normal trajectory surface hopping scheme.

   Parameter *aampl*, *adir*, *afreq*, and *apara1* are read from the next line.
   The amplitude of the vector potential is provided in *aampl* and its
   polarization is given in *adir* (1 = x-polarized, 2 = y-polarized, 3 =
   z-polarized, 4 = all components). The keyword *afreq* gives the frequency of
   the field and *apara1* is a free parameter for a specific user-specified pulse.

   Important points: the applied electromagnetic field needs to be hard coded in
   the subroutine sh_tddft_utils.mod.F90, in the subroutine SH_EXTPOT. The vector
   potential is used for the coupling with the amplitudes equations. Be careful to
   use a time step small enough for a correct description of the pulse. The pulse
   is printed in the file SH_EXTPT.dat (step, A(t), E(t)).

.. keyword:: FORCE STATE
   :section: &TDDFT

   The state for which the forces are calculated is read from the next line.
   Default is for state 1.

.. keyword:: ISC-SHTDDFT
   :section: &TDDFT

   Non-adiabatic surface hopping dynamics including intersystem crossing
   (ISC) through spin-orbit coupling (SOC). Two integers are read from the
   next line: ``nskip socskip``. ``nskip`` is the interval (in MD steps) at
   which the ISC step (``do_tsh_isc``, Landau-Zener type hopping between
   spin manifolds) is performed, ``socskip`` the interval at which the
   spin-orbit couplings are recomputed; both are also performed whenever a
   surface hop occurs. **Defaults** are 1 and 1. The keyword implies
   :kw:`T-SHTDDFT` (Tully surface hopping) and requires singlet and triplet
   response states to be computed.

   Example::

      ISC-SHTDDFT
       1 1

   .. todo::

      MISSING: required companion keywords (e.g. :kw:`STATES` with SINGLET
      and TRIPLET, SOC settings in ``&CPMD``), physical reference.

.. source: src/lr_in_utils.mod.F90:795

.. keyword:: LOCALIZATION
   :section: &TDDFT

   Use localized orbitals in the TDDFT calculation. Default is to use canonical
   orbitals.

.. keyword:: LR-TDDFT
   :section: &TDDFT

   Use full linear response version of TDDFT. Default is to use
   :kw:`TAMM-DANCOFF` approximation.

.. keyword:: LZ-SHTDDFT
   :section: &TDDFT

   Non adiabatic (nonadiabatic, non-adiabatic) trajectory surface hopping scheme
   based on Landau-Zener transition probabilities. Excited state dynamics is
   performed on TDDFT potential energy surfaces. To be used together with the
   keywords :kw:`MOLECULAR DYNAMICS` BO and
   :kw:`TDDFT` in the &CPMD section (see
   section :ref:`sec:TDDFTdynamics`). Do NOT use this keyword together with
   :kw:`SURFACE HOPPING` in &CPMD, which invokes the surface
   hopping scheme based on :kw:`ROKS` (see :kw:`SURFACE HOPPING`). See also the related approach based on
   Tully’s hopping probabilities :kw:`T-SHTDDFT`.

.. keyword:: MOLECULAR STATES
   :section: &TDDFT

   Calculate and group Kohn–Sham orbitals into molecular states for a TDDFT
   calculation.

.. keyword:: NACV_FULL
   :section: &TDDFT

   For the calculation of non-adiabatic coupling vectors (``NACV`` in
   ``&CPMD``) include the linear-response contribution. By default only the
   direct term is computed. No further input is read.

   .. todo::

      MISSING: formula/reference for the two contributions and the cost
      implications.

.. source: src/lr_in_utils.mod.F90:793

.. keyword:: OSCILLATOR STRENGTH
   :section: &TDDFT
   :options: {QP, BF, BFH, BERRY}

   Requests oscillator strengths. Options on the same line: ``QP`` (only
   Kohn-Sham oscillator strengths, no TDDFT), ``BF`` (transition dipoles by
   Bessel-function integration), ``BFH`` ("Hutter-Bessel" variant; switches
   off ``BF``) and ``BERRY`` (generalized Berry-phase formula; switches off
   the other options). The selection is echoed in the TDDFT summary of the
   output.

   In the current code the flags set by this keyword are only used to print
   that summary; they are not referenced anywhere else.

   .. todo::

      MISSING: whether this keyword has any effect on the calculation; how
      oscillator strengths are actually computed and printed.

.. source: src/lr_in_utils.mod.F90:767

.. keyword:: OUTPUT
   :target: OUTPUT TDDFT
   :section: &TDDFT
   :alternatives: [NORMAL, FULL, MINIMAL]

   Print level of the TDDFT routines used for excited-state forces
   (optimization, molecular dynamics, finite-difference checks):
   ``NORMAL`` (**Default**), ``FULL`` or ``MINIMAL``.

.. source: src/lr_in_utils.mod.F90:781

.. keyword:: PCG PARAMETER
   :section: &TDDFT

   The parameters for the PCG diagonalization are read from the next line. If
   *MINIMIZE* was used in the :kw:`DIAGONALIZER` then the total
   number of steps (default 100) and the convergence criteria (default
   :math:`10^{-8}`) are read from the next line. Without minimization in addition
   the step length (default 0.5) has also to be given.

.. keyword:: PROPERTY
   :section: &TDDFT
   :options: { STATE }

   Calculate properties of excited states at the end of an :kw:`ELECTRONIC SPECTRA` calculations. default is to calculate
   properties for all states. Adding the keyword **STATE** allows to restrict the
   calculation to only one state. The number of the state is read from the next
   line.

.. keyword:: RANDOMIZE
   :target: RANDOMIZE TDDFT
   :section: &TDDFT

   Randomize the initial vectors for the diagonalization in a TDDFT calculation.
   The amplitude is read from the next line. Default is not to randomize the
   vectors.

.. keyword:: REORDER
   :section: &TDDFT

   Reorder the canonical Kohn–Sham orbitals prior to a TDDFT calculation. The
   number of states to be reordered is read from the next line. On the following
   line the final rank of each states has to be given. The first number given
   corresponds to the HOMO, the next to the HOMO - 1 and so on. All states down to
   the last one changed have to be specified, no holes are allowed. This keyword
   can be used together with the *SUBSPACE* option in a
   :kw:`TAMM-DANCOFF` calculation to select arbitrary states.
   Default is to use the ordering of states according to the Kohn–Sham
   eigenvalues.

.. keyword:: REORDER LOCAL
   :section: &TDDFT

   Reorder the localized states according to a distance criteria. The number of
   reference atoms is read from the next line. On the following line the position
   of the reference atoms within the set of all atoms has to be given. The keyword
   :kw:`LOCALIZE` is automatically set. The minimum distance of the
   center of charge of each state to the reference atoms is calculated and the
   states are ordered with respect to decreasing distance. Together with the
   *SUBSPACE* option in a :kw:`TAMM-DANCOFF` calculation this can
   be used to select specific states for a calculation.

.. keyword:: ROTATION PARAMETER
   :section: &TDDFT

   The parameters for the orbital rotations in an optimized subspace calculation
   (see :kw:`TAMM-DANCOFF`) are read from the next line. The total
   number of iterations (default 50), the convergence criteria (default
   :math:`10^{-6}`) and the step size (default 0.5) have to be given.

.. keyword:: SH_LCT
   :section: &TDDFT

   Local control theory (LCT) applied to TDDFT surface hopping: an external
   field is constructed on the fly to drive the population towards a target
   state. Five values are read from the next line:
   ``seed_t seed_m lct_l adir lct_state``. Up to the time ``seed_t`` a
   constant seed field of amplitude ``seed_m`` is applied; afterwards the
   field amplitude is ``-lct_l`` times the transition dipole of the target
   state times the imaginary part of the product of the amplitudes of the
   ground and target states. ``adir`` selects the polarization (1 = x,
   2 = y, 3 = z, 4 = all components) and ``lct_state`` the target excited
   state. **Defaults** are 100, 1, 100, (none) and 1. The field is written
   to ``SH_LCON.dat``.

   Example::

      SH_LCT
       100.0 0.001 100.0 3 1

   .. todo::

      MISSING: units of the parameters (time presumably in a.u.), whether
      :kw:`T-SHTDDFT` must be given explicitly, literature reference.

.. source: src/lr_in_utils.mod.F90:813

.. keyword:: STATES
   :target: STATES TDDFT
   :section: &TDDFT
   :options: {MIXED,SINGLET,TRIPLET}

   The number of states to be calculated is read from the next line. The type of
   state *SINGLET, TRIPLET* can be given for non-spinpolarized calculations.
   Default is to calculate one singlet state for LDA and 1 mixed state for LSD
   calculations.

.. keyword:: STORE_V_AB
   :section: &TDDFT

   For TDDFT with hybrid functionals, store the exchange pair potentials
   :math:`v_{ab}` computed for the response kernel in memory (an array of
   size grid points :math:`\times` states :math:`\times` states) so that they
   are not recomputed. No further input is read. See also :kw:`WRITE_V_AB`.

.. source: src/lr_in_utils.mod.F90:696

.. keyword:: T-SHTDDFT
   :section: &TDDFT

   Non adiabatic (nonadiabatic, non-adiabatic) Tully’s trajectory surface hopping
   dynamics using TDDFT energies and forces. To be used together with the keywords
   :kw:`MOLECULAR DYNAMICS` BO and :kw:`TDDFT` in
   the &CPMD section (see section :ref:`sec:TDDFTdynamics`). Do NOT use the
   keyword :kw:`T-SHTDDFT` together with the keyword :kw:`SURFACE HOPPING` in &CPMD, which invokes the SH scheme based on
   :kw:`ROKS` [surfhop]_ (see :kw:`SURFACE HOPPING`).

   For a given initial configuration, the run produces a trajectory that undergoes
   surface hopping according to the algorithm by Tully adapted to
   TDDFT [TDDFT-SH]_. The forces on the excited state surfaces are
   computed using TDDFT as for the adiabatic case. A sufficiently large number of
   excited states must be declared using the keyword :kw:`STATES` in the
   section &TDDFT. The initial running surface is specified with the keyword
   :kw:`FORCE STATE` in the section &TDDFT. This can change during
   the dynamics when a surface hop occurs. After a restart the value of the
   running state is taken from the file SH_STATE.dat (see below). The run produces
   a series of additional files: SH_COEFA.dat (absolute value of the state
   amplitudes), SH_COEFC.dat (their complex values), SH_COUPL.dat (the coupling
   strength per state), SH_ENERG.dat (the energy of the different states: setp
   number, ground state energy, first excited state energy, :math:`\dots`, highest
   excited state energy, energy of the running state), SH_PROBS.dat (transition
   probabilities between running state and all other states), SH_STATE.dat (the
   running state at each step). All these files (in addition to SH_WAVEFUNCTIONS
   and SH_LRWAVEFUNCTIONS) are needed to restart the SH dynamics. Note that each
   run produces a single SH trajectory. Several subsequent runs starting from
   different initial coordinates and velocities are required to collect
   statistics.

.. keyword:: TAMM-DANCOFF
   :section: &TDDFT
   :alternatives: [SUBSPACE,OPTIMIZE]

   Use the Tamm–Dancoff approximation. This is the default for TDDFT calculations.
   Optionally, only a *SUBSPACE* of the occupied orbitals can be included in the
   calculation. The subspace can be optimized at each step (not recommended).
   Default is to use all states.

.. keyword:: TD_METHOD_A
   :target: TD METHOD A
   :section: &TDDFT
   :alternatives: [{\em functionals}]

   Use a different potential for the eigenvalue difference part of the response
   equations than was used to generate the ground state orbitals. The potential
   generating functional has to be given after the keyword. For possible
   functionals see the code. Most likely you want to use the **SAOP** functional.

   This functional does not affect the choice of functional used in the TDDFT
   kernel. The kernel functional is set in the &DFT section. It is either the
   standard functional or the functional defined by the keyword :kw:`LR KERNEL`.

.. keyword:: USE_LIN_LIN
   :section: &TDDFT

   Use the Lin Lin scheme (adaptively compressed exchange) to speed up the
   exact-exchange part of the TDDFT response kernel. No further input is
   read. Not implemented for :kw:`LSD` (the program stops). During the
   Davidson part of the calculation the scheme is temporarily switched off.

   .. todo::

      MISSING: literature reference, accuracy/performance considerations.

.. source: src/lr_in_utils.mod.F90:703

.. keyword:: WRITE_V_AB
   :section: &TDDFT

   Like :kw:`STORE_V_AB`, but the exchange pair potentials are written to
   (and read back from) an unformatted scratch file per MPI task (Fortran
   unit ``2000+rank``) instead of being kept in memory. No further input is
   read.

   .. todo::

      MISSING: the actual file names (compiler-dependent default names for
      unit 2000+rank) and disk-space requirements.

.. source: src/lr_in_utils.mod.F90:699
