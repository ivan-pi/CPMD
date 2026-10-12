.. _input-resp:

``&RESP`` section
=================

Response calculations. This section is always evaluated, even if it is not used. See also :doc:`/hints/linear-response`.

The section starts with ``&RESP`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &RESP

Keywords
--------

.. keyword:: ATOMICRHO
   :section: &RESP

   Option of :kw:`INTERACTION`: use an atomic density instead of the Wannier
   density for the reference density :math:`\rho_0`. Only recognized after
   the :kw:`INTERACTION` keyword and only if preceded by a blank.

   .. todo::

      MISSING: how the atomic density is built.

.. source: src/respin_p_utils.mod.F90:999

.. keyword:: ATOMICWAVEFUNCTION
   :section: &RESP

   Option of :kw:`INTERACTION`: use atomic wavefunctions instead of Wannier
   functions for the unperturbed orbitals :math:`W_0`. Implies that no
   ground-state optimization is done and no restart file is read. Only
   recognized after the :kw:`INTERACTION` keyword and only if preceded by a
   blank.

   .. todo::

      MISSING: source of the atomic wavefunctions.

.. source: src/respin_p_utils.mod.F90:996

.. keyword:: ATOMLIST
   :section: &RESP

   Restricts a :kw:`VOA` calculation in ``AT`` mode to a subset of atoms. The
   atom indices are given on the same line as four-digit, zero-padded
   numbers (e.g. ``0001 0012``); an atom is selected if its formatted index
   appears anywhere in the line. Only recognized if preceded by a blank
   (indent the line or put it on the :kw:`VOA` line). Because the line is
   limited to 80 characters, only about 14 atoms can be listed.

   **Note:** if ``ATOMLIST`` is placed on the :kw:`VOA` line, the substring
   test for the ``AT`` option also matches; with ``MD`` this leads to the
   error "AT and MD are not compatible".

   Example::

      VOA AT
       ATOMLIST 0001 0002 0005

.. source: src/respin_p_utils.mod.F90:762

.. keyword:: CG-ANALYTIC
   :section: &RESP

   | The number of steps for which the step length in the conjugate gradient
     optimization is calculated assuming a quadratic functional E(2) (quadratic in
     the linear response vectors). No accuracy impact, pure convergence speed
     tuning.
   | **Default** value is **3** for NMR and **99** otherwise.

.. keyword:: CG-FACTOR
   :section: &RESP

   | The analytic length calculation of the conjugate-gradient step lengths yields
     in general a result that is slightly too large. This factor is used to
     correct for that deficiency. No accuracy impact, pure convergence speed
     tuning.
   | **Default** is **0.8**.

.. keyword:: CONVERGENCE
   :target: CONVERGENCE RESP
   :section: &RESP

   Convergence criterion on the gradient :math:`\delta E/\delta \psi^*`
   **Default** value is **10\ :math:`^{-5}`**.

.. keyword:: DIIS
   :section: &RESP

   Use (preconditioned) DIIS instead of conjugate gradients to optimize the
   linear-response wavefunctions. The size of the DIIS buffer is set with
   :kw:`MDIIS` (**Default** 10). See also :kw:`PCG_FIRST`. Any line containing
   ``DIIS`` as a substring selects this optimizer.

.. source: src/respin_p_utils.mod.F90:1060

.. keyword:: DISCARD
   :section: &RESP
   :alternatives: [OFF, PARTIAL, TOTAL, LINEAR]

   Request to discard trivial modes in vibrational analysis from linear response
   (both :kw:`PHONON` and :kw:`LANCZOS`).

   **OFF** = argument for performing no projection.

   **PARTIAL** = argument for projecting out only translations (this is the
   default).

   **TOTAL** = argument for projecting both rotations and translations.

   **LINEAR** = argument for projecting rotations around the :math:`C - \infty`
   axis in a linear molecule (not implemented yet).

.. keyword:: DX
   :section: &RESP

   Option of :kw:`INTERACTION`: read a displacement vector of the Wannier
   functions for each molecule. :kw:`NMOL` lines with three reals each are
   read. Only recognized after :kw:`INTERACTION`; the number of molecules
   must be set before with :kw:`NMOL`. Since the test is a substring test,
   any line containing ``DX`` triggers this keyword in that context.

   .. todo::

      MISSING: units (presumably Bohr) and how the displacement is applied.

.. source: src/respin_p_utils.mod.F90:916

.. keyword:: EIGENSYSTEM
   :section: &RESP

   Not documented.

.. keyword:: EPR
   :section: &RESP
   :notes: options

   Calculate the EPR :math:`g` tensor for the system. This routine accepts most,
   if not all, of the options available in the NMR routine (RESTART, NOSMOOTH,
   NOVIRTUAL, PSI0, RHO0, OVERLAP and FULL). Most important new options are:

   **FULL SMART**: does a calculation with improved accuracy. A threshold value
   (between 0 and 1) must be present on the next line. The higher the threshold
   value, the lower the computational cost, but this will also reduce the accuracy
   (a bit). Typically, a value of 0.05 should be fine.

   **OWNOPT**: for the calculation of the :math:`g` tensor, an effective potential
   is needed. By default, the EPR routine uses the local potential
   (:math:`V_{LOC} = V_{PP,LOC} + V_{HARTREE} + V_{XC}`). This works well with
   Goedecker pseudopotentials, but rather poor with Troullier-Martins
   pseudopotentials. When using this option, the following potential is used
   instead:

   .. math:: V_{EFF} = -\frac{Z}{r}\mathrm{erf}(r/r_c) + V_{HARTREE} + V_{XC}

   and :math:`r_c` (greater than 0) is read on the next line.

   **HYP**: calculates the hyperfine tensors. See epr_hyp_utils.mod.F90 for
   details.

   Contact Reinout.Declerck@UGent.be should you require further information.

.. keyword:: FASTDEXC
   :section: &RESP

   Option of :kw:`INTERACTION`: use the LDA to compute the second derivative
   of the exchange-correlation energy :math:`d^2E_{xc}/dn_1 dn_2` in
   perturbation calculations with gradient-corrected functionals. Only
   recognized after :kw:`INTERACTION` and only if preceded by a blank.

.. source: src/respin_p_utils.mod.F90:989

.. keyword:: FUKUI
   :section: &RESP
   :alternatives: [N=nf]

   Calculates the response to a change of occupation number of chosen orbitals.
   The indices of these orbitals are read from the following nf lines (**default
   nf=1**). The orbitals themselves are not read from any :kw:`RESTART`
   file but from WAVEFUNCTION.\* files generated with :kw:`RHOOUT` in
   the &CPMD section; to recall this the orbital numbers have to be negative, just
   like for the :kw:`RHOOUT` keyword.

   A weight can be associated with each orbital if given just after the orbital
   number, on the same line. It corresponds to saying how many electrons are put
   in or taken from the orbital. For example;

   ::

         FUKUI N=2
         -i 1.0
         -j -1.0

   corresponds to the response to taking one electron from orbital i and put it in
   orbital j.

.. keyword:: HAMILTONIAN CUTOFF
   :target: HAMILTONIAN CUTOFF RESP
   :section: &RESP

   The value where the preconditioner (the approximate Greens function
   :math:`(V_{kin}+V_{coul}-\epsilon_{KS})^{-1})` is truncated. HTHRS can also be
   used. Default value is 0.5.

.. keyword:: HARDNESS
   :section: &RESP

   Not documented.

.. keyword:: INTERACTION
   :section: &RESP

   Not documented.

.. keyword:: ISOLATED
   :section: &RESP

   For :kw:`NMR` / :kw:`EPR`: assign overlapping orbitals to the same
   virtual cell using an overlap threshold of 0.1 (equivalent to
   ``OVERLAP`` with 0.1) and switch off the ``FULL`` calculation.
   Intended for isolated systems. Detected as the substring ``ISOLAT``.

.. source: src/respin_p_utils.mod.F90:729

.. keyword:: KEEPREALSPACE
   :section: &RESP

   Like the standard CPMD option, this keeps the C0 ground state wavefunctions in
   the direct space representation during the calculation. Can save a lot of time,
   but is incredibly memory intensive.

.. keyword:: KPERT
   :section: &RESP
   :alternatives: [MONKHORSTPACK,SCALE]

   Calculation of total energy and electronic density of states with an arbitrary
   number of k-points (at almost no additional computational effort). The method
   is based on a :math:`{\bf k}\cdot {\bf p}-`\ like approximation developed in
   the framework of the density functional perturbation theory
   [mimp]_. For a sampling of the BZ determined by the Monkhorst-Pack
   algorithm, the option **MONKHORSTPACK** has to be specified, followed by the
   dimension of the mesh along the 3 reciprocal space axis
   :math:`(NK_{1} , NK_{2} , NK_{3})`. If omitted, the individual absolute
   coordinates of the k-points have to be given one by one in the following lines.
   The :kw:`SCALE` option allows to specify them in units of the
   reciprocal cell vectors.

   The line after **KPERT** has to contain the total number of k-points
   :math:`(NKPTS)`, which have then to be given by their coordinates and the
   associated weights :math:`(RK,WK)` in the format:

   :math:`NKPTS`

   :math:`RK_{x1} \; RK_{y1} \; RK_{z1} \; WK_{1}`

   :math:`\dots`

   :math:`RK_{x NKPTS} \; RK_{y NKPTS} \; RK_{z NKPTS} \; WK_{NKPTS}`.

   Three response wavefunctions are calculated, corresponding to the three
   independent orientations of the k basis vectors in reciprocal space. Therefore,
   3 independent optimization loops are started (:math:`x,y` and :math:`z`), and
   the 3 sets of wfns are stored (you need 4 times the memory required for a
   standard wavefunction optimization). The second order correction to the
   :math:`\Gamma`-point total energy is calculated for the requested k-point mesh.

   Further options are (each in a new line of the input file ):

   WRITE_C1
      the 3 sets of response wfns are stored in three separate restart files.

   HAMILTONIAN
      the k-dependent Hamiltonian is constructed via the second order perturbation
      theory approximation, and the corresponding KS energies are calculated. Due
      to technical reasons, for each k-point :math:`2*NSTATE` KS energies are
      calculated, however only those corresponding to occupied orbitals are
      reliable.

   READ_C1
      the response wfns are read from RESTART.P\_{xyz}.

   BUILD_C00
      the set of k-dependent wfns (first order correction) is calculated from the
      unperturbed :math:`\Gamma`-point wfns together with the response orbitals.
      They are then written in a standard :kw:`RESTART` file. From this
      restart file one can perform a calculation of the Hamiltonian matrix for
      each kpoint and calculate the KS energies (use **LANCZOS DIAGO** in &CPMD
      and the **KPOINT** option **ONLYDIAG** in &SYSTEM. The k-point mesh must be
      the same used in the linear response calculation. set also **NOSPHERICAL
      CUTOFF** in &SYSTEM).

   NORESTART
      no RESTART file is written.

.. keyword:: LANCZOS
   :section: &RESP
   :alternatives: [CONTINUE,DETAILS]

   lanczos_dim iterations conv_threshold lanczos_dim= dimension of the vibrational
   d.o.f. iterations = no.of iterations desired for this run conv_threshold =
   threshold for convergence on eigenvectors CONTINUE = argument for continuing
   Lanczos diagonalization from a previous run (reads file LANCZOS_CONTINUE)
   DETAILS = argument for verbosity. prints a lot of stuff

.. keyword:: LINEARSCALING
   :section: &RESP

   Option of :kw:`INTERACTION`: the response orbitals are only required to
   be orthogonal to "their own" ground-state orbital,
   :math:`\langle\psi^{(1)}_i|\psi^{(0)}_i\rangle = 0`, and not to the other
   ground-state orbitals. Only recognized after :kw:`INTERACTION` and only if
   preceded by a blank.

.. source: src/respin_p_utils.mod.F90:837

.. keyword:: LMAX
   :section: &RESP

   Option of :kw:`INTERACTION`: maximum angular momentum used for the
   rotations of Wannier functions. The integer is read from the next line.
   **Default** is 4.

.. source: src/respin_p_utils.mod.F90:850

.. keyword:: MANNOPT
   :section: &RESP

   Option of :kw:`INTERACTION`: use the old serial implementation of the
   perturbative interaction code instead of the parallel one. In parallel
   runs the option is switched off with a warning. Only if preceded by a
   blank.

.. source: src/respin_p_utils.mod.F90:1008

.. keyword:: MDIIS
   :section: &RESP

   Use DIIS for the linear-response optimization (as :kw:`DIIS`) and read
   the size of the DIIS buffer from the next line. **Default** is 10.

   Example::

      MDIIS
       6

.. source: src/respin_p_utils.mod.F90:1065

.. keyword:: MULTIREF
   :section: &RESP

   Option of :kw:`INTERACTION`: allow a separate reference Wannier-function
   file for each molecule.

   .. todo::

      MISSING: names of the reference files.

.. source: src/respin_p_utils.mod.F90:960

.. keyword:: NMOL
   :section: &RESP

   Option of :kw:`INTERACTION`: number of molecules in the system, read from
   the next line. **Default** is 1. Must precede :kw:`THETA`, :kw:`VEC` and
   :kw:`DX`, which read one line per molecule.

.. source: src/respin_p_utils.mod.F90:856

.. keyword:: NMR
   :section: &RESP
   :notes: options

   Calculate the NMR chemical shielding tensors for the system. Most important
   option: FULL, does a calculation with improved accuracy for periodic systems
   but takes a lot of time. Isolated systems: Use OVERLAP and 0.1 (on next line)
   for the same effect. *Be careful for non-hydrogen nuclei.* The shielding is
   calculated without contribution from the core electrons. Contact
   sebastia@mpip-mainz.mpg.de for further details.

.. keyword:: NOLOC
   :section: &RESP

   Option of :kw:`INTERACTION` (on its own line): disables the localization
   of the optimized wavefunctions (same flag as the ``NOLOC`` option of
   :kw:`NMR`). Only recognized after :kw:`INTERACTION` and only if preceded
   by a blank.

.. source: src/respin_p_utils.mod.F90:976

.. keyword:: NOMIN
   :section: &RESP

   Switch off the analytic line minimization in the conjugate-gradient
   optimization of the response wavefunctions, i.e. use a fixed step length
   (the output shows "using FIXED STEP LENGTH").

   .. todo::

      MISSING: which step length is used.

.. source: src/respin_p_utils.mod.F90:1053

.. keyword:: NONORTHOGWANNIER
   :section: &RESP

   Option of :kw:`INTERACTION`: use a non-orthogonal combination of Wannier
   functions (by default the Wannier orbitals are orthogonalized). Only if
   preceded by a blank.

.. source: src/respin_p_utils.mod.F90:1002

.. keyword:: NOOPT
   :section: &RESP

   Do not perform a ground state wfn optimization. Be sure the restarted wfn is at
   the BO-surface.

.. keyword:: NOROT
   :section: &RESP

   Option of :kw:`INTERACTION`: a rotation by zero degrees is performed
   without projection (by default all rotations of the Wannier functions are
   carried out).

.. source: src/respin_p_utils.mod.F90:950

.. keyword:: NSAVESTATE
   :section: &RESP

   Option of :kw:`INTERACTION`: number of states displayed and saved at the
   end of the perturbative interaction calculation, read from the next line.
   **Default** is 10.

.. source: src/respin_p_utils.mod.F90:862

.. keyword:: OACP
   :section: &RESP
   :alternatives: [DENSITY, REF_DENSITY, FORCE]

   Not documented.

.. keyword:: ORTHOGONAL
   :section: &RESP

   Option of :kw:`INTERACTION`: enforce orthogonality of the perturbation
   orbitals among themselves. Only if preceded by a blank.

.. source: src/respin_p_utils.mod.F90:844

.. keyword:: OWNPOT
   :section: &RESP

   For :kw:`EPR`: use the effective potential
   :math:`V_{EFF} = -\frac{Z}{r}\mathrm{erf}(r/r_c) + V_{HARTREE} + V_{XC}`
   instead of the local potential for the calculation of the :math:`g`
   tensor. The smoothing radius :math:`r_c` (must be greater than zero) is
   read from the next line. (Documented in the :kw:`EPR` entry under the
   misspelled name OWNOPT.)

   Example::

      OWNPOT
       0.5

.. source: src/respin_p_utils.mod.F90:800

.. keyword:: PCG_FIRST
   :section: &RESP

   With :kw:`DIIS`, perform the first optimization step with preconditioned
   conjugate gradients before switching to DIIS.

.. source: src/respin_p_utils.mod.F90:1056

.. keyword:: PCG_P
   :section: &RESP

   Use preconditioned conjugate gradients for the optimization of the
   response wavefunctions (this is the **default**). Overrides an earlier
   :kw:`DIIS` or :kw:`MDIIS`.

.. source: src/respin_p_utils.mod.F90:1083

.. keyword:: PHONON
   :section: &RESP

   Calculate the harmonic frequencies from perturbation theory.

.. keyword:: POLAK
   :section: &RESP

   Uses the Polak-Ribiere formula for the conjugate gradient algorithm. Can be
   safer in the convergence.

.. keyword:: RAMAN
   :section: &RESP

   Calculate the polarizability (also in periodic systems) as well as Born-charges
   and dipole moment.

.. keyword:: RESTART
   :target: RESTART RESP
   :section: &RESP

   Restart a response calculation from previously written response files
   (:kw:`NMR`: ``RESTART.NMR``, ``RESTART.p_x`` ...; also used by :kw:`EPR`).
   The keyword is detected as a substring anywhere in a line of ``&RESP``
   and also implies :kw:`NOOPT` (no ground-state optimization; the
   wavefunction is read from the restart file).

.. source: src/respin_p_utils.mod.F90:441

.. keyword:: RESTRAINED
   :section: &RESP

   Option of :kw:`INTERACTION`: the plane-wave cutoff is lowered for the
   perturbation calculation (``cutoff_restr``).

   .. todo::

      MISSING: by how much the cutoff is lowered and how this differs from
      :kw:`TRICK`.

.. source: src/respin_p_utils.mod.F90:970

.. keyword:: SAVELOCW0
   :section: &RESP

   Option of :kw:`INTERACTION`: save the orthogonalized unperturbed Wannier
   functions :math:`W_0` in a RESTART file and quit. Only if preceded by a
   blank.

.. source: src/respin_p_utils.mod.F90:992

.. keyword:: SCFPTTOL
   :section: &RESP

   Option of :kw:`INTERACTION`: convergence tolerance of the self-consistent
   perturbation-theory cycles (see :kw:`SELFCONSISTENT`), read from the next
   line. **Default** is :math:`10^{-6}`.

.. source: src/respin_p_utils.mod.F90:873

.. keyword:: SELFCONSISTENT
   :section: &RESP

   Option of :kw:`INTERACTION`: number of self-consistent iterations for
   the optimization of the first-order Wannier functions :math:`W_1`, read
   from the next line. **Default** is 1.

.. source: src/respin_p_utils.mod.F90:868

.. keyword:: SIGMA
   :section: &RESP

   For :kw:`NMR` / :kw:`EPR`: spread of the Gaussian used to determine the
   virtual cells. The lower left corner of the virtual cell of each orbital
   is placed where the integral of the orbital density over the Gaussian
   :math:`\exp(-(x-x_0)^2/2\sigma^2)` is smallest. The value is read from the
   next line. **Default** is 0.1. Should be used with care, results may
   depend on it.

   .. todo::

      MISSING: units of sigma.

.. source: src/respin_p_utils.mod.F90:1104

.. keyword:: SIMPLEMODEL
   :section: &RESP

   Option of :kw:`INTERACTION`: use a simplified interaction model (a
   simplified orthogonalization of the response orbitals is used in the
   optimizer). Only if preceded by a blank.

   .. todo::

      MISSING: description of the simplified model.

.. source: src/respin_p_utils.mod.F90:1005

.. keyword:: SMART
   :section: &RESP

   For :kw:`EPR` with ``FULL``: perform the full calculation with a
   threshold that reduces its cost. The threshold is read from the next
   line and must lie in the interval [0, 0.1] (the program stops
   otherwise). A higher value lowers the cost and the accuracy; the manual
   recommends 0.05. (Documented in the :kw:`EPR` entry as ``FULL SMART``.)

   Example::

      FULL SMART
       0.05

.. source: src/respin_p_utils.mod.F90:812

.. keyword:: STATEPREC
   :section: &RESP

   Use a state-dependent preconditioner for the response optimization: for
   each state :math:`i` the diagonal :math:`\eta_{GG} = H_{GG}-\epsilon_i`
   is smoothed as :math:`(\eta^2+\delta^2)^{-1/2}` (cf. :kw:`HAMILTONIAN CUTOFF`).
   Mutually exclusive with :kw:`TIGHTPREC` (which takes precedence).

.. source: src/respin_p_utils.mod.F90:1047

.. keyword:: STEEPEST DESCENT
   :target: STEEPEST DESCENT RESP
   :section: &RESP

   A line containing both ``STEE`` and ``DESC`` selects steepest descent for
   the response optimization, which is no longer supported: the program
   stops immediately with an error message.

.. source: src/respin_p_utils.mod.F90:1075

.. keyword:: THETA
   :section: &RESP

   Option of :kw:`INTERACTION`: rotation angle of the Wannier functions for
   each molecule; :kw:`NMOL` lines with one real each are read.

   .. todo::

      MISSING: units (degrees or radians).

.. source: src/respin_p_utils.mod.F90:879

.. keyword:: TIGHTPREC
   :section: &RESP

   Uses a harder preconditioner. For experts: The Hamiltonian is approximated by
   the kinetic energy, the G-diagonal Coulomb potential and the KS-energies. The
   number obtained this way must not be close to zero. This is achieved by
   smoothing it with This is achieved by smoothing it with

   .. math:: x \to f(x) = \sqrt{x^2 + \epsilon^2} \; \; [{\rm default}]

   or

   .. math:: x \to f(x) = (x^2 + \epsilon ^2)/x \; \; [{\rm this \; option}]

   The HARD option conserves the sign of the approximate Hamiltonian whereas the
   default formula does never diverge.

.. keyword:: TIGHTSTATEPREC
   :section: &RESP

   Intended to select a state-dependent preconditioner with the
   :math:`\eta/(\eta^2+\delta^2)` smoothing (cf. :kw:`TIGHTPREC`).

   **Warning:** because the pattern ``STATEPREC`` is tested first and is a
   substring of ``TIGHTSTATEPREC``, this keyword currently selects
   :kw:`STATEPREC` instead; the tight state-dependent variant is not
   reachable from the input.

.. source: src/respin_p_utils.mod.F90:1050

.. keyword:: TRICK
   :section: &RESP

   Option of :kw:`INTERACTION`: "cutoff trick", a lower cutoff is used for
   the perturbation wavefunctions.

   .. todo::

      MISSING: the cutoff used.

.. source: src/respin_p_utils.mod.F90:965

.. keyword:: VEC
   :section: &RESP

   Option of :kw:`INTERACTION`: rotation axis of the Wannier functions for
   each molecule; :kw:`NMOL` lines with three reals each are read. Detected
   as a substring, so any line containing ``VEC`` triggers it in that
   context.

.. source: src/respin_p_utils.mod.F90:894

.. keyword:: VOA
   :section: &RESP
   :options: {AT, MD, POSITION_FORM, AMATRIX, CURRENT, DENSITY, VERBOSE, DEBUG}

   Vibrational optical activity from nuclear velocity perturbation theory.
   Exactly one of the options ``AT`` and ``MD`` must be given on the same
   line:

   * ``AT``: atomic axial and polar tensors (AAT/APT) in velocity form for
     all atoms (or those selected with :kw:`ATOMLIST`).
   * ``MD``: dipole moments along a molecular dynamics run.

   Further options: ``POSITION_FORM`` (also the APT in position form, only
   with ``AT``), ``AMATRIX`` (A-matrix, only with ``AT``), ``CURRENT``
   (electronic probability current density, only with ``MD``), ``DENSITY``
   (write ``DENSITY-<step>.cube`` files, only with ``CURRENT``),
   ``VERBOSE`` and ``DEBUG`` (implies ``VERBOSE``). Only cells with
   ``IBRAV`` 0, 1 or 8 are accepted. The keyword is recognized as ``VOA``
   followed by a blank; the options must be preceded by a blank.

   Example::

      VOA AT POSITION_FORM

   .. todo::

      MISSING: description of the output files; literature: Scherrer,
      Sebastiani, Vuilleumier, J. Chem. Theory Comput. 9, 5305 (2013) is
      printed by the code.

.. source: src/respin_p_utils.mod.F90:749

.. keyword:: WCPREDICT
   :section: &RESP

   Option of :kw:`INTERACTION`, for Born-Oppenheimer dynamics only: switches
   on the prediction of the Wannier centers. The number of MD steps between
   localization cycles is read from the next line. Only if preceded by a
   blank.

.. source: src/respin_p_utils.mod.F90:979

.. keyword:: WGAUSSIAN
   :section: &RESP

   Option of :kw:`INTERACTION`: use a Gaussian representation for the input
   and output of Wannier functions (with :kw:`WLOAD` / :kw:`WSAVE`).

.. source: src/respin_p_utils.mod.F90:955

.. keyword:: WLOAD
   :section: &RESP

   Option of :kw:`INTERACTION`: read the Wannier functions from file. Implies
   that no ground-state optimization is done and no restart file is read.
   Mutually exclusive with :kw:`WSAVE` (the last one given wins).

   .. todo::

      MISSING: file name(s).

.. source: src/respin_p_utils.mod.F90:938

.. keyword:: WSAVE
   :section: &RESP

   Option of :kw:`INTERACTION`: save the Wannier functions to a file.
   Mutually exclusive with :kw:`WLOAD` (the last one given wins).

   .. todo::

      MISSING: file name(s).

.. source: src/respin_p_utils.mod.F90:944
