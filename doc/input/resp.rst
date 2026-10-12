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

.. keyword:: NMR
   :section: &RESP
   :notes: options

   Calculate the NMR chemical shielding tensors for the system. Most important
   option: FULL, does a calculation with improved accuracy for periodic systems
   but takes a lot of time. Isolated systems: Use OVERLAP and 0.1 (on next line)
   for the same effect. *Be careful for non-hydrogen nuclei.* The shielding is
   calculated without contribution from the core electrons. Contact
   sebastia@mpip-mainz.mpg.de for further details.

.. keyword:: NOOPT
   :section: &RESP

   Do not perform a ground state wfn optimization. Be sure the restarted wfn is at
   the BO-surface.

.. keyword:: OACP
   :section: &RESP
   :alternatives: [DENSITY, REF_DENSITY, FORCE]

   Not documented.

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
