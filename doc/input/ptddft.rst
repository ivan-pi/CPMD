.. _input-ptddft:

``&PTDDFT`` section
===================

Propagation TDDFT for Ehrenfest dynamics and spectra. This section contains
specific information for Ehrenfest dynamics and for spectra computed by the
propagation of the perturbed Kohn-Sham orbitals (Fourier transform of the
induced dipole fluctuation) [taver_eh]_ [taver1]_.

The section starts with ``&PTDDFT`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &PTDDFT

Keywords
--------

.. keyword:: ACCURACY
   :section: &PTDDFT

   Specifies the accuracy to be reached in the Cayley propagation scheme used in
   Ehrenfest type of dynamics and spectra calculation.

.. keyword:: ADD_MASK
   :section: &PTDDFT

   Applies an absorbing-boundary mask to the propagated Kohn-Sham orbitals
   after each Cayley propagation step (with :kw:`RUNGE_KUTTA` only in the
   second half step). Three real numbers are read from the next line:
   ``p1 p2 p3``. Along each Cartesian direction the mask is 1 inside
   :math:`|x| < p_1 L`, decays as
   :math:`\cos^2\left(\frac{\pi}{2}\frac{|x|-p_1 L}{(p_2-p_1)L}\right)` for
   :math:`p_1 L \le |x| < p_2 L` and is 0 beyond :math:`p_2 L`, where
   :math:`L` is the box length in that direction and :math:`x` is measured
   from the mask center. If ``p1`` equals ``p2`` the mask is a sharp cut-off.
   The third value is a scaling factor for both boundaries that is only
   applied in a momentum-distribution mode of the mask routine (comments in
   the code suggest 0.42, 0.48 and 0.8 as typical values).

   Using :kw:`ADD_MASK` automatically switches on :kw:`DO_NOT_NORMALIZE`, so
   that the norm lost through the boundary is not restored.

   Example::

      ADD_MASK
       0.42 0.48 0.8

   .. todo::

      MISSING: which centering (center of charge of the ions or center of the
      box) is used in production runs, and when the momentum-distribution
      mode (third parameter) is active; recommended values.

.. source: src/read_prop_utils.mod.F90:224

.. keyword:: BERRY
   :section: &PTDDFT

   Sets an internal flag (``tberry``). No further input is read.

   In the current code the flag is only broadcast and is not used anywhere
   else, i.e. the keyword has no effect.

   .. todo::

      MISSING: intended purpose (presumably a Berry-phase dipole evaluation
      during the propagation); the feature is not implemented.

.. source: src/read_prop_utils.mod.F90:216

.. keyword:: BOHMIANTRAJ
   :section: &PTDDFT

   Sets the flag ``bohm%tbohmtraj`` and reads the number of Bohmian
   trajectories (integer) from the next line. In an Ehrenfest run
   (:kw:`MOLECULAR DYNAMICS` EH) combined with :kw:`PLANE_WAVE`, the number
   must be a perfect square, otherwise the program stops.

   Apart from this consistency check the values are not used outside the
   input parser in the current code, i.e. the keyword has no further effect.

   Example::

      BOHMIANTRAJ
       16

   .. todo::

      MISSING: meaning and use of the Bohmian trajectories; the feature does
      not appear to be implemented in this version.

.. source: src/read_prop_utils.mod.F90:239

.. keyword:: BOHM_R
   :section: &PTDDFT

   Sets the flag ``bohm%restart_bohm``. No further input is read.

   In the current code the flag is not used outside the input parser, i.e.
   the keyword has no effect.

   .. todo::

      MISSING: intended purpose (presumably restart of Bohmian trajectories,
      cf. :kw:`BOHMIANTRAJ`); the feature is not implemented in this version.

.. source: src/read_prop_utils.mod.F90:243

.. keyword:: CAL_CURRENT
   :section: &PTDDFT

   Computes the electronic current density during the Cayley propagation
   and writes its Cartesian components as cube files named
   ``currx.cube_<step>``, ``curry.cube_<step>``, ... An integer is read
   from the next line (``curr_step``, **Default** is 100).

   In the current code the value read is not used: the current is written
   at every call of the propagator.

   Example::

      CAL_CURRENT
       100

   .. todo::

      MISSING: intended meaning of the integer (presumably the output stride
      in steps), which is currently ignored.

.. source: src/read_prop_utils.mod.F90:219

.. keyword:: DENS_R
   :section: &PTDDFT

   Sets the flag ``bohm%eh_restart_dens`` and reads an integer
   (``bohm%eh_restart_n``) from the next line.

   In the current code neither value is used outside the input parser, i.e.
   the keyword has no effect.

   .. todo::

      MISSING: intended purpose (presumably a density restart of Ehrenfest
      dynamics); the feature is not implemented in this version.

.. source: src/read_prop_utils.mod.F90:246

.. keyword:: DO_NOT_NORMALIZE
   :section: &PTDDFT

   Do not renormalize the propagated Kohn-Sham orbitals. By default the
   orbitals are renormalized after each Cayley propagation step and in the
   orthogonalization step of Ehrenfest dynamics / propagation spectra. With
   this keyword the norm of each orbital is allowed to change (e.g. when
   charge is absorbed). No further input is read.

   It is switched on automatically by :kw:`ADD_MASK` and
   :kw:`IONIZE_FROM_RDM`.

.. source: src/read_prop_utils.mod.F90:263

.. keyword:: EXT_POTENTIAL
   :target: EXT-POTENTIAL
   :section: &PTDDFT

   An external potential is applied to the KS orbitals before propagation. This
   keyword is used with the principal keyword **MOLECULAR DYNAMICS EH** or
   :kw:`PROPAGATION SPECTRA`. The type of the perturbation
   is specified with the keyword **PERT_TYPE** (in &PTDDFT) and its intensity with
   **PERT_AMPLI** (in &PTDDFT).

.. keyword:: EXT_PULSE
   :target: EXT-PULSE
   :section: &PTDDFT

   A Dirac type pulse is applied to the KS orbitals before propagation. This
   keyword is used with the principal keyword :kw:`PROPAGATION SPECTRA`. The intensity of the perturbing field is
   read from the next line.

.. keyword:: FIX_SPIN_DENS
   :section: &PTDDFT

   Constrains the spin density during the Cayley propagation "assuming a
   doublet" (code comment): after each step the beta orbitals are paired with
   the alpha orbitals in order, and both orbitals of a pair are replaced by
   their average. The pairing stops at the alpha orbital whose index equals
   the state given with :kw:`IONIZE_FROM_STATE` (as coded, all remaining
   beta orbitals are then left unchanged). No further input is read. The
   keyword also sets the start from real wavefunctions (``cntl%start_real``).

   .. todo::

      MISSING: physical justification, intended use cases and whether it
      requires :kw:`LSD`; the pairing of alpha and beta orbitals should be
      checked by an expert.

.. source: src/read_prop_utils.mod.F90:259

.. keyword:: GAUGEFIELD PARA
   :section: &PTDDFT

   Intended to read three parameters of the laser pulse used with the gauge
   field (``GAUGEFIELD`` in ``&CPMD``) from the next line:
   ``gswtt gdfreq gimax`` (**Defaults** 250, 0.00025 and 0). In the pulse
   :math:`E(t) = E_0 \sin^2(\pi t/T)\sin(\omega t)` the code uses
   ``gdfreq`` as :math:`1/T`; in the circularly polarized variant ``gswtt``
   is the width and ``gdfreq`` the center of a Gaussian envelope. ``gimax``
   is not used.

   **Warning:** the keyword line in ``&PTDDFT`` is read list-directed, so
   only the first word (``GAUGEFIELD``) is kept and the two-word pattern
   ``GAUGEFIELD PARA`` can never be matched. The keyword is therefore
   currently ineffective unless the whole line is enclosed in quotes
   (``'GAUGEFIELD PARA'``).

   .. todo::

      MISSING: units and recommended values; confirm the intended keyword
      spelling (and fix the parser).

.. source: src/read_prop_utils.mod.F90:200

.. keyword:: GAUGEPULSE PARA
   :section: &PTDDFT

   Intended to read the number of periods of the gauge-field pulse
   (``gaugefield%nperiod``, **Default** 5) from the next line. The value is
   only referenced in commented-out code and has no effect.

   **Warning:** as for :kw:`GAUGEFIELD PARA`, only the first word of the
   keyword line is kept by the parser, so this two-word keyword can never be
   matched (unless the line is enclosed in quotes).

   .. todo::

      MISSING: intended purpose; the feature is not implemented.

.. source: src/read_prop_utils.mod.F90:204

.. keyword:: GAUGE_VECTOR
   :section: &PTDDFT

   Sets the flag ``gaugefield%tkgauge`` and reads an integer
   (``gaugefield%kgauge``) from the next line.

   The integer is not used. If the gauge-field propagation is used
   (``GAUGEPULSE`` / ``GAUGEFIELD`` in ``&CPMD``) the program stops with
   the message "tkgauge not implemented".

   .. todo::

      MISSING: intended purpose; the feature is not implemented.

.. source: src/read_prop_utils.mod.F90:208

.. keyword:: GAUSSIAN_WP
   :section: &PTDDFT

   Sets the flag ``bohm%tgaussian_wp`` and reads three lines:
   two reals (``gw_par1 gw_par2``), a position (``gw_r1 gw_r2 gw_r3``) and a
   momentum (``gw_p1 gw_p2 gw_p3``).

   In the current code these values are not used outside the input parser,
   i.e. the keyword has no effect.

   .. todo::

      MISSING: intended purpose (presumably an initial Gaussian wave packet)
      and meaning of the parameters; the feature is not implemented in this
      version.

.. source: src/read_prop_utils.mod.F90:228

.. keyword:: IONIZE_FROM_RDM
   :section: &PTDDFT

   Start the propagation from an ionized state constructed from the
   one-particle reduced density matrix (``restart_from_ionic_state``) of the
   orbitals read from the restart file. Requires :kw:`LSD` (the program
   stops otherwise). No further input is read. The keyword sets the start
   from real wavefunctions (``cntl%start_real``) and switches on
   :kw:`DO_NOT_NORMALIZE`.

   .. todo::

      MISSING: description of how the ionic state is constructed and of the
      required restart files.

.. source: src/read_prop_utils.mod.F90:255

.. keyword:: IONIZE_FROM_STATE
   :section: &PTDDFT

   Start the propagation from an ionized configuration: the orbitals are
   taken from the restart file and the orbital whose index (within the
   beta-spin orbitals) is read from the next line is removed. The keyword
   sets the start from real wavefunctions (``cntl%start_real``). The same
   index is excluded from the spin constraint of :kw:`FIX_SPIN_DENS`.

   Example::

      IONIZE_FROM_STATE
       1

   .. todo::

      MISSING: exact orbital indexing convention and occupation handling
      (the code also moves the first beta orbital to the end of the list);
      requirement of :kw:`LSD`.

.. source: src/read_prop_utils.mod.F90:250

.. keyword:: NUMBER_ZEROS
   :section: &PTDDFT

   Number of zeros used in the Chebyshev propagator (printed as "number of
   zeros") or the number of subintervals in the Cayley propagator (printed
   as "number of subintervalls"). The integer is read from the next line.
   **Default** is 6.

   Example::

      NUMBER_ZEROS
       6

.. source: src/read_prop_utils.mod.F90:143

.. keyword:: N_CYCLES
   :target: N-CYCLES
   :section: &PTDDFT

   Defines the number of cycles (time steps) used in the propagation of the
   perturbed KS orbitals in a :kw:`PROPAGATION SPECTRA`
   calculation. The number of cycles is read from the next line.

.. keyword:: PERT_AMPLI
   :target: PERT-AMPLI
   :section: &PTDDFT

   The amplitude of the perturbation used with :kw:`MOLECULAR DYNAMICS` EH or :kw:`PROPAGATION SPECTRA` is read from the next line.

.. keyword:: PERT_DIRECTION
   :target: PERT-DIRECTION
   :section: &PTDDFT

   Use with MOLECULAR DYNAMICS EH or PROPAGATION SPECTRA to specify the direction
   (x=1,y=2,z=3) of the perturbing field. The direction code (1,2,3) is read from
   the next line.

.. keyword:: PERT_TYPE
   :target: PERT-TYPE
   :section: &PTDDFT

   The type of the static perturbation used with :kw:`MOLECULAR DYNAMICS` EH or :kw:`PROPAGATION SPECTRA` is read from the next line. 1: Dirac pulse,
   2: Heaviside step function (constant value after :math:`t=0`).

.. keyword:: PIPULSE
   :section: &PTDDFT

   Specifies a time dependent pi-pulse to be used with MOLECULAR DYNAMICS EH. Use
   PIPULSE together with TD_POTENTIAL. The pulse strength is read from the next
   line (see subroutine gaugepot_laser in td_utils.mod.F90 for further details).

.. keyword:: PLANE_WAVE
   :section: &PTDDFT

   Sets the flag ``bohm%tplanewave_wp`` and reads two lines: a position
   (``gw_r1 gw_r2 gw_r3``) and a momentum (``gw_p1 gw_p2 gw_p3``).

   Apart from the perfect-square check of :kw:`BOHMIANTRAJ`, the values are
   not used outside the input parser in the current code.

   .. todo::

      MISSING: intended purpose (presumably a plane-wave initial wave packet);
      the feature is not implemented in this version.

.. source: src/read_prop_utils.mod.F90:234

.. keyword:: POINT_CHARGE
   :section: &PTDDFT

   Adds a time-dependent external potential generated by a point charge.
   Seven values are read from the next line:
   ``t0 dt q fx fy fz iwrite``. The potential :math:`-q/|\mathbf{r}-\mathbf{r}_q|`
   is switched on during the time window
   :math:`t_0 - \Delta t/2 < t < t_0 + \Delta t/2`. The position
   :math:`\mathbf{r}_q` is given in fractions ``fx fy fz`` of the cell
   edges (diagonal elements of the cell vectors, shifted by half a mesh
   spacing). If ``iwrite`` is 0, the initial density is written to the file
   ``DENSITY_REF`` for later analysis. **Defaults** are 0 for all values.

   This perturbation is mandatory for ``PROPAGATION DISTRUB`` (sic) in
   ``&CPMD`` (the program stops otherwise). The field enters the total
   potential but its interaction energy is not added to the total energy.

   Example::

      POINT_CHARGE
       5.0 2.0 1.0 0.5 0.5 0.8 0

   .. todo::

      MISSING: units of ``t0``/``dt`` (presumably atomic units of time) and of
      ``q``; meaning of other values of ``iwrite``.

.. source: src/read_prop_utils.mod.F90:173

.. keyword:: PPCORR
   :section: &PTDDFT

   Requests the gauge correction of the non-local pseudopotentials for the
   gauge-field propagation (the output shows "GAUGEING OF THE NON-LOCAL
   PSEUDOPOTENTIALS"). No further input is read.

   In the current code the correction term is computed only if an internal
   flag (``update_p``) is set, which is hard-coded to false, so the keyword
   has no effect on the results.

   .. todo::

      MISSING: confirm that the feature is disabled on purpose.

.. source: src/read_prop_utils.mod.F90:213

.. keyword:: PROP_TSTEP
   :target: PROP-TSTEP
   :section: &PTDDFT

   Propagation timestep used in Ehrenfest dynamics. It is used in the spectra
   calculation (**PROPAGATION SPECTRA** ) to specify the time step for the
   propagation of the KS orbitals.

.. keyword:: RESTFILE
   :target: RESTFILE PTDDFT
   :section: &PTDDFT

   Defines a restart code for the restart of the Ehrenfest dynamics (:kw:`MOLECULAR DYNAMICS` EH) and the propagation spectra
   (:kw:`PROPAGATION SPECTRA`). The restart option is read
   from the next line: 0(=default) restart from the (complex)wavefunctions in the
   file wavefunctions. This option is used in case of a continuation run; 1.
   restart from the orbital files WAVEFUNCTION.n, where :math:`n` is the index of
   the KS orbital and runs from :math:`1` to the number of states (This states a
   prepare in a previous run using the KOHN-SHAM ENERGIES principal keyword), 2;
   restart from the orbitals stored in RESTART (obtained from a optimization run
   with tight convergence (at least 1.D-7)).

.. keyword:: TD_POTENTIAL
   :target: TD-POTENTIAL
   :section: &PTDDFT

   Defines a time dependent external potential to be used in Ehrenfest dynamics
   (:kw:`MOLECULAR DYNAMICS` EH). Can be used with the
   keyword :kw:`PIPULSE`. The frequency of the external field is read
   from the next line (in atomic units). For more information see subroutine
   gaugepot_laser in td_utils.mod.F90 .
