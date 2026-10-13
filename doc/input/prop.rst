.. _input-prop:

``&PROP`` section
=================

Calculation of properties. This section is only fully evaluated if the :kw:`PROPERTIES` keyword is given in the ``&CPMD`` section.

The section starts with ``&PROP`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &PROP

Keywords
--------

.. keyword:: AVERAGED POTENTIAL
   :section: &PROP

   Calculate averaged electrostatic potential in spheres of radius Rcut around the
   atomic positions.

   Parameter Rcut is read in from next line.

.. keyword:: CHARGES
   :section: &PROP
   :alternatives: [PARA, HYPER]

   Calculate atomic charges. Charges are calculated according to the method of
   Hirshfeld [Hirshfeld77]_ and charges derived from the
   electrostatic potential [Cox84]_.

   The ESP charges are obtained by a least-squares fit to the electrostatic
   potential with the constraint of the total charge. Two restrained
   variants are available:

   - ``PARA``: harmonic restraint towards reference charges, which are read
     from the file ``CZERO`` (one value per atom, in input order) if it
     exists, otherwise zero. A reference value of 100 or more restrains the
     corresponding charge strongly towards zero.
   - ``HYPER``: hyperbolic restraint towards zero
     (:math:`\sqrt{q^2+b^2}` with :math:`b=0.1`), solved iteratively (at most
     5 iterations, until the relative RMS deviation changes by less than
     :math:`10^{-4}`).

   The restraint weight is :math:`5.8\cdot 10^{-6}` times the number of
   grid points used in the fit.

   .. todo::

      MISSING: literature reference for the restrained fits (RESP-like).

.. keyword:: CONDUCTIVITY
   :section: &PROP

   Computes the optical conductivity according to the Kubo-Greenwod formula

   .. math::
      :label: condu

      \sigma(\omega) =
           \frac{2 \pi e^2}{3m^2 V_{\rm cell}} \frac{1}{\omega }
           \sum_{i,j} (f_i-f_j)
           |\langle \psi _i| \hat{\bf p} |\psi _j \rangle |^2
           \delta(\epsilon _i -\epsilon_j - \hbar \omega)

   where :math:`\psi _i` are the Kohn-Sham eigenstates, :math:`\epsilon _i` their
   corresponding eigenvalues, :math:`f_i` the occupation number and the difference
   :math:`f_i-f_j` takes care of the fermionic occupancy. This calculation is
   executed when the keyword PROPERTIES is used in the section &CPMD ... &END. In
   the section &PROP ... &END the keyword CONDUCTIVITY must be present; the
   next line contains an integer ``n`` and the optional entries ``STEP=s``
   and ``N=m``. ``n`` is the interval (in MD steps) at which the
   conductivity is computed when the keyword is used in an MD run (it is
   ignored in a PROPERTIES run). ``s`` is the bin width of the
   :math:`\sigma(\omega)` histogram in eV (**default** 0.5) and ``m`` the
   number of bins (**default** 40); non-positive values reset the defaults.
   Example::

      CONDUCTIVITY
        1 STEP=0.14 N=100

   Note that, since this is a "PROPERTIES" calculation,
   *you must have previously computed the electronic structure of your system and
   have a consistent* :kw:`RESTART` *file ready to use*.
   A file MATRIX.DAT is written in your working directory,
   where all the non-zero transition amplitudes and related information are
   reported (see the header of MATRIX.DAT). An example of application is given in
   Refs. [solve]_ [solve2]_.

   .. todo::

      MISSING: confirm that ``n`` has no meaning in PROPERTIES runs.

.. keyword:: CORE SPECTRA
   :section: &PROP

   Computes the X-ray adsorption spectrum and related transition matrix elements
   according to Ref. [xray]_. This calculation is executed when the
   keyword PROPERTIES is used in the section &CPMD ... &END. In the section &PROP
   ... &END the keyword CORE SPECTRA must be present and the core atom number
   (e.g. 10 if it is the 10\ :math:`th` atom in your list) and core level energy
   (in au) are read from the next line, while in the following line the :math:`n`
   and :math:`l` quantum numbers of the selected core level, along with the
   exponential factor :math:`a` of the STO orbital for the core level must be
   provided. In the case of :math:`1s` states, the core orbital is reconstructed
   as

   .. math::
      :label: 1s

      \psi _{1s}(r) = 2 a^{\frac{3}{2}} r \cdot \exp (-a\cdot r)

   and it is this :math:`a` value in au that must be supplied in input. As a
   general rule, first-row elements in the neutral case have the following
   :math:`a` values: B (4.64), C (5.63), N (6.62), O (7.62). For an excited atom
   these values would be of course a bit larger; e.g. for O it is 7.74453, i.e.
   1.6 % larger. Since this is a "PROPERTIES" calculation, *you must have
   previously computed the electronic structure of your system and have a
   consistent* :kw:`RESTART` *file ready to use*. A file XRAYSPEC.DAT
   is written in your working directory, containing all the square transition
   amplitudes and related information, part of which are also written in the
   standard output. Warning: in order to use this keyword you need special
   pseudopotentials. These are provided, at least for some elements, in the PP
   library of CPMD and are named as \*_HOLE.psp

.. keyword:: CUBECENTER
   :section: &PROP

   Sets the center of the cubefiles produced by the :kw:`CUBEFILE`
   flag. The next line has to contain the coordinates of the center in Bohr or
   Angstrom, depending on whether the :kw:`ANGSTROM` keyword was
   given. **Default** is the geometric center of the system.

.. keyword:: CUBEFILE
   :section: &PROP
   :options: {ORBITALS, DENSITY, POTENTIAL, HALFMESH}

   Plots the requested objects in Gaussian cube format. Several of the
   options can be combined on one line.

   - ``DENSITY``: total density to ``RHO_TOT.cube``; with :kw:`LSD` also the
     spin density to ``RHO_SPIN.cube``.
   - ``ORBITALS`` (aliases ``ORBITAL`` and ``WAVEFUNCTION``): the number of
     orbitals is read from the next line and their indices from the line
     after; one file ``PSI.<n>.cube`` per orbital, or ``PSI_A.<n>.cube`` and
     ``PSI_B.<n>.cube`` with :kw:`LSD`. The orbitals are first rotated to
     canonical (Kohn-Sham) form.
   - ``POTENTIAL``: the potential computed by ``rscpot`` from the current
     density is written to ``V_Hxc.cube``.
   - ``HALFMESH``: halves the number of grid points per direction (file size
     reduced by a factor of 8).

   The centre of the cube can be set with :kw:`CUBECENTER`.

   .. todo::

      MISSING: exact content of ``V_Hxc.cube`` (the file name suggests the
      Hartree plus exchange-correlation potential), and the orbital
      numbering (index ``i`` selects the stored state ``nstate-i+1``).

.. keyword:: DIPOLE MOMENT
   :section: &PROP
   :alternatives: [BERRY, RS]

   Calculate the dipole moment.

   Without additional keywords the dipole is obtained by integrating the
   charge density in reciprocal space over a sphere around the centre of
   the core charge; this is only implemented for simple cubic, fcc and
   orthorhombic (SYMMETRY 8) supercells. The keyword **RS** requests the
   real-space algorithm (a warning is printed for periodic systems). The
   keyword **BERRY** requests the Berry phase algorithm. ``DIPOLE BERRY``
   (without MOMENT) is equivalent to ``DIPOLE MOMENT BERRY``. If
   ``EXCITED`` is also present on the line, :kw:`EXCITED DIPOLE` is
   activated in addition.

   .. todo::

      MISSING: confirm the identification of SYMMETRY 8 as orthorhombic for
      the integration radius :math:`0.5\,\Omega^{1/3}`.

.. keyword:: EFG
   :section: &PROP

   Calculates electric field gradient tensors at the nuclei. The total
   valence density (stored after the SCF) plus the Gaussian ionic charges,
   minus the ionic charge of the atom itself, is converted into a potential
   via :math:`4\pi/G^2`, and its second derivatives at each nucleus are
   evaluated in reciprocal space; the tensor is made traceless. The
   elements (1,1) (2,2) (3,3) (1,2) (1,3) (2,3) and the electric field
   are appended to the file ``EFG`` and the tensor eigenvalues to
   ``EV_EFG``. The calculation is performed at the end of a wavefunction
   optimization and at every step of a Car-Parrinello MD run (not with
   broken-symmetry runs); it is not part of a :kw:`PROPERTIES` run.

   .. todo::

      MISSING: units of the printed tensors and fields, literature
      reference, and whether core (PAW-like) corrections are needed for
      quantitative results.

.. source: src/propin_utils.mod.F90:540

.. keyword:: EPR
   :target: EPR PROP
   :section: &PROP

   Calculates the anisotropic spin-dipole coupling tensors at the nuclei
   from the spin density :math:`\rho_\alpha-\rho_\beta`. The keyword is
   only active together with :kw:`LSD`; otherwise it is ignored silently.
   The spin density (stored after the SCF) is converted into a potential
   via :math:`4\pi/G^2`, and its second derivatives at each nucleus are
   evaluated in reciprocal space. The elements (1,1) (2,2) (3,3) (1,2)
   (1,3) (2,3) and the spin density at the nucleus are appended to the
   file ``EPR`` and the tensor eigenvalues to ``EV_EPR``. The calculation
   is performed at the end of a wavefunction optimization and at every step
   of a Car-Parrinello MD run (not with broken-symmetry runs); it is not
   part of a :kw:`PROPERTIES` run. Not to be confused with the EPR option
   of the ``&RESP`` section.

   .. todo::

      MISSING: units of the printed tensors, conversion to hyperfine
      coupling constants, and literature reference.

.. source: src/propin_utils.mod.F90:537

.. keyword:: EXCITED DIPOLE
   :section: &PROP

   Calculate the difference of dipole moments between the ground state density and
   a density generated by differently occupied Kohn-Sham orbitals.

   On the next line the number of dipole moments to calculate and the total number
   orbitals has to be given. On the following lines the occupation of the states
   for each calculation has to be given. By default the dipoles are calculated by
   the method used for the **DIPOLE MOMENT** option and the same restrictions
   apply. If the **LOCAL DIPOLE** option is specified the dipole moment
   differences are calculated within the same boxes.

   ``EXCITED`` can be combined with :kw:`DIPOLE MOMENT` or
   :kw:`LOCAL DIPOLE` on the same line.

.. keyword:: G_SPACE_LOC
   :section: &PROP
   :alternatives: [ALL, PART, LIST]

   Localizes the orbitals in reciprocal space (routines ``g_loc_dr*``)
   within a :kw:`PROPERTIES` run. The states to be localized are chosen by
   the option:

   - ``ALL`` (**default**, also used when no option is given): all states;
   - ``PART``: the first and last state index of a contiguous range are
     read from the next line;
   - ``LIST``: the number of states is read from the next line and the
     state indices from the following line(s).

   The localization functional and its parameters are controlled by
   ``&CPMD`` keywords such as ``GLOCALIZATION PARAMETERS``,
   ``GFUNCTIONAL TYPE``, ``G_REAL``/``G_COMPLEX`` and ``GLOC WFNOUT``.

   Example::

      G_SPACE_LOC PART
        5 12

   .. todo::

      MISSING: description of the localization functional and literature
      reference; output files produced.

.. source: src/propin_utils.mod.F90:543

.. keyword:: LDOS
   :section: &PROP

   Calculate the layer projected density of states. The number of layers
   and the direction perpendicular to the layers (1, 2 or 3 for the first,
   second or third real-space grid direction) are read from the next line.

   To use the LDOS keyword, the user must first have performed a wavefunction
   optimization and then restart with the :kw:`PROPERTIES` and
   :kw:`LANCZOS DIAGONALIZATION` keywords in the &CPMD
   section (and LDOS in the &PROP section).

   **Warning:** If you use special k-points for a special structure you need to
   symmetrize charge density for which you must specify the :kw:`POINT GROUP`.

.. keyword:: LOCAL DIPOLE
   :section: &PROP

   Calculate :math:`numloc` local dipole moments.

   | :math:`numloc` is read from the next line followed by two numloc lines with
     the format:
   | :math:`xmin` :math:`ymin` :math:`zmin`
   | :math:`xmax` :math:`ymax` :math:`zmax`

   Any line containing both ``DIPOLE`` and ``LOCAL`` (without ``MOMENT`` or
   ``BERRY``) is accepted. The six box limits of each region are read with
   one list-directed READ, so they may be given on one or two lines. They
   are in Bohr, or in Angstrom if :kw:`ANGSTROM` is set in ``&SYSTEM``. At
   most 20 regions are allowed.

.. keyword:: LOCALIZE
   :section: &PROP
   :options: {WAN_SPREAD}

   Localize the molecular orbitals [Hutter94b]_ as defined through
   the atomic basis set.

   The same localization transformation is then applied also to the wavefunctions
   in the plane wave basis. These wavefunction can be printed with the keyword
   **RHOOUT** specified in the section &CPMD section.

   With ``WAN_SPREAD`` a spread threshold is read from the next line. It is
   used by :kw:`WANNIER DOS` (``&CPMD``): the Wannier-function Hamiltonian is
   reordered by spread and written to ``WANNIER_HAM_SPREAD_ORDERED``, and the
   states with spread below and above the threshold are diagonalized
   separately. A value that is zero or negative additionally prepares a
   harmonic analysis (``PREPARATION FOR HARMONIC ANALYSIS``) and its
   absolute value is used as threshold. The run stops if no state has a
   larger spread.

   .. todo::

      MISSING: units of the threshold (spread as stored in the Wannier
      centres array) and the purpose of the harmonic analysis.

.. keyword:: n-CENTER CUTOFF
   :section: &PROP

   The cutoff for printing the 3-center (``3-CENTER CUTOFF``) or 4-center
   (``4-CENTER CUTOFF``) shared electron numbers of the Davidson
   :kw:`POPULATION ANALYSIS` is read from the next line. **Default** is 0.05
   for 3-center and 0.01 for 4-center terms. All one- and two-center terms
   are printed; ``2-CENTER CUTOFF`` is not recognized.

.. keyword:: NOPRINT ORBITALS
   :section: &PROP

   Do not print the wavefunctions in the atomic basis set. The singular form
   ``ORBITAL`` is also accepted.

.. keyword:: OPTIMIZE SLATER EXPONENTS
   :section: &PROP

   Sets an internal flag (``prop1%opts``) that is not used anywhere in the
   code; the keyword has no effect. Any line containing ``OPTIMIZE`` and
   ``SLATER`` is accepted.

.. keyword:: POLARIZABILITY
   :section: &PROP

   Computes the polarisability of a system, intended as dipole moment per unit
   volume. An integer ``n`` is read from the next line: during MD the
   polarisability is computed every ``n`` steps (and at the last step).
   The spellings ``POLARIZABILITIES``, ``POLARISABILITY`` and
   ``POLARISABILITIES`` are also accepted.

.. keyword:: POPULATION ANALYSIS
   :section: &PROP
   :alternatives: [MULLIKEN, DAVIDSON],[2-CENTER, 3-CENTER, 4-CENTER]

   The type of population analysis that is performed with the projected
   wavefunctions.

   Löwdin charges are given with both options. For the Davidson
   analysis [Davidson67]_ the maximum complexity is selected with one of
   the words ``2-CENTER`` (default), ``3-CENTER`` or ``4-CENTER``. The
   keyword implies :kw:`PROJECT WAVEFUNCTION`.

   For the Davidson option one has to specify the number of atomic orbitals
   that are used in the analysis, one value per species in ``&ATOMS``
   order. They are read with one list-directed READ, so they may be given
   on one or several lines. An input example
   for a water molecule is given in the hints section :ref:`hints:pop`.

.. keyword:: PROJECT HARMONICS
   :section: &PROP
   :options: {RADIAL, SPREAD}

   Projects the Kohn-Sham orbitals (after orthogonalization) on real
   spherical harmonics centred at a given point. The next line contains the
   number of spherical harmonics ``ny`` and the three coordinates of the
   centre (used without unit conversion, i.e. in Bohr); the maximum angular
   momentum is :math:`l_{max}=\mathrm{INT}(\sqrt{ny-1})`, so ``ny`` should be
   :math:`(l_{max}+1)^2`. The overlap coefficients are written to the file
   ``YLMCOEF``.

   With ``RADIAL`` a further line with ``rmax nmax`` is read: the radius of
   integration and the number of equidistant radial grid points (increased
   by one if odd); the radial parts are written to the unformatted file
   ``RADIAL``. With ``SPREAD`` an additional line ``spr_min spr_max`` is
   read; it only takes effect together with ``RADIAL``: the spreads are read
   from the file ``SPREAD`` (fifth column) and only states whose spread lies
   strictly between ``spr_min`` and ``spr_max`` are analysed (without
   ``SPREAD`` the window -1.0 to 10.0 is used with all spreads set to 0).

   Example::

      PROJECT HARMONICS RADIAL
        16 0.0 0.0 0.0
        8.0 100

   .. todo::

      MISSING: the physical purpose of the analysis, the format of
      ``YLMCOEF`` and ``RADIAL``, and the expected format of ``SPREAD``
      (written by which keyword).

.. source: src/propin_utils.mod.F90:243

.. keyword:: PROJECT WAVEFUNCTION
   :section: &PROP

   The wavefunctions are projected on atomic orbitals.

   The projected wavefunctions are then used to calculate atomic populations and
   bond orders. The atomic orbitals to project on are taken from the &BASIS
   section. If there is no &BASIS section in the input a minimal Slater basis is
   used. See section :ref:`input:basis` for more details.

   The plural form ``PROJECT WAVEFUNCTIONS`` is also accepted. Note that
   :kw:`POPULATION ANALYSIS` switches this option on automatically, and that
   PROJECT WAVEFUNCTION cannot be combined with :kw:`LOCALIZE` (in ``&PROP``)
   or Wannier localization (the run stops with ``PROJECT WAVEFUNCTION AND
   LOCALIZATION ARE INCOMPATIBLE``).

.. keyword:: TRANSITION MOMENT
   :section: &PROP

   Calculate the dipole transition matrix element.

   On the following lines, the number of transitions and the involved orbitals are
   given. Example:

   ::

      TRANSITION MOMENT
      2
      6 7
      6 8


   This calculates the dipole transition matrix elements between KS states 6 and
   7, and between 6 and 8.
