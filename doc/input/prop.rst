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

   Calculate atomic charges. Charges are calculated according to the method of
   Hirshfeld [Hirshfeld77]_ and charges derived from the
   electrostatic potential [Cox84]_.

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
   the section &PROP ... &END the keyword CONDUCTIVITY must be present and the
   interval interval :math:`\Delta \omega` for the calculation of the spectrum is
   read from the next line. Note that, since this is a "PROPERTIES" calculation,
   *you must have previously computed the electronic structure of your system and
   have a consistent* :kw:`RESTART` *file ready to use*. Further
   keyword: ``STEP=0.14``, where (e.g.) 0.14 is the bin width in eV of the
   :math:`\sigma(\omega)` histogram if you want it to be different from
   :math:`\Delta \omega`. A file MATRIX.DAT is written in your working directory,
   where all the non-zero transition amplitudes and related information are
   reported (see the header of MATRIX.DAT). An example of application is given in
   Refs. [solve]_ [solve2]_.

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
   :options: ORBITALS,DENSITY
   :alternatives: HALFMESH

   Plots the requested objects in .CUBE file format. If ORBITALS are demanded, the
   total number as well as the indices have to be given on the next and second
   next line. HALFMESH reduces the number of grid points per direction by 2, thus
   reducing the file size by a factor of 8.

.. keyword:: DIPOLE MOMENT
   :section: &PROP
   :alternatives: [BERRY, RS]

   Calculate the dipole moment.

   Without the additional keywords **BERRY** or **RS** this is only implemented
   for simple cubic and fcc supercells. The keyword **RS** requests the use of the
   real-space algorithm. The keyword **BERRY** requests the use of the Berry phase
   algorithm.

   **Default** is to use the real-space algorithm.

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

   Calculate the layer projected density of states. The number of layers is read
   from the next line.

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

.. keyword:: LOCALIZE
   :section: &PROP

   Localize the molecular orbitals [Hutter94b]_ as defined through
   the atomic basis set.

   The same localization transformation is then applied also to the wavefunctions
   in the plane wave basis. These wavefunction can be printed with the keyword
   **RHOOUT** specified in the section &CPMD section.

.. keyword:: n-CENTER CUTOFF
   :section: &PROP

   The cutoff for printing the n-center shared electron numbers is read from the
   next line. All one and two center terms are printed.

.. keyword:: NOPRINT ORBITALS
   :section: &PROP

   Do not print the wavefunctions in the atomic basis set.

.. keyword:: OPTIMIZE SLATER EXPONENTS
   :section: &PROP

   Not documented

.. keyword:: POLARIZABILITY
   :section: &PROP

   Computes the polarisability of a system, intended as dipole moment per unit
   volume.

.. keyword:: POPULATION ANALYSIS
   :section: &PROP
   :alternatives: [MULLIKEN, DAVIDSON],[n-CENTER]

   The type of population analysis that is performed with the projected
   wavefunctions.

   Löwdin charges are given with both options. For the Davidson
   analysis [Davidson67]_ the maximum complexity can be specified
   with the keyword **n-CENTER**.

   Default for n is 2, terms up to 4 are programmed. For the Davidson option one
   has to specify the number of atomic orbitals that are used in the analysis. For
   each species one has to give this number in a separate line. An input example
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
