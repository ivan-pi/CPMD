.. _input-classic:

``&CLASSIC`` section
====================

Simple classical code interface.

The section starts with ``&CLASSIC`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &CLASSIC

Keywords
--------

.. keyword:: DEFORMATION
   :section: &CLASSIC
   :options: {SAMPLE}

   Periodically applies an affine scaling about a fixed centre to atomic
   coordinates during a classical-force-field MD run (atoms with the
   ``CLASSIC`` label in the ``&ATOMS`` section). Any line containing ``DEFORM`` is
   accepted.

   If ``SAMPLE`` is given on the same line, the number of MD steps
   ``ndefo`` between two deformations is read from the next line
   (**Default** is 1, i.e. every step). Then two further lines are read:

   - ``expx expy expz``: the scaling factors along :math:`x`, :math:`y`
     and :math:`z`;
   - ``xcen ycen zcen``: the centre of the scaling.

   At every step whose step number is a multiple of ``ndefo`` each
   coordinate is replaced by :math:`x \to x_{cen} + e_x (x - x_{cen})`
   (and likewise for :math:`y` and :math:`z`).

   Note that the routine ``deformation`` (``src/mdclas_utils.mod.F90``)
   only loops over classical atom types whose ``is_qm`` index is non-zero,
   i.e. over the atoms of the quantum species stored in the classical
   coordinate array.

   Example::

      DEFORMATION SAMPLE
        100
        1.001 1.001 1.0
        0.0 0.0 0.0

   .. todo::

      MISSING: units of the centre coordinates (no conversion is applied,
      so presumably Bohr), the physical purpose of the option, and whether
      restricting the deformation to the quantum species is intended.

.. source: src/cl_init_utils.mod.F90:125

.. keyword:: FORCE FIELD ... END FORCE FIELD
   :section: &CLASSIC

   Block with the classical force-field parameters, terminated by a line
   containing ``END`` and ``FORCE``. Each line starts with a tag:

   - ``R12 lab1 lab2 a rc``: purely repulsive :math:`(a/R)^{12}` pair
     potential between classical atom types with labels ``lab1`` and
     ``lab2`` (the ``LABEL=`` values of the ``CLASSIC`` species in
     ``&ATOMS``, at most 4 characters), parameter ``a`` and cutoff radius
     ``rc``. The parameters are stored symmetrically. Pairs in which both
     atoms are quantum atoms are skipped.
   - ``LJ``: recognized but ignored (not implemented).

   Lines with other tags are ignored.

   Example::

      FORCE FIELD
        R12 Ar Ar 6.0 15.0
      END FORCE FIELD

   .. todo::

      MISSING: units of ``a`` and ``rc`` (presumably Bohr) and the exact
      energy expression (``sigmasq`` = :math:`a^2` is used in
      ``clas_force``).

.. source: src/readff_utils.mod.F90:24

.. keyword:: FREEZE QUANTUM
   :section: &CLASSIC

   Freeze the quantum atoms and performs a classical MD on the others (in QMMM
   mode only !).

.. keyword:: FULL TRAJECTORY
   :section: &CLASSIC
   :options: {SAMPLE}

   Writes the coordinates and velocities of all classical atoms to the
   unformatted file ``CLASSIC`` (six single-precision numbers per atom) during
   a classical-force-field MD run. With ``SAMPLE`` the output interval ``n``
   in MD steps is read from the next line (**default** 1); a frame is written
   whenever :math:`\mathrm{mod}(n_{fi}-1, n) = 0`. Any line containing
   ``FULL`` and ``TRAJ`` is accepted.

.. source: src/cl_init_utils.mod.F90:107

.. keyword:: MOVIE
   :target: MOVIE CLASSIC
   :section: &CLASSIC

   Writes the coordinates of all classical atoms in movie format to the
   file ``CLASSIC_MOVIE`` during a classical-force-field MD run. The output
   interval ``imovc`` (in MD steps) is read from the next line; a frame is
   written whenever :math:`\mathrm{mod}(n_{fi}-1, imovc) = 0`. Each line of
   a frame contains the Cartesian coordinates in Angstrom, the atomic number
   (6 is used for atom types without an atomic number) and an integer type
   code (array ``imtyp``).
   Any line containing ``MOVIE`` is accepted.

   Example::

      MOVIE
        10

.. source: src/cl_init_utils.mod.F90:118

.. keyword:: PRINT COORDINATES
   :section: &CLASSIC

   Prints the coordinates of all classical atoms (index, label, QM species
   index, x, y, z) to the standard output at the start and at the end of a
   classical-force-field MD run. Any line containing ``PRINT`` and ``COORD``
   is accepted.

.. keyword:: PRINT FF
   :section: &CLASSIC

   Prints the parameters of the force field read in
   :kw:`FORCE FIELD ... END FORCE FIELD` (for the :math:`(a/R)^{12}` potential:
   the pair of atom labels, the parameter and the cutoff RC) to the standard
   output. Any line containing ``PRINT`` and ``FF`` is accepted.
