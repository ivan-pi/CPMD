.. _input-atoms:

``&ATOMS`` section
==================

Atoms, pseudopotentials and related parameters (**required**). The pseudopotential input is explained in more detail in :ref:`S_Pseudopotentials`.

The section starts with ``&ATOMS`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &ATOMS

Keywords
--------

.. keyword:: ATOMIC CHARGES
   :section: &ATOMS

   Changes the default charge (0) of the atoms for the initial guess to the values
   read from the next line. One value per atomic species has to be given.

.. keyword:: CHANGE BONDS
   :section: &ATOMS

   The buildup of the empirical Hessian can be affected.

   You can either add or delete bonds. The number of changed bonds is read from
   the next line. This line is followed by the description of the bonds. The
   format is

   *{ ATOM1   ATOM2   FLAG}* .

   *ATOM1* and *ATOM2* are the numbers of the atoms involved in the bond. A *FLAG*
   of :math:`-1` causes a bond to be deleted and a *FLAG* of :math:`1` a bond to
   be added.

   Example:

   ::

      CHANGE BONDS
      2
      1 2 +1
      6 8 -1


.. keyword:: CONFINEMENT POTENTIAL
   :section: &ATOMS

   The use of this label activates a spherical Gaussian confinement potential in
   the calculation of the form factor of pseudopotentials. In the next line(s) two
   parameters for each atomic species must be supplied: the amplitude
   :math:`\alpha` and the cut off radius :math:`r_c`. The Gaussian spherical
   amplitude is computed as :math:`A=\pi ^{3/2}r_c^3\cdot \alpha` and the Gaussian
   confinement potential reads

   .. math::
      :label: pconf

      V({\bf G}) = \sum_{\bf G} A \cdot |{\bf G}|\cdot e^{-G^2r_c^2/4}

   being **G** the G-vectors, although in practice the loop runs only on the
   G-shells :math:`G=|{\bf G}|`.

.. keyword:: CONSTRAINTS ... END CONSTRAINTS
   :section: &ATOMS

   With this option you can specify several constraints and restraints on the
   atoms. (see section :ref:`sec:cnstr` for more information on the
   available options and the input format).

.. keyword:: DUMMY ATOMS
   :section: &ATOMS

   The definition of dummy atoms follows this keyword.

   Three different kinds of dummy atoms are implemented. Type 1 is fixed in space,
   type 2 lies at the arithmetic mean, type 3 at the center of mass of the
   coordinates of real atoms.

   The first line contains the total number of dummy atoms. The following lines
   start with the type label **TYPE1, TYPE2, TYPE3, TYPE4**.

   For type 1 dummy atoms the label is followed by the Cartesian coordinates.

   For type 2 and type 3 dummy atoms the first number specifies the total number
   of atoms involved in the definition of the dummy atom. Then the number of these
   atoms has to be specified on the same line. A negative number of atoms stands
   for all atoms. For type 4, the dummy atom is defined as a weighed average of
   coordinates of real atoms with user-supplied weights. This feature is useful
   e. g.in constrained dynamics, because allows to modify positions and weights of
   dummy atoms according to some relevant quantity such as forces on selected
   atoms. Example:

   ::

      DUMMY ATOMS
      3
      TYPE1 0.0 0.0 0.0
      TYPE2 2 1 4
      TYPE3 -1


   Note: Indices of dummy atoms always start with total-number-of-atoms plus 1. In
   the case of a Gromos-QM/MM interface simulations with dummy hydrogen atoms for
   capping, it is total-number-of-atoms plus number-of-dummy-hydrogens plus 1

.. keyword:: GENERATE COORDINATES
   :section: &ATOMS

   The number of generator atoms for each species are read from the next line.

   These atoms are used together with the point group information to generate all
   other atomic positions. The input still has to have entries for all atoms but
   their coordinates are overwritten. Also the total number of atoms per species
   has to be correct.

.. keyword:: ISOTOPE
   :section: &ATOMS

   Changes the default masses of the atoms.

   This keyword has to be followed by *NSP* lines (number of atom types). In each
   line the new mass (in a.m.u.) of the respective species has to be specified (in
   order of their definition).

.. keyword:: METADYNAMICS ... END METADYNAMICS
   :section: &ATOMS

   Initiate Metadynamics (see section :ref:`sec:meta` for more information
   on the available options and the input format).

.. keyword:: MOVIE TYPE
   :section: &ATOMS

   Assign special movie atom types to the species.

   The types are read from the next line. Values from 0 to 5 were allowed in the
   original MOVIE format.

.. keyword:: VELOCITIES ... END VELOCITIES
   :section: &ATOMS

   Sets an **initial velocity** for specified atoms.

   The first line contains first the total number of specified atomic velocities
   followed **on the same line** by the list of atomic numbers for which the
   velocities will be read. On each of the following lines the x, y and z
   coordinates of the velocities of an atom have to be specified. These values
   will ignored in case of starting with :kw:`RESTART` VELOCITIES..

   NOTE: these velocities are rescaled to produce the initial temperature as
   specified by :kw:`TEMPERATURE`. The default temperature,
   however, is 0K, so you **have** to set the matching temperature or your initial
   velocities will be useless.
