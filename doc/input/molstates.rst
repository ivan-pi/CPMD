.. _input-molstates:

``&MOLSTATES`` section
======================

The ``&MOLSTATES`` section defines a partitioning of the system into
molecules. It is read only if :kw:`MOLECULAR STATES` is given in the
``&TDDFT`` section. The occupied Kohn-Sham orbitals are localized (Wannier
functions), each Wannier function is assigned to the molecule of the atom
closest to its center, and the Kohn-Sham Hamiltonian is diagonalized within
the block of each molecule ("molecular states"). The molecular eigenvalues
are printed and, together with :kw:`TD_METHOD_A`, the highest eigenvalue of
each molecule is used as reference energy for the statistical potential
(e.g. SAOP) of that molecule in an :kw:`ELECTRONIC SPECTRA` calculation.
Not implemented for :kw:`LSD` and QM/MM, and not available for TDDFT forces
(the program stops).

Syntax::

   &MOLSTATES
     MOLECULES
       m(1) m(2) ... m(NAT)
     STATES
       n(1) n(2) ... n(NMOL)
   &END

:kw:`MOLECULES` must be given before :kw:`STATES <STATES MOLSTATES>`.

Summary
-------

.. keyword-summary:: &MOLSTATES

Keywords
--------

.. keyword:: MOLECULES
   :section: &MOLSTATES

   Assigns every atom to a molecule. The next line(s) contain one integer
   molecule index per atom, in the order of the atoms in the ``&ATOMS``
   section (list-directed read, so the values may span several lines). The
   number of molecules is the largest index given; indices should run from 1
   to that number.

   Example (two water molecules, atoms ordered O O H H H H)::

      MOLECULES
       1 2 1 1 2 2

.. source: src/lr_in_utils.mod.F90:1216

.. keyword:: STATES
   :target: STATES MOLSTATES
   :section: &MOLSTATES

   Number of occupied orbitals belonging to each molecule, one integer per
   molecule, read from the next line(s). Must follow :kw:`MOLECULES` (the
   program stops with "Molecules not defined" otherwise). If the number of
   Wannier functions assigned to a molecule differs from the value given, a
   warning "NUMBER OF STATES IS INCONSISTENT" is printed.

   Example::

      STATES
       4 4

.. source: src/lr_in_utils.mod.F90:1226
