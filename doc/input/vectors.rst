.. _input-vectors:

``&VECTORS`` section
====================

The ``&VECTORS ... &END`` section supplies user-defined search directions
(nuclear displacement vectors) for the iterative linear-response calculation
of non-adiabatic couplings. It is only read when the option ``SPECIFY`` is
given together with ``NVECT=n`` on the :kw:`COUPLINGS LINRES` line of the
``&SYSTEM`` section, i.e. ``COUPLINGS LINRES NVECT=n SPECIFY``. It is read at
each analytic coupling calculation; if it is missing the run stops with
``END OF FILE``.

Syntax::

   &VECTORS
   [FIXED]
   v(1,1) v(2,1) v(3,1)  v(1,2) ...  v(3,NAT)      (vector 1)
   ...                                             (vector n)
   &END

The first line after ``&VECTORS`` is a control line and is always consumed;
only its first word is examined (see :kw:`FIXED`). Then ``n`` vectors are
read (``n-1`` with FIXED), each with one list-directed READ, so a vector must
start on a new line but may span several lines. A vector contains the
:math:`3 N_{at}` Cartesian components for all atoms in the internal atom
order (species by species, as in ``&ATOMS``). Each vector is normalized after
reading. The vectors replace the automatically generated directions of the
iterative scheme (the vectors are still orthogonalized against the previous
ones).

Summary
-------

.. keyword-summary:: &VECTORS

Keywords
--------

.. keyword:: FIXED
   :section: &VECTORS

   If the first word of the first line of the ``&VECTORS`` section is
   ``FIXED``, the first search vector is not read but taken as the
   normalized constant-density coupling vector computed by the program,
   and only vectors 2 to ``n`` are read from the following lines. If the
   first line contains anything else, all ``n`` vectors are read starting
   with the second line.

   Example (``NVECT=3``, 2 atoms)::

      &VECTORS
      FIXED
        1.0 0.0 0.0   -1.0 0.0 0.0
        0.0 1.0 0.0    0.0 -1.0 0.0
      &END

   .. todo::

      MISSING: guidance on how to choose the vectors and literature
      reference (Ref. nonadiabatic of COUPLINGS).

.. source: src/cplngs_utils.mod.F90:1636
