.. _input-hardness:

``&HARDNESS`` section
=====================

Input for ``HARDNESS`` calculations. This section is only evaluated if the :kw:`ORBITAL HARDNESS` ``LR`` keyword is given in the ``&CPMD`` section.

The section starts with ``&HARDNESS`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &HARDNESS

Keywords
--------

.. keyword:: DIAGONAL
   :section: &HARDNESS
   :alternatives: [OFF]

   Selects the diagonal orbital hardness: one linear-response calculation is
   performed per reference orbital and the matrix :math:`\eta_{ij}` is printed
   and written to the file ``HARDNESS``. With ``OFF`` the generalized
   (extended) orbital hardness is computed, which needs one response
   calculation per pair of reference orbitals (:math:`n(n+1)/2`
   calculations). **Default** is DIAGONAL.

.. source: src/orbhard_utils.mod.F90:214

.. keyword:: LOCALIZE
   :target: LOCALIZE HARDNESS
   :section: &HARDNESS

   Use localized orbitals in an orbital hardness calculation

.. keyword:: ORBITALS
   :section: &HARDNESS

   Specify the number of orbitals to be used in a hardness calculation on the next
   line. Any line containing ``ORBITAL`` is accepted. Without :kw:`REFATOM` (or
   without :kw:`LOCALIZE <LOCALIZE HARDNESS>`) the last ``n`` occupied orbitals
   are used; **default** is all orbitals.

.. keyword:: REFATOM
   :section: &HARDNESS

   Specify the reference atom to be used in a hardness calculation on the next
   line. This option is to be used together with the :kw:`ORBITALS`
   and :kw:`LOCALIZE <LOCALIZE HARDNESS>`. With LOCALIZE the ``n`` orbitals
   whose Wannier centres (read from the file ``WANNIER_CENTER``) are closest
   to the reference atom are selected.
