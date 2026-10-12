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

   Not documented

.. keyword:: LOCALIZE
   :target: LOCALIZE HARDNESS
   :section: &HARDNESS

   Use localized orbitals in an orbital hardness calculation

.. keyword:: ORBITALS
   :section: &HARDNESS

   Specify the number of orbitals to be used in a hardness calculation on the next
   line.

.. keyword:: REFATOM
   :section: &HARDNESS

   Specify the reference atom to be used in a hardness calculation on the next
   line. This option is to be used together with the :kw:`ORBITALS`
   and :kw:`LOCALIZE`.
