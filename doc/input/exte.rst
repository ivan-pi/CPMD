.. _input-exte:

``&EXTE`` section
=================

The ``&EXTE ... &END`` section defines the external electrostatic field of
the classical environment for the file-based QM/MM interfaces selected with
:kw:`INTERFACE` ``EGO`` or ``INTERFACE GMX`` in the ``&CPMD`` section. In
these runs CPMD waits for a file ``QMCONTINUE`` (in the :kw:`FILEPATH`
directory); each time it appears, CPMD re-reads the input file given as the
first command-line argument: the atomic coordinates are taken from the
``&ATOMS`` section (for each ``*`` species line: one line is skipped, then the
number of atoms and the coordinates are read) and the external field from
``&EXTE``. If ``&EXTE`` is absent, no external field is applied. Since the
input file is re-read in every cycle, the classical MD program is expected to
update these sections before creating ``QMCONTINUE``.

The section has a fixed format. Every data block is preceded by one title
line, which is read and echoed but not interpreted. Coordinates are in Bohr.

Common head::

   <title>
   xmin xmax ymin ymax zmin zmax          (box of the real-space grid)
   <title>
   xmin xmax ymin ymax zmin zmax          (QM box)

``INTERFACE EGO`` then continues with the local expansions::

   <title>
   nloc
   (nloc times:)
   id  x  y  z                             (position of the expansion)
   c1
   c2 c3 c4
   c5 c6 c7 c8 c9 c10
   nincl
   <near-inclusion list, 10 entries of 7 characters per line>

The potential at a grid point :math:`\mathbf r` is evaluated from the
nearest expansion centre :math:`\mathbf R` as
:math:`-(c_1 + c_2\Delta x + c_3\Delta y + c_4\Delta z + \tfrac12 c_5\Delta x^2
+ \tfrac12 c_6\Delta y^2 + \tfrac12 c_7\Delta z^2 + c_8\Delta x\Delta y
+ c_9\Delta x\Delta z + c_{10}\Delta y\Delta z)`, plus the contribution of
the near atoms in its inclusion list. The inclusion list contains 0-based
indices into the near-atom list below; ``INT(nincl/10)+1`` lines are always
read (one line even for ``nincl = 0``).

``INTERFACE GMX`` instead continues with::

   <title>
   stpno                                   (current step number)
   <title>
   pwlmax
   <title>
   intmf                                   (update interval, intermediate layer)
   <title>
   outmf                                   (update interval, outermost layer)

Both interfaces then read the list of near MM atoms::

   <title>
   nnear
   (nnear times:)
   id  q  cutoff_width  cutoff_value  core_charge  sigma  x  y  z

For ``INTERFACE EGO`` a non-positive ``sigma`` sets the core charge to zero
and uses :math:`|\sigma|`. ``INTERFACE GMX`` finally reads two more layers::

   <title>
   nintm
   (nintm lines, only if mod(stpno,intmf) = 0:)  id  q  sigma  x  y  z
   <title>
   noutr
   (noutr lines, only if mod(stpno,outmf) = 0:)  id  q  sigma  x  y  z

The array limits ``maxind``, ``maxnear``, ``maxintm`` and ``maxoutr``
(module ``epot_types``) bound ``nloc``, ``nnear``/``nincl``, ``nintm`` and
``noutr``.

.. todo::

   MISSING: units of the charges, of ``cutoff_width``/``cutoff_value`` and
   of the expansion coefficients; meaning of ``pwlmax`` (described in the
   source only as "LMAX for partial wave expansion of MM-charge"); meaning
   of the 2nd to 4th columns of the near-atom list beyond the source
   comments; reference to the EGO/GROMACS interface papers.

.. source: src/egointer_utils.mod.F90:417
