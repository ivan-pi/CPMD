.. _input-path:

``&PATH`` section
=================

Mean free energy path calculation (MFEP). This section is only evaluated if the :kw:`PATH MINIMIZATION` keyword is given in the ``&CPMD`` section.

This section contains specific information for mean free energy path
searches [Eijnden06]_. The space of collective variables in which the search
is performed has to be defined using **restraints** in the ``&ATOMS`` section
(see :ref:`sec:cnstr`). The initial string in this collective variable space
is read from the external file ``STRING.n``, where ``n`` is the value of
:kw:`NPREVIOUS` (``STRING.0`` by default). This file contains one line per
replica, each giving first the replica index and then the list of collective
variable values. The base name of the directories in which each replica is run
should be specified with :kw:`FILEPATH`, or else with the environment variable
``CPMD_FILEPATH``. The code writes the files ``CONSTRAINT_p.x`` and
``METRIC_p.x`` during the execution, where ``p`` is the replica index and
``x`` the string number, and the free energy profile along the path to the
file ``PATHMIN``. After every string iteration the new string is written to
``STRING.x`` and its number to ``LATEST_STRING``.

The section starts with ``&PATH`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &PATH

Keywords
--------

.. keyword:: ALPHA
   :section: &PATH

   Smoothing parameter for iterating the string (see [Eijnden06]_).

   **Default** value is **0.2**

.. keyword:: CLASSICAL TEST
   :target: CLASSICAL TEST PATH
   :section: &PATH

   Test option inherited from the path-integral code. It sets the flag
   ``testpi``, which forces :kw:`REPLICA NUMBER` to 1 and prints
   ``THIS IS A TEST WITH NP=1 FOR COMPARISON WITH CLASSICAL RUNS``. Any
   line containing both ``CLASS`` and ``TEST`` is accepted.

   .. todo::

      MISSING: whether this option has any further use for mean free energy
      path runs (the only other place where it is tested is inside the
      GENERATE REPLICAS branch, which stops as not implemented).

.. source: src/pm_cntl_utils.mod.F90:135

.. keyword:: COLLECTIVE VARIABLES
   :section: &PATH

   Restricts the string to a subset of the restraints defined in the
   ``&ATOMS`` section (see :kw:`CONSTRAINTS ... END CONSTRAINTS`). The number
   ``ncvsp`` of collective variables is read from the next line, followed
   by ``ncvsp`` restraint indices (list-directed, may span several lines).
   Only these collective variables are read from and written to the
   ``STRING.n`` files, used in the convergence test (see
   :kw:`CONVERGENCE` in ``&PATH``) and printed. ``ncvsp`` must lie between 1
   and the number of restraints, and each index must be a valid restraint
   number. **Default** is to use all restraints in their input order. Any
   line containing both ``COLLE`` and ``VAR`` is accepted.

   Example::

      COLLECTIVE VARIABLES
        2
        1 3

.. source: src/pm_cntl_utils.mod.F90:221

.. keyword:: CONVERGENCE
   :target: CONVERGENCE PATH
   :section: &PATH

   Convergence criterion for the string iterations of a mean free energy
   path search. The tolerance is read from the next line. After each string
   update the maximum absolute change of the string positions (over all
   replicas and the selected collective variables, in internal units:
   lengths in Bohr, angles in radians) is compared with this value, and the
   loop over :kw:`NLOOP` strings is terminated when it is smaller.
   **Default** is 0.0, i.e. all :kw:`NLOOP` iterations are performed. Any
   line containing ``CONVERG`` is accepted.

   Example::

      CONVERGENCE
        1.0D-3

.. source: src/pm_cntl_utils.mod.F90:217

.. keyword:: FACTOR
   :section: &PATH

   Step for propagating string (see [Eijnden06]_).

   **Default** value is **1.0**

.. keyword:: GENERATE REPLICAS
   :target: GENERATE REPLICAS PATH
   :section: &PATH

   Requests that the replica coordinates be generated from the ``&ATOMS``
   configuration when :kw:`INITIALIZATION` is active. For mean free energy
   path runs this is **not implemented**: the run stops with
   ``GENERATE REPLICA NOT IMPLEMENTED``. It is also selected automatically
   for wavefunction optimizations without :kw:`RESTART` unless
   :kw:`READ REPLICAS` is given, so such runs need READ REPLICAS. Any line
   containing ``GENERATE`` is accepted.

.. source: src/pm_cntl_utils.mod.F90:147

.. keyword:: INITIALIZATION
   :target: INITIALIZATION PATH
   :section: &PATH

   Takes the starting coordinates of all replicas from
   :kw:`READ REPLICAS` (or :kw:`GENERATE REPLICAS`, which is not
   implemented for path runs) instead of the restart files. For
   wavefunction optimizations it is switched on automatically if no
   :kw:`RESTART` is requested; a wavefunction optimization that runs
   without INITIALIZATION reads the replica coordinates from the restart
   files (RESTART COORDINATES is implied). Any line containing
   ``INITIALIZA`` is accepted.

   .. todo::

      MISSING: behaviour for MD or geometry-optimization path runs when
      INITIALIZATION is given without READ REPLICAS (the replica array is
      then allocated but not filled; this looks like a missing check).

.. source: src/pm_cntl_utils.mod.F90:139

.. keyword:: NEQUI
   :section: &PATH

   Number of equilibration steps discarded to calculate the mean force.

.. keyword:: NLOOP
   :section: &PATH

   Maximum number of string searches for Mean Free Energy Path searches.

.. keyword:: NPREVIOUS
   :section: &PATH
   :options: {LATEST}

   Index ``n`` of the last string of a previous run, read from the next line
   (**default** 0). The initial string is read from the file ``STRING.n``
   (one line per replica: the replica index followed by the values of the
   collective variables; lengths in Bohr, or in Angstrom if :kw:`ANGSTROM`
   is set in ``&SYSTEM``; angles in degrees), and the new strings are
   numbered from ``n+1``. With ``LATEST`` no line is read and ``n`` is taken
   from the file ``LATEST_STRING``, which is written after every string
   iteration.

.. source: src/pm_cntl_utils.mod.F90:202

.. keyword:: OUTPUT
   :target: OUTPUT PATH
   :section: &PATH
   :alternatives: [ALL, GROUPS, PARENT]

   Idem as above, here for Mean Free Energy Path runs. ``GROUP`` is
   accepted as well as ``GROUPS``. **Default** is GROUPS.

.. keyword:: PRINT LEVEL
   :target: PRINT PATH
   :section: &PATH

   Idem as above, here for Mean Free Energy Path searches.

.. keyword:: PROCESSOR GROUPS
   :target: PROCESSOR GROUPS PATH
   :section: &PATH

   Idem as above, here for mean free energy path search.

.. keyword:: PROJOUT
   :section: &PATH

   When updating the string, project out the component of the (metric
   weighted) mean force that is parallel to the local string tangent
   before the step :math:`s_{new} = s_{old} - \mathrm{FACTOR}\cdot F`
   (see :kw:`FACTOR`). The tangent at replica :math:`i` is built from the
   neighbouring replica selected by the sign of the force along
   :math:`s_{i+1}-s_i`; the end points are kept fixed in both cases.
   **Default** is not to project out the forces.

   .. todo::

      MISSING: literature reference (presumably the string method of
      Ref. Eijnden06, but this is not stated in the code).

.. source: src/pm_cntl_utils.mod.F90:212

.. keyword:: READ REPLICAS
   :target: READ REPLICAS PATH
   :section: &PATH

   Reads the initial coordinates of all replicas from an external file
   whose name is given on the next line; used together with
   :kw:`INITIALIZATION`. The file contains, for each replica in turn, one
   header line followed by the Cartesian coordinates (three numbers per
   line) of all atoms in the order of the ``&ATOMS`` section. If the header
   line contains ``SCALE`` the coordinates are interpreted in units of the
   lattice vectors; if it contains ``ANGST`` they are converted from
   Angstrom to Bohr; otherwise they are taken in Bohr. Any line containing
   both ``READ`` and ``REPL`` is accepted.

   Example::

      READ REPLICAS
      replicas.xyz

.. source: src/pm_cntl_utils.mod.F90:152

.. keyword:: REPLICA NUMBER
   :section: &PATH

   Number of replicas along the string.
