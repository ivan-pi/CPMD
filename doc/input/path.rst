.. _input-path:

``&PATH`` section
=================

Mean free energy path calculation (MFEP). This section is only evaluated if the :kw:`PATH MINIMIZATION` keyword is given in the ``&CPMD`` section.

This section contains specific information for mean free energy path
searches [Eijnden06]_. The space of collective variables in which the search
is performed has to be defined using **restraints** in the ``&ATOMS`` section
(see :ref:`sec:cnstr`). The initial string in this collective variable space
is read from the external file ``STRING.0``. This file contains one line per
replica, each giving first the replica index and then the list of collective
variable values. The base name of the directories in which each replica is run
should be specified with :kw:`FILEPATH`, or else with the environment variable
``CPMD_FILEPATH``. The code writes the files ``CONSTRAINT.x`` and
``METRIC.x`` during the execution, where ``x`` is the string number.

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

.. keyword:: FACTOR
   :section: &PATH

   Step for propagating string (see [Eijnden06]_).

   **Default** value is **1.0**

.. keyword:: NEQUI
   :section: &PATH

   Number of equilibration steps discarded to calculate the mean force.

.. keyword:: NLOOP
   :section: &PATH

   Maximum number of string searches for Mean Free Energy Path searches.

.. keyword:: NPREVIOUS
   :section: &PATH

   String index to restart from. Note that this is just for numbering files, the
   initial path in collective variables for the search is always *string.inp*.

.. keyword:: OUTPUT
   :target: OUTPUT PATH
   :section: &PATH
   :alternatives: [ALL, GROUPS, PARENT]

   Idem as above, here for Mean Free Energy Path runs.

.. keyword:: PRINT LEVEL
   :target: PRINT PATH
   :section: &PATH

   Idem as above, here for Mean Free Energy Path searches.

.. keyword:: PROCESSOR GROUPS
   :target: PROCESSOR GROUPS PATH
   :section: &PATH

   Idem as above, here for mean free energy path search.

.. keyword:: REPLICA NUMBER
   :section: &PATH

   Number of replicas along the string.
