.. _inputfile:

Input File Reference
====================

The following sections **try** to explain the various keywords and syntax of a
CPMD input file. It is not meant to teach how to create good CPMD input files,
but as a reference manual.

Basic rules
-----------

-  **Warning:** Do not expect the input to be logical. The programmers logic
   may be different from yours.

-  **Warning:** This input description may not refer to the actual version of
   the program you are using. Therefore the ultimate and authoritative input
   guide is the source code. Most of the input file is read via the code in the
   files

   ``control_utils.mod.F90``, ``sysin_utils.mod.F90``, ``dftin_utils.mod.F90``,
   ``ratom_utils.mod.F90``, ``recpnew_utils.mod.F90``, ``detsp_utils.mod.F90``,
   ``proppt_utils.mod.F90``, ``setbasis_utils.mod.F90``,
   ``vdwin_utils.mod.F90``, ``respin_p_utils.mod.F90``,
   ``lr_in_utils.mod.F90``, ``orbhard_utils.mod.F90``,
   ``egointer_utils.mod.F90``, ``pi_cntl_utils.mod.F90``,
   ``cl_init_utils.mod.F90``, ``cplngs_utils.mod.F90``, ``mts_utils.mod.F90``

-  The input is free format except when especially stated

-  In most cases, only the first 80 characters of a line are read (exceptions
   are lists, that have to be on one line).

-  Lines that do not match a keyword are treated as comments and thus ignored.

   **Warning:** For most sections there will be a report of unrecognized
   keywords. For the ``&ATOMS`` this is not possible, so please check the atom
   coordinates in the output with particular care.

-  Take warnings seriously. There are a few warnings, that can be safely
   ignored under specific circumstances, but usually warnings are added to a
   program for a good reason.

-  The order of the keywords is arbitrary unless it is explicitly stated. For
   keywords that select one of many alternatives (e.g. the algorithm for
   wavefunction optimization), the last one ‘wins’.

-  Only keywords with capital letters match

-  Lists enclosed in **{ }** imply that you have to choose **exactly one** of
   the items

-  Lists enclosed in **[ ]** imply that you can choose **any number** of items
   on the same line

-  Arguments to a keyword are given on the following line(s) if not explicitly
   stated otherwise

-  The full keyword/input line has to be within columns 1 to 80

-  There are exceptions to those rules.

.. _sections:

Input sections
--------------

The input file is composed of different sections. Each section is started by
``&SECTIONNAME`` and ended by ``&END``. All input outside the sections is
ignored.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Section
     - Purpose
   * - ``&INFO``
     - A place to put comments about the job. The contents of this section
       are copied to the output file at the beginning of the calculation.
   * - :doc:`&CPMD <cpmd>`
     - General control parameters for the calculation (**required**).
   * - :doc:`&SYSTEM <system>`
     - Simulation cell and plane wave parameters (**required**).
   * - :doc:`&ATOMS <atoms>`
     - Atoms, pseudopotentials and related parameters (**required**).
       See also :ref:`S_Pseudopotentials`.
   * - :doc:`&DFT <dft>`
     - Exchange and correlation functional and related parameters.
   * - :doc:`&PROP <prop>`
     - Calculation of properties. Only fully evaluated if the
       :kw:`PROPERTIES` keyword is given in the ``&CPMD`` section.
   * - ``&BASIS``
     - Atomic basis sets for properties or the initial guess; see
       :ref:`input:basis`.
   * - :doc:`&RESP <resp>`
     - Response calculations. This section is always evaluated, even if it is
       not used.
   * - :doc:`&PTDDFT <ptddft>`
     - Propagation TDDFT for Ehrenfest dynamics and spectra.
   * - :doc:`&LINRES <linres>`
     - General input for ``HARDNESS`` and ``TDDFT`` calculations.
   * - :doc:`&HARDNESS <hardness>`
     - Input for ``HARDNESS`` calculations. Only evaluated if the
       :kw:`ORBITAL HARDNESS` ``LR`` keyword is given in the ``&CPMD`` section.
   * - :doc:`&TDDFT <tddft>`
     - Input for TDDFT calculations.
   * - :doc:`&QMMM <qmmm>`
     - Input for the Gromos QM/MM interface (see :ref:`sec:qmmm`).
       **Required** if the :kw:`QMMM` keyword is given in the ``&CPMD`` section.
   * - :doc:`&CLASSIC <classic>`
     - Simple classical code interface.
   * - ``&EXTE``
     - External field definition for the EGO QM/MM interface.
   * - :doc:`&VDW <vdw>`
     - Settings associated with van der Waals correction schemes. Only
       evaluated if one of :kw:`DCACP`, :kw:`VDW CORRECTION` or
       :kw:`VDW WANNIER` is given in the ``&CPMD`` section.
   * - :doc:`&MTS <mts>`
     - Parameters for the multiple time step MD scheme.
   * - :doc:`&PIMD <pimd>`
     - Path integral molecular dynamics (PIMD). Only evaluated if the
       :kw:`PATH INTEGRAL` keyword is given in the ``&CPMD`` section.
   * - :doc:`&PATH <path>`
     - Mean free energy path calculation (MFEP). Only evaluated if the
       :kw:`PATH MINIMIZATION` keyword is given in the ``&CPMD`` section.

All keywords are also listed in the :ref:`keyword index <cpmd-kwindex>`.
Additional components of CPMD input files that do not fit into the keyword
lists are explained in :doc:`details`.

.. toctree::
   :hidden:

   cpmd
   system
   atoms
   dft
   prop
   resp
   ptddft
   linres
   hardness
   tddft
   qmmm
   classic
   vdw
   mts
   pimd
   path
