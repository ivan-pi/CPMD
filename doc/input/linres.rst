.. _input-linres:

``&LINRES`` section
===================

General input for ``HARDNESS`` and ``TDDFT`` calculations.

The section starts with ``&LINRES`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &LINRES

Keywords
--------

.. keyword:: CONVERGENCE
   :target: CONVERGENCE LINRES
   :section: &LINRES

   Convergence criterion for linear response calculations.

   **Default** value is **10\ :math:`^{-5}`**.

.. keyword:: DIFF FORMULA
   :section: &LINRES

   Number of points used in finite difference formula for second derivatives of
   exchange–correlation functionals. Default is two point central differences.

.. keyword:: GAUGE
   :section: &LINRES
   :options: {PARA,GEN,ALL}

   Gauge of the linear-response wavefunctions. Default is the parallel-transport
   gauge (PARA) for closed-shell calculations and a sensible combination of the
   parallel-transport gauge and the full-rotation gauge (GEN) for all other cases.
   The full-rotation gauge can be enforced for all states by selecting ALL. See
   [lsets]_.

.. keyword:: HTHRS
   :section: &LINRES

   Threshold for Hessian in preconditioner for linear response optimizations.
   Default is 0.5.

.. keyword:: MAXSTEP
   :target: MAXSTEP LINRES
   :section: &LINRES

   Maximum number of optimization steps for linear response optimizations. Default
   is 1000.

.. keyword:: OPTIMIZER
   :section: &LINRES
   :alternatives: [SD,DIIS,PCG,AUTO]

   Optimizer to be used for linear response equations. Default is “AUTO” which
   will first use PCG, then switch to DIIS and finally switch to DIIS with full
   storage and state dependent preconditioner. :kw:`THAUTO` sets the two
   tolerances for when to do the switch.

.. keyword:: QS_LIMIT
   :target: QS-LIMIT
   :section: &LINRES

   Tolerance above which we use quadratic search algorithm in linear response
   calculations.

.. keyword:: STEPLENGTH
   :section: &LINRES

   Step length for steepest descent and preconditioned conjugate gradient methods
   used in linear response calculations. Default is 0.1.

.. keyword:: THAUTO
   :section: &LINRES

   The two values read from the next line control the switch to different
   optimizers for an automatic selection of optimizers during a linear response
   calculation. This also applies to the Z-vector optimization for TDDFT forces.
   The first value is the threshold for switching from conjugate gradients to DIIS
   (with compressed storage and averged preconditioner, subspace size defined with
   :kw:`ODIIS`). The second value is the threshold for switching to DIIS
   with full storage and state dependent preconditioner. See also
   :kw:`ZDIIS` for specification of the subspace size.

.. keyword:: XC_ANALYTIC
   :target: XC-ANALYTIC
   :section: &LINRES

   Use analytic second derivatives of the XC functional (only available for some
   LDA functionals)

.. keyword:: XC_DD_ANALYTIC
   :target: XC-DD-ANALYTIC
   :section: &LINRES

   Use analytic second derivatives of the XC functional, see
   Ref. [xcder]_ (only available for some LDA and gradient-corrected
   functionals). For the analytic third derivatives of some LDA XC functionals,
   :kw:`XC_ANALYTIC <XC-ANALYTIC>` can be combined with this keyword

.. keyword:: XC_EPS
   :target: XC-EPS
   :section: &LINRES

   Finite difference parameter for XC derivative. Default is
   :math:`5 \cdot 10^{-4}`.

.. keyword:: ZDIIS
   :section: &LINRES

   The subspace size for the optimizer is read from the next line.
