Geometry Optimization
=====================

Any combination of methods for geometry optimization and wavefunction
optimization is allowed. Possible options for geometry optimization are GDIIS,
LBFGS, PRFO, RFO, BFGS and steepest descent. If you choose steepest descent for
both, geometry variables and the wavefunction, a combined method is used. For
all other combinations a full wavefunction optimization is performed between
changes of the ionic coordinates. The convergence criteria for the wavefunction
optimization can be adapted to the requirements of the geometry optimization
(CONVERGENCE ADAPT and CONVERGENCE ENERGY). The default options are GDIIS and
ODIIS. Some quasi-Newton methods (GDIIS, RFO and BFGS) are using the BFGS
method to update an approximate Hessian. At the beginning of a run the Hessian
can either be initialized as a unit matrix **HESSIAN UNIT** or with an
empirical force field. Two force fields are implemented: The **DISCO** and the
**SCHLEGEL** force field. The algorithm for the empirical force fields has to
identify bonds in the system. For unusual geometries this may fail and the
Hessian becomes singular. To prevent this you can add or delete bonds with the
keyword **CHANGE BONDS**.

The linear-scaling geometry optimizers (options LBFGS and PRFO) do not require
an approximate Hessian. To achieve linear scaling with the system size, the
L-BFGS optimizer starts from a unit Hessian and applies the BFGS update on the
fly using the history of the optimization. The P-RFO method can find transition
states by following eigenmodes of the Hessian. The mode to be followed does not
necessarily have to be the lowest eigenvalue initially (**PRFO MODE**). For
larger systems, only the reaction core should be handled by the P-RFO optimizer
(**PRFO NVAR** and **PRFO CORE**), and the environment is variationally
decoupled using the L-BFGS optimizer. The normal LBFGS options can be used for
the environment. The Hessian used for transition-state search therefore spans
only a subset of all degrees of freedom, is separate from the Hessian for the
other optimizers and the vibrational analysis, but it can be transferred into
the appropriate degrees of freedom of the regular Hessian (**HESSIAN
PARTIAL**). In order to allow negative eigenvalues, the Powell update instead
of BFGS is used for transition-state search.

Although tailored for large systems, the linear-scaling geometry optimizers are
suitable for smaller systems as well.
