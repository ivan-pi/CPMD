.. _hints:symm0:

Cell Size for Calculations with **SYMMETRY 0**
==============================================

Calculations of isolated systems (i.e. decoupling of the electrostatic images
in the Poisson solver) are initialized with:

::

   SYMMETRY
     0

(see :kw:`SYMMETRY`).

The box is assumed to be orthorhombic. With the additional options
:kw:`SURFACE` or :kw:`POLYMER` periodic boundary
conditions in two or one dimensions, respectively, are assumed.

Three different kinds of :kw:`POISSON SOLVER`
{HOCKNEY,TUCKERMAN,MORTENSEN} are available. All methods require that the
charge density is zero at the border of the box. For normal systems and the
Hockney solver, this means that about 3 Angstrom space between the outermost
atoms and the box should be enough. For large molecules, Tuckerman will require
a considerably higher margin, see below. However, for some systems and for high
accuracy these rules of thumb may not be enough. Some methods have additional
requirements (see below).

Note that the :kw:`ISOLATED MOLECULE` keyword has only an
effect on the calculation of the degrees of freedom (3N-6 vs. 3N-3 for periodic
systems). The main purpose of :kw:`CENTER MOLECULE` ON/OFF
is to center the molecule (center of mass) in the box. This is needed for the
HOCKNEY Poisson solver. This solver gives wrong results if the charge density
is not centered in the computational box. All other solvers behave like the
periodic counterpart, i.e. the relative position of the charge density and the
box are not important.

Further requirements on the Poisson solvers:

HOCKNEY Method:

-  molecule has to be in the center of the box

-  box size molecule + 3 Å border

-  expensive for very small systems

-  not available for some response calculations

-  :kw:`POLYMER` is available but gives (currently, Version 3.9)
   wrong results.

-  :kw:`SURFACE` is available and works.

TUCKERMAN-MARTYNA Method:

-  box size : molecule + 3 Å border **and** - imperatively! - *twice* the size
   of the electron charge distribution (check the density after a single-point
   wavefunction optimization).

-  expensive for large systems, smaller boxes might be used without loosing too
   much accuracy

-  :kw:`SURFACE` or :kw:`POLYMER` are not available

MORTENSEN Method:

-  same as TUCKERMAN, but using analytic formula made possible by using special
   boundary conditions (sphere, rod)

-  :kw:`SURFACE` and :kw:`POLYMER` are available and
   should be safe to use (MORTENSEN is default for SURFACE and POLYMER)

-  If you do an isolated system calculation, your cell has to be cubic, if you
   use :kw:`POLYMER` cell dimensions b and c have to be equal.

Finally, for many systems using a large enough cell and periodic boundary
conditions is also an option. In general, the computed properties of molecules
should be independent of the scheme used (either pbc or isolated box) except in
difficult cases such as charged molecules, where the calculation in an isolated
box is recommended. The PBC calculation is always cheaper for a box of the same
size, so for a neutral molecule such as water molecule you would save time and
memory by not using :kw:`SYMMETRY` 0.
