.. _further_input:

Further details of the input
----------------------------

.. _S_Pseudopotentials:

Pseudopotentials
~~~~~~~~~~~~~~~~

The general format for entering the pseudo potentials in the input file is:

-  The input for a **new atom type** is started with a in the first column.
   This line further contains:

   -  the **file name** (*ECPNAME*) where to find the **pseudopotential**
      information starting in column 2

   -  and several **labels**:

      -  The **first label**

         **[GAUSS-HERMITE, KLEINMAN-BYLANDER]**

         specifies the method to be used for the calculation of the **nonlocal
         parts** of the **pseudopotential** [KB]_. It can be
         omitted for Vanderbilt pseudopotentials and Stefan Goedecker’s
         pseudopotentials. For semi-local pseudopotentials the default is
         Gauss-Hermite integration with 20 special points. The number of
         integration points can be changed using **GAUSS-HERMITE=xx**.

      -  It is further possible to specify **nonlinear core correction**
         [NLCC]_ **[NLCC]** and the width of the **ionic charge
         distribution [RAGGIO]**. (**Default** is **no NLCC** and the default
         value for **RAGGIO** is **1.2**.)

      -  The label **UPF** indicates that the pseudopotential was stored using
         the Universal Pseudopotential Format.

      -  For **Vanderbilt ultrasoft pseudopotentials** one of the following
         options has to be specified: **BINARY** indicates the binary version
         of the output file from Vanderbilt’s atomic code.

         **FORMATTED** indicates the formatted version of the Vanderbilt
         pseudopotential files after a conversion with the program ‘reform.f’
         from the Vanderbilt atomic code package (see section
         :ref:`intro:further`)

         For Vanderbilt pseudopotentials the option NLCC is ignored. The
         nonlinear core correction will always be used if the pseudopotential
         was generated with a partial core.

         It is strongly recommended to use only Vanderbilt pseudopotentials
         that were generated with a new version of Vanderbilt’s atomic code
         (version 7.3.0 or newer).

      -  The label **CLASSIC** indicates that the following atoms are to be
         treated with classical force fields only. See section &CLASSIC for
         more information.

      -  The label **\_EAM\_** indicates that the following atoms are treated
         using the EAM approach.

      -  The label **FRAC** indicates that the core charge of a pseudopotential
         should not be rounded for the calculation of the number of electrons
         (for pseudopotentials with fractional core charge).

      -  The label **ADD_H** indicates that the potential should used to
         saturate dangling bonds or “hydrogenize” united atom potentials in a
         CPMD/Gromos-QM/MM calculation (see section :ref:`sec:qmmm` for
         more details).

-  The next line contains information on the **nonlocality** of the
   **pseudopotential** (*LMAX, LOC, SKIP*) [BHS]_ [SGS]_ [TM]_.

-  On the following lines the **coordinates** for this **atomic species** have
   to be given.

   The first line gives the number of atoms (*NATOMS*) of the current type.

   Afterwards the coordinates of the atoms are listed (in Cartesian coordinates
   by default). For CPMD/Gromos-QM/MM calculation, however, the Gromos atom
   numbers have to be given instead of coordinates (see section
   :ref:`sec:qmmm` for more details).

The information on the **nonlocal part** of the pseudopotential
[BHS]_ [SGS]_ [TM]_ can be given in two different styles:

-  You can specify the maximum :math:`l` - quantum number with
   “**LMAX**\ =\ :math:`l`” where :math:`l` is S, P or D.

   If this is the only input, the program assumes that LMAX is the local
   potential (LOC).

   You can use another local function by specifying “**LOC**\ =\ :math:`l`”.

   In addition it is possible to assign the local potential to a further
   potential with “**SKIP**\ =\ :math:`l`”.

-  Alternatively you can specify these three angular quantum numbers by their
   numerical values (S=0, P=1, D=2) in the order “LMAX LOC SKIP”.

   If values for LOC and SKIP are provided outside the range 0 - LMAX the
   program uses the default.

**Examples:** The following lines are **equivalent**::

   LMAX=P
   LMAX=P LOC=P
   1   1   2
   1   2   2

**Note:**

Also for Vanderbilt and Goedecker pseudopotentials [GTH]_ this
line has to be in a **valid format**, but the actual values are **not** used.

.. _sec:cnstr:

Constraints and Restraints
~~~~~~~~~~~~~~~~~~~~~~~~~~

:kw:`CONSTRAINTS ... END CONSTRAINTS`

Within this input block you can specify several **constraints** and
**restraints** on the atoms.

Please note, that for calculations using the Gromos QM/MM-interface (see
section :ref:`sec:qmmm`) the atom indices refer to the ordering of the
atoms as it appears in the respective Gromos coordinate file. In all cases the
indices of dummy atoms start sequentially from total-number-of-atoms plus one.

The following suboptions are possible:

FIX ALL
   **All coordinates** of **all atoms** are kept fixed.

   For wavefunction optimization via simulated annealing.

FIX QM
   **All coordinates** of **all QM atoms** are kept fixed.

   This is the same as above unless you are running a QM/MM calculation with
   the Gromos interface code.

FIX MM
   **All coordinates** of **all MM atoms** are kept fixed.

   This is ignored unless you are running a QM/MM calculation with the Gromos
   interface code.

FIX SOLUTE
   **All coordinates** of **all solute atoms** are kept fixed.

   This is ignored unless you are running a QM/MM calculation with the Gromos
   interface code. The definition of what is a solute is taken from the
   respective GROMOS topology file.

FIX SEQUENCE
   **All coordinates** of a series of atoms are kept fixed.

   This keyword is followed by the index numbers of the first and the last
   atoms to be fixed in the next line. Example:

   ::

      FIX SEQUENCE
      5 25

   This example: all coordinates of atoms no. 5 to 25 are kept fixed.

FIX ELEMENT [SEQUENCE]
   **All coordinates** of all atoms belonging to the same **element** are kept
   fixed. This works across pseudopotential types or QM and MM atoms in case of
   a QM/MM calculation. The keyword is followed by the core charge of the
   respective element. With the optional SEQUENCE modifier two more numbers are
   read in, specifying the first and the last index of a sequence of atoms to
   which this keyword will be applied. Example:

   ::

      FIX ELEMENT
      8

   This example: all coordinates of oxygen atoms are kept fixed.

FIX PPTYPE [SEQUENCE]
   **All coordinates** of all atoms belonging to the same potential type are
   kept fixed. The keyword is followed by the atom type index number on the
   next line, corresponds to the sequence of how the atom types are specified
   in the &ATOMS section of the CPMD input. In case of a QM/MM calculation this
   is expanded to respective classical atom types. In this case the QM atom
   types come first followed by the GROMOS atom types. With the optional
   SEQUENCE modifier two more numbers are read in, specifying the first and the
   last index of a sequence of atoms to which this keyword will be applied.
   Example:

   ::

      FIX PPTYPE SEQUENCE
        2   5 25

   This example: atoms corresponding to the second atom type with an index between 5 and 25 are kept fixed.

FIX ATOMS
   **All coordinates** of certain atoms can be fixed.

   This keyword is followed by the number of atoms to be fixed and a list of
   these atoms specifying them by the number of their position in the input
   file (NOTE: in the file GEOMETRY.xyz the atoms have the same ordering).
   Example:

   ::

      FIX ATOMS
      5    2  5 20 21 23

   This example: all coordinates of atoms 2, 5, 20, 21, and 23 are kept fixed.

FIX COORDINATES
   **Certain coordinates** of atoms are fixed.

   This keyword is followed by the number of atoms with fixed coordinates and a
   list of these atoms together with flags indicating which coordinates are
   fixed. A zero indicates a fixed coordinate. Example:

   ::

      FIX COORDINATES
      2
      1   1 1 0
      4   0 1 0

   Two atoms have fixed coordinates: for atom 1, :math:`z` is fixed; for atom
   4, :math:`x` and :math:`z` are fixed.

FIX COM
   Fix the center of mass.

   **NOTE:** This currently works only for :kw:`OPTIMIZE GEOMETRY` and not for the :kw:`LBFGS`
   optimizer.

FIX STRUCTURE [SHOVE]
   This keyword starts a group of individual constraints where whole
   **structural units** can be fixed. The keyword is followed by the number of
   individual constraints on the next line.

   DIST
      :math:`n1` :math:`n2` :math:`R`

      Fixes the distance :math:`R` between the atoms :math:`n1` and :math:`n2`.

   STRETCH
      :math:`n1` :math:`n2` :math:`R`

      Fixes :math:`R^2` defined by the atoms :math:`n1` and :math:`n2`.

   DIFFER
      :math:`n1` :math:`n2` :math:`n3` :math:`R`

      Fixes :math:`R_{12}-R_{23}` defined by the atoms :math:`n1`, :math:`n2`,
      and :math:`n3`, where :math:`R_{ab}` is the distance between atoms a and
      b.

   DISAXIS
      :math:`n1` :math:`n2` :math:`k` :math:`R_0` [+1,-1,0] GROWTH aa DEST Rd

      Distance between two atoms :math:`n1` and :math:`n2` along :math:`x` or
      :math:`y` or :math:`z` direction. :math:`n1` :math:`n2` :math:`k` and the
      value :math:`R_0` are read next on the same line. Here :math:`k=1` means
      :math:`x`, :math:`k=2` means :math:`y` and :math:`k=3` means :math:`z`
      coordinate.

   BEND
      :math:`n1` :math:`n2` :math:`n3` :math:`\theta`

      Fixes the bending angle defined by the atoms :math:`n1`, :math:`n2` and
      :math:`n3`.

   TORSION
      :math:`n1` :math:`n2` :math:`n3` :math:`n4` :math:`\Theta`

      Fixes the torsion angle defined by the atoms :math:`n1`, :math:`n2`,
      :math:`n3` and :math:`n4`.

   OUTP
      :math:`n1` :math:`n2` :math:`n3` :math:`n4` :math:`\Theta`

      “Out of Plane”; Angle between plane (:math:`n1`, :math:`n2`, :math:`n3`)
      and atom :math:`n4` is fixed.

   RIGID
      :math:`nr` :math:`n1` :math:`n2` ... :math:`nx`

      Keeps the structure formed by the :math:`nr` atoms :math:`n1`,
      :math:`n2`, …

      You can put your atom index in several lines. The number of constraints
      **nfix** is equal to :math:`3nr-6` for :math:`nr>2` (:math:`nfix=1` for
      :math:`nr=2`).

   COORD
      :math:`n1` :math:`\kappa` :math:`Rc` :math:`d^0`

      Constraint on the coordination (number of atoms around a selected one
      within a specific spherical range of radius :math:`\sim Rc`) for atom
      :math:`n1`. The parameters :math:`\kappa` and :math:`Rc` for the Fermi
      function are given in Bohr (:math:`Rc`) and 1/Bohr (:math:`\kappa`), or
      in Å(:math:`Rc`) and 1/Å(:math:`\kappa`), if the keyword
      :kw:`ANGSTROM` is set. See Ref. [Sprik98a]_.

   COORSP
      :math:`n1` :math:`jsp` :math:`\kappa` :math:`Rc` :math:`d^0`

      Fixes the coordination number (CN) of one selected atom :math:`i` with
      respect to only one selected species :math:`jsp`. The CN is defined by a
      Fermi-like function as for :math:`COORD`, but in this case :math:`j` runs
      only over the atoms belonging to the selected species :math:`jsp`.

   COOR_RF
      :math:`n1` :math:`jsp` :math:`p` :math:`q` :math:`Rc` :math:`d^0`

      CN of one selected atom :math:`i` with respect to one selected species,
      :math:`jsp`. The CN value is calculated as the sum of rational functions

      .. math::

         CN_i = \sum_{j \neq i}^{n_{list}} \frac{1-\left(\frac{d_{ij}}{d^0}\right)^{p}}
                       {1-\left(\frac{d_{ij}}{d^0}\right)^{p+q}},

      where j runs over the indexes of the atoms belonging to :math:`jsp` or
      over the indexes given in the list :math:`j1 \cdots jn_{list}`.

   BNSWT
      :math:`n1` :math:`n2` :math:`p` :math:`q` :math:`Rc` :math:`d^0`

      Reciprocal CN between 2 selected atoms, defined with the same functional
      form as the one described for :math:`COOR\_RF`. This coordinate states
      the presence of the bond between the two atoms :math:`i` and :math:`j`.

   TOT_COOR
      :math:`isp` :math:`jsp` :math:`p` :math:`q` :math:`Rc` :math:`d^0`

      Average CN of the atoms belonging to a selected species :math:`isp` with
      respect to a second selected species, :math:`jsp`, or with respect to a
      given list of atoms, :math:`j1 \cdots jn_{list}`. The same functional
      forms and input options are used, as those described for
      :math:`COOR\_RF`, but the index of one selected species :math:`isp` is
      read in place of the index of one atom.

   :math:`n1`, ... are the atom numbers, :math:`R` distances and :math:`\Theta`
   angles. A function value of -999. for :math:`R` or :math:`\Theta` refers to
   the current value to be fixed. The constraint is linearly added to the CP
   Lagrangian according to the *Blue Moon* ensemble
   prescription\ [Sprik98b]_. The values of the Lagrange
   multipliers and of the actual constraint are printed in the file CONSTRAINT.

   The options **DIST, STRETCH, BEND, TORSION, OUTP, DIFFER, COORD, COORSP,
   COOR_RF, TOT_COOR** can have an optional additional keyword at the end of
   the line of the form

   ``DIST 1 2 -999. GROWTH 0.001``

   The keyword **GROWTH** indicates that the constraint value should be changed
   at each time step. The rate of change is given after the keyword in units
   per atomic time unit, i.e. **independent** from the current length of a time
   step.

   **Note:** In MD runs **only** the actual initial value (-999.) can be fixed.

   The **SHOVE** option requires an additional entry at the end of each
   constraint line. This entry has to be either :math:`-1`, :math:`0`, or
   :math:`1`. The constraint is then either fixed (:math:`0`) or allowed to
   shrink (:math:`-1`) or grow (:math:`1`).

RESTRAINTS [HYPERPLANE [K=scal
   ]]

   Defines restraints.

   nres
      Number of restraints.

   DIST
      :math:`n1` :math:`n2` :math:`R` :math:`kval`

      Restrains the distance :math:`R` between the atoms :math:`n1` and
      :math:`n2` by a harmonic potential having spring constant :math:`kval`.

   STRETCH
      :math:`n1` :math:`n2` :math:`R` :math:`kval`

      Restrains :math:`R^2` defined by the atoms :math:`n1` and :math:`n2` by a
      harmonic potential having spring constant :math:`kval`.

   DIFFER
      :math:`n1` :math:`n2` :math:`n3` :math:`R` :math:`kval`

      Restrains :math:`R_{12}-R_{23}` defined by the atoms :math:`n1`,
      :math:`n2`, and :math:`n3`, where :math:`R_{ab}` is the distance between
      atoms a and b by a harmonic potential having spring constant
      :math:`kval`.

   DISAXIS
      :math:`n1` :math:`n2` :math:`k` :math:`R_0` :math:`kval`

      Restraints the distance between two atoms :math:`n1` and :math:`n2` along
      :math:`x` or :math:`y` or :math:`z` direction. :math:`n1` :math:`n2`
      :math:`k` and the value :math:`R_0` are read next on the same line. Here
      :math:`k=1` means :math:`x`, :math:`k=2` means :math:`y` and :math:`k=3`
      means :math:`z` coordinate.

   BEND
      :math:`n1` :math:`n2` :math:`n3` :math:`\theta` :math:`kval`

      Restrains the bending angle defined by the atoms :math:`n1`, :math:`n2`
      and :math:`n3` by a harmonic potential having spring constant
      :math:`kval`.

   TORSION
      :math:`n1` :math:`n2` :math:`n3` :math:`n4` :math:`\Theta` :math:`kval`

      Restrains the torsion angle defined by the atoms :math:`n1`, :math:`n2`,
      :math:`n3` and :math:`n4` by a harmonic potential having spring constant
      :math:`kval`.

   OUTP
      :math:`n1` :math:`n2` :math:`n3` :math:`n4` :math:`\Theta` :math:`kval`

      “Out of Plane”; Angle between plane (:math:`n1`, :math:`n2`, :math:`n3`)
      and atom :math:`n4` is restrained by a harmonic potential having spring
      constant :math:`kval`.

   COORD
      :math:`n1` :math:`K` :math:`RC` :math:`C0` :math:`kval`

      Coordination restraint for atom :math:`n1`. The parameters :math:`K` and
      :math:`RC` for the Fermi-like function are given in Bohr (:math:`RC`) and
      1/Bohr (:math:`K`), or in Å(:math:`RC`) and 1/Å(:math:`K`), if the
      keyword :kw:`ANGSTROM` is set. See
      Ref. [Sprik98a]_. The harmonic potential spring constant is
      :math:`kval`.

   COORSP
      :math:`n1` :math:`is` :math:`K` :math:`RC` :math:`C0` :math:`kval`

      Retraint on the coordination number (CN) of one selected atom :math:`n1`
      with respect to a single selected species :math:`is`. The CN is defined
      by a Fermi-like function as for :math:`COORD`. As in all the above cases,
      the harmonic potential spring constant is :math:`kval`.

   COOR_RF
      :math:`n1` :math:`is` :math:`N` :math:`M` :math:`RC` :math:`C0`
      :math:`kval`

      Restraint on CN of one selected atom :math:`n1` with respect to one
      selected species, :math:`is`. The CN value is calculated as the sum of
      rational functions

      .. math::

         CN_i = \sum_{j \neq i}^{n_{list}} \frac{1-\left(\frac{d_{ij}}{d^0}\right)^{p}}
                       {1-\left(\frac{d_{ij}}{d^0}\right)^{p+q}},

      where j runs over the indexes of the atoms belonging to :math:`is` or
      over the indexes given in the list :math:`j1 \cdots jn_{list}`. As in all
      the above cases, the harmonic potential spring constant is :math:`kval`.

   BNSWT
      :math:`n1` :math:`n2` :math:`N` :math:`M` :math:`RC` :math:`C0`
      :math:`kval`

      Reciprocal CN between 2 selected atoms, defined with the same functional
      form as the one described for :math:`COOR\_RF`. This coordinate states
      the presence of the bond between the two atoms :math:`n1` and :math:`n2`.
      As in all the above cases, the harmonic potential spring constant is
      :math:`kval`.

   TOT_COOR
      :math:`is1` :math:`is2` :math:`N` :math:`M` :math:`RC` :math:`C0`
      :math:`kval`

      Average CN of the atoms belonging to a selected species :math:`is1` with
      respect to a second selected species, :math:`is2`, or with respect to a
      given list of atoms, :math:`j1 \cdots jn_{list}`. The same functional
      forms and input options are used, as those described for
      :math:`COOR\_RF`, but the index of one selected species :math:`isp` is
      read in place of the index of one atom. As in all the above cases, the
      harmonic potential spring constant is :math:`kval`.

   RESPOS
      :math:`n1`, :math:`x_0`, :math:`y_0`, :math:`z_0` :math:`d_0`
      :math:`kval`

      Restrains the position :math:`{\bf R}=(x,y,z)` of atom :math:`n1` to
      oscillate around :math:`{\bf R}_0=(x_0,y_0,z_0)` with a constraint
      harmonic potential :math:`V_c=(kval/2)(|{\bf R}-{\bf R}_0|-d_0)^2`
      [cco]_. The limits :math:`kval=0` and
      :math:`kval \to \infty` correspond to free and fixed atomic positions,
      respectively. The keyword GROWTH is not supposed to be used for this
      restraint. For the sake of clarity and consistency with the atomic units
      used through the code, coordinates and distances are expected to be in
      atomic units (not Å).

   :math:`n1`, ... are the atom numbers, :math:`R` distances and :math:`\Theta`
   angles. A function value of -999. for :math:`R` or :math:`\Theta` refers to
   the current value. The restraining potential is harmonic with the force
   constant :math:`kval`. The options can have an optional additional keyword
   at the end of the line of the form

   ``DIST 1 2 -999. 0.1 GROWTH 0.001``

   The keyword **GROWTH** indicates that the constraint value should be changed
   at each time step. The rate of change is given after the keyword in units
   per atomic time unit.

   If the keyword **HYPEPLANE** is set, the system is not restrained around a
   point in the collective variable space but in an hyperplane. This hyperplane
   is defined as going through a point in the collective variable space,
   defined from the :math:`R` and :math:`\Theta` above, and by a vector defined
   from the :math:`kval` values. **K=scal** applies a scaling to the vector
   defining the hyperplane so as to modulate the strength of the restraint.

   The energy formula for an hyperplane restraint is then:

   :math:`E_r=\frac{1}{2}\left((\vec{c}-\vec{c}_0)\cdot \vec{n}\right)^2`,

   where the vectors are vectors in the collective variable space.

   If a file **RESVAL** is found after parsing the input, the current restraint
   target values will be replaced by the values found in this file.

PENALTY
   The weight factors for the penalty function for stretches, bends and
   torsions are read from the next line.

.. _input:basis:

Atomic Basis Set
~~~~~~~~~~~~~~~~

The &BASIS section is used in CPMD only to provide an atomic basis set for
generating the initial guess and for analyzing orbitals. If the input file
contains no &BASIS section, a minimal Slater basis is used.

There have to be *number of species* different entries in this section.

The order of the basis sets has to correspond with the order of the atom types
in the section **&ATOMS …&END**.

With the keyword **SKIP** the species is skipped and the default minimal Slater
function basis is used.

Basis sets are either specified as Slater functions or given on an additional
input file.

The respective input formats are given below:

Slater type basis

::

   SLATER   nshell  [OCCUPATION]
     n1   l1    exp1
     ..   ..    ....
     nx   lx    expx
     [f1 f2 ... ]

Pseudo atomic orbitals

::

   PSEUDO AO nshell  [OCCUPATION]
     l1   l2   ..   lx     !a function with l=-1 is skipped
     [f1 f2 ... ]

Numerical functions

::

   *filename  nshell  FORMAT=n  [OCCUPATION]
     l1   l2   ..   lx
     [f1 f2 ... ]

Gaussian basis functions

::

   *filename  nshell  GAUSSIAN  [OCCUPATION]
     l1   l2   ..   lx
     [f1 f2 ... ]

Skip atom type and use default minimal slater function

::

   SKIP

``nshell`` is the number L-values  ``l1 l2 .. lx``  to be used.

``[f1 f2 ... ]`` is their occupation. If :kw:`CDFT` with the ``HDA``
option is active, a second line with occupations is read after the first one
(for every basis type that uses OCCUPATION); these occupations are used to
reset the initial occupations when the CDFT optimizers are re-initialized.

The format **PSEUDO AO** refers to the &WAVEFUNCTION section on the
corresponding pseudopotential file. For Vanderbilt ultrasoft
pseudopotentials the pseudo-atomic wavefunctions stored in the Vanderbilt
file are used instead of a ``&WAVEFUNCTION`` section.

With a L-value of -1 a specific function can be skipped.

A line is taken as a numerical (or Gaussian) basis entry if it contains a
\* anywhere (it does not have to be in the first column); a line containing
``SLATER`` is, however, always taken as a Slater basis entry. The file is
searched in the pseudopotential library directory. The default
format is 1, other possible formats are 2 and 3. The numbers correspond to the
format numbers in the old pseudopotential definitions for the atomic
wavefunctions:

- FORMAT=1 (default): first line ``mesh clog``, then ``mesh`` lines with
  the radius and the ``nshell`` radial functions.
- FORMAT=2: unformatted file, one record per shell containing ``mesh clog``
  and the pairs (r, f(r)).
- FORMAT=3: one comment line, then for each shell a line ``mesh clog``
  followed by ``mesh`` lines with r and f(r).

A value of ``clog`` larger than 1 is interpreted as the grid ratio and its
logarithm is taken. The functions are interpolated to a 256-point
logarithmic grid.

The format **GAUSSIAN** allows to use any linear combination of Gaussian
functions. The format of the file is as follows:

::

   Comment line
   number of l values (Lmax+1)
   (for each l = 0 ... Lmax)
     Comment line
     # of functions; # of exponents
     Comment line
     exp1 exp2 ... expn
     c11  c21      cn1
     c12  c22      cn2
     ...  ...      ...
     c1m  c2m      cnm

Van der Waals potential
~~~~~~~~~~~~~~~~~~~~~~~

This section ( &VDW …&END) contains information about the van der Waals
correction scheme. Currently, three schemes are available:

-  A Wannier-center based approach (:kw:`VDW WANNIER`)

-  A non-local potential-based approach
   (:kw:`DCACP` [opt-ecp04]_)

-  A parametrized force-field-like potential (:kw:`VDW CORRECTION`)

| The first two approaches both depend on the actual electron density of the
  system, the third is a post-SCF correction that is solely dependent on
  nuclear coordinates.
|

| **DCACP:** When using DCACP, the electron density is updated
  self-consistently according to the vdW potential that the system experiences.
  See the description of the keywords :kw:`DCACP Z=z`,
  :kw:`NO_CONTRIBUTION <NO CONTRIBUTION>` and
  :kw:`INCLUDE_METALS <INCLUDE METALS>` for customisation options.
|

**VDW CORRECTION:** Two major types of correction types are implemented in
connection with the :kw:`VDW CORRECTION` keyword: The one
described by M. Elstner *et al.*\ [Elstner]_ requires parameter
sets designed to use only in conjunction with a specific corresponding density
functional, the alternate parametrization by
S. Grimme\ [Grimme06]_ is less specific and therefore easier to
use and independent of the chose functional. Parameters for elements up to Xe
have been directly coded into CPMD. See the description of the keyword :kw:`VDW PARAMETERS` for more details regarding custom input.
