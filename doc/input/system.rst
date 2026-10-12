.. _input-system:

``&SYSTEM`` section
===================

Simulation cell and plane wave parameters (**required**).

The section starts with ``&SYSTEM`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &SYSTEM

Keywords
--------

.. keyword:: ACCEPTOR
   :section: &SYSTEM
   :alternatives: [HDASINGLE,WMULT]

   Set the :kw:`CDFT` acceptor atoms. Parameter NACCR must be specified
   next to the keyword. NACCR :math:`\in [1,2,...,N]` is the number of acceptor
   Atoms (:math:`N` being the total number of atoms). The indices of NACCR atoms
   separated by whitespaces are read from the next line.

   **HDASINGLE** **off** if set together with CDFT HDA, CPMD performs a
   constrained HDA calculation with only an ACCEPTOR group weight but different
   constraint values :math:`N_\text{c}`. **WMULT** **off** if set together with
   CDFT HDA, CPMD performs a constrained HDA calculation with two different an
   ACCEPTOR group weights for the two states.

   **HDASINGLE** and **WMULT** are mutually exclusive.

.. keyword:: ANGSTROM
   :section: &SYSTEM

   The atomic coordinates and the supercell parameters and several other
   parameters are read in Ångstroms.

   **Default** is **atomic units** which are always used internally. Not supported
   for :kw:`QMMM` calculations.

.. keyword:: CELL
   :section: &SYSTEM
   :options: {ABSOLUTE, DEGREE, VECTORS}

   The parameters specifying the super cell are read from the next line. Six
   numbers in the following order have to be provided: :math:`a`, :math:`b/a`,
   :math:`c/a`, :math:`\cos \alpha`, :math:`\cos \beta`, :math:`\cos \gamma`. For
   cubic phases, :math:`a` is the lattice parameter. CPMD will check those values,
   unless you turn off the test via :kw:`CHECK SYMMETRY`. With
   the keyword **ABSOLUTE**, you give :math:`a`, :math:`b` and :math:`c`. With the
   keyword **DEGREE**, you provide :math:`\alpha`, :math:`\beta` and
   :math:`\gamma` in degrees instead of their cosine. With the keyword
   **VECTORS**, the lattice vectors :math:`a1`, :math:`a2`, :math:`a3` are read
   from the next line instead of the 6 numbers. In this case the **SYMMETRY**
   keyword is not used.

.. keyword:: CHARGE
   :section: &SYSTEM

   The total charge of the system is read from the next line.

   **Default** is **0**.

.. keyword:: CHECK SYMMETRY
   :section: &SYSTEM
   :options: [OFF]

   The precision with which the conformance of the :kw:`CELL` parameters
   are checked against the (supercell) :kw:`SYMMETRY` is read from the
   next line. With older versions of CPMD, redundant variables could be set to
   arbitrary values; now **all** values have to conform. If you want the old
   behavior back, you can turn the check off by adding the keyword **OFF** or by
   providing a negative precision. **Default** value is: **1.0e-4**

.. keyword:: CLASSICAL CELL
   :section: &SYSTEM
   :alternatives: [ABSOLUTE, DEGREE]

   Not documented.

.. keyword:: CLUSTER
   :section: &SYSTEM

   Isolated system such as a molecule or a cluster. Same effect as
   :kw:`SYMMETRY` 0, but allows a non-orthorhombic cell. Only rarely
   useful.

.. keyword:: CONSTANT CUTOFF
   :section: &SYSTEM

   Apply a cutoff function to the kinetic energy term [bernasconi95]_
   in order to simulate constant cutoff dynamics. The parameters :math:`A`,
   :math:`\sigma` and :math:`E_o` are read from the next line (all quantities have
   to be given in Rydberg).

   .. math::

      G^2 \to G^2 + A \left[ 1 + \mbox{erf}
            \left( {\frac{1}{2} G^2 -  \frac{E_o}{\sigma}} \right) \right]

.. keyword:: COUPLINGS
   :section: &SYSTEM
   :options: {FD=\epsilon,PROD=\epsilon}
   :alternatives: [NAT]

   Calculate non-adiabatic couplings [nonadiabatic]_ using finite
   differences (FD and PROD are two different finite-difference approximations).
   The displacement :math:`\epsilon` is expected in atomic units. If
   NAT=\ :math:`n` is given, the coupling vector acting on only a subset of
   :math:`n` atoms is calculated. In this case, a line containing :math:`n` atom
   sequence numbers is expected. See :kw:`COUPLINGS NSURF`.

.. keyword:: COUPLINGS LINRES
   :section: &SYSTEM
   :options: {BRUTE FORCE,NVECT=n}
   :alternatives: [THR,TOL]

   Calculate non-adiabatic couplings [nonadiabatic]_ using
   linear-response theory. With BRUTE FORCE, the linear response to the nuclear
   displacements along all Cartesian coordinates is calculated. With
   NVECT=\ :math:`n`, at most :math:`n` cycles of the iterative scheme in
   [nonadiabatic]_ are performed. However, the iterative calculation
   is also stopped earlier if its contribution to the non-adiabatic coupling
   vector is smaller a given tolerance (TOL=\ :math:`C_{\mathrm{tol}}`). In the
   case of the iterative scheme, also the option THR can be given, followed by
   three lines each containing a pair of a threshold contribution to the
   non-adiabatic coupling vector and a tolerance for the linear-response
   wavefunction (see [nonadiabatic]_). Do not forget to include a
   &LINRES section in the input, even if the defaults are used. See :kw:`COUPLINGS NSURF`.

.. keyword:: COUPLINGS NSURF
   :section: &SYSTEM

   Required for non-adiabatic couplings: the Kohn-Sham states involved in the
   transition. For the moment, only one pair of states makes sense, NSURF=1. On
   the following line, the orbital numbers of the two Kohn-Sham states and a
   weight of 1.0 are expected. For singlet-singlet transitions, the ROKS-based
   Slater transition-state density (:kw:`LOW SPIN EXCITATION LSETS`) should be used. For doublet-doublet
   transitions, the local spin-density approximation (:kw:`LSD`) with the
   occupation numbers (:kw:`OCCUPATION`, :kw:`NSUP`,
   :kw:`STATES`) of the corresponding Slater transition-state density
   should be used.

.. keyword:: CUTOFF
   :section: &SYSTEM
   :alternatives: [SPHERICAL,NOSPHERICAL]

   The **cutoff** for the plane wave basis in **Rydberg** is read from the next
   line. The keyword **SPHERICAL** is used with k points in order to have
   :math:`|g + k|^2 < E_{cut}` instead of :math:`|g|^2 < E_{cut}`. This is the
   default.

.. keyword:: DENSITY CUTOFF
   :section: &SYSTEM
   :alternatives: [NUMBER]

   Set the plane wave energy cutoff for the density. The value is read from the
   next line. The density cutoff is usually automatically determined from the
   wavefunction :kw:`CUTOFF` via the :kw:`DUAL` factor.

   With the additional flag **NUMBER** the number of plane waves can be specified
   directly. This is useful to calculate bulk modulus or properties depending on
   the volume. The given energy cutoff has to be bigger than the one to have the
   required plane wave density number.

.. keyword:: DONOR
   :section: &SYSTEM

   Set the :kw:`CDFT` donor atoms. Parameter NACCR must be specified next
   to the keyword. NDON :math:`\in \mathbb{R}_+` is the number of Donor Atoms
   (:math:`N` being the total number of atoms). If NDON\ :math:`>0` the indices of
   NDON atoms separated by whitespace are read from the next line else only use an
   Acceptor group in the CDFT weight.

.. keyword:: DUAL
   :section: &SYSTEM

   The ratio between the wavefunction energy :kw:`CUTOFF` and the
   :kw:`DENSITY CUTOFF` is read from the next line.

   **Default** is **4**.

   There is little need to change this parameter, except when using ultra-soft
   pseudopotentials, where the wavefunction cutoff is very low and the
   corresponding density cutoff is too low to represent the augmentation charges
   accurately. In order to maintain good energy conservation and have good
   convergence of wavefunctions and related parameters, **DUAL** needs to be
   increased to values of 6–10.

   Warning: You can have some trouble if you use the **DUAL** option with the
   symmetrization of the electronic density.

.. keyword:: ENERGY PROFILE
   :section: &SYSTEM

   Perform an energy profile calculation at the end of a wavefunction optimization
   using the ROKS or ROSS methods.

.. keyword:: EXTERNAL FIELD
   :section: &SYSTEM

   Applies an external electric field to the system using the Berry phase. The
   electric field vector in AU is read from the next line.

.. keyword:: HFX CUTOFF
   :section: &SYSTEM

   Set an additional cutoff for wavefunction and density to be used in the
   calculation of exact exchange. Cutoffs for wavefunctions and densities are read
   from the next line in Rydberg units. Defaults are the same cutoffs as for the
   normal calculation. Only lower cutoffs than the defaults can be specified.

.. keyword:: ISOTROPIC CELL
   :section: &SYSTEM

   Specifies a constraint on the super cell in constant pressure dynamics or
   geometry optimization. The shape of the cell is held fixed, only the volume
   changes.

.. keyword:: KPOINTS
   :section: &SYSTEM
   :notes: options

   With no option, read in the next line with the number of k-points and for each
   k-point, read the components in the Cartesian coordinates
   (units :math:`2\pi/a`) and the weight.

   MONKHORST-PACK
      Read in the next line three numbers for the Monkhorst-Pack mesh. The program
      calculates then the special k-points. With the keyword **SHIFT=kx ky kz** in
      the same line, you can precise the constant vector shift.

   SYMMETRIZED
      Symmetrized special k-points mesh (useful if you use a constant vector
      shift).

   FULL
      Construct full Monkhorst-Pack mesh with only inversion symmetry. Useful for
      molecular dynamics simulation The keywords **SYMMETRIZED FULL** preserves
      all symmetry of Bravais lattice so there is no need to symmetrize density
      and forces.

   SCALED
      You can give k-points in reciprocal space coordinates.

   BANDS
      This option is to calculate the band structure.

      For each line you have to specify the number of k-points for the band, the
      initial and the final k-point. To finish the input, put:

      0 0. 0. 0. 0. 0. 0.

   BLOCK=n [OPTIONS]
      The block option, specifies the number of k-points in the memory. The
      program uses a swap file to store the wavefunctions only by default. With
      the following options, you can change this behaviour:

      ALL
         Three swap files are used to store wavefunctions and others arrays
         related to k-points. Swap files are in the current directory or the
         temporary directory given by environment variable TMPDIR. The use of
         memory is smaller than with the above option.

      CALCULATED
         One swap file is used to store only wavefunctions. The other arrays
         related to k-points are calculated each time if needed.

      NOSWAP
         The wavefunctions are not swapped. This is useful to calculate
         eigenvalues for each k point with few memory used.

         **Warning:** The wavefunctions calculated are irrelevant. You have to
         specify explicitly some other options to use it:

         | MAXSTEP 1 and
         | STORE OFF WAVEFUNCTIONS DENSITY POTENTIAL.

.. keyword:: LOW SPIN EXCITATION
   :section: &SYSTEM
   :alternatives: [ROKS,ROSS,ROOTHAAN,CAS22]

   Use the low spin excited state functional [Frank98]_. For ROKS
   calculations, see also the :kw:`ROKS` keyword in the &CPMD-section.

.. keyword:: LOW SPIN EXCITATION LSETS
   :section: &SYSTEM

   Slater transition-state density with restricted open-shell Kohn-Sham (low spin
   excited state). Currently works only with ROKS but not with ROSS, ROOTHAAN, or
   CAS22. See Ref. [lsets]_.

.. keyword:: LSE PARAMETERS
   :section: &SYSTEM

   Determines the energy expression used in LSE calculations. The two parameters
   LSEA and LSEB are read from the next line.

   .. math:: E = \mbox{LSEA} \cdot E(Mixed) + \mbox{LSEB} \cdot E(Triplet)

   The default (LSEA :math:`= 2` and LSEB :math:`= 1`) corresponds to singlet
   symmetry. For the lowest triplet state, the :kw:`LSE PARAMETERS` must be set to 0 and 1 (zero times mixed
   state plus triplet). See ref [Frank98]_ for a description of the
   method.

.. keyword:: MESH
   :section: &SYSTEM

   The number of **real space mesh** points in :math:`x-`, :math:`y-` and
   :math:`z-`\ direction is read from the next line.

   If the values provided by the user are not compatible with the plane-wave
   cutoff or the requirements of the FFT routines the program chooses the next
   bigger valid numbers.

   **Default** are the **minimal values** compatible with the energy cutoff and
   the **FFT** requirements.

.. keyword:: MULTIPLICITY
   :section: &SYSTEM

   This keyword only applies to LSD calculations.

   The multiplicity (2\ :math:`S`\ +1) is read from the next line.

   **Default** is the **smallest possible** multiplicity.

.. keyword:: NSUP
   :section: &SYSTEM

   The number of states of the same spin as the first state is read from the next
   line.

   This keyword makes only sense in spin-polarized calculations (keyword
   :kw:`LSD`).

.. keyword:: OCCUPATION
   :section: &SYSTEM
   :alternatives: [FIXED]

   The occupation numbers are read from the next line.

   This keyword must be preceded by :kw:`STATES`. The FIXED option fixes
   the occupation numbers for the diagonalization scheme, otherwise this option is
   meaningless.

.. keyword:: POINT GROUP
   :section: &SYSTEM
   :alternatives: [MOLECULE], [AUTO], [DELTA=delta]

   The point group symmetry of the system can be specified in the next line. With
   the keyword *AUTO* in the next line, the space group is determined
   automatically. This affects the calculation of nuclear forces and ionic
   positions. The electronic density and nuclear forces are symmetrized in
   function of point group symmetry. The group number is read from the next line.

   Crystal symmetry groups:

   ::

                     1  1 (c1)     9   3m (c3v)   17 4/mmm (d4h)   25  222 (d2)
                     2 <1>(ci)    10  <3>m(d3d)   18   6   (c6)    26  mm2 (c2v)
                     3  2 (c2)    11   4  (c4)    19  <6>  (c3h)   27  mmm (d2h)
                     4  m (c1h)   12  <4> (s4)    20   6/m (c6h)   28  23  (t)
                     5 2/m(c2h)   13  4/m (c4h)   21   622 (d6)    29  m3  (th)
                     6  3 (c3)    14  422 (d4)    22   6mm (c6v)   30  432 (o)
                     7 <3>(c3i)   15  4mm (c4v)   23  <6>m2(d3h)   31 <4>3m(td)
                     8 32 (d3)    16 <4>2m(d2d)   24  6/mmm(d6h)   32  m3m (oh)

   You can specify the point group by its name using the keyword *NAME=* followed
   by the name of the point group (one of both notations).

   For molecular point groups the additional keyword *MOLECULE* has to be
   specified. The Schönflies symbol of the group is read in the following format
   from the next line:

   *Group symbol; order of principle axis*

   Possible group symbols are any Schönflies symbol with the axis number replaced
   by :math:`n` (e.g. DNH). For molecular point groups a special orientation is
   assumed. The principle axis is along :math:`z` and vertical symmetry planes are
   orthogonal to :math:`x`.

   *DELTA=* specifies the required accuracy (default=\ :math:`10^{-6}`).

   With the keyword **AUTO**, the point group is determined automatically.

.. keyword:: POISSON SOLVER
   :section: &SYSTEM
   :options: {HOCKNEY, TUCKERMAN, MORTENSEN}
   :alternatives: [PARAMETER]

   This keyword determines the method for the solution of the Poisson equation for
   isolated systems. Either Hockney’s method [Hockney70]_ or Martyna
   and Tuckerman’s method [Martyna99]_ is used. The smoothing
   parameter (for Hockney’s method) or :math:`L \times \alpha` for Tuckerman’s
   method can be read from the next line using the **PARAMETER** keyword.

   For more information about the usage of this parameter see also section
   :ref:`hints:symm0`.

.. keyword:: POLYMER
   :section: &SYSTEM

   Assume **periodic boundary** condition in **:math:`x`-direction**.

.. keyword:: PRESSURE
   :section: &SYSTEM

   The **external pressure** on the system is read from the next line (in
   **kbar**).

.. keyword:: REFERENCE CELL
   :section: &SYSTEM
   :alternatives: [ABSOLUTE, DEGREE, VECTORS]

   This cell is used to calculate the Miller indices in a constant pressure
   simulation. This keyword is only active together with the option
   **PARRINELLO-RAHMAN**.

   The parameters specifying the reference (super) cell are read from the next
   line.

   Six numbers in the following order have to be provided: :math:`a`, :math:`b/a`,
   :math:`c/a`, :math:`\cos \alpha`, :math:`\cos \beta`, :math:`\cos \gamma`.

   The keywords **ABSOLUTE** and **DEGREE** are described in **CELL** option.

.. keyword:: SCALE
   :section: &SYSTEM
   :alternatives: [CARTESIAN]
   :notes: [S=sascale] [SX=sxscale] [SY=syscale] [SZ=szscale]

   **Scale atomic coordinates** of the system with the lattice constants (see
   **CELL**). You can indicate an additional scale for each axis with the options
   **SX**, **SY** and **SZ**. For instance, if you indicate SX=sxscale, you give
   your x-coordinates between :math:`0.` and sxscale (by default :math:`1.`). This
   is useful when you use many primitive cells. With the keyword **CARTESIAN**,
   you specify that the given coordinates are in Cartesian basis, otherwise the
   default with the **SCALE** option is in direct lattice basis. In all cases, the
   coordinates are multiplied by the lattice constants. If this keyword is present
   an output file GEOMETRY.scale is written. This file contains the lattice
   vectors in Åand atomic units together with the atomic coordinates in the direct
   lattice basis.

.. keyword:: STATES
   :section: &SYSTEM

   The number of states used in the calculation is read from the next line.

   This keyword has to precede the keyword **OCCUPATION**.

.. keyword:: STRESS TENSOR
   :target: STRESS TENSOR SYSTEM
   :section: &SYSTEM

   In extension to the keyword PRESSURE the complete **stress tensor** in kbar can
   be specified. The **stress** on the system is read in the form:

   .. container:: center

      | :math:`t_{11}\ t_{12}\ t_{13}`
      | :math:`t_{21}\ t_{22}\ t_{23}`
      | :math:`t_{31}\ t_{32}\ t_{33}`

.. keyword:: SURFACE
   :section: &SYSTEM
   :alternatives: [XY, YZ, ZX]

   By default, if nothing is specified, assume **periodic boundary** condition in
   **:math:`x`- and :math:`y`-direction**. With the extra keywords *XY*, *YZ* or
   *ZX*, the periodicity of the systems is assumed to be along :math:`(x,y)`,
   :math:`(y,z)` or :math:`(z,x)`, respectively.

.. keyword:: SYMMETRIZE COORDINATES
   :section: &SYSTEM

   **Input coordinates** are **symmetrized** according to the **point group**
   specified.

   This only makes sense when the structure already is close to the symmetric one.

.. keyword:: SYMMETRY
   :section: &SYSTEM

   The **supercell symmetry type** is read from the next line.

   You can put a number or a keyword.

   .. container:: description

      **ISOLATED** system in a cubic/orthorhombic
      box [Hockney70]_ [Landman]_ with ISOLATED MOLECULE option
      activated. By default the Hockney method (see :kw:`POISSON SOLVER`) is used for solving the Poisson equations.
      You can use this option in combination with :kw:`POLYMER` or
      :kw:`SURFACE` for systems that are periodic in only 1 or 2
      dimensions. The default Poisson solver is MORTENSEN in this case. See the
      Hints and Tricks section for some additional requirements when calculating
      isolated system.

      Simple **CUBIC**

      **FACE CENTERED CUBIC** (**FCC**)

      **BODY CENTERED CUBIC** (**BCC**)

      **HEXAGONAL**

      **TRIGONAL** or **RHOMBOHEDRAL**

      **TETRAGONAL**

      **BODY CENTRED TETRAGONAL** (**BCT**)

      **ORTHORHOMBIC**

      **MONOCLINIC**

      **TRICLINIC**

   Warning: This keyword should not be used with the keyword **CELL VECTORS**.

.. keyword:: TESR
   :section: &SYSTEM

   The number of additional supercells included in the real space sum for the
   Ewald term is read from the next line. Default is 0, for small unit cells
   larger values (up to 8) have to be used.

.. keyword:: WCUT
   :section: &SYSTEM
   :notes: CUT

   Set the radial :kw:`CDFT` weight cutoff for all atom species to CUT,
   which is specified next to the keyword. Default is a species specific cutoff at
   the distance where the magnitude of the respective promolecular density is
   smaller than :math:`10^{-6}`.

.. keyword:: WGAUSS
   :section: &SYSTEM
   :notes: NWG

   Use Gaussian weight functions instead of Hirshfeld promolecular orbitals in the
   :kw:`CDFT` weight. Parameter NWG is specified next to the keyword and
   has to be equal to the number of different atom species in the calculation. The
   Gaussian widths :math:`\sigma_i` of the species :math:`i` are read from
   subsequent lines.

.. keyword:: ZFLEXIBLE CELL
   :section: &SYSTEM

   Specifies a constraint on the super cell in constant pressure dynamics or
   geometry optimizations. The supercell may only shrink or grow in z-direction.
   Should be very useful for “dense slab” configurations, e.g. a water layer
   between solid slabs.

   **Please note:** this is by no means intended to give a statistically
   meaningful ensemble, but merely to provide a tool for efficient equilibration
   of a specific class of system.
