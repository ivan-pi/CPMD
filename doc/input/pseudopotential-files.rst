.. _pp-file-format:

Pseudopotential file format
===========================

This page describes the sections of the pseudopotential files in the CPMD
format, which are read for every species of the ``&ATOMS`` section that does
not use another format (e.g. UPF). The general usage of pseudopotentials in the
input file is described in :ref:`S_Pseudopotentials`.

``&ATOM`` section
-----------------

Pseudopotential files in the CPMD format (species without the ``UPF``,
``BINARY``/``FORMATTED`` (Vanderbilt) or ``_EAM_`` labels) are read by
routine ``recpnew``. The file is searched first in the pseudopotential
library directory and then in the current directory. Sections are located by
substring search for their name; except ``&INFO`` they need no ``&END``.
The ``&ATOM`` section is mandatory and must start with four lines::

   &ATOM
    Z  = <atomic number>
    ZV = <valence charge>
    XC = <functional code> <Slater alpha>
    <fourth line, read and ignored (conventionally TYPE = ...)>
   &END

The values are read after the first ``=`` of each line. The last two digits
of the functional code (``code mod 100``) are compared with the gradient
correction code of the functional selected in ``&DFT`` (or the code returned
by the XC driver), and the Slater alpha with the alpha of the ``&DFT``
section; a mismatch only produces the warning
``XC FUNCTIONALS INCONSISTENT``.

.. todo::

   MISSING: table of functional codes; documentation of the UPF and
   Vanderbilt formats (read by other routines).

.. source: src/recpnew_utils.mod.F90:198

``&INFO`` section
-----------------

Free-text section of the pseudopotential file. Up to 60 lines (first 66
characters of each) until ``&END`` are stored and reproduced in the output.
Not to be confused with the ``&INFO`` section of the input file.

.. source: src/recpnew_utils.mod.F90:242

``&POTENTIAL`` section
----------------------

Mandatory section with the (semi-)local pseudopotential. The first line
selects the format: ``BHS``, ``CAR`` or ``GOEDECKER`` for analytic forms;
otherwise the numerical format is assumed (optionally preceded by
``LOGSCALE``). Below, ``LMAX`` is the value from the LMAX line of the species
in ``&ATOMS``.

.. keyword:: BHS
   :section: &POTENTIAL

   Analytic pseudopotential of Bachelet-Hamann-Schlueter type. The next
   line contains ``wrc1 rc1 wrc2 rc2`` (weights and radii of the two
   Gaussians of the local core part), followed, for each angular momentum
   :math:`l = 0 \ldots \mathrm{LMAX}`, by three lines ``rcl al bl``.

   .. todo::

      MISSING: explicit functional form and units.

.. source: src/recpnew_utils.mod.F90:264

.. keyword:: CAR
   :section: &POTENTIAL

   Analytic pseudopotential of von Barth-Car type. The next line contains
   the core radius ``rc``, followed, for each angular momentum
   :math:`l = 0 \ldots \mathrm{LMAX}`, by one line ``rcl al bl``.

   .. todo::

      MISSING: explicit functional form and units.

.. source: src/recpnew_utils.mod.F90:276

.. keyword:: GOEDECKER
   :section: &POTENTIAL

   Analytic separable pseudopotential of Goedecker type (the LMAX/LOC
   values in ``&ATOMS`` are ignored). Following lines:

   - ``nl``: number of nonlocal angular momentum channels
     (:math:`l = 0 \ldots nl-1`);
   - ``rloc``: radius of the local part;
   - ``nc c1 ... cnc``: number of coefficients and the coefficients of the
     local part;
   - for each channel :math:`l`: ``rl np h11 h12 ... h1np h22 ... hnpnp``,
     i.e. the radius, the number of projectors and the upper triangle of
     the :math:`h^l` matrix, row by row (list-directed, may continue on the
     next lines).

   Any first line containing ``GOEDECK`` is accepted.

   .. todo::

      MISSING: reference to the GTH papers and the maximum numbers of
      coefficients/projectors (``mcfun``, ``mpro``).

.. source: src/recpnew_utils.mod.F90:286

.. keyword:: LOGSCALE
   :section: &POTENTIAL

   Numerical pseudopotential. If the first line contains ``LOGSCALE`` (or
   both ``ECHELLE`` and ``LOG``), the mesh given in the file is used as is
   (it must be logarithmic) and the number of mesh points is read from the
   following line; in this case the ``&WAVEFUNCTION`` section must use the
   same mesh. Without ``LOGSCALE`` the first line holds the number of mesh
   points and the data are interpolated to a new logarithmic grid. Then one
   line per mesh point follows with the radius and ``LMAX+1`` potential
   values (:math:`l = 0 \ldots \mathrm{LMAX}`). The number of points must
   not exceed the spline-point limit (``MESH TOO BIG`` otherwise).

   Example::

      &POTENTIAL
         421
        0.0000E+00  -0.123E+02  -0.456E+01
        ...

.. source: src/recpnew_utils.mod.F90:333

``&WAVEFUNCTION`` section
-------------------------

Radial atomic pseudo-wavefunctions. The section is required for numerical
pseudopotentials used with ``KLEINMAN-BYLANDER`` when at least one nonlocal
channel exists. It is also read by the ``PSEUDO AO`` option of the
``&BASIS`` section and by the Hubbard-U code. The first line contains the
number of mesh points ``n``, followed by ``n`` lines with the radius and the
wavefunctions :math:`l = 0, 1, \ldots`. For ``recpnew``, ``LMAX+1`` columns
are read, or ``LMAX`` columns if LOC equals LMAX (the default); for
``PSEUDO AO`` the number of columns equals the number of shells given there.
With ``LOGSCALE`` in ``&POTENTIAL`` the number of mesh points must be the same
as in ``&POTENTIAL``.

.. todo::

   MISSING: normalization convention of the tabulated functions
   (presumably :math:`r R_l(r)`).

.. source: src/recpnew_utils.mod.F90:389

``&NLCC`` section
-----------------

Core charge for the nonlinear core correction. It is read only for species
with the ``NLCC`` label in ``&ATOMS``; if the label is missing but the file
contains ``&NLCC``, the warning ``PSEUDOPOTENTIAL GENERATED WITH NLCC`` is
printed. The first line selects the format:

.. keyword:: ANALYTIC
   :section: &NLCC

   The next line contains three numbers :math:`A`, :math:`B`, :math:`E`
   of the analytic core charge. The reciprocal-space expression in
   ``nlccset`` corresponds to
   :math:`\rho_c(r) = (A + B r^2)\, e^{-E r^2}`. Any first line containing
   ``ANAL`` is accepted.

.. source: src/recpnew_utils.mod.F90:475

.. keyword:: NUMERIC
   :section: &NLCC

   The next line contains the number of mesh points ``n``, followed by
   ``n`` lines with the radius and the core charge density. The data are
   interpolated to a logarithmic grid. Any first line containing ``NUME``
   is accepted.

   .. todo::

      MISSING: normalization convention of the tabulated core density
      (with or without :math:`4\pi r^2`).

.. source: src/recpnew_utils.mod.F90:479
