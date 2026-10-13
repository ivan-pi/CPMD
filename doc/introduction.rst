.. _about:

About this manual
=================

Many members of the CPMD consortium contributed to this manual. This version of the manual is based on a compilation by Barbara
Kirchner, Ari P. Seitsonen and Jürg Hutter working at the Physical Chemistry
Institute of the University of Zurich. Recent updates by Mauro Boero,
Alessandro Curioni, Jürg Hutter, Axel Kohlmeyer, Nisanth Nair and Wolfram
Quester.

CPMD is now developed as open source by the OpenCPMD community at
https://github.com/OpenCPMD. The original web site of the CPMD consortium,
www.cpmd.org, together with its mailing list (cpmd-list@cpmd.org) and contact
address (cpmd@cpmd.org), is no longer available. If you want to contribute or
have constructive criticism, please open an issue or a pull request in the
repository of the code, or use the discussions of the OpenCPMD organization on
GitHub.

.. _citation:

Citation
========

Publications of results obtained with CPMD should acknowledge its use by an
appropriate citation of the following kind:

*CPMD, https://github.com/OpenCPMD/CPMD,
Copyright IBM Corp 1990-2023,
Copyright MPI für Festkörperforschung Stuttgart 1997-2001*.

Important constants and conversion factors
==========================================

Input and output are in Hartree atomic units (a.u.), unless otherwise
explicitly mentioned.

IMPORTANT NOTICE:
~~~~~~~~~~~~~~~~~

As of CPMD version 3.15.1 all constants and conversion factors have been
consolidated and updated to the CODATA 2006 data
set\ [codata2006]_. For details see the file ``cnst.inc`` and
https://physics.nist.gov/constants.

+-------------------+---------------------------------------------------------+
| *quantity*        | *conversion factor*                                     |
+-------------------+---------------------------------------------------------+
| time step         | 1 a.u. = 0.02418884326505 fs                            |
+-------------------+---------------------------------------------------------+
| coordinates       | 1 Bohr = 1 a\ :math:`_0` = 0.52917720859 Å              |
+-------------------+---------------------------------------------------------+
| velocity          | 1 a.u. = 1 Bohr / 1 a.t.u. = 2187691.2541 m/s           |
+-------------------+---------------------------------------------------------+
| energy            | 1 Ha = 27.21138386 eV = 627.5094706 kcal/mol =          |
|                   | 2625.4996251 kJ/mol                                     |
+-------------------+---------------------------------------------------------+
| plane wave cutoff | 1 Ry = 1/2 Ha = 13.60569193 eV                          |
+-------------------+---------------------------------------------------------+
| dipole moment     | 1 a.u. = 2.5417462289 Debye                             |
+-------------------+---------------------------------------------------------+
| atomic mass       | 1 a.u. = 0.00054857990943 a.m.u                         |
+-------------------+---------------------------------------------------------+

.. _intro:further:

Recommendations for further reading
===================================

-  **General Introduction to Theory and Methods**

   Jorge Kohanoff, "Electronic Structure Calculation for Solids and Molecules",

   Cambridge University Press, 2006, ISBN-13 978-0-521-81591-8

   https://doi.org/10.1017/CBO9780511755613

-  **General introduction to Car-Parrinello simulation**

   D. Marx and J. Hutter, "Ab Initio Molecular Dynamics - Basic Theory and
   Advanced Methods”,

   Cambridge University Press, 2009

   D. Marx and J. Hutter, "Modern Methods and Algorithms of Quantum Chemistry",

   Forschungszentrum Jülich, NIC Series, Vol. 1 (2000), 301-449

   W. Andreoni and A. Curioni,

   "New Advances in Chemistry and Material Science with CPMD and Parallel
   Computing",

   *Parallel Computing*, 26 (2000) 819.

-  **Electronic Structure Theory**

   Richard M. Martin, "Electronic Structure: Basic Theory and Practical
   Methods’",

   Cambridge University Press, 2004, ISBN-13: 978-0-521-78285-2


-  **General overview about quantum simulation techniques**

   J. Grotendorst, D. Marx, and A. Muramatsu, *Quantum Simulations of Complex
   Many–Body Systems: From Theory to Algorithms*, (John von Neumann Institute
   for Computing, Forschungszentrum Jülich 2002); Printed Version: ISBN
   3–00–009057–6. Electronic version:
   https://juser.fz-juelich.de/record/24560/files/NIC-Band-10.pdf

-  **Molecular dynamics simulation**

   M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids* (Clarendon
   Press, Oxford, 1987; reprinted 1990).

   D. Frenkel and B. Smit, *Understanding Molecular Simulation – From
   Algorithms to Applications* (Academic Press, San Diego, 1996).

   M. E. Tuckerman and G. J. Martyna, J. Phys. Chem. B **104** 159 (2000).

-  **Pseudopotentials**

   https://github.com/OpenCPMD/pseudopotentials-extended-library (benchmarked
   pseudopotentials for CPMD)

   https://github.com/OpenCPMD/GTH-pseudopotentials (Goedecker-Teter-Hutter
   pseudopotentials for all elements)

   https://github.com/cp2k/cp2k-data/tree/master/potentials/Goedecker/cpmd
   (Goedecker pseudopotentials in CPMD format from CP2K)

   https://doi.org/10.17632/kygjf563w3.1 (the fhi98PP pseudopotential
   generator, CPC Program Library)

   https://www.physics.rutgers.edu/~dhv/uspp/ (Vanderbilt ultrasoft
   pseudopotentials)

   https://pseudopotentials.quantum-espresso.org/ (pseudopotential library of
   Quantum ESPRESSO)

-  **Parallelization & Performance**

   J. Hutter and A. Curioni, Parallel Computing **31**, 1 (2005).

   J. Hutter and A. Curioni, ChemPhysChem **6**, 1788-1793 (2005).

   C. Bekas and A.Curioni, Comp. Phys. Comm.\ **181**, 1057 (2010).
