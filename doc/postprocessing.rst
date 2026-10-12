.. _post:

Postprocessing
==============

The given output from a calculation with the CPMD code can be used for
postprocessing. There are several types of output. The most typical types of
output are density-like files, trajectories and/or xyz-files. These can be
visualized or analyzed with a number of different programs. Some of them are
(in no specific order):

- Molden (homepage: https://www.theochem.ru.nl/molden/)

- gOpenMol (no longer maintained)

- Molekel (homepage: https://github.com/ugovaretto/molekel)

- VMD (homepage: https://www.ks.uiuc.edu/Research/vmd/)

Starting with version 1.8.2 VMD does fully support the CPMD trajectory format,
xyz-movie format and Gaussian Cube files. Since version 1.8.3 VMD supports
periodic display of non-orthogonal supercells. A tutorial on Visualization and
Analysis of CPMD data with VMD was available at
http://www.theochem.ruhr-uni-bochum.de/go/cpmd-vmd.html; the tutorials of
Axel Kohlmeyer are now archived at
https://www.theochem.ruhr-uni-bochum.de/~legacy.akohlmey/.

Density files
-------------

List
~~~~

DENSITY.x, ELF, LSD_ELF, SPINDEN.x, WANNIER_1.x ...

.. _postprocessing-1:

Postprocessing
~~~~~~~~~~~~~~

These files are created in a binary format, they have to be transformed to a
Gaussian cube-File format to be readable by visualization programs. The
``cpmd2cube.x`` to convert the output is part of the CPMD add-ons
(https://github.com/OpenCPMD/Addons, directory ``cpmd2cube``) and is
used in the following way:

::

   cpmd2cube: Convert CPMD's Wannier-function files to cube
   usage is: cpmd2cube [options] Wannier_file [Wannier_file...]
      If you specify more than one Wannier file, they MUST have the
      same g-vectors and (for the moment) atom positions
      The program will create one cube file for each Wannier file
      and one pdb file with the atom positions
   Example:
      cpmd2cube.x WANNIER_1.*
   possible options are
      -v <verbosity>:
         <verbosity> is 0-2 (default is 1)
      -double:
         Read the density in double precision (default is single)
      -halfmesh:
         leave out half the grid points in each direction.
         Reduces the file size by 1/8th (on by default).
      -fullmesh:
         use the full real space grid.
      -n <n1> <n2> <n3>:
         change the REAL-space mesh. Default is to take the same mesh as CPMD
      -o <prefix>:
         specify the prefix of the name used for the cube and pdb-files
      -rep <n1> <n2> <n3>:
         replicate the cell n<j> times along the <j>-th direction by periodicity
      -shift <r1> <r2> <r3>:
         shift cube density by r1*a1+r2*a2+r3*a3
      -centre:
      -center:
         centre density around centre of mass of system.
      -inbox:
         put atoms inside unit cell centred around origin
      -rho:
      -dens:
         store the density instead of the wavefunction into the cube file.
      -psi:
      -wave:
         store the wavefunction instead of the density into the cube file.
      --:
         last option. Useful if you have a file with the same name as an option
      -h  or  -?  or  -help  or  --help  or no files:
         write this help

xyz-files
---------

.. _list-1:

List
~~~~

GEOMETRY.xyz, ION+CENTERS.xyz

.. _postprocessing-2:

Postprocessing
~~~~~~~~~~~~~~

These files can be directly read by a visualisation program. Note, that at the
current status it is useful to include a reference point (WANNIER REFERENCE)
for the ION+CENTERS.xyz which has to be put at the middle of the box.

TRAJECTORY-File
---------------

.. _list-2:

List
~~~~

TRAJECTORY

.. _postprocessing-3:

Postprocessing
~~~~~~~~~~~~~~

1. Looking at a movie

The TRAJECTORY files contains the coordinates and the velocities. To create a
movie file in the xyz-format, you have to transfer the coordinates from Bohr to
Å and you have to add the symbols of the atoms in the first position. Two lines
have to be at the beginning of each time step, from which the first line gives
the number of the total atoms. An .xyz file can also be recorded directly
during the simulation Using the :kw:`TRAJECTORY` keyword with the
option **XYZ**.

Please note, that CPMD does not apply the minimum image convention to these
trajectory files, i.e. atoms are **not** replaced by their images if they leave
the supercell.

2. Calculating radial pair distribution functions

The simplest analysis of the structure is given by the radial pair distribution
function g(r). This quantity is a to unity normalized function and describes
the probability of finding two atoms separated by a distance r relative to the
probability expected for a completely random distribution at the same density.
It is formally defined as:

.. math::

   \begin{aligned}
   g(r)=\rho^{-2}\left < \sum_i\sum_{i\neq j} \delta({\bf r}_i)
   \delta({\bf r}_j-{\bf r}) \right > \nonumber\\
   = \frac{V}{N^2}\left < \sum_i\sum_{i\neq j} \delta({\bf r}-{\bf r}_{ij})
    \right >
   \end{aligned}

with r being the atomic separation, :math:`\rho` the number density, N the
number of atoms, V the volume, **r** the atomic position, and
**r**\ :math:`_{ij}` the position of the atom i relative to the atom j. The
average :math:`<>` is taken over particles and time. Examples of code can be
found in the following references [Allen87]_ [Frenkel02]_.

The MOVIE format
----------------

Besides the TRAJECTORY file CPMD also produces a specially formatted trajectory
in the MOVIE file. This file contains the position of all atoms of the system
in Angstrom units at a rather low precision (:math:`10^{-4}`). The sampling of
the positions can be done independently from the trajectory file. The sets of
coordinates are following each other contiguously. The format of each line is

x-coordinate, y-coordinate, z-coordinate, atomic number, type .
