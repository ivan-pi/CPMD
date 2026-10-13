.. _FILES:

Files
=====

Incomplete list of the files used or created by CPMD:

+--------------------------------------+--------------------------------------+
| *File*                               | *Contents*                           |
+--------------------------------------+--------------------------------------+
| RESTART                              | Restart file                         |
+--------------------------------------+--------------------------------------+
| RESTART.x                            | Old/New restart files                |
+--------------------------------------+--------------------------------------+
| LATEST                               | Info file on the last restart file   |
+--------------------------------------+--------------------------------------+
| GEOMETRY                             | Current ionic positions and          |
|                                      | velocities                           |
+--------------------------------------+--------------------------------------+
| GEOMETRY.xyz                         | Current ionic positions and          |
|                                      | velocities in Å                      |
+--------------------------------------+--------------------------------------+
| GEOMETRY.scale                       | Current unit cell vectors in Åand    |
|                                      | atomic units                         |
+--------------------------------------+--------------------------------------+
| and ionic positions in scaled        |                                      |
| coordinates                          |                                      |
+--------------------------------------+--------------------------------------+
| GSHELL                               | **G**\ :math:`^2` (NOT normalized),  |
|                                      | G-shells :math:`|{\bf G}|` in a.u.   |
+--------------------------------------+--------------------------------------+
| and related shell index              |                                      |
+--------------------------------------+--------------------------------------+
| TRAJECTORY                           | All ionic positions and velocities   |
|                                      | along the trajectory                 |
+--------------------------------------+--------------------------------------+
| TRAJEC.xyz                           | All ionic positions along the        |
|                                      | trajectory in xyz-format             |
+--------------------------------------+--------------------------------------+
| TRAJEC.dcd                           | All ionic positions along the        |
|                                      | trajectory in dcd-format             |
+--------------------------------------+--------------------------------------+
| GEO_OPT.xyz                          | All ionic positions along the        |
|                                      | geometry optimization                |
+--------------------------------------+--------------------------------------+
| in xyz-format                        |                                      |
+--------------------------------------+--------------------------------------+
| ENERGIES                             | All energies along the trajectory    |
+--------------------------------------+--------------------------------------+
| MOVIE                                | Atomic coordinates in Movie format   |
+--------------------------------------+--------------------------------------+
| STRESS                               | The "trajectory" of stress tensors   |
+--------------------------------------+--------------------------------------+
| CELL                                 | The "trajectory" of the unit cell    |
+--------------------------------------+--------------------------------------+
| NOSE_ENERGY                          | The kinetic, potential and total     |
|                                      | energies of the Nose-Hoover          |
|                                      | thermostat(s)                        |
+--------------------------------------+--------------------------------------+
| NOSE_TRAJEC                          | The "trajectory" of the Nose-Hoover  |
|                                      | variables and velocities             |
+--------------------------------------+--------------------------------------+
| CONSTRAINT                           | The "trajectory" of                  |
|                                      | constraint/restraint forces          |
+--------------------------------------+--------------------------------------+
| METRIC                               | The "trajectory" of collective       |
|                                      | variable metric (restraints only)    |
+--------------------------------------+--------------------------------------+
| DIPOLE                               | The "trajectory" of dipole moments   |
+--------------------------------------+--------------------------------------+
| dipole.dat                           | The dipole file produced in the      |
|                                      | spectra calculation using the        |
|                                      | propagation TDDFT scheme             |
+--------------------------------------+--------------------------------------+
| DENSITY.x                            | Charge density in Fourier space      |
+--------------------------------------+--------------------------------------+
| SPINDEN.x                            | Spin density in Fourier space        |
+--------------------------------------+--------------------------------------+
| ELF                                  | Electron localization function in    |
|                                      | Fourier space                        |
+--------------------------------------+--------------------------------------+
| LSD_ELF                              | Spin polarized electron localisation |
|                                      | function                             |
+--------------------------------------+--------------------------------------+
| ELF_ALPHA                            | Electron localisation function of    |
|                                      | the alpha density                    |
+--------------------------------------+--------------------------------------+
| ELF_BETA                             | Electron localisation function of    |
|                                      | the beta density                     |
+--------------------------------------+--------------------------------------+
| WAVEFUNCTION                         | Wavefunction instead of density is   |
|                                      | stored                               |
+--------------------------------------+--------------------------------------+
| HESSIAN                              | Approximate Hessian used in geometry |
|                                      | optimization                         |
+--------------------------------------+--------------------------------------+
| FINDIF                               | Positions and gradients for finite   |
+--------------------------------------+--------------------------------------+
| difference calculations              |                                      |
+--------------------------------------+--------------------------------------+
| VIBEIGVEC                            | Eigenvectors of Hessian              |
+--------------------------------------+--------------------------------------+
| MOLVIB                               | The matrix of second derivatives,    |
+--------------------------------------+--------------------------------------+
| as used by the program MOLVIB        |                                      |
+--------------------------------------+--------------------------------------+
| VIB1.log                             | Contains the modes 4 – 3N in a style |
+--------------------------------------+--------------------------------------+
| similar to the gaussian output for   |                                      |
+--------------------------------------+--------------------------------------+
| visualization with MOLDEN,           |                                      |
| MOLEKEL,...                          |                                      |
+--------------------------------------+--------------------------------------+
| VIB2.log                             | Contains the modes 1 – 3N-3          |
+--------------------------------------+--------------------------------------+
| ENERGYBANDS                          | Eigenvalues for each k points        |
+--------------------------------------+--------------------------------------+
| KPTS_GENERATION                      | Output of k points generation        |
+--------------------------------------+--------------------------------------+
| WANNIER_CENTER                       | Centers of the Wannier functions     |
+--------------------------------------+--------------------------------------+
| WC_SPREAD                            | Spread of the Wannier functions      |
+--------------------------------------+--------------------------------------+
| WC_QUAD                              | Second moments of the Wannier        |
|                                      | functions as :                       |
+--------------------------------------+--------------------------------------+
| ISTEP QXX QXY QXZ QYY QYZ QZZ        |                                      |
+--------------------------------------+--------------------------------------+
| IONS+CENTERS.xyz                     | Trajectory of ionic positions and    |
|                                      | WFCs in Å                            |
+--------------------------------------+--------------------------------------+
| WANNIER_DOS                          | Projection of the Wannier functions  |
+--------------------------------------+--------------------------------------+
| onto the Kohn-Sham states            |                                      |
+--------------------------------------+--------------------------------------+
| WANNIER_HAM                          | KS Hamiltonian in the Wannier states |
|                                      | representation                       |
+--------------------------------------+--------------------------------------+
| WANNIER_1.x                          | Wannier orbitals or orbital          |
|                                      | densities                            |
+--------------------------------------+--------------------------------------+
| wavefunctions                        | Contains the complex wavefunctions   |
|                                      | for Ehrenfest dynamics               |
+--------------------------------------+--------------------------------------+
| HARDNESS                             | Orbital hardness matrix              |
+--------------------------------------+--------------------------------------+
| RESTART.NMR                          | Files needed to restart a NMR/EPR    |
|                                      | calculation                          |
+--------------------------------------+--------------------------------------+
| RESTART.EPR                          |                                      |
+--------------------------------------+--------------------------------------+
| RESTART.L_x                          |                                      |
+--------------------------------------+--------------------------------------+
| RESTART.L_y                          |                                      |
+--------------------------------------+--------------------------------------+
| RESTART.L_z                          |                                      |
+--------------------------------------+--------------------------------------+
| RESTART.p_x                          |                                      |
+--------------------------------------+--------------------------------------+
| RESTART.p_y                          |                                      |
+--------------------------------------+--------------------------------------+
| RESTART.p_z                          |                                      |
+--------------------------------------+--------------------------------------+
| j\ :ma                               | Files in .CUBE format that contain   |
| th:`\alpha`\ \_B\ :math:`\beta`.cube | the                                  |
+--------------------------------------+--------------------------------------+
| B\ :ma                               | induced current densities and        |
| th:`\alpha`\ \_B\ :math:`\beta`.cube | magnetic fields                      |
+--------------------------------------+--------------------------------------+
| in an NMR calculation, respectively  |                                      |
| (:math:`\alpha,\beta`\ =x,y,z)       |                                      |
+--------------------------------------+--------------------------------------+
| PSI_A.\ :math:`i`.cube               | Files in .CUBE format that contain   |
|                                      | (spinpolarized)                      |
+--------------------------------------+--------------------------------------+
| PSI_B.\ :math:`i`.cube               | orbitals and densities.              |
+--------------------------------------+--------------------------------------+
| RHO_TOT.cube                         |                                      |
+--------------------------------------+--------------------------------------+
| RHO_SPIN.cube                        |                                      |
+--------------------------------------+--------------------------------------+
| QMMM_ORDER                           | Relation between the various         |
|                                      | internal atom order lists in QM/MM   |
|                                      | runs.                                |
+--------------------------------------+--------------------------------------+
| QM_TEMP                              | “local” temperatures of the QM/MM    |
|                                      | subsystems.                          |
+--------------------------------------+--------------------------------------+
| CRD_INI.grm                          | Positions of all atoms of first step |
|                                      | in Gromos format (QM/MM only).       |
+--------------------------------------+--------------------------------------+
| CRD_FIN.grm                          | Positions of all atoms of last step  |
|                                      | in Gromos format (QM/MM only).       |
+--------------------------------------+--------------------------------------+
| MM_TOPOLOGY                          | New Gromos format topology file      |
|                                      | (QM/MM only).                        |
+--------------------------------------+--------------------------------------+
| ESP                                  | Contains the ESP charges of the QM   |
|                                      | atoms (QM/MM only).                  |
+--------------------------------------+--------------------------------------+
| EL_ENERGY                            | Contains the electrostatic           |
|                                      | interaction energy (QM/MM only).     |
+--------------------------------------+--------------------------------------+
| MULTIPOLE                            | Contains the dipole- and the         |
|                                      | quadrupole-moment of the quantum     |
|                                      | system.                              |
+--------------------------------------+--------------------------------------+
| MM_CELL_TRANS                        | Contains the trajectory of the       |
|                                      | offset of the QM cell (QM/MM only).  |
+--------------------------------------+--------------------------------------+
| OCCMAT                               | The "trajectory" of the Occupation   |
|                                      | Matrix (DFT+U only)                  |
+--------------------------------------+--------------------------------------+

In case of path integral runs every replica :math:`s=\{1, \dots , P\}` gets its
own ``RESTART_``\ :math:`s`, ``RESTART_``\ :math:`s`\ ``.x``,
``DENSITY_``\ :math:`s`\ ``.x``, ``ELF_``\ :math:`s`\ ``.x``, and
``GEOMETRY_``\ :math:`s` file.

In contrast, in case of mean free energy path searches, each replica uses its
own directory for nearly all files. Exception are the file ``LATEST`` and the
geometry files named as for path integrals: ``GEOMETRY_``\ :math:`s`. The
directories are built from the :kw:`FILEPATH` path (or else from
the environment variable **CPMD_FILEPATH**) with the suffixes ``_``\ :math:`s`,
:math:`s=\{1, \dots , P\}`.

In case of bicanonical Ensemble MD

+-------------------+---------------------------------------------------------+
| *File*            | *Contents*                                              |
+-------------------+---------------------------------------------------------+
| OUTPUT_CNF2       | standard output of the smaller canonical system         |
+-------------------+---------------------------------------------------------+
| RESTART_CNF1      | Restart file of the larger canonical system             |
+-------------------+---------------------------------------------------------+
| RESTART_CNF2      | Restart file of the smaller canonical system            |
+-------------------+---------------------------------------------------------+
| RESTART_CNF1.x    | Old/New restart files                                   |
+-------------------+---------------------------------------------------------+
| RESTART_CNF2.x    | Old/New restart files                                   |
+-------------------+---------------------------------------------------------+
| LATEST_CNF1       | Info file on the last restart file                      |
+-------------------+---------------------------------------------------------+
| LATEST_CNF2       | Info file on the last restart file                      |
+-------------------+---------------------------------------------------------+
| ENERGIES          | Chemical potential, bicanonical weight, Kohn-Sham       |
|                   | energy difference along the trajectory                  |
+-------------------+---------------------------------------------------------+
| ENERGIES_CNF1     | All energies along the trajectory larger canonical      |
|                   | system                                                  |
+-------------------+---------------------------------------------------------+
| ENERGIES_CNF2     | All energies along the trajectory smaller canonical     |
|                   | system                                                  |
+-------------------+---------------------------------------------------------+
| GEOMETRY_CNF1     | Current ionic positions and velocities larger canonical |
|                   | system                                                  |
+-------------------+---------------------------------------------------------+
| GEOMETRY_CNF2     | Current ionic positions and velocities smaller          |
|                   | canonical system                                        |
+-------------------+---------------------------------------------------------+
| GEOMETRY_CNF1.xyz | Current ionic positions and velocities in Ålarger       |
|                   | canonical system                                        |
+-------------------+---------------------------------------------------------+
| GEOMETRY_CNF2.xyz | Current ionic positions and velocities in Åsmaller      |
|                   | canonical system                                        |
+-------------------+---------------------------------------------------------+
| TRAJECTORY_CNF1   | All ionic positions and velocities along the trajectory |
|                   | larger canonical system                                 |
+-------------------+---------------------------------------------------------+
| TRAJECTORY_CNF2   | All ionic positions and velocities along the trajectory |
|                   | smaller canonical system                                |
+-------------------+---------------------------------------------------------+
| TRAJEC_CNF1.xyz   | All ionic positions along the trajectory in xyz-format  |
|                   | larger canonical system                                 |
+-------------------+---------------------------------------------------------+
| TRAJEC_CNF2.xyz   | All ionic positions along the trajectory in xyz-format  |
|                   | smaller canonical system                                |
+-------------------+---------------------------------------------------------+
| TRAJEC_CNF1.dcd   | All ionic positions along the trajectory in dcd-format  |
|                   | larger canonical system                                 |
+-------------------+---------------------------------------------------------+
| TRAJEC_CNF2.dcd   | All ionic positions along the trajectory in dcd-format  |
|                   | smaller canonical system                                |
+-------------------+---------------------------------------------------------+

In general, existing files are **overwritten**!

Exceptions are “trajectory” type files (``TRAJECTORY``, ``ENERGIES``,
``MOVIE``, ``STRESS``, ...), in them data are **appended**.
