TDDFT
=====

The TDDFT part of CPMD is rather new. Therefore it hasn’t yet reached the
stability of other parts of the code. It has to be used with special care.

There are four different type of calculations that can be performed using the
TDDFT module; calculation of the electronic spectra, geometry optimization and
vibrational analysis, and molecular dynamics in excited states.

All options (spectra and forces, spin polarized and unpolarized) are
implemented for the Tamm–Dancoff approximation to TDDFT. Only part of these
options are available for the full TDDFT response calculation.

.. _sec:electronicspectra:

Electronic spectra
------------------

Electronic excitation energies can be calculated using the keyword :kw:`ELECTRONIC SPECTRA` in the &CPMD section. This calculation is
performed in three parts. First, the ground state wavefunctions are optimized,
then a limited set of unoccupied orbitals is determined and finally the TDDFT
response equations are solved. A typical input for such a calculation would
look like

::

   &CPMD
     ELECTRONIC SPECTRA
     DIAGONALIZATION LANCZOS
     COMPRESS WRITE32
   &END

   &TDDFT
     STATES SINGLET
       5
     TAMM-DANCOFF
     DAVIDSON PARAMETER
     150 1.D-7 50
   &END

For this calculation of the electronic spectra defaults are used for the ground
state optimization (ODIIS and :math:`10^{-5}` convergence). The calculation of
the empty states is performed using the Lanczos diagonalizer with default
settings. The final wavefunction will be stored in the restart file using 32
bit precision.

Five single states with the Tamm–Dancoff approximation have to be calculated.
The parameters for the Davidson diagonalization have been changed to 50 for the
Davidson subspace and a convergence criteria of :math:`10^{-7}` is used.

Restarting this type of calculation has to be done with care. At the end of
each phase of the calculation a new restart file is written. If the defaults
are used, each time the file ``RESTART.1`` is overwritten. For a restart from
converged ground state wavefunctions and canonical Kohn–Sham orbitals a restart
with

``RESTART WAVEFUNCTION COORDINATES``

will be used. A restart also including the linear response orbitals will use

``RESTART WAVEFUNCTION COORDINATES LINRES``.

In this case only restarts from the file ``RESTART`` are possible as after
phase one and two the file ``RESTART.1`` would be overwritten and the
information on the linear response orbitals, read only in phase three, would be
lost.

.. _sec:TDDFTdynamics:

Geometry optimizations and molecular dynamics
---------------------------------------------

Geometry optimizations and molecular dynamics simulations can only be performed
after an electronic spectra calculation. A typical input file would contain the
sections

::

    &CPMD
      OPTIMIZE GEOMETRY
      TDDFT
      RESTART WAVEFUNCTION COORDINATES LINRES
    &END

    &TDDFT
     STATES SINGLET
       1
     TAMM-DANCOFF
     DAVIDSON PARAMETER
     150 1.D-7 50
     FORCE STATE
      1
    &END

The keywords in section &CPMD are all mandatory. The section &TDDFT specifies
that the optimization should be performed for the first excited singlet state.
Replacing :kw:`OPTIMIZE GEOMETRY` by :kw:`MOLECULAR DYNAMICS` BO would result in a molecular dynamics
simulation. In this case further input specifying the time step, maximal number
of steps, thermostats, etc. would also be supplied.
