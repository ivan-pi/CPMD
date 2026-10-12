Restarts
========

General information
-------------------

All restart information for CPMD simulations are stored within one binary file.
There are very few exceptions we will discuss later. The name of the restart
files is ``RESTART`` or ``RESTART.n``, where :math:`n` stands for an integer
number. If the keyword :kw:`RESTART` is found the program processes
the file with the name ``RESTART``. Using suboptions to the
:kw:`RESTART` option, the information retained from the file can be
specified. For example the suboptions ``COORDINATES WAVEFUNCTION`` will force
the program to use the geometry and orbitals from the ``RESTART`` file.

At the end of a simulation or at regular intervals (using the keyword
:kw:`STORE`) a restart file with the default name ``RESTART.1`` is
written. If this happens more than once (e.g. during a molecular dynamics run)
the restart file is being overwritten. Using the keyword
:kw:`RESTFILE` it can be specified that more than one restart file
should be used for writing. If the :kw:`RESTFILE` parameter was set
to 4, then 4 restart files with the names ``RESTART.1``, :math:`\ldots`,
``RESTART.4`` will be written. If more than 4 restart files are needed the
first one will be overwritten. This option is useful if there is the
possibility that restart files get corrupted (e.g. on instable systems), or if
simulations are performed that might lead to unphysical results. In this case
it might be possible to go back to a restart file which contains still intact
information.

The name of the last restart file written is stored in the file ``LATEST``.
Using the suboption LATEST to the keyword RESTART changes the default name of
the file to be read from RESTART to the name found in the file ``LATEST``. The
danger of using this option is that the file from which the simulation is
started gets overwritten during the simulation. Using the default (starting
from RESTART) ensures that the original file stays intact. However, it requires
the renaming of the final file of a simulation from ``RESTART.n`` to
``RESTART``.

Typical restart scenarios
-------------------------

Wavefunction optimizations
~~~~~~~~~~~~~~~~~~~~~~~~~~

| The restart options used in wavefunction optimizations are
| **RESTART WAVEFUNCTION COORDINATES**. The suboption COORDINATES is not really
  necessary but it is advised to use it anyway, as in this way the
  correspondence of wavefunction and ionic geometry is assured.

Geometry optimizations
~~~~~~~~~~~~~~~~~~~~~~

| Typical suboptions used in a geometry optimizations are
| **RESTART WAVEFUNCTION COORDINATES HESSIAN**. With the suboption HESSIAN the
  information from previous runs stored in the updated approximate HESSIAN can
  be reused.

.. _molecular-dynamics-1:

Molecular dynamics
~~~~~~~~~~~~~~~~~~

| Molecular dynamics simulations use restart options of the kind
|  **RESTART WAVEFUNCTION COORDINATES VELOCITIES**. These are the minimal
  options needed for a smooth continuation of a Car–Parrinello molecular
  dynamics simulation. Use of the suboption ACCUMULATORS ensures that the
  calculated means (e.g. temperature) are correct for the whole simulation, not
  just the current run. If Nosé thermostats are used it is important also the
  restart the thermostat variables. This is achieved by adding the
  corresponding keywords to the RESTART (NOSEE, NOSEP, NOSEC).

Kohn–Sham energies
~~~~~~~~~~~~~~~~~~

The calculation of canonical Kohn–Sham orbitals requires a restart. In general,
this will be a restart from converged orbitals from a wavefunction
optimization. There is no way that the program can check this. However, if the
same convergence criteria are used, the number of occupied states orbitals
should converge in the first iteration of the diagonalization.

Some special cases
------------------

The suboption VELOCITIES will result in a restart from both, ionic and
wavefunction velocities. In special cases, this is not the desired behavior.
Using the additional keyword :kw:`QUENCH` the read velocities can be
set back to zero. This will be most likely used for wavefunctions with QUENCH
ELECTRONS. Another possibility is to reoptimize the wavefunction at the start
of a molecular dynamics simulation. This is achieved with the keywords QUENCH
BO.

For performance reasons the writing of the restart file should be done only
occasionally. This might cause problems if the simulation was terminated
incorrectly. Several hundreds or thousands of simulation steps might be lost.
For this reason CPMD writes a special output file ``GEOMETRY`` after each
molecular dynamics step. Together with a normal restart file this allows to
start the simulation form the last ionic configuration and velocities. To
achieve this another suboption GEOFILE has to be added to the RESTART keyword.
After reading the positions and velocities of the ions from the restart file,
they are also read from the GEOMETRY file and overwritten.

Special restarts to be used with the keywords :kw:`TDDFT` and
:kw:`VIBRATIONAL ANALYSIS` are discussed in the
sections covering that type of simulations.
