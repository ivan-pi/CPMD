.. _hints:wfconv:

Wavefunction Convergence
========================

Some general comments on wavefunction optimizations:

Any optimization that takes more than 100 steps should be considered slow.

Optimizations using ODIIS that have repeated resets (more than a few) will
probably never converge. Convergence for LSD is normally slower and more
difficult than for unpolarized cases.

If the ODIIS convergence gets stuck (more than one reset) stop the run and
restart using

::

     PCG MINIMIZE
     TIMESTEP
      20

The conjugate gradient minimizer with line search is much more robust. For LSD
and larger systems it should be used from the start.

A typical behavior will be that after the restart the energy goes down and the
gradient increases. This means that we are in a region where there are negative
curvatures. In such regions the DIIS minimizer moves in the wrong direction.
After some iterations we will be back to normal behavior, energy and gradient
get smaller. At this point it may be save to switch back to ODIIS.

Sometimes, it can also be helpful to wait longer for a DIIS reset and to
diagonalize after repeated resets to get out of this region. This can be
accomplished using

::

     ODIIS NO_RESET=20
       5
     LANCZOS DIAGONALIZATION RESET=2

Starting a Car-Parrinello MD from a random wavefunction with all atom positions
fixed, a comparatively high electron mass and using
:kw:`ANNEALING` ELECTRONS is another alternative to get to a
reasonably converged wavefunction. Due to the exponential convergence of the
annealing procedure, one should switch to a different optimizer as soon as the
fictitious kinetic energy of the electrons drops below the typical range for a
(normal) MD run.

Wavefunction optimizations for geometries that are far from equilibrium are
often difficult. If you are not really interested in this geometry (e.g. at the
beginning of a geometry optimization or this is just the start of a MD) you can
relax the convergence criteria to :math:`10^{-3}` or :math:`10^{-4}` and do
some geometry steps. After that optimization will be easier.

Some general remarks on comparing the final energies:

Converge the wavefunction very well, i.e. set :kw:`CONVERGENCE`
ORBITALS to :math:`10^{-6}` or better.

Make sure that all parameters are the same:

- same geometry,

- same functional,

- same number of grid points (this may differ if you use different FFT
  libraries)

- same number of spline points for the PP

(IMPORTANT: the default for :kw:`SPLINE` **POINTS** has changed
between different CPMD versions, :math:`500 \to 3000 \to 5000`). A very good
test is to start always from the same RESTART file and only do one single step.
This way ALL energies have to be exactly the same and no problems with
different convergence rates occur.
