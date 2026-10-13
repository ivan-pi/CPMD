Molecular Dynamics
==================

.. _hints:nose:

Choosing the Nosé-Hoover chain thermostat parameters
----------------------------------------------------

The Nosé-Hoover chain thermostat is defined by specifying three parameters: A
target kinetic energy, a frequency and a chain length. For the ions, given the
target temperature :math:`T_W`, the target kinetic energy is just
:math:`gkT_W`, where :math:`g` is the number of degrees of freedom involved in
a common thermostat. For example, if there is one thermostat on the entire
ionic system, then :math:`g=3N_{AT}-N_{const}`, where :math:`N_{const}` is the
number of constraints to which the atoms are subject. The frequency for the
ionic thermostat should be chosen to be some characteristic frequency of the
ionic system for which one wishes to insure equilibration. In water, for
example, one could choose the O-H bond vibrational frequency. (Having a precise
value for this frequency is not important, as one only wishes to insure that
the thermostat will couple to the mode of interest.) The choice of chain length
is not terribly important as it only determines how many extra thermostats
there will be to absorb energy from the system. Usually a chain length of 4 is
sufficient to insure effective equilibration. Longer chains may be used in
situations where heating or cooling effects are more dramatic.

For the electrons, the target kinetic energy is not usually known *a priori* as
it is for the ions. However, by performing a short run without thermostats, one
can determine a value about which the electron kinetic energy ‘naturally’
fluctuates and take this as the target value. While the precise value is not
important, a little experience goes a long way, as a choice that is either too
small or too large can cause spurious damping of the ions or departures from
the Born-Oppenheimer surface, respectively. A good choice for the frequency of
the electron thermostat can be made based on :math:`\Omega_I^{\rm
max}`, the maximum frequency in the phonon spectrum. The frequency of the
electron thermostat should be at least 2-3 times this value to avoid coupling
between the ions and the electron thermostats. As an example, for silicon, the
highest frequency in the phonon spectrum is 0.003 a.u., so a good choice for
the electron thermostat frequency is 0.01 a.u. The chain length of the electron
thermostat can be chosen in the same way as for the ions. 4 is usually
sufficient, however longer chains may be used if serious heating is expected.
In addition, the electron thermostats have an extra parameter that scales the
number of dynamical degrees of freedom for the electrons. (:math:`1/\beta_e =
2E_e/N_e`, where :math:`E_e` is the desired electron kinetic energy and
:math:`N_e` is the number of dynamical degrees of freedom for the electrons –
see Eq. (3.4) in Ref.\ [Tuckerman94]_). The default value is the
true number of dynamical degrees of freedom
:math:`N_e = (2*N_{GW}-1)*N_{ST} - N_{ST}^p`, where :math:`p=2` for
orthonormality constraints and :math:`p=1` for norm constraints. When this
number is very large, it may not be possible to integrate the electron chain
thermostats stably using a frequency above that top of the phonon spectrum.
Should this be the case in your problem, then the number of dynamical degrees
of freedom should be scaled to some smaller number such that the system can
once again be integrated stably. This parameter has no other effect that to
change the relative time scales between the first element of the electron
thermostat chain and the other elements of the chain.

In addition to the basic parameters defining the chains themselves, one needs
to specify two more parameters related to the integration of the thermostated
equations of motion. The first is the order :math:`M_{SUZ}` of the Suzuki
integrator. Experience shows that the choice :math:`M_{SUZ}=3` is sufficient
for most applications. Finally, one must specify the number of times the Suzuki
integrator will be applied in a given update. This is the parameter
:math:`N_{SUZ}` which determines the basic Suzuki time step :math:`\delta t` =
:math:`\Delta t/N_{SUZ}`, where :math:`\Delta t` is the time step being used in
the MD run. :math:`N_{SUZ} =2` or 3 is usually large enough to give stable
integration. If more stable integration is required, try :math:`M_{SUZ}=4` or
make :math:`N_{SUZ}` larger.

Fictitious electron mass, time step and energy conservation
-----------------------------------------------------------

In Car-Parrinello molecular dynamics the fictitious electron mass
(:kw:`EMASS`) and the time step (:kw:`TIMESTEP`) have to be chosen together.
The fictitious kinetic energy of the electrons (``EKINC`` in the output)
should remain small compared with the kinetic energy of the ions and should
not drift during the run; a larger electron mass increases ``EKINC`` and
allows a longer time step, a smaller electron mass keeps the electrons closer
to the Born-Oppenheimer surface but requires a shorter time step. The quality
of a run is judged from the conservation of the constant of motion
(``ECLASSIC``, or ``EHAM`` when thermostats are used) and from the behavior
of ``EKINC``.

.. figure:: ../fig/md-ekinc.*
   :width: 80%
   :alt: Fictitious electronic kinetic energy for different electron masses and time steps

   Fictitious kinetic energy of the electrons, ``EKINC``, during
   Car-Parrinello runs with Troullier-Martins (MT, 70 Ry cutoff) and
   Vanderbilt (VDB, 35 Ry cutoff) pseudopotentials for different electron
   masses :math:`\mu` (a.u.) and time steps :math:`\Delta t` (a.u.). ``EKINC``
   increases with the electron mass; for the largest masses it also drifts
   upwards during the run.

.. figure:: ../fig/md-timestep-eham.*
   :width: 80%
   :alt: Conserved energy for time steps from 1 to 5 a.u. and the corresponding CPU time

   Constant of motion ``EHAM`` (shifted and scaled) during 1 ps of
   Car-Parrinello dynamics with time steps from 1 to 5 a.u., and the CPU time
   needed for each run (bars). Longer time steps reduce the cost considerably
   while the fluctuations of the conserved energy stay small.

.. figure:: ../fig/md-econs.*
   :width: 60%
   :alt: Conserved energy as a function of time for three pairs of simulations

   Conserved energy :math:`E_\mathrm{cons}` as a function of time for three
   pairs of simulations (open and filled symbols in each panel). In each
   panel one of the two runs shows a steady drift of the conserved energy.

   .. todo::

      MISSING: the simulations compared in this figure (system, method and
      parameters of the runs shown with open and filled symbols) are not
      documented, and neither are the details (system, functional,
      pseudopotentials) of the two figures above.
