Hints on using FEMD
===================

There are several parameters which crucially affect the speed, accuracy and
robustness of the FEMD method. These are related to: LANCZOS PARAMETERS, STATES
and ANDERSON MIXING. Less crucially, the ELECTRON TEMPERATURE.

.. _hints:lanczos:

Lanczos Parameters
------------------

Several parameters related to the Lanczos ( Friesner-Pollard) method are given.
Generically:

::

                LANCZOS PARAMETER [N=n]
                   ncycle  nkrylov  nblock  tolerance
                   drhomax(2)    tolerance(n)
                     .....
                     .....
                   drhomax(n)    tolerance(n)

:math:`Ncycle` can always be safely set to 50. Similarly, :math:`Nkrylov=8` is
almost always a good choice. Exceptionally, for certain d-metallic systems,
increasing :math:`nkrylov=16` may be more efficient. :math:`Nblock` is the
dimension of the blocking in the evaluation of
:math:`H[\psi_1,...,\psi_{nblock}]`. :math:`Nblock` should be a divisor of
:math:`NSTATE` and recommended values lie in the range of 20-100. The
:math:`tolerance` specifies the accuracy to be achieved in the Lanczos method.
States are considered converged if

.. math:: |H\psi-\epsilon\psi|^2 < tolerance

For efficient calculations, the tolerance should vary according to closeness to
self-consistency ( as measured by DRHOMAX). During initial stages of the SC
cycle, the tolerance can be loose, gradually tightening until close to SC it is
high. An example of this might be:

::

                LANCZOS PARAMETER N=5
                   50  8  20 1.D-9
                   0.05      1.D-11
                   0.01      1.D-13
                   0.0025    1.D-16
                   0.001     1.D-18

For accurate forces, a final tolerance of at least 1.D-16 is recommended,
although accurate energies can be got using a lower tolerance. It is worth
experimenting how best to tighten the tolerance - it could be system dependent.

Other important FEMD parameters
-------------------------------

The keyword STATES defines the dimension of the subspace used in the
diagonalization. STATES must be greater than or equal to :math:`N_{el}/2`, but
it is generally good to allow for a number of more or less empty bands (usually
10% or so). Finally, ANDERSON MIXING determines the rate of convergence to
self-consistency. Properly chosen the convergence can be very fast. Typically
for bulk systems we use values between 0.2-0.5, smaller values being necessary
for large systems. For metallic surfaces, small values are necessary (
typically 0.03-0.05).

If using k-points, then it is usually a good idea ( and this is done by default
if using MONKHORST PACK k-points) to exploit symmetries. In this case, however,
beware of including the POINT GROUP keyword to symmetrise properly the density.
Finally, if starting from a high-symmetry structure, you may nevertheless want
to use the full k-point mesh ( apart from time-inversion symmetry related
k-points), and in this case specify the keyword FULL.
