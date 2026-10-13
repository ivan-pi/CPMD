.. _hints:roks:

Restricted Open-Shell Calculations
==================================

Molecular dynamics simulations in the first excited state can be performed
using Restricted Open-Shell Kohn-Sham (ROKS) theory [Frank98]_.
The keyword :kw:`ROKS` in the &CPMD section defaults to the first
excited singlet state. Solving open-shell equations is not simple unless

#. a high-spin state is computed.

#. the two singly occupied molecular orbitals (SOMOs) have different spatial
   symmetry.

In these two cases the Goedecker-Umrigar-Algorithm (GOEDECKER) may be used
which shows the best convergence properties and is applicable in connection
with Car-Parrinello molecular dynamics. Otherwise it is necessary to use a
modified variant of the Goedecker-Umrigar-Algorithm and to do Born-Oppenheimer
molecular dynamics (unless you know what you are doing). In almost all cases,
the default algorithm (DELOCALIZED) is applicable, whereas for example some
dissociation reactions require the localized variant to enable localization of
the orbitals on the fragments.

::

         ROKS LOCALIZED

In order to make sure that the chosen algorithm works for a certain system, the
conservation of energy during a molecular dynamics simulation and the shape of
the orbitals should always be checked. One of the SOMOs should have the same
nodal structure as the HOMO obtained by a ground state calculation. If using
the unmodified Goedecker-Umrigar scheme (GOEDECKER), the energy of the singlet
may collapse to approximately the triplet energy if the two SOMOs do not have
different symmetries. The triplet energy can be calculated by specifying

::

          ROKS TRIPLET

or also

::

          ROKS TRIPLET GOEDECKER

See the description of the keywords :kw:`LOW SPIN EXCITATION`, :kw:`LSE PARAMETERS`
and :kw:`MODIFIED GOEDECKER` for a description of how to
do ROKS calculations using the older input LOW SPIN EXCITATION ROKS. ROKS
GOEDECKER corresponds to LOW SPIN EXCITATION ROKS whereas ROKS DELOCALIZED
corresponds to LOW SPIN EXCITATION ROKS with MODIFIED GOEDECKER. Do not use LOW
SPIN EXCITATION in the &SYSTEM section and ROKS in the &CPMD section at the
same time.

ROKS is not implemented with Vanderbilt pseudopotentials.

A Slater transition-state density between a singlet ground state and the first
excited singlet state (or any pair of states described with ROKS) can be useful
whenever one set of Kohn-Sham states is required which is equally well suited
for each of the states involved in a transition, e.g., to calculate the
couplings between the electronic transition and an external influence. This
method is analogous to state-averaged multiconfigurational SCF methods and
shares many of their benefits with them. In CPMD, it can be used to calculate
non-adiabatic couplings between singlet
states [nonadiabatic]_ [lsets]_, see options
:kw:`COUPLINGS`.
