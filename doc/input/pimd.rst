.. _input-pimd:

``&PIMD`` section
=================

Path integral molecular dynamics (PIMD). This section is only evaluated if the :kw:`PATH INTEGRAL` keyword is given in the ``&CPMD`` section.

The section starts with ``&PIMD`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &PIMD

Keywords
--------

.. keyword:: CENTROID DYNAMICS
   :section: &PIMD

   Adiabatic centroid molecular dynamics, see
   Ref. [Cao93]_ [Martyna96]_ [aicmd]_ for theory and details of our
   implementation, which yields quasiclassical dynamics of the nuclear centroids
   at a specified temperature of the non–centroid modes. This keyword makes only
   sense if used in conjunction with the normal mode propagator via the keyword
   NORMAL MODES *and* FACSTAGE :math:`>1.0` *and* WMASS :math:`=1.0`. The centroid
   adiabaticity control parameter FACSTAGE, which makes the non-centroid modes
   artificially fast in order to sample adiabatically the quantum fluctuations,
   has to be chosen carefully; note that FACSTAGE :math:`= 1/\gamma` as introduced
   in Ref. [aicmd]_ in eq. (2.51). The abbreviation ``CMD`` is also
   accepted.

.. keyword:: CLASSICAL TEST
   :section: &PIMD

   Test option to reduce the path integral branch to the classical code for the
   special case :math:`P=1` in order to allow for a one-to-one comparison to a run
   using the standard branch of CPMD. It works only with primitive propagator,
   i.e. not together with NORMAL MODES, STAGING and/or
   :kw:`DEBROGLIE` CENTROID.

.. keyword:: DEBROGLIE
   :section: &PIMD
   :alternatives: [CENTROID]

   An initial configuration assuming quantum free particle behavior is generated
   for each individual atom according to its physical mass at the temperature
   given in Kelvin on the following input line. Using DEBROGLIE each nuclear
   position obtained from the &ATOMS … &END section serves as the starting point
   for a Gaussian Lévy walk of length :math:`P` in three dimensions, see
   e.g. Ref. [Fosdick66]_. Using DEBROGLIE CENTROID each nuclear
   position obtained from the &ATOMS … &END section serves as the centroid (center
   of geometry) for obtaining the centroid (center of geometry) for obtaining the
   :math:`P` normal modes in three dimensions, see
   e.g. Ref. [Tuckerman96]_. This option does only specify the
   generation of the initial configuration if INITIALIZATION and GENERATE REPLICAS
   are active. Default is DEBROGLIE CENTROID and 500 Kelvin. ``CENTROIDS``
   is accepted as an alias of ``CENTROID``.

.. keyword:: FACMASS
   :section: &PIMD

   Obtain the fictitious nuclear masses :math:`M_I^\prime` within path integral
   molecular dynamics from the real physical atomic masses :math:`M_I` (as
   tabulated in the DATA ATWT / …/ statement in atoms_utils.mod.F90) by
   *multiplying* them with the dimensionless factor WMASS that is read from the
   following line; if the NORMAL MODES or STAGING propagator is used obtain
   :math:`M_I^{\prime (s)}= \mbox{WMASS} \cdot M_I^{(s)}` for *all* replicas
   :math:`s=1, \dots , P`; see e.g. Ref. [aicmd]_ eq. (2.37) for
   nomenclature.

   **Default** value of WMASS is **1.0**

.. keyword:: GENERATE REPLICAS
   :section: &PIMD

   Generate quantum free particle replicas from scratch given a classical input
   configuration according to the keyword :kw:`DEBROGLIE`
   specification. This is the default if :kw:`INITIALIZATION` is
   active.

.. keyword:: GLE_LAMBDA
   :target: GLE LAMBDA
   :section: &PIMD

   Set the scaling factor :math:`\lambda` of the generalized Langevin thermostat
   read from the next line for removing the resonances between the vibrations of
   the system and the harmonic potential representing the quantum kinetic energy
   term in the description of the path integral. The default value is 0.5 as
   suggested in Ref. [Rossi2014]_.

.. keyword:: INITIALIZATION
   :section: &PIMD

   Provide an initial configuration for all replicas as specified either by
   :kw:`GENERATE REPLICAS` or by :kw:`READ REPLICAS`. This option is automatically activated if
   :kw:`RESTART` COORDINATES is not specified. It is defaulted to
   GENERATE REPLICAS together with :kw:`DEBROGLIE` CENTROID and a
   temperature of 500 Kelvin.

.. keyword:: NORMAL MODES
   :section: &PIMD

   Use the normal mode representation [Tuckerman96]_ of the path
   integral propagator. It is possible to impose a mass disparity between centroid
   and non–centroid coordinates by dividing the fictitious masses of only the
   *non*–centroid :math:`s=2, \dots ,P` replicas by the adiabaticity control
   factor FACSTAGE. This dimensionless factor *must always* be specified in the
   following line. Note: the eigen–*frequencies* of the :math:`s>1` replicas are
   changed by only :math:`\sqrt{\mbox{FACSTAGE}}`, see
   Ref. [Martyna96]_\ (b). Using FACSTAGE :math:`\not= 1.0` makes
   only sense in conjunction with CENTROID DYNAMICS where WMASS=1.0 has to be used
   as well.

.. keyword:: OUTPUT
   :section: &PIMD
   :alternatives: [ALL, GROUPS, PARENT]

   Output files for each processor, processor group, or only grandparent. Default
   is PARENT to standard output file (Note: some information such as messages for
   correct reading / writing of restart files is lost); GROUPS and ALL write to
   the files OUTPUT\_\ :math:`n` where :math:`n` is the group and bead number,
   respectively. ``GROUP`` is accepted as an alias of ``GROUPS``.

.. keyword:: PRINT LEVEL
   :section: &PIMD

   The detail of printing information is read as an integer number from the next
   line. Currently there is only minimal output for :math:`<5` and maximal output
   for :math:`\geq 5`.

.. keyword:: PROCESSOR GROUPS
   :section: &PIMD

   This is only needed for *fine*–tuning load balancing in case of path integral
   runs *iff* two level parallelization is used. The default optimizes the
   combined load balancing of the parallelization over replicas and g–vectors. The
   default load distribution is usually optimal. Separate the total number of
   processors into a certain number of processor groups that is read from the
   following line; only 2\ :math:`^N` = 2, 4, 8, 16, :math:`\dots` groups are
   allowed and the maximum number of groups is the number of replicas. Every
   processor group is headed by one PARENT and has several CHILDREN that work
   together on a single replica at one time; the processor groups work
   sequentially on replicas if there is more than one replica assigned to one
   processor group. Note: if the resulting number of processor groups is much
   smaller than the number of replicas (which occurs in “odd” cases) specifying
   the number of processor groups to be equal to the number of replicas might be
   more efficient. This keyword is only active in parallel mode.

.. keyword:: READ REPLICAS
   :section: &PIMD

   Read all :math:`P` replicas from a file with a name to be specified in the
   following line, for the input format see subroutine rreadf_utils.mod.F90.

.. keyword:: RING-POLYMER DYNAMICS
   :section: &PIMD

   Perform ring-polymer molecular dynamics (RPMD) [Craig2004]_. In
   the present implementation the fictitious mass of the nuclei for all the
   replicas is set to the real physical mass divided by the number of replicas (P)
   in the primitive representation and the physical mass in the normal mode
   representation for compatibility with the formulation of the conventional path
   integral. According to our choice of the fictitious masses, the target
   temperature you should set is not TP in the original formulation of
   Ref. [Craig2004]_ but just T as in the standard path integral and
   centroid MD. Note that RPMD in the staging mode representation is not supported
   and the generalized Langevin thermostat is available in the normal mode
   representation only. The abbreviation ``RPMD`` is also accepted.

.. keyword:: STAGING
   :section: &PIMD

   Use the staging representation [Tuckerman96]_ of the path integral
   propagator. It is possible to impose a mass disparity between centroid and
   non–centroid coordinates by dividing the fictitious masses of only the
   *non*–centroid :math:`s=2, \dots ,P` replicas by the adiabaticity control
   factor FACSTAGE. This dimensionless factor *must always* be specified in the
   following line. Note: the eigen–*frequencies* of the :math:`s>1` replicas are
   changed by only :math:`\sqrt{\mbox{FACSTAGE}}`, see
   Ref. [Martyna96]_\ (b). Note: using FACSTAGE :math:`\not= 1.0`
   essentially makes no sense within the STAGING scheme, but see its use within
   CENTROID DYNAMICS and NORMAL MODES.

.. keyword:: TROTTER DIMENSION
   :section: &PIMD

   The Trotter number :math:`P`, i.e. the number of “replicas”, “beads”, or
   “imaginary time slices” which are used in order to discretize the Feynman–Kac
   path integral of the nuclei, is read from the next line. If NORMAL MODES or
   STAGING is not activated the path integral is discretized in cartesian
   coordinates in real space (so–called “primitive coordinates”). A discussion
   about controlling discretization errors and on estimating :math:`P` in advance
   is given in Ref. [knoll-marx-00]_.
