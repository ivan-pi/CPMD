.. _sec:meta:

Metadynamics
============

These are some notes about the use of the metadynamics (MTD) machinery within
CPMD. It is just a first version of a manual that I hope will be improved by
the comment and possibly the contributions of the users of this method.

The metadynamics can run in a standard NVE/NVT MD run or in a NPE/NPT run
(variable cell). In order to apply the MTD algorithms in CPMD some (not few)
lines have to be added in the input file. These lines are to be in the &ATOMS
.. &END section and they provide information about the kind of MTD to be
performed, the choice of collective variables (CV), some parameters required to
determine the time dependent potential and some other options. All the MTD
input must be between an initial and a final line which are:

``METADYNAMICS``
:math:`\cdots`
``END METADYNAMICS``

If the initial line contains also the keyword :math:`COLLECTIVE`
:math:`VARIABLES`, the standard MTD, with one single set of CV, is initialized.
If, instead, the keyword :math:`MULTI` is found, more than one MTD are
performed simultaneously on the same system; therefore, the whole set of CV is
constituted by :math:`NSUBSYS` subsets, which are independent one from each
other. The number of subsystems is given on the same line by writing
:math:`NS=` followed by an integer number (default: 1). Alternatively, if
:math:`MULTIPLE` :math:`WALKERS` is present in the same line multiple walker
metadynamics is performed using the extended Lagrangian metadynamics; number of
walkers is read in the same line immediately after :math:`NW=` (see
:ref:`sect:mw`). Instead, if :math:`CELL  FULL` is the keyword, the CV
are the 6 cell parameters (3 side lengths and 3 angles), and the MTD is
performed without extended Lagrangian, i.e. the contribution coming from
:math:`V(t)` is directly added into the stress tensor (see below in **MTD
Algorithm**).

For almost all the input parameters there is a reasonable default value, but,
since the range of applications of MTD is quite wide, it is likely that the
default values do not fit your problem. Therefore some effort is required to
choose the optimal conditions for your run. Of course, it is important to know
something about MTD before using it. There are some references about the method
[alaio]_ [iannuzzi]_ [micheletti]_ [gerlaio]_ [michele]_, and about some
successful applications, as e.g.
[andras]_ [gervasio]_ [iannuzzi2]_ [sergey]_ [scwo]_ [ikeda]_ [rna]_ [Nair-jacs-08]_ [mb11]_.
It can be of great help to read about the choices and results obtained by other
users. But I remark that there are very few general rules that can be held
valid for different problems and systems.

The method is based on the definition of a manifold of CV as functions of the
degrees of freedom characterizing your system,
:math:`{\bf S} = \{S_{\alpha}({\bf R},{\bf \phi},{\bf h})\}`, where
:math:`{\bf R}` are the ionic degrees of freedom, :math:`{\bf \phi}` are the
electronic wavefunctions, and :math:`{\bf h}` defines the cell box. The CV
which are implemented in the code, have been chosen according to the needs of
those who used the method up to now. Of course they do not exhaust all the
problems, and many more CV might be needed in the future. To implement them,
once the analytical formula and its derivatives are available, is not
complicated at all. In principle, the implementation should be easy for anybody
who knows a bit the CPMD code.

MTD Algorithm
-------------

Once the CV have bee chosen, the MTD method can be applied in two different
fashions.

**Direct MTD:** The simplest approach is to define the time dependent potential
as function of :math:`{\bf S}`, :math:`V(t,{\bf S})`, and apply it directly
onto the involved degrees of freedom. In this case, the equations of motion of
the dynamic variables of the system, :math:`{\bf R},{\bf \phi},{\bf h}`, will
include an additional term in the total forces, due to the contribution of
:math:`V(t,{\bf S})`. The disadvantage of this simplified version is that there
is scarce control on the dynamics in the space defined by the CV (CV-space),
which is a projection of the space of all the possible configurations. In
general, we would like to span thoroughly the CV-space, and to acquire
information about the underlying potential. Often, this means that we need a
slow dynamics in this space, where, for each set of values of the CV, we allow
the system to equilibrate and to choose the configuration with the highest
probability. Only in this way we will be able to construct a reasonable
probability distribution in the configurational space that has been explored
and consequently we will be able to reproduce the Free Energy surface.

**Lagrangian MTD:** This formulation is based on the method of the extended
Lagrangian. In addition to the dynamic variables that characterize your system,
a new set of variables :math:`{\bf s}=\{s_{\alpha}\}` is introduced. Each
:math:`s_{\alpha}` is associated to one of the selected :math:`S_{\alpha}`, it
has a fictitious mass :math:`M_{\alpha}` and velocity :math:`\dot{s}_{\alpha}`.
The equations of motion for the :math:`s_{\alpha}` variables are derived by a
properly extended Lagrangian, where we add the fictitious kinetic energy and
the potential energy as a function of :math:`{\bf s}`. Therefore the total
potential energy includes two new terms, a sum of harmonic potentials, which
couple the :math:`s_\alpha` to the respective
:math:`S_{\alpha}({\bf R},{\bf \phi},{\bf h})`,
:math:`\sum_{\alpha}k_{\alpha}(S_{\alpha}(\cdots)-s_{\alpha})^{2}`, and the
time dependent potential, which now is a function of :math:`{\bf s}`,
:math:`V(t,{\bf s})`. The coupling constants :math:`\{k_{\alpha}\}` and the
fictitious masses :math:`\{M_{\alpha}\}` are the parameters that determine the
dynamics of the :math:`\{s_{\alpha}\}` in the CV-space. Please notice that the
units of :math:`k` are Hartree divided by the square power of u.s., the
characteristic units of the specific CV (if CV is a distance it will be
:math:`a.u.^2`, if an angle :math:`radiants^2`, etc.). In analogy, the units of
the fictitious mass are :math:`Hartree ((t)/(u.s.))^2`, where :math:`t`
indicates the unit of time. Some guide lines on the choice of these parameters
will be given in the following paragraphs. By choosing the temperature
:math:`T_{\bf s}`, the velocities of the components of :math:`\bf s` can be
initialized giving via a Boltzmann distribution. Moreover, the velocities can
be kept in a desired range by the activation of a temperature control algorithm
(at the moment only the rescaling of velocity is implemented).

The Shape of :math:`V(t)`
-------------------------

Several shapes have been tested (and more might be proposed in the future). The
default choice is the construction of :math:`V(t)` by the accumulation of
Gaussian-like hills, i.e. (within the Lagrangian formulation, but the
expressions are the same for the direct MTD approach, providing to exchange
:math:`\bf s` with :math:`{\bf S(\cdots)}`)

.. math::
   :label: eq-hills

   \begin{aligned}
   V(t,{\bf s}) &=  \sum_{t_{i} < t } \Bigg[W_{i}\exp \left \{-
   \frac{({\bf s}-{\bf s}^{i})^{2}}{2 (\Delta s^{\perp})^{2}} \right \} \nonumber \\
    &  \exp \left
   \{-\frac{\left( ({\bf s}^{i+1}-{\bf s}^{i}) \cdot ({\bf s}-{\bf s}^{i})\right)^{2}}
   {2 (\Delta s^{||}_{i})^4} \right\}\Bigg],
   \end{aligned}

Here, :math:`t` indicates the actual simulation time, :math:`i` counts the
metadynamics steps, the first exponential gives the hill’s shape in the
direction perpendicular to the trajectory, whereas the second exponential tunes
the shape along the trajectory. In this form, the width of the hill along the
trajectory is determined by the displacement in the CV-space, walked between
two consecutive metadynamics steps,
:math:`\Delta s^{||}_{i} = f_{b}\sqrt{\Big[\sum_{\alpha}(s_{\alpha}^{i+1}-(s_{\alpha}^{i})^2\Big]}`.
:math:`f_b` is a factor, which is read in input and can be used to change the
size of the hills along the trajectory, by default it is 1. The height
:math:`W` and the width :math:`\Delta s^{\perp}` are input parameters that can
also be tuned during the MTD, in order to better fit the hill shape to the
curvature of the underlying energy surface (in the CV-space). As a rule of
thumb, :math:`\Delta s^{\perp}` should have roughly the size of the
fluctuations of CV at equilibrium (half the amplitude of the well) and
:math:`W` should not exceed few percents of the barrier’s height. These
information can be obtained by some short MD runs at equilibrium (without MTD)
and from some insight in the chemical/physical problem at hand. Since, in
general, different CV fluctuate in wells of different size, it is important to
define one scaling factor :math:`scf_{\alpha}` for each component
:math:`s_{\alpha}`, so that
:math:`\langle \delta s_{\alpha}\rangle /scf_{\alpha} = \Delta s^{\perp} \, \forall \alpha`.

Other implemented shapes of :math:`V(t)` are:

**Shift:** the tails of the Gaussians are cutoff, by setting to zero the
Gaussian at a distance :math:`R_{cutoff}\Delta s^{\perp}` from its center. In
this way the problem of the overlap of the tails in regions far from the actual
trajectory is reduced.

**Rational:** instead of Gaussian-like hills, some kind of rational functions
are used,

.. math::
   :label: eq-hills2

   \begin{aligned}
   V(t,{\bf s}) &=  \sum_{t_{i} < t } \Bigg[W_{i}
   \frac{1- \left(\frac{\sqrt{({\bf s}-{\bf s}^{i})^{2}}}{\Delta s^{\perp}}\right)^{n}}
   {1- \left(\frac{\sqrt{({\bf s}-{\bf s}^{i})^{2}}}{\Delta s^{\perp}}\right)^{m}}
    \nonumber \\
    &  \exp \left
   \{-\frac{\left( ({\bf s}^{i+1}-{\bf s}^{i}) \cdot ({\bf s}-{\bf s}^{i})\right)^{2}}
   {2 (\Delta s^{||}_{i})^4} \right\}\Bigg],
   \end{aligned}

where the exponents :math:`n` and :math:`m` determine the decay.

**Lorentzian:** Lorentzian functions are used in place of Gaussians.

In all the cases, a new hill is added at each step of MTD, where
:math:`\Delta t _{meta}=t_{i+1}-t_{i}` is usually chosen equal to
:math:`10\div500` steps of CP-MD (it depends on the relaxation time of the
system and the size of the hills). The center of the new hill at time
:math:`t_{i+1}` is positioned along the vector :math:`{\bf s} - {\bf s^{i}}`.

Metadynamics Keywords
---------------------

Now let’s start with the explanation of the keywords. First, the definition of
the CV is required. The selected CV are read from the input subsection enclosed
between the initial and final lines:

:math:`DEFINE` :math:`VARIABLES`

:math:`\cdots`

:math:`END` :math:`DEFINE`

Between these two lines the first input line states how many CV are used,
:math:`NCOLVAR`. In the following, each variable is described by one or more
lines, according to its type. In general, each line must start with the name of
the CV, :math:`type-name`, followed by some indexes or parameters that are
needed to specify its analytical function and the kind of atoms or species that
are involved. At the end, always on the same line of the :math:`type-name`, the
scaling factor :math:`scf` and, if the extended Lagrangian is used, :math:`k`
and :math:`M` can be given. If not specified :math:`scf`, :math:`k` and
:math:`M` take some default values.

**scf:** by default, :math:`scf=1` and it is fixed during the whole run.
Otherwise, you can write :math:`{\bf SCF}` followed by the value, or
:math:`{\bf SCA}` followed by the value, a lower bound and an upper bound. In
the latter case, the :math:`scf` is tuned along the MTD run. In practice, the
average amplitude of the CV fluctuation is checked every time to time, and, if
:math:`scf_{\alpha}\cdot\delta s_{\alpha}` is far from
:math:`\Delta s^{\perp}`, the :math:`scf_{\alpha}` is changed accordingly.

**M:** it determines how fast the :math:`s` variable spans the entire well.
Given the width of the well :math:`scf\cdot\Delta s^{\perp}` and the
temperature :math:`T_s`, it is possible to choose :math:`M` by stating the
number of complete fluctuations per ps. The default value is taken for 10
fluctuation per ps. Otherwise, you can write :math:`{\bf MCV}` followed by the
desired value for :math:`M` in
:math:`Hartree ((t)/(u.s.))^2/1822 = a.m.u. (a.u./a.s.)^2`.

**k:** it determines the dynamics of :math:`s` with respect to the dynamics of
the physical CV. If :math:`S(\cdots)` is dominated by fast modes, it is
recommended that :math:`s` be slower and its fluctuations span the entire well.
Given the characteristic frequencies of the normal modes :math:`\omega_{0}`,
:math:`k` can be chosen such that :math:`\sqrt{k/M}< \omega_{0}`. On the other
hand, we want :math:`k` big enough, so that :math:`s` and :math:`S` stay close,
and :math:`S` fluctuates many times at each position in the configuration
space. By satisfying the latter condition, the average forces due to the
underlying potential can be accurately estimated, and the trajectory lays on
the minimum energy path. Therefore, also for :math:`k`, the default value is
chosen in terms of :math:`T_s` and :math:`scf_{\alpha}\cdot\delta s_{\alpha}`
Otherwise, you can write :math:`{\bf KCV}` followed by the desired value.

On the same line, by writing **WALL+** or **WALL-**, some fixed upper and lower
boundaries for the CV can be determined. After the keyword the position of the
boundary and the value of the constant repulsive force have to be specified.

*Warning*: if even only one :math:`k` or one :math:`M` is read from input, the
Lagrangian formulation of the MTD is initialized.

The Implemented Types of CV
---------------------------

Please note, that for calculations using the Gromos QM/MM-interface (see
section :ref:`sec:qmmm`) the atom indices refer to the ordering of the
atoms as it appears in the respective GROMOS coordinate file.

-  STRETCH: Bond stretch: give the indexes of the 2 atoms :math:`i1`
   :math:`i2`, :math:`s = (d_{i1,i2})^{2}`

-  BEND: Bond angle: give the indexes of the 3 atoms defining the angle,
   :math:`i1` :math:`i2` :math:`i3`.

-  TORSION: Torsion angle: give the indexes of the 4 atoms defining the torsion
   angle, :math:`i1` :math:`i2` :math:`i3` :math:`i4`.

-  DIST: Distance between two atoms: give the indexes of the 2 atoms :math:`i1`
   :math:`i2`, :math:`s = d_{i1,i2}`.

-  DISAXIS: Distance between two atoms :math:`i1` and :math:`i2` along
   :math:`x` or :math:`y` or :math:`z` direction. :math:`i1` :math:`i2`
   :math:`n` are read next on the same line. Here :math:`n=1` means :math:`x`,
   :math:`n=2` means :math:`y` and :math:`n=3` means :math:`z` coordinate.

-  OUTP: Angle out of plane: give the indexes of the 3 atoms defining the plane
   and a fourth index of the atom for which the angle out of plane is computed,
   :math:`i1` :math:`i2` :math:`i3` :math:`i4`.

-  COORD: Coordination number (CN) of one atom with respect to all the other
   atoms in the system. The CN is defined by a Fermi-like function

   .. math:: CN_i = \sum_{j \neq i}^{NATOM}\frac{1}{1+e^{k(d_{ij}-d^0)}}

   where :math:`i` is the index of the selected atom, :math:`j` runs over all
   the other atoms in the system, :math:`k` is the parameter which determines
   the steepness of the decay and :math:`d^0` is the reference distance. After
   the type-name, in the same line, give :math:`i` :math:`k` :math:`d^0`.

-  DIFFER: Difference between two distances, give the indexes of the 3 atoms
   defining the 2 vectors, :math:`i1` :math:`i2` :math:`i3`,
   :math:`s = d_{i1i2}-d_{i2i3}`.

-  COORSP: CN of one selected atom :math:`i` with respect to only one selected
   species :math:`jsp`. The CN is defined by a Fermi like function as for
   :math:`COORD`, but in this case :math:`j` runs only over the atoms belonging
   to the selected species :math:`jsp`. After the type-name, in the same line,
   give :math:`i` :math:`jsp` :math:`k` :math:`d^0`.

-  COORGROUP: Sum of the CN of a group of atoms :math:`A` with respect to
   individual group of atoms (:math:`B`). CN is estimated using the Fermi
   function. Different cutoff distances are allowed for each type of :math:`A`
   atoms.

   .. math:: CN = \sum_{i}^{N_A}  \sum_{j}^{N_B(i)} \frac{1}{1+e^{k\left[d_{ij}-d^0(i)\right]}}

   After the keyword ``COORGROUP``, :math:`N_A` and :math:`k` should be
   specified. In the next lines should be:

   :math:`i`, :math:`d^0(i)`, :math:`N_B(i)`

   :math:`j(1) \cdots j(N_B(i))`

   This has to be done for all :math:`i` in list of :math:`A` type atoms.

-  COOR_RF: CN of one selected atom :math:`i` with respect to one selected
   species, :math:`jsp`, or a list of atoms, :math:`j1 \cdots jn_{list}`. The
   CN value is calculated as the sum of rational functions

   .. math::

      CN_i = \sum_{j \neq i}^{n_{list}} \frac{1-\left(\frac{d_{ij}}{d^0}\right)^{p}}
      {1-\left(\frac{d_{ij}}{d^0}\right)^{p+q}},

   where j runs over the indexes of the atoms belonging to :math:`jsp` or over
   the indexes given in the list :math:`j1 \cdots jn_{list}`. For the option of
   the species, you should provide, after the type-name, the indexes :math:`i`
   and :math:`jsp`, the exponents :math:`p` and :math:`q`, and the reference
   distance :math:`d^0` are read. If, instead, the list option is your choice,
   write immediately after the type-name the keyword :math:`INDAT`, and next
   the values of :math:`i`, :math:`n_{list}`, :math:`p`, :math:`q`, and
   :math:`d^0`. The indexes of the atoms belonging to the list are read from
   the next line.

   If the keyword :math:`2SHELL` is found, in the same line as
   :math:`COOR\_RF`, the first real number after this keyword is a second
   reference distance :math:`d_{2sh}`. In this case, the functional form of CN
   is modified, in order to take into account only the neighbors belonging to
   one farther shell, and :math:`d_{2sh}` is the average distance of these
   atoms from :math:`i`:

   .. math::

      CN_i^{2sh} = \sum_{j \neq i}^{n_{list}} \frac{1-\left(\frac{(d_{ij}-d_{2s})}{d^0}\right)^{p}}
      {1-\left(\frac{(d_{ij}-d_{2s})}{d^0}\right)^{p+q}}.

   For the modified CN the exponents must be even.

-  BNSWT: Reciprocal CN between 2 selected atoms, defined with the same
   functional form as the one described for :math:`COOR\_RF`. This coordinate
   states the presence of the bond between the two atoms :math:`i` and
   :math:`j`. After the type-name, give :math:`i`, :math:`j`, :math:`p`,
   :math:`q`, and :math:`d^0`.

-  TOT_COOR: Average CN of the atoms belonging to a selected species
   :math:`isp` with respect to a second selected species, :math:`jsp`, or with
   respect to a given list of atoms, :math:`j1 \cdots jn_{list}`. The same
   functional forms and input options are used, as those described for
   :math:`COOR\_RF`, but the index of one selected species :math:`isp` is read
   in place of the index of one atom.

-  DISPL1: Average displacement of one group of species with respect to a
   second group of species, computed along one specified direction in space
   (lattice axis in crystals). This CV is useful to study diffusion processes
   in condensed matter. If the keyword :math:`MIL` is found, the 3 Miller
   indexes, which define the direction in space, are read immediately after
   (default: :math:`{\bf v}=(hkl)=(100)`).

-  COORDIS:

-  PLNANG: Angle between two planes. Each plane is defined by the coordinates
   of 3 atoms; after the type-name, give the indexes of the 3 atoms defining
   the first plane, :math:`i1` :math:`i2` :math:`i3`, and the indexes of the
   atoms defining the second plane, :math:`j1` :math:`j2` :math:`j3`.

-  HBONDCH:

-  DIFCOOR: Difference between the CN of two different atoms, :math:`i1` and
   :math:`i2`, with respect to the same species :math:`jsp`, or the same list
   of atoms,\ :math:`j1 \cdots jn_{list}`. The same functional forms and input
   options are used, as those described for :math:`COOR\_RF`, but the index of
   two selected atoms are read, :math:`i1` and :math:`i2`, rather than one.

-  RMSD_AB: Given two atomic configurations :math:`A` and :math:`B`, the root
   mean square displacements (RMSD) of the actual configuration from :math:`A`,
   :math:`rmsdA`, and from :math:`B`, :math:`rmsdB`, are calculated (global
   translation and rotation are subtracted by the method of quaternions). The
   RMSD can be calculated on selected group of species: after the type name
   give the number of species (:math:`NUMSPEC`) and the indexes of the selected
   species (:math:`IS_{1} \cdots IS_{NUMSPEC}`). If :math:`NUMSPEC = 0` all the
   species are included. If in the same line the keyword :math:`{\bf FILEAB}`
   is found, next the file name is read, where the atomic positions of the
   configurations :math:`A` and :math:`B` are given. Otherwise the file name is
   by default :math:`{\bf STRUCTURE\_AB}`. File format: 2 consecutive blocks of
   :math:`1+NATOM` lines. In each block, the first line is a title (Character)
   and it is followed by the list of atomic coordinates in a.u.
   (:math:`element\_name \, x \, y \, z`).

-  COOR_CHAIN: Conditioned CN. Given three species :math:`isp1`, :math:`isp2`,
   and :math:`isp3`, the following average CN is calculated

   .. math::

      CN = \frac{1}{N_{sp1}}\sum_{i1=1}^{N_{sp1}}\left[ \sum_{i2=1}^{N_{sp2}}\left(
      \frac{1-\left(\frac{d_{i1i2}}{d_{12}^0}\right)^{p}}{1-\left(
      \frac{d_{i1i2}}{d_{12}^0}\right)^{p+q}}\times \sum_{i3=1}^{N_{sp3}}
      \frac{1-\left(\frac{d_{i2i3}}{d_{23}^0}\right)^{p}}{1-\left(
      \frac{d_{i2i3}}{d_{23}^{0}}\right)^{p+q}}\right)\right].

   After the type-name, the parameters :math:`isp1`, :math:`isp2`,
   :math:`isp3`, :math:`p`, :math:`q`, :math:`d_{12}^0`, and :math:`d_{23}^0`
   are read.

-  HYDRONIUM:

-  DIS_HYD:

-  SPIN: Distance between a selected atom and the center of the spin
   polarization :math:`(\rho_{\uparrow} - \rho_{\downarrow})`, where
   :math:`\rho_{}` indicate the polarized density. The center is located where
   the difference is maximum, and this kind of variable is useful only when
   some spin polarization is present. The position of the center in systems
   with PBC can be calculated by the definition proposed by Resta
   [resta]_ [berghold]_. Obviously, this CV can be used only together
   with :math:`LSD`. After the type-name, give the index of the selected atom.

-  VOLVAR: Volume of the cell. It can be used only with NPE/NPT MD.

-  CELLSIDE: Length of one cell’s side: give the cell-side’s index :math:`i`
   (:math:`i_a=1`,\ :math:`i_b=2`,\ :math:`i_c=3`).It can be used only with
   NPE/NPT MD.

-  CELLANGLE: Cell-angle: give the cell-angle’s index :math:`i`
   (:math:`i_{\alpha}=1`,\ :math:`i_{\beta}=2`,\ :math:`i_{\gamma}=3`)). It can
   be used only with NPE/NPT MD.

-  VCOORS This CV represents the coordination of one point (V) with respect to
   a selected species of atoms :math:`jspec` in the system:

   .. math::

      \begin{aligned}
      CN^{jspec}_{V}=\sum^{N_{jspec}}_{i\in jspec}{1\over 1 + e^{k(d_{iV}-d_0)}}
      \end{aligned}

   After the keyword the parameters :math:`jspec,k,d_0` are read. In the next
   line the coordinates of the point :math:`V` are read in a.u.

-  DIPOLE The dipole of the atoms :math:`i_1,...,i_N` with respect to the atom
   :math:`j` is defined as:

   .. math::

      \begin{aligned}
      \vec D_j={1\over Q}\sum_{i=i_1,..,i_N}q_i(\vec r_i -\vec r_j);\;\;\; Q=\sum_{i=i_1,..,i_N} q_i
      \end{aligned}

   The three spherical coordinates of :math:`\vec D_j`, that is
   :math:`(\rho_j,\theta_j,\phi_j)`, can be used independently as CV. The
   keywords are **DIPOLERHO**, **DIPOLETHA**, **DIPOLEPHI**. In the same line
   after the keywords are read the index of the atom :math:`j` and the number
   :math:`N` of atoms which constitute the dipole. In the next two lines are
   read the indexes of the :math:`N` atoms and the corresponding charge
   :math:`q_i`.

If :math:`CELL` :math:`FULL` is defined in the first line of the MTD input,
none of the CV defined above is used. The CV are the 6 cell parameters. In the
section :math:`DEFINE` :math:`VARIABLE`, the number of CV is 6 and in the
following 6 lines the scaling factors are given: for each line write the index
:math:`i` of the corresponding CV
(:math:`i_a=1`,\ :math:`i_b=2`,\ :math:`i_c=3`,\ :math:`i_{\alpha}=4`,\ :math:`i_{\beta}=5`,\ :math:`i_{\gamma}=6`)
followed by :math:`SCF` or :math:`SCA` and the desired values (see the
description at the beginning of this subsection).

Other Keywords
--------------

-  ANALYSIS: A standard MD run is performed, where the equations of motions are
   not affected by the hills-potential or the coupling potential. The selected
   CV are monitored and the values are reported in the output file, after every
   10 MD steps. This option is useful in order to observe the behavior of the
   selected CV in equilibrium conditions. With this option only two output
   files are written: *istvar_mtd*, and *enevar_mtd* (see section
   :ref:`sub:outmtd`). The former file contains the values of the
   :math:`S_{\alpha}` and their averages in time.

-  | METASTEPNUM: The maximum number of MTD steps is read from the next line,
   | :math:`I\_META\_MAX` (default: 100).

-  META_RESTART: To restart a metadynamic’s run from where the former run has
   stopprd, one can use this keyword, and write in the following line, the
   number of meta-steps completed already :math:`I\_META\_RES` (default: 0).
   Beware that for restarting the MTD in this way, the output files of the
   previous run are to be available in the run’s directory and must contain a
   number of lines at least equal to the :math:`I\_META\_RES`. From this files
   the previous history of the MTD is read and the MTD is initialized
   accordingly.

   | Otherwise, it is possible to restart from the restart file of the MTD,
     :math:`MTD\_RESTART`. This is an unformatted file, which is written
     whenever the standard CPMD restart file, :math:`RESTART.1`, is also
     written. It contains the number of meta-steps already performed, the
     number of CV used in the previous run and the information about the
     position and the size of the hills which have been already located. To
     restart the MTD from this file, the same keyword is used, and the keyword
     :math:`RFILE` is added in the same line, providing that the unformatted
     restart file is available in the run’s directory. In this second case the
     number of performed meta-step is not read from the input file but from the
     restart file
   | . Obviously, when a run is restarted, the same number and the same kinds
     of CV must be used. However masses, force constants, scaling factors, and
     the width and height of the hills can be changed.

-  MINSTEPNUM INTERMETA: The minimum number of MD steps between two MTD steps
   (in general, the MTD step is characterized by the positioning of a new hill
   in the CV-space) is read from the next line, :math:`INTER\_HILL` (default:
   100). This is a lower bound because, before the construction of a new hill,
   the displacement in the CV-space is checked, and the new step is accepted
   only if the calculated displacement is above a given tolerance.

-  MOVEMENT CHECK: The tolerance for the acceptance of a new MTD step is read
   from the next line (default: 0.001D0)

-  CHECK DELAY: The number of MD steps to be run, before a new check of the
   displacement is done, is read from the following line (default: 20).

-  MAXSTEPNUM INTERMETA: The maximum number of MD steps that can be run, before
   a new MTD step is accepted anyway, is read from the following line (default:
   300).

-  METASTORE [NO TRAJECTORY]: In the following line, three integer numbers are
   given, which indicate respectively how often (in terms of MTD steps) the
   :math:`RESTART.1` and the :math:`MTD\_RESTART` files are over-written
   (default: 50), how often the trajectory files are appended (default: 1), and
   how often a quench of the electronic wavefunctions onto the BO surface is
   performed (default: 10000). With the additional flag NO TRAJECTORY, the
   trajectories are still written according to the settings as indicated by the
   :kw:`TRAJECTORY` keyword in the &CPMD section. The selection
   of the files (e.g. turning on TRAJEC.xyz via the XYZ flag and turning
   TRAJECOTORY off via a negative value of NTRAJ in combination with the SAMPLE
   flag) is always honored.

-  MAXKINEN: From the following line, the maximum electronic kinetic energy is
   given, above which a quench of the electronic wavefunctions onto the BO
   surface is performed anyway (by default no quench is done whatever is the
   kinetic energy).

-  LAGRANGE TEMPERATURE: The temperature :math:`T_s` used to initialize the
   velocities of the :math:`{\bf s}` CV is read from the next line. By default
   :math:`T_s` is chosen equal to the temperature of the ions, the units are
   Kelvin. Notice that this keyword causes the initialization of the Lagrangian
   formulation of MTD.

-  LAGRANGE TEMPCONTROL: The control of the :math:`T_s` is activated, and the
   rescaling velocities algorithm is used. The average temperature and the
   permitted range of variation are read from the next line. By default
   :math:`T_s` is not controlled. Notice that this keyword causes the
   initialization of the extended degrees of freedom of the Lagrangian
   formulation of MTD.

-  LAGRANGE LANGEVIN: Performs Langevin dynamics for the Lagrangian formulation
   of MTD. In the next line the Temperature (Kelvin) and the friction
   :math:`\gamma` (a.u.) are read. The Langevin equation in its standard form
   writes (:math:`z` is a CV):

   .. math::

      \begin{aligned}
      M\ddot z=F(z)-\gamma M \dot z+\sqrt{2k_BT\gamma M}\eta(t)
      \end{aligned}

   where :math:`M` is the CV mass, :math:`F` a generic force field,
   :math:`\gamma` the friction coefficient, :math:`T` the temperature and
   :math:`\eta` is a white noise. The integration algorithm is a
   Velocity-Verlet which can be written as:

   .. math::

      \begin{aligned}
      &  \dot z_{n+1/2}=\dot z_n+{1\over 2}dt\Big[{F(z_n)\over M}-\gamma\dot z_n\Big]+{1\over2}\sqrt{2k_BT\gamma dt\over M}\xi_n\nonumber \\
      &  z_{n+1}=z_n+\dot z_{n+1}dt\nonumber \\
      &  \dot z_{n+1}=\dot z_{n+1/2}+{1\over 2}dt\Big[{F(z_{n+1})\over M}-\gamma\dot z_{n+1/2}\Big]+{1\over2}\sqrt{2k_BT\gamma dt\over M}\xi_n
      \end{aligned}

   the :math:`\xi_n` are independent Gaussian random numbers with mean zero and
   variance one.

-  HILLS: With this keyword it is defined the shape of :math:`V(t)`. If
   :math:`OFF` is read in the same line, no hill potential is used. The default
   hills’ shape is the Gaussian-like one described above. If :math:`SPHERE` is
   read from the same line, the second exponential term in equation
   :eq:`eq-hills` is not applied. i.e., a normal Gaussian function
   rather than a Gaussian tube formalism is used. If, instead,
   :math:`LORENZIAN` is read in the same line as :math:`HILLS`, Lorentzian
   functions are used in place of Gaussians. If, instead, :math:`RATIONAL` is
   read in the same line as :math:`HILLS`, the rational function described in
   the previous section is used; in this case, if :math:`POWER` is read, the
   exponents :math:`n` and :math:`m`, and the boosting factor :math:`f_b` are
   also read immediately after (defaults: :math:`n=2`, :math:`m=16`,
   :math:`f_b=1`). If, instead, :math:`SHIFT` is read on the same line as
   :math:`HILLS`, the shifted Gaussians are used, where the tails after a given
   cutoff are set equal to zero; in this case, if :math:`RCUT` is read, the
   cutoff :math:`rshift` and the boosting factor :math:`f_b` are also read
   immediately after (defaults: :math:`rshift = 2`,\ :math:`f_b=1`).

   In all this cases, if the symbol :math:`=` is read in the same line as
   :math:`HILLS`, the perpendicular width :math:`\Delta s^{\perp}` and the
   height :math:`W` are read immediately after (defaults:
   :math:`\Delta s^{\perp} = 0.1 \, u.s.`, :math:`W = 0.001 \, Hartree`).

-  TUNING HHEIGHT: With this keyword the height of the hills is tuned according
   to the curvature of the underlying potential. If the symbol :math:`=` is
   read in the same line as :math:`TUNING` :math:`HHEIGHT`, the lower and upper
   bounds of :math:`W` are read immediately after (defaults:
   :math:`W_{min}=0.0001 \, Hartree`, :math:`W_{max}=0.016 \, Hartree`).

-  HILLVOLUME: With this keyword the volume of the hills,
   :math:`\sim W\cdot (\Delta s^{\perp})^{NCOLVAR}`, is kept constant during
   the MTD run, i.e. when the height changes due to the tuning (see the
   previous keyword), the width is changed accordingly. This option is can be
   used only if the tuning of the hills’ height is active.

-  CVSPACE BOUNDARIES: With this keyword the confinement of the CV-space is
   required, in the direction of a selected group of CV. The number of
   dimensions of the CV-space, in which the confinement is applied, is read
   from the next line, :math:`NUMB` (default: 0). The following :math:`NUMB`
   lines describe the type of confinement. For each line the following
   parameters are required: index of the CV (as from the list given in the
   definition of the CV), the strength of the confinement potential
   :math:`V_{0}` in Hartree, two real numbers, :math:`C1` and :math:`C2`, that
   determine from which value of the CV the confining potential is active.
   Finally, if the keyword :math:`EPS` is read in the same line, the real
   number, which is read immediately after, determines how smoothly the
   confining potential is switched on (default :math:`\varepsilon = 0.01`). The
   confining potential can be used with the following CV:

   DIST: :math:`V_{conf} = V_{0}\left(\frac{s_{\alpha}}{C1}\right)^{4}` and it
   becomes active only if :math:`s_{\alpha} > C2`.

   DIFFER:
   :math:`V_{conf} = V_{0}\left(\frac{\vert s_{\alpha}\vert}{C1}\right)^{4}`
   and it becomes active only if :math:`s_{\alpha} > C2` or
   :math:`s_{\alpha} < -C2` .

   Coordination numbers: if :math:`(CN-\varepsilon) < C1`, :math:`V_{conf} =
   V_{0} (\frac{C1}{CN})^{10}`, if :math:`(CN+\varepsilon) > C2`,
   :math:`V_{conf} = V_{0} (\frac{CN}{C2})^{10}`.

-  RESTRAIN VOLUME: With this keyword a confining potential is applied to the
   volume variations. The option can be used only in combination with the
   NPE/NPT MD. From the next line, the following parameters are required:
   :math:`f_{min}`, :math:`f_{max}`, and :math:`V_0`. The factors
   :math:`f_{min}` and :math:`f_{max}` multiplied by the initial volume of the
   cell, give, respectively, the lower and upper bounds for the volume, whereas
   :math:`V_0` gives the strength of the confining potential.

-  MULTI NUM: his keyword should be used when a multiple set of MTD runs are
   performed simultaneously on the same system. Here the number of separated
   sets CV for each subsystem has to be given in the following line,
   :math:`NCVSYS_1 \cdots NCVSTS_{NSUBSYS}`. This means that the first
   :math:`NCVSYS_1` CV, in the list of those defined in the :math:`DEFINE`
   :math:`VARIABLES` section, will belong to the first subset, the next
   :math:`NCVSYS_2` to the second, and so on. This option is implemented only
   together with the extended Lagrangian formulation.

-  MONITOR: This keyword requires that an additional monitoring of the values
   of the CV is performed along the MD trajectory. This means that the values
   are written on an output file every *WCV_FREQ* MD steps, even if no hill is
   added at that step. The frequency for updating the file is read from the
   following line, the file name is *cvmdck_mtd*, and it is not created if this
   option is not activated.

-  RANDWALK: In the case of multiple walker metadynamics, collective variables
   of all walkers are initialized with different random velocities.

.. _sub:outmtd:

Output files
------------

During a run of MTD, several output file are updated at each MTD step, which
are characterized by the extension :math:`\_mtd`. These files contain the
history of the MTD, the parameters giving the additional potential, and other
information that can be useful during the analysis of results. The first column
for all these files is the CPMD step number (NFI) that corresponds to this MTD
step. In the case of :math:`MULTI` MTD some of the files have a further
extension :math:`\_s\#`, which indicates the related subsystem.

Extended Lagrangian MTD:

*istvar_mtd*: the first :math:`NCOLVAR` columns are the
:math:`S_{\alpha}(\cdots)`, the next :math:`NCOLVAR` columns are differences
:math:`S_{\alpha}(\cdots) - s_{\alpha}`.

*colvar_mtd*: the first :math:`NCOLVAR` columns are the :math:`s_{\alpha}`, the
next :math:`NCOLVAR` columns are the corresponding scaling factors
:math:`scf_{\alpha}`.

*parvar_mtd*: norm of the total displacement in the CV-space,
:math:`\Delta s^{||}`, hill’s width, :math:`\Delta s^{\perp}`, hill’s height,
:math:`W`.

*disvar_mtd*: the first :math:`NCOLVAR` columns are the displacements of the
:math:`s_{\alpha}`, the next :math:`NCOLVAR` columns are their diffusivities,
the final :math:`NCOLVAR` columns are the coupling constants, :math:`k_\alpha`.

*velvar_mtd*: the velocities of the :math:`s_{\alpha}`.

*forfac_mtd*: the first :math:`NCOLVAR` columns are the forces coming from the
coupling potential (sum of harmonic terms, :math:`V_{harm}`), the next
:math:`NCOLVAR` columns are the forces coming from :math:`V(t)`, the last
:math:`NCOLVAR` are the forces coming from the confining potential.

*enevar_mtd*: ions temperature, electrons kinetic energy, :math:`{\bf s}` CV
temperature :math:`2 K_{\bf s}/(NCOLVAR k_{B})`, :math:`V_{harm}`, :math:`V(t)`
at the actual position in CV-space, KS energy :math:`E_{KS}`,
:math:`E_{tot}+ K_{s} + V_{harm}`, :math:`E_{tot}+ K_{s} + V_{harm} + V(t)`.
*cvmdck_mtd*: monitoring o the CV along the MD trajectory. This file is updated
every it WCV_FREQ MD steps (see previous section *MONITOR*)

Direct MTD:

*colvar_mtd*: the first :math:`NCOLVAR` columns are the
:math:`s(\cdots)_{\alpha}`, the next :math:`NCOLVAR` columns are the
corresponding scaling factors :math:`scf_{\alpha}`.

*sclvar_mtd*: the first :math:`NCOLVAR` columns are the scaled
:math:`s(\cdots)_{\alpha}`, the next :math:`NCOLVAR` columns are the
corresponding scaled diffusivities.

*parvar_mtd*: norm of the total displacement in the CV-space,
:math:`\Delta s^{||}`, hill’s width, :math:`\Delta s^{\perp}`, hill’s height,
:math:`W`.

*disvar_mtd*: the first :math:`NCOLVAR` columns are the displacements of the
:math:`s_{\alpha}`, the next :math:`NCOLVAR` columns are their diffusivities,
the final :math:`NCOLVAR` columns are the coupling constants, :math:`k_\alpha`.

*forfac_mtd*: the first :math:`NCOLVAR` columns are the forces coming from
:math:`V(t)`, the last :math:`NCOLVAR` are the forces coming from the confining
potential.

*enevar_mtd*: ion temperature, electrons kinetic energy, :math:`V(t)` at the
actual position in CV-space, KS energy :math:`E_{KS}`, :math:`E_{tot}+ V(t)`.

.. _sect:mw:

Using multiple walker metadynamics
----------------------------------

Multiple walker metadynamics is activated using the :math:`MULTIPLE`
:math:`WALKER` keyword in the initial input line of metadynamics (i. e., after
:math:`METADYNAMICS` keyword). Multiple walker using the extended Lagrangian
metadynamics in combination with the Car–Parrinello type of dynamics is only
implemented at the moment. From the same line of :math:`MULTIPLE`
:math:`WALKER` keyword, the number of walkers (NWALK) is read as
:math:`NW=NWALK` (without any space in between).

:math:`NWALK` replicas are created, and for each DFT forces and energy
calculations are done independently. However, all the replicas fill the same
free energy surface encompassed by the set of reaction coordinates
specified [Raiteri05]_. Implementation is done in such a way that
each replica belongs to a different processor group,and each processor group is
able to perform independent DFT calculations in parallel; if :math:`NPROC`
number of processors are used, each replica is using :math:`NPROC/NWALK` number
of processors [Nair-inside-08]_ for computations. See the output
file for the details on the division of processors in to corresponding
processor groups. Note that output of all the walkers are currently dumped in
to the (same) standard output file. Trajectory, geometry and energy files of
each walkers are separately written out in files having their usual names
augmented with \_\ :math:`IWALK`, where :math:`IWALK` is the walker ID.

If a multiple walker run has to be started from a single restart file, copy or
link it :math:`NWALK` times as RESTART_1, :math:`\cdots`
RESTART\_\ :math:`NWALK` (similarly the MTD_RESTART file, if is also
restarted). In the procedure of creating new walkers, like above from one
restart or increasing walker numbers during the run, it is advised to initially
run with zero hill height, still keeping all the biasing potentials accumulated
upto then (i. e. restarting from the MTD_RESTART file if any previous
metadynamics runs have been done), and use :math:`RANDWALK` keyword until the
(all) walkers are far apart (at least 1.5 x hill width) from each other. In the
following (restart) run, use the required hill height and remove the
:math:`RANDWALK` keyword. Note that the frequency of adding hills will be
nearly :math:`NMTD/NWALK` if :math:`NMTD` is the number of MD steps required to
add a hill using 1 walker. Thus consider decreasing the :math:`MINSTEPNUM`
:math:`INTERMETA` appropriately. However, it is highly recommended to use the
adaptive metadynamics time step using :math:`MOVEMENT` :math:`CHECK` keyword,
and the tolerance is typically 1.5 times the hill width
parameter [Nair-jacs-08]_; better set the :math:`CHECK`
:math:`DELAY` to 1. Displacement tolerance are also forced to satisfy between
the walkers when :math:`MOVEMENT` :math:`CHECK` is used. Note that it is
possible to decrease the number of walkers, in a straight forward manner during
a restart run.

Shooting from a Saddle
----------------------

Once one reactive trajectory has been found, one may want to determine more
precisely the position of the transition state region. A standard way to do
this is to select some points along the trajectory, and, by shooting with
random velocities a new MD from this point, measure the probability to reach
the surrounding basins of attraction [chandler]_. The different
basins of attraction can be identified by different values of a selected set of
CV. One can say that the trajectory has fallen in one of the known basins when
all the actual CV values satisfy the values characterizing that basin, within a
certain tolerance. Given a set of coordinates, one can start a CPMD run where
this check is iterated as many times as you like, in order to establish the
commitor distribution.The search for the saddle point region is initialized
when in the section &ATOMS &END of the input file, the keyword :math:`SADDLE`
:math:`POINT` is found. In what follows, a subsection for the description of
the selected CV is required. It has the same format as the one used for the MTD
run.

Keywords
--------

The list of keyword regarding the shooting needs to be ended by the line

:math:`END` :math:`SADDLE`

Other keywords are:

-  KNOWN_MINIMA The values of the CV, which characterize the known basin of
   attraction, are read from the following lines. The first line after the
   keyword contains the number of the known minima :math:`NCVMIN`. The next
   :math:`NCVMIN` lines contain the set of values for each of these minima. The
   list of the values must keep the same order used in the definition of the
   CV. If on the same line as :math:`KNOWN\_MINIMA`, the keyword
   :math:`EXTENDED` is also found, each line contains :math:`NCOLVAR` more
   entries, which are the tolerances for the acceptance of the corresponding
   minimum configuration (the order of the tolerances must be the same as the
   one for the CV values). By using the :math:`EXTENDED` keyword, each minimum
   configuration can be accepted with different tolerances.

-  SADDLE TOLERANCES If the :math:`EXTENDED` keyword is not used, one single
   set of tolerances (one for each CV) can be given by using this keyword. The
   tolerances are read from the next line in the same order used for the CV
   definition. Otherwise, default values are assigned at each tolerance.

-  MAXSEARCH The maximum number of trials, where a new MD trajectory is
   generated, is read from the next line. At each trial, the MD starts from the
   same initial coordinates, whereas the initial velocities are randomly
   generated at every new restart. During one trials, every :math:`NSTEP` the
   actual values of the CV are checked and compared to the values given for the
   known minima. If all the values of one of these minima are satisfied within
   the given tolerances, the MD is stopped and restarted for the next trial.

-  STEPCHECK The number of MD steps between two consecutive checks is read from
   the next.

-  MAXCHECKS The maximum number of checks for each trial is read from the next
   line.
