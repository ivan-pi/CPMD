Perturbation Theory / Linear Response
=====================================

General
-------

Ref: [apdsmp]_.

Perturbation theory describes the reaction of a system onto an external
perturbation. At the time when the perturbation occurs, the system is in its
ground state (or unperturbed state). The perturbation then changes slightly the
potential energy surface and therefore also the state where the system’s energy
is minimum. As a consequence, the system tries to move towards that state of
minimum energy. This movement or the new state often have properties which can
be accessed experimentally. Example: An external electric field will slightly
deform the electronic cloud, creating a dipole. That dipole can then be
measured.

Assume that the magnitude of the perturbation is small compared to the strength
of the forces acting in the unperturbed system. Then, the change in the minimum
energy state will be small as well and perturbation theory can be applied to
compute how the system reacts onto the perturbation. Generally, the Schrödinger
equation is expanded in powers of the perturbation parameter (ex: the strength
of the electric field), and the equations obtained for those powers are solved
individually. At power zero, one finds the equation of the unperturbed system:

.. math::

   \begin{aligned}
     (H^{(0)}- \varepsilon_k) \vert\varphi_{k}^{(0)}\rangle &=  0.
   \end{aligned}

For the part which is linear in the perturbation, the general format of the
resulting equation is

.. math::
   :label: eqn-sternheimer

   \begin{aligned}
     (H^{(0)}- \varepsilon_k) \vert\varphi_{k}^{(1)}\rangle  &=  - H^{(1)}\vert\varphi_{k}^{(0)}\rangle.
   \end{aligned}

Grosso modo, this equation is solved during a linear response calculation
through a wavefunction optimization process for
:math:`\vert\varphi_{k}^{(1)}\rangle`.

The presence of a first order perturbation correction for the wavefunctions,
:math:`\vert\varphi_{k}^{(\mathsf{tot})}\rangle = \vert\varphi_{k}^{(0)}\rangle + \vert\varphi_{k}^{(1)}\rangle`
implies that the total density of the perturbed system is no longer equal to
the unperturbed one, :math:`n^{(0)}`, but also contains a first order
perturbation correction, :math:`n^{(1)}`. That density is given by

.. math::

   \begin{aligned}
   n^{(1)}(\mathbf{r}) &=  \sum_k
              \langle \varphi^{(1)}_k \vert \mathbf{r}\rangle \;
              \langle \mathbf{r}\vert \varphi^{(0)}_k \rangle
              + \text{c.c.}
   \end{aligned}

The Hamiltonian depends on the electronic density. Therefore, the first order
density correction implies automatically an additional indirect perturbation
Hamiltonian coming from the expansion of the unperturbed Hamiltonian in the
density. It has to be added to the explicit perturbation Hamiltonian determined
by the type of the (external) perturbation. The contribution is

.. math::
   :label: eqn-implicit-hamiltonian

   \begin{aligned}
   H^{(1)}_{\mathsf{indirect}}(\mathbf{r}) &=  \int d^3r^\prime\;
      \frac{\partial H^{(0)}(\mathbf{r})}{\partial n^{(0)}(\mathbf{r}^\prime)}\;
      n^{(1)}(\mathbf{r}^\prime)
   \end{aligned}

The calculation of this indirect Hamiltonian represents almost 50% of the
computational cost of the response calculation, especially in connection with
xc-functionals. After several unsuccessful trials with analytic expressions for
the derivative of the xc-potential with respect to the density, this is done
numerically. That means that at each step, the xc-potential is calculated for
the density :math:`n^{(0)}+\epsilon n^{(1)}` and for
:math:`n^{(0)}-\epsilon n^{(1)}` (with an :math:`\epsilon` empirically set to
0.005), and the derivative needed in (:eq:`eqn-implicit-hamiltonian`)
is calculated as

.. math::

   \begin{aligned}
     \int d^3r^\prime\; n^{(1)}(\mathbf{r}^\prime)\;
     \frac{\partial v_{\text{xc}}}{\partial n^{(0)}(\mathbf{r}^\prime)}
     &=  \frac{v_{\text{xc}}[n^{(0)}+\epsilon n^{(1)}]
         - v_{\text{xc}}[n^{(0)}-\epsilon n^{(1)}]}
         {2\epsilon}.
   \end{aligned}

In the case of the local density approximation, the derivative can be done
analytically, in which case it only needs to be done once. This improves the
performance of the optimization.

.. _sec:resp-section:

&RESP section input
-------------------

Generally, the keyword :kw:`LINEAR RESPONSE` in the &CPMD
input section initiates the calculation. In the section &RESP, the type of the
perturbation needs to be specified. Either one of the following keywords must
appear:

   :kw:`PHONON` :kw:`LANCZOS` :kw:`RAMAN`
   :kw:`FUKUI` :kw:`KPERT` :kw:`NMR` :kw:`EPR`
   :kw:`HARDNESS` :kw:`EIGENSYSTEM`
   :kw:`INTERACTION` :kw:`OACP`

The first six types are discussed in detail in the following. An overview is
also contained in the file respin_p_utils.mod.F90. In addition to the specific
keywords of every option, there are several keywords which are common to all
perturbation types. They determine fine-tuning parameters of the wavefunction
optimization process, and usually you do not need to change them. A # indicates
that a command takes an argument which is read from the next line. All other
keywords toggle between two states and do not require any argument. Those
keywords can be put together, and they can also figure in the main keyword line
(example: **NMR NOOPT FAST WANNIERCENTERS**)

NB: The linear response code works with all cell symmetries [1]_, but it is
**not implemented** for :math:`k`-points.

-  **CG-ANALYTIC:** The wavefunction optimization uses a preconditioned
   conjugate gradient technique. The optimum length of the “time step” can be
   calculated analytically assuming a purely linear equation, according to the
   Numerical Recipes Eq. 10.6.4. However, this is somewhat expensive, and
   experience shows that the time step is almost constant except at the very
   beginning. Therefore, it is only calculated a few times, and later on, the
   last calculated value is used. This option controls the number of times the
   step length is calculated analytically. Default is 3 for NMR and 99 for all
   other perturbations.

-  **CG-FACTOR:** The analytic formula for the time step assumes that the
   equation to be solved is purely linear. However, this is not the case, since
   the right hand side can still depend on the first order wavefunctions
   through the dependence of the perturbation Hamiltonian :math:`H^{(1)}` on
   the perturbation density :math:`n^{(1)}`. Therefore, the analytic formula
   has a tendency to overshoot. This is corrected by an empirical prefactor
   which is controlled by this option. Default is 0.7.

-  **CONVERGENCE:** The criterion which determines when convergence is reached
   is that the maximum element of the gradient of the energy with respect to
   the wavefunction coefficients be below a certain threshold. This value is
   read from the next line. Default is 0.00001. Experience shows that often, it
   is more crucial to use a strict convergence criterion on the ground state
   wavefunctions than for the response. A rule of thumb is that good results
   are obtained with a 10 times stricter convergence on the ground state
   orbitals compared to that of the response orbitals.

-  **HTHRS or HAMILTONIAN CUTOFF:** The preconditioning calculates the diagonal
   :math:`(\mathbf{G},\mathbf{G})` matrix elements of
   :math:`\eta=H^{(0)}-\frac 1N \sum_k\varepsilon_k` to do an approximate
   inversion of Eq. (:eq:`eqn-sternheimer`). However, these diagonal
   values can become very small, yielding numerical instabilities. Therefore, a
   smoothing is applied instead of simply taking the reciprocal values:

   .. math::
      :label: eqn-precond-1

      \eta^{-1} \mapsto \left(\eta^2 + \delta^2\right)^{-1/2}

   .. math::
      :label: eqn-precond-2

      \eta^{-1} \mapsto \frac {\eta} {\eta^2 + \delta^2}

   The value of the parameter :math:`\delta` in a.u. is read from the line
   after **HTHRS**, default is 0.5. By default,
   Eq. (:eq:`eqn-precond-1`) is used. **TIGHTPREC** switches to
   Eq. (:eq:`eqn-precond-2`).

-  **NOOPT:** In order for the wavefunction optimization to work properly, the
   ground state wavefunction must be converged. For this reason, a ground state
   optimization is performed by default prior to computing the response. When
   restarting from an already converged wavefunction, this step can be skipped
   through this keyword and the computer time for initializing the ground state
   optimization routine is saved. However, the use of this option is strongly
   discouraged.

-  **POLAK:** There are several variants of the conjugate gradient algorithm.
   This keyword switches to the Polak-Ribiere formulation (see the Numerical
   Recipes, Eq. 10.6.7) which is usually significantly slower but safer in the
   convergence. By default, the Fletcher-Reeves formula is used.

-  **TIGHTPREC:** Switches to another preconditioning formula. See **HTHRS**.

Response output
---------------

While the calculations are being done, the program prints the progress of the
optimization process:

-  The “scaling input” prints the number by which the right-hand-side of
   Eq. (:eq:`eqn-sternheimer`) is multiplied in order to deal with
   reasonable numbers during the optimization. The number is determined through
   the condition

   .. math::

      \begin{aligned}
        \vert\vert\;H^{(1)}\vert\varphi_{k}^{(0)}\rangle\;\vert\vert_2 &=  1.
      \end{aligned}

   When leaving the optimization routine, the inverse scaling is applied to the
   :math:`\vert\varphi_{k}^{(1)}\rangle`, of course.

-  The standard output shows (A) the maximum gradient of the second order
   energy with respect to the first order perturbation wavefunction, (B) the
   norm of that gradient, (C) the second order energy, (D) its difference with
   respect to the last value, and finally (E) the CPU time needed for one step.
   The last value decreases by somewhat after the number of steps determined by
   the **CG-ANALYTIC** keyword, because the analytic line search is no longer
   performed, as discussed above.

-  A line full of tildes (:math:`\tilde{\ }`) indicates that the energy has
   increased. In that case, the conjugate direction is erased and the conjugate
   gradient routine is restarted. Also the time step is calculated again using
   the analytic quadratic line search formula.

-  The warning “line search instable” indicates that the length which has been
   calculated using the analytic quadratic line search approximation (Numerical
   Recipes Eq. 10.6.4) has given a numerical values larger than 3. This does
   not happen in normal cases, and it can yield to program crashes due to
   floating point overflows (:math:`\mapsto` ``NaN``, not a number). Thus, a
   safe value, 1, is used instead.

-  The warning “gradient norm increased by more than 200%” indicates that the
   quadratic approximation is not valid in the momentary position of the
   optimization. If it was, the gradient (of the energy with respect to the
   first order perturbation wavefunctions) would only decrease and finally
   reach zero. If this situation occurs, the conjugate direction is erased and
   the conjugate gradient algorithm is restarted.

Phonons
-------

**Theory**

A phonon corresponds to small displacements of the ionic positions with respect
to their equilibrium positions. The electrons principally follow them, in order
to minimize again the energy of the system.

The expansion of the Hamiltonian in powers of the displacement
:math:`u^R_\alpha` of the ion (labeled by its position :math:`R`) in the
Cartesian direction :math:`\alpha=1,2,3` consists of two parts [2]_:

.. math::

   \begin{aligned}
   H^{(1)}&=  H^{(1)}_C + H^{(1)}_{PP}
   \end{aligned}

The contribution :math:`H^{(1)}_C` comes from the Coulomb term, the
electrostatic potential:

.. math::

   \begin{aligned}
   H^{(1)}_C &=  u^R_\alpha
       \frac \partial {\partial R_\alpha} \; \frac{Z_R}{|\mathbf{r}- \mathbf{R}|}.
   \end{aligned}

The second is due to the pseudopotential which is rigidly attached to the ions
and which must be moved simultaneously. In particular, the nonlocal
pseudopotential projectors must be taken into account as well:

.. math::
   :label: eq-phon

   \begin{aligned}
   H^{(1)}_{PP} &=  u^R_\alpha
       \frac \partial {\partial R_\alpha} \;
       \left[\sum_i
          \vert\text{\tt P}^R_i\rangle \langle\text{\tt P}^R_i\vert
       \right]\\
      &=  u^R_\alpha\sum_i\;\left[
       \left[\frac\partial{\partial R_\alpha}\;
             \vert\text{\tt P}^R_i\rangle\right]
             \langle\text{\tt P}^R_i\vert
       +
       \vert\text{\tt P}^R_i\rangle
             \left[\frac\partial{\partial R_\alpha}\;
             \langle\text{\tt P}^R_i\vert\right]
       \right]
   \end{aligned}

where :math:`\vert\text{\tt P}_i^R\rangle` designates the projectors, whatever
type they are. The index :math:`i` comprises the :math:`l` and :math:`m`
quantum numbers, for example. The superscript :math:`R` just says that of
course only the projectors of the pseudopotential of the displaced atom at
:math:`R` are considered in this equation.

In CPMD 4.3, these projectors are stored in G-space, and only one copy is
stored (that is the one for a fictitious ion at the origin, :math:`R=0`). The
projectors for an ion at its true coordinates is then obtained as

.. math::
   :label: eqn-translation-formula

   \begin{aligned}
     \langle \mathbf{G}\vert\text{\tt P}_i^R\rangle &=
     \text{e}^{i\mathbf{G}\cdot\mathbf{R}}\;\langle \mathbf{G}\vert\text{\tt P}_i^{R=0}\rangle.
   \end{aligned}

This makes the derivative :math:`\frac{\partial}{\partial R_\alpha}`
particularly simple, as only the :math:`i\mathbf{G}` comes down, and the
translation formula (:eq:`eqn-translation-formula`) remains valid for
:math:`\vert\text{\tt P}_i^R\rangle`. Thus, there is only an additional
nonlocal term appearing which can be treated almost in the same way as the
unperturbed pseudopotential projectors.

A perturbative displacement in Cartesian coordinates can have components of
trivial eigenmodes, that is translations and rotations. They can be written *a
priori* (in mass weighted coordinates in this case) and thus projected out from
the Hessian matrix;

.. math::

   \begin{aligned}
   \mathbf{t}_j=\left(\begin{array}{c}
               \sqrt{m_1}(\mathbf{e}_j) \\
               \sqrt{m_2}(\mathbf{e}_j) \\
               \vdots          \\
               \sqrt{m_N}(\mathbf{e}_j) \\
               \end{array} \right) \qquad
   \mathbf{s}_j=\left(\begin{array}{c}
               \sqrt{m_1}(\mathbf{e}_j\times\mathbf{R}_1) \\
               \sqrt{m_2}(\mathbf{e}_j\times\mathbf{R}_2) \\
               \vdots          \\
               \sqrt{m_N}(\mathbf{e}_j\times\mathbf{R}_N) \\
               \end{array} \right)
   \end{aligned}

where j=x,y,z, :math:`\mathbf{e}_j` = unit vectors in the Cartesian directions
and :math:`\mathbf{R}_k` = positions of the atoms (referred to the COM). The
rotations :math:`\mathbf{s}_j` constructed in this way are not orthogonal;
before being used they are orthogonalized with respect to the translations and
with each other.

The projection on the internal mode subspace is done this way: first one
constructs a projector of the form,

.. math::

   \begin{aligned}
   \mathbf{P}=\sum_j(\mathbf{t}_j\cdot\mathbf{t}_j^T+\mathbf{s}_j\cdot\mathbf{s}_j^T)
   \end{aligned}

and then applies the projector to the Hessian matrix,

.. math::

   \begin{aligned}
   \mathbf{A}=(\mathbf{I}-\mathbf{P})\cdot\mathbf{A}\cdot(\mathbf{I}-\mathbf{P})
   \end{aligned}

with  being the unit matrix. The projection is controlled by the keyword
:kw:`DISCARD`, *vide infra*.

**Phonon input**

The input for the phonon section is particularly simple due to the absence of
any special keywords. Only the word :kw:`PHONON` should appear in the
&RESP section.

**Phonon output**

In total analogy to CPMD 4.3’s :kw:`VIBRATIONAL ANALYSIS`, the displacement of all atoms in all
Cartesian directions are performed successively. The difference is, of course,
that there is no real displacement but a differential one, calculated in
perturbation theory. Thus, only one run is necessary per ion/direction [3]_.

At the end, the harmonic frequencies are printed like in the :kw:`VIBRATIONAL ANALYSIS`. They should coincide to a few percent.

Lanczos
-------

**Theory**

Ref: [ffmp]_

A different way of diagonalizing the Hessian matrix comes from the Lanczos
procedure. It is easy to generalize eq.\ :eq:`eq-phon` to a collective
displacement of atoms,

.. math::

   \begin{aligned}
   H^{(1)}_{PP} &=  \sum_{R,\alpha}u^R_\alpha
       \frac \partial {\partial R_\alpha} \;
       \left[\sum_i
          \vert\text{\tt P}^R_i\rangle \langle\text{\tt P}^R_i\vert
       \right]\\
      &=   \sum_{R,\alpha} u^R_\alpha\sum_i\;\left[
       \left[\frac\partial{\partial R_\alpha}\;
             \vert\text{\tt P}^R_i\rangle\right]
             \langle\text{\tt P}^R_i\vert
       +
       \vert\text{\tt P}^R_i\rangle
             \left[\frac\partial{\partial R_\alpha}\;
             \langle\text{\tt P}^R_i\vert\right]
       \right]
   \end{aligned}

In this way the information contained in the Hessian matrix, **A**, can be used
to compute terms of the type

.. math::
   :label: eq-lr-12

   \begin{aligned}
   \mathbf{A}\cdot\mathbf{w}=
   \left( \begin{array}{c}
          \frac{\partial^2{{E}}}{\partial{{R_{1x}}}\partial{{\mathbf{w}}}}\\
          \frac{\partial^2{{E}}}{\partial{{R_{1y}}}\partial{{\mathbf{w}}}}\\
          \vdots                 \\
          \frac{\partial^2{{E}}}{\partial{{R_{Nz}}}\partial{{\mathbf{w}}}}
          \end{array}\right)
   \end{aligned}

where **w** is the collective displacement. This is the building block for
Lanczos diagonalization, that is performed by iterative application of the
symmetric matrix to be diagonalized, over a set of vectors, known as Lanczos
vectors, according to the scheme

-  :math:`\mathbf{r}_0=\mathbf{q}_1; \beta_0=1; \mathbf{q}_0=0; k=0`

-  :math:`\beta_k\neq 0`

-  :math:`k=k+1`

-  :math:`\alpha_k=\mathbf{q}_k^T\mathbf{A}\mathbf{q}_k`

-  :math:`\mathbf{r}_k=\mathbf{A}\mathbf{q}_k-\alpha_k\mathbf{q}_k-\beta_{k-1}\mathbf{q}_{k-1}`

-  :math:`\beta_k=\parallel\mathbf{r}_k\parallel_2`

-  :math:`\mathbf{q}_{k+1}=\mathbf{r}_k/\beta_k`

-  Orthogonalization,
   :math:`\mathbf{q}_{k+1}\perp\left\{\mathbf{q}_1,\cdots,\mathbf{q}_k\right\}`

-  Diagonalization of :math:`\mathbf{T}_k`

-

The matrix  is thus projected onto a k\ :math:`\times`\ k subspace, in a
tridiagonal form, :math:`\mathbf{T}_k`. ’s eigenvalues are approximations to
’s, while the eigenvectors are brought in the n\ :math:`\times`\ n space by
means of the orthogonal matrix
:math:`\mathbf{Q}=[\mathbf{q}_1,\mathbf{q}_2,\dots,\mathbf{q}_k]`. The
advantage with respect to the usual diagonalization schemes is that the highest
and the lowest end of the spectrum tend to converge before the ionic degrees of
freedom are fully explored.

The projection of the trivial eigenmodes is performed by eliminating their
contribution from every Lanczos vector at the orthogonalization step reported
in the algorithm above,

.. math::

   \begin{aligned}
   \mathbf{q}_k=\mathbf{q}_k-\sum_j\left((\mathbf{q}_k^T\cdot\mathbf{t}_j)\mathbf{t}_j +
   (\mathbf{q}_k^T\cdot\mathbf{s}_j)\mathbf{s}_j\right)
   \end{aligned}

This procedure seems slightly different from that of the :kw:`PHONON`
case, but they are perfectly equivalent. It is controlled by the same keyword
:kw:`DISCARD` with the same arguments.

**Lanczos input**

The input for the Lanczos routine is given by the keyword
:kw:`LANCZOS` plus some optional arguments and a mandatory following
line with numerical data. Please note that this keyword is analogous to that
for the Lanczos diagonalization of the electronic degrees of freedom. Given its
presence in the &RESP section there should be no ambiguity.

-  **LANCZOS [CONTINUE,DETAILS];** the keyword :kw:`LANCZOS` simply
   activates the diagonalization procedure. At every cycle the program writes a
   file, ``LANCZOS_CONTINUE.1`` which contains the information about the
   dimension of the calculation, the number of iteration, the elements of the
   Lanczos vectors and the diagonal and subdiagonal elements of the  matrix.

   -  With the option **CONTINUE** one restart a previous job, asking the
      program to read the needed information from the file
      ``LANCZOS_CONTINUE``; if the program does not find this file, it stops.
      This option is to be used when the calculation that one wants to perform
      does not fit in the time of a queue slot.

   -  The argument **DETAILS** prints a lot of information about the procedure
      at the end of the output. It is intended for debugging purposes.

   The subsequent line is mandatory, and contains three numerical parameters;

   -  ``Lanczos_dimension``; it is the dimension of the vibrational degrees of
      freedom to be explored. Normally it is
      3\ :math:`\times`\ N\ :math:`_{\rm at}`, -3 if you are eliminating the
      translations or -6 if you are eliminating also the rotations(*vide
      infra*). It is possible to use lower values. Higher are non sense.

   -  ``no. of iterations``; it is the number of cycles that are to be
      performed during the actual calculation. It must be :math:`\leq` to
      ``Lanczos_dimension``; the program checks and stops if it is higher. In
      case of a continued job, the program checks if the given number of
      iterations + the index of the iteration from which it restarts is within
      the limit of ``Lanczos_dimension``; if not it resets the number to fit
      the dimension, printing a warning on std. output.

   -  ``conv_threshold``; it is the threshold under which an eigenvector is to
      be considered as converged. It is the component of the eigenvector over
      the last Lanczos vector. The lower it is, the lower is the information
      about this mode which has still to be explored.

-  **DISCARD {PARTIAL,TOTAL,OFF,LINEAR};** this keyword controls the
   elimination of the trivial eigenmodes (*i.e.* translations and rotations)
   from the eigenvector calculations. It works both for :kw:`PHONON`
   and :kw:`LANCZOS`. Omitting it is equivalent to **DISCARD
   PARTIAL**. When using it, the choice of one of the arguments is mandatory;
   the program stops otherwise. They are;

   -  **PARTIAL**; only the translational degrees of freedom are eliminated
      (useful for crystals). This is the default.

   -  **TOTAL**; both translational and rotational degrees of freedom are
      eliminated.

   -  **OFF**; the option is disactivated and no projection is performed.

**Lanczos output** In the output, all the information about the perturbed
wavefunction optimization are reported, just like in the :kw:`PHONON`
case. The differences are in the form of the perturbation and in the eigenmodes
information coming from the diagonalization of :math:`\mathbf{T}_k`.

The report of the perturbation for a water dimer reads like;

::

    **********************   perturbations    **********************

    **** atom=    1     O    .04171821      .07086763      .07833475
    **** atom=    2     O    .03615521      .08499767      .07082773
    **** atom=    3     H    .29222111      .13357862      .09032954
    **** atom=    4     H    .33764012      .15750912      .25491923
    **** atom=    5     H    .06426954      .26020430      .01822161
    **** atom=    6     H    .15765937      .27370013      .29183999
    cpu time for wavefunction initialization:          89.22 seconds

where at every atom corresponds the x,y,z displacements applied to calculate
the perturbed wavefunctions.

At every iteration information about the elements of the :math:`\mathbf{T}_k`
matrix, alpha’s and beta’s, are reported. Here we are at the end of the
calculation. Note the numerical zero in the final +1 beta value. Then the
spectrum in a.u. is reported, together with the convergence information. At the
last iteration there are vectors which are not “converged”. But this comes only
from the definition, since some of the eigenvectors *must* have a component
over the last Lanczos vectors. Following there are the familiar eigenvalues in
cm\ :math:`^{-1}`.

::

    ****************************************************************
     *+* L2 norm[n[i+1]]       =  0.100262974102936016E-02
     *=*   overlap: alpha[ 12 ] =  0.982356847531824437E-03
     *=*   off-diag: beta[ 13 ] =  0.300011777832112837E-19
     *=*   norm:               =  0.300011777832112837E-19
    ****************************************************************
    *** SPECTRUM, run 12 :
    *** eigenvalue  12 =     .4855525       (converged: .000000).
    *** eigenvalue  11 =     .4785309       (converged: .000000).
    *** eigenvalue  10 =     .4605248       (converged: .000000).
    *** eigenvalue   9 =     .4303958       (converged: .000000).
    *** eigenvalue   8 =     .0973733       (converged: .000000).
    *** eigenvalue   7 =     .0943757       (converged: .000000).
    *** eigenvalue   6 =     .0141633       (converged: .000004).
    *** eigenvalue   5 =     .0046807   (NOT converged: .004079).
    *** eigenvalue   4 =     .0020023   (NOT converged: .010100).
    *** eigenvalue   3 =     .0010310   (NOT converged: .313417).
    *** eigenvalue   2 =     .0009886   (NOT converged: .933851).
    *** eigenvalue   1 =     .0006303   (NOT converged: .171971).

    ****************************************************************
     harmonic frequencies [cm**-1]:

           129.0591        161.6295        165.0552        230.0241
           351.6915        611.7668       1579.1898       1604.0732
          3372.3938       3488.4362       3555.9794       3581.9733

    ****************************************************************

Raman
-----

Ref: [apmp]_

Upon the specification of the RAMAN keyword, the polarizabilities of the system
are calculated by evaluating the response to an applied electrical field. The
calculation is done by means of the Berry phase approach, which is also suited
for periodic systems.

Nuclear Magnetic Resonance
--------------------------

**Theory**

Ref: [dsmp]_

A magnetic field :math:`\mathbf{B}` is applied to the system, which reacts by
induced electronic ring currents. These currents produce an additional magnetic
field by themselves, which is not homogeneous in space. Therefore, the actual
magnetic field at the ionic positions is different for all atoms in the cell.
This field determines the resonance frequency of the nuclear spin, and this
resonance can be measured with a very high accuracy.

The perturbation Hamiltonian is given by

.. math::
   :label: eqn-nmr

   \begin{aligned}
     H^{(1)}&=  \frac 12 \frac em \mathbf{p}\times \mathbf{r}\cdot \mathbf{B}.
   \end{aligned}

The difficult part of this Hamiltonian lies in the position operator which is
ill defined in a periodic system. To get around this, the wavefunctions are
localized and for each localized orbital, Eq. (:eq:`eqn-nmr`) is
applied individually assuming the orbital being isolated in space. Around each
orbital, a *virtual cell* is placed such that the wavefunction vanishes at the
borders of that virtual cell.

The perturbation and therefore also the response are purely imaginary, so that
there is no first order response density. This simplifies the equations and
speeds up convergence.

**NMR input**

The options which control the specific NMR features are discussed below. None
of them requires an argument, they can all be put in the same line.

-  **RESTART:** The run is restarted from a previous stop. The user has to take
   care that the required files RESTART.\ *xxx* (where *xxx* are NMR, p_x, p_y,
   p_z) exist.

-  **CURRENT:** Three density files containing current density values are
   written to disk. Further, the nucleus-independent chemical shift fields are
   written to disk. All in cube format.

-  **NOSMOOTH:** At the border of each virtual cell, the position operator is
   smoothened through an :math:`\exp -r^2`. This option turns smoothing off.
   Can be useful if your cell is too small so that this smoothing would already
   occur in a region where the orbital density has not yet vanished.

-  **NOVIRTUAL:** If this keyword is specified, no individual virtual cells are
   used at all. All orbitals will have the same virtual cell, equal to the
   standard unit cell.

-  **PSI0:** With this option, the Wannier orbitals are plotted as cube files.

-  **RHO0:** With this option, the orbital densities are plotted as cube files.

-  **OVERLAP:** Overlapping orbitals (overlap criterion read from next line)
   are assigned the same virtual cells. Useful for isolated systems.

-  **FULL** A full calculation of the :math:`\Delta j` term in the NMR scheme
   is done. Relatively expensive, but quite important for accurate results in
   periodic systems.

**NMR output**

At the end of the six perturbation runs, the results are printed. This includes
the magnetic susceptibility and the chemical shieldings. For the chemical
shieldings, there are two values: the raw and the net ones. The raw shieldings
correspond to a molecule in vacuum, without the susceptibility correction,
whereas the net shieldings contains that correction for a spherical sample. It
consists of an additive constant to all eigenvalues, which is printed out at
the end of the net shieldings.

In more detail, the results are:

-  The magnetic susceptibility tensor in both SI and cgs units. This quantity
   is an extensive quantity, i.e. two molecules in the simulation box will give
   twice the susceptibility of one.

-  The raw shielding matrices of all atoms and its eigenvalues. Generally, the
   shielding tensor is not symmetric. To obtain unique eigenvalues, it is
   symmetrized (:math:`A_{ij} \mapsto (A_{ij}+A_{ji})/2`) before the
   diagonalization.

   All values are given in ppm, parts per million. Chemical shieldings are
   dimensionless quantities.

-  The principal values (eigenvalues) of the raw and net shielding tensors. As
   mentioned above, they only differ by an additive constant, the
   susceptibility correction for a spherical sample, which is printed out at
   the end of the list. The numbers are the shielding eigenvalues from most
   negative to most positive, the isotropic shielding (the average of the
   eigenvalues), and the anisotropy (the difference between the most positive
   and the average of the other two).

**What to look at?** If you search for the values which are peaked in a
spectrum, you have to take the isotropic shieldings (the **iso** column of the
output). If the system is a gas, take the raw shielding, if it is condensed
matter, take the net shielding. If your system is a molecule in vacuo, but the
experimentalists measure it in a solvent, add the susceptibility correction to
the raw shieldings by yourself.

**Why are my numbers so strange in absolute value?** One more point shall be
mentioned: For all nuclei except hydrogen, pseudopotentials screen the core
electrons. The chemical shielding, however, is very sensitive to core and
semi-core electrons. This can be corrected through a semi-empirical additive
constant in many cases. This constant still needs to be added to the values
given from the program. It depends on the nucleus, on the pseudopotential, and
on the xc-functional.

In other cases, the calculated chemical shieldings are completely meaningless
due to this deficiency. Then, you have to use a pseudopotential which leaves
out the semicore states such that they are correctly taken into account.
Example: carbon shieldings can be corrected very well through a constant
number, silicon shieldings cannot. For Si, you have to take the :math:`n=3`
shell completely into the valence band, requiring a cutoff larger than 300Ry.

**How to compare to experiment?** Usually, experimentalists measure the
*difference* between the resonance frequency of the desired system and that of
a reference system, and they call it :math:`\delta` (the **shift**) instead of
:math:`\sigma` (the **shielding**). To make life more complicated, they usually
define the shift of nucleus A of molecule X with respect to reference molecule
ref as
:math:`\delta^{\mathsf{A}}_{\mathsf{ref}}(\mathsf{X}) = \sigma^{\mathsf{A}}(\mathsf{ref}) - \sigma^{\mathsf{A}}(\mathsf{X})`.
Example: To calculate :math:`\delta^{\mathsf{H}}_{\mathsf{TMS}}(\mathsf{CH}_4)`, where TMS=tetramethylsilane, the standard reference molecule for
H-shifts, one would have to calculate the H-shielding of TMS and of
CH\ :math:`_4` and subtract them. Unfortunately, TMS is nontrivial to
calculate, because it is a large molecule and the geometry is complicated (and
the shieldings probably must be calculated taking into account vibrational and
rotational averaging). Thus, in most cases it is better to take for instance
the CH\ :math:`_4` shielding as a (computational) reference, and transform the
shieldings relative to CH\ :math:`_4` to those relative to TMS through the
experimental shielding of CH\ :math:`_4` with respect to TMS.

While doing so, you should not forget that the shielding is a property which is
usually not yet fully converged when energies and bonding are. Therefore, the
reference molecule should be calculated with the same computational parameters
as the desired system (to reproduce the same convergence error). In particular,
computational parameters include the type of the pseudopotential and its
projectors, the xc-functional and the cutoff.

**What accuracy can I expect?** This is a difficult question, and there is no
overall answer. First, one has to consider that on the DFT-pseudopotential
level of theory, one will never reach the accuracy of the quantum chemistry
community. However, for “normal” systems, the absolute accuracy is typically
0.5-1 ppm for hydrogen and about 5-20ppm for Carbon. The spread between extreme
regions of the carbon spectrum is not reached: instead of 200ppm, one only
reaches 150ppm between aliphatic and aromatic atoms, for instance. The
anisotropy and the principal values of the shielding tensor can be expected to
be about 10-20% too small. For hydrogen shieldings, these values are usually
better, the error remains in the region of a single ppm.

FUKUI
-----

Compute, within the linear response perturbation theory, the nuclear Fukui
function :math:`\phi _I` [fukui]_ formally identified as a
reactivity index of the density functional
theory [fukui2]_ [fukui3]_, according to the postulated criterion for
:math:`\delta \mu`. The quantity :math:`\delta \mu` is the change in the
electronic chemical potential :math:`\mu` and is given by

.. math::

   \delta \mu =\int f({\bf r}) \delta V^{ext}({\bf r}) d^3r
   = -\sum_I \phi_I \delta {\bf R}_I

where :math:`\phi _I =(\partial {\bf F}_I/\partial N)_{V^{ext}}=
-(\partial \mu /\partial R_I)_N`, :math:`N` is the number of electrons,
:math:`f({\bf r})` the electronic Fukui function [fukui]_ [fukui2]_,
:math:`V^{ext}({\bf r})` the external potential at **r**, :math:`{\bf R}_I` the
Cartesian coordinate of the :math:`I^{th}` nucleus and :math:`{\bf F}_I` the
force on the :math:`I^{th}` nucleus.

KPERT: kdp k-point calculations
-------------------------------

Described in the entry of the :kw:`KPERT` keyword. Ref: [mimp]_

.. [1]
   except for isolated systems, **SYMMETRY=0**, where only the NMR part is
   adapted to.

.. [2]
   The reader might wonder whether Einstein summations over repeated indices
   like :math:`\alpha` and :math:`R` are done or not. The answer is that *it
   depends*. If only one atom is displaced along one of the axes, there is no
   summation needed. The generalization to a simultaneous displacement of all
   atoms in arbitrary directions is trivially done by assuming the summation
   over :math:`\alpha` and :math:`R`.

.. [3]
   In contrast to this, the :kw:`VIBRATIONAL ANALYSIS`
   displaces each atom in each direction first by +0.01 and then by -0.01.
