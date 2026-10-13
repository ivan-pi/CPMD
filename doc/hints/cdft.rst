.. _sec:CDFT_FO-DFT:

CDFT Theory
===========

We implemented the constrained DFT (CDFT) method as developed by Wu and van
Voorhis [Wu05]_ [Wu06jcp]_ [Wu06jctc]_ [Wu06jpca]_ [Wu07]_.

The density constraint
----------------------

Here, one imposes a scleronomic constraint on the electron density to reproduce
a certain charge distribution on the atoms,

.. math:: \int w(\mathbf{r}) \rho(\mathbf{r})\; d\mathbf{r}-N_\text{c}=0,

with :math:`\rho(r)` being the electron density :math:`w(r)` the weight and
:math:`\text{N}_\text{c}` the constraint value. The new energy functional:

.. math:: W[\rho,V]=E[\rho]+V (\int w(\mathbf{r}) \rho(\mathbf{r})\; d\mathbf{r}-N_\text{c}),

with :math:`E[\rho]` being the normal Kohn-Sham energy functional.

The weight
----------

We chose the Hirshfeld partitioning scheme [Hirshfeld77]_ in order
to define charges. The weight for imposing a charge difference between a donor
:math:`D` and an acceptor :math:`A` is then given by:

.. math:: w(\mathbf{r})=\frac{\sum_{i \in D} \rho_i(\mathbf{r}-\mathbf{R}_i)-\sum_{i \in A} \rho_i(\mathbf{r}-\mathbf{R}_i)}{\sum_{i=1}^{N} \rho_i(\mathbf{r}-\mathbf{R}_i)},

where the sums in the numerator range over all donor and acceptor atoms,
respectively, whereas the sum in the denominator ranges over all :math:`N`
atoms. :math:`\rho_i(r)` is the unperturbed electron density of atom :math:`i`
which is given by

.. math:: \rho_i(r)=\sum_j n_j \frac{\vert \psi_i^j(r)\vert^2}{r^2},

where the sum ranges over all orbitals and :math:`\psi_i^j(r)` and :math:`n_j`
denote the reference orbitals and :math:`n_j` their occupation number,
respectively.

Constraint forces
-----------------

Using the Hellmann-Feynman theorem we can express forces on atom :math:`k` at
position :math:`R_k`

.. math:: \frac{\partial W}{\partial \mathbf{R}_k}=\langle \psi\vert \frac{\partial \hat{H}}{\partial \mathbf{R}_k}\vert \psi\rangle.

Applying this relation we can calculate the additional forces due to the bias
potential:

.. math::

   \begin{aligned}
   \mathbf{F}_{\text{Bias},k}= - V \int \rho(\mathbf{r})\frac{\partial w(\mathbf{r}-\mathbf{R}_k)}{\partial \mathbf{R}_k}\;d\mathbf{r},
   \end{aligned}

with

.. math:: \frac{\partial w(\mathbf{r}-\mathbf{R}_k)}{\partial \mathbf{R}_k}=-\frac{\rho_k^\prime(\vert\mathbf{r}-\mathbf{R}_k\vert)}{\sum \rho_i(\vert\mathbf{r}-\mathbf{R}_k\vert)} G_k(\mathbf{r}-\mathbf{R}_k)

.. math::

   \begin{aligned}
   G_k(\mathbf{r}-\mathbf{R}_k) = \left\lbrace
   \begin{array}{l c r}
   w(\mathbf{r}-\mathbf{R}_k)-1 &  k \in D\\
   w(\mathbf{r}-\mathbf{R}_k)+1 &  k \in A\\
   w(\mathbf{r}-\mathbf{R}_k) &  k \not\in D\cup A\\
   \end{array}\right.
   \end{aligned}

Finally the derivative of :math:`\rho_i` is given by

.. math::

   \begin{aligned}
   \rho_k^\prime(\vert\mathbf{r}-\mathbf{R}_k\vert)&= \frac{\partial\rho_k(\vert\mathbf{r}-\mathbf{R}_k\vert)}{\partial\mathbf{R}_k}= \frac{\partial\rho_k(\vert\mathbf{r}-\mathbf{R}_k\vert)}{\partial \vert\mathbf{r}-\mathbf{R}_k\vert}\frac{\partial\vert\mathbf{r}-\mathbf{R}_k\vert}{\partial\mathbf{R}_k},\nonumber\\
   &=  \frac{\partial\rho_k(\vert\mathbf{r}-\mathbf{R}_k\vert)}{\partial \vert\mathbf{r}-\mathbf{R}_k\vert}\frac{\mathbf{r}-\mathbf{R}_k}{\vert\mathbf{r}-\mathbf{R}_k\vert}.
   \end{aligned}

The radial partial derivative :math:`\rho_i` of is finally calculated
numerically using splines.

For a more thorough treatment of the topic of constrained DFT MD please consult
reference [Oberhofer09]_.

Cutoff correction
~~~~~~~~~~~~~~~~~

In order to avoid integrating over every real space gridpoint we introduced a
cutoff :math:`R_c` in the generation of the weights. :math:`R_c` is chosen for
each atom species such that the total reference density is smaller than
:math:`10^{-6}`.

The action of the cutoff can be described by a Heaviside function
:math:`\theta`

.. math:: \rho_k(\vert\mathbf{r}-\mathbf{R}_k\vert)\rightarrow\rho_k(\vert\mathbf{r}-\mathbf{R}_k\vert)\theta(R_c-\vert\mathbf{r}-\mathbf{R}_k\vert).

Therefore the derivative of the reference density becomes

.. math:: \rho_k^\prime(\vert\mathbf{r}-\mathbf{R}_k\vert)\rightarrow\rho_k^\prime(\vert\mathbf{r}-\mathbf{R}_k\vert)\theta(R_c-\vert\mathbf{r}-\mathbf{R}_k\vert)+\rho_k(\vert\mathbf{r}-\mathbf{R}_k\vert)\delta(R_c-\vert\mathbf{r}-\mathbf{R}_k\vert),

with the Dirac :math:`\delta` function. Thus the full force splits up into

.. math:: \mathbf{F}_{\text{Full},k}=\mathbf{F}_{\text{Bias},k}+\mathbf{F}_{\text{Bound},k}.

Here :math:`\mathbf{F}_{\text{Bound},k}` denotes forces due to the finite
cutoff. They can be expressed by the following surface integral in spherical
coordinates

.. math::

   \begin{aligned}
   \mathbf{F}_{\text{Bound},k}=-V_c \rho_k(R_c) R_c \int \frac{\rho(R_c,\vartheta,\varphi)G(R_c,\vartheta,\varphi)}{\sum \rho_i(R_c,\vartheta,\varphi)}R_c
   \left(
   \begin{array}{c}
   \sin\vartheta\cos\varphi\\
   \sin\vartheta\sin\varphi\\
   \cos\vartheta
   \end{array}
   \right)\sin\vartheta \;d\vartheta d\varphi
   \end{aligned}

As we do have a Cartesian grid we need to perform the integral as an integral
over a thin shell in real space. :math:`R_c` times the vector is just the
position vector of a point on the surface :math:`(x,y,z)` and the surface
element can be expressed in Cartesian coordinates as

.. math:: \sin\vartheta\;d\vartheta d\varphi = \mathrm{sgn}(z) \frac{y\;dx dz-x\;dy dz}{R_c(R_c^2-z^2)}

Transition Matrix Element calculation
-------------------------------------

The following technical description of CDFT matrix element calculations can
also be found, together with some test calculation and an investigation of the
influence of the choice of the weight function on the results, in
ref. [Oberhofer10acie]_.

In order to calculate the electronic transition matrix element we first need to
calculate two constrained wave functions :math:`\phi_A` and :math:`\phi_B` and
write down the Hamiltonian in the non-orthogonal constrained basis spanned by
these two states.

Hamiltonian in the non-orthogonal constrained basis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: \mathbb{H}_\text{no}=\left(\begin{array}{c c} H_{BB} & H_{BA}\\ H_{AB} & H_{AA}\end{array}\right),

where :math:`A` and :math:`B` denote the two states, respectively. Here the
diagonal elements are just given by the respective energies of donor and
acceptor due to:

.. math::

   \begin{aligned}
   H_{BB}&= \langle \phi_B \vert H \vert \phi_B \rangle=E_B \\
   H_{AA}&= \langle \phi_A \vert H \vert \phi_A \rangle=E_A,
   \end{aligned}

while the off-diagonal elements are determined by:

.. math::

   \begin{aligned}
   H_{BA}&= F_B \langle \phi_B \vert \phi_A \rangle -V_c^B \langle \phi_B \vert W\vert \phi_A \rangle\\
   &=  F_B S_{BA} -V_c^B W_{BA}\\
   H_{AB}&= F_B S_{AB} -V_c^A W_{AB}.
   \end{aligned}

Here we also introduced the overlap matrix element :math:`S_{AB}` and the
weight matrix element :math:`W_{AB}`. From now on we will only consider one of
the off-diagonal elements as the calculation of the other one is analogous.
Approximating the wave functions of the two states by slater determinants built
from Kohn-Sham orbitals we can write in a plane wave basis:

.. math::

   \begin{aligned}
   S_{AB}=\langle \phi_B \vert \phi_A \rangle= \det \left[\sum_G (c_i^A(G))^\star c_j^B(G)\right]=\det \Phi_{ij},
   \end{aligned}

where the :math:`c_{i,j}^{A,B}(G)` denote the plane wave coefficients of the
respective wave function for their respective electronic states :math:`i` and
:math:`j` and :math:`\Phi` denotes the full overlap matrix. Note that the
diagonal elements of :math:`\mathbb{S}` are equal to one due to the
normalisation of the two states. The off-diagonal :math:`\mathbb{W}` matrix
element is given by:

.. math::
   :label: eq-wab

   W_{AB}=\langle \phi_B \vert W \vert \phi_A \rangle=N\sum_{i=1}^N \sum_{j=1}^N \langle \phi^i_B \vert W \vert \phi^j_A \rangle  (-1)^{i+j} \det \Phi_{(i,j)},

with :math:`\Phi_{(i,j)}` being the :math:`i,j`\ th minor of :math:`\Phi` and
the integral given by

.. math::
   :label: eq-gspacetransme

   \langle \phi^i_B \vert W \vert \phi^j_A \rangle=\sum_{G,G^\prime}(c_i^A(G))^\star c_j^B(G) \tilde{W}(G-G^\prime).

Note that Equ. :eq:`eq-gspacetransme` actually constitutes a
convolution in G-space. Therefore the most efficient way to calculate the
elements :math:`\langle \phi^i_B \vert W \vert \phi^j_A \rangle` is by
transforming Equ. :eq:`eq-gspacetransme` back into real-space where
the convolution becomes a simple multiplication. Additionally CPMD already
stores the wavefunctions in real-space (REAL SPACE WFN KEEP which is activated
automatically by the CDFT HDA calculation) and thus saves us the necessary
FFTs. Calculating every minor determinant of the overlap matrix in
Equ. (:eq:`eq-wab`) would be rather costly, therefore we note that the
second part of Equ. (:eq:`eq-wab`) is actually the cofactor matrix
:math:`C` of :math:`\Phi`.

.. math:: (-1)^{i+j} \det \Phi_{(i,j)}=C_{ij},

which can be calculated using the Laplace Expansion of the inverse of
:math:`\Phi`:

.. math:: C^T=\Phi^{-1}. [\det(\Phi) I ]

Thus we only have to calculate the inverse and the determinant of :math:`\Phi`
instead of all the minors. The diagonal elements of :math:`\mathbb{W}` are just
the constraint values :math:`N_c` for the two states. Having performed this
calculation we symmetrise :math:`\mathbb{H}_\text{no}` in order to correct for
inaccuracies of the approximations we made.

Full diabatic non-orthogonal Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

However, the matrix :math:`\mathbb{H}_\text{no}` is not the full Hamiltonian as
it neglects the fact that the two diabatic wavefunctions form a non-orthogonal
basis (:math:`S_{AB}\neq 0`). Nevertheless, we can construct the full
Hamiltonian using the overlap matrix element :math:`S_{AB}`:

.. math:: \mathbb{H}_\text{full}=\frac{1}{1-S_{AB}^2}\left(\begin{array}{c c} H_{BB}-S_{AB}H_{AB} & H_{BA}-S_{AB}H_{AA}\\ H_{AB}-S_{AB}H_{BB} & H_{AA}-S_{AB}H_{BA}\end{array}\right),

Diabatic orthogonal Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to compare our matrix elements with other methods an experiments we
follow the procedure of Wu and Van Voorhis and first solve the two-dimensional
generalised eigenvalue problem

.. math:: \mathbb{W} \mathbb{V}=\mathbb{S}\mathbb{V}\mathbb{L},

where :math:`\mathbb{V}` is the matrix of generalised eigenstates and
:math:`\mathbb{L}` is the diagonal matrix of generalised eigenvalues. To
calculate the diabatic orthogonal Hamiltonian we then perform the similarity
transformation

.. math:: \mathbb{H}_\text{diab}=\mathbb{V}^{-1} \mathbb{H}_\text{full} \mathbb{V}.

The off-diagonal elements of :math:`\mathbb{H}_\text{diab}` –or their average
if the two states have been very different– are then the desired transition
matrix element.

Adiabatic Hamiltonian
~~~~~~~~~~~~~~~~~~~~~

The adiabatic energies of ground and excited state :math:`\varepsilon_{g,e}`
and their mixing coefficients :math:`\mathbf{x}_{g,e}` can be calculated by
simply diagonalising the diabatic Hamiltonian (in principle either one as a
similarity transform does not change the eigenvalues of a matrix, here we use
the full Hamiltonian)

.. math:: \mathbb{H}_\text{full}\mathbf{x}_{g,e}=\varepsilon_{g,e}\mathbf{x}_{g,e}

With the matrix :math:`\mathbb{X}` of mixing coefficients we could then in
principle also calculate the adiabatic states themselves via

.. math::
   :label: eq-adiabS

   \left( \begin{array}{c} \psi^\text{ad}_a\\ \psi^\text{ad}_b\end{array}\right)=\mathbb{X} \left( \begin{array}{c} \psi_A\\ \psi_B\end{array}\right)

Projection on reference states
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unfortunately, CDFT only constrains a total charge density, thus sometimes one
may get a spurious delocalization of the orbitals over both donor and acceptor
leading to a too large overlap :math:`S_{AB}` of the two diabatic states. One
possible way to overcome this is to project the wavefunctions on
non-interacting reference states and transform the diabatic Hamiltonian into
this basis.

The reference states are constructed from four wavefunctions
:math:`\Phi_{D^-},\Phi_{D^+},` :math:`\Phi_{A^+},\Phi_{A^-}` describing donor
and acceptor in the charge states corresponding to the two diabatic states

.. math::

   \begin{aligned}
   \psi_a^0 &=  \Phi_{D^-}\times \Phi_{A^+} \nonumber\\
   \psi_b^0 &= \Phi_{D^+}\times \Phi_{A^-}.
   \end{aligned}

Then we project the adiabatic states - Equ. :eq:`eq-adiabS` - onto the
reference states

.. math::
   :label: eq-projection

   \left( \begin{array}{c} \psi^\prime_a\\ \psi^\prime_b\end{array}\right)=\left(
   \begin{array}{c c}
    \langle \psi^\text{ad}_a \vert \psi^0_a\rangle & \langle \psi^\text{ad}_b \vert \psi^0_a\rangle\\
    \langle \psi^\text{ad}_a \vert \psi^0_b\rangle & \langle \psi^\text{ad}_b \vert \psi^0_b\rangle
   \end{array} \right) \left( \begin{array}{c} \psi^\text{ad}_a\\ \psi^\text{ad}_b\end{array}\right) =
   \mathbb{K} \left( \begin{array}{c} \psi^\text{ad}_a\\ \psi^\text{ad}_b\end{array}\right).

These, so called “dressed”, states are not yet orthogonal which we rectify
using the Löwdin scheme

.. math:: \left( \begin{array}{c} \psi_1\\ \psi_2\end{array}\right)=\mathbb{S}^{-1/2} \left( \begin{array}{c} \psi^\prime_a\\ \psi^\prime_b\end{array}\right),

with the elements of the Löwdin matrix :math:`\mathbb{S}` defined as
:math:`\mathbb{S}_{j,k}=\langle
\psi_j^\prime \vert \psi_k^\prime \rangle` where :math:`j,k \in {a,b}`.
Combining this with Equ. :eq:`eq-projection` we can then express the
new projected diabatic states in terms of the adiabatic constrained states:

.. math:: \left( \begin{array}{c} \psi_1\\ \psi_2\end{array}\right)=\mathbb{S}^{-1/2} \mathbb{K} \left( \begin{array}{c} \psi^\text{ad}_a\\ \psi^\text{ad}_b\end{array}\right).

The similarity transformation from the adiabatic Hamiltonian matrix
:math:`\mathbb{H}_\text{ad}` to the projected diabatic Hamiltonian
:math:`\mathbb{H}^\text{p}_\text{diab}` is then simply

.. math:: \mathbb{H}^\text{p}_\text{diab}= (\mathbb{S}^{-1/2} \mathbb{K}) \mathbb{H}_\text{ad} (\mathbb{S}^{-1/2} \mathbb{K})^\text{T}.
