Fragment Orbital DFT (FO-DFT)
=============================

The basic idea of FO-DFT, as formulated by Senthilkumar, Grozema, Bickelhaupt,
and Siebbeles [Senthilkumar03]_, is to construct the two reference
states from calculations of the non-interacting isolated donor and acceptor
groups, respectively. This way one avoids spurious delocalization of electrons
over both groups, but also neglects polarisation effects of the donor on the
acceptor and vice versa. Thus the diabatic states - and therefore the matrix
element :math:`H_{AB}` - will in general be slightly different from those
calculated with CDFT.

With these diabatic reference states the transition matrix element is then in
principle given by:

.. math:: H_{AB}= \left\langle \psi^D_\text{HOMO}\left\vert \mathcal{H}_\text{KS} \right\vert\psi^A_\text{HOMO}\right\rangle,

where :math:`\mathcal{H}_\text{KS}` denotes the Kohn-Sham Hamiltonian and
:math:`\psi^{D,A}_\text{HOMO}` the highest occupied molecular orbitals (HOMO’s)
of donor and acceptor, respectively. Here we make the assumption that we can -
to a sufficient degree of accuracy - approximate the acceptors LUMO (lowest
unoccupied molecular orbital) by its HOMO in a system with one more electron.
Successful applications of our implementation of FO-DFT can be found in
ref.\ [Oberhofer10acie]_.

FODFT with CPMD
---------------

An FO-DFT calculation in CPMD consists of six different calls to CPMD:

1 and 2:
   the respective wavefunction optimisation of donor and acceptor, on positions
   as they would be if they were combined (**CENTER MOLECULE OFF**). Note that
   here the calculated acceptor state should be the desired acceptor state plus
   one electron, which is necessary because CPMD does not calculate Kohn-Sham
   matrix elements for unoccupied orbitals. These two calls are the only time
   consuming ones in the procedure.

3 and 4:
   the diagonalization of the two states in order to get their respective
   Kohn-Sham Orbitals (KOHN-SHAM ENERGIES).

5:
   the combination and orthogonalisation of both wavefunctions with keywords
   COMBINE WAVEFUNCTION and one of the orthogonalisation schemes
   (ORTHOGONALIZATION LOWDIN of GRAM-SCHMIDT).

6:
   the calculation of the Kohn-Sham matrix for the combined and orthogonalised
   states (keyword KSHAM).

If both donor and acceptor atoms are set the constraint value
:math:`N_\text{c}` denotes as usual the charge difference between donor and
acceptor. If NDON\ :math:`=0`, :math:`N_\text{c}` is the desired charge of the
acceptor.

TODO’s and usage WARNINGS
~~~~~~~~~~~~~~~~~~~~~~~~~

-  WARNING: don’t use diagonalization schemes (e.g. Lanczos) during the MD
   because it breaks the force calculation.

-  WARNING: don’t use non-orthogonal supercells, only SYMMETRY 0, 1 or 8 should
   be used.

-  WARNING: CDFT and FODFT are not implemented for QM/MM calculations, yet.
