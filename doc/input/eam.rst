.. _input-eam:

``&EAM`` section
================

The ``&EAM ... &END`` section contains the parameters of the embedded-atom
model (EAM) potentials used for species that are labelled ``_EAM_`` on their
``*`` line in the ``&ATOMS`` section. The section is read (routine ``eamin``)
during the initialization of the pseudopotentials whenever at least one
species carries the ``_EAM_`` label; if the section is missing the run stops
with ``No &EAM section``. EAM species carry no electrons; for them the
pseudopotential file only needs ``&ATOM`` (``Z=`` and ``ZV=`` lines),
``&INFO`` and a section spelled exactly ``&Potential`` (mixed case) that
contains one number, read as the core radius ``rc`` of an analytic local
potential of von Barth-Car type (if ``rc`` < 1e-10 the internal flag
``igau`` is set to -1). The LMAX line in ``&ATOMS`` must be present but is
ignored.

Syntax::

   &EAM
     re  fe  phie  alpha  beta  gamma  rhoe  PHIe  Ec  rc  rs   (species 1)
     ...                                                        (species NSP)
     i  j                                                       (pair list)
     ...
   &END

The first ``NSP`` lines (one for *every* species of ``&ATOMS``, in the same
order, including non-EAM species) contain eleven numbers each, read in this
order:

========== ==================================================================
``re``     equilibrium distance :math:`r_e` (Bohr)
``fe``     density prefactor :math:`f_e` (dimensionless)
``phie``   pair-potential prefactor :math:`\phi_e` (eV)
``alpha``  exponent :math:`\alpha` (dimensionless)
``beta``   exponent :math:`\beta` (dimensionless)
``gamma``  exponent :math:`\gamma` (dimensionless)
``rhoe``   reference density :math:`\rho_e` (dimensionless)
``PHIe``   embedding prefactor :math:`\Phi_e` (eV)
``Ec``     cohesive energy :math:`E_c` (eV)
``rc``     cutoff radius :math:`r_c` (Bohr)
``rs``     smoothing width :math:`r_s` of the cutoff function (Bohr)
========== ==================================================================

The energies are converted internally from eV to Hartree. After these lines
up to ``NSP*NSP`` lines with two species indices ``i j`` are read; reading
stops at the first line that cannot be read as two integers (e.g. ``&END``).
Each pair is symmetrized, and the EAM density and pair interaction between
atoms of species ``i`` and ``j`` are only computed for listed pairs whose two
species are both EAM species. Self pairs (e.g. ``1 1``) must be listed
explicitly; otherwise atoms of that species do not interact at all.

As implemented in ``eam_pot`` (``src/eam_pot_utils.mod.F90``), with
:math:`f_c(r) = 1/(\exp[(r-\bar r_c)/\bar r_s]+1)` (:math:`\bar r_c,\bar r_s`
being the averages over the two species) the density at atom :math:`i` is
:math:`\rho_i = \sum_j f_c(r_{ij})\, f_{e,j} \exp[-\beta_j (r_{ij}/r_{e,j}-1)]`,
the embedding energy is
:math:`F(\rho) = -E_c (1-\ln x)\,x - \Phi_e\, y` with
:math:`x=(\rho/\rho_e)^{\alpha/\beta}` and :math:`y=(\rho/\rho_e)^{\gamma/\beta}`,
and the pair term is :math:`\phi_e \exp[-\gamma (r/r_e-1)]` for equal species
and :math:`\tfrac12\,[f_b(r)/f_a(r)\,\phi_{aa}(r) + f_a(r)/f_b(r)\,\phi_{bb}(r)]`
for unlike species :math:`a, b`, where :math:`f_a(r)` and :math:`\phi_{aa}(r)`
are the density and pair functions of species :math:`a`; all terms are
multiplied by :math:`f_c(r)`.

Example layout (species 1 is an EAM species, species 2 is not and only
needs a placeholder line; the values are symbolic, not a parameter set)::

   &EAM
     <re> <fe> <phie> <alpha> <beta> <gamma> <rhoe> <PHIe> <Ec> <rc> <rs>
     0 0 0 0 1 0 1 0 0 0 1
     1 1
   &END

.. todo::

   MISSING: literature reference for this EAM functional form, recommended
   parameter sets, and confirmation of the units of ``rs`` (printed as a.u.
   by the code). The source comment in ``eamin`` lists the parameters in a
   different order ("rc re fe ...") than the actual READ statement; the
   order given above follows the READ statement.

.. source: src/eam_pot_utils.mod.F90:85
