.. _input-vdw:

``&VDW`` section
================

Settings associated with van der Waals correction schemes. This section is only evaluated if one of the keywords :kw:`DCACP`, :kw:`VDW CORRECTION` or :kw:`VDW WANNIER` is given in the ``&CPMD`` section.

The section starts with ``&VDW`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &VDW

Keywords
--------

.. keyword:: DCACP Z=z
   :section: &VDW

   Set custom DCACP parameters for the element with core charge :math:`Z=`\ ``z``.
   The values of the DCACP projectors :math:`\sigma_1` and :math:`\sigma_2` are
   read from the next line.

   The keyword can be repeated multiple times, once for every element. Useful if
   the desired parameters are not available in the internally coded library,
   *e.g.* for an exotic xc functional.

.. keyword:: INCLUDE_METALS
   :target: INCLUDE METALS
   :section: &VDW

   When using the :kw:`DCACP` dispersion-correction scheme, include
   contributions from metal centres, too (where available). This should only
   rarely be needed. In the vast majority of cases (metal complexes,
   metallo-organic compounds), contributions from the metal centres are
   vanishingly small.

.. keyword:: NO_CONTRIBUTION
   :target: NO CONTRIBUTION
   :section: &VDW

   If the :kw:`DCACP` van der Waals-correction scheme is adpoted, the
   indeces of atomic species (in the same order as specified in the &ATOMS
   section) that shall not contribute to the dispersion forces (capping atoms,
   metals...) are read from the next line.

.. keyword:: VDW DCACP
   :section: &VDW

   Synonymous to :kw:`DCACP`.

.. keyword:: VDW PARAMETERS
   :section: &VDW

   Parameters for empirical van der Waals correction schemes are set with the
   keyword. This requires the :kw:`VDW CORRECTION` keyword to be
   set in the &CPMD section. For Grimme’s **DFT-D2** and **DFT-D3** type (see
   below) an automatic assignment of the parameters can be requested by putting
   **ALL DFT-D2** or **ALL DFT-D3** on the next line. Otherwise the number of
   pairs *NVDW* is read from the next line and followed by *NVDW* lines of
   parameters: *TYPE*, :math:`\alpha`, :math:`\beta`, :math:`C_6^{\alpha\beta}`,
   :math:`R_0^{\alpha\beta}`, and :math:`d` for each pair of atom types
   :math:`\alpha` and :math:`\beta`, where :math:`\alpha` and :math:`\beta` are
   the indexes of pseudopotentials (and their associated groups of atoms) in the
   order they are listed in the &ATOMS section. For type **DFT-D2** and **DFT-D3**
   only :math:`\alpha` and :math:`\beta` are required. If the other parameters are
   omitted the internal table of parameters is used.

   A presently implemented damped dispersion model, described by M. Elstner *et
   al.*\ [Elstner]_, having the same form as that constructed by
   Mooij *et al.*\ [mooij:99]_, is activated by specifying **C6** as
   *TYPE*. This model is expressed as

   .. math::
      :label: elstner-damping-function

      E_{vdW} = \sum_{ij}
          \frac{C_6^{\alpha\beta}}{{R^{\alpha\beta}_{ij}}^6}
      \left(1 - \exp{
      \left[-d
      \left(\frac{R^{\alpha\beta}_{ij}}{R^{\alpha\beta}_0}
      \right)^7
      \right]}
      \right)^4.

   A table of parameters appropriate for this particular model, using the PBE and
   BLYP functionals, is available [williams-vdw:06]_.

   Alternatively Van der Waals correction according to Grimme can be used
   [Grimme06]_ by selecting *TYPE* **DFT-D2** and its D3 version
   [GrimmeD3]_ with **DFT-D3**.

   .. math::

      E_{disp} = - s_6 \sum_{i=1}^{N_{at} -1} \sum_{j=i+1}^{N_{at}}
                     \frac{C_6^{ij}}{R_{ij}^6} f_{dmp} (R_{ij})

   The values of :math:`C_6` and :math:`R_0` are not specific that are used by
   this method are taken from [Grimme06]_ and stored internally (see
   above for details). Namely, all elements from H (:math:`Z=1`) to Rn
   (:math:`Z=86`) are available, whereas elements beyond Rn give by default a zero
   contribution. Note that the parameter :math:`s_6` depends on the functional
   used and has to be provided consistently with the DFT one chosen for the
   calculation. The following line has to be added S6GRIMME and the type of
   functional is read from the next line. One of the following labels has to be
   provided: BP86, BLYP, B3LYP, PBE, TPSS, REVPBE, PBE0. Note that Grimme vdW does
   not support other functionals. In the **DFT-D3** version
   [GrimmeD3]_ there is no need to specify manually :math:`s_6`

.. keyword:: VDW-CELL
   :section: &VDW

   The number of additional supercells to be included in the sum of van der Waals
   correction.

.. keyword:: VDW-CUTOFF
   :section: &VDW

   On the next line the short range cutoff of van der Waals correction has to be
   specified. The default value is :math:`10^{-2}`.

.. keyword:: WANNIER CORRECTION ... END WANNIER CORRECTION
   :section: &VDW

   Between these opening and ending keywords, the partitioning of the system and
   the calculation procedure must be selected. Three implementations are available
   for partitioning the system: (1) choosing a *zlevel*, namely a z coordinate
   separating the first fragment form the second (this is appropriate for cases
   where there are only two fragments such as, for instance two graphene layers or
   adsorption of molecules on surfaces); in this case the keyword FRAGMENT ZLEVEL
   must be used. (2) give reference ion and a cut-off radius around which WFCs are
   supposed to belong to the given atom or fragment; in this case the keyword
   FRAGMENT RADIUS must be used. (3) the system is subdivided into fragments
   automatically detected by using predefined covalent bond radii. in this case
   the keyword FRAGMENT BOND must be used. This is also the default in case no
   specification is done.

   The syntax for the different options is:

   VERSION

   iswitchvdw (method 1 [psil1]_ or 2 [psil2]_)

   FRAGMENT ZLEVEL

   zlevel (in a.u.)

   FRAGMENT RADIUS

   multifrag

   i radius(i)

   ...

   FRAGMENT BOND

   tollength

   DAMPING [DIPOLE]

   a6

   RESTART WANNIER

   ENERGY MONOMER

   enmonomer

   TOLERANCE WANNIER

   tolwann

   TOLERANCE REFERENCE

   tolref

   CHANGE BONDS

   nboadwf

   i j :math:`\pm` 1

   CELL

   nxvdw nyvdw nzvdw

   PRINT :math:`[`\ INFO,FRAGMENT,C6,FORCES\ :math:`]`

   Note that the total number of WFCs in your system depends on the spin
   description you use (1 for LSD, 2 for LDA). The coefficient a6 is the smoothing
   parameter and the reference total energy intended as a sum of all the total
   energies of your fragments (e.g. the ETOT you get by a standard calculation not
   including vdW corrections). For a6 the suggested parameter is 20.0
   [molphy]_. An alternative method avoiding the use of the empirical
   parameter a6 is the one exploiting the replacement of the short-range damping
   function by an estimation of the Pauli exchange repulsion given in terms of
   solely Wannier function centers and spreads [psil3]_. To activate
   this option the keyword DAMPING has to be given along with DIPOLE on the same
   line. The following line (a6) has to be obviously removed, since anyhow is
   neither read nor used in this (DAMPING DIPOLE) case. Note that the two possible
   vdW options, EMPIRICAL CORRECTION and WANNIER CORRECTION are mutually
   exclusive.
