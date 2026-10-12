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

.. keyword:: EMPIRICAL CORRECTION ... END EMPIRICAL CORRECTION
   :section: &VDW

   Block with the parameters of the empirical van der Waals correction. It is
   only read if :kw:`VDW CORRECTION` is set in the ``&CPMD`` section. The block
   ends at the first line containing the word ``END``. Unknown lines inside the
   block are silently ignored. The following sub-keywords are recognised:

   ``ALL DFT-D2``
      Use Grimme's DFT-D2 correction with the built-in :math:`C_6` and
      :math:`R_0` tables for all pairs of species (elements H to Rn; heavier
      elements up to Z=99 give a zero contribution with a warning).
   ``ALL DFT-D3``
      Use Grimme's DFT-D3 correction (see the ``DFT-D3 ...`` sub-keywords).
   ``CUSTOM``
      The number of pairs ``NVDW`` is read from the next line, followed by
      ``NVDW`` lines ``TYPE ia ib C6 R0 d`` where ``ia``, ``ib`` are species
      indices (order of the ``&ATOMS`` section) and ``TYPE`` is ``C6`` (damped
      dispersion of Elstner et al.), ``TANH`` or ``DFT-D2``. For ``DFT-D2`` the
      three numbers may be omitted to use the built-in values. ``DFT-D3`` is not
      allowed here.
   ``CUTOFF`` (alias ``VDW-CUTOFF``)
      Short-range cutoff, read from the next line. **Default** :math:`10^{-2}`.
   ``CELL`` (alias ``VDW-CELL``)
      Three integers ``NX NY NZ`` (number of additional cells in each direction
      for the real-space sum) are read from the next line. **Default** 0 0 0.
   ``S6GRIMME``
      Functional label for the DFT-D2 global scaling factor :math:`s_6`, read
      from the next line: ``BLYP`` (1.20), ``B3LYP`` (1.05), ``BP86`` (1.05),
      ``TPSS`` (1.00), ``PBE`` (0.75), ``REVPBE`` (1.25), ``PBE0`` (0.60). If not
      given, :math:`s_6` is taken from the functional of the ``&DFT`` section;
      the program stops if this is not one of the labels above.
   ``DFT-D3 PARAMETERS [THREEBODY] [NUMGRAD]``
      The real-space cutoffs for the dispersion energy and for the coordination
      numbers are read from the next line (Bohr, or Angstrom with
      :kw:`ANGSTROM`). **Defaults** :math:`\sqrt{9000}` and :math:`\sqrt{1600}`
      Bohr. ``THREEBODY`` adds the three-body term, ``NUMGRAD`` uses numerical
      gradients.
   ``DFT-D3 FUNCTIONAL [VERSION=n] [TZ]``
      The functional name as understood by the DFT-D3 library (e.g. ``b-lyp``,
      ``b3-lyp``, ``pbe``, ``pbe0``) is read from the next line. ``VERSION``
      selects the damping: 2 (Fermi-type), 3 or 5 (zero damping), 4 or 6
      (Becke-Johnson damping). ``TZ`` selects the special TZ parameter set.
      **Default** is ``VERSION=4`` and the functional derived from the
      ``&DFT`` section (BLYP, B3LYP, BP86, TPSS, PBE, REVPBE, PBE0, HSE06,
      REVPBE0, CAM-B3LYP, HCTH, PW91, OLYP, PBES, OPBE; otherwise the program
      stops).
   ``DFT-D3 NO_CONTRIBUTION``
      Species indices that are excluded from the DFT-D3 correction are read from
      the next line.

   Example::

      &VDW
       EMPIRICAL CORRECTION
        ALL DFT-D3
        DFT-D3 FUNCTIONAL VERSION=4
         pbe
       END EMPIRICAL CORRECTION
      &END

   .. todo::

      MISSING: functional form of the TANH type; meaning of the TZ parameter
      set and of VERSION 5 and 6 beyond the damping family printed by the code.

.. source: src/vdwin_utils.mod.F90:143 (sub-parser empvdwin, 410-760;
.. s6 table set_s6 1285; D3 mapping dftd3func 1325; damping report 316-330)

.. keyword:: GRIMME CORRECTION
   :section: &VDW

   Read the input for the external Grimme van der Waals library. Only available
   if CPMD was built with that library (CMake option ``CPMD_ENABLE_VDWLIB``,
   ``configure.sh`` option for vdW, preprocessor flag ``__HAS_LIBGRIMMEVDW``) and
   only read if :kw:`VDW CORRECTION` is set in the ``&CPMD`` section. The lines
   following the keyword are parsed by the library routine
   ``vdw_grimme_read_input``; the functional passed to the library is ``PBE`` or
   ``revPBE`` if the corresponding exchange gradient correction is used, and the
   functional label of the ``&DFT`` section otherwise. When present, this input
   replaces the built-in :kw:`EMPIRICAL CORRECTION ... END EMPIRICAL CORRECTION`
   scheme.

   .. todo::

      MISSING: input format of the library (vdw_grimme_read_input is not part of
      the CPMD source tree), terminating keyword, literature reference.

.. source: src/vdwin_utils.mod.F90:155 (guarded by #ifdef _HAS_LIBGRIMMEVDW;
.. macro mapping src/cpmd_global.h:89)

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

   If the :kw:`DCACP` van der Waals-correction scheme is adopted, the
   indices of atomic species (in the same order as specified in the &ATOMS
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
