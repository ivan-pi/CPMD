.. _input-dft:

``&DFT`` section
================

Exchange and correlation functional and related parameters.

The section starts with ``&DFT`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &DFT

Keywords
--------

.. keyword:: ACM0
   :section: &DFT

   Add exact exchange according to the adiabatic connection method
   0. [acm0]_ [adamo2000]_ The fraction of exact exchange is set to 0.25 and
   the fractions of the LDA and gradient-corrected exchange to 0.75.

   .. warning::

      The ACM keywords only take effect with the legacy drivers when the
      functional is assembled with :kw:`GRADIENT CORRECTION` (and
      :kw:`SLATER`/:kw:`LDA CORRELATION`). Together with :kw:`FUNCTIONAL` the
      settings are overwritten by the functional definition, which is applied
      after the whole section has been read; in particular no exact exchange is
      computed for a non-hybrid functional such as ``PBE``. With
      :kw:`XC_DRIVER` the program stops; use :kw:`SCALES` and
      :kw:`HFX_SCALE` instead.

.. source: src/dftin_utils.mod.F90:896

.. keyword:: ACM1
   :section: &DFT

   Add exact exchange according to the adiabatic connection method
   1. [adamo2000]_ [acm1]_ The parameter :math:`a` is read from the next line;
   the fraction of exact exchange is set to :math:`a` and the fractions of the
   LDA and gradient-corrected exchange to :math:`1-a`. The same restrictions
   as for :kw:`ACM0` apply.

.. keyword:: ACM3
   :section: &DFT

   Add exact exchange according to the adiabatic connection method
   3. [adamo2000]_ [acm3]_ The three parameters :math:`a_1`, :math:`a_2` and
   :math:`a_3` are read from the next line; the fraction of exact exchange is
   :math:`a_1`, the fractions of the LDA exchange, gradient-corrected exchange
   and gradient-corrected correlation are :math:`1-a_1`, :math:`1-a_2` and
   :math:`1-a_3`. The same restrictions as for :kw:`ACM0` apply.

.. keyword:: ANALYTICAL_DIV
   :target: ANALYTICAL DIV
   :section: &DFT

   Calculate the :math:`\mathbf{G}=0` term for the Coulomb-attenuated exact
   exchange using an analytical description of the integral term. This is the
   default in combination with all internal hybrid xc functionals.

.. keyword:: BECKE BETA
   :section: &DFT

   Change the :math:`\beta` parameter in Becke’s exchange
   functional [Becke88]_ to the value given on the next line.

   **Default** is 0.0042.

.. keyword:: COULOMB ATTENUATION
   :section: &DFT

   Activates the long-range correction/Coulomb attenuation method (**LC, CAM**).
   With this keyword, all selected *internal* LDA or GGA exchange functionals will
   be Coulomb-attenuated, and the corresponding exact exchange contribution is
   added by default. The three CAM parameters are read from the next line in
   the order :math:`\mu` (range-separation parameter; printed by the code as
   "MU [BOHR]"), :math:`\alpha` and :math:`\beta`. Setting :math:`\alpha=0`
   and :math:`\beta=1.0` is equivalent to using the long-range correction (LC)
   method; :math:`\alpha+\beta` must not exceed 1.

   The keyword switches on the xc driver (:kw:`XC_DRIVER <XC DRIVER>`) by
   itself; the functionals must therefore be given with xc-driver names. The
   CAM parameters also apply to the functionals of :kw:`LR KERNEL` and
   :kw:`MTS_LOW_FUNC`. The aliases ``CAM``, ``ATTENUATION`` and
   ``RANGE SEPARATION`` are also accepted (but not on a line that contains
   ``SCREENED``, see :kw:`SCREENED EXCHANGE`).

   Example (CAM-B3LYP parameters)::

      COULOMB ATTENUATION
       0.33  0.19  0.46

.. source: src/dftin_utils.mod.F90:686 (parameter order: hfx_func_init in
.. src/cp_xc_utils.mod.F90:1881)

.. keyword:: CP_LIBRARY_ONLY
   :target: CP LIBRARY ONLY
   :section: &DFT

   Use the new xc driver and take all functionals of :kw:`FUNCTIONAL` from the
   internal CP library (libxc is not used).

   The keyword switches on the xc driver by itself, i.e. it implies
   :kw:`XC_DRIVER <XC DRIVER>`; a separate ``XC_DRIVER`` line is not needed.
   Consequently:

   * the functionals must be given with :kw:`FUNCTIONAL` using the xc-driver
     names (e.g. ``GGA_XC_PBE``); the legacy names (``PBE``, ``BLYP``, ...)
     are not valid;
   * :kw:`GRADIENT CORRECTION`, :kw:`LDA CORRELATION`, :kw:`SLATER` (except
     ``SLATER ALPHA``), :kw:`ACM0`, :kw:`ACM1` and :kw:`ACM3` stop the program;
   * :kw:`OLDCODE` and :kw:`NEWCODE` have no effect.

   The keyword sets the library of the functionals of :kw:`FUNCTIONAL` only,
   overriding a :kw:`LIBRARY` line given before it. If :kw:`LR KERNEL` is not
   given, the kernel is a copy of the functional and therefore also uses the
   internal library. A :kw:`LIBRARY`, :kw:`KERNEL_LIBRARY` or
   :kw:`MTS_LOW_FUNC_LIBRARY` entry ``LIBXC`` (or ``LXC``) given *after* this
   keyword stops the program. If both :kw:`CP_LIBRARY_ONLY` and
   :kw:`LIBXC_ONLY` are given, the last one is used.

   The aliases ``USE_CP`` and ``INTERNAL_ONLY`` are also accepted.

   Example::

      CP_LIBRARY_ONLY
      FUNCTIONAL GGA_XC_PBE

.. source: src/dftin_utils.mod.F90:319

.. keyword:: EXCHANGE CORRELATION TABLE
   :section: &DFT
   :alternatives: [NO]

   Specifies the range and the granularity of the lookup table for the local
   exchange-correlation energy and potential. The number of table entries and the
   maximum density have to be given on the next line.

   Note that this keyword is only relevant when using :kw:`OLDCODE` and
   even then it is set to **NO** by default. Previous default values were 30000
   and 2.0. The table is always switched off for hybrid functionals. Only the
   words ``EXCHANGE`` and ``TABLE`` are required on the keyword line.

.. source: src/dftin_utils.mod.F90:801

.. keyword:: FORCE_NEW_HFX
   :section: &DFT

   Force the use of the new Hartree-Fock exchange (HFX) module even when
   :kw:`HFX SCREENING` is given. By default the new HFX module is used, but any
   :kw:`HFX SCREENING` line switches back to the old module; with
   ``FORCE_NEW_HFX`` the new module is kept, so that thresholding of the
   integrals is done with the new module as well. The keyword takes no argument.
   The settings :kw:`HFX_BLOCK_SIZE` and :kw:`HFX_DISTRIBUTION` only apply to the
   new module.

   **Default** is to use the new HFX module unless :kw:`HFX SCREENING` is present.

   Example::

      HFX SCREENING EPS_INT
       1.0E-6
      FORCE_NEW_HFX

.. source: src/dftin_utils.mod.F90:988 (flag applied in check_options, line 2155;
.. module selection in src/hfx_drivers.mod.F90:124)

.. keyword:: FUNCTIONAL
   :section: &DFT
   :notes: functional_1, functional_2, ...

   | The xc functional(s) are specified on the same line, separated by a space.
   |

   There is no default functional: if neither :kw:`FUNCTIONAL` nor one of
   :kw:`GRADIENT CORRECTION`, :kw:`SLATER`, :kw:`LDA CORRELATION`,
   :kw:`HARTREE` or :kw:`HARTREE-FOCK` is given, the program stops with "NO
   FUNCTIONAL HAS BEEN DEFINED - DEFAULTS NO LONGER APPLY". The ``&DFT``
   section is therefore mandatory. The alias ``MTS_HIGH_FUNC`` is also
   accepted (see :kw:`MTS_HIGH_FUNC`).

   **XC_DRIVER**:

   Both functionals from the internal CP library and those linked from libxc are
   available within the same run and can be freely mixed when using the
   :kw:`XC_DRIVER <XC DRIVER>`. Please note that the internal CP functionals
   offer a considerable performance advantage over libxc. While the code will
   default to the internal implementations, use of libxc for single functionals
   can be requested by the keyword :kw:`LIBRARY`. The functional code
   names for the functionals available internally follow the nomenclature
   introduced by libxc: The functional name must be preceded by a tag specifying
   its level (hybrid, MGGA, GGA...). B3LYP would therefore be specified as
   HYB_GGA_XC_B3LYP.

   | Fixed (hybrid) exchange-correlation functionals available in the internal CP
     library include:
   |

   **GGA_XC** \_BLYP, \_BP86, \_HCTH_93, \_HCTH_120, \_HCTH_147, \_HCTH_407,
   \_OLYP, \_OPBE, \_PBE, \_PBE_SOL, \_REVPBE

   **MGGA_XC** \_TPSS, \_M06_L, \_REVM06_L, \_M11_L, \_MN12_L, \_MN15_L

   **HYB_GGA_XC** \_B3LYP, \_CAM_B3LYP, \_O3LYP, \_PBE0, \_PBE0_SOL, \_REVPBE0

   | **HYB_MGGA_XC** \_M05, \_M05_2X, \_M06, \_M06_2X, \_M06_HF, \_M08_SO,
     \_M08_HX, \_M11, \_MN12_SX
   |

   The following aliases are also accepted:

   .. list-table:: Aliases of xc-driver names (internal library)
      :header-rows: 1
      :widths: 30 70

      * - Name
        - Meaning
      * - ``GGA_XC_BP``
        - Alias of ``GGA_XC_BP86``.
      * - ``GGA_XC_PBE_R``
        - Alias of ``GGA_XC_REVPBE``.
      * - ``HYB_GGA_XC_PBE0_R``
        - Alias of ``HYB_GGA_XC_REVPBE0``.
      * - ``GGA_X_PBE_R``
        - Alias of ``GGA_X_REVPBE``.

   Any other name ``X`` taken from the internal library is looked up as the
   internal routine ``CP_X``; unknown names stop the program with "Could not
   set functional pointer; functional CP_X n/a". The range-separation
   parameters of ``MGGA_XC_M11_L``, ``HYB_MGGA_XC_M11`` (:math:`\mu=0.25`) and
   ``HYB_MGGA_XC_MN12_SX`` (:math:`\mu=0.11`) are hard-coded. For the libxc
   functionals ``HYB_GGA_XC_B3LYP5`` and ``HYB_GGA_XC_PBEH`` the exact-exchange
   fractions 0.20 and 0.25 are set by default.

   **Important note on Minnesota functionals:**

   | At the usual plane wave cutoff values, the Minnesota functionals require
     unusually dense integration grids. In order to avoid introducing noise, this
     should be achieved by increasing the value of the :kw:`DUAL`. For
     usage recommendations and more details: Bircher, Lopez-Tarifa and
     Rothlisberger\ [bircher_rothlisberger19]_, DOI:
     10.1021/acs.jctc.8b00897.
   |

   The Minnesota functionals of the internal library stop the program ("Invalid
   setup for Minnesota functional. Please set DUAL to 8-12") unless
   :kw:`DUAL` or the density cutoff (see :kw:`CUTOFF`) is given explicitly in
   the ``&SYSTEM`` section.

   | The following exchange and correlation functionals are available as separate
     keywords within the internal CP library, both in spin-restricted and
     spin-unrestricted form, and can be freely combined:
   |

   **LDA_X**

   **LDA_C** \_PZ, \_PW, \_OB_PW, \_VWN

   **GGA_X** \_B88, \_OPTX, \_PBE, \_PBE_SOL, \_REVPBE, \_PBE_FLEX

   **GGA_C** \_LYP, \_P86, \_PBE, \_PBE_SOL, \_PBE_FLEX

   **MGGA_X** \_TPSS

   | **MGGA_C** \_TPSS
   |

   The parameters for VWN used within the CP library correspond to what is
   sometimes denoted ‘VWN5’. The version of HCTH available in prior releases
   corresponds to GGA_XC_HCTH_120. The definitions for the LDA part of the P86,
   PBE and TPSS correlation have been updated to comply with the standard
   definitions of most codes, for certain functionals, the definitions used in
   CPMD versions prior to 4.3 may be accessed by using
   :kw:`OLD_DEFINITIONS <OLD DEFINITIONS>`. A custom value of the parameter
   :math:`\beta` used in GGA_X_B88 can be input using the keyword :kw:`BECKE BETA`.

   | A completely customised parameterisation of PBE can be set up by using the
     \_PBE_FLEX routines in combination with the PBE_FLEX\_…keywords in the &DFT
     section. Since many different but highly specific parametrisations of PBE
     exist, this allows access to virtually all declinations based on the PBE xc
     equations.
   |

   **NEWCODE** and **OLDCODE**:

   Single keyword for setting up XC-functionals. Only *one single* functional
   keyword can be specified when using :kw:`OLDCODE` or
   :kw:`NEWCODE`.

   Available functionals are NONE, SONLY, LDA (in PADE form), BONLY, BP, BLYP,
   XLYP, GGA (=PW91), PBE, PBES, REVPBE, HCTH, OPTX, OLYP, TPSS, PBE0, B1LYP,
   B3LYP, X3LYP, HSE06.

   The following additional names are accepted. Except for the two aliases
   and the model potentials ``SAOP``, ``LB94`` and ``GLLB``, they all switch to
   :kw:`OLDCODE` internally.

   .. list-table:: Additional functional names of the legacy drivers
      :header-rows: 1
      :widths: 15 85

      * - Name
        - Definition in the code
      * - ``BP86``
        - Alias of ``BP``.
      * - ``PW91``
        - Alias of ``GGA``.
      * - ``PBE1W``
        - Slater exchange + Perdew-Zunger LDA correlation + PBE
          gradient-corrected exchange + 0.74 times the PBE gradient correction
          for correlation.
      * - ``REVPBE0``
        - Hybrid: 0.75 Slater + 0.75 revPBE gradient-corrected exchange + VWN
          LDA correlation + PBE gradient-corrected correlation + 25 % exact
          exchange.
      * - ``ZPBE0``
        - Hybrid with Ashcroft-screened exact exchange (25 %) built on the ZPBE
          exchange and correlation routines of DFRepository. Accepted by the
          parser, but these routines are not part of the distributed source:
          the run stops with "you need the DFRepository.F".
      * - ``XPBE0``
        - As ``ZPBE0`` with the XPBE hybrid exchange; also stops with "you need
          the DFRepository.F".
      * - ``PBEoriginal``
        - Removed. The program stops with "PBEoriginal isnt supported anymore"
          (only matched with exactly this mixed-case spelling).
      * - ``OLDX3LYP``
        - X3LYP as defined in CPMD 3.17 and earlier: the LDA part of LYP is
          used instead of VWN (0.782 Slater, 1.0 X exchange (0.542 B88 + 0.167
          PW91), 1.0 LDA-LYP, 0.871 LYP gradient correction, 21.8 % exact
          exchange).
      * - ``MODX3LYP``
        - X3LYP with 100 % Perdew-Zunger LDA correlation instead of 19 % VWN
          (reference printed by the code: J. Phys. Chem. B 110, 3685 (2006)).
      * - ``OLDB3LYP``
        - B3LYP as defined in CPMD 3.17 and earlier: the LDA part of LYP is
          used instead of VWN (0.8 Slater, 0.72 B88, 1.0 LDA-LYP, 0.81 LYP
          gradient correction, 20 % exact exchange).
      * - ``MODB3LYP``
        - B3LYP with 100 % Perdew-Zunger LDA correlation instead of 19 % VWN
          (reference printed by the code: J. Phys. Chem. B 110, 3685 (2006)).
      * - ``OLDB1LYP``
        - B1LYP as defined in CPMD 3.17 and earlier (0.75 Slater, 0.75 B88,
          LYP, 25 % exact exchange). In the present code the parameters are
          identical to ``B1LYP``; only the printed message differs.
      * - ``MODB1LYP``
        - B1LYP with Perdew-Zunger LDA correlation instead of the LDA part of
          LYP (reference printed by the code: J. Phys. Chem. B 110, 3685
          (2006)).
      * - ``SAOP``
        - Model exchange-correlation potential for TDDFT (no energy functional
          is available; total energies are not valid). Output label "SAOP
          POTENTIAL (modified)".
      * - ``LB94``
        - Model potential for TDDFT, as ``SAOP`` ("LB94 POTENTIAL
          (modified)").
      * - ``GLLB``
        - Model potential for TDDFT, as ``SAOP`` ("GLLB POTENTIAL
          (modified)").

   With the legacy drivers, VWN and PW LDA correlation are not available for
   spin-polarized (LSD) runs: the program stops with "LSD_VWN N/A, you may
   want to switch to USE_XC_DRIVER" (resp. ``LSD_PW``). This affects e.g.
   ``B3LYP``, ``X3LYP``, ``PBE0``, ``REVPBE0`` and :kw:`LDA CORRELATION`
   ``VWN``/``PW``. The fractions of the chosen legacy functional can be
   changed with :kw:`PXLDA`, :kw:`PXGC`, :kw:`PCLDA`, :kw:`PCGC` and
   :kw:`PHFX`.

.. source: src/dftin_utils.mod.F90:935 (legacy names: set_xc_functional,
.. line 1386; xc-driver names: src/cp_xc_utils.mod.F90:586)

.. keyword:: GC-CUTOFF
   :section: &DFT

   On the next line the density cutoff for the calculation of the gradient
   correction has to be specified. The default value is :math:`10^{-8}`.
   Experience showed that for a small :kw:`CUTOFF` value (e.g. when
   using Vanderbilt pseudopotentials) a bigger values have to be used. A
   reasonable value for a 25 ryd cutoff calculation is :math:`5 \cdot 10^{-6}`.

   **Warning**: for the HCTH functional, since it includes both the :math:`xc`
   part and the gradient correction in a unique functional, a GC-CUTOFF too high
   (e.g. :math:`\geq 5 \cdot 10^{-5}`) could result in not including any
   :math:`xc` part with uncontrolled related consequences.

.. keyword:: GRADIENT CORRECTION
   :section: &DFT
   :notes: [functionals]

   Individual components of gradient corrected functionals can be selected. Rarely
   needed anymore, use the :kw:`FUNCTIONAL` keyword instead. Only
   available for **OLDCODE** and **NEWCODE**; together with the xc driver
   (:kw:`XC_DRIVER <XC DRIVER>`, :kw:`CP_LIBRARY_ONLY`, :kw:`LIBXC_ONLY`) the
   program stops. The new xc driver offers similar flexibility in combination
   with the **FUNCTIONAL** keyword.

   Functionals implemented are for the exchange energy:

   **BECKE88** [Becke88]_, **GGAX** [Perdew92]_ (alias **PW91X**),
   **PBEX** [Perdew96]_, **REVPBEX** [Zhang98]_,
   **HCTH** [Handy98]_,
   **OPTX** [Optx]_,\ **PBESX** [Perdew07]_

   and for the correlation part:

   **PERDEW86** [Perdew86]_, **LYP** [Lee88]_,
   **GGAC** [Perdew92]_ (alias **PW91C**), **PBEC** [Perdew96]_,
   **REVPBEC** [Zhang98]_, **HCTH** [Handy98]_,
   **PBESC** [Perdew07]_.

   Note that for HCTH, exchange and correlation are treated as a unique
   functional. ``HCTH`` and ``OPTX`` switch to :kw:`OLDCODE`. ``REVPBEC`` is
   accepted as an alias of ``PBEC`` and selects exactly the same PBE
   correlation routine. The exchange part of the OLYP functional is ``OPTX``
   (OLYP itself is a :kw:`FUNCTIONAL` name).

   The words ``EXCH`` and ``CORREL`` select only the default exchange
   (BECKE88) or correlation (PERDEW86) gradient correction, respectively; the
   full words ``EXCHANGE`` and ``CORRELATION`` are not recognized. If no
   functionals are specified, the default functionals for both exchange and
   correlation are used.

   Example::

      GRADIENT CORRECTION BECKE88 LYP

.. source: src/dftin_utils.mod.F90:526

.. keyword:: HARTREE
   :section: &DFT

   Do a Hartree calculation. Only of use for testing purposes.

.. keyword:: HARTREE-FOCK
   :section: &DFT

   Do a Hartree–Fock calculation. This only works correctly for isolated systems.
   It should be used with care, it needs enormous amounts of CPU time. The
   spelling ``HARTREE FOCK`` (two words) is also accepted.

.. keyword:: HFX SCREENING
   :section: &DFT
   :alternatives: [RADIUS, DIAG, KEEP_RADIUS_FIXED, EPS_INT, RECOMPUTE_TWO_INT_LIST_EVERY]

   Perform the calculation of exact exchange using Wannier functions with
   screening of orbital pairs. Exactly one option is evaluated per line; give
   the keyword several times to combine options. Without a recognized option
   the program stops ("need to specify something...").

   RADIUS
      Orbital pairs are screened according to the distance of the Wannier
      centres; the radius is read from the next line.

   KEEP_RADIUS_FIXED
      Keep the screening radius fixed. No value is read.

   EPS_INT
      Orbital pairs are screened according to the value of the integrals; the
      threshold is read from the next line.

   RECOMPUTE_TWO_INT_LIST_EVERY
      The number of steps after which the integral list is recomputed is read
      from the next line (by default the list is never recomputed).

   DIAG
      Only the diagonal terms are included. No value is read; this disables
      ``RADIUS`` and ``EPS_INT``.

   Any ``HFX SCREENING`` line selects the old HFX module unless
   :kw:`FORCE_NEW_HFX` is given. The keyword is ignored if the functional
   contains no exact exchange. The alias ``HFX-SCREENING`` is also accepted.

   Example::

      HFX SCREENING EPS_INT
       1.0E-6
      HFX SCREENING RECOMPUTE_TWO_INT_LIST_EVERY
       10

.. source: src/dftin_utils.mod.F90:990

.. keyword:: HFX_BLOCK_SIZE
   :target: HFX BLOCK SIZE
   :section: &DFT

   Set the block size for the HFX calculations. This keyword is active only when
   no thresholding of the integrals is used. Default is 2. The block size is read
   from the next line.

.. keyword:: HFX_DISTRIBUTION
   :target: HFX DISTRIBUTION
   :section: &DFT

   Set the block distribution for the parallel HFX calculations. This keyword is
   active only when no thresholding of the integrals is used. The distribution
   (``BLOCK_CYCLIC`` or ``DYNAMIC``) is given on the same line, e.g.
   ``HFX_DISTRIBUTION DYNAMIC``. **Default** is ``BLOCK_CYCLIC``.

.. source: src/dftin_utils.mod.F90:663

.. keyword:: HFX_NEW_SCDM
   :section: &DFT

   Select the modified SCDM ("selected column density matrix") localization of
   the orbitals used in the :kw:`HFX_SCDM` exact-exchange scheme: the columns are
   pre-selected from the grid points where the electron density is larger and its
   gradient smaller than their averages over the grid (routine ``loc_scdm_new``,
   described in the source as a modified version of the SCDM implementation of
   Quantum ESPRESSO). Without this keyword the original routine
   ``LOCALIZE_SCDM`` is used. The keyword takes no argument and only has an
   effect together with :kw:`HFX_SCDM`; the SCDM localization itself is only
   performed when ACE is active (``USE_ACE`` in the ``&MTS`` section).

   **Default** is the original SCDM localization.

   .. todo::

      MISSING: literature reference and recommended use of the modified SCDM
      variant; whether HFX_NEW_SCDM is meaningful without USE_ACE.

.. source: src/dftin_utils.mod.F90:682 (used in src/hfx_utils.mod.F90:2146)

.. keyword:: HFX_SCALE
   :target: HFX SCALE
   :section: &DFT
   :notes: scaling_for_exact_exchange

   Only available in combination with :kw:`XC_DRIVER <XC DRIVER>`. On the same
   line, a scaling value (:math:`\in [0,1]`) for the exact exchange contribution
   to the overall xc functional is specified. If used with a built-in hybrid
   functional, the default value is changed; if used with any other functional,
   the keyword activates the inclusion of exact exchange. The alias
   ``HFX_SCALES`` is also accepted; :kw:`PHFX` has the same effect.

   Do not use this keyword with Coulomb-attenuated functionals
   (:kw:`COULOMB ATTENUATION`, ``SCREENED EXCHANGE CAM``,
   ``HYB_GGA_XC_CAM_B3LYP``): for the CAM operator the value given here
   replaces the parameter :math:`\alpha` (whichever of ``HFX_SCALE`` and the
   CAM parameter line is read last is used), and with ``HYB_GGA_XC_CAM_B3LYP``
   it disables the built-in attenuation. Change :math:`\alpha` and
   :math:`\beta` via the CAM parameter line instead.

.. source: src/dftin_utils.mod.F90:818

.. keyword:: HFX_SCDM
   :section: &DFT

   Compute the exact exchange with the SCDM ("selected column density matrix")
   scheme (output banner: "HFX WITH SCDM APPROACH", S. Mandal, IIT Kanpur).
   The pair cutoff ``SCDM_CUTOFF`` is read from the next line. An orbital pair
   :math:`(i,j)` is neglected in the exchange sum if
   :math:`\frac{1}{N}\sum_{\mathbf{r}}|\psi_i(\mathbf{r})\psi_j(\mathbf{r})| <`
   ``SCDM_CUTOFF``, where :math:`N` is the number of real-space grid points.
   When ACE is active (``USE_ACE`` in the ``&MTS`` section) the orbitals are first
   localized with the SCDM method (see :kw:`HFX_NEW_SCDM`). The scheme is only
   used with the new HFX module (i.e. not after :kw:`HFX SCREENING` unless
   :kw:`FORCE_NEW_HFX` is given) and only if exact exchange is part of the
   functional.

   **Default** for ``SCDM_CUTOFF`` is :math:`10^{-8}`; the scheme is off by default.

   Example::

      HFX_SCDM
       1.0E-6

   .. todo::

      MISSING: literature reference for the SCDM-based HFX implementation and
      recommended values of SCDM_CUTOFF. The commented-out parameters
      ``de_cutoff`` and ``n_loop`` (dftin_utils.mod.F90:674) are not read.

.. source: src/dftin_utils.mod.F90:673 (used in src/hfx_drivers.mod.F90:124-138,
.. src/hfx_utils.mod.F90:3134 SELECT_PAIR)

.. keyword:: HUBBARD
   :section: &DFT
   :options: [NORM,ORTHO,NUATM=nuatm,OCCMAT=printfreq,VERB]

   Use the **HUBBARD-U** correction to partly correct the self interaction error
   of DFT. Specifying NORM and/or ORTHO causes the orbitals for the projections to
   be normalized and/or orthogonalized, respectively. VERB causes the occupation
   matrix to be printed in the output, otherwise it is printed to the file OCCMAT.
   The keyword OCCMAT=\ :math:`printfreq` specifies the frequency in which the
   OCCMAT is written. Default is 1. :math:`nuatm` is the number of atoms which are
   treated with the Hubbard U correction (U-atom). Default is 1.

   For each U-atom the atom number, the U parameter [eV], the linear response
   parameter :math:`\alpha` [eV] and the number of different angular momenta for
   projectors is read from a new line. For each of those, the shell number of the
   projector and angular momentum :math:`l` is read from the next lines.

   The options are matched as substrings of the keyword line. If ``NUATM=`` is
   missing, a warning is printed and :math:`nuatm=1` is used. The keyword
   cannot be combined with broken-symmetry calculations (``BSYMM``); the
   program stops.

.. source: src/dftin_utils.mod.F90:1028

.. keyword:: KERNEL_HFX_SCALE
   :target: KERNEL HFX SCALE
   :section: &DFT
   :notes: scaling_for_exact_exchange_in_kernel

   Identical functionality as described in :kw:`HFX_SCALE <HFX SCALE>`, but
   applies to the :kw:`LR KERNEL` rather than the
   :kw:`FUNCTIONAL`. The alias ``KERNEL_HFX_SCALES`` is also accepted.

.. keyword:: KERNEL_LIBRARY
   :target: KERNEL LIBRARY
   :section: &DFT
   :notes: library_for_kernel_1, library_for_kernel_2, ...

   Same functionality as :kw:`LIBRARY`, but applies to the functionals of
   :kw:`LR KERNEL` rather than to those of :kw:`FUNCTIONAL`. The number of
   entries must equal the number of kernel functionals. **Default** is the
   internal library; :kw:`CP_LIBRARY_ONLY` and :kw:`LIBXC_ONLY` do not change
   this default, but an entry contradicting them stops the program if it is
   given after them.

.. source: src/dftin_utils.mod.F90:390

.. keyword:: KERNEL_SCALES
   :target: KERNEL SCALES
   :section: &DFT
   :notes: scaling_for_kernel_1, scaling_for_kernel_2, ...

   Identical functionality as described in :kw:`SCALES`, but applies to
   the :kw:`LR KERNEL` rather than the
   :kw:`FUNCTIONAL`. The alias ``KERNEL_SCALE`` is also accepted; the line
   must not contain the word ``HFX``.

.. keyword:: LDA CORRELATION
   :section: &DFT
   :notes: [functional]

   The LDA correlation functional is specified. Only available for **OLDCODE** and
   **NEWCODE**; together with the xc driver the program stops. The keyword is
   only effective together with :kw:`GRADIENT CORRECTION` (or alone); it is
   overwritten by :kw:`FUNCTIONAL`.

   Possible functionals are **NO** (no correlation functional),
   **PZ** [Perdew81]_, the Perdew and Zunger fit to the data of Ceperley and
   Alder [Ceperley80]_, **VWN** [Vosko80]_,
   **LYP** [Lee88]_ and **PW** [Perdew91]_. Further options are:

   .. list-table:: Additional LDA CORRELATION options
      :header-rows: 1
      :widths: 15 85

      * - Option
        - Meaning
      * - ``WIGNER``
        - Wigner correlation.
      * - ``HEDIN``
        - Hedin-Lundqvist correlation.
      * - ``OBPZ``
        - Ortiz-Ballone parametrization (PRB 50, 1391 (1994)), Perdew-Zunger
          form.
      * - ``OBPW``
        - Ortiz-Ballone parametrization (PRB 50, 1391 (1994)), Perdew-Wang
          form.
      * - ``TETER``
        - Pade approximation of Goedecker, Teter and Hutter for exchange and
          correlation together (Slater exchange is switched off).
      * - ``HCTH``
        - HCTH exchange-correlation functional (switches to :kw:`OLDCODE`;
          LDA exchange and correlation are switched off and the HCTH gradient
          correction is used).

   If the keyword is not given, the Pade form is used with :kw:`NEWCODE` and
   PZ with :kw:`OLDCODE`.

.. source: src/dftin_utils.mod.F90:614 (defaults: line 2086)

.. keyword:: LIBRARY
   :section: &DFT
   :notes: library_for_functional_1, library_for_functional_2, ...

   In combination with :kw:`XC_DRIVER <XC DRIVER>`, both functionals from a new,
   internal CP library and from libxc are available. On the same line, for every
   functional specified in **FUNCTIONAL**, the library has to be indicated by
   **INTERNAL** or **CP** for the internal CP library, or by **LIBXC** (or
   **LXC**) for libxc. The number of entries must equal the number of
   functionals. If the keyword is omitted, the functionals will default to the
   internal CP library.

   :kw:`CP_LIBRARY_ONLY` and :kw:`LIBXC_ONLY` set the library of all
   functionals and override a ``LIBRARY`` line given before them; a
   ``LIBRARY`` line given after them must not contain an entry of the other
   library (the program stops).

   The use of LIBXC is usually associated to a considerable overhead with respect
   to the internal implementation. It is therefore discouraged for functionals
   that are also available in the internal library, *cf.*
   :kw:`FUNCTIONAL`.

.. keyword:: LIBXC_ONLY
   :target: LIBXC ONLY
   :section: &DFT

   Use the new xc driver and take all functionals of :kw:`FUNCTIONAL` from
   libxc (the internal CP library is not used). CPMD must have been compiled
   with libxc; otherwise the program stops with "no libxc available".

   The keyword switches on the xc driver by itself, i.e. it implies
   :kw:`XC_DRIVER <XC DRIVER>`; a separate ``XC_DRIVER`` line is not needed.
   Consequently:

   * the functionals must be given with :kw:`FUNCTIONAL` using libxc names;
     the legacy names (``PBE``, ``BLYP``, ...) are not valid;
   * :kw:`GRADIENT CORRECTION`, :kw:`LDA CORRELATION`, :kw:`SLATER` (except
     ``SLATER ALPHA``), :kw:`ACM0`, :kw:`ACM1` and :kw:`ACM3` stop the program;
   * :kw:`OLDCODE` and :kw:`NEWCODE` have no effect.

   The keyword sets the library of the functionals of :kw:`FUNCTIONAL` only,
   overriding a :kw:`LIBRARY` line given before it. If :kw:`LR KERNEL` is not
   given, the kernel is a copy of the functional and therefore also uses
   libxc. Functionals given with :kw:`LR KERNEL` or :kw:`MTS_LOW_FUNC`,
   however, still default to the internal library and need
   :kw:`KERNEL_LIBRARY` or :kw:`MTS_LOW_FUNC_LIBRARY` with ``LIBXC`` entries.
   A :kw:`LIBRARY`, :kw:`KERNEL_LIBRARY` or :kw:`MTS_LOW_FUNC_LIBRARY` entry
   ``CP`` (or ``INTERNAL``) given *after* this keyword stops the program. If
   both :kw:`CP_LIBRARY_ONLY` and :kw:`LIBXC_ONLY` are given, the last one is
   used.

   The alias ``USE_LIBXC`` is also accepted.

   Example::

      LIBXC_ONLY
      FUNCTIONAL GGA_X_B88 GGA_C_LYP

.. source: src/dftin_utils.mod.F90:334

.. keyword:: LR KERNEL
   :section: &DFT
   :notes: functional

   Use another functional for the linear-response kernel. With
   :kw:`XC_DRIVER <XC DRIVER>` the functionals are given as for
   :kw:`FUNCTIONAL`; with the legacy drivers only ``NONE``, ``LDA``, ``PBE``
   and ``PBE0`` are available (other names stop the program with "TDDFT kernel
   functional not available"). If the keyword is not given, the kernel is the
   same as the functional. The aliases ``LR_KERNEL`` and ``XC_KERNEL`` are also
   accepted.

.. source: src/dftin_utils.mod.F90:984 (kernel names: line 1989)

.. keyword:: MTS_HIGH_FUNC
   :target: MTS HIGH FUNC
   :section: &DFT
   :notes: functionals

   Alias for :kw:`FUNCTIONAL` that makes sense when used with the
   MTS scheme since it corresponds to the high level functional in that context
   (see &MTS section). **Note:** The functionals in combination with the MTS
   scheme have to be set with :kw:`XC_DRIVER <XC DRIVER>`.

.. keyword:: MTS_LOW_FUNC
   :target: MTS LOW FUNC
   :section: &DFT
   :notes: functionals

   Select the low level functional in the MTS scheme. To be used like
   :kw:`FUNCTIONAL`. **Note:** The functionals in combination with
   the MTS scheme have to be set with :kw:`XC_DRIVER <XC DRIVER>`.

.. keyword:: MTS_LOW_FUNC_HFX_SCALE
   :section: &DFT
   :notes: scaling_for_exact_exchange

   Same as :kw:`HFX_SCALE`, but applies to the low-level functional of the
   multiple-time-step (MTS) scheme, :kw:`MTS_LOW_FUNC`, rather than to
   :kw:`FUNCTIONAL`. The value is read from the same line. Giving the keyword
   switches on exact exchange for the low-level functional and disables the
   built-in hybrid defaults of the functionals listed in :kw:`MTS_LOW_FUNC`.
   Only meaningful with :kw:`XC_DRIVER`. The alias ``MTS_LOW_FUNC_HFX_SCALES`` is
   also accepted.

   Example::

      MTS_LOW_FUNC          GGA_X_PBE GGA_C_PBE
      MTS_LOW_FUNC_SCALES   0.75      1.00
      MTS_LOW_FUNC_HFX_SCALE 0.25

.. source: src/dftin_utils.mod.F90:840

.. keyword:: MTS_LOW_FUNC_LIBRARY
   :section: &DFT
   :notes: library_for_functional_1, library_for_functional_2, ...

   Same as :kw:`LIBRARY`, but applies to the functionals given with
   :kw:`MTS_LOW_FUNC`. On the same line, one library tag per functional is
   expected, in the same order: ``CP`` or ``INTERNAL`` for the internal library,
   ``LIBXC`` or ``LXC`` for libxc. The number of tags must equal the number of
   functionals. Choosing ``LIBXC`` after :kw:`CP_LIBRARY_ONLY`, or ``CP`` after
   :kw:`LIBXC_ONLY`, stops the program.

   **Default** is the internal library for all low-level functionals; note that
   :kw:`LIBXC_ONLY` does not change this default.

   Example::

      MTS_LOW_FUNC          LDA_X LDA_C_VWN_3
      MTS_LOW_FUNC_LIBRARY  CP    LIBXC

.. source: src/dftin_utils.mod.F90:418

.. keyword:: MTS_LOW_FUNC_SCALES
   :section: &DFT
   :notes: scaling_for_functional_1, scaling_for_functional_2, ...

   Same as :kw:`SCALES`, but applies to the functionals given with
   :kw:`MTS_LOW_FUNC`. The scaling factors are read from the same line in the
   order of the functionals; their number must equal the number of functionals.
   Only meaningful with :kw:`XC_DRIVER`. The alias ``MTS_LOW_FUNC_SCALE`` is also
   accepted; the line must not contain the word ``HFX``.

   **Default** is 1.0 for every functional.

.. source: src/dftin_utils.mod.F90:483

.. keyword:: NEWCODE
   :section: &DFT

   | This keyword will be deprecated in a future release. Both **NEWCODE** and
     **OLDCODE** have been replaced by the new :kw:`XC_DRIVER <XC DRIVER>`.
   |

   The following description applies up to CPMD version prior to 4.1:

   **NEWCODE** is a switch to select one out of two versions of code to calculate
   exchange-correlation functionals.

   NEWCODE is the default, but not all functionals are available with NEWCODE, if
   you select one of these, :kw:`OLDCODE` is used automatically.
   NEWCODE is highly recommended for all new projects and especially for vector
   computers, also some of the newer functionality is untested or non-functional
   with OLDCODE.

.. keyword:: NUMERICAL_DIV
   :target: NUMERICAL DIV
   :section: &DFT

   Calculate the :math:`\mathbf{G}=0` term for the Coulomb-attenuated exact
   exchange using a numerical description of the integral term. This is the
   default only when CAM-B3LYP is used with LIBXC. In all other situations,
   :kw:`ANALYTICAL_DIV <ANALYTICAL DIV>` is used.

.. keyword:: OLDCODE
   :section: &DFT

   This keyword will be deprecated in a future release. Please use the new
   :kw:`XC_DRIVER <XC DRIVER>`. *Cf.* :kw:`NEWCODE`.

.. keyword:: OLD_DEFINITIONS
   :target: OLD DEFINITIONS
   :section: &DFT

   For certain functionals, older versions of CPMD did not use standard
   definitions. When using the new :kw:`XC_DRIVER <XC DRIVER>`, those definitions
   can be invoked by using **OLD_DEFINITIONS**. The option affects the PBE
   correlation in ``GGA_C_PBE``, ``GGA_C_PBE_SOL``, ``GGA_XC_PBE``,
   ``GGA_XC_REVPBE`` (``GGA_XC_PBE_R``), ``GGA_XC_PBE_SOL`` and
   ``HYB_GGA_XC_PBE0``, and the TPSS correlation in ``MGGA_C_TPSS`` and
   ``MGGA_XC_TPSS``. ``HYB_GGA_XC_REVPBE0`` is not affected, and with
   ``HYB_GGA_XC_PBE0_SOL`` the program stops ("Old definition n/a; Pade xc not
   coded"). The alias ``COMPATIBILITY`` is also accepted.

.. source: src/dftin_utils.mod.F90:501 (src/cp_xc_utils.mod.F90:596-1360)

.. keyword:: PBE_FLEX_BETA
   :target: PBE FLEX BETA
   :section: &DFT

   In combination with the PBE_FLEX correlation functional, the parameter
   :math:`\beta` can be specified on the next line.

.. keyword:: PBE_FLEX_GAMMA
   :target: PBE FLEX GAMMA
   :section: &DFT

   In combination with the PBE_FLEX correlation functional, the parameter
   :math:`\gamma` can be specified on the next line.

.. keyword:: PBE_FLEX_KAPPA
   :target: PBE FLEX KAPPA
   :section: &DFT

   In combination with the PBE_FLEX exchange functional (``GGA_X_PBE_FLEX``),
   the parameter :math:`\kappa` is read from the next line.

.. keyword:: PBE_FLEX_MU
   :target: PBE FLEX MU
   :section: &DFT

   In combination with the PBE_FLEX exchange functional (``GGA_X_PBE_FLEX``),
   the parameter :math:`\mu` is read from the next line.

.. keyword:: PBE_FLEX_UEG_CORRELATION
   :target: PBE FLEX UEG CORRELATION
   :section: &DFT
   :notes: functional

   In combination with the PBE_FLEX correlation functional, the LDA correlation
   functional can be specified on the same line. Available options are
   ``LDA_C_PW`` (the default), ``LDA_C_PZ``, ``LDA_C_VWN`` and
   ``LDA_C_OB_PW``. Any other value silently falls back to ``LDA_C_PW``.

.. source: src/dftin_utils.mod.F90:521 (src/cp_gga_correlation_utils.mod.F90:1127)

.. keyword:: PCGC
   :section: &DFT

   Overwrite the scaling factor of the gradient-correction part of the
   correlation functional (internal variable ``func3%pcgc``). The new value is
   read from the next line. Only used by the legacy drivers
   (:kw:`OLDCODE`/:kw:`NEWCODE`) and only when the functional is set with
   :kw:`FUNCTIONAL`; the value replaces the one implied by the chosen functional
   (e.g. 0.81 for B3LYP). It has no effect with :kw:`XC_DRIVER` (use
   :kw:`SCALES` instead) or when the functional is assembled with
   :kw:`GRADIENT CORRECTION`. In the gradient-correction routines the factor is
   applied only if a hybrid functional is active; for hybrids containing LYP
   it also scales the LYP contribution in the LDA routines.

   **Default** is the value set by the functional (1.0 for non-hybrid
   functionals).

   Example::

      FUNCTIONAL B3LYP
      PCGC
       0.81

   .. todo::

      MISSING: intended use cases; confirm whether the factor is used by the
      NEWCODE LDA routines for non-hybrid functionals.

.. source: src/dftin_utils.mod.F90:872 (applied at 1865; used in
.. src/gcener_utils.mod.F90:677, src/functionals_utils.mod.F90:219)

.. keyword:: PCLDA
   :section: &DFT

   Overwrite the scaling factor of the LDA correlation functional (internal
   variable ``func3%pclda``). The new value is read from the next line. Same
   restrictions as :kw:`PCGC`: only used with :kw:`FUNCTIONAL` in the legacy
   drivers, ignored with :kw:`XC_DRIVER`.

   **Default** is the value set by the functional (e.g. 0.19 for B3LYP and
   X3LYP, 1.0 otherwise).

   .. todo::

      MISSING: confirm whether the factor is used by the NEWCODE LDA routines
      for non-hybrid functionals.

.. source: src/dftin_utils.mod.F90:888 (applied at 1867; used in
.. src/functionals_utils.mod.F90:155, src/lsd_func_utils.mod.F90:114)

.. keyword:: PHFX
   :section: &DFT

   Set the fraction of exact (Hartree-Fock) exchange. The value is read from
   the next line and must lie between 0 (no exact exchange) and 1 (100 % exact
   exchange).

   With the legacy drivers (:kw:`OLDCODE`/:kw:`NEWCODE`) the value replaces the
   fraction built into the hybrid functional selected with :kw:`FUNCTIONAL`
   (e.g. 0.25 for PBE0 and HSE06, 0.20 for B3LYP). The fractions of the
   semi-local exchange are not adjusted (see :kw:`PXLDA` and :kw:`PXGC`), and
   the keyword has no effect for non-hybrid functionals (exact exchange is not
   switched on) or when the functional is assembled with
   :kw:`GRADIENT CORRECTION`. For HSE06 this keyword is the only way to change
   the default fraction of 0.25 (:kw:`ACM1` is overwritten by the functional
   definition).

   With :kw:`XC_DRIVER <XC DRIVER>`, PHFX has the same effect as
   :kw:`HFX_SCALE` (which reads the value from the same line); the use of
   :kw:`HFX_SCALE` is recommended there.

   Example::

      FUNCTIONAL HSE06
      PHFX
       0.30

.. source: src/dftin_utils.mod.F90:851 (applied at line 1874)

.. keyword:: PXGC
   :section: &DFT

   Overwrite the scaling factor of the gradient-correction part of the
   exchange functional (internal variable ``func3%pxgc``). The new value is read
   from the next line. Same restrictions as :kw:`PCGC`: only used with
   :kw:`FUNCTIONAL` in the legacy drivers, ignored with :kw:`XC_DRIVER`; in the
   gradient-correction routines the factor is applied only if a hybrid
   functional is active.

   **Default** is the value set by the functional (e.g. 0.75 for PBE0, 0.72 for
   B3LYP, 1.0 for non-hybrid functionals).

   Example::

      FUNCTIONAL PBE0
      PHFX
       0.30
      PXGC
       0.70
      PXLDA
       0.70

.. source: src/dftin_utils.mod.F90:864 (applied at 1864; used in
.. src/gcener_utils.mod.F90:667-676)

.. keyword:: PXLDA
   :section: &DFT

   Overwrite the scaling factor of the LDA (Slater) exchange functional
   (internal variable ``func3%pxlda``). The new value is read from the next
   line. Same restrictions as :kw:`PCGC`: only used with :kw:`FUNCTIONAL` in the
   legacy drivers, ignored with :kw:`XC_DRIVER`.

   **Default** is the value set by the functional (e.g. 0.75 for PBE0, 0.8 for
   B3LYP, 1.0 for non-hybrid functionals).

   .. todo::

      MISSING: confirm whether the factor is used by the NEWCODE LDA routines
      for non-hybrid functionals.

.. source: src/dftin_utils.mod.F90:880 (applied at 1866; used in
.. src/functionals_utils.mod.F90:92, src/lsd_func_utils.mod.F90:73)

.. keyword:: REFUNCT
   :section: &DFT
   :notes: functionals

   Use a special reference functional in a calculation. This option is not
   active: the functional is parsed, but the values are not used anywhere in
   the code. Accepted names are ``LDA``, ``PBE``, ``BLYP``, ``OLYP`` and
   ``BP`` (alias ``BP86``); other names stop the program ("Reference functional
   not available").

.. source: src/dftin_utils.mod.F90:943

.. keyword:: SCALED EXCHANGE
   :section: &DFT

   Alias for :kw:`SCEX`.

.. keyword:: SCALES
   :section: &DFT
   :notes: scaling_for_functional_1, scaling_for_functional_2, ...

   Only available in combination with :kw:`XC_DRIVER <XC DRIVER>`. On the same
   line, rescaling values (:math:`\in [0,1]`) for the different
   exchange-correlation functionals specified in :kw:`FUNCTIONAL`
   are read in the same order in which the functionals were specified. The
   number of values must equal the number of functionals. The alias ``SCALE``
   is also accepted; the line must not contain the word ``HFX``.

.. keyword:: SCEX
   :section: &DFT

   Activate the use of the coordinate-scaled exact exchange scheme by Bircher and
   Rothlisberger\ [bircher_rothlisberger18b]_ for *isolated systems*.
   Speedups of up to one order of magnitude can be achieved. Please note that the
   simulation supercell must imperatively span twice the size of the charge
   density, and that the molecule has to remain centred in the box.

.. keyword:: SCREENED EXCHANGE
   :section: &DFT
   :options: {ASHCROFT,CAM,ERFC,EXP}

   Screening / range separation applied to the Hartree-Fock exchange term *only*.
   Without option, ``EXP`` is used. Possible options are:

   ASHCROFT
      Ashcroft exchange. :math:`r_\text{cut}` is read from the next line.

   CAM
      Coulomb-attenuation method. The parameters are read from the next line
      in the order :math:`\mu`, :math:`\alpha`, :math:`\beta` (as for
      :kw:`COULOMB ATTENUATION`). Please note that in case attenuation of the
      (semi-)local part is desired, too, the keyword
      :kw:`COULOMB ATTENUATION` has to be used. ``LC`` is accepted as an alias
      of ``CAM``.

   ERFC
      Screened exchange using the complementary error function. The switching
      parameter :math:`\gamma` is read from the next line.

   EXP
      Screened exchange using an exponential function. The switching parameter
      :math:`\gamma` is read from the next line.

   With :kw:`XC_DRIVER <XC DRIVER>`, this keyword replaces the built-in
   exact-exchange settings of hybrid functionals: for ``EXP``, ``ERFC`` and
   ``ASHCROFT`` the fraction of exact exchange must then be given with
   :kw:`HFX_SCALE` (otherwise it is zero); for ``CAM`` the fractions are the
   :math:`\alpha` and :math:`\beta` read from the input. The settings also
   apply to the functionals of :kw:`LR KERNEL` and :kw:`MTS_LOW_FUNC`. With the
   legacy drivers, the settings are overwritten by ``FUNCTIONAL HSE06``
   (:math:`\omega` is fixed to 0.106).

.. source: src/dftin_utils.mod.F90:718

.. keyword:: SLATER
   :section: &DFT
   :alternatives: [NO]

   The :math:`\alpha` value for the Slater exchange
   functional [Slater51]_ is read from the next line and Slater exchange is
   selected. With NO the exchange functional is switched off. With ``ALPHA``
   (``SLATER ALPHA``) only the value of :math:`\alpha` is read from the next
   line, without selecting Slater exchange; this is the only form that may be
   used together with the xc driver.

   Default is a value of 2/3.

   This option together with no correlation functional, allows for :math:`X\alpha`
   theory.

.. keyword:: SMOOTH
   :section: &DFT

   A smoothening function is applied to the density [Laasonen93]_. The alias
   ``SMOOTHING`` is also accepted.

   The function is of the Fermi type.

   .. math::

      f(G) = \frac{1}{
          \displaystyle{1 + e^{\frac{\scriptstyle{G - G_{\scriptstyle cut}}}
                                    {\scriptstyle\Delta}}}}

   G is the wavevector, :math:`G_{cut} = \alpha\,G_{max}` and
   :math:`\Delta = \beta\,G_{max}`. Values for :math:`\alpha` and :math:`\beta`
   have to be given on the next line.

.. keyword:: XC_DRIVER
   :target: XC DRIVER
   :section: &DFT

   Use the new CP xc driver. This keyword replaces both :kw:`OLDCODE`
   and :kw:`NEWCODE` and offers more ample customisation options (*cf.*
   :kw:`COULOMB ATTENUATION`, :kw:`LIBRARY` and
   :kw:`SCALES`). All functionals are now consistently available for
   both spin-restricted and open-shell systems. Additionally, functionals from
   libxc can be used and freely mixed with the internal functionals, *cf.*
   :kw:`LIBRARY`. Please note that the new
   :kw:`FUNCTIONAL` abbreviations differ from **OLDCODE** and
   **NEWCODE**, *cf.* the corresponding entry for a list of available functionals.
   The new driver uses the standard definitions for all functionals; where those
   definitions differ from older versions of CPMD, the old definitions (where
   implemented) can be accessed by using the keyword
   :kw:`OLD_DEFINITIONS <OLD DEFINITIONS>`.

   The aliases ``USE_XC_DRIVER`` and ``USE_DRIVER`` are also accepted. The xc
   driver is also switched on implicitly by :kw:`CP_LIBRARY_ONLY`,
   :kw:`LIBXC_ONLY` and :kw:`COULOMB ATTENUATION`.

.. source: src/dftin_utils.mod.F90:349

.. keyword:: XC_KERNEL
   :target: XC KERNEL
   :section: &DFT
   :notes: functional

   Alias for :kw:`LR KERNEL`.
