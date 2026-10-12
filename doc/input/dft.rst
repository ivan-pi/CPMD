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

   Add exact exchange to the specified :kw:`FUNCTIONAL` according to
   the adiabatic connection method 0. [acm0]_ [adamo2000]_

.. keyword:: ACM1
   :section: &DFT

   Add exact exchange to the specified :kw:`FUNCTIONAL` according to
   the adiabatic connection method 1. [adamo2000]_ [acm1]_ The parameter
   is read from the next line.

.. keyword:: ACM3
   :section: &DFT

   Add exact exchange to the specified :kw:`FUNCTIONAL` according to
   the adiabatic connection method 3. [adamo2000]_ [acm3]_ The three
   needed parameters are read from the next line.

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

.. keyword:: COULOMB ATTENUATION
   :section: &DFT

   Activates the long-range correction/Coulomb attenuation method (**LC, CAM**).
   With this keyword, all selected *internal* LDA or GGA exchange functionals will
   be Coulomb-attenuated, and the corresponding exact exchange contribution is
   added by default. The three CAM parameters :math:`\alpha`, :math:`\beta` and
   :math:`\mu` are read from the next line. Setting :math:`\alpha=0` and
   :math:`\beta=1.0` is equivalent to using the long-range correction (LC) method.
   This option is only available in combination with the new
   :kw:`XC_DRIVER <XC DRIVER>`. The keyword is aliased and may be abbreviated by
   **CAM**.

.. keyword:: CP_LIBRARY_ONLY
   :target: CP LIBRARY ONLY
   :section: &DFT

   Use the new xc driver, but use only functionals that are available internally
   (this excludes libxc). *Cf.* :kw:`XC_DRIVER <XC DRIVER>`.

.. keyword:: EXCHANGE CORRELATION TABLE
   :section: &DFT
   :alternatives: [NO]

   Specifies the range and the granularity of the lookup table for the local
   exchange-correlation energy and potential. The number of table entries and the
   maximum density have to be given on the next line.

   Note that this keyword is only relevant when using :kw:`OLDCODE` and
   even then it is set to **NO** be default. Previous default values were 30000
   and 2.0.

.. keyword:: FUNCTIONAL
   :section: &DFT
   :notes: functional_1, functional_2, ...

   | The xc functional(s) are specified on the same line, separated by a space.
   |

   **XC_DRIVER**:

   Both functionals from the internal CP library and those linked from libxc are
   available within the same run and can be freely mixed when using the
   :kw:`XC_DRIVER <XC DRIVER>`. Please note that the internal CP functionals
   offer a considerable performance advantage over libxc. While the code will
   default to the internal implementations, use of libxc for single functionals
   can be requested by the keyword :kw:`LIBRARY`. The functional code
   names for the functionals available internally follow the nomenclature
   introduced by libxc: The functional name must be preceeded by a tag specifying
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

   **Important note on Minnesota functionals:**

   | At the usual plane wave cutoff values, the Minnesota functionals require
     unusually dense integration grids. In order to avoid introducing noise, this
     should be achieved by increasing the value of the :kw:`DUAL`. For
     usage recommendations and more details: Bircher, Lopez-Tarifa and
     Rothlisberger\ [bircher_rothlisberger19]_, DOI:
     10.1021/acs.jctc.8b00897.
   |

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
   keyword can be specified in when using :kw:`OLDCODE` or
   :kw:`NEWCODE`.

   Available functionals are NONE, SONLY, LDA (in PADE form), BONLY, BP, BLYP,
   XLYP, GGA (=PW91), PBE, PBES, REVPBE, HCTH, OPTX, OLYP, TPSS, PBE0, B1LYP,
   B3LYP, X3LYP,PBES, HSE06

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
   available for **OLDCODE** and **NEWCODE**; but the new
   :kw:`XC_DRIVER <XC DRIVER>` offers similar flexibility in combination with the
   **FUNCTIONAL** keyword.

   Functionals implemented are for the exchange energy:

   **BECKE88** [Becke88]_, **GGAX** [Perdew92]_
   **PBEX** [Perdew96]_, **REVPBEX** [Zhang98]_,
   **HCTH** [Handy98]_,
   **OPTX** [Optx]_,\ **PBESX** [Perdew07]_

   and for the correlation part:

   **PERDEW86** [Perdew86]_, **LYP** [Lee88]_,
   **GGAC** [Perdew92]_, **PBEC** [Perdew96]_,
   **REVPBEC** [Zhang98]_, **HCTH** [Handy98]_
   **OLYP** [Optx]_,\ **PBESC** [Perdew07]_.

   Note that for HCTH, exchange and correlation are treated as a unique
   functional.

   The keywords **EXCHANGE** and **CORRELATION** can be used for the default
   functionals (currently BECKE88 and PERDEW86). If no functionals are specified
   the default functionals for exchange and correlation are used.

.. keyword:: HARTREE
   :section: &DFT

   Do a Hartree calculation. Only of use for testing purposes.

.. keyword:: HARTREE-FOCK
   :section: &DFT

   Do a Hartree–Fock calculation. This only works correctly for isolated systems.
   It should be used with care, it needs enormous amounts of CPU time.

.. keyword:: HFX SCREENING
   :section: &DFT
   :options: {WFC,DIAG,EPS_INT,RECOMPUTE_TWO_INT_LIST_EVERY}

   Read value from the next line.

   Perform the calculation of exact exchange using Wannier functions. Orbital
   pairs are screened according to the distance of the Wannier centers *WFC*, the
   value of the integrals *EPS_INT*, or only the diagonal terms are included
   (*DIAG*). *RECOMPUTE_TWO_INT_LIST_EVERY* allows to set how often the integral
   list is recomputed.

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
   active only when no thresholding of the integrals is used. Default is
   BLOCK_CYCLIC. The distribution is read from the next line.

.. keyword:: HFX_SCALE
   :target: HFX SCALE
   :section: &DFT
   :notes: scaling_for_exact_exchange

   Only available in combination with :kw:`XC_DRIVER <XC DRIVER>`. On the same
   line, a scaling value (:math:`\in [0,1]`) for the exact exchange contribution
   to the overall xc functional is specified. If used with a built-in hybrid
   functional, the default value is changed; if used with any other functional,
   the keyword activates the inclusion of exact exchange. Note that its value
   should always be set to one if combined with a Coulomb-attenuated functional
   (in which case the global scaling can be modified by changing the CAM
   parameters :math:`\alpha` and :math:`\beta`, instead).

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

.. keyword:: KERNEL_HFX_SCALE
   :target: KERNEL HFX SCALE
   :section: &DFT
   :notes: scaling_for_exact_exchange_in_kernel

   Identical functionality as described in :kw:`HFX_SCALE <HFX SCALE>`, but
   applies to the :kw:`LR KERNEL` rather than the
   :kw:`FUNCTIONAL`.

.. keyword:: KERNEL_LIBRARY
   :target: KERNEL LIBRARY
   :section: &DFT
   :notes: library_for_kernel_1, library_for_kernel_2, ...

   Same functionalty as :kw:`LIBRARY`, but applies to the
   :kw:`FUNCTIONAL` rather than the :kw:`LR KERNEL`.

.. keyword:: KERNEL_SCALES
   :target: KERNEL SCALES
   :section: &DFT
   :notes: scaling_for_kernel_1, scaling_for_kernel_2, ...

   Identical functionality as described in :kw:`SCALES`, but applies to
   the :kw:`LR KERNEL` rather than the
   :kw:`FUNCTIONAL`.

.. keyword:: LDA CORRELATION
   :section: &DFT
   :notes: [functional]

   The LDA correlation functional is specified. Only available for **OLDCODE** and
   **NEWCODE**.

   Possible functionals are **NO** (no correlation functional),
   **PZ** [Perdew81]_, **VWN** [Vosko80]_,
   **LYP** [Lee88]_ and **PW** [Perdew91]_.

   Default is the **PZ**, the Perdew and Zunger fit to the data of Ceperley and
   Alder [Ceperley80]_.

.. keyword:: LIBRARY
   :section: &DFT
   :notes: library_for_functional_1, library_for_functional_2, ...

   In combination with :kw:`XC_DRIVER <XC DRIVER>`, both functionals from a new,
   internal CP library and from libxc are available. On the same line, for every
   functional specified in **FUNCTIONAL**, the library has to be indicated by
   **INTERNAL** or **CP** for the internal CP library, or by **LIBXC** for libxc.
   If the keyword is omitted, the functionals will default to the internal CP
   library.

   The use of LIBXC is usually associated to a considerable overhead with respect
   to the internal implementation. It is therefore discouraged for functionals
   that are also available in the internal library, *cf.*
   :kw:`FUNCTIONAL`.

.. keyword:: LIBXC_ONLY
   :target: LIBXC ONLY
   :section: &DFT

   Use the new xc driver, but use only functionals from libxc (this excludes all
   internal functionals). *Cf.* :kw:`XC_DRIVER <XC DRIVER>`.

.. keyword:: LR KERNEL
   :section: &DFT
   :notes: functional

   Use another functional for the linear response kernel. To be used like
   :kw:`FUNCTIONAL`.

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
   can be invoked by using **OLD_DEFINITIONS**. The option is available in
   combination with any derivative of PBE and TPSS correlation.

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

   In combination with the PBE_FLEX exchange or correlation functionals, the
   parameter :math:`\kappa` can be specified on the next line.

.. keyword:: PBE_FLEX_MU
   :target: PBE FLEX MU
   :section: &DFT

   In combination with the PBE_FLEX exchange or correlation functionals, the
   parameter :math:`\mu` can be specified on the next line.

.. keyword:: PBE_FLEX_UEG_CORRELATION
   :target: PBE FLEX UEG CORRELATION
   :section: &DFT
   :notes: functional

   In combination with the PBE_FLEX correlation functional, the LDA correlation
   functional can be specified on the same line. Available options are LDA_C_PZ,
   LDA_C_PW (the default), LDA_C_OBPW and LDA_C_PZ.

.. keyword:: PHFX
   :section: &DFT

   Set the amount of HFX in the next line. The numerical value must be between 0
   (no HFX) and 1 (100% of HFX. For HSE06 the standard valus is 0.25 if this
   keyword and related numerical value ios not is not selected. Warning: PHFX is
   incompatible with XC_DRIVER. Use NEWCODE instead.

.. keyword:: REFUNCT
   :section: &DFT
   :notes: functionals

   Use a special reference functional in a calculation. This option is not active.

.. keyword:: SCALED EXCHANGE
   :section: &DFT

   Alias for :kw:`SCEX`.

.. keyword:: SCALES
   :section: &DFT
   :notes: scaling_for_functional_1, scaling_for_functional_2, ...

   Only available in combination with :kw:`XC_DRIVER <XC DRIVER>`. On the same
   line, rescaling values (:math:`\in [0,1]`) for the different
   exchange-correlation functionals specified in :kw:`FUNCTIONAL`
   are read in the same order in which the functionals were specified.

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
   Possible options are:

   ASHCROFT
      Ashcroft exchange. :math:`r_\text{cut}` is read from the next line.

   CAM
      Coulomb-attenuation method. :math:`\alpha`, :math:`\beta` and :math:`\mu`
      are read from the next line. Please note that in case attenuation of the
      (semi-)local part is desired, too, the keyword :kw:`COULOMB ATTENUATION` has to be used.

   ERFC
      Screened exchange using the complementary error function. The switching
      parameter :math:`\gamma` is read from the next line.

   EXP
      Screened exchange using an exponential function. The switching parameter
      :math:`\gamma` is read from the next line.

.. keyword:: SLATER
   :section: &DFT
   :alternatives: [NO]

   The :math:`\alpha` value for the Slater exchange
   functional [Slater51]_ is read from the next line. With NO the
   exchange functional is switched off.

   Default is a value of 2/3.

   This option together with no correlation functional, allows for :math:`X\alpha`
   theory.

.. keyword:: SMOOTH
   :section: &DFT

   A smoothening function is applied to the density [Laasonen93]_.

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

.. keyword:: XC_KERNEL
   :target: XC KERNEL
   :section: &DFT
   :notes: functional

   Alias for :kw:`LR KERNEL`.
