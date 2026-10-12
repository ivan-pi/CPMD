Using the new xc driver and assembling functionals
==================================================

The new :kw:`XC_DRIVER <XC DRIVER>` offers more flexibility compared to
**OLDCODE** and **NEWCODE**. The minimal setup for the use of, *e.g.*, the BLYP
xc functional:

::

   &DFT
    XC_DRIVER
    FUNCTIONAL GGA_XC_BLYP
   &END

is equivalent to specifying

::

   &DFT
    XC_DRIVER
    FUNCTIONAL GGA_X_BLYP GGA_C_LYP
   &END

or, in a more verbose way,

::

   &DFT
    XC_DRIVER
    FUNCTIONAL GGA_X_BLYP GGA_C_LYP
    SCALES     1.0        1.0
    LIBRARY    CP         CP
   &END

If not explicitly stated, the xc driver will always use the internal CP
functional library. The use of libxc functionals can be requested by using the
keyword :kw:`LIBRARY`, and both libxc and internal CP functionals
may be used simultaneously. The following example combines a particular
parametrisation of VWN available in libxc with Slater exchange computed within
the internal CP library.

::

   &DFT
    XC_DRIVER
    FUNCTIONAL LDA_X LDA_C_VWN_3
    LIBRARY    CP    LIBXC
   &END

If not specified, all functionals will contribute equally to the total xc
energy and the potential. A scaling factor for every functional can be
introduced by using :kw:`SCALES`. A manual setup equivalent to
HYB_GGA_XC_CAM_B3LYP, but using different attenuation parameters
:math:`\alpha=0.19, \beta=0.66, \mu=0.330` would read:

::

   &DFT
    XC_DRIVER
    FUNCTIONAL GGA_X_B88 GGA_C_LYP LDA_C_VWN
    SCALES     1.00      0.81      0.19
    COULOMB ATTENUATION
    0.190  0.660  0.330
   &END

It is also possible to use a linear response kernel different from the one used
to generate the ground-state orbitals. In this case, duplicates of the
aforementioned exist with the prefix **KERNEL\_**. Note that is currently not
possible to combine different forms of the Coulomb operator between the kernel
and the ground-state functional (*i.e.* it is not possible to use the
combination of B3LYP and CAM-B3LYP; however, combinations between GGA and any
(screened) hybrid are possible).

::

   &DFT
    XC_DRIVER
    FUNCTIONAL       GGA_X_B88 GGA_C_LYP LDA_C_VWN_3
    LIBRARY          CP        CP        LIBXC
    SCALES           0.80      0.81      0.19
    HFX_SCALE        0.20
    LR_KERNEL        GGA_X_B88 GGA_C_LYP LDA_C_VWN
    KERNEL_LIBRARY   CP        CP        CP
    KERNEL_SCALES    0.75      0.81      0.19
    KERNEL_HFX_SCALE 0.25
   &END
