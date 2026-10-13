.. _input-mts:

``&MTS`` section
================

Parameters for the multiple time step MD scheme.

The section starts with ``&MTS`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &MTS

Keywords
--------

.. keyword:: HIGH_LEVEL_FORCES
   :target: HIGH LEVEL FORCES
   :section: &MTS
   :options: {DFT, EXTERNAL}

   Set the computational model for the calculation of high level forces in the MTS
   scheme. For now two options exists: (i) The forces are obtained from the DFT
   code internal to CPMD. In that case, the high level functional has to be set in
   the &DFT section with the :kw:`MTS_HIGH_FUNC <MTS HIGH FUNC>` keyword. (ii)
   The forces are calculated by an external program which called via a script
   which should be named ``EXT_HIGH_FORCES``. When the forces are needed, CPMD
   will write a file ``geometry.xyz`` containing the current geometry in xyz
   format (cartesian coordinates in Angstroms). Then the script will be called and
   CPMD will wait until a file ``forces.xyz`` becomes available with the new
   forces in xyz format (atomic units are expected). See the
   ``get_external_forces`` routine in the ``interface_utils.mod.F90`` file.

   The keyword :kw:`LOW_LEVEL_FORCES <LOW LEVEL FORCES>` can be used in exactly
   the same way for the selection of the computational model used to calculate the
   low level forces in the MTS scheme. Similarly, if the DFT code internal to CPMD
   is chosen as low level, the keyword :kw:`MTS_LOW_FUNC <MTS LOW FUNC>` (&DFT
   section), can be used to select the low level density functional.

   **Default:** Forces are calculated with the DFT code internal to CPMD.

.. keyword:: LOW_LEVEL_FORCES
   :target: LOW LEVEL FORCES
   :section: &MTS
   :options: {DFT, EXTERNAL}

   Analogue to :kw:`HIGH_LEVEL_FORCES <HIGH LEVEL FORCES>`.

.. keyword:: PRINT_FORCES
   :target: PRINT FORCES
   :section: &MTS
   :options: [OFF]

   Turn on/off the printing of the high and low level forces along the MTS-MD
   trajectory. The forces are printed to two distinct trajectory files in the
   ``.xyz`` format, ``MTS_LOW_FORCES.xyz`` and ``MTS_HIGH_FORCES.xyz``.

   **Default**: OFF.

.. keyword:: TIMESTEP_FACTOR
   :target: TIMESTEP FACTOR
   :section: &MTS

   The time step factor (:math:`n`) used in the MTS scheme is read from the next
   line.

   For instance, if :kw:`TIMESTEP IONS` is set to **15 a.u.**,
   and (:math:`n=4`) then the high level correction to the forces in the MTS
   scheme will be calculated every 4 steps, i.e. every **60 a.u.** **Default** is
   a factor of 1, i.e. the MD will be identical to a Velocity-Verlet MD with the
   high level forces.

.. keyword:: USE_ACE
   :section: &MTS

   Switches on the ACE formulation of the Hartree-Fock exchange (projector
   operator built in ``build_ace_proj``, ``src/hfx_utils.mod.F90``) within
   the multiple-time-step scheme. It sets the flag ``use_ace`` in module
   ``ace_hfx``; during the MTS loop the low-level DFT force evaluation is
   then flagged with ``status_ace=.TRUE.`` and the high-level one with
   ``status_ace=.FALSE.``, and the exchange operator is applied through the
   ACE projectors. With ``use_ace`` and :kw:`EXTRAPOLATE WFN` active, the
   extrapolation of the high-level wavefunction is replaced by a copy of
   the current wavefunction. The section is only read when :kw:`USE_MTS` is given in
   ``&CPMD``. **Default** is not to use ACE. The value is printed as
   ``USING ACE APPROACH`` in the MTS summary.

   Example::

      &MTS
        TIMESTEP_FACTOR
          4
        USE_ACE
      &END

   .. todo::

      MISSING: meaning of the acronym ACE, literature reference, which
      combinations of :kw:`LOW_LEVEL_FORCES`/:kw:`HIGH_LEVEL_FORCES` and
      hybrid functionals are supported, and whether the SCDM localization
      used in ``hfx_utils`` must be enabled separately.

.. source: src/mts_utils.mod.F90:130
