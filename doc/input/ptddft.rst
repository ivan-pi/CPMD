.. _input-ptddft:

``&PTDDFT`` section
===================

Propagation TDDFT for Ehrenfest dynamics and spectra. This section contains
specific information for Ehrenfest dynamics and for spectra computed by the
propagation of the perturbed Kohn-Sham orbitals (Fourier transform of the
induced dipole fluctuation) [taver_eh]_ [taver1]_.

The section starts with ``&PTDDFT`` and ends with ``&END``.

Summary
-------

.. keyword-summary:: &PTDDFT

Keywords
--------

.. keyword:: ACCURACY
   :section: &PTDDFT

   Specifies the accuracy to be reached in the Cayley propagation scheme used in
   Ehrenfest type of dynamics and spectra calculation.

.. keyword:: EXT_POTENTIAL
   :target: EXT-POTENTIAL
   :section: &PTDDFT

   An external potential is applied to the KS orbitals before propagation. This
   keyword is used with the principal keyword **MOLECULAR DYNAMICS EH** or
   :kw:`PROPAGATION SPECTRA`. The type of the perturbation
   is specified with the keyword **PERT_TYPE** (in &PTDDFT) and its intensity with
   **PERT_AMPLI** (in &PTDDFT).

.. keyword:: EXT_PULSE
   :target: EXT-PULSE
   :section: &PTDDFT

   A Dirac type pulse is applied to the KS orbitals before propagation. This
   keyword is used with the principal keyword :kw:`PROPAGATION SPECTRA`. The intensity of the perturbing field is
   read from the next line.

.. keyword:: N_CYCLES
   :target: N-CYCLES
   :section: &PTDDFT

   Defines the number of cycles (time steps) used in the propagation of the
   perturbed KS orbitals in a :kw:`PROPAGATION SPECTRA`
   calculation. The number of cycles is read from the next line.

.. keyword:: PERT_AMPLI
   :target: PERT-AMPLI
   :section: &PTDDFT

   The amplitude of the perturbation used with :kw:`MOLECULAR DYNAMICS` EH or :kw:`PROPAGATION SPECTRA` is read from the next line.

.. keyword:: PERT_DIRECTION
   :target: PERT-DIRECTION
   :section: &PTDDFT

   Use with MOLECULAR DYNAMICS EH or PROPAGATION SPECTRA to specify the direction
   (x=1,y=2,z=3) of the perturbing field. The direction code (1,2,3) is read from
   the next line.

.. keyword:: PERT_TYPE
   :target: PERT-TYPE
   :section: &PTDDFT

   The type of the static perturbation used with :kw:`MOLECULAR DYNAMICS` EH or :kw:`PROPAGATION SPECTRA` is read from the next line. 1: Dirac pulse,
   2: Heaviside step function (constant value after :math:`t=0`).

.. keyword:: PIPULSE
   :section: &PTDDFT

   Specifies a time dependent pi-pulse to be used with MOLECULAR DYNAMICS EH. Use
   PIPULSE together with TD_POTENTIAL. The pulse strength is read from the next
   line (see subroutine gaugepot_laser in td_utils.mod.F90 for further details).

.. keyword:: PROP_TSTEP
   :target: PROP-TSTEP
   :section: &PTDDFT

   Propagation timestep used in Ehrenfest dynamics. It is used in the spectra
   calculation (**PROPAGATION SPECTRA** ) to specify the time step for the
   propagation of the KS orbitals.

.. keyword:: RESTFILE
   :target: RESTFILE PTDDFT
   :section: &PTDDFT

   Defines a restart code for the restart of the Ehrenfest dynamics (:kw:`MOLECULAR DYNAMICS` EH) and the propagation spectra
   (:kw:`PROPAGATION SPECTRA`). The restart option is read
   from the next line: 0(=default) restart from the (complex)wavefunctions in the
   file wavefunctions. This option is used in case of a continuation run; 1.
   restart from the orbital files WAVEFUNCTION.n, where :math:`n` is the index of
   the KS orbital and runs from :math:`1` to the number of states (This states a
   prepare in a previous run using the KOHN-SHAM ENERGIES principal keyword), 2;
   restart from the orbitals stored in RESTART (obtained from a optimization run
   with tight convergence (at least 1.D-7)).

.. keyword:: TD_POTENTIAL
   :target: TD-POTENTIAL
   :section: &PTDDFT

   Defines a time dependent external potential to be used in Ehrenfest dynamics
   (:kw:`MOLECULAR DYNAMICS` EH). Can be used with the
   keyword :kw:`PIPULSE`. The frequency of the external field is read
   from the next line (in atomic units). For more information see subroutine
   gaugepot_laser in td_utils.mod.F90 .
