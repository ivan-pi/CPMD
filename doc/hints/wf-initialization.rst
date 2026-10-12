Wavefunction Initialization
===========================

The default initial guess for the wavefunctions is calculated from the atomic
pseudo-wavefunctions and gives usually excellent results. Good results can also
be obtained by using wavefunctions from other calculations with a different
cutoff or slightly different geometry. The other initialization available,
starts from random coefficients and should only be used as a last resort. Cases
where the default method does not work are when the molecule has less occupied
states than one of the atoms (in this case add some empty states for the
molecule) or when the additional memory required for the atomic calculation is
not available.

Using Vanderbilt Ultrasoft Pseudopotentials
-------------------------------------------

When using Vanderbilt ultrasoft pseudopotentials (USPPs)
[Vanderbilt]_ and starting from atomic pseudo-wavefunctions, the
calculations often do not converge or converge to a wrong state, especially if
3d-elements are involved. Convergence is generally much better when assigning
(partial) charges via the :kw:`ATOMIC CHARGES` keyword in the
&SYSTEM …&END section. Values from a classical MD forcefield or an NBO
calculation are usually good values. Alternatively a random initialization of
the wave functions (via :kw:`INITIALIZE WAVEFUNCTION` RANDOM) can be used.

Also, due to the comparatively small plane wave cutoffs, you will have small
but significant modulations of the density in especially in regions with little
electron density. These lead to "strange" effects with gradient corrected
functionals, causing the optimization to fail. To avoid this, you can skip the
calculation of the gradient correction for low electron density areas using
:kw:`GC-CUTOFF` with a value between 1.D-6 and 1.D-5 in the &DFT
section.

In case of geometry optimizations, also the accurate calculation of the forces
due to the augmentation charges may need a higher density cutoff and/or a
tighter real space grid. This can be achieved by either using a higher plane
wave cutoff or via increasing :kw:`DUAL` to 5.0 or even 6.0 and/or
setting the real space grid explicitly via the :kw:`MESH` keyword in
the &SYSTEM section. For the same reason, these options may be needed to
increase energy conservation during molecular dynamics runs. Use these options
with care, as they will increase the cpu time and memory requirements
significantly and thus can easily take away one of the major advantages of
ultra-soft pseudopotentials.
