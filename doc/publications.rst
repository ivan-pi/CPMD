Publications
============

CPMD has been used in thousands of scientific studies since the early 1990s,
in fields ranging from liquid water and aqueous chemistry, through catalysis
and materials science, to enzymes and photobiology. This page lists the
papers that describe the code and its core methods, followed by a curated
selection of application papers that used CPMD.

The list is not complete. It is a representative sample, weighted towards
frequently cited work. The bibliographic data of every entry (title, journal,
volume, pages and DOI) was checked against the publisher's record. Note that
the Car-Parrinello method is also implemented in other codes such as CP2K and
Quantum ESPRESSO (``cp.x``); for a few application papers the use of the CPMD
program is stated in the abstract, for most it was inferred from the authors'
long-standing use and development of CPMD.

.. todo::

   MISSING: confirm the use of the CPMD program for each entry of Section B
   from the full text of the papers, and extend the list, in particular with
   recent work (see "Updating this list" below). Corrections and additions
   are welcome as pull requests.

Publications that report results obtained with CPMD should cite the code as
described in :ref:`citation`. Where appropriate, they
should also cite the method papers in Section A.


A. Papers describing CPMD and its methods
-----------------------------------------

The code and its parallel implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- R. Car, M. Parrinello, "Unified approach for molecular dynamics and
  density-functional theory", *Phys. Rev. Lett.* **55**, 2471 (1985).
  doi:`10.1103/PhysRevLett.55.2471 <https://doi.org/10.1103/PhysRevLett.55.2471>`__
- W. Andreoni, A. Curioni, "New advances in chemistry and materials science
  with CPMD and parallel computing", *Parallel Comput.* **26**, 819 (2000).
  doi:`10.1016/S0167-8191(00)00014-4 <https://doi.org/10.1016/S0167-8191(00)00014-4>`__
- J. Hutter, A. Curioni, "Dual-level parallelism for ab initio molecular
  dynamics: Reaching teraflop performance with the CPMD code",
  *Parallel Comput.* **31**, 1 (2005).
  doi:`10.1016/j.parco.2004.12.004 <https://doi.org/10.1016/j.parco.2004.12.004>`__
- J. Hutter, A. Curioni, "Car-Parrinello molecular dynamics on massively
  parallel computers", *ChemPhysChem* **6**, 1788 (2005).
  doi:`10.1002/cphc.200500059 <https://doi.org/10.1002/cphc.200500059>`__
- D. Marx, J. Hutter, *Ab Initio Molecular Dynamics: Basic Theory and Advanced
  Methods*, Cambridge University Press, Cambridge (2009).
  doi:`10.1017/CBO9780511609633 <https://doi.org/10.1017/CBO9780511609633>`__
- J. Hutter, "Car-Parrinello molecular dynamics", *WIREs Comput. Mol. Sci.*
  **2**, 604 (2012).
  doi:`10.1002/wcms.90 <https://doi.org/10.1002/wcms.90>`__
- T. Klöffel, G. Mathias, B. Meyer, "Integrating state of the art compute,
  communication, and autotuning strategies to multiply the performance of ab
  initio molecular dynamics on massively parallel multi-core supercomputers",
  *Comput. Phys. Commun.* **260**, 107745 (2021).
  doi:`10.1016/j.cpc.2020.107745 <https://doi.org/10.1016/j.cpc.2020.107745>`__

QM/MM and multiscale simulations (CPMD QM/MM interface, MiMiC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- A. Laio, J. VandeVondele, U. Rothlisberger, "A Hamiltonian electrostatic
  coupling scheme for hybrid Car-Parrinello molecular dynamics simulations",
  *J. Chem. Phys.* **116**, 6941 (2002).
  doi:`10.1063/1.1462041 <https://doi.org/10.1063/1.1462041>`__
- J. M. H. Olsen, V. Bolnykh, S. Meloni, et al., "MiMiC: A novel framework for
  multiscale modeling in computational chemistry", *J. Chem. Theory Comput.*
  **15**, 3810 (2019).
  doi:`10.1021/acs.jctc.9b00093 <https://doi.org/10.1021/acs.jctc.9b00093>`__
- V. Bolnykh, J. M. H. Olsen, S. Meloni, et al., "Extreme scalability of
  DFT-based QM/MM MD simulations using MiMiC", *J. Chem. Theory Comput.*
  **15**, 5601 (2019).
  doi:`10.1021/acs.jctc.9b00424 <https://doi.org/10.1021/acs.jctc.9b00424>`__
- A. Antalík, A. Levy, S. Kvedaravičiūtė, et al., "MiMiC: A high-performance
  framework for multiscale molecular dynamics simulations",
  *J. Chem. Phys.* **161**, 022501 (2024).
  doi:`10.1063/5.0211053 <https://doi.org/10.1063/5.0211053>`__

Methods implemented in CPMD
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- D. Marx, M. Parrinello, "Ab initio path integral molecular dynamics: Basic
  ideas", *J. Chem. Phys.* **104**, 4077 (1996).
  doi:`10.1063/1.471221 <https://doi.org/10.1063/1.471221>`__
- M. E. Tuckerman, D. Marx, M. L. Klein, et al., "Efficient and general
  algorithms for path integral Car-Parrinello molecular dynamics",
  *J. Chem. Phys.* **104**, 5579 (1996).
  doi:`10.1063/1.471771 <https://doi.org/10.1063/1.471771>`__
- A. Putrino, D. Sebastiani, M. Parrinello, "Generalized variational density
  functional perturbation theory", *J. Chem. Phys.* **113**, 7102 (2000).
  doi:`10.1063/1.1312830 <https://doi.org/10.1063/1.1312830>`__
- D. Sebastiani, M. Parrinello, "A new ab-initio approach for NMR chemical
  shifts in periodic systems", *J. Phys. Chem. A* **105**, 1951 (2001).
  doi:`10.1021/jp002807j <https://doi.org/10.1021/jp002807j>`__
- N. L. Doltsinis, D. Marx, "Nonadiabatic Car-Parrinello molecular dynamics",
  *Phys. Rev. Lett.* **88**, 166402 (2002).
  doi:`10.1103/PhysRevLett.88.166402 <https://doi.org/10.1103/PhysRevLett.88.166402>`__
- M. Iannuzzi, A. Laio, M. Parrinello, "Efficient exploration of reactive
  potential energy surfaces using Car-Parrinello molecular dynamics",
  *Phys. Rev. Lett.* **90**, 238302 (2003).
  doi:`10.1103/PhysRevLett.90.238302 <https://doi.org/10.1103/PhysRevLett.90.238302>`__
- O. A. von Lilienfeld, I. Tavernelli, U. Rothlisberger, et al., "Optimization
  of effective atom centered potentials for London dispersion forces in
  density functional theory", *Phys. Rev. Lett.* **93**, 153004 (2004).
  doi:`10.1103/PhysRevLett.93.153004 <https://doi.org/10.1103/PhysRevLett.93.153004>`__
- I. Tavernelli, U. F. Röhrig, U. Rothlisberger, "Molecular dynamics in
  electronically excited states using time-dependent density functional
  theory", *Mol. Phys.* **103**, 963 (2005).
  doi:`10.1080/00268970512331339378 <https://doi.org/10.1080/00268970512331339378>`__
- E. Tapavicza, I. Tavernelli, U. Rothlisberger, "Trajectory surface hopping
  within linear response time-dependent density-functional theory",
  *Phys. Rev. Lett.* **98**, 023001 (2007).
  doi:`10.1103/PhysRevLett.98.023001 <https://doi.org/10.1103/PhysRevLett.98.023001>`__
- M. P. Bircher, U. Rothlisberger, "Plane-wave implementation and performance
  of à-la-carte Coulomb-attenuated exchange-correlation functionals for
  predicting optical excitation energies in some notorious cases",
  *J. Chem. Theory Comput.* **14**, 3184 (2018).
  doi:`10.1021/acs.jctc.8b00069 <https://doi.org/10.1021/acs.jctc.8b00069>`__


B. Selected publications using CPMD
-----------------------------------

Liquid water, aqueous ions and proton transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- M. Tuckerman, K. Laasonen, M. Sprik, et al., "Ab initio molecular dynamics
  simulation of the solvation and transport of hydronium and hydroxyl ions in
  water", *J. Chem. Phys.* **103**, 150 (1995).
  doi:`10.1063/1.469654 <https://doi.org/10.1063/1.469654>`__
- M. Sprik, J. Hutter, M. Parrinello, "Ab initio molecular dynamics simulation
  of liquid water: Comparison of three gradient-corrected density
  functionals", *J. Chem. Phys.* **105**, 1142 (1996).
  doi:`10.1063/1.471957 <https://doi.org/10.1063/1.471957>`__
- D. Marx, M. E. Tuckerman, J. Hutter, et al., "The nature of the hydrated
  excess proton in water", *Nature* **397**, 601 (1999).
  doi:`10.1038/17579 <https://doi.org/10.1038/17579>`__
- P. L. Silvestrelli, M. Parrinello, "Water molecule dipole in the gas and in
  the liquid phase", *Phys. Rev. Lett.* **82**, 3308 (1999).
  doi:`10.1103/PhysRevLett.82.3308 <https://doi.org/10.1103/PhysRevLett.82.3308>`__
- P. L. Silvestrelli, M. Parrinello, "Structural, electronic, and bonding
  properties of liquid water from first principles", *J. Chem. Phys.*
  **111**, 3572 (1999).
  doi:`10.1063/1.479638 <https://doi.org/10.1063/1.479638>`__
- M. Boero, K. Terakura, T. Ikeshoji, et al., "Hydrogen bonding and dipole
  moment of water at supercritical conditions: A first-principles molecular
  dynamics study", *Phys. Rev. Lett.* **85**, 3245 (2000).
  doi:`10.1103/PhysRevLett.85.3245 <https://doi.org/10.1103/PhysRevLett.85.3245>`__
- P. L. Geissler, C. Dellago, D. Chandler, et al., "Autoionization in liquid
  water", *Science* **291**, 2121 (2001).
  doi:`10.1126/science.1056991 <https://doi.org/10.1126/science.1056991>`__
- R. Vuilleumier, M. Sprik, "Electronic properties of hard and soft ions in
  solution: Aqueous Na+ and Ag+ compared", *J. Chem. Phys.* **115**, 3454
  (2001).
  doi:`10.1063/1.1388901 <https://doi.org/10.1063/1.1388901>`__
- M. E. Tuckerman, D. Marx, M. Parrinello, "The nature and transport mechanism
  of hydrated hydroxide ions in aqueous solution", *Nature* **417**, 925
  (2002).
  doi:`10.1038/nature00797 <https://doi.org/10.1038/nature00797>`__
- J. Blumberger, L. Bernasconi, I. Tavernelli, et al., "Electronic structure
  and solvation of copper and silver ions: A theoretical picture of a model
  aqueous redox reaction", *J. Am. Chem. Soc.* **126**, 3928 (2004).
  doi:`10.1021/ja0390754 <https://doi.org/10.1021/ja0390754>`__
- A. Chandra, M. E. Tuckerman, D. Marx, "Connecting solvation shell structure
  to proton transport kinetics in hydrogen-bonded networks via population
  correlation functions", *Phys. Rev. Lett.* **99**, 145901 (2007).
  doi:`10.1103/PhysRevLett.99.145901 <https://doi.org/10.1103/PhysRevLett.99.145901>`__

Nuclear quantum effects and matter under high pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- D. Marx, M. Parrinello, "Structural quantum effects and three-centre
  two-electron bonding in CH5+", *Nature* **375**, 216 (1995).
  doi:`10.1038/375216a0 <https://doi.org/10.1038/375216a0>`__
- M. Benoit, M. Bernasconi, P. Focher, et al., "New high-pressure phase of
  ice", *Phys. Rev. Lett.* **76**, 2934 (1996).
  doi:`10.1103/PhysRevLett.76.2934 <https://doi.org/10.1103/PhysRevLett.76.2934>`__
- M. E. Tuckerman, D. Marx, M. L. Klein, et al., "On the quantum nature of the
  shared proton in hydrogen bonds", *Science* **275**, 817 (1997).
  doi:`10.1126/science.275.5301.817 <https://doi.org/10.1126/science.275.5301.817>`__
- M. Benoit, D. Marx, M. Parrinello, "Tunnelling and zero-point motion in
  high-pressure ice", *Nature* **392**, 258 (1998).
  doi:`10.1038/32609 <https://doi.org/10.1038/32609>`__
- M. Bernasconi, P. L. Silvestrelli, M. Parrinello, "Ab initio infrared
  absorption study of the hydrogen-bond symmetrization in ice",
  *Phys. Rev. Lett.* **81**, 1235 (1998).
  doi:`10.1103/PhysRevLett.81.1235 <https://doi.org/10.1103/PhysRevLett.81.1235>`__

Vibrational spectroscopy and NMR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- P. L. Silvestrelli, M. Bernasconi, M. Parrinello, "Ab initio infrared
  spectrum of liquid water", *Chem. Phys. Lett.* **277**, 478 (1997).
  doi:`10.1016/S0009-2614(97)00930-5 <https://doi.org/10.1016/S0009-2614(97)00930-5>`__
- M. Bühl, M. Parrinello, "Medium effects on 51V NMR chemical shifts: A density
  functional study", *Chem. Eur. J.* **7**, 4487 (2001).
  doi:`10.1002/1521-3765(20011015)7:20\<4487::AID-CHEM4487\>3.0.CO;2-G <https://doi.org/10.1002/1521-3765(20011015)7:20%3C4487::AID-CHEM4487%3E3.0.CO;2-G>`__
- M. Pagliai, G. Cardini, R. Righini, et al., "Hydrogen bond dynamics in
  liquid methanol", *J. Chem. Phys.* **119**, 6655 (2003).
  doi:`10.1063/1.1605093 <https://doi.org/10.1063/1.1605093>`__
- M.-P. Gaigeot, M. Sprik, "Ab initio molecular dynamics computation of the
  infrared spectrum of aqueous uracil", *J. Phys. Chem. B* **107**, 10344
  (2003).
  doi:`10.1021/jp034788u <https://doi.org/10.1021/jp034788u>`__
- M.-P. Gaigeot, R. Vuilleumier, M. Sprik, et al., "Infrared spectroscopy of
  N-methylacetamide revisited by ab initio molecular dynamics simulations",
  *J. Chem. Theory Comput.* **1**, 772 (2005).
  doi:`10.1021/ct050029z <https://doi.org/10.1021/ct050029z>`__

Chemical reactions, catalysis and mechanochemistry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- A. Curioni, M. Sprik, W. Andreoni, et al., "Density functional theory-based
  molecular dynamics simulation of acid-catalyzed chemical reactions in liquid
  trioxane", *J. Am. Chem. Soc.* **119**, 7218 (1997).
  doi:`10.1021/ja970935o <https://doi.org/10.1021/ja970935o>`__
- M. Boero, M. Parrinello, K. Terakura, "First principles molecular dynamics
  study of Ziegler-Natta heterogeneous catalysis", *J. Am. Chem. Soc.*
  **120**, 2746 (1998).
  doi:`10.1021/ja972367i <https://doi.org/10.1021/ja972367i>`__
- L. Guidoni, K. Spiegel, M. Zumstein, et al., "Green oxidation catalysts:
  Computational design of high-efficiency models of galactose oxidase",
  *Angew. Chem. Int. Ed.* **43**, 3286 (2004).
  doi:`10.1002/anie.200454081 <https://doi.org/10.1002/anie.200454081>`__
- M. Bühl, A. Chaumont, R. Schurhammer, et al., "Ab initio molecular dynamics
  of liquid 1,3-dimethylimidazolium chloride", *J. Phys. Chem. B* **109**,
  18591 (2005).
  doi:`10.1021/jp0518299 <https://doi.org/10.1021/jp0518299>`__
- N. N. Nair, E. Schreiner, D. Marx, "Glycine at the pyrite-water interface:
  The role of surface defects", *J. Am. Chem. Soc.* **128**, 13815 (2006).
  doi:`10.1021/ja063295a <https://doi.org/10.1021/ja063295a>`__
- E. Schreiner, N. N. Nair, D. Marx, "Influence of extreme thermodynamic
  conditions and pyrite surfaces on peptide synthesis in aqueous media",
  *J. Am. Chem. Soc.* **130**, 2768 (2008).
  doi:`10.1021/ja7108085 <https://doi.org/10.1021/ja7108085>`__
- P. Dopieralski, J. Ribas-Arino, P. Anjukandi, et al., "The Janus-faced role
  of external forces in mechanochemical disulfide bond cleavage",
  *Nat. Chem.* **5**, 685 (2013).
  doi:`10.1038/nchem.1676 <https://doi.org/10.1038/nchem.1676>`__
- P. Dopieralski, J. Ribas-Arino, P. Anjukandi, et al., "Unexpected
  mechanochemical complexity in the mechanistic scenarios of disulfide bond
  reduction in alkaline solution", *Nat. Chem.* **9**, 164 (2017).
  doi:`10.1038/nchem.2632 <https://doi.org/10.1038/nchem.2632>`__

Materials, surfaces and nanostructures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- W. Andreoni, A. Curioni, "Freedom and constraints of a metal atom
  encapsulated in fullerene cages", *Phys. Rev. Lett.* **77**, 834 (1996).
  doi:`10.1103/PhysRevLett.77.834 <https://doi.org/10.1103/PhysRevLett.77.834>`__
- P. Carloni, W. Andreoni, M. Parrinello, "Self-assembled peptide nanotubes
  from first principles", *Phys. Rev. Lett.* **79**, 761 (1997).
  doi:`10.1103/PhysRevLett.79.761 <https://doi.org/10.1103/PhysRevLett.79.761>`__
- K. C. Hass, W. F. Schneider, A. Curioni, et al., "The chemistry of water on
  alumina surfaces: Reaction dynamics from first principles", *Science*
  **282**, 265 (1998).
  doi:`10.1126/science.282.5387.265 <https://doi.org/10.1126/science.282.5387.265>`__
- A. Curioni, M. Boero, W. Andreoni, "Alq3: ab initio calculations of its
  structural and electronic properties in neutral and charged states",
  *Chem. Phys. Lett.* **294**, 263 (1998).
  doi:`10.1016/S0009-2614(98)00829-X <https://doi.org/10.1016/S0009-2614(98)00829-X>`__
- H. Grönbeck, A. Curioni, W. Andreoni, "Thiols and disulfides on the Au(111)
  surface: The headgroup-gold interaction", *J. Am. Chem. Soc.* **122**, 3839
  (2000).
  doi:`10.1021/ja993622x <https://doi.org/10.1021/ja993622x>`__
- B. Meyer, D. Marx, "Density-functional study of the structure and stability
  of ZnO surfaces", *Phys. Rev. B* **67**, 035403 (2003).
  doi:`10.1103/PhysRevB.67.035403 <https://doi.org/10.1103/PhysRevB.67.035403>`__

Biomolecules and enzymes (including QM/MM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- C. Rovira, K. Kunc, J. Hutter, et al., "Equilibrium geometries and
  electronic structure of iron-porphyrin complexes: A density functional
  study", *J. Phys. Chem. A* **101**, 8914 (1997).
  doi:`10.1021/jp9722115 <https://doi.org/10.1021/jp9722115>`__
- P. Carloni, M. Sprik, W. Andreoni, "Key steps of the cis-platin-DNA
  interaction: Density functional theory-based molecular dynamics
  simulations", *J. Phys. Chem. B* **104**, 823 (2000).
  doi:`10.1021/jp992590x <https://doi.org/10.1021/jp992590x>`__
- S. Piana, P. Carloni, M. Parrinello, "Role of conformational fluctuations in
  the enzymatic reaction of HIV-1 protease", *J. Mol. Biol.* **319**, 567
  (2002).
  doi:`10.1016/S0022-2836(02)00301-7 <https://doi.org/10.1016/S0022-2836(02)00301-7>`__
- M. Sulpizi, A. Laio, J. VandeVondele, et al., "Reaction mechanism of
  caspases: Insights from QM/MM Car-Parrinello simulations", *Proteins*
  **52**, 212 (2003).
  doi:`10.1002/prot.10275 <https://doi.org/10.1002/prot.10275>`__
- K. Spiegel, U. Rothlisberger, P. Carloni, "Cisplatin binding to DNA
  oligomers from hybrid Car-Parrinello/molecular dynamics simulations",
  *J. Phys. Chem. B* **108**, 2699 (2004).
  doi:`10.1021/jp036230s <https://doi.org/10.1021/jp036230s>`__
- M. De Vivo, B. Ensing, M. L. Klein, "Computational study of phosphatase
  activity in soluble epoxide hydrolase: High efficiency through a water
  bridge mediated proton shuttle", *J. Am. Chem. Soc.* **127**, 11226 (2005).
  doi:`10.1021/ja053049j <https://doi.org/10.1021/ja053049j>`__
- M. Dal Peraro, A. J. Vila, P. Carloni, et al., "Role of zinc content on the
  catalytic efficiency of B1 metallo beta-lactamases", *J. Am. Chem. Soc.*
  **129**, 2808 (2007).
  doi:`10.1021/ja0657556 <https://doi.org/10.1021/ja0657556>`__
- X. Biarnés, A. Ardèvol, A. Planas, et al., "The conformational free energy
  landscape of beta-D-glucopyranose. Implications for substrate preactivation
  in beta-glucoside hydrolases", *J. Am. Chem. Soc.* **129**, 10686 (2007).
  doi:`10.1021/ja068411o <https://doi.org/10.1021/ja068411o>`__
- M. De Vivo, M. Dal Peraro, M. L. Klein, "Phosphodiester cleavage in
  ribonuclease H occurs via an associative two-metal-aided catalytic
  mechanism", *J. Am. Chem. Soc.* **130**, 10955 (2008).
  doi:`10.1021/ja8005786 <https://doi.org/10.1021/ja8005786>`__
- A. Ardèvol, C. Rovira, "The molecular mechanism of enzymatic glycosyl
  transfer with retention of configuration: Evidence for a short-lived
  oxocarbenium-like species", *Angew. Chem. Int. Ed.* **50**, 10897 (2011).
  doi:`10.1002/anie.201104623 <https://doi.org/10.1002/anie.201104623>`__
- A. Ardèvol, C. Rovira, "Reaction mechanisms in carbohydrate-active enzymes:
  Glycoside hydrolases and glycosyltransferases. Insights from ab initio
  quantum mechanics/molecular mechanics dynamic simulations",
  *J. Am. Chem. Soc.* **137**, 7528 (2015).
  doi:`10.1021/jacs.5b01156 <https://doi.org/10.1021/jacs.5b01156>`__
- M. G. Chiariello, V. Bolnykh, E. Ippoliti, et al., "Molecular basis of CLC
  antiporter inhibition by fluoride", *J. Am. Chem. Soc.* **142**, 7254
  (2020).
  doi:`10.1021/jacs.9b13588 <https://doi.org/10.1021/jacs.9b13588>`__
- B. Raghavan, M. Paulikat, K. Ahmad, et al., "Drug design in the exascale
  era: A perspective from massively parallel QM/MM simulations",
  *J. Chem. Inf. Model.* **63**, 3647 (2023).
  doi:`10.1021/acs.jcim.3c00557 <https://doi.org/10.1021/acs.jcim.3c00557>`__
- G. L. Hoang, M. Röck, A. Tancredi, et al., "Refining ligand poses in
  RNA/ligand complexes of pharmaceutical relevance: A perspective by QM/MM
  simulations and NMR measurements", *J. Phys. Chem. Lett.* **16**, 1702
  (2025).
  doi:`10.1021/acs.jpclett.4c03456 <https://doi.org/10.1021/acs.jpclett.4c03456>`__

Excited states, photochemistry and photobiology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- U. F. Röhrig, I. Frank, J. Hutter, et al., "QM/MM Car-Parrinello molecular
  dynamics study of the solvent effects on the ground state and on the first
  excited singlet state of acetone in water", *ChemPhysChem* **4**, 1177
  (2003).
  doi:`10.1002/cphc.200300650 <https://doi.org/10.1002/cphc.200300650>`__
- U. F. Röhrig, L. Guidoni, A. Laio, et al., "A molecular spring for vision",
  *J. Am. Chem. Soc.* **126**, 15328 (2004).
  doi:`10.1021/ja048265r <https://doi.org/10.1021/ja048265r>`__
- O. Valsson, P. Campomanes, I. Tavernelli, et al., "Rhodopsin absorption from
  first principles: Bypassing common pitfalls", *J. Chem. Theory Comput.*
  **9**, 2441 (2013).
  doi:`10.1021/ct3010408 <https://doi.org/10.1021/ct3010408>`__
- P. Campomanes, M. Neri, B. A. C. Horta, et al., "Origin of the spectral
  shifts among the early intermediates of the rhodopsin photocycle",
  *J. Am. Chem. Soc.* **136**, 3842 (2014).
  doi:`10.1021/ja411303v <https://doi.org/10.1021/ja411303v>`__


Updating this list
------------------

Candidate papers for this page can be found with the OpenAlex API, for
example:

- Look up the OpenAlex record (and its ``W...`` identifier) of the Andreoni and
  Curioni (2000) CPMD paper:
  https://api.openalex.org/works/doi:10.1016/S0167-8191(00)00014-4
- List the works that cite it, most cited first (replace ``W...`` with the
  identifier from the previous query):
  ``https://api.openalex.org/works?filter=cites:W...&sort=cited_by_count:desc&per-page=50``
- Search titles and abstracts that mention the code explicitly:
  ``https://api.openalex.org/works?filter=title_and_abstract.search:%22CPMD%20code%22&sort=cited_by_count:desc&per-page=50``
  (and likewise ``%22CPMD%20program%22`` and ``%22CPMD%20package%22``)

Search hits must be checked by hand before they are added. "CPMD" is also
widely used as an abbreviation of the Car-Parrinello method, which CP2K,
Quantum ESPRESSO and other codes implement.
