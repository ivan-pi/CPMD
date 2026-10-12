CPMD on parallel computers
==========================

There are three different parallel strategies implemented in CPMD. The actual
strategy has to be chosen at compile time.

-  Shared memory parallelization

   This strategy uses the OpenMP library. The code is compiled *without* the
   **PARALLEL** preprocessor flag and compilation and linking need the
   corresponding OpenMP flags (dependent on compiler).

   Depending on the overhead of the OpenMP system implementation good speedups
   can be achieved for small numbers of processors (typically 4 to 8). The
   advantages of this version of the code are small additional memory usage and
   it can be used in non-dedicated CPU environments.

-  Distributed memory parallelization

   This is the standard parallelization scheme used in CPMD based on MPI
   message passing library.

   The single processor version of this code typically shows an overhead of ca.
   10% with respect to the optimal serial code. This overhead is due to
   additional copy and sort operations during the FFTs.

   All the basic system data and many matrices of the size of the number of
   electrons are replicated on all processors. This leads to considerable
   additional memory usage (calculated as the sum over the memory of all
   processors compared to the memory needed on a single processor). For large
   systems distributed over many processors the replicated data can dominate
   the memory usage.

   The efficiency of the parallelization depends on the calculated system (e.g.
   cutoff and number of electrons) and the hardware platform, mostly latency
   and bandwidth of the communication system. The most important bottleneck in
   the distributed memory parallelization of CPMD is the load-balancing problem
   in the FFT. The real space grids are distributed over the first dimension
   alone (see line ``REAL SPACE MESH:`` in the output. As the mesh sizes only
   vary between 20 (very small systems, low cutoffs) and 300 (large systems,
   high cutoff) we have a rather coarse grain parallelization. To avoid load
   unbalance the number of processors should be a divisor of the mesh size. It
   is therefore clear that even for large systems no speedup can be achieved
   beyond 300 processors. A solution to this problem is provided with the
   keyword :kw:`CP GROUPS`. This technique, together with optimal mapping,
   allow to scale from thousands to millions of cores on modern supercomputers
   such as IBM BG/Q (in particular when running HFX calculations). To learn
   more about the distributed memory parallelization of CPMD consult D. Marx
   and J. Hutter, "Modern Methods and Algorithms of Quantum Chemistry",
   Forschungszentrum Jülich, NIC Series, Vol. 1 (2000), 301-449. For recent
   developments see https://github.com/OpenCPMD.

   .. figure:: ../fig/fft-flowchart.png
      :width: 80%
      :alt: Flowcharts of the density calculation and of the application of the local potential

      The two basic operations that dominate the use of FFTs in a plane wave
      code. Left: the orbitals :math:`c_i(\mathbf{G})` are transformed to real
      space and their squares are accumulated to the density
      :math:`n(\mathbf{R})`. Right: the local potential is applied to an
      orbital in real space, and the product is transformed back to
      reciprocal space.

   When selecting **NSTBLK** for :kw:`BLOCKSIZE STATES` it
   is important to take into account the granularity of the problem at hand.
   For example, in cases where the number of :kw:`STATES` is smaller
   than the total number of the available processors, one must choose a value
   for **NSTBLK** such that only a subgroup of the processors participate in
   the distributed linear algebra calculations. The same argument is also
   relevant when the number of :kw:`STATES` is only moderately larger
   than the number of processors.

-  Mixed shared/distributed memory parallelization

   The two parallelization schemes described above are implemented in such a
   way that they don’t interfere. Therefore it is easy to combine them if this
   is supported by hardware (shared/distributed memory architecture) and
   software (libraries). Since the MPI parallelization is very efficient for a
   small to medium number of nodes and all modern MPI libraries are able to
   take advantage from shared memory communication, using the mixed
   shared/distributed memory parallelization is of most use if you run a job on
   a large number of SMP nodes, when the distributed memory parallelization has
   reached its scalability limit (see above). To learn more about the mixed
   parallelization scheme of CPMD consult [mixed]_.

As with all general statements, these are only guidelines. The only way to get
reliable information is to run benchmarks with the system you want to
calculate.
