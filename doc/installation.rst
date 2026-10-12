.. _installation:

Installation
============

CPMD is written in Fortran (with a few C and, optionally, CUDA sources). It
can be built in two ways that use the same sources and the same machine
configuration files in the ``configure/`` directory:

* with **CMake** (3.18 or newer), described in :ref:`install-cmake`, or
* with the classic **configure.sh** script, which writes a ``Makefile`` for
  ``make``, described in :ref:`install-configure`.

Prerequisites
-------------

* A recent Fortran compiler (e.g. GNU gfortran, Intel
  ifx/ifort, NVIDIA/PGI, Cray, Fujitsu, ...) and a C compiler.
* BLAS and LAPACK (e.g. OpenBLAS, Intel MKL, ESSL).
* An MPI library with the ``mpi_f08`` module for parallel runs. Without MPI
  (no ``-D__PARALLEL``) a serial version is built.
* Optional: FFTW3 (otherwise the built-in FFT is used), `libxc
  <https://libxc.gitlab.io/>`_ for additional exchange-correlation
  functionals, CUDA for the GPU code, MINPACK, and the Gromos sources for the
  QM/MM interface (see :ref:`install-qmmm`).

The machine configuration files in ``configure/`` (for example
``LINUX-X86_64-GFORTRAN-MPI``, ``LINUX-X86_64-INTEL-MPI``,
``CRAY-HPE-GFORTRAN-AMDEPYC``) set compilers, flags, preprocessor
definitions and libraries for a given platform. Most likely a configuration
will **not** match your machine exactly and you have to adapt the compiler
names, optimization flags and library locations; copy the closest file and
edit it, or pass variables on the command line as shown below.

.. _install-cmake:

Building with CMake
-------------------

Quick start::

   cmake -B build                      # MPI, LAPACK/BLAS found by CMake
   cmake --build build -j 8
   ./build/bin/cpmd.x input.inp

``cmake --install build --prefix /opt/cpmd`` installs ``bin/cpmd.x`` and
``lib/libcpmd.a``.

Without a machine configuration CMake picks the compilers (set ``FC``/``CC``
or ``CMAKE_Fortran_COMPILER``/``CMAKE_C_COMPILER`` to change them), finds MPI,
OpenMP and LAPACK with ``find_package()``, and adds the preprocessor
definitions from the options below. The build type defaults to ``Release``.
``-DBLA_VENDOR=OpenBLAS`` (or ``Intel10_64lp``, ...) selects a BLAS/LAPACK
vendor.

Using a machine configuration from ``configure/``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The machine configurations used by ``configure.sh`` can be used unchanged::

   cmake -B build -DCPMD_CONFIG=LINUX-X86_64-GFORTRAN-MPI
   cmake -B build -DCPMD_CONFIG=CRAY-HPE-GFORTRAN-AMDEPYC -DCPMD_ENABLE_OPENMP=ON
   cmake -B build -DCPMD_CONFIG=/path/to/MY-MACHINE        # a file outside configure/

The configuration can also be given through the environment, which is
convenient in job scripts or module files::

   export CPMD_CONFIG=LINUX-X86_64-GFORTRAN-MPI
   cmake -B build

CMake evaluates the configuration file with
``scripts/cmake/import_config.sh``, exactly as ``configure.sh`` does, so the
``if [ $omp ]``, ``if [ $debug ]`` and ``if [ $coverage ]`` branches inside
the file keep working. The switches of ``configure.sh`` correspond to:

=================  ===================================
``configure.sh``   CMake
=================  ===================================
``-omp``           ``-DCPMD_ENABLE_OPENMP=ON``
``-debug``         ``-DCMAKE_BUILD_TYPE=Debug``
``-coverage``      ``-DCPMD_ENABLE_COVERAGE=ON``
``-qmmm``          ``-DCPMD_ENABLE_QMMM=ON``
``-iphigenie``     ``-DCPMD_ENABLE_IPHIGENIE=ON``
``-scr``           ``-DCPMD_ENABLE_SCRATCHLIB=ON``
``-vdw``           ``-DCPMD_ENABLE_VDWLIB=ON``
``-minpack=L``     ``-DCPMD_MINPACK_LIBRARY=L``
``-DEST=dir``      ``-B dir``
=================  ===================================

To see what CMake will take from a configuration, run the helper by hand::

   scripts/cmake/import_config.sh -omp LINUX-X86_64-GFORTRAN-MPI

The variables of a configuration file are used as follows:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Use in the CMake build
   * - ``FC``, ``CC``, ``NVCC``
     - Compilers; extra words are treated as flags. They can only be set in a
       fresh build directory and are ignored if ``CMAKE_<LANG>_COMPILER`` is
       given.
   * - ``FFLAGS``
     - Fortran compile flags, also passed when linking ``cpmd.x``.
   * - ``CFLAGS``, ``NVCCFLAGS``
     - C and CUDA compile flags.
   * - ``CPPFLAGS``
     - ``-D`` become compile definitions, ``-I`` include directories.
   * - ``LFLAGS`` (``LIBS``)
     - Link libraries and flags of ``cpmd.x``.
   * - ``FFLAGS_GROMOS[_MODULES]``
     - Extra flags for the QM/MM sources.
   * - ``NOOPT_OBJS``, ``NOOPT_FLAG``
     - ``NOOPT_FLAG`` is appended to the listed sources.
   * - ``OMP3_DISABLED``
     - Adds ``-D__HASNT_OMP_COLLAPSE`` with OpenMP.
   * - ``CPP``, ``LD``, ``AR``, ``RANLIB``, ``IRAT``
     - Not used.

Notes:

* CMake always lets the Fortran compiler preprocess the sources
  (``CMAKE_Fortran_PREPROCESS``); an external ``CPP`` is not used.
* With a configuration, CMake does not add optimization flags of its own
  unless ``CMAKE_BUILD_TYPE`` is set; ``Debug`` selects the ``$debug`` branch
  of the configuration file.
* References to make or environment variables in a configuration, such as
  ``$(MKLROOT)``, are expanded from the environment.
* Changing ``CPMD_CONFIG`` in an existing build directory updates the flags,
  but not the compilers. Use a new build directory for a different compiler.

CMake options
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Default
     - Meaning
   * - ``CPMD_CONFIG``
     - ``$CPMD_CONFIG``
     - Machine configuration (name in ``configure/`` or path).
   * - ``CPMD_ENABLE_MPI``
     - ``ON``
     - Parallel version (``-D__PARALLEL``); taken from ``CPPFLAGS`` with a
       configuration. With ``OFF`` CPMD is built without MPI.
   * - ``CPMD_ENABLE_OPENMP``
     - ``OFF``
     - OpenMP.
   * - ``CPMD_FFT``
     - ``DEFAULT``
     - ``DEFAULT`` (built-in, ``-D__HAS_FFT_DEFAULT``) or ``FFTW3``
       (``-D__HAS_FFT_FFTW3``); without a configuration. See
       ``src/cpmd_global.h`` for how the FFT library is selected.
   * - ``CPMD_ENABLE_LIBXC``
     - ``OFF``
     - libxc functionals (``-D__HAS_LIBXC``); without a configuration.
   * - ``CPMD_ENABLE_CUDA``
     - ``OFF``
     - CUDA code (``-D__HAS_CUDA``); without a configuration, else taken from
       ``NVCC``/``CPPFLAGS``.
   * - ``CPMD_ENABLE_QMMM``
     - ``OFF``
     - QM/MM, needs the Gromos sources in ``modules/Gromos``.
   * - ``CPMD_ENABLE_IPHIGENIE``
     - ``OFF``
     - Builds ``libiffiinter.a`` for IPHIGENIE.
   * - ``CPMD_ENABLE_SCRATCHLIB``
     - ``OFF``
     - Scratch library from ``modules/scratchmodule_lib``.
   * - ``CPMD_ENABLE_VDWLIB``
     - ``OFF``
     - Grimme vdW library from ``modules/vdw_lib``.
   * - ``CPMD_ENABLE_COVERAGE``
     - ``OFF``
     - Coverage instrumentation (GNU).
   * - ``CPMD_MINPACK_LIBRARY``
     -
     - MINPACK library (``-D__MINPACK``).
   * - ``CPMD_BUILD_TESTS``
     - ``ON``
     - Add the tests in ``tests/`` and ``.github/tests/``.
   * - ``CPMD_TEST_NPROCS``
     - ``2``
     - MPI processes used by the tests.
   * - ``CPMD_TEST_NTHREADS``
     - ``1``
     - OpenMP threads per process used by the tests.

Tests
~~~~~

``ctest --test-dir build`` runs a few short calculations and compares the
last total energy of each with a reference value: a smoke test in ``tests/``
and the regression tests used by CI in ``.github/tests/`` (registered in
``.github/tests/tests.cmake``, pseudopotentials in ``.github/tests/pseudo/``).
They cover wavefunction optimisation (ODIIS, PCG, LSD, PBE, k-points, CP
groups), restarting (also with MPI-IO), geometry optimisation, BO and CP
molecular dynamics and, with ``CPMD_ENABLE_VDWLIB``, the Grimme D3
correction. MPI builds run them with ``mpiexec -n 2``; set
``MPIEXEC_EXECUTABLE``, ``MPIEXEC_PREFLAGS`` or ``CPMD_TEST_NPROCS`` to change
that, and ``CPMD_TEST_NTHREADS`` for the number of OpenMP threads per process
(default 1). These tests check that the build works and that results do not
change; they are not a replacement for the CPMD regression test suite
(https://github.com/CPMD-code/Regtests).

FFTW3 builds take much longer to run the tests (about a minute each) because
CPMD creates its FFTW plans with ``FFTW_EXHAUSTIVE``.

The GitHub Actions workflow ``.github/workflows/cmake.yml`` builds CPMD with
the gfortran versions of the Ubuntu runner (native configuration with MPI,
with FFTW3, with the scratch and vdW libraries and two OpenMP threads,
without MPI, and one ``configure/`` file) and runs the tests. A separate job
checks that MPI is only used through the ``mpi_f08`` module.

Layout of the build
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Target
     - File
     - Contents
   * - ``cpmd``
     - ``bin/cpmd.x``
     - ``cpmd.F90``, ``timetag.F90``
   * - ``cpmd_lib``
     - ``lib/libcpmd.a``
     - all other sources in ``src/``
   * - ``cpmd_gromos``
     - ``lib/libgromos.a``
     - ``modules/Gromos`` (QM/MM)
   * - ``cpmd_mm_interface``
     - ``lib/libinterface.a``
     - ``modules/MM_Interface`` (QM/MM)
   * - ``cpmd_iffiinter``
     - ``lib/libiffiinter.a``
     - ``modules/IPhigenie_Interface``
   * - ``cpmd_scratch``
     - ``lib/libscratch_module.a``
     - ``modules/scratchmodule_lib``
   * - ``cpmd_grimme``
     - ``lib/libgrimme.a``
     - ``modules/vdw_lib``

The source files are listed in ``src/CMakeLists.txt`` and
``modules/CMakeLists.txt``. In ``src/CMakeLists.txt`` they are grouped by
topic, following the ``SRC_*`` groups of ``src/SOURCES`` (wavefunction
optimization, geometry optimization, molecular dynamics, metadynamics, ...).
New files must be added to one of these groups; as long as ``configure.sh``
is maintained they must also be added to ``src/SOURCES``. CMake warns when
the two lists differ. The compilation order of the Fortran modules is
determined by CMake; no dependency list has to be maintained.

.. _install-configure:

Building with ``configure.sh`` and make
---------------------------------------

Run ``configure.sh`` from the top-level directory of the source tree. Without
arguments it lists the available machine configurations::

   ./configure.sh

Choose the configuration closest to your machine and run the script again
with its name. By default the ``Makefile`` and the ``obj/``, ``lib/`` and
``bin/`` directories are created in the top-level directory; with ``-DEST``
the build goes into a separate directory, which is recommended when you
build several variants (for example serial and OpenMP, or different
compilers) from the same sources::

   ./configure.sh -DEST=build/gfortran-mpi LINUX-X86_64-GFORTRAN-MPI
   cd build/gfortran-mpi
   make -j 8

The executable is ``bin/cpmd.x`` inside the build directory. Additional
information about a configuration, and the full list of options, are printed
with::

   ./configure.sh -info LINUX-X86_64-GFORTRAN-MPI
   ./configure.sh -help

The options of ``configure.sh`` are:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Meaning
   * - ``-help``, ``-h``
     - Print the help message.
   * - ``-info``, ``-i``
     - Print additional information about a machine configuration.
   * - ``-debug``, ``-d``
     - Use the debugging branch of the configuration file.
   * - ``-verbose``, ``-v``
     - Display the names of files for the dependencies.
   * - ``-coverage``, ``-c``
     - GNU compilers only: build with code coverage / profiling
       instrumentation (if the configuration file supports it).
   * - ``-omp``, ``-o``
     - Enable OpenMP (if the configuration file allows it). Whether OpenMP 3
       ``collapse`` clauses are used is controlled by ``OMP3_DISABLED`` in the
       configuration file.
   * - ``-disable_omp3``
     - Disable OpenMP 3 ``collapse`` clauses regardless of the configuration.
   * - ``-qmmm``, ``-q``
     - Build the Gromos QM/MM interface (needs the Gromos sources, see
       :ref:`install-qmmm`).
   * - ``-iphigenie``
     - Build the external interface to IPHIGENIE (see
       :doc:`hints/iphigenie`).
   * - ``-scr``
     - Use the scratch library (``modules/scratchmodule_lib``).
   * - ``-vdw``
     - Compile and include the Grimme van der Waals library
       (``modules/vdw_lib``).
   * - ``-minpack=<PATH>``
     - Enable MINPACK and link with the given library.
   * - ``-VAR=value``
     - Override a variable of the generated ``Makefile``, e.g.
       ``-DEST=dir`` (build directory), ``-FC=...``, ``-FFLAGS=...``,
       ``-LFLAGS=...``, ``-CFLAGS=...``, ``-CPPFLAGS=...``, ``-CPP=...``,
       ``-LD=...``, ``-AR=...``, ``-RANLIB=...``.

If you change the preprocessor flags (the ``CPPFLAGS`` entry), e.g. going from
a serial to a parallel compilation, delete all preprocessed, module and
object files first::

   make clean

To create a new configuration, add a new file to the ``configure/``
directory, using one of the existing files as a template.

.. _install-qmmm:

QM/MM interfaces
----------------

The Gromos-based QM/MM interface (:ref:`sec:qmmm`) needs the Gromos sources,
which are proprietary (ETH Zurich) and cannot be distributed with CPMD under
the MIT license; ``modules/Gromos`` therefore only contains a README. If you
have a license for the Gromos sources, place them in ``modules/Gromos`` and
build with ``-qmmm`` (``-DCPMD_ENABLE_QMMM=ON``). The interface to IPHIGENIE
(:doc:`hints/iphigenie`) is distributed with CPMD.

Building this manual
--------------------

The manual is written in reStructuredText and built with `Sphinx
<https://www.sphinx-doc.org/>`_::

   pip install -r doc/requirements.txt
   sphinx-build -b html doc doc/_build/html          # web pages
   sphinx-build -M latexpdf doc doc/_build           # PDF (needs LaTeX)

The GitHub Actions workflow ``.github/workflows/docs.yml`` builds both and
publishes the web pages on GitHub Pages.
