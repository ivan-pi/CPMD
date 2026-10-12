# Building CPMD with CMake

CPMD can be built with CMake (3.18 or newer) as an alternative to
`configure.sh` + `make`. Both build systems use the same sources and the same
machine configurations in `configure/`.

## Quick start

```sh
cmake -B build                      # MPI, LAPACK/BLAS found by CMake
cmake --build build -j 8
./build/bin/cpmd.x input.inp
```

`cmake --install build --prefix /opt/cpmd` installs `bin/cpmd.x` and
`lib/libcpmd.a`.

Without a machine configuration CMake picks the compilers (set `FC`/`CC` or
`CMAKE_Fortran_COMPILER`/`CMAKE_C_COMPILER` to change them), finds MPI,
OpenMP and LAPACK with `find_package()`, and adds the preprocessor
definitions from the options below. The build type defaults to `Release`.
`-DBLA_VENDOR=OpenBLAS` (or `Intel10_64lp`, ...) selects a BLAS/LAPACK
vendor.

## Using a machine configuration from `configure/`

The machine configurations used by `configure.sh` can be used unchanged:

```sh
cmake -B build -DCPMD_CONFIG=LINUX-X86_64-GFORTRAN-MPI
cmake -B build -DCPMD_CONFIG=CRAY-HPE-GFORTRAN-AMDEPYC -DCPMD_ENABLE_OPENMP=ON
cmake -B build -DCPMD_CONFIG=/path/to/MY-MACHINE        # a file outside configure/
```

The configuration can also be given through the environment, which is
convenient in job scripts or module files:

```sh
export CPMD_CONFIG=LINUX-X86_64-GFORTRAN-MPI
cmake -B build
```

CMake evaluates the configuration file with
`scripts/cmake/import_config.sh`, exactly as `configure.sh` does, so the
`if [ $omp ]`, `if [ $debug ]` and `if [ $coverage ]` branches inside the file
keep working. The switches of `configure.sh` correspond to:

| configure.sh  | CMake                                   |
|---------------|-----------------------------------------|
| `-omp`        | `-DCPMD_ENABLE_OPENMP=ON`               |
| `-debug`      | `-DCMAKE_BUILD_TYPE=Debug`              |
| `-coverage`   | `-DCPMD_ENABLE_COVERAGE=ON`             |
| `-qmmm`       | `-DCPMD_ENABLE_QMMM=ON`                 |
| `-iphigenie`  | `-DCPMD_ENABLE_IPHIGENIE=ON`            |
| `-scr`        | `-DCPMD_ENABLE_SCRATCHLIB=ON`           |
| `-vdw`        | `-DCPMD_ENABLE_VDWLIB=ON`               |
| `-minpack=L`  | `-DCPMD_MINPACK_LIBRARY=L`              |
| `-DEST=dir`   | `-B dir`                                |

To see what CMake will take from a configuration, run the helper by hand:

```sh
scripts/cmake/import_config.sh -omp LINUX-X86_64-GFORTRAN-MPI
```

The variables of a configuration file are used as follows:

| Variable                  | Use in the CMake build                                    |
|---------------------------|-----------------------------------------------------------|
| `FC`, `CC`, `NVCC`        | compilers; extra words are treated as flags. They can only be set in a fresh build directory and are ignored if `CMAKE_<LANG>_COMPILER` is given |
| `FFLAGS`                  | Fortran compile flags, also passed when linking `cpmd.x`  |
| `CFLAGS`, `NVCCFLAGS`     | C and CUDA compile flags                                  |
| `CPPFLAGS`                | `-D` become compile definitions, `-I` include directories |
| `LFLAGS` (`LIBS`)         | link libraries and flags of `cpmd.x`                      |
| `FFLAGS_GROMOS[_MODULES]` | extra flags for the QM/MM sources                         |
| `NOOPT_OBJS`, `NOOPT_FLAG`| `NOOPT_FLAG` is appended to the listed sources            |
| `OMP3_DISABLED`           | adds `-D__HASNT_OMP_COLLAPSE` with OpenMP                 |
| `CPP`, `LD`, `AR`, `RANLIB`, `IRAT` | not used                                        |

Notes:

* CMake always lets the Fortran compiler preprocess the sources
  (`CMAKE_Fortran_PREPROCESS`), an external `CPP` is not used.
* With a configuration, CMake does not add optimisation flags of its own
  unless `CMAKE_BUILD_TYPE` is set; `Debug` selects the `$debug` branch of the
  configuration file.
* References to make or environment variables in a configuration, such as
  `$(MKLROOT)`, are expanded from the environment.
* Changing `CPMD_CONFIG` in an existing build directory updates the flags,
  but not the compilers. Use a new build directory for a different compiler.

## Options

| Option                     | Default   | Meaning                                                   |
|----------------------------|-----------|-----------------------------------------------------------|
| `CPMD_CONFIG`              | `$CPMD_CONFIG` | machine configuration (name in `configure/` or path) |
| `CPMD_ENABLE_MPI`          | `ON`      | parallel version (`-D__PARALLEL`); taken from `CPPFLAGS` with a configuration |
| `CPMD_ENABLE_OPENMP`       | `OFF`     | OpenMP                                                    |
| `CPMD_FFT`                 | `DEFAULT` | `DEFAULT` (built-in) or `FFTW3`; without a configuration  |
| `CPMD_ENABLE_LIBXC`        | `OFF`     | libxc functionals (`-D__HAS_LIBXC`); without a configuration |
| `CPMD_ENABLE_CUDA`         | `OFF`     | CUDA code (`-D__HAS_CUDA`); without a configuration, else from `NVCC`/`CPPFLAGS` |
| `CPMD_ENABLE_QMMM`         | `OFF`     | QM/MM, needs the Gromos sources in `modules/Gromos`       |
| `CPMD_ENABLE_IPHIGENIE`    | `OFF`     | builds `libiffiinter.a` for IPHIGENIE                     |
| `CPMD_ENABLE_SCRATCHLIB`   | `OFF`     | scratch library from `modules/scratchmodule_lib`          |
| `CPMD_ENABLE_VDWLIB`       | `OFF`     | Grimme vdW library from `modules/vdw_lib`                 |
| `CPMD_ENABLE_COVERAGE`     | `OFF`     | coverage instrumentation (GNU)                            |
| `CPMD_MINPACK_LIBRARY`     |           | MINPACK library (`-D__MINPACK`)                           |
| `CPMD_BUILD_TESTS`         | `ON`      | add the smoke tests in `tests/`                           |
| `CPMD_TEST_NPROCS`         | `2`       | MPI processes used by the tests                           |

## Tests

`ctest --test-dir build` runs a few short calculations (`tests/`) and
compares the final total energy with a reference value. MPI builds run them
with `mpiexec -n 2`; set `MPIEXEC_EXECUTABLE`, `MPIEXEC_PREFLAGS` or
`CPMD_TEST_NPROCS` to change that. These are smoke tests for the build, not a
replacement for the CPMD regression test suite.

The GitHub Actions workflow `.github/workflows/cmake.yml` builds CPMD with
the gfortran versions of the Ubuntu runner (native configuration with and
without MPI, and one `configure/` file) and runs the tests.

## Layout of the build

| Target              | File                    | Contents                                   |
|---------------------|-------------------------|--------------------------------------------|
| `cpmd`              | `bin/cpmd.x`            | `cpmd.F90`, `timetag.F90`                  |
| `cpmd_lib`          | `lib/libcpmd.a`         | all other sources in `src/`                |
| `cpmd_gromos`       | `lib/libgromos.a`       | `modules/Gromos` (QM/MM)                   |
| `cpmd_mm_interface` | `lib/libinterface.a`    | `modules/MM_Interface` (QM/MM)             |
| `cpmd_iffiinter`    | `lib/libiffiinter.a`    | `modules/IPhigenie_Interface`              |
| `cpmd_scratch`      | `lib/libscratch_module.a` | `modules/scratchmodule_lib`              |
| `cpmd_grimme`       | `lib/libgrimme.a`       | `modules/vdw_lib`                          |

The source files are listed in `src/CMakeLists.txt` and `modules/CMakeLists.txt`.
In `src/CMakeLists.txt` they are grouped by topic, following the `SRC_*`
groups of `src/SOURCES` (wavefunction optimisation, geometry optimisation,
molecular dynamics, metadynamics, ...). New files must be added to one of
these groups; as long as `configure.sh` is maintained they must also be
added to `src/SOURCES`. CMake warns when the two lists differ.

The compilation order of the Fortran modules is determined by CMake, no
dependency list has to be maintained.
