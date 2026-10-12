The Code
==============

The CPMD code is a parallelized **plane wave / pseudopotential** implementation of Density Functional Theory, particularly designed for ab-initio molecular dynamics.


## Copyright Notice

The CPMD program is © 1990-2023 by IBM Corp. and © 1994-2001 by Max Planck Institute, Stuttgart. 


## License

The CPMD is freely distributed under the MIT License.


## Building

CPMD can be built either with the classic configuration script

    ./configure.sh [-omp] [-qmmm] ... CONFIGURATION
    make

or with CMake, which can use the same machine configurations from `configure/`:

    cmake -B build [-DCPMD_CONFIG=CONFIGURATION] [-DCPMD_ENABLE_OPENMP=ON]
    cmake --build build

See [doc/cmake.md](doc/cmake.md) for the CMake options.
