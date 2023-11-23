The Code
==============

The CPMD code is a parallelized **plane wave / pseudopotential** implementation of Density Functional Theory, particularly designed for ab-initio molecular dynamics.


## Copyright Notice

The CPMD program is © 1990-2023 by IBM Corp. and © 1994-2001 by Max Planck Institute, Stuttgart. 


## License

The CPMD is freely distributed under the MIT License.


## Building on SNG Phase2

``` 
cd scratchmodule_lib/; make
cd ..
cd vdw_lib
make
cd .. 
./configure.sh INTEL-XHOST-IFX-MPI-SCRATCHLIBRARY-VDWLIBRARY-OMPOFFLOAD -omp -omp_offload -DEST=build/base_omp; 
cd build/base_omp/
make -j
```
