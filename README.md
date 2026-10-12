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

See [doc/installation.rst](doc/installation.rst) for the build options. The manual
is built with Sphinx from the sources in [doc/](doc/), see below.


## Documentation

The manual is written in reStructuredText in [doc/](doc/) and built with
[Sphinx](https://www.sphinx-doc.org/):

    pip install -r doc/requirements.txt
    sphinx-build -b html doc doc/_build/html     # web pages
    sphinx-build -M latexpdf doc doc/_build      # PDF manual (needs LaTeX)

The GitHub Actions workflow `.github/workflows/docs.yml` publishes the web
pages on GitHub Pages.
