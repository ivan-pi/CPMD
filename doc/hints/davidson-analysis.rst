.. _hints:pop:

The Davidson analysis and the shared electron number
====================================================

The calculation of the shared electron number can have the following input
section:

::

    &PROPERTIES
      PROJECT WAVEFUNCTION
      POPULATION ANALYSIS MULLIKEN DAVIDSON 2-CENT 3-CENT
        4
        1
      WAVEFUNCTION LATEST
    &END

Note, that for the hydrogen it is enough to specify one atomic orbital to
project on, for the elements Li to Ne it is sufficient to specify 4 atomic
orbitals.
