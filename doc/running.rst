.. _runcpmd:

Running CPMD
============

The CPMD program is started with the following command:

::

   # cpmd.x {file.in} [{PP_path}] > {file.out}

Running ``cpmd.x`` requires the following files:

-  an **input file**  *file.in* (see section :ref:`inputfile`)

-  **pseudopotential files** for all atomic species specified in *file.in* (see
   section :ref:`S_Pseudopotentials`).

The path to the pseudopotential library can be given in different ways:

-  The second command line argument  [*PP_path*]  is set.

-  | If the second command line argument  [*PP_path*]  is not given, the
     program checks the environment variables
   |  **CPMD_PP_LIBRARY_PATH and PP_LIBRARY_PATH**.

-  If neither the environment variables nor the second command line argument
   are set, the program assumes that the pseudopotential files are in the
   current directory

During the run cpmd.x creates different outputs:

-  Various status messages to monitor the correct operation of the program is
   printed to standard output (in our case redirected to the file *file.out*).

-  Detailed data is written into different files (depending on the keywords
   specified in the input *file.in*). An overview on them is given in section
   :ref:`FILES`. Large files are written either to the current
   directory, the directory specified by the environment variable
   **CPMD_FILEPATH**, or the directory specified in the input file using the
   keyword :kw:`FILEPATH`.

In case CPMD quits with a non-zero exit status, an error message is written to
a file called *LocalError-X-X-X.log* by every task, indicating the procedure in
which the error occurred and the call stack. In certain cases (mostly when
processing the input), more ample error information is written to the output
unit as well.

Jobs can be stopped at the next breakpoint by creating a file:

.. container:: center

   *EXIT*

in the run-directory.
