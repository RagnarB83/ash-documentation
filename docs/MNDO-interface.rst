MNDO interface
======================================

`MNDO <https://mndo.kofo.mpg.de>`_  is a semiempirical quantum chemistry program written by Walter Thiel at the MPI für Kohlenforschung in Mülheim, Germany.

ASH features a basic interface to MNDO.
Energies and gradients are available in the interface so MNDOTheory in ASH can be used for single-point energies, geometry optimizations, 
numerical frequencies, surface scans, NEB and molecular dynamics within ASH. 
Full QM/MM pointcharge-support is available so the interface can be used fully in QM/MM jobs.

MNDO requires a valid license. See `MNDO licensing <https://mndo.kofo.mpg.de/license.php>`_  .
See also `MNDO input documentation <https://mndo.kofo.mpg.de/input.php>`_ .

**MNDOTheory class:**

.. code-block:: python
    
    class MNDOTheory:
        def __init__(self, mndodir=None, filename='mndo', method=None, printlevel=2, label="MNDO",
                    numcores=1, restart_option=True, diis=False, guess_option=0):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``mndodir``
     - string
     - None
     - Directory where MNDO binaries are. ASH will search for the 'mndo2020' binary in this dir.
   * - ``filename``
     - string
     - 'mndo'
     - Filename used for MNDO input/output files.
   * - ``mndo``
     - string
     - None
     - Name of method to use. Options: 'MNDO', 'AM1', 'PM3', 'OM2', 'OM3', 'ODM2', 'ODM3'
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``label``
     - string
     - 'MNDO'
     - Optional label for MNDOTheory object.
   * - ``numcores``
     - integer
     - 1
     - Number of cores. MNDO parallelization is currently not active.
   * - ``restart_option``
     - Boolean
     - True
     - Whether MNDO will read/write density-matrix to disk for restarts. Can affect speed.
   * - ``guess_option``
     - integer
     - 0
     - What type of guess_option to use. Controls ktrial keyword in MNDO. See `MNDO input documentation <https://mndo.kofo.mpg.de/input.php>`_
   * - ``diis``
     - Boolean
     - False
     - Whether DIIS is on or off. Default mode does not use DIIS unless convergence problem.
   * - ``numcores``
     - integer
     - 1
     - Number of cores. MNDO parallelization is currently not active.

################################################################################
MNDO installation
################################################################################

MNDO requires a valid license. See `MNDO licensing <https://mndo.kofo.mpg.de/license.php>`_ 

Once MNDO has been downloaded, to use it with ASH you need to either :
i) provide the *mndodir* keyword option to MNDOTHeory to provide the location of the directory where the mndo2020 binary resides.
or 
ii) Make sure the mndo2020 binary is in the PATH variable of your shell environment. ASH will try to find an executable the *mndo2020* binary automatically when creating the MNDOTheory object.

################################################################################
Parallelization
################################################################################

MNDO by default does not run in parallel. There are some MPI parallelization options but they are not yet compatible with ASH.
The *numcores* keyword in MNDOTheory is not used at the moment.

################################################################################
Examples
################################################################################

**MNDO-AM1 example:**

.. code-block:: python

    from ash import *

    #Variable defining number of cores
    numcores=1

    #Fragment
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)

    #Define MNDOTheory object
    mndo_object = MNDOTheory(method="AM1")

    Singlepoint(theory=mndo_object, fragment=frag)
