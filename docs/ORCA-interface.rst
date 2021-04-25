ORCA interface
======================================


To be done: ORCATheory Class definition and arguments



The ORCA interface is quite flexible. It currently requires the path to the ORCA installation to be passed on as a keyword
argument when creating object. orcasimpleinput and orcablocks keyword arguments (accepts single or multi-line strings) have to be provided
and these keywords define what the ORCA-inputfile looks like. The geometry block would be added to the inputfile by Ash.
Functionality such as adding to inputfile what orbitals to read and parallelization block would be handled by Ash as well.


.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    #ORCA
    orcadir='/opt/orca_4.2.1'
    input="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    blocks="""
    %scf
    maxiter 200
    end
    %basis
    newgto F "ma-def2-SVP" end
    end
    """

    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput=input, orcablocks=blocks, nprocs=8)

    #Run a single-point energy job
    Singlepoint(theory=ORCAcalc, fragment=HF_frag)
    #An Energy+Gradient calculation
    Singlepoint(theory=ORCAcalc, fragment=HF_frag, Grad=True)



Here a fragment (here called HF_frag) is defined (from XYZ file) and passed to the Singlepoint function along with an
ORCAtheory object (called ORCAcalc). The orcadir, input, and blocks string variables are defined and passed onto the ORCA object via keywords, as
are charge and spin multiplicity. By default, the ORCA autostart feature is active, meaning that if an inputfile with name "orca-input.inp" is run, ORCA will
try to read orbitals from "orca-input.gbw" file if present. This is utilized automatically during geometry optimizations, numerical frequencies as well
as multiple single-point calculations sequentially. It is possible to turn this off by adding "!Noautostart" in the simple-inputline of the orcasimpleinput variable.
It is also possible to have each ORCA-calculation read in orbitals from another source by using the: moreadfile keyword argument option:

.. code-block:: python

    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput=input,
                        orcablocks=blocks, nprocs=8, moreadfile="orbitals.gbw")


Note: For parallel-ASH calculations (ASH in parallel, ORCA in serial). The full path to the moreadfile may be required.


The ORCA object is then used by passing it to a function: e.g. Singlepoint, an optimizer, a QM/MM object, NumFreq function etc.
When the ORCA object is run (e.g. by the Singlepoint function, an optimizer etc.) it will create an ORCA inputfile
that will always be called orca-input.inp. This inputfile will look familiar to any ORCA user as it will contain a "Simpleinput line", Block-input
a coordinate block etc. (cordinates in Å). ASH will then tell ORCA to run the inputfile and an outputfile called orca-input.out will be created.
Once the ORCA calculation is done the outputfile (or other files) is read for information (usually the energy and gradient) by ASH
and ASH will continue. The ORCA inputfile , "orca-input.inp" may be replaced later (e.g. if an optimization job" and ORCA
will be run again.





################################################################################
Parallelization
################################################################################

ORCA parallelization is handled by OpenMPI. By specifying the nprocs=X, a *%pal nprocs X end block* will be added to the
ORCA inputfile created by Ash. ORCA will then call the OpenMPI mpirun binary when needed and this requires the
correct OpenMPI version to be available.
Make sure the recommended OpenMPI version for the ORCA version you are using is available. This typically requires
setting:

export PATH=/path/to/openmpi/bin:$PATH and export LD_LIBRARY_PATH=/path/to/openmpi/lib:$LD_LIBRARY_PATH

or alternatively loading the appropriate module. Set these variables in the job-script (see :doc:`basics`) that you are using

