==========================
QM Interfaces
==========================

QM interfaces currently supported in Yggdrasill are:
ORCA, PySCF, Psi4 and xTB.
To use the interfaces you define a QMtheory object of the appropriate Class.
The QM-interface Classes are:
ORCATheory, PySCFTheory, Psi4Theory, xTBTheory

When defining a QMtheory object you are creating an instance of one of the QMTheory classes.
When defining the object, a few keyword arguments are required, that differs between classes.
It is not necessary to define a fragment (Yggdrasill fragment) as part of the QMTheory (e.g. not done for QM/MM),
but is necessary for a single-point energy QM job. For a geometry optimization the fragment does not have to be part
of the QMTheory object (instead passed to the optimizer object/function).

###########################
ORCATheory
###########################
To be done: ORCATheory Class definition and arguments

The ORCA interface is quite flexible. It currently requires the path to the ORCA installation to be passed on as a keyword
argument when creating object. The charge, multiplicity are also necessary keyword arguments (integers).
As are orcasimpleinput and orcablocks keyword arguments (accepts single or mult-line strings).


.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    #ORCA
    orcadir='/opt/orca_4.2.1'
    input="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    blocks="%scf maxiter 200 end"
    ORCAcalc = ORCATheory(orcadir=orcadir, fragment=HF_frag, charge=0, mult=1,
                                orcasimpleinput=input, orcablocks=blocks)

    #Run a single-point energy job
    ORCAcalc.run()
    #An Energy+Gradient calculation running on 8 cores
    ORCASP.run(Grad=True, nprocs=8)



Here a fragment (here called HF_frag) is defined (from XYZ file) and passed to the ORCAtheory object (called ORCAcalc).
Additionally, orcadir, input, and blocks string variables are defined and passed onto the ORCA object via keywords, as
are charge and spin multiplicity.
The object can then be run via the object function (ORCAcalc.run) which is the equivalent of a single-point energy job.
By using Grad=True keyword argument ORCAcalc.run, a gradient is also requested and by setting the nprocs=8 argument,
an 8-core run ORCA calculation is requested (handled via OpenMPI, requires OpenMPI variables to be set outside Python).

###########################
Psi4Theory
###########################


###########################
xTBTheory
###########################
The xTB interface comes in two forms, a shared-library interface and a file-based interface.
The shared-library interface is recommended as no disk I/O is required while running xTB. Yggdrasill and xTB communicate via a Python C-API.
As no files are written to disk, this makes the interface faster than the file-based interface, useful for e.g. fast MD.
The file-based interface writes an XYZ-file to disk, calls an xTB executable which reads the XYZ file, runs the job and writes the output to disk which is then read by Yggdrasill.
For regular jobs, e.g. geometry optimizations, the speed-difference between interface will probably not matter.

To use either interface is quite simple, when an xTB object is created, charge and multiplicity keywords should be provided
as well as the xtbmethod keyword argument that takes values: "GFN2", "GFN1" for the GFN2-xTB and GFN1-xTB Hamiltonians, respectively.
An optional fragment object can also be associated with the xTB-object (makes only sense for single-point jobs).
An optional runmode argument is also available: runmode='library' or runmode='inputfile'.

The runmode='library' option is used by default and requires the shell environment variable LD_LIBRARY_PATH to include the xtb library dir.
e.g. export LD_LIBRARY_PATH=/path/to/your/xtb_6_2_3/lib64:$LD_LIBRARY_PATH

The runmode='inputfile' option requires an additional xtbdir variable to be set that points to the dir containing the xtb executable, e.g. xtbdir=/path/to/xtb_6_2_3/bin .

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    xTBcalc = xTBTheory(fragment=HF_frag, charge=0, mult=1, xtbmethod='GFN2')
    #xTBcalc = xTBTheory(fragment=HF_frag, charge=0, mult=1, xtbmethod='GFN2', runmode='inputfile', xtbdir='/path/to/xtb_6_2_3/bin')

    #Run a single-point energy job on the fragment associated with the xtb-object
    xTBcalc.run()
    #An Energy+Gradient calculation running on 8 cores
    xTBcalc.run(Grad=True, nprocs=8)




Troubleshooting:
==================

- If the library-interface is not working, the reason is likely that something is missing from the LD_LIBRARY_PATH environment variable,  Make sure the lib64 dir of xtb is part of the LD_LIBRARY_PATH in the shell from which you are running (or in the jobscript you are submitting).
e.g. export LD_LIBRARY_PATH=/path/to/your/xtb_6_2_3/lib64:$LD_LIBRARY_PATH

- Fortran libraries may also be missing for xTB. Make sure to load the necessary libraries (e.g. loading a module or  sourcing the Intel compilervars.sh script)


- If the problem is not resolved, try to load the Yggdrasill xtb-interface directly in a script:

.. code-block:: python

    import xtb_interface_library
    test = xtb_interface_library.XTBLibrary()

That should reveal what libraries are not found.