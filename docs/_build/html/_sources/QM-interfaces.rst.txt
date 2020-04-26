==========================
QM Interfaces
==========================

Quantum chemistry codes that you can currently use with Yggdrasill are: ORCA, PySCF, Psi4 and xTB.
To use the interfaces you define a QMtheory object of the appropriate Class.
The QM-interface Classes available are: ORCATheory, PySCFTheory, Psi4Theory, xTBTheory

When defining a QMtheory object you are creating an instance of one of the QMTheory classes.
When defining the object, a few keyword arguments are required, that differs between classes.
It is not necessary to define a fragment (Yggdrasill fragment) as part of the QMTheory (e.g. not done for QM/MM),
but is necessary for a single-point energy QM job. For a geometry optimization the fragment does not have to be part
of the QMTheory object (instead passed to the optimizer object/function).
Parallelization of the QM codes differs behind the scenes but is controlled by a nprocs=X keyword for all interfaces.

Example:

Below, we create a dummy QMcalc object of the dummy class QMTheory. We would always set the charge, mult and nprocs keyword (available for all QM theories).
nprocs=1 is the default and the keyword can be skipped if one wants a serial calculation.
Here we also add a fragment (not necessary for an optimization but necessary for a single-point job).
One would also add other keywords that are specific to the QMtheory used (that define the QM method and basis etc.).
We can then run the object  created using the object run function (single-point energy job), with or without a Grad
keyword (request calculation of gradient). It is possible to manipulate some behaviour of the QM object by providing
other arguments to the run-function (here Grad=True and nprocs=12 is provided) but typically this should be kept to a minimum
and should be done when object is created instead.

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    # Defining an object of the (dummy) class QMTheory
    QMcalc = QMTheory(fragment=HF_frag, charge=0, mult=1, nprocs=8)

    #Run a single-point energy job
    QMcalc.run()
    #An Energy+Gradient calculation where we change the number of cores to 12
    QMcalc.run(Grad=True, nprocs=12)


###########################
ORCATheory
###########################
To be done: ORCATheory Class definition and arguments

The ORCA interface is quite flexible. It currently requires the path to the ORCA installation to be passed on as a keyword
argument when creating object. orcasimpleinput and orcablocks keyword arguments (accepts single or multi-line strings) have to be provided
and these keywords define what the ORCA-inputfile looks like. The geometry block would be added to the inputfile by Yggdrasill.
Functionality such as adding to inputfile what orbitals to read and parallelization block would be handled by Yggdrasill as well.


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

    ORCAcalc = ORCATheory(orcadir=orcadir, fragment=HF_frag, charge=0, mult=1,
                                orcasimpleinput=input, orcablocks=blocks, nprocs=8)

    #Run a single-point energy job
    ORCAcalc.run()
    #An Energy+Gradient calculation
    ORCASP.run(Grad=True)



Here a fragment (here called HF_frag) is defined (from XYZ file) and passed to the ORCAtheory object (called ORCAcalc).
Additionally, orcadir, input, and blocks string variables are defined and passed onto the ORCA object via keywords, as
are charge and spin multiplicity.
The object can then be run via the object function (ORCAcalc.run) which is the equivalent of a single-point energy job.
By using Grad=True keyword argument ORCAcalc.run, a gradient is also requested.

**Parallelization**

ORCA parallelization is handled by OpenMPI. By specifying the nprocs=X, a *%pal nprocs X end block* will be added to the
ORCA inputfile created by Yggdrasill. ORCA will then call the OpenMPI mpirun binary when needed and this requires the
correct OpenMPI version to be available.
Make sure the recommended OpenMPI version for the ORCA version you are using is available. This typically requires
setting:

export PATH=/path/to/openmpi/bin:$PATH and export LD_LIBRARY_PATH=/path/to/openmpi/lib:$LD_LIBRARY_PATH

or alternatively loading the appropriate module. Set these variables in the job-script (see :doc:`basics`) that you are using


###########################
Psi4Theory
###########################
The Psi4 interface comes in two versions, a library-based interface and an inputfile-based interface.
The library interface means that Yggdrasill will load Psi4 Python libraries that have to be part of the same Python installation.
In the inputfile-based interface (Psithon), Yggdrasill will create a Psi4 inputfile in Psithon syntax and will then call
a separate Psi4 executable (can be a separate Python installation) via the psi4dir variable (or will find psi4 in shell PATH).

Both interfaces are quite flexible. Most Psi4 settings are controlled by setting the psi4settings dictionary.

Todo:
- Allow to pass dictionaries for other modules
- Enable e.g. coupled-cluster E+G calculation

Polarizable Embedding via Psi4 and the CPPE library is possible (described later).
Set pe=True and give path to potfile to use.

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    #Psi4 variables defined as a dictionary:
    psi4settings={'scf_type': 'pk', 'soscf': True, 'basis' : 'def2-SVP' }
    psi4functional='b3lyp'

    #Psi4: Input-file based interface: using psi4dir to set path
    psi4dir='/path/to/psi4_install/bin/psi4'
    Psi4calc = Psi4Theory(fragment=HF_frag, charge=0, mult=1, psi4settings, psi4functional, runmode='psithon'
                                psi4dir=psi4dir, pe=False, outputname='psi4output.dat', label='psi4input',
                                 psi4memory=3000, prinsetting=False)
    #Psi4: Library-based interface
    Psi4calc = Psi4Theory(fragment=HF_frag, charge=0, mult=1, psi4settings, psi4functional, runmode='library'
                                pe=False, outputname='psi4output.dat', label='psi4input', psi4memory=3000)

    #Run a single-point energy job
    Psi4calc.run()
    #An Energy+Gradient calculation
    Psi4calc.run(Grad=True)

**Parallelization**

The Psi4 parallelization is thread-based. The nprocs keyword provided to the Psi4-interface is used to specify the number
of threads available to Psi4 when the job is run (command-line argument for Psithon and environment variable for library).

###########################
PySCFTheory
###########################
The PySCF interface is library-based and requires a PySCF installation via Pip (pip install pyscf).
At the moment, the interface is not very flexible and only allows for simple DFT calculations with a specific basis set.

Valid keywords are: pyscfbasis, pyscffunctional, fragment, charge, mult, pyscfmemory, nprocs, outputname and printsetting.
Printsetting controls whether to write pyscf-output to a file (False) or to stdout (True).

The interface will become more flexible in the future.

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    #PySCF
    PySCFcalc = PySCFTheory(pyscfbasis="def2-SVP", pyscffunctional="B3LYP", nprocs=2
    fragment=HF_frag, charge=0, mult=1, pyscfmemory=3000, outputname='pyscf.out', printsetting=False)

    #Run a single-point energy job
    PySCFcalc.run()
    #An Energy+Gradient calculation
    PySCFcalc.run(Grad=True)


**Parallelization**

The PySCF parallelization is OpenMP thread-based. The nprocs keyword is used to specify the number of threads available
to PySCF.

###########################
xTBTheory
###########################
The xTB interface comes in two forms, a shared-library interface and a file-based interface.
The shared-library interface is recommended as no disk I/O is required while running xTB. Yggdrasill and xTB communicate via a Python C-API.
As no files are written to disk, this makes the interface faster than the file-based interface, useful for e.g. fast MD.
The file-based interface writes an XYZ-file to disk, calls an xTB executable which reads the XYZ file, runs the job and writes the output to disk which is then read by Yggdrasill.
For regular jobs, e.g. geometry optimizations, the speed-difference between interfaces will probably not matter.

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


**Parallelization**

The xTB parallelization is OpenMP or MKL thread-based and can be controlled via the nprocs keyword.
Currently OMP threads are set equal to nprocs and MKL threads are set equal to 1.
Todo: confirm that this actually works

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