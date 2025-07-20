QM Interfaces
==========================

Quantum chemistry codes that you can currently use with ASH are: ORCA, xTB, PySCF, Psi4, Dalton, CFour, ccpy, MRCC, NWChem, TeraChem, QUICK, CP2K.
Additionally codes like Dice and Block2 have interfaces in ASH, which requires PySCF installed as well.

To use the interfaces you define a QMtheory object of the appropriate Class.
The QM-interface Classes are: 

- ORCATheory (:doc:`ORCA-interface`)
- xTBTheory (:doc:`xTB-interface`)
- PySCFTheory (:doc:`PySCF-interface`)
- Psi4Theory (:doc:`Psi4-interface`)
- DaltonTheory (:doc:`Dalton-interface`)
- CFourTheory (:doc:`CFour-interface`)
- MRCCTheory (:doc:`MRCC-interface`)
- QUICKTheory (:doc:`QUICK-interface`)
- NWChemTheory (:doc:`NWChem-interface`)
- TeraChemTheory (:doc:`TeraChem-interface`)
- CP2KTheory (:doc:`CP2K-interface`)
- DiceTheory (:doc:`Dice-interface`)
- BlockTheory (:doc:`Block-interface`)
- ccpyTheory (:doc:`ccpy-interface`)


When defining a QMtheory object you are creating an instance of one of the QMTheory classes.
When defining the object, a few keyword arguments are required, that differs between classes.
Parallelization of the QM codes also differs behind the scenes but is all controlled by ASH with the numcores=X keyword for all interfaces.


**Example**

Below, we create a dummy QMcalc object of the dummy class **QMTheory**. 
The arguments that **QMTheory** takes will depend on the interface but typically we would at least set the numcores keyword (available for all QM theories) that defines how many CPU cores the QM program is allowed to use.
numcores=1 is always the default and the keyword can be skipped if one wants a serial calculation.
One would then add other keywords that are specific to the **QMtheory** used that will define the QM method and basis etc.).
Note that it is a design choice in ASH to not define general variables such as functional, basis set (that could in theory be used for many interfaces) as this would heavily restrict the flexibility of the interface.
Instead, each interface differs in how the electronic structure details are defined, with the aim of allowing you to select any method you prefer to use in the respective QM code and with the options in the program you like to use.

Once both a **Fragment** object and a **Theory** object has been created we can run a job, e.g. a single-point energy calculation or a geometry optimization.
This is done by calling a job-function (here **Singlepoint**) that usually requires both theory and fragment keyword arguments.


.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz', charge=0, mult=1)
    # Defining an object of the (dummy) class QMTheory
    QMcalc = QMTheory(numcores=8)

    #Run a single-point energy job
    Singlepoint(theory=QMcalc, fragment=HF_frag)



#############################################################
Attributes and methods available to all QM interfaces
#############################################################

**Attributes**

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``printlevel``
     - integer
     - 2
     - The level of printing to use when QMTheory is defined or run.
   * - ``numcores``
     - integer
     - 1
     - The number of CPU cores that the QM program will use (parallelization may be MPI or thread-based).
   * - ``label``
     - string
     - None
     - A string-label that can be useful to distinguish different QMTheory objects.
   * - ``filename``
     - string
     - None
     - A string that may be used to name inputfiles for the QMTheory.


**Methods**

- run(self, current_coords=None, charge=None, mult=None, current_MM_coords=None, MMcharges=None, qm_elems=None, elems=None, Grad=False, Hessian=False, PC=False, numcores=None, label=None).

- cleanup(self)

Each QMTheory class has a run method that will be called by a jobtype function (e.g. Singlepoint or geomeTRICOptimizer) and the current coordinates will be provided.
However, it is recommended to instead use the job-function **Singlepoint** for running a simple energy or energy+gradient calculation.

The cleanup method removes temporary files created by the QM-program (or ASH) that may interfer with the next calculation.


