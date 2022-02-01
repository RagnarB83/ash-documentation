==========================
QM Interfaces
==========================

Quantum chemistry codes that you can currently use with ASH are: ORCA, xTB, PySCF, Psi4, Dalton, CFour, MRCC.
To use the interfaces you define a QMtheory object of the appropriate Class.
The QM-interface Classes available are: ORCATheory, xTBTheory, PySCFTheory, Psi4Theory, DaltonTheory, CFourTheory, MRCCTheory.

When defining a QMtheory object you are creating an instance of one of the QMTheory classes.
When defining the object, a few keyword arguments are required, that differs between classes.
Parallelization of the QM codes differs behind the scenes but is controlled by a numcores=X keyword for all interfaces.

Example:

Below, we create a dummy QMcalc object of the dummy class QMTheory. 
The arguments that QMTheory takes will depend on the interface but typically we would at least set the numcores keyword (available for all QM theories) that defines how many CPU cores the QM program is allowed to use.
numcores=1 is the default and the keyword can be skipped if one wants a serial calculation.
One would then add other keywords that are specific to the QMtheory used that will define the QM method and basis etc.).
Note that it is a design choice in ASH to not define general variables such as functional, basis set (that could in theory be used for all interfaces) as this would heavily restrict the flexibility of the interface.
Instead, each interface differs in how the electronic structure details are defined, with the aim of allowing you to select any method you prefer to use in the respective QM code.

Once both a fragment object and a theory object has been created we can run a basic single-point energy calculation.
This is done using the Singlepoint function that requires both theory and fragment keyword arguments. Additionally charge and multiplicity information needs to be provided.

Charge and multiplicity information is either provided to the fragment object:

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz', charge=0, mult=1)
    # Defining an object of the (dummy) class QMTheory
    QMcalc = QMTheory(numcores=8)

    #Run a single-point energy job
    Singlepoint(theory=QMcalc, fragment=HF_frag)

or it can be provided to the jobtype function (here Singlepoint):

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    # Defining an object of the (dummy) class QMTheory
    QMcalc = QMTheory(numcores=8)

    #Run a single-point energy job
    Singlepoint(theory=QMcalc, fragment=HF_frag, charge=0, mult=1)



When running an ASH script involving multiple fragments (that differ w.r.t. charge/mult), it is best to associate charge/mult attributes to the fragment.
The charge/mult information from the fragment will then be provided to the theory object when run via the respective job-function (Singlepoint, geomeTRICOptimizer etc.) instead.

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

**Methods**

- run(self, current_coords=None, charge=None, mult=None, current_MM_coords=None, MMcharges=None, qm_elems=None, elems=None, Grad=False, Hessian=False, PC=False, numcores=None, label=None).

- cleanup(self)

Each QMTheory class has a run method that will be called by a jobtype function (e.g. Singlepoint or geomeTRICOptimizer) and the current coordinates will be provided.

The cleanup method removes temporary files created by the QM-program (or ASH) that may interfer with the next calculation.


###########################
ORCATheory
###########################

See :doc:`ORCA-interface`

###########################
xTBTheory
###########################

See :doc:`xTB-interface`


###########################
Psi4Theory
###########################

See :doc:`Psi4-interface`


###########################
PySCFTheory
###########################

See :doc:`PySCF-interface`


###########################
DaltonTheory
###########################

See :doc:`Dalton-interface`


###########################
CFourTheory
###########################

See :doc:`CFour-interface`


###########################
MRCCTheory
###########################

See :doc:`MRCC-interface`

