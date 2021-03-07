==========================
QM Interfaces
==========================

Quantum chemistry codes that you can currently use with Ash are: ORCA, PySCF, Psi4 and xTB.
To use the interfaces you define a QMtheory object of the appropriate Class.
The QM-interface Classes available are: ORCATheory, PySCFTheory, Psi4Theory, xTBTheory

When defining a QMtheory object you are creating an instance of one of the QMTheory classes.
When defining the object, a few keyword arguments are required, that differs between classes.
Parallelization of the QM codes differs behind the scenes but is controlled by a nprocs=X keyword for all interfaces.

Example:

Below, we create a dummy QMcalc object of the dummy class QMTheory. We would always set the charge, mult and nprocs keyword (available for all QM theories).
nprocs=1 is the default and the keyword can be skipped if one wants a serial calculation.
One would also add other keywords that are specific to the QMtheory used (that define the QM method and basis etc.).
We can then run a single-point calculation of a fragment using the object.
This is done using the Singlepoint function that requires both theory and fragment keyword arguments.
Additionally an Energy+Gradient calculation can be requested.

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    # Defining an object of the (dummy) class QMTheory
    QMcalc = QMTheory(charge=0, mult=1, nprocs=8)

    #Run a single-point energy job
    Singlepoint(theory=QMcalc, fragment=HF_frag)
    #An Energy+Gradient singlepoint calculation.
    Singlepoint(theory=QMcalc, fragment=HF_frag, Grad=True)

    #Cleanup
    QMcalc.cleanup()

Since some QM-interfaces may generate a lot of files it can be a good idea to request a cleanup after the job is run
as sometimes old files may interfer with new calculations.
This would be done via the internal cleanup function of the object as shown.


###########################
ORCATheory
###########################

See :doc:`ORCA-interface`



###########################
Psi4Theory
###########################

See :doc:`Psi4-interface`


###########################
PySCFTheory
###########################

See :doc:`PySCF-interface`



###########################
xTBTheory
###########################

See :doc:`xTB-interface`


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
