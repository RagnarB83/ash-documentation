PySCF interface
======================================

The PySCF interface is library-based and requires a PySCF installation via Pip (pip install pyscf).
At the moment, the interface is not very flexible and only allows for simple DFT calculations with a specific basis set.

Valid keywords are: pyscfbasis, pyscffunctional, fragment, charge, mult, pyscfmemory, numcores, outputname and printsetting.
Printsetting controls whether to write pyscf-output to a file (False) or to stdout (True).

The interface will become more flexible in the future.

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    #PySCF
    PySCFcalc = PySCFTheory(pyscfbasis="def2-SVP", pyscffunctional="B3LYP", numcores=2,
    charge=0, mult=1, pyscfmemory=3000, outputname='pyscf.out', printsetting=False)

    #Run a single-point energy job
    Singlepoint(theory=PySCFcalc, fragment=HF_frag)
    #An Energy+Gradient calculation
    Singlepoint(theory=PySCFcalc, fragment=HF_frag, Grad=True)



**Parallelization**

The PySCF parallelization is OpenMP thread-based. The numcores keyword is used to specify the number of threads available
to PySCF.