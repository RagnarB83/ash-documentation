PySCF interface
======================================


PySCFTheory class:

.. code-block:: python
    
    class PySCFTheory:
        def __init__(self, printsetting=False, printlevel=2, pyscfbasis='', pyscffunctional='',
                    pe=False, potfile='', filename='pyscf', pyscfmemory=3100, numcores=1):



**PySCFTheory** options:

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``pyscfbasis``
     - string
     - ''
     - Name of PySCF basis set
   * - ``pyscffunctional``
     - string
     - ''
     - Name of PySCF DFT functional
   * - ``pe``
     - Boolean
     - False
     - Whether to use polarizable embedding in PySCF via CPPE library.
   * - ``potfile``
     - string
     - ''
     - Name of potential file for in PySCF CPPE polarizable embedding
   * - ``filename``
     - string
     - 'pyscf'
     - Filename used for PySCF output
   * - ``pyscfmemory``
     - integer
     - 3100
     - Memory (in MB) used by PySCF .
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores used by PySCF
   * - ``printsetting``
     - string
     - None
     - Printsetting. if True: printing to standard output otherwise disk.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel

The PySCF interface is library-based and requires a PySCF installation via Pip (pip install pyscf).
At the moment, the interface is not very flexible and only allows for simple DFT calculations with a specific basis set.

Valid keywords are: pyscfbasis, pyscffunctional, fragment, charge, mult, pyscfmemory, numcores, outputname and printsetting.
Printsetting controls whether to write pyscf-output to a file (False) or to stdout (True).

The interface will become more flexible in the future.

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz', charge=0, mult=1)
    #PySCF
    PySCFcalc = PySCFTheory(pyscfbasis="def2-SVP", pyscffunctional="B3LYP", numcores=2,
    pyscfmemory=3000, outputname='pyscf.out', printsetting=False)

    #Run a single-point energy job
    Singlepoint(theory=PySCFcalc, fragment=HF_frag)
    #An Energy+Gradient calculation
    Singlepoint(theory=PySCFcalc, fragment=HF_frag, Grad=True)



**Parallelization**

The PySCF parallelization is OpenMP thread-based. The numcores keyword is used to specify the number of threads available
to PySCF.