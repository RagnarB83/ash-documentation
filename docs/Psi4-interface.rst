Psi4 interface
======================================

`Psi4 <https://psicode.org>`_  is a full-featured open-source quantum chemistry program written in both Python and C++.
ASH has an interface both to the Psithon inputfile-based interface of Psi4 but also to the Psi4 Python library.

**Psi4Theory class:**

.. code-block:: python
    
    class Psi4Theory:
        def __init__(self, psi4dir=None, runmode='library', printsetting='False', psi4settings=None, 
                    psi4method=None, pe=False, potfile='', filename='psi4_', label='psi4input',
                    psi4memory=3000, numcores=1, printlevel=2, fchkwrite=False):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``psi4dir``
     - string
     - None
     - Path to Psi4 directory. Used by runmode='psithon'
   * - ``runmode``
     - string
     - 'library'
     - Psi4 runmode. Options: runmode='library' or runmode='psithon'
   * - ``printsetting``
     - Boolean
     - False
     - Printsetting. If True: Psi4 output is printed to stdout. If False, printing to disk.
   * - ``psi4settings``
     - dict
     - None
     - Psi4 settings dictionary
   * - ``psi4method``
     - string
     - None
     - Psi4 method 
   * - ``pe``
     - Boolean
     - False
     - Whether to use polarizable embedding in Psi4 via CPPE library.
   * - ``potfile``
     - string
     - ''
     - Name of potential file for in Psi4 CPPE polarizable embedding
   * - ``filename``
     - string
     - 'psi4'
     - Name of Psi4 inputfile for psithon runmode.
   * - ``label``
     - string
     - None
     - Optional label of object
   * - ``psi4memory``
     - integer
     - 3000
     - Memory in MB used by Psi4
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores Psi4 will use
   * - ``printlevel``
     - integer
     - 2
     - The printlevel.
   * - ``fchkwrite``
     - Boolean
     - False
     - XXXX

The Psi4 interface comes in two versions, a library-based interface and an inputfile-based interface.
The library interface means that ASH will load Psi4 Python libraries that have to be part of the same Python installation.
In the inputfile-based interface (Psithon), ASH will create a Psi4 inputfile in Psithon syntax and will then call
a separate Psi4 executable (can be a separate Python installation) via the psi4dir variable (or will find psi4 in shell PATH).

Both interfaces are quite flexible. Most Psi4 settings are controlled by setting the psi4settings dictionary.
The Psi4 method is controlled by the psi4method argument. It can be set to a functional name : e.g. 'B3LYP', 'B3LYP-D3BJ'
or even 'CCSD'  and 'CCSD(T)'. Analytical gradients are available in Psi4 for both CCSD and CCSD(T).

Todo:
- Allow to pass dictionaries for other modules

Polarizable Embedding via Psi4 and the CPPE library is possible (described later).
Set pe=True and give path to potfile to use.

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz', charge=0, mult=1)
    #Psi4 variables defined as a dictionary:
    psi4settings={'scf_type': 'pk', 'soscf': True, 'basis' : 'def2-SVP' }
    psi4method='b3lyp'

    #Psi4: Input-file based interface: using psi4dir to set path
    psi4dir='/path/to/psi4_install/bin/psi4'
    Psi4calc = Psi4Theory(psi4settings=psi4settings, psi4method=psi4method, runmode='psithon',
                                psi4dir=psi4dir, pe=False, outputname='psi4output.dat', label='psi4input',
                                 psi4memory=3000, prinsetting=False)
    #Psi4: Library-based interface
    Psi4calc = Psi4Theory(psi4settings=psi4settings, psi4method=psi4method, runmode='library',
                                pe=False, outputname='psi4output.dat', label='psi4input', psi4memory=3000)

    #Run a single-point energy job
    Singlepoint(theory=Psi4calc, fragment=HF_frag)
    #An Energy+Gradient calculation
    Singlepoint(theory=Psi4calc, fragment=HF_frag, Grad=True)

**Parallelization**

The Psi4 parallelization is thread-based. The numcores keyword provided to the Psi4-interface is used to specify the number
of threads available to Psi4 when the job is run (command-line argument for Psithon and environment variable for library).