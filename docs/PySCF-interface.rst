PySCF interface
======================================

PySCF is a very powerful open-source quantum chemistry code. ASH currently features a limited interface that will be extended
in the future.

**PySCFTheory class:**

.. code-block:: python
    
  class PySCFTheory:
      def __init__(self, printsetting=False, printlevel=2, numcores=1, 
                  scf_type=None, basis=None, functional=None, gridlevel=5, 
                  pe=False, potfile='', filename='pyscf', memory=3100, conv_tol=1e-8, verbose_setting=4, 
                  CC=False, CCmethod=None, CC_direct=False, frozen_core_setting='Auto', 
                  frozen_virtuals=None, FNO=False, FNO_thresh=None):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``printsetting``
     - string
     - None
     - Printsetting. if True: printing to standard output otherwise disk.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores used by PySCF
   * - ``scf_type``
     - string
     - None
     - Type of SCF-determinant. Options: 'RHF','UHF','RKS','UKS'.
   * - ``basis``
     - string
     - None
     - Name of basis set (must be valid PySCF basis-name).
   * - ``functional``
     - string
     - None
     - Name of DFT functional (must be valid PySCF functional-name).
   * - ``gridlevel``
     - string
     - ''
     - Name of DFT functional (must be valid PySCF functional-name).
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
   * - ``memory``
     - integer
     - 3100
     - Memory (in MB) used by PySCF .
   * - ``conv_tol``
     - float
     - 1e-8
     - Convergence tolerance in Eh .
   * - ``verbose_setting``
     - int
     - 4
     - How verbose PySCF output is.
   * - ``CC``
     - Boolean
     - False
     - Whether to do coupled-cluster on top of SCF or not.
   * - ``CCmethod``
     - string
     - None
     - Type of CCSD-method. Options:'CCSD', 'CCSD(T)'. More options will be available.
   * - ``CC_direct``
     - Boolean
     - False
     - Whether to use integral-direct CC or not.
   * - ``frozen_core_setting``
     - string
     - 'Auto'
     - How frozen core is handled. The ASH-default option is 'Auto' which means that frozen core settings are chosen by ASH (mimics ORCA-settings).
   * - ``frozen_virtuals``
     - list
     - None
     - Optionally freeze selected virtual orbitals in CC calculation.
   * - ``FNO``
     - Boolean
     - False
     - Do frozen natural orbital coupled cluster using MP2 natural orbitals.
   * - ``FNO_thresh``
     - float
     - None
     - Optional threshold to choose virtual natural orbitals to be skipped, based on natural occupation (from MP2 occupations).


The PySCF interface is library-based and requires a PySCF installation via Pip (pip install pyscf).
The interface is currently not very flexible and only supports SCF and coupled cluster calculations at the moment.
The interface will be extended in the future.


**DFT-SCF example:**

.. code-block:: python

  from ash import *

  n2_singlet= Fragment(diatomic="N2", diatomic_bondlength=1.09, charge=0, mult=1)

  #Define PySCF theory: RKS-PBE0 hybrid-DFT calculation
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RKS', functional="PBE0", gridlevel=6,
    numcores=2, memory=3000, filename='pyscf.out', printsetting=False)

  Singlepoint(theory=PySCFcalc, fragment=n2_singlet)



**Unrestricted CCSD(T) example:**

.. code-block:: python

  from ash import *

  o2_triplet= Fragment(diatomic="O2", diatomic_bondlength=1.2075, charge=0, mult=3)

  #PySCF with UHF SCF and CCSD(T) on top
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", numcores=2, scf_type="UHF", CC=True,
    CCmethod='CCSD(T)', memory=3000, filename='pyscf.out', printsetting=False)


  Singlepoint(theory=PySCFcalc, fragment=o2_triplet)




**Parallelization**

The PySCF parallelization is OpenMP thread-based. The numcores keyword is used to specify the number of threads available
to PySCF.