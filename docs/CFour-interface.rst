CFour interface
======================================

**CFourTheory class:**

.. code-block:: python

    class CFourTheory:
        def __init__(self, cfourdir=None, printlevel=2, cfouroptions=None, numcores=1,
                    filename='cfourjob',specialbasis=None, ash_basisfile='def2-SVP'):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``cfourdir``
     - string
     - None
     - Path to CFour directory.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``cfouroptions``
     - dict
     - None
     - CFour keywords as a Python dictionary 
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores CFour will use
   * - ``filename``
     - string
     - 'cfourjob'
     - Name of CFour inputfile
   * - ``specialbasis``
     - dictionary
     - None
     - Optional specialbasis option.
   * - ``ash_basisfile``
     - string
     - 'def2-SVP'
     - ASH-internal basis set file for CFour. Options: 'def2-SVP'



**Example:**

.. code-block:: python

    from ash import *

    #Add coordinates to fragment
    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1)

    cfouroptions = {
    'CALC':'CCSD',
    'BASIS':'PVTZ',
    'REF':'RHF',
    'FROZEN_CORE':'ON',
    'MEM_UNIT':'MB',
    'MEMORY':3100,
    'PROP':'FIRST_ORDER',
    'CC_PROG':'ECC',
    'SCF_CONV':10,
    'LINEQ_CONV':10,
    'CC_MAXCYC':300,
    'SYMMETRY':'OFF',
    'HFSTABILITY':'OFF'
    }

    cfourcalc = CFourTheory(cfourdir='/path/to/cfour', cfouroptions=cfouroptions)

    #Simple Energy SP calc
    energy = Singlepoint(theory=cfourcalc, fragment=HF_frag)