MRCC interface
======================================

**MRCCTheory class:**

.. code-block:: python

    class MRCCTheory:
        def __init__(self, mrccdir=None, filename='mrcc', printlevel=2,
                    mrccinput=None, numcores=1):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``mrccdir``
     - string
     - None
     - Path to MRCC directory.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``mrccinput``
     - string
     - None
     - MRCC input as a multi-line string 
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores MRCC will use
   * - ``filename``
     - string
     - 'mrcc'
     - Name of MRCC inputfile



**Example:**


.. code-block:: python

    from ash import *

    #Add coordinates to fragment
    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1)

    #Defining MRCCTheory object
    mrccinput="""
    basis=def2-SVP
    calc=CCSDT
    mem=9000MB
    scftype=UHF
    ccmaxit=150
    core=frozen
    """
    MRCCcalc = MRCCTheory(mrccdir="/path/to/mrccdir", mrccinput=mrccinput, numcores=1)
    
    result=Singlepoint(theory=MRCCcalc,fragment=frag)