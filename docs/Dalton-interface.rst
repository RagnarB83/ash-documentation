Dalton interface
======================================

**DaltonTheory class:**

.. code-block:: python

    class DaltonTheory:
        def __init__(self, mrccdir=None, filename='mrcc', printlevel=2,
                    mrccinput=None, numcores=1):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``daltondir``
     - string
     - None
     - Path to Dalton directory.
   * - ``dalton_input``
     - string
     - None
     - Dalton input as a multi-line string 
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores Dalton will use
   * - ``filename``
     - string
     - 'dalton'
     - Name of Dalton inputfile
   * - ``basis_name``
     - string
     - None
     - Name of Dalton basis set.
   * - ``basis_dir``
     - string
     - None
     - Name of basis set directory
   * - ``pe``
     - Boolean
     - False
     - Whether to use polarizable embedding in Dalton.
   * - ``potfile``
     - string
     - ''
     - Name of potential file for polarizable embedding option.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``label``
     - string
     - None
     - Optional label for object

**Example:**

.. code-block:: python

    from ash import *

    #Add coordinates to fragment
    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1)

    #Defining DaltonTheory object
    input_string="""
    **DALTON
    .RUN RESPONSE
    .DIRECT
    **WAVE FUNCTIONS
    .DFT
    CAMB3LYP
    **RESPONSE
    .TDA
    *LINEAR
    .SINGLE RESIDUE
    .ROOTS
    5
    **END OF
    """
    Daltonobject = DaltonTheory(daltondir="/path/to/daltondir",dalton_input=input_string)
    
    #Running single-point job
    result=Singlepoint(theory=Daltonobject,fragment=frag)