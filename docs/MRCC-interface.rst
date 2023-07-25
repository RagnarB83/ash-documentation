MRCC interface
======================================

MRCC is a very powerful wavefunction-theory QM code, with special support for higher-order coupled cluster methods (CCSDT, CCSDT(Q), CCSDTQ etc.) 
as well as local correlation methods (LNO-CCSD(T)). Compiled binaries are available that makes the program easy to set up.
ASH can be used to drive MRCC calculations in a convenient way by using the MRCCTheory class.
The interface is flexible as you provide a multi-line string containing a MRCC-inputfile.


**MRCCTheory class:**

.. code-block:: python

  class MRCCTheory:
      def __init__(self, mrccdir=None, filename='mrcc', printlevel=2,
                  mrccinput=None, numcores=1, parallelization='OMP-and-MKL'):

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
   * - ``parallelization``
     - string
     - 'OMP-and-MKL'
     - Type of parallelization used. Options: 'OMP', 'OMP-and-MKL' or 'MPI'
   * - ``filename``
     - string
     - 'mrcc'
     - Name of MRCC inputfile



################################
Finding the MRCC program
################################

To use the interface you need to have MRCC-binaries already installed. 
You can then either provide the path to the MRCC binaries using the mrccdir keyword or alternatively ASH will try to find the binaries (searches for the dmrcc binary) automatically in PATH.


################################
Using the interface
################################

You need to provide the mrccinput keyword which should be a multi-line string containing the MRCC syntax. Other keywords are optional.

.. note:: Do **NOT** provide geometry, charge or spin information in the mrccinput string. This information is handled automatically by ASH.


################################
Parallelization
################################

Parallelization of MRCC can be based on OpenMP, BLAS MKL libraries or MPI. 
The MPI parallelization is currently not supported by ASH.
The OpenMP and MKL parallelization is controlled by setting the numcores keyword.

.. note:: **Expert note:** ASH sets ccsdthreads and ptthreads options to numcores. 


################################
Examples
################################


*CCSDT(Q) singlepoint calculation:*

.. code-block:: python

    from ash import *

    #Add coordinates to fragment
    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1)

    #Defining MRCCTheory object
    mrccinput="""
    basis=def2-SVP
    calc=CCSDT(Q)
    mem=9000MB
    scftype=UHF
    ccmaxit=150
    core=frozen
    """
    MRCCcalc = MRCCTheory(mrccinput=mrccinput, numcores=2)
    
    result=Singlepoint(theory=MRCCcalc,fragment=frag)

.. note:: If you are running MRCC on multiple fragments that might be either closed-shell or open-shell it might be convenient to leave out scftype=RHF/scftype=UHF and let MRCC choose this. MRCC may otherwise complain.
