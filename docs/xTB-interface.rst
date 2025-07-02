xTB interface
======================================

`xTB <https://xtb-docs.readthedocs.io>`_  is a semiempirical tightbinding DFT quantum chemistry code from the Grimme group.
It features the extended tight-binding (xTB) methods (GFN0-xTB,GFN1-xTB and GFN2-xTB) that feature a fairly accurate 
electrostatic interaction term and built-in dispersion correction and are parameterized in a general element parameterization 
fashion that avoids pair-potentails, unlike standard DFTB-based methods.
GFN-xTB methods give overall fairly accurate geometries, frequencies and noncovalent interactions, 
often even for transition metal complex while energies are further from regular DFT. 
The speed of these methods is 100-1000 x compared to regular DFT.

More recently, a new method, g-xTB, has been described that gives much more accurate energies, much closer to regular DFT:
`g-xTB preprint <https://chemrxiv.org/engage/chemrxiv/article-details/685434533ba0887c335fc974>`_

**gxTBTheory class:**

The g-xTB method is available as a preliminary implementation in the gxtb binary, 
see `g-xTB Github repository <https://github.com/grimme-lab/g-xtb>`_ .
As this implementation features only a numerical gradient, 
geometry optimizations will be slow and will suffer from some numerical noise.
A future implementation is expected in the tblite library that ASH will support once available.

ASH features a very basic interface to the gxtb binary that allows for energies and slow geometry optimizations.

.. code-block:: python

  class gxTBTheory(Theory):
      def __init__(self, printlevel=2, numcores=1):

No QM/MM is supported in gxTBTheory yet as pointcharge-support is not available in the gxtb program. For xTB-based QM/MM, see xTBTheory class below.


**xTBTheory class:**

*xTBTheory* is the interface to the general xTB program that supports the GFN-xTB methods (GFN0, GFN1 and GFN2) with analytic
gradients, implicit solvation and pointcharges.
xTBTheory can be used for QM/MM and ONIOM embedding etc.

.. code-block:: python

    class xTBTheory:
        def __init__(self, xtbdir=None, xtbmethod='GFN1', runmode='inputfile', numcores=1, printlevel=2, filename='xtb_',
                    maxiter=500, electronic_temp=300, label=None, accuracy=0.1, hardness_PC=1000, solvent=None):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``xtbdir``
     - string
     - None
     - Path to xTB directory. Optional. See "Finding the xTB program"
   * - ``xtbmethod``
     - string
     - 'GFN1'
     - The xTB Hamiltonian to use. Options: 'GFN2', 'GFN1', 'GFN0', 'GFNFF', 'IPEA'
   * - ``runmode``
     - string
     - inputfile
     - Whether to run xTB using 'inputfile' mode (read/write files on disk) or using 'library' (share data in memory).
   * - ``numcores``
     - integer
     - 1
     - The number of cores used in xTB parallelization (thread-based).
   * - ``printlevel``
     - integer
     - 2
     - The level of printing desired.
   * - ``filename``
     - string
     - xtb\_
     - The basename used to create the xTB XYZ-file read by xTB.
   * - ``maxiter``
     - integer
     - 500
     - Max number of SCC iterations
   * - ``electronic_temp``
     - integer
     - 300
     - The electronic temperature used by xTB in fractional occupation smearing.
   * - ``label``
     - string
     - None
     - Optional label to define for xTB object.
   * - ``accuracy``
     - integer
     - 0.1
     - The xtB accuracy parameter in the SCC calculation. See `xtb-docs <https://xtb-docs.readthedocs.io/en/latest/sp.html?highlight=accuracy#accuracy-and-iterations>`_ 
   * - ``hardness_PC``
     - integer
     - 1000
     - What hardness to use for pointcharges in electrostatic QM/MM embedding. 1000 is recommended.
   * - ``solvent``
     - string
     - None
     - Name of solvent to use with the internal ALPB solvent model. Only available for inputfile interface. See `xtbdocs <https://xtb-docs.readthedocs.io/en/latest/gbsa.html#implicit-solvation>`_ 

The xTB interface comes in two forms, a shared-library interface and a file-based interface.

The shared-library interface is generally recommended (limitation: no QM/MM yet) as no disk I/O is performed while running xTB. ASH and xTB communicate via a Python C-API.
It requires first installation of xtb-python (https://xtb-python.readthedocs.io/en/latest/installation.html), e.g. via conda: conda install -c conda-forge xtb-python

The file-based interface writes an XYZ-file to disk, calls an xTB executable which reads the XYZ file, runs the job and writes the output to disk which is then read by ASH.
For regular jobs, e.g. geometry optimizations, the speed-difference between interfaces will probably not matter.

################################
Finding the xTB program
################################

ASH can find the xTB program in a few different ways. For the inputfile-based runmode:

- ASH will first check if the xtbdir argument has been set which should be a string that points to the directory where the xtb program is located, e.g. "xtbdir=/path/to/xtb_6_4_0/bin". This option takes precedence.
- If the xtbdir argument has not been provided ASH will next see if xtbdir has been provided in the ASH settings: :doc:`basics`
- If xtbdir has also not been defined there, ASH will search the PATH environment variable for an executable "xtb" and if found will set the xtbdir accordingly. This can be a convenient option if you make sure to define your shell environments carefully in your jobscript or shell-startup file. Be careful, however, if you have multiple version of the program available.

runmode='library' on the other hand requires the installation of the xtb-python interface to your Python environment. ASH will check and complain if it does not find the library. 

################################
Examples
################################

To use either interface is quite simple, when an xTB object is created, the xtbmethod keyword is used to select xTB method: "GFN2", "GFN1" for the GFN2-xTB and GFN1-xTB Hamiltonians, respectively.
The optional runmode argument is also available: runmode='library' or runmode='inputfile'. Default runmode: "inputfile"


.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz', charge=0, mult=1)
    xTBcalc = xTBTheory(xtbmethod='GFN2', runmode='library')

    #Run a single-point energy job on the fragment associated with the xtb-object
    Singlepoint(theory=xTBcalc, fragment=HF_frag)
    #An Energy+Gradient calculation running on 8 cores
    Singlepoint(theory=xTBcalc, fragment=HF_frag, Grad=True)



################################
Parallelization
################################
The xTB parallelization is OpenMP or MKL thread-based and can be controlled via the numcores keyword.
Currently OMP threads are set equal to numcores and MKL threads are set equal to 1.