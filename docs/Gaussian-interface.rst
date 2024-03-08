Gaussian interface
======================================

`Gaussian <https://gaussian.com`_  is one of the oldest and widely used quantum chemistry programs.

ASH features a basic interface to Gaussian, allowing you to use it for DFT and WFT based methods within the program.
Energies and gradients are available in the interface so GaussianTheory in ASH can be used for single-point energies, geometry optimizations, 
numerical frequencies, surface scans, NEB and molecular dynamics. 
Warnign: Pointcharge-support is available, however, as the pointcharge gradient is not available, 
electrostatic embedding QM/MM is limited to single-point energy calculations, for now.

**GaussianTheory class:**

.. code-block:: python
    
    class GaussianTheory:
        def __init__(self, gaussiandir=None, gauss_executable='g16', filename='gaussian', 
                    file_extension='.com', printlevel=2, label="Gaussian",
                    gaussianinput=None, memory='800MB', numcores=1):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``gaussiandir``
     - string
     - None
     - Directory where Gaussian binaries are.
   * - ``gauss_executable``
     - string
     - 'g16'
     - Name of the Gaussian executable to be used (e.g. 'g16', 'g09', 'g03').
   * - ``filename``
     - string
     - 'gaussian'
     - Filename used for Gaussian input/output files.
   * - ``file_extension``
     - string
     - '.com'
     - File extension used for Gaussian input/output files created by ASH, e.g. '.com' or '.gjf'.
   * - ``gaussianinput``
     - string
     - None
     - Multi-line string containing the Gaussian input-lines.
   * - ``memory``
     - string
     - '800MB'
     - A string specifying the memory to be given to Gaussian.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``numcores``
     - integer
     - 1
     - The number of CPU cores used for OpenMP parallelization of Gaussian.
   * - ``label``
     - string
     - 'Gaussian'
     - Optional label for GaussianTheory object.

################################################################################
Gaussian installation
################################################################################

Gaussian needs to be installed separately.
As we do not have access to the Gaussian source code, we cannot offer any installation instructions.

Once Gaussian has been installed you need to either :
i) use the gaussiandir keyword option to GaussianTheory to provide the location of the directory where the Gaussian binary resides 
or 
ii) Make sure the g16/g09/g03 binary is in the PATH variable of your shell environment. ASH will try to find an executable the *gauss_executable* (e.g. 'g16') automatically when creating the GaussianTheory object.

################################################################################
Parallelization
################################################################################

Gaussian can be parallelized using the built-in OpenMP. This parallelization will only work for single node parallelization.
The parallelization is controlled via the *numcores* keyword.

Note that you should not add any parallelization options to the Gaussian input string, as ASH will handle this.


################################################################################
Examples
################################################################################

**DFT-SCF example:**

.. code-block:: python

    from ash import *

    #Variable defining number of cores
    numcores=4

    #Fragment
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)

    #Multi-line inputstring to define the main Gaussian options
    #No geometry-input, job-keywords (opt, freq etc.) or parallelization options should be provided here (handled by ASH)
    input_string=""""
    # b3lyp cc-pVDZ scf=tight int=ultrafine
    """"

    gaussian_object = GaussianTheory(gauss_executable='g09', gaussianinput=input_string, numcores=numcores)

    Singlepoint(theory=gaussian_object, fragment=frag)
