QUICK interface
======================================

`QUICK <https://quick-docs.readthedocs.io/en/latest/about.html>`_  is an open-source HF/DFT code that runs on the GPU. 
Basis functions are supported only up to d angular momentum.

ASH features a simple interface to QUICK that allows QUICK energy+gradient calculations on the GPU.
The QUICK interface can be used in a QM/MM setting, thus allowing both the QM and MM steps to run on the GPU (QUICK for QM and OpenMM for MM). 

The QUICK interface is a bit limited at the moment.

**QUICKTheory class:**

.. code-block:: python
    
  class QUICKTheory:
      def __init__(self, quickdir=None, filename='quick', printlevel=2,
                  quickinput=None, numcores=1, quickbinary="quick.cuda"):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``quickdir``
     - string
     - None
     - Directory where QUICK binaries are.
   * - ``quickinput``
     - string
     - None
     - Single-line string containing the QUICK input-keywords (without job-keyword).
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores used (only for CPU-version of QUICK).
   * - ``filename``
     - string
     - 'quick'
     - Filename used for QUICK input/output files.
   * - ``quickbinary``
     - string
     - 'quick.cuda'
     - Name of the QUICK binary to use. Options: 'quick.cuda' (GPU-CUDA), 'quick.cuda.MPI' (for multiple CUDA GPUs) or 'quick' (CPU)


################################################################################
QUICK installation
################################################################################

QUICK needs to be compiled and installed separately, paying attention to different GPU architectures available.
See `QUICK installation guide <https://quick-docs.readthedocs.io/en/latest/installation-guide.html>`_.
In our tests we have successfully compiled the CUDA version for a workstation with an NVIDIA Geforce RTX 2080 Ti like below:

.. code-block:: text

  QUICK_HOME=/path/to/quick-install-dir
  mkdir ${QUICK_HOME}/builddir
  cd ${QUICK_HOME}/builddir
  cmake .. -DCOMPILER=GNU -DMPI=TRUE -DCUDA=TRUE -DQUICK_USER_ARCH="volta turing" \
    -DCMAKE_INSTALL_PREFIX=${QUICK_INSTALL}
  make
  make install


################################################################################
QUICK performance control and parallelization
################################################################################

The QUICK binary (quickbinary keyword in QUICKTheory) should generally chosen to be quick.cuda when running on the GPU.
The quick.cuda.mpi binary should be used when multiple GPU cards are available (requires MPI installation).



################################################################################
Examples
################################################################################

**DFT-SCF example:**

.. code-block:: python

  from ash import *

  n2_singlet= Fragment(diatomic="N2", diatomic_bondlength=1.09, charge=0, mult=1)

  #Define QUICK theory : B3LYP/6-31G* calculation (note: do not put job-keyword in string)
  quick = QUICKTheory(quickinput="B3LYP BASIS=6-31G* cutoff=1.0e-8")

  Singlepoint(theory=quick, fragment=n2_singlet)
