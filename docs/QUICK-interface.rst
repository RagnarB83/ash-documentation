QUICK interface
======================================

`QUICK <https://quick-docs.readthedocs.io/en/latest/about.html>`_  is an open-source HF/DFT code that runs on the GPU. 
Basis functions only up to d angular momentum are supported.

ASH features a simple interface to it that allows QUICK energy+gradient calculations on the GPU.
The QUICK interface be used in a QM/MM setting allowing both the QM and MM steps to run on the GPU (QUICK for QM and OpenMM for MM). 
Interface is a bit limited at the moment.

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
     - Name of the QUICK binary to use. Options: 'quick.cuda' (GPU-CUDA) or 'quick' (CPU)



QUICK needs to be compiled and installed separately.

**DFT-SCF example:**

.. code-block:: python

  from ash import *

  n2_singlet= Fragment(diatomic="N2", diatomic_bondlength=1.09, charge=0, mult=1)

  #Define QUICK theory : B3LYP/6-31G* calculation (note: do not put job-keyword in string)
  quick = QUICKTheory(quickinput="B3LYP BASIS=6-31G* cutoff=1.0e-8")

  Singlepoint(theory=quick, fragment=n2_singlet)
