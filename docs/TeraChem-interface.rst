TeraChem interface
======================================

`TeraChem <http://www.petachem.com/products.html>`_  is a quantum chemistry code written exclusively for the Nvidia GPU.
The program is commercial, however, the demo version can be used for up to 15 minutes per session.

ASH features a simple interface to it that allows TeraChem energy+gradient calculations on the GPU.
Can be used in a QM/MM setting allowing both the QM and MM steps to run on the GPU via TeraChem and OpenMM. 
Interface is a bit limited at the moment.

**TeraChemTheory class:**

.. code-block:: python
    
  class TeraChemTheory:
      def __init__(self, terachemdir=None, filename='terachem', printlevel=2,
                  teracheminput=None, numcores=1):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``terachemdir``
     - string
     - None
     - Directory where TeraChem binaries are.
   * - ``teracheminput``
     - Python dict
     - None
     - Python dictionary containing string key-value pairs of TeraChem input options.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``filename``
     - string
     - 'pyscf'
     - Filename used for TeraChem input/output files.


TeraChem needs to be installed separately.

**DFT-SCF example:**

.. code-block:: python

    from ash import *
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)

    #Python dictionary with TeraChem input options. See manual: http://www.petachem.com/doc/userguide.pdf
    teracheminput_dict={
    'method':'ublyp',
    'basis':'6-31g*',
    'sphericalbasis':'false',
    'scf':'diis+a',
    'timings':'yes'
    }

    terachem = TeraChemTheory(teracheminput=teracheminput_dict)

    Singlepoint(theory=terachem, fragment=frag)

**GPU control:**

How Terachem runs on the GPU/GPUs is controlled outside of ASH via environment variables.
The most important variable to control is what GPU device to use (1 or multiple).

.. code-block:: text

    #Use GPU device 0 and 1
    export CUDA_VISIBLE_DEVICES=0,1

