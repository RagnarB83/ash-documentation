Fairchem interface
======================================

`Fairchem <https://fair-chem.github.io>`_ is a library containing the machine-learning models from the META Fair project.
The `UMA models <https://arxiv.org/abs/2506.23971>`_ in particular, trained on the large Open Molecules 2025 dataset (`OMol25 <https://arxiv.org/abs/2505.08762>`_ ) 

ASH features a simple interface to the Fairchem library that allows easy use of the models.

Once a FairchemTheory has been created, specifying a model or previously downloaded file, it
can be used to calculate energies and gradients (just like a regular QM or MM theory).
A FairchemTheory object can be used for single-point energies, geometry optimizations, numerical frequencies, surface scans, NEB, molecular dynamics etc. within ASH.
Even hybrid ONIOM and QM/MM calculations are possible (with some limitations).

Periodic boundary conditions are also supported.


**FairchemTheory class:**

.. code-block:: python
    
    class FairchemTheory():
        def __init__(self, model_name=None, model_file=None, task_name=None, platform="cpu", seed=41, numcores=1,
                    periodic=False, periodic_cell_vectors=None, periodic_cell_dimensions=None):

.. list-table::
    :widths: 15 15 15 60
    :header-rows: 1

    * - Keyword
      - Type
      - Default value
      - Details
    * - ``model_name``
      - string
      - None
      - Name of model. Options: 'uma-s-1', 'uma-s-1p1', 'uma-m-1p1', 'esen-md-direct-all-omol', 'esen-sm-conserving-all-omol', 'esen-sm-direct-all-omol'
    * - ``model_file``
      - string
      - None
      - Name (or path) of file available used to load model.
    * - ``task_name``
      - string
      - None
      - Name task. Options: 'oc20', 'omol', 'omat', 'odac', 'omc'
    * - ``platform``
      - string
      - 'cuda'
      - Platform/device used for using Fairchem model within PyTorch. Options are: 'cuda', 'opencl', 'mps', 'cpu'.
    * - ``numcores``
      - integer
      - 1
      - Number of CPU cores used (if platform is CPU)
    * - ``periodic``
      - Boolean
      - False
      - Whether to use PBCs or not
    * - ``periodic_cell_vectors``
      - numpy array
      - None
      - Cell vectors as 3x3 numpy array in Angstrom.
    * - ``periodic_cell_vectors``
      - list
      - None
      - Cell dimensions as list of cell lengths and angles in Angstrom and degrees.

################################################################################
Fairchem installation
################################################################################

See Fairchem installation instructions: https://fair-chem.github.io/install/
Most likely the installation can be simply carried out using pip: 

.. code-block:: python

    pip install fairchem-core


################################################################################
Fairchem Examples with UMA models
################################################################################

To specify a model one should specify either *model_name* or *model_file* and also a *task_name* .

Currently available models include e.g. : 'uma-s-1p1' and 'uma-m-1p1'

*task_name* specifies a specific parameterization within each model. 

Options: 'oc20', 'omol', 'omat', 'odac', 'omc'
For molecular chemistry use: *task_name* ='omol'


See `Fairchem UMA models documentation <https://fair-chem.github.io/uma/>`_ for more information.

**Use a UMA model by model_name**

This option will automatically download the UMA model from Hugging-Face and use directly.
It requires a Hugging-Face account and one must have been granted access to the `UMA-model <https://huggingface.co/facebook/UMA>`_
One must also create a Hugging-Face token and set a shell environment variable. `Read here <https://fair-chem.github.io/install/#access-to-gated-models-on-huggingface>`_

Available models are documented at: https://fair-chem.github.io/uma

.. code-block:: python

    from ash import *

    frag = Fragment(databasefile="h2o.xyz")
    # Here selecting the recommended UMA-model uma-s-1p1
    theory = FairchemTheory(model_name="uma-s-1p1", task_name="omol")
    Singlepoint(theory=theory, fragment=frag)


**Use a UMA model by downloaded file**

If one has access to the `UMA-models <https://huggingface.co/facebook/UMA>`_ one can also download
the file, upload to the relevant computer/cluster and load them via the FairchemTheory *model_file* option.

.. code-block:: python

    from ash import *

    frag = Fragment(databasefile="h2o.xyz")
    # Here selecting the pre-downloaded uma-s-1p1.pt file 
    theory = FairchemTheory(model_file="/data/ML-models/uma-s-1p1.pt", task_name="omol")
    Singlepoint(theory=theory, fragment=frag)

**Use a UMA model with PBCs**

Here is an example of how to use periodic boundary conditions.

.. code-block:: python

  from ash import *

  numcores=1
  frag = Fragment(xyzfile="ammonia.xyz", charge=0, mult=1)

  #Cell
  cell_vectors = np.array([[5.01336,0.0,0.0],[0.0,5.01336,0.0],[0.0,0.0,5.01336]])

  #Periodic CP2KTheory definition with specified cell dimensions
  theory = FairchemTheory(model_file="/path/to/uma-s-1p1.pt",
                  periodic=True, periodic_cell_vectors=cell_vectors)

  # Here using the geomeTRIC Optimizer that recognizes that PBCs are active in FairChemTheory and optimizes atoms and cell positions
  Optimizer(theory=theory, fragment=frag, coordsystem="hdlc")