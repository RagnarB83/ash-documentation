Torch interface
======================================

PyTorch is a very general deep-learning framework within the Python ecosystem and one of the most popular approches to training neural networks.
TorchANI is a PyTorch-based implementation of the ANI neural network potential for describing potential energy surfaces of molecular systems (and other properties).
AIMNet2 is another even newer neural network potential, also implemented using PyTorch, capable of unprecedented accuracy for up to 14 chemical elements, 
support for charged systems and is capable of describing long-range interactions and dispersion.
AIMNet2 paper: https://doi.org/10.26434/chemrxiv-2023-296ch-v2


ASH features a basic interface to PyTorch that furthermore supports both TorchANI and AIMNet2 neural network models,
that allows easy use of pre-trained neural network potentials and even training of new models. 
This allows the use of such neural network potentials as any other theory-level within ASH.

Energies and gradients can be requested (just like a regular QM or MM theory) and so a valid TorchTheory object can be used
for single-point energies, geometry optimizations, numerical frequencies, surface scans, NEB, molecular dynamics etc. within ASH.
Even hybrid ONIOM and QM/MM calculations are possible (with some limitations).


**TorchTheory class:**

.. code-block:: python
    
    class TorchTheory():
        def __init__(self, filename="torch.pt", model_name=None, model_object=None,
                    model_file=None, printlevel=2, label="TorchTheory", numcores=1, 
                    platform=None, train=False):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``filename``
     - string
     - 'torch.pt'
     - Filename used for PyTorch models when saving.
   * - ``model_name``
     - string
     - None
     - Name of pretrained model to use. Options: 'AIMNet2','ANI1x', 'ANI2x', 'ANI1ccx' (requires TorchANI or AIMNet2)
   * - ``model_object``
     - PyTorch NN model object
     - None
     - Read in a PyTorch model object directly.
   * - ``model_file``
     - Name of Pytorch model-file to read in
     - None
     - Read in a PyTorch model from a file.
   * - ``platform``
     - string
     - None
     - Type of platform to use when training/running. Options: 'cpu', 'cuda', 'mps' (Apple Silicon)
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``label``
     - string
     - 'TorchTheory'
     - Label of TorchTheory object
   * - ``numcores``
     - integer
     - 1
     - Number of cores.
   * - ``train``
     - Boolean
     - False
     - Whether to create a TorchTheory object that will be used for training.

################################################################################
Torch/TorchANI/AIMNet2 installation
################################################################################

*PyTorch*

 `PyTorch <pytorch.org>`_  needs to be installed on your system to use TorchTheory. The easiest way to install PyTorch is via pip:

.. code-block:: python

    pip install torch

*TorchANI*

To use pre-trained ANI neural network potentials you need to install TorchANI.
See `TorchANI documentation <https://aiqm.github.io/torchani/>`_ and `TorchANI repository <https://github.com/aiqm/torchani>`_
Torchani can be installed via pip:

.. code-block:: python

    pip install torchani

*AIMNet2*

To use AIMNet2 follow the installation instructions at https://github.com/isayevlab/AIMNet2
We installed it like this in a conda environment where Pytorch and ASH had already been installed:

.. code-block:: shell

  git clone https://github.com/isayevlab/AIMNet2.git
  cd AIMNet2
  python setup.py install

################################################################################
AIMNet2 Examples
################################################################################

*Basic AIMNet2 example*

It is easy to use the AIMNet2 neural network potential with TorchTheory.
The available models are: 'AIMNet2', and it is available for elements: 'H', 'C', 'N', 'O', 'F', 'Cl', 'S', 'Si', 'B', 'P', 'Br', 'As', 'I', 'Se'

.. code-block:: python

    from ash import *

    #H2O fragment
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
  
    # Create a TorchTheory object using the AIMNet2 neural network potential
    theory = TorchTheory(model_name="AIMNet2", platform="cpu")
    
    #Run a single-point energy+gradient calculation
    #Optimizer,NumFreq, MolecularDynamics etc. should also work
    result = Singlepoint(theory=theory, fragment=frag, Grad=True)

    print(result.energy)
    print(result.gradient)

################################################################################
TorchANI Examples
################################################################################

*Basic TorchANI example*

A pretrained ANI-based model using the TorchANI library can easily be used as well.
The available models are: 'ANI1x', 'ANI1ccx', 'ANI2x' and they are available for elements: 'H', 'C', 'N', 'O'

.. code-block:: python

    from ash import *

    #H2O fragment
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
    # Create a TorchTheory object using the ANI2x neural network potential
    theory = TorchTheory(model_name="ANI2x", platform="cpu")
    
    #Run a single-point energy+gradient calculation
    #Optimizer,NumFreq, MolecularDynamics etc. should also work
    result = Singlepoint(theory=theory, fragment=frag, Grad=True)

    print(result.energy)
    print(result.gradient)

*Loading a pretrained model from file*

It's also possible to load a neural-network potential from file.
Here we show an example of this by loading the ANI1x model from a file (Pytorch .pt format).
This file was generated by first creating a TorchTheory object like above and then calling the *save_model* method.

.. code-block:: python

    from ash import *

    #H2O fragment
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
    # Create a TorchTheory object using the ANI1x neural network potential from a saved-file
    theory = TorchTheory(model_file="savedANI1x.pt")
    #Run a single-point energy+gradient calculation
    result = Singlepoint(theory=theory, fragment=frag, Grad=True)

*Defining a new PyTorch model from scratch*

Not yet ready