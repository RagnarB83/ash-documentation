MACE interface
======================================

`MACE <https://mace-docs.readthedocs.io/>`_ is an equivariant message passing neural network potential and software that looks very promising.

ASH features a simple interface to MACE via the PyTorch implementation that allows both the direct use of pre-trained MACE models
to be used in ASH (just like any other theory) as well as direct training of a MACE-potential.

Energies and gradients can be requested (just like a regular QM or MM theory) and so a valid 
MACETheory object can be used for single-point energies, geometry optimizations, numerical frequencies, surface scans, NEB, molecular dynamics etc. within ASH.
Even hybrid ONIOM and QM/MM calculations are possible (with some limitations).


**MACETheory class:**

.. code-block:: python
    
    class MACETheory():
        def __init__(self, config_filename="config.yml", 
                    model_file=None, printlevel=2, 
                    label="MACETheory", numcores=1, device="cpu"):

.. list-table::
    :widths: 15 15 15 60
    :header-rows: 1

    * - Keyword
      - Type
      - Default value
      - Details
    * - ``config_filename``
      - string
      - 'config.yml'
      - Filename used for for MACE training configuration file.
    * - ``model_file``
      - string
      - 'mace.model'
      - Name of file used to load model (when using) or file created when training new model.
    * - ``printlevel``
      - integer
      - 2
      - Printlevel
    * - ``label``
      - string
      - 'MACETheory'
      - Label used for object
    * - ``device``
      - string
      - 'cpu'
      - Device used for training or using of MACE model within PyTorch. Options are: 'cuda', 'opencl', 'mps', 'cpu'.
    * - ``numcores``
      - integer
      - 1
      - Number of CPU cores used (if device is CPU)



**MACETheory train method:**

.. code-block:: python

    def train(self, config_file="config.yml", name="model",model="MACE", device='cpu',
                      valid_fraction=0.1, train_file="train_data_mace.xyz",E0s=None,
                      energy_key='energy_REF', forces_key='forces_REF',        
                      energy_weight=1, forces_weight=100,
                      max_num_epochs=500, swa=True, batch_size=10,
                      max_L = 0, r_max = 5.0, num_channels=128,  
                      results_dir= "MACE_models", checkpoints_dir = "MACE_models", 
                      log_dir ="MACE_models", model_dir="MACE_models"):

################################################################################
MACE installation
################################################################################

See MACE installation instructions: https://mace-docs.readthedocs.io/en/latest/guide/installation.html#installations
Most likely the installation can be simply carried out using pip: 

.. code-block:: python

    pip install mace-torch


################################################################################
MACE training
################################################################################

Training a MACE model is performed by first creating a MACETheory object and then calling the train method of that object.

The train method features a lot of options with the following defaults.

.. code-block:: python

    def train(self, config_file="config.yml", name="model",model="MACE", device='cpu',
                      valid_fraction=0.1, train_file="train_data_mace.xyz",E0s=None,
                      energy_key='energy_REF', forces_key='forces_REF',        
                      energy_weight=1, forces_weight=100,
                      max_num_epochs=500, swa=True, batch_size=10,
                      max_L = 0, r_max = 5.0, num_channels=128,  
                      results_dir= "MACE_models", checkpoints_dir = "MACE_models", 
                      log_dir ="MACE_models", model_dir="MACE_models"):

To train a MACE model one needs a file containing the training data (specify by *train_file* keyword).
The file is in a modified XYZ format that contains coordinates, energies and forces.
It is best to use the *create_ML_training_data* function (see :doc:`Machine_learning_in_ASH`)
to generate this file (named train_data_mace.xyz).
One can then train a model in the simple manner below, which uses default training parameters (see above).

.. code-block:: python

    theory = MACETheory()
    theory.train(train_file="train_data_mace.xyz")

Calling the train method calls MACE and begins training which can take a while depending on the size of the training set, size of molecules, 
number of epochs requested (500 by default) etc. The training uses the CPU by default (device='cpu') but can run much faster on the GPU using e.g. device='cuda'.


################################################################################
MACE Examples
################################################################################

**Train a MACE model**

Here we train a MACE potential from scratch.
We first need to create some basic training data for a small molecule in the form of geometries and energies and forces.
Here we generate a high-temperature Wigner ensemble from a Hessian calculated at our target level of theory. 
This is an easy way of generating some relevant geometries (though unknown how good this choice is for any given molecule).

We then call the **create_ML_training_data** function to generate energies and forces for these geometries or possibly a subset (this is our training data).

After that we can create a MACETheory object and then call the *train* method of the object using the train_data_mace.xyz file, 
a multigeometry XYZ-file that contains coordinates, energies and forces (generated by **create_ML_training_data**).
We can choose to either train on the CPU (device="cpu") or alternatively using a GPU (use device="cuda" or device="opencl").
Here we modify the weight on forces vs. energy to be 100:1 (forces_weight=100).
The code below will create multiple model files but the file named model_stagetwo_compiled.model would typically be the one to use.

.. code-block:: python

    from ash import *

    # Define initial fragment
    frag = Fragment(xyzfile="3fgaba.xyz", charge=0, mult=1)

    # Define a QM theory
    xtb = xTBTheory(xtbmethod="GFN2-xTB")

    ###########################################
    # Generate training data (Wigner ensemble)
    ###########################################
    #Run Hessian
    result_freq = NumFreq(theory=xtb, fragment=frag)
    # Generate Wigner ensemble of 100 samples at elevated temperature
    wigner_frags = wigner_distribution(fragment=frag, hessian=result_freq.hessian,
        temperature=500, num_samples=100)
    # This creates file Wigner_traj.xyz

    # Read the Wigner ensemle snapshots and generate energies and forces training data 
    create_ML_training_data(xyz_trajectory="Wigner_traj.xyz", theory_1=xtb, Grad=True)
    # This creates file train_data_mace.xyz
    
    ###################################
    # TRAIN MACE model
    ###################################
    #Create MACE model object
    macetheory = MACETheory()
    #Train model
    macetheory.train(train_file="train_data_mace.xyz", device="cpu", max_num_epochs="3", forces_weight=100)



**Use a previously trained model**

Here it is assumed that you have already either trained a model yourself previously (see above) or that you want to use a pre-trained MACE model
from the literature. 

.. code-block:: python

    from ash import *

    frag = Fragment(xyzfile="3fgaba.xyz", charge=0, mult=1)

    #Create MACE model object and point to modelfile
    macetheory = MACETheory(model_file="model_stagetwo_compiled.model")
    
    #Run a job with Theory
    res = Optimize(theory=macetheory, fragment=frag)