MACE interface
======================================

`MACE <https://mace-docs.readthedocs.io/>`_ is an equivariant message passing neural network potential and software that looks very promising.

ASH features an interface to MACE via the PyTorch implementation that allows both the direct use of pre-trained MACE models
to be used in ASH (just like any other theory) as well as direct training of a MACE-potential.

Energies and gradients can be requested (just like a regular QM or MM theory) and so a valid 
MACETheory object can be used for single-point energies, geometry optimizations, numerical frequencies, surface scans, NEB, molecular dynamics etc. within ASH.
Even hybrid ONIOM and QM/MM calculations are possible (with some limitations).

Periodic boundary conditions are now also supported for running pretrained models.

**MACETheory class:**

.. code-block:: python
    
    def __init__(self, config_filename="config.yml",
                 model_name=None, model_name_subtype=None, model_name_head=None,
                 model_file=None, printlevel=2, mace_load_dispersion=False,
                 label="MACETheory", numcores=1, platform="cpu", device=None, return_zero_gradient=False, 
                 polarmace=False, default_dtype="float64",
                 energy_weight=None, forces_weight=None, max_num_epochs=None, valid_fraction=None,
                 periodic=False, periodic_cell_vectors=None, periodic_cell_dimensions=None):

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
      - None
      - Name of file used to load model (when using) or file created when training new model.
    * - ``model_name``
      - string
      - None
      - Name of foundational model used to load.
    * - ``model_name_subtype``
      - string
      - None
      - Name of foundational model subtype to choose (when appropriate).
    * - ``model_name_head``
      - string
      - None
      - For multihead models, which head to choose.
    * - ``printlevel``
      - integer
      - 2
      - Printlevel
    * - ``label``
      - string
      - 'MACETheory'
      - Label used for object
    * - ``platform``
      - string
      - 'cpu'
      - Platform (device) used for training or running of MACE model within PyTorch. Options are: 'cuda', 'opencl', 'mps', 'cpu'.
    * - ``default_dtype``
      - string
      - 'float64'
      - Whether to use float64 or float32 in Pytorch.
    * - ``numcores``
      - integer
      - 1
      - Number of CPU cores used (if device is CPU)
    * - ``energy_weight``
      - float
      - None
      - Weight on energy for training
    * - ``forces_weight``
      - float
      - None
      - Weight on forces for training
    * - ``max_num_epochs``
      - integer
      - None
      - Max number of epochs to use for training
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
Loading pre-trained or foundational MACE models
################################################################################
Pre-trained MACE models can either be previously trained models by yourself or a general foundational model.
If you simply want to load a MACE model in MACETheory, for the purpose of using it as a Theory object to run calculations,
you can either load the model via file (*model_file* keyword) or by name in case of a foundational model (*model_name* keyword).
For option *model_name* the model is looked up inside the mace-torch library and automatically downloaded.
Sometimes it may be necessary also to specify a subtype (*model_name_head*) or a specific head in the case of multi-head models (*model_name_head*).

For information on available foundational models, see:

- `MACE-docs Foundational models <https://mace-docs.readthedocs.io/en/latest/guide/foundation_models.html>`_
- `MACE-foundations repository <https://github.com/ACEsuit/mace-foundations>`_
- `MACE modeles on Hugging Face <https://huggingface.co/mace-foundations>`_

Examples on different ways to load models into MACETheory:

.. code-block:: python

  frag = Fragment(xyzfile="acetone.xyz", charge=0, mult=1)

  # Load a MACE foundational model via previously downloaded file (here MACE-OFF23)
  theory = MACETheory(model_file="/path/to/MACE-foundational-models/MACE-OFF23_medium.model")
  # Load a MACE foundational model via name (here MACE-OFF23)
  theory = MACETheory(model_name="mace-off23")
  # Load a MACE foundational model via name (here MACE-OFF23) and specify subtype (size)
  theory = MACETheory(model_name="mace-off23", model_name_subtype="large")
  # Load a MACE MH-1 foundational model via name and specify a specific head 
  theory = MACETheory(model_file="/path/to/MACE-foundational-models/mace-mh-1.model", model_name_head="omat_pbe")


Recognized *model_name* options by MACETheory in ASH:

- mace-off23 (*model_name_subtype* options: 'small', 'medium', 'large')
- mace-mp or mace-mh (*model_name_subtype* options: 'small', 'medium', 'large' or 'medium-mpa-0', 'mh-1', 'mh-0')
- mace-polar  (*model_name_subtype* options: 'polar-1-s', 'polar-1-m', 'polar-1-l')
- mace-ani-cc
- mace_omol

Note that for mace-mp / mace-mh models, D3 dispersion can also be activated via the *mace_load_dispersion* =True keyword and the *mace_dispersion_xc* keyword can be used to specify
the name of the DFT functional to use for evaluating the dispersion formula (e.g. *mace_dispersion_xc*="PBE" for PBE-D3 correction).
Dispersion is often recommended for use with those models (see literature). Note that the torch-dftd library is needed for this (pip install torch-dftd)
Code block below shows how to perform a mace-mh-1-omat-D3 model calculation (see preprint: https://arxiv.org/pdf/2510.25380):

.. code-block:: python

  # Choosing MACE-MH-1 model by name, subtype and head. And activating D3 dispersion
  qm = MACETheory(model_name="mace-mh-1", model_name_subtype="mh-1", model_name_head="omat_pbe",
                  mace_load_dispersion=True, mace_dispersion_xc="pbe")
  # Choosing MACE-MH-1 model by file and head.  And activating D3 dispersion
  qm = MACETheory(model_file="/path/to/MACE-foundational-models/mace-mh-1.model", model_name_head="omat_pbe",
                  mace_load_dispersion=True, mace_dispersion_xc="pbe")


Note that for other MACE foundational models it is easiest to download them separately and select via *model_file* keyword.

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

Here it is assumed that you have already either trained a model yourself previously (see above) or that you want to use a foundational MACE model
from the literature. 

.. code-block:: python

    from ash import *

    frag = Fragment(xyzfile="3fgaba.xyz", charge=0, mult=1)

    #Create MACE model object and point to modelfile
    macetheory = MACETheory(model_file="model_stagetwo_compiled.model")
    
    #Run a job with Theory
    res = Optimize(theory=macetheory, fragment=frag)

**Run PBC calculations with MACE foundational models**

MACE foundational models can be used for periodic calculations.
Example below activates PBCs in MACETheory object and then uses the geomeTRICOptimizer to 
optimize atom and cell positions.

.. code-block:: python

  from ash import *

  numcores=1
  frag = Fragment(xyzfile="ammonia.xyz", charge=0, mult=1)

  #Cell
  cell_vectors = np.array([[5.01336,0.0,0.0],[0.0,5.01336,0.0],[0.0,0.0,5.01336]])

  #Periodic CP2KTheory definition with specified cell dimensions
  theory = MACETheory(model_file="/path/to/MACE-foundational-models/MACE-POLAR-1-M.model",
                  periodic=True, periodic_cell_vectors=cell_vectors)

  # Here using the geomeTRIC Optimizer that recognizes that PBCs are active in MACETheory and optimizes atoms and cell positions
  Optimizer(theory=theory, fragment=frag, coordsystem="hdlc")