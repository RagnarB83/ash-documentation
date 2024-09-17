Machine learning in ASH
=========================================================

ASH is well-suited for utilizing machine-learning within computational chemistry,
being a scriptable Python library with interfaces to various quantum chemistry codes and the OpenMM molecular mechanics code.
This makes ASH convenient for generating training data for machine-learning potentials.
Additionally, as almost all ASH job-types are theory-agnostic, ML-potentials are just as valid as input to computational chemistry jobs within ASH (as long as an interface to that ML-potential is available).
Thanks to interfaces to PyTorch, MLatom and other machine learning libraries, 
the training of ML-potentials can also be performed within the ASH environment and pre-trained models can be utilized
directly in ASH just as easily as any other Theory object.

Machine learning capabilities of ASH are expected to grow in the future.

################################################################################
Creating training data for machine-learning potentials or Δ-learning
################################################################################

ASH features a function **create_ML_training_data** that can be used to generate energy or energy+gradient data,
suitable for training machine-learning potentials or Δ-learning corrections. 
Other properties are currently not supported, a manual workflow is needed for those.

.. code-block:: python
        
    # Function to create ML training data given XYZ-files and 2 ASH theories
    def create_ML_training_data(xyzdir=None, xyz_trajectory=None, theory_1=None, theory_2=None, 
                                num_snapshots=100, random_snapshots=False, charge=0, mult=1, Grad=True):

The user should either specify a directory with XYZ-files (*xyzdir* keyword) or a single XYZ trajectory file (*xyz_trajectory* keyword, file should contain multiple XYZ geometries in Xmol format).
The number of snapshots (geometries) should be given (*num_snapshots*, defaults to 100), and the snapshots can be randomized or not (*random_snapshots*, defaults to False).
Charge and multiplicity of the molecule should be provided (*charge* and *mult* keywords).
For training a potential using both energy and gradient data, set *Grad* to True (default), this is generally preferable (more accurate training). For energy-only training, set *Grad* to False. 

Theory levels should be provided via keywords *theory_1* and *theory_2*. Any ASHTheory object is in principle suitable.
If only *theory_1* is provided, the function will generate energy data for that level alone (suitable for training machine-learning potentials).
If both *theory_1* and *theory_2* are provided, Δ-learning mode is activated and the energy-difference (and gradient-difference if Grad=True) will be outputted (suitable for training Δ-learning corrections).

**Examples**

The example below assumes that you have already created either:

i) a directory containing XYZ-files of the molecules you want to train on, or 

ii) a single XYZ trajectory file containing multiple snapshots of the molecule.

The latter can e.g. come from a molecular dynamics simulation. 

.. code-block:: python

    from ash import *

    numcores=1

    #Variables
    method="HF-3c" #String that defines and ORCA-keyword
    num_snaps=100 #Number of snapshots to use

    #Training data directory
    #xyzdir="/Users/rb269145/ash-tests/ML-deltacorrection-3fgaba/individual-molecules"
    xyztraj = "/Users/rb269145/ash-tests/ML-deltacorrection-3fgaba/MD-data/walker0_trajectory.xyz"
    #Theory levels for delta_learning
    #Theory 1 (low-level)
    theory_gas=ORCATheory(orcasimpleinput=f"! {method} tightscf", numcores=numcores)
    #Theory 2 (high-level)
    theory_solv=ORCATheory(orcasimpleinput=f"! {method} CPCM(water) tightscf", numcores=numcores)

    #Call create_ML_training_data using 2 theory levels (delta-learning)
    #create_ML_training_data(xyz_dir=xyzdir, num_snapshots=num_snaps, random_snapshots=True,
    #    theory_1=theory_gas, theory_2=theory_solv, Grad=True)
    create_ML_training_data(xyz_trajectory=xyztraj, num_snapshots=num_snaps, random_snapshots=True,
        theory_1=theory_gas, theory_2=theory_solv, Grad=True)
    #produces files: train_data.xyz, train_data.energies, train_data.gradients

Now that the training data has been created it can be used as input to a machine-learning training library.
Below we show as an example how to use the ASH interface to the MLatom library to train a Δ-learning correction potential:

.. code-block:: python

    #Create MLatomTheory model and train
    ml_theory = MLatomTheory(ml_model="ANI", model_file=f"ANI-3fgaba_delta_snap{num_snaps}_{method}.pt")
    ml_theory.train(molDB_xyzfile="train_data.xyz", molDB_scalarproperty_file="train_data.energies",
                molDB_xyzvecproperty_file="train_data.gradients")



################################################################################
Interface to Torch and TorchANI
################################################################################

The interface to  `PyTorch <pytorch.org>`_ is documented at :doc:`torch_interface` that can be used for both ANI-style and AIMNet2 potentials.
This interface allows easy use of Torch-based machine-learning potentials in ASH.

.. code-block:: python

    from ash import *

    #H2O fragment
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
    # Create a TorchTheory object using the ANI1x neural network potential
    theory = TorchTheory(model_name="ANI1x", platform="cpu") #built-in
    #theory = TorchTheory(model_file="savedANI1x.pt") #from saved file
    
    #Run a geometry optimization
    Optimizer(theory=theory, fragment=frag)

################################################################################
Interface to MLatom
################################################################################

MLatom is a library for training and using ML potentials in computational chemistry.
The ASH interface can be used for both training and using ML-atom potentials.
See :doc:`MLatom-interface` for more.

################################################################################
Using machine-learning potentials in OpenMMTheory
################################################################################

A trained machine learning potential can be used directly by OpenMM thanks to 
the `OpenMM_Torch <https://github.com/openmm/openmm-torch>`_ and `OpenMM-ML <https://github.com/openmm/openmm-ml>`_ 
additions to OpenMM (need to be separately installed).
The advantage of using machine-learning potentials with OpenMM is that the simulation may run faster 
than other options requiring additional interfaces, as OpenMM is then responsible for propagating the system with 
optimized C++ or CUDA/OpenCL code. OpenMM can also be used for mixed systems where part is described by MM and part by ML.

If OpenMM-Torch is installed then a ML-force can be loaded and added to an OpenMMTheory object like this:

.. code-block:: python

    from ash import *
    from openmmtorch import TorchForce

    #Fragment
    frag = Fragment(pdbfile="file.pdb")

    #Load a Torch model from file using OpenMM-Torch to get an OpenMM-compatible force
    force = TorchForce('model.pt')

    #Create the ASH OpenMMTheory object without any force
    omm = OpenMMTheory(fragment=fragment, dummysystem=True)
    #Add ML force
    omm.add_force(mlforce)

    #Run a simulation e.g.
    MolecularDynamics(theory=omm, fragment=frag, simulation_steps=1000, timestep=0.001)


`OpenMM-ML <https://github.com/openmm/openmm-ml>`_ is a higher-level API that allows even easier use of 
pretrained built-in ML models together with OpenMM.
The most useful feature is to be able to easily create a mixed OpenMM system that uses both
MM forces and ML forces. 
The ASH interface allows easy creation of a mixed system like this:


.. code-block:: python

    from ash import *
    from openmmml import MLPotential

    pdbfile="relaxbox_NPT_lastframe.pdb"
    frag = Fragment(pdbfile=pdbfile)

    #Creating OpenMM object
    omm = OpenMMTheory(xmlfiles=["MOL_F57D69.xml"], pdbfile=pdbfile)

    mlpot = MLPotential('ani2x')  #Load the ANI2x ML potential
    mlatoms=[3069,3070,3071, 3072, 3073,3074] #Specify which atoms are ML
    omm.create_mixed_ML_system(mlpot,mlatoms) #Create the mixed ML/MM system

    # Run a simulation
    MolecularDynamics(theory=omm, fragment=frag, simulation_steps=1000, timestep=0.001)
