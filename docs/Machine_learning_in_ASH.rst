Machine learning in ASH
=========================================================

ASH is a well-suited program for utilizing machine-learning within computational chemistry,
being a scriptable Python library with interfaces to various quantum chemistry codes and the OpenMM molecular mechanics code.
This makes ASH convenient for generating training data for machine-learning potentials.
Additionally, as the ASH job-types are theory-agnostic, ML-potentials are just as valid as input to computational chemistry jobs within ASH.
Thanks to interfaces to PyTorch, MLatom and other machine learning libraries, 
the training can also be performed within the ASH environment and pre-trained models can be utilized
directly in ASH just as easily as any other Theory object.

Machine learning capabilities of ASH are expected to grow in the future.

################################################################################
Interface to Torch and TorchANI
################################################################################

The interface to  `PyTorch <pytorch.org>`_ and `TorchANI <https://aiqm.github.io/torchani/>`_ 
is documented at :doc:`torch_interface`.

It allows easy use of Torch-based machine-learning potentials in ASH.

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

    #Load a Torch model using OpenMM-Torch to get an OpenMM-compatible force
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

    from ash import *#
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


################################################################################
Interface to ML-atom
################################################################################

ML-atom is a library for training and using ML potentials in computational chemistry.
The ASH interface can be used for both training and using ML-atom potentials.
See :doc:`MLatom-interface` for more.