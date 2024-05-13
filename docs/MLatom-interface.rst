MLatom interface
======================================

`The MLatom program <http://mlatom.com>`_ (see also `XACS_MLatom documentation <https://xacs.xmu.edu.cn/docs/mlatom/>`_) is a program for AI/ML enchanced quantum chemistry, developed by Pavlo Dral's research group at Xiamen university.

Mlatom features various interfaces to machine-learning potentials and QM software and can be used for both training and prediction of molecular energies and properties.
ASH features a basic interface to the MLatom Python API that allows direct use of various pre-trained models using the ML potentials supported by MLatom. Some basic training is also supported.

If a valid *MLatomTheory* object is created using a pretrained model or a model is correctly trained, an *MLatomTheory* object will behave like any other Theory level within ASH.
That is: energies and gradients can be requested (just like a regular QM or MM theory) and so *MLatomTheory* can be used for single-point energies, geometry optimizations, 
numerical frequencies, surface scans, NEB , molecular dynamics etc. within ASH. 

It is also possible to use ASH to provide training data (usually energies and gradients) to MLatom and to train a new ML model (even directly within ASH).

WARNING: As the interface to MLatom is new and MLatom is under rapid development, the interface may regularly change.


**MLatomTheory class:**

.. code-block:: python
    
    class MLatomTheory:
        def __init__(self, method=None, ml_model=None, model_file=None, qm_program=None, ml_program=None,
                    printlevel=2, numcores=1, label="mlatom"):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``method``
     - string
     - None
     - Name of pretrained method to. Options: 'AIQM1', 'AIQM1\@DFT', 'AIQM1\@DFT*', 'ANI-1ccx', 'ANI-1x', 'ANI-1x-D4', 'ANI-2x', 'ANI-2x-D4'
   * - ``ml_model``
     - string
     - None
     - Name of ML model to use. Options: 'ani', 'dpmd', 'gap', 'physnet', 'sgdml', 'mace'.
   * - ``ml_program``
     - string
     - None
     - Name of helper ML-program. Used by ml_model='kreg' (Options: 'KREG_API' and 'MLatomF')
   * - ``qm_program``
     - string
     - None
     - Name of QM-program that MLatom may use as part of the method (e.g. 'mndo' or 'sparrow').
   * - ``model_file``
     - string
     - None
     - Read in a model from a file.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``label``
     - string
     - 'mlatom'
     - Label of MLatomTheory object
   * - ``numcores``
     - integer
     - 1
     - Number of cores. 

################################################################################
MLatom installation
################################################################################


 `MLatom <mlatom.com>`_  needs to be installed on your system. 
 Installation instructions can be found at `XACS Docs <https://xacs.xmu.edu.cn/docs/mlatom/installation.html>`_ and `MLatom docs <http://mlatom.com/docs/installation.html>`_.
 The easiest way to install mlatom is via pip:

.. code-block:: python

    pip install mlatom

Additionally, mlatom has a few optional (but mosty needed) dependencies:

.. code-block:: python

    pip install sgdml rmsd openbabel xgboost scikit-learn pyscf rmsd rdkit pandas ase fortranformat tensorflow geometric

Additional dependencies may be needed depending on the specific ML-model form to be used.
See `MLatom documentation <http://mlatom.com/docs/installation.html>`_ and `XACS-MLatom documentation <https://xacs.xmu.edu.cn/docs/mlatom/installation.html>`_ for more information.


################################################################################
Examples
################################################################################

To use MLatom we need to choose to use either a method or a model (MLatom syntax).

- A MLatom-method is a general pretrained ML model and is designed to be general and work outside the box (just like a DFT method). It is specified by the *method* keyword in the ASH interface (a string). Examples of methods are: 'AIQM1' and 'ANI-1x'.
- A MLatom-model is a ML-model that needs to have a specified form can be kernel-based, neural-network based etc.) and needs to be trained or parameters loaded. It is specified by the *ml_model keyword* in the ASH interface.

MLatomTheory requires you to specify either a *method* or a *ml_model* when defining the object.

*Pretrained AIQM1 method example*

Since the `AIQM1 model <https://www.nature.com/articles/s41467-021-27340-2>`_ is built on top of a semiempirical QM method (ODM2), 
we also need to specify the semiempirical QM program that MLatom will use to define the ODM2 Hamiltonian. The options are: 'mndo' and 'sparrow' and these programs need to be separately installed on your system (and available in PATH).

.. code-block:: python

    from ash import *

    frag = Fragment(databasefile="glycine.xyz")
    theory = MLatomTheory(method="AIQM1", qm_program="mndo")
    Singlepoint(theory=theory, fragment=frag, Grad=True)

*Pretrained ANI-1x method example*

The ANI models (ANI-1ccx, ANI-1x, ANI-1x-D4, ANI-2x, ANI-2x-D4), based on the ANI neural network potentials are available in MLatom.
They require pytorch and torchani to be installed.
See also :doc:`torch_interface` for direct use of TorchANI/PyTorch (without MLatom).

.. code-block:: python

    from ash import *

    frag = Fragment(databasefile="glycine.xyz")
    theory = MLatomTheory(method="ANI-1x")
    Singlepoint(theory=theory, fragment=frag)

*Loading and running pretrained model from file*

We next show how to use a ML-model (*ml_model* keyword). If the training has already been performed and available as a file, can we load it.
First we have to choose what type of ML-model potential we want to use. The options are: 'ani', 'dpmd', 'gap', 'physnet', 'sgdml', 'mace'.
Next we must choose the file containing the model. This file often has a .pt suffix (for pytorch models) or a .pkl suffix (for scikit-learn models) or various other extensions.

.. code-block:: python

    from ash import *

    #Here defining a MACE ML-model (requires installing MACE separately) 
    #And downloading init.xyz and mace.pt from here: https://xacs.xmu.edu.cn/docs/mlatom/tutorial_geomopt.html
    theory = MLatomTheory(ml_model="mace", model_file="mace.pt")
    #theory = MLatomTheory(ml_model="ani", model_file="ani_model.pt")
    #theory = MLatomTheory(ml_model="kreg", model_file="kreg_model.unf")


    #Defining a molecule Fragment. NOTE: This must match the training data used to train the model (same molecule, same atom-order etc.)
    #See https://xacs.xmu.edu.cn/docs/mlatom/tutorial_geomopt.html for the init.xyz file
    frag = Fragment(xyzfile="init.xyz")

    Singlepoint(theory=theory, fragment=frag)

*Training a new model using MLatomTheory*

ASH features a very basic way to train a new ML model using the MLatom API.
It should be noted that training a new ML model can be a labororious, complicated process and it may be better to use MLatom directly (either the PythonAPI or the command-line interface) to have more control over the training process.
ASH and it's interfaces to various QM programs can still be used to generate the training data.
See `MLatom training documentation <https://xacs.xmu.edu.cn/docs/mlatom/tutorial_mlp.html#training>`_

Currently ASH can be used to train very basic ML-model potentials based on energies and gradients like the following examples.

See also `MLatom Machine learning potentials tutorial <https://xacs.xmu.edu.cn/docs/mlatom/tutorial_mlp.html>`_ for a tutorial on training machine learning potentials in general,
as well as links to download training data used below (H2.xyz, H2_HF.en, H2_HF.grad).

**ANI-example**

.. code-block:: python

    from ash import *

    #Create MLatomTheory model
    theory = MLatomTheory(ml_model="ANI")
    #Train model using 3 databasefiles containing XYZ-coords, energies and gradients
    #Download from; https://xacs.xmu.edu.cn/docs/mlatom/tutorial_mlp.html
    theory.train(molDB_xyzfile="H2.xyz", molDB_scalarproperty_file="H2_HF.en",
                molDB_xyzvecproperty_file="H2_HF.grad")
    #Model is now trained and can be used directly,


    #Molecule Fragment to use for simulation (needs to be compatible with training data)
    frag = Fragment(diatomic="H2", bondlength=1.0, charge=0, mult=1)

    result = Singlepoint(theory=theory, fragment=frag, Grad=True)

    print("Energy:", result.energy)
    print("Gradient:", result.gradient)

    result = Optimizer(theory=theory, fragment=frag, Grad=True)



**KREG-example**

.. code-block:: python

    from ash import *

    #Create MLatomTheory model
    theory = MLatomTheory(ml_model="kreg", ml_program='MLatomF')
    #Train model using 3 databasefiles containing XYZ-coords, energies and gradients
    #Download from; https://xacs.xmu.edu.cn/docs/mlatom/tutorial_mlp.html
    theory.train(molDB_xyzfile="H2.xyz", molDB_scalarproperty_file="H2_HF.en",
                molDB_xyzvecproperty_file="H2_HF.grad")
    #Model is now trained and can be used directly,


    #Molecule Fragment to use for simulation (needs to be compatible with training data)
    frag = Fragment(diatomic="H2", bondlength=1.0, charge=0, mult=1)

    result = Singlepoint(theory=theory, fragment=frag, Grad=True)

    print("Energy:", result.energy)
    print("Gradient:", result.gradient)

    result = Optimizer(theory=theory, fragment=frag)