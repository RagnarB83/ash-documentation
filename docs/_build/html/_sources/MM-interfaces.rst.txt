==========================
MM Interfaces
==========================

Molecular mechanics in ASH is possible via either the internal NonBondedTheory or via an interface to the external OpenMM package: :doc:`OpenMM-interface`

The internal NonBondedTheory can only perform nonbonded interactions: electrostatics and short-range Lennard-Jones interactions.
The Coulomb+LJ interaction is quite fast as it is written in Julia.
Nonbonded Theory can be used in geometry optimizations provided that the MM atoms are always frozen.
It can be used as MM code in QM/MM (:doc:`module_QM-MM`) theory objects.

The interface to the external OpenMM code (:doc:`OpenMM-interface`) has recently become available. OpenMM is a fast C++ code with a Python API that ASH is interfaced to. OpenMM can calculate both bonded and nonbonded interactions, both periodic and nonperiodic systems and can read in multiple types of forcefield files.


###########################
NonBondedTheory
###########################

The NonBondedTheory class:

.. code-block:: python


    class NonBondedTheory:
        def __init__(self, atomtypes=None, forcefield=None, charges = None, LJcombrule='geometric',
                     codeversion='julia', printlevel=2):


Defining a NonBondedTheory object is easy and can be accomplished in a few different ways.
A simple way is to define all information in the script itself. This requires defining the MM forcefield as a dictionary
and then provides a list of atomtypes of the system as a minimum. In the forcefield we then define the Coulomb and Lennard-Jones parameters
associated with the atomtypes.

Simple way (forcefield_dict and atomtypes):

.. code-block:: python

    from ash import *

    HF_frag=Fragment(xyzfile="hf.xyz")
    #Defining atomtypes 'HT' for hydrogen and atomtype 'FX' for fluorine. These atomtypes can be named anything.
    atomtypes =['HT', 'FX']

    MM_forcefield = {}
    MM_forcefield[atomtype]=AtomMMobject()
    MM_forcefield[atomtype].add_charge(atomcharge=charge)
    MM_forcefield[atomtype].add_LJparameters(LJparameters=[sigma_i,eps_i])
    HF_MM_forcefield= {'HT' : LJparameters}

    MMobject = NonBondedTheory(forcefield=MM_forcefield, atomtypes=atomtypes, LJcombrule='geometric')

    Singlepoint(fragment=HF_frag,theory=MMobject)


Alternative is to define the forcefield in a forcefieldfile that is read-in.

.. code-block:: python

    from ash import *

    HF_frag=Fragment(xyzfile="hf.xyz")
    #Defining atomtypes 'HT' for hydrogen and atomtype 'FX' for fluorine. These atomtypes can be named anything.
    atomtypes =['HT', 'FX']
    MM_forcefield=MMforcefield_read('forcefield.ff')

    MMobject = NonBondedTheory(forcefield=MM_forcefield, atomtypes=atomtypes, LJcombrule='geometric')

    Singlepoint(fragment=HF_frag,theory=MMobject)

###########################
OpenMM interface
###########################


See :doc:`OpenMM-interface`

