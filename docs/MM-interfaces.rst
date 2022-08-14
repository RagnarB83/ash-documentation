MM Interfaces
==========================

Molecular mechanics in ASH is possible via either the internal **NonBondedTheory** or via an interface to the external OpenMM package: :doc:`OpenMM-interface`

The internal **NonBondedTheory** is only capable of calculating nonbonded interactions: Coulomb electrostatics and short-range Lennard-Jones interactions.
It will use routines written in Julia and thus requires activating the Python-Julia interface.
**NonbondedTheory** can be used in geometry optimizations provided that the MM atoms are always frozen (should be specified by the **Optimizer**).
It can be used as MM code in QM/MM (:doc:`module_QM-MM`) theory objects.

The interface to the external OpenMM code (:doc:`OpenMM-interface`) allows for full-fledged molecular mechanics.
OpenMM is a fast C++ library with a Python API that ASH is interfaced to. 
OpenMM is capable of bonded and nonbonded interactions, can treat periodic and nonperiodic systems and can read in multiple types of forcefield files.
It can also run both on the CPU and the GPU.
It can be used as MM code in QM/MM (:doc:`module_QM-MM`) theory objects.

###########################
NonBondedTheory
###########################

The NonBondedTheory class:

.. code-block:: python


    class NonBondedTheory:
        def __init__(self, atomtypes=None, forcefield=None, charges = None, LJcombrule='geometric',
                     codeversion='julia', printlevel=2):


Defining a **NonBondedTheory** object is easy and can be accomplished in a few different ways.
A simple way is to define all information in the script itself. This requires defining the MM forcefield as a dictionary
and then provides a list of atomtypes of the system as a minimum. In the forcefield we then define the Coulomb and Lennard-Jones parameters
associated with the atomtypes.

*Simple way (forcefield_dict and atomtypes):*

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


*An alternative is to define the forcefield in a forcefieldfile that is read-in.*

.. code-block:: python

    from ash import *

    HF_frag=Fragment(xyzfile="hf.xyz")
    #Defining atomtypes 'HT' for hydrogen and atomtype 'FX' for fluorine. These atomtypes can be named anything.
    atomtypes =['HT', 'FX']
    MM_forcefield=MMforcefield_read('forcefield.ff')

    MMobject = NonBondedTheory(forcefield=MM_forcefield, atomtypes=atomtypes, LJcombrule='geometric')

    Singlepoint(fragment=HF_frag,theory=MMobject)

where forcefield.ff contains e.g.:

.. code-block:: text

    LennardJones_i_sigma FX 3.150574227 -0.1521
    LennardJones_i_sigma HT 3.550053212    -0.07
    charge FX -0.9
    charge HT 0.9

The forcefield file will read and parse lines like:


.. code-block:: text

    LennardJones_i_sigma <atomtype> <sigma> <epsilon> # Specify atomtype and sigma-value (Å) and epsilon value (kcal/mol)
    LennardJones_i_R0 <atomtype> <R0> <epsilon>  # Specify atomtype and R0-value (Å) and epsilon value (kcal/mol)
    LennardJones_ij <atomtype1> <atomtype2>  <R0_ij> <epsilon pair parameter> #Specify pair-potential. Currently inactive option
    charge <atomtype> <chargevalue> # Specify atomtype and charge-value


or other options:

.. code-block:: text

    combination_rule <combrule option>  #For LJ potential. Can be geometric, arithmetic, mixed_geoepsilon, mixed_geosigma
    XX_atomtypes <list of atomtypes> #Atomtypes for residue XX. List is space-separated.
    XX_charges <list of charges> #Charges for residue XX. List is space-separated.
    XX_elements <list of elements> #Elements for residue XX. List is space-separated.


###########################
OpenMM interface
###########################


See :doc:`OpenMM-interface`

