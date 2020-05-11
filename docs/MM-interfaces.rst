==========================
MM Interfaces
==========================

Molecular mechanics in Ash is currently quite rudimentary.
There is no full-fledged internal forcefield code available (i.e. that handles both bonded and nonbonded terms).
There is, however, a flexible nonbonded forcefield code (see NonBondedTheory below) that allows for Coulomb and Lennard-Jones
energy+gradient evaluations. This allows for rigid (MM atoms frozen) MM and QM/MM (:doc:`QM-MM`) theory objects that can be used in geometry optimizations
and molecular dynamics simulations (:doc:`job-types`).

In addition, an interface to OpenMM is available.


###########################
NonBondedTheory
###########################

Defining a NonBondedTheory object is easy and can be accomplished in a few different ways.


A simple way is to define all information in the script itself. This requires defining the MM forcefield as a dictionary
and then provides a list of atomtypes of the system as a minimum. In the forcefield we then define the Coulomb and Lennard-Jones parameters
associated with the atomtypes.

Simple way (forcefield_dict and atomtypes):

.. code-block:: python

    from ash import *
    import sys
    settings_ash.init() #initialize

    HF_frag=Fragment(xyzfiles="hf.xyz")
    #Defining atomtypes 'HT' for hydrogen and atomtype 'FX' for fluorine. These atomtypes can be named anything.
    atomtypes =['HT', 'FX']

    #MM_forcefield = {}
    #MM_forcefield[atomtype]=AtomMMobject()
    #MM_forcefield[atomtype].add_charge(atomcharge=charge)
    #MM_forcefield[atomtype].add_LJparameters(LJparameters=[sigma_i,eps_i])
    #HF_MM_forcefield= {'HT' : LJparameters}


    MMobject = NonBondedTheory(forcefield=MM_forcefield, atomtypes=atomtypes, LJcombrule='geometric')
    MMobject.run()


Alternative is to define the forcefield in a forcefieldfile that is read-in.

###########################
OpenMM interface
###########################
TODO...