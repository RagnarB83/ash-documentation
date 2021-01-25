==========================
MM Interfaces
==========================

Molecular mechanics in ASH is possible via either the internal NonBondedTheory or via an interface to the external OpenMM package.

The internal NonBondedTheory can only perform nonbonded interactions: electrostatics and short-range Lennard-Jones interactions.
The Coulomb+LJ interaction is quite fast as it is written in Julia.
Nonbonded Theory can be used in geometry optimizations provided that the MM atoms are always frozen.
It can be used as MM code in QM/MM (:doc:`QM-MM`) theory objects.

The interface to the external OpenMM code has recently become available. OpenMM is a fast C++ code with a Python API
that ASH is interfaced to. OpenMM can calculate both bonded and nonbonded interactions, both periodic and nonperiodic systems
and can read in multiple types of forcefield files.


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
    settings_ash.init() #initialize

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
    settings_ash.init() #initialize

    HF_frag=Fragment(xyzfile="hf.xyz")
    #Defining atomtypes 'HT' for hydrogen and atomtype 'FX' for fluorine. These atomtypes can be named anything.
    atomtypes =['HT', 'FX']
    MM_forcefield=MMforcefield_read('forcefield.ff')

    MMobject = NonBondedTheory(forcefield=MM_forcefield, atomtypes=atomtypes, LJcombrule='geometric')

    Singlepoint(fragment=HF_frag,theory=MMobject)

###########################
OpenMM interface
###########################

The OpenMMTheory class:

.. code-block:: python

    class OpenMMTheory:
        def __init__(self, pdbfile=None, platform='CPU', active_atoms=None, frozen_atoms=None,
                     CHARMMfiles=False, psffile=None, charmmtopfile=None, charmmprmfile=None,
                     GROMACSfiles=False, gromacstopfile=None, grofile=None, gromacstopdir=None,
                     Amberfiles=False, amberprmtopfile=None, printlevel=2, nprocs=1,
                     xmlfile=None, periodic=False, periodic_cell_dimensions=None, customnonbondedforce=False):


Example creation of an OpenMMtheory object with CHARMM-files:

.. code-block:: python

    forcefielddir="/path/to/dir"
    topfile=forcefielddir+"top_all36_prot.rtf"
    parfile=forcefielddir+"par_all36_prot.prm"
    psffile=forcefielddir+"new-XPLOR-psffile.psf"
    openmmobject = OpenMMTheory(CHARMMfiles=True, psffile=psffile, charmmtopfile=topfile,
                               charmmprmfile=parfile)

Example creation of an OpenMMtheory object with GROMACS-files:

.. code-block:: python

    openmmobject = OpenMMTheory(GROMACSfiles=True, gromacstopdir="/path/to/gromacstopdir",
                    gromacstopfile="gromacstopfile.top", grofile="grofile.gro")

Example creation of an OpenMMtheory object with AMBER files:

.. code-block:: python

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile="/path/to/amberprmtopfile")

Example creation of an OpenMMtheory object with OpenMM XML file:

.. code-block:: python

    openmmobject = OpenMMTheory(xmlfile="exampl.xml")


An openmmtheory object can then be used to create a QM/MM theory object. See :doc:`QM-MM` page.