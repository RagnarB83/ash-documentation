==========================
Basic examples
==========================

This page shows the basics of using ASH via examples.

#########################################
Python basics
#########################################

An ASH inputfile is a Python script/program with ASH loaded as a library that gives access to specific ASH functionality (the Fragment, ORCATheory, Singlepoint, geomeTRICOptimizer etc. names are ASH object, only available once ASH has been imported).
This means that ASH does not behave like a regular QM program with a specific inputfile syntax for a specific job.
Instead, an ASH inputfile is its own program (in Python), written by the user that when executed will carry out all instructions in the file.
This means that when a user's ASH script is executed, the Python interpreter will carry out all of the instructions defined in the script but nothing else.
This gives a lot of flexiblity as it allows the user to create simple or complex, specific or general, single-job or multi-job calculations depending on the need or taste.

The only drawback is that it requires the user to learn some Python programming basics in order to use the program.
There are many many resources available for learning Python.

Examples:

- https://wiki.python.org/moin/BeginnersGuide/NonProgrammers
- https://python-course.eu/python-tutorial/
- https://blog.finxter.com/python-cheat-sheet/
- https://www.freecodecamp.org/news/the-python-guide-for-beginners/


#########################################
Defining coordinates in different ways
#########################################

Molecular coordinates can be read into ASH in many different ways.
The key thing is to create an ASH fragment object as this is how ASH reads molecular information.

*Defining an atom:*

We can define an atom in ASH in a very simple way by using the atom keyword and providing an element symbol when creating an ASH Fragment.
ASH will by default create Cartesian coordinates at position: 0.0 0.0 0.0 (assuming Å).

.. code-block:: python

    from ash import *

    #Create oxygen atom fragment
    O=Fragment(atom="O")


When we run this script ASH will print output indicating that the fragment has been created:

.. code-block:: text

    --------------------------------------------------------------------------------
                                    New ASH fragment
    --------------------------------------------------------------------------------

    Creating Atom Fragment
    Creating/Updating fragment attributes...
    Number of Atoms in fragment: 1
    Formula: O1
    Label: None
    Charge: None Mult: None

    --------------------------------------------------------------------------------


Typically we want to also define the charge and spin multiplicity of the fragment. Below we define a neutral oxygen atom in the triplet ground state:

.. code-block:: python

    from ash import *

    #Create oxygen atom fragment
    O=Fragment(atom="O", charge=0, mult=3)

*Defining a diatomic molecule:*

Similarly we can create a diatomic molecule using the diatomic keyword (where a 2-element symbol string is expected) and by providing the diatomic_bondlength keyword as well (bondlength in Å).
ASH will by default put the first atom at position 0.0 0.0 0.0 and the next at coordinate 0.0 0.0 X (where X is equal to the diatomic_bondlength value) 


.. code-block:: python

    from ash import *

    #Create HCl diatomic fragment
    HCl=Fragment(diatomic="HCl", diatomic_bondlength="1.3", charge=0, mult=1)

*Defining a water molecule from a multi-line string of Cartesian coordinates (Å):*

Here we define first a Python multi-line string (the 3 quotation marks are necessary) and then use the coordsstring keyword of Fragment to point to this string. 

.. code-block:: python

    from ash import *

    #Create H2O fragment from a multi-line string
    coords="""
    O       -1.377626260      0.000000000     -1.740199718
    H       -1.377626260      0.759337000     -1.144156718
    H       -1.377626260     -0.759337000     -1.144156718
    """
    H2O=Fragment(coordsstring=coords)

*Defining a water molecule from an Xmol xyz-file containing Cartesian coordinates (Å):*


.. code-block:: python

    from ash import *

    #Create H2O fragment from xyz-file
    H2O=Fragment(xyzfile="h2o.xyz", charge=0, mult=1)

where h2o.xyz must be present in working directory and should look like (a 2-line header is always necessary containing the number of atoms in the first line):

.. code-block:: text

    2
    h2o title line
    O       -1.377626260      0.000000000     -1.740199718
    H       -1.377626260      0.759337000     -1.144156718
    H       -1.377626260     -0.759337000     -1.144156718

*Defining a protein fragment from a PDB-file:*

.. code-block:: python

    from ash import *

    #Create a protein fragment from PDB-file
    protein=Fragment(pdbfile="lysozyme.pdb")

where lysozyme.pdb must be present in working directory and be a regular PDB-file.

.. note:: When ASH creates a Fragment from a PDB-file, it will only extract element and coordinate information from the file, not atom-type or topology information. OpenMMTheory is needed for reading topology from a PDB-file.


#########################################
Defining theories 
#########################################

We can define theory levels using any theory level defined in ASH that has a valid interface to an external QM or MM program.
See :doc:`QM-interfaces`, :doc:`MM-interfaces` and :doc:`module_QM-MM`

The syntax can very different for different theory levels.

*Defining an ORCATheory level:*

.. code-block:: python

    from ash import *

    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=8)

When a Theory object is created, ASH by default prints out information on the object and may also check whether it can find the external program.
For the ORCATheory object created above, ASH would print out:

.. code-block:: text

                    #####################################
                    #                                   #
                #     ORCATheory initialization     #
                    #                                   #
                    #####################################


    Checking for ORCA location
    No orcadir argument passed to ORCATheory. Attempting to find orcadir variable in ASH settings file (~/ash_user_settings.ini)
    Found no orcadir variable in ASH settings file either.
    Checking for ORCA in PATH environment variable.
    Found orca binary in PATH. Using the following directory: /Applications/orca_5_0_3_macosx_arm64_openmpi411
    Checking if ORCA binary works... yes
    ORCA parallel job requested. Make sure that the correct OpenMPI version (for the ORCA version) is available in your environment
    OpenMPI binary directory found: /Users/bjornsson/miniconda/bin
    Testing that mpirun is executable... yes
    OpenMPI version: 4.1.1

    Creating ORCA object
    ORCA dir: /Applications/orca_5_0_3_macosx_arm64_openmpi411
    ! BP86 def2-SVP def2/J


    ORCATheory object created!

Note, however, that defining a Theory object will not result in a calculation to be carried out.

*Defining an xTBTheory level:*

.. code-block:: python

    from ash import *

    xTBcalc = xTBTheory(xtbmethod='GFN2', runmode='library')

*Defining an OpenMMTheory level:*

.. code-block:: python

    from ash import *

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile="/path/to/amberprmtopfile")



#########################################
A few different job examples on H2O
#########################################

*Single-point calculation at the DFT-level (BP86/def2-SVP) using ORCA where the charge/mult is defined as part of the fragment:*

Here is a very simple script that defines an H\ :sub:`2`\O\  fragment (called H2O) from an available h2o.xyz file, defining charge and spin multiplicity as well, next
creating the ORCATheory object (called ORCAcalc) and then calling the **Singlepoint** function that takes as input argument the ASH fragment (here H2O) and an ASH theory object (here ORCAcalc).

.. code-block:: python

    from ash import *

    #Create H2O fragment
    H2O=Fragment(xyzfile="h2o.xyz", charge=0, mult=1)

    #Defining ORCA-related variables
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput)

    #Single-point energy job on H2O with ORCAcalc theory
    energy = Singlepoint(fragment=H2O, theory=ORCAcalc)

    print("Final energy:", energy)


ASH will print information related to the creation of the H2O Fragment object and the creation of the ORCATheory object and will then run and print output related to the Singlepoint function:

.. code-block:: text

                    ################################
                    #                              #
                    #     Singlepoint function     #
                    #                              #
                    ################################


    Warning: Charge/mult was not provided to Singlepoint
    Fragment contains charge/mult information: Charge: 0 Mult: 1 Using this instead
    Make sure this is what you want!
    Doing single-point Energy job on fragment. Formula: H2O1 Label: OHH
    ------------RUNNING ORCA INTERFACE-------------
    Running ORCA object with 1 cores available
    Job label: None
    Creating inputfile: orca.inp
    ORCA input:
    ! BP86 def2-SVP def2/J tightscf



    Charge: 0  Mult: 1
    ORCA Calculation started.
    ORCA Calculation done.
    ORCA converged in 11 iterations

    ORCA energy: -76.360561445295
    Single-point ORCA energy: -76.360561445295
    ------------ENDING ORCA-INTERFACE-------------

    ------------------------------------------------------------
    Time to calculate step (ORCA run): 0.4 seconds, 0.0 minutes.
    ------------------------------------------------------------
    Energy:  -76.360561445295

    ---------------------------------------------------------------
    Time to calculate step (Singlepoint): 0.4 seconds, 0.0 minutes.
    ---------------------------------------------------------------
    Final energy: -76.360561445295







*Single-point calculation where charge/mult is given as input to the jobtype:*

If you don't define charge and multiplicity as part of the fragment (generally recommended) it is also possible to provide this information to the job-type function.

.. code-block:: python

    from ash import *

    #Create H2O fragment
    H2O=Fragment(xyzfile="h2o.xyz")

    #Defining a temporary string that will become part of the ORCA inputfile
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    #Defining the ORCATheory
    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput)

    #Single-point energy job on H2O with ORCAcalc theory with charge/mult provided
    energy = Singlepoint(fragment=H2O, theory=ORCAcalc, charge=0, mult=1)



*Geometry optimization at the DFT-level (BP86/def2-SVP) using ORCA:*

Instead of a single-point energy calculation we can run a geometry optimization instead.
To use the recommended geomeTRICOptimizer function, the geomeTRIC Python library needs to have been installed.

.. code-block:: python

    from ash import *

    #Create H2O fragment with charge/mult information
    H2O=Fragment(xyzfile="h2o.xyz", charge=0, mult=1)

    #Defining ORCA-related variables
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput)

    #Geometry optimization on H2O with ORCAcalc theory
    geomeTRICOptimizer(fragment=H2O, theory=ORCAcalc)

*Numerical frequency calculation at the DFT-level (BP86/def2-SVP) using ORCA:*

Or we can run a numerical frequency job instead.

.. code-block:: python

    from ash import *

    #Create H2O fragment with charge/mult information
    H2O=Fragment(xyzfile="h2o.xyz", charge=0, mult=1)

    #Defining ORCA-related variables
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput)

    #Numerical frequencies
    NumFreq(fragment=H2O, theory=ORCAcalc)
