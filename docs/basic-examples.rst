==========================
Basic examples
==========================

This page shows the basics of ASH in a simple way.

#########################################
Defining coordinates in different ways
#########################################

Molecular coordinates can be read into ASH in many different way.
The key thing is to create an ASH fragment object as this is how ASH 

*Defining an atom:*

We can define an atom in ASH in a very simple way by using the atom keyword and providing an element symbol when creating an ASH Fragment.
ASH will by default create Cartesian coordinates at position: 0.0 0.0 0.0 (assuming Å).

.. code-block:: python

    from ash import *

    #Create oxygen atom fragment
    O=Fragment(atom="O")

Typically we want to also define the charge and spin multiplicity of the fragment. Here defining a neutral oxygen atom in the triplet ground state.

.. code-block:: python

    from ash import *

    #Create oxygen atom fragment
    O=Fragment(atom="O", charge=0, mult=3)

*Defining a diatomic molecule:*

Similarly we can create a diatomic molecule using the diatomic keyword (expecting the two element symbols) and by providing the diatomic_bondlength keyword as well (bondlength in Å).
ASH will by default put the first atom at position 0.0 0.0 0.0 and the next at coordinate 0.0 0.0 X (where X is equal to the diatomic_bondlength value) 



.. code-block:: python

    from ash import *

    #Create HCl diatomic fragment
    HCl=Fragment(diatomic="HCl", diatomic_bondlength="1.3", charge=0, mult=1)

*Defining a water molecule from a multi-line string of Cartesian coordinates (Å):*

.. code-block:: python

    from ash import *

    #Create H2O fragment from string
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
A few different job examples on H2O
#########################################

*Single-point calculation at the DFT-level (BP86/def2-SVP) using ORCA where the charge/mult was defined as part of the fragment:*

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
