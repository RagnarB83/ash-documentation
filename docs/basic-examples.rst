==========================
Basic examples
==========================


#########################################
Defining coordinates in different ways
#########################################

*Defining an atom:*

.. code-block:: python

    from ash import *

    #Create oxygen atom fragment
    O=Fragment(atom="O")

*Defining a diatomic molecule:*

.. code-block:: python

    from ash import *

    #Create HCl diatomic fragment
    HCl=Fragment(diatomic="HCl", diatomic_bondlength="1.3")

*Defining an H2O molecule from a multi-line string of Cartesian coordinates (Å):*

.. code-block:: python

    from ash import *

    #Create H2O fragment from string
    coords="""
    O       -1.377626260      0.000000000     -1.740199718
    H       -1.377626260      0.759337000     -1.144156718
    H       -1.377626260     -0.759337000     -1.144156718
    """
    H2O=Fragment(coordsstring=coords)

*Defining an H2O molecule from an Xmol xyz-file containing Cartesian coordinates (Å):*

.. code-block:: python

    from ash import *

    #Create H2O fragment from xyz-file
    H2O=Fragment(xyzfile="h2o.xyz)

where h2o.xyz must be present in working directory and contain:

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

Single-point calculation at the DFT-level (BP86/def2-SVP) using ORCA:

.. code-block:: python

    from ash import *

    #Create H2O fragment
    coords="""
    O       -1.377626260      0.000000000     -1.740199718
    H       -1.377626260      0.759337000     -1.144156718
    H       -1.377626260     -0.759337000     -1.144156718
    """
    H2O=Fragment(coordsstring=coords)

    #Defining ORCA-related variables
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    ORCAcalc = ORCATheory(orcadir='/opt/orca_5.0.2', charge=0, mult=1, 
                            orcasimpleinput=orcasimpleinput)

    #Single-point energy job on H2O with ORCAcalc theory
    energy = Singlepoint(fragment=H2O, theory=ORCAcalc)



Geometry optimization at the DFT-level (BP86/def2-SVP) using ORCA:

.. code-block:: python

    from ash import *

    #Create H2O fragment
    H2O=Fragment(xyzfile="h2o.xyz")

    #Defining ORCA-related variables
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    ORCAcalc = ORCATheory(orcadir='/opt/orca_5.0.2', charge=0, mult=1, 
                            orcasimpleinput=orcasimpleinput)

    #Geometry optimization on H2O with ORCAcalc theory
    geomeTRICOptimizer(fragment=H2O, theory=ORCAcalc)
