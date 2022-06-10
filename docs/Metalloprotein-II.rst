Metalloprotein tutorial II: Ferredoxin
======================================

Here we encounter a slightly more complex metalloprotein example, ferredoxin.
This protein setup require slightly more work due to errors arising parsing the more complicated PDB-file with a slightly more complicated metal cluster.

Ferredoxin contains 2 copies of the [2Fe-2S] cofactor bound to 4 deprotonated cysteines. 
Additionally, the X-ray structure contains additional molecules: SCN\ :sup:`-` \ and SO\ :sub:`4`:sup:`2-` \ anions and a benzamidine molecule. 


######################################################
1. OpenMM_Modeller parsing errors
######################################################

If we download the 6lk1.pdb file and read into OpenMM_Modeller:

.. code-block:: python

    from ash import *
    #Original raw PDB-file (no hydrogens, nosolvent)
    pdbfile="6lk1.pdb"
    # Setting up system via Modeller
    OpenMM_Modeller(pdbfile=pdbfile,forcefield="CHARMM36")

we will get this error:

.. code-block:: text

    ValueError: No template found for residue 191 (FES).  This might mean your input topology is missing some atoms or bonds, 
    or possibly that you are using the wrong force field.

This arises simply because the FES residue, i.e. the [2Fe-2S] cluster is not present in the CHARMM36 forcefield.
In the PDB-file the residue is defined like this:

.. code-block:: text

    HETATM 1457 FE1  FES A 101       5.508  -6.531  22.895  1.00  7.55          FE
    HETATM 1458 FE2  FES A 101       4.734  -4.618  21.086  1.00  6.50          FE
    HETATM 1459  S1  FES A 101       5.095  -4.352  23.250  1.00  7.17           S
    HETATM 1460  S2  FES A 101       4.997  -6.752  20.748  1.00  7.46           S


If we tell OpenMM_Modeller about the extra residue: 

.. code-block:: python

    from ash import *

    #Original raw PDB-file (no hydrogens, nosolvent)
    pdbfile="6lk1.pdb"
    #XML-file to deal with cofactor
    extraxmlfile="specialresidue.xml"
    # Setting up system via Modeller
    OpenMM_Modeller(pdbfile=pdbfile,forcefield="CHARMM36", extraxmlfile="specialresidue.xml")

where specialresidue.xml contains:

.. code-block:: python

    <ForceField>
    <AtomTypes>
    <Type name="FEX" class="Fe" element="Fe" mass="55.84700"/>
    <Type name="SXM" class="S" element="S" mass="32.065"/>
    </AtomTypes>
    <Residues>
    <Residue name="FES">
    <Atom name="FE1" type="FEX"/>
    <Atom name="FE2" type="FEX"/>
    <Atom name="S1" type="SXM"/>
    <Atom name="S2" type="SXM"/>
    </Residue>
    </Residues>
    <NonbondedForce coulomb14scale="1.0" lj14scale="1.0">
    <Atom type="FEX" charge="0.0" sigma="1.3" epsilon="0.0"/>
    <Atom type="SXM" charge="0.0" sigma="1.3" epsilon="0.0"/>
    </NonbondedForce>
    <LennardJonesForce lj14scale="1.0">
    <Atom type="FEX" sigma="0.3" epsilon="0.00000"/>
    <Atom type="SXM" sigma="0.3" epsilon="0.00000"/>
    </LennardJonesForce>
    </ForceField>

.. warning:: For OpenMM to correctly parse the specialresidue.xml file, it is important that the PDB-file contains element definitions (column 77-78) for
    each element of the special residue and the atom names in the XML file much match the atom names in the PDB-file.

and run this script we get instead another error:

.. code-block:: text

    ValueError: No template found for residue 191 (FES).  The set of atoms matches FES, but the bonds are different.  
    Perhaps the chain is missing a terminal group?

This error messages suggest that OpenMM Modeller recognizes our residue definition for FES but is confused about bonding. This arises due to the presence of bonding information in the bottom of the PDB-file
in the form of CONE lines (each line indicates which atom indices should be considered to have a bond between them):

.. code-block:: text

    CONECT  299 1457
    CONECT  338 1457
    CONECT  356 1458
    CONECT  575 1458
    CONECT 1033 1523
    CONECT 1064 1523
    CONECT 1082 1524
    CONECT 1298 1524
    CONECT 1457  299  338 1459 1460
    CONECT 1458  356  575 1459 1460
    CONECT 1459 1457 1458


While we could add bonding information to specialresidue.xml and try to match the connectivity in the PDB-file an easier solution is to remove the connectivity information by creating a modified version
of the PDB-file. This should probably always be an acceptable solution since we will constrain our metal-cluster residue anyway, as we don't have forcefield parameters available.
Thus we make a modified version, called 6lk1-mod.pdb, that does not contain the CONE lines and we also remove most header lines of the PDB-file (all lines before ATOM/HETATM section begins).

Next we run our script again (now using 6lk1-mod.pdb as input PDB file):

.. code-block:: python

    from ash import *

    #Original raw PDB-file (no hydrogens, nosolvent)
    pdbfile="6lk1-mod.pdb"
    #XML-file to deal with cofactor
    extraxmlfile="specialresidue.xml"
    # Setting up system via Modeller
    OpenMM_Modeller(pdbfile=pdbfile,forcefield="CHARMM36", extraxmlfile="specialresidue.xml")


When we run this script we instead get a different error message:

.. code-block:: text

    ValueError: No template found for residue 192 (BEN).  The set of atoms is similar to INDA, but it is missing 6 hydrogen atoms.

This message refers to the fact that the PDB-file contains BEN residues (benzamidine) that are neither present in the CHARMM36 protein forcefield or in our specialresidue.xml file.
As benzamidine is simply a crystallized molecule from the buffer solution and has little do with the protein, we here make the choice to simply remove the BEN residues from 6lk1-mod.pdb.

We then get similar messages associated with missing residue definitions for SCN and SO4

.. code-block:: text

    ValueError: No template found for residue 191 (SCN).  The set of atoms is similar to THAZ, but it is missing 5 atoms.

    ValueError: No template found for residue 191 (SO4).  The set of atoms is similar to MSO4, but it is missing 4 atoms.

and again we make the choice to remove these crystallized contaminants from 6lk1-mod.pdb.
Once we have done this, OpenMM_Modeller proceeds without problems but this does not mean of course that the system is correctly set up.

######################################################
2. OpenMM residue variants: protonation states
######################################################

As previously occurred for rubredoxin, OpenMM Modeller protonates the cysteine residues that are coordinated to the Fe ions.
Since we want to avoid this, we again define a dictionary with information about abnormal residues and pass this on to OpenMM_Modeller.
Since the protein contains two protein chains (named 'A' and 'B' in the PDB-file) with the [2Fe-2S] cofactor coordinates to 4 cysteines in each chain,
we need to define these cysteines as deprotonated ('CYX' label)

.. code-block:: python

    from ash import *

    #Original raw PDB-file (no hydrogens, nosolvent)
    pdbfile="6lk1-mod.pdb"
    #XML-file to deal with cofactor
    extraxmlfile="specialresidue.xml"

    #Defining deptonated cysteine residues
    residue_variants={'A':{5:'CYX',8:'CYX',38:'CYX',41:'CYX'},'B':{5:'CYX',8:'CYX',38:'CYX',41:'CYX'}}
    # Setting up system via Modeller
    OpenMM_Modeller(pdbfile=pdbfile,forcefield="CHARMM36", extraxmlfile="specialresidue.xml", residue_variants=residue_variants)


The printed table shows what Cys residues we selected to deprotonate:

.. code-block:: text

    User defined residue variants per chain:
    Chain A : {6: 'CYX', 9: 'CYX', 39: 'CYX', 42: 'CYX'}

    MODELLER TOPOLOGY - RESIDUES TABLE

    ASH-resid   Resname      Chain-index  Chain-name   ResID-in-chain       User-modification
    ----------------------------------------------------------------------------------------------------
    0           MET          0            A            1
    1           ASP          0            A            2
    2           ILE          0            A            3
    3           TYR          0            A            4
    4           VAL          0            A            5
    5           CYS          0            A            6                   -- This residue will be changed to: CYX --
    6           THR          0            A            7
    7           VAL          0            A            8
    8           CYS          0            A            9                   -- This residue will be changed to: CYX --
    9           GLY          0            A            10
    10          TYR          0            A            11
    11          GLU          0            A            12
    12          TYR          0            A            13
    ...



###########################################################
3. A more realistic nonbonded model for the [2Fe-2S] 
###########################################################

While a pragmatic solution to dealing with simple inorganic residues like the [2Fe-2S] cluster is to simply create 
dummy forcefield parameters as in the specialresidue.xml file above, this will not always work.
If the charges of the Fe and S atoms in [2Fe-2S] are zero, then this means no electrostatic interaction is present between
these atoms and the rest of the protein+solvent. Furthermore, with epsilon and/or sigma parameters being 0.0 no repulsion (or attractive dispersion)
forces are present between [2Fe-2S] and other atoms, meaning that other atoms could occupy the same space as the [2Fe-2S] cluster.

.. note:: In electrostatically embedded QM/MM the metal cluster will most often be in the QM-region and any atom charges defined for the cluster will not be used.
    Note, however, that the LJ interactions between QM and MM atoms are calculated and the LJ parameters may be important.


Thus a more realistic scenario is to come up with a proper nonbonded model for the [2Fe-2S] cluster: i.e. charges and Lennard-Jones parameters.
There are two main choices here:
1. Search the literature for a study using nonbonded MM parameters for the same/similar residue. Ideally with the same protein forcefield.
2. Derive the parameters using similar residues already present in the forcefield.
3. Derive the parameters from a DFT calculation and a population analysis.

Option 3 is the most general solution. 
----THIS IS NOT YET FINISHED----