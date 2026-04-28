Force Fields in ASH
======================================

ASH uses the OpenMM library for almost all forcefield functionality.
OpenMM can use many different forcefields and can even read in forcefields in older Amber and CHARMM file formats.
For setting up a new system it is usually recommended, however, to define the forcefields via the XML format
as this allows the most flexibility.

This page lists the various options for how to go about selecting, creating or defining forcefield using XML-files.

################################################################################################
Creating an OpenMMTheory object from XML-files
################################################################################################

When using the XML-file option one defines the OpenMMTheory ASH object by providing a PDB-file and 1 or more XML-file.
The PDB-file serves the purpose of defining the topology of the system.
It's important to note that coordinates present in the topology PDB-file will not be used by ASH unless requested separately
(coordinates always come from the Fragment object).

**Basic example:**

.. code-block:: python

    # File example.xml and file.pdb needs to be present in dir
    openmmobject = OpenMMTheory(pdbfile="file.pdb", xmlfiles=["example.xml"])

One can provide multiple XML-files to the *xmlfiles* keyword option. Because OpenMM has some built-in forceields as XML-files (e.g. charmm36.xml) it is possible to
specify them directly by name. No need to give their absolute paths or copy them to directory. See more below on built-in forcefield options.

.. code-block:: python

    # Here only file.pdb and specialresidue.xml need to be present in dir
    openmmobject = OpenMMTheory(pdbfile="file.pdb", xmlfiles=["charmm36.xml", "charmm36/water.xml", "specialresidue.xml"]) 

################################################################################################
The XML format
################################################################################################

The XML forcefield file can be defined entirely from scratch. 
See `OpenMM Creating Force Fields documentation <https://docs.openmm.org/7.6.0/userguide/application/05_creating_ffs.html>`_ 
for official information on the format and all the options.

For a brand new forcefield of a system, the file should contain a definition of the atom types, 
definition of residues including bonding information and finally the bonded and nonbonded parameters.
The forcefield information can be split into different files.

The XML-file below defines a simple forcefield for water.

.. code-block:: text

    <?xml version="1.0" encoding="utf-8"?>
    <ForceField>
    <AtomTypes>
    <Type name="OT" class="OT" element="O" mass="18.00000"/>
    <Type name="HT" class="HT" element="H" mass="1.007947"/>
    </AtomTypes>
    <Residues>

    <Residue name="H2O">
        <Atom name="O_1" type="OT"/>
        <Atom name="H_1" type="HT"/>
        <Atom name="H_2" type="HT"/>
        <Bond atomName1="O_1" atomName2="H_1"/>
        <Bond atomName1="O_1" atomName2="H_2"/>
    </Residue>
    </Residues>
    # Bonded forces
    <HarmonicBondForce>
    <Bond type1="OT" type2="HT" length="0.09572" k="462750.4"/>
    </HarmonicBondForce>
    # Nonbonded forces
    <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
        <Atom type="OT" charge="-0.834" sigma="0.31507524065751241" epsilon="0.635968"/>
        <Atom type="HT" charge="0.417" sigma="1" epsilon="0"/>
    </NonbondedForce>
    </ForceField>

To use this forcefield for water we would just need to create a PDB-file that defines the topology (what atoms are present and how are they bonded):

water.pdb:

.. code-block:: text

    REMARK   1 CREATED WITH OPENMM 8.3, 2026-04-28
    HETATM    1  O   MOL A   1       0.225  -0.836  -1.476  1.00  0.00           O
    HETATM    2  H1  MOL A   1       0.225  -0.065  -0.907  1.00  0.00           H
    HETATM    3  H2  MOL A   1       0.225  -1.607  -0.907  1.00  0.00           H
    TER       4      MOL A   1
    CONECT    1    2    3    2    3
    CONECT    2    1    1
    CONECT    3    1    1
    END

This file can be conveniently created by an ASH script like this:

.. code-block:: python

    from ash import*
    # Create fragment from XYZ-file
    frag = Fragment(xyzfile="h2o.xyz")
    # Write PDBfile of system using OpenMM
    pdbfile =frag.write_pdbfile_openmm()



We can then define the OpenMMTheory object and run an MM geometry optimization like this:

.. code-block:: python

    from ash import*

    frag = Fragment(databasefile="h2o.xyz")
    theory = OpenMMTheory(xmlfiles=["h2o.xml"], pdbfile="water.pdb", 
                autoconstraints=None, rigidwater=False)

    Optimizer(theory=theory, fragment=frag)


################################################################################################
Matching topology and forceield
################################################################################################

In order for OpenMM to match the forcefield against the system, it needs both the forcefield (via the XML-file)
and the PDB topology which must match.

Most common reasons for OpenMM not matching topology and forcefield are:

- **M**issing atoms:** if the residue definition in the XML-file is missing e.g. an H-atom present in PDB-file or other way around.
- **Missing CONECT lines in PDB-file:** If residue definition contains bond-information (e.g. <Bond atomName1="O_1" atomName2="H_1"/>)  then CONECT lines must be present in PDB-file. It is easiest to use the frag.write_pdbfile_openmm() approach above to get this right.
- **Wrong CONECT lines in PDB-file:** If you define a purely nonbonded forcefield in XML-file (i.e. no <Bond> lines) then no CONECT lines should be present in PDB-file. Make sure <Bond> lines in XML-file and CONECT lines in PDB-file match in terms of connectivity.
- **Missing element information:** PDB-file usually has to  contain an element column in order to match correctly.




################################################################################################
OpenMM built-in forcefields
################################################################################################

OpenMM has built-in various forcefields in XML format, for proteins, nucleic acids, solvents, ions etc.
Inspecting these files can be helpful when creating your own XML-files for small molecules or when you want to know what parameters are available or what OpenMM actually is using.
The OpenMM built-in Amber and CHARMM XML forcefield files are available in your environment and will be automatically found by ASH (if OpenMM was installed correctly) when you select e.g.:

.. code-block:: python

  OpenMMTheory(...,xmlfiles=["amber14/tip3p_standard.xml"])
  OpenMMTheory(...,xmlfiles=["charmm36/water.xml"])
  OpenMMTheory(...,xmlfiles=["charmm36.xml", "charmm36/water.xml"])
  OpenMMTheory(...,xmlfiles=["amber14-all.xml", "implicit/obc2.xml", "gaff_ligand.xml"])

All of these files and others can also be inspected on your system like this:

.. code-block:: shell

    # Go into OpenMM data directory within conda environment (make sure the environment is loaded)
    cd $(dirname $(which test-openmm-platforms))/../lib/python3.11/site-packages/openmm/app/data/
    ls
    # Result:
    DLPC.pdb              amber03.xml           amber99Test.xml       amoeba2013.xml        glycam-hydrogens.xml  spce.pdb              tip4pew.pdb
    DLPE.pdb              amber03_obc.xml       amber99_obc.xml       amoeba2013_gk.xml     hydrogens.xml         spce.xml              tip4pew.xml
    DMPC.pdb              amber10.xml           amber99sb.xml         amoeba2018.xml        iamoeba.xml           swm4ndp.pdb           tip4pfb.xml
    DOPC.pdb              amber10_obc.xml       amber99sbildn.xml     amoeba2018_gk.xml     implicit              swm4ndp.xml           tip5p.pdb
    DPPC.pdb              amber14               amber99sbnmr.xml      charmm36              opc.xml               test.pdb              tip5p.xml
    POPC.pdb              amber14-all.xml       amberfb15.xml         charmm36.xml          opc3.xml              tip3p.pdb
    POPE.pdb              amber96.xml           amoeba2009.xml        charmm_polar_2013.xml pdbNames.xml          tip3p.xml
    absinth.xml           amber96_obc.xml       amoeba2009_gk.xml     charmm_polar_2019.xml residues.xml          tip3pfb.xml



################################################################################################
UFF
################################################################################################

ASH features an implementation of the Universal Force Field (UFF) as an XML-file.
The XML-file is present in : /YOUR/ASH/PATH/ash/databases/forcefields/uff_mod.xml

To use UFF for your system, you have to match the topology and forcefield.
This requires defining the residue in XML-format and choose the desired UFF atomtypes.
The UFF atomtypes as defined in ASH are in: /YOUR/ASH/PATH/ash/databases/forcefields/uff_mod.xml


WARNING: The UFF implementation in ASH is not identical to other programs...

--------------------
Manual definition
--------------------

For a small molecule it is fairly straightforward to set things up. Here we will set up a UFF calculation of methane as an example.
One needs a PDB-file, created manually or perhaps by another program, that will serve the primary purpose of topology.
In the PDB-file then one needs to be sure that the molecule is part of the same residue (here  MOL), same segment (here A) and same residue-ID (her e 1) and that a
column with element names (here C, H etc.) is present. Finally, one needs to specify what atoms are connected to each other by adding CONECT statements at the end.
Note that the atomnames (C_3, H_1 etc. ; these are not atomtypes) or the residue-name (MOL) are actually not important since OpenMM
matches based on connectivity.

methane_uff.pdb:

.. code-block:: text

    REMARK  Methane with UFF atom type names
    HETATM    1  C_3 MOL A   1       0.000   0.000   0.000  1.00  0.00           C
    HETATM    2  H_1 MOL A   1       0.630   0.630   0.630  1.00  0.00           H
    HETATM    3  H_2 MOL A   1      -0.330  -0.630   0.630  1.00  0.00           H
    HETATM    4  H_3 MOL A   1      -0.630   0.830  -0.630  1.00  0.00           H
    HETATM    5  H_4 MOL A   1       0.630  -0.630  -0.930  1.00  0.00           H
    CONECT    1    2    3    4    5
    CONECT    2    1
    CONECT    3    1
    CONECT    4    1
    CONECT    5    1
    END

Next we need to define the residue in a separate XML-file, here named uff_residues.xml.
This file can contain many residue definitions if desired.

Generally, this file needs to contain a definition of all the atoms in the residue, unique atom names and the desired UFF atom types.
Finally, we must specify what atoms are bonded to each other and this information must match the information in the CONECT lines of the PDB-file.


uff_residues.xml:

.. code-block:: text

    <?xml version="1.0" encoding="utf-8"?>
    <ForceField>
    <Residues>

    <!-- ── Methane ───────────────────────────────────────────────────────── -->
    <Residue name="CH4">
        <Atom name="C" type="C_3"/>
        <Atom name="H_1" type="H_"/>
        <Atom name="H_2" type="H_"/>
        <Atom name="H_3" type="H_"/>
        <Atom name="H_4" type="H_"/>
        <Bond atomName1="C_3" atomName2="H_1"/>
        <Bond atomName1="C_3" atomName2="H_2"/>
        <Bond atomName1="C_3" atomName2="H_3"/>
        <Bond atomName1="C_3" atomName2="H_4"/>
    </Residue>

    </Residues>
    </ForceField>

Once we have defined the residue in the file we can run:

.. code-block:: python

    from ash import *

    # Defining PDB file. Needed for topology
    pdbfile="methane_uff.pdb"
    # Here reading coordinates from PDB-file into ASH Fragment as well.
    # Note: could alternatively read in coordinates of methane from an XYZ-file instead.
    fragment=Fragment(pdbfile=pdbfile, charge=0, mult=1)


    # Create list of the XML-files (one for the UFF parameters, one for the UFF residue definition)
    uff_xml_definition = [f"{ashpath}/databases/forcefields/uff_mod.xml", "uff_residues.xml"]
    print("uff_xml_definition:", uff_xml_definition)
    # Define OpenMMTheory object specifying the system topology and FF via PDB-file and XML-files.
    theory = OpenMMTheory(pdbfile=pdbfile, xmlfiles=uff_xml_definition,
                autoconstraints=None, rigidwater=False)

    # Performing a geometry optimization to test things out
    Optimizer(theory=theory, fragment=fragment)

-----------------------
Automatic definition
-----------------------


Since only a PDB-file and residue definition is needed for using UFF on any molecule, one could try to automate this.
ASH contains a function **define_uff** that, given an ASH fragment will create both the UFF-residue XML file and PDB-file topology.
It will guess the atom types based on a bond-order matrix provided by an xTB-calculation (requires tblite library to be installed).

This information can then be directly fed to OpenMMTheory.
WARNING: this automatic residue definition should be checked carefully and the residue-file (uff_residues.xml) 
has to be manually modified if the guess is wrong.


.. code-block:: python

    from ash import *

    # Define ASH fragment from e.g. an XYZ-file
    fragment=Fragment(xyzfile="file.xyz, charge=0, mult=1)

    # Automatically create the UFF XML definition and PDB-file topology
    uff_xml_definition, pdbfile = define_uff(fragment)
    print("uff_xml_definition:", uff_xml_definition)
    print("pdbfile:", pdbfile)
    # Define OpenMMTheory object specifying the system topology and FF via PDB-file and XML-files.
    theory = OpenMMTheory(pdbfile=pdbfile, xmlfiles=uff_xml_definition,
                autoconstraints=None, rigidwater=False)

    # Performing a geometry optimization to test things out
    Optimizer(theory=theory, fragment=fragment)


#######################################################
Automatic forcefield fitting (Seminario method)
#######################################################

SOON TO BECOME AVAILABLE


#######################################################
Defining forcefield for ligand / small molecule
#######################################################

Often one wants to perform a classical or QM/MM simulation of a small molecule in solution (either as part of a biomolecular system or on its own)
but one lacks forcefield parameters to do so. One has a feew options for how to proceed in this case:

- OpenMM built-in forcefields (see above). Mostly limited to biomolecular systems.
- UFF (see above).
- Create a nonbonded forcefield (charges and Lennard-Jones parameters) for the small molecule.
- Create a full forcefield for the small molecule (bonded and nonbonded parameters).

The nonbonded option is sufficient if one primarily intends to perform QM/MM simulations where the molecule will always be in the QM-region.
This may also be the only easy option if the molecule is inorganic (e.g. a metal complex) where forcefield parameterization is less straightforward. 
The nonbonded forcefield can also be used in classical simulation if one makes sure the ligand is rigid (all bonds constrained, possibly angles and dihedrals as well).
See next section below: **write_nonbonded_FF_for_ligand**

The second option (full forcefield) is generally better and is required if one wants to perform classical MM simulations where the molecule is flexible.
ASH features a function (**small_molecule_parameterizer**) that allows one to expedite this process with the help of the `openmm-forcefields <https://github.com/openmm/openmmforcefields>`_, 
that provides a convenient way of getting forcefield parameters from the `GAFF <https://ambermd.org/antechamber/gaff.html>`_ and `OpenFF <https://openforcefield.org>`_ projects. 
The limitation is that this option is primarily available for organic or drug-like molecules.
Additionally these small-molecule forcefields are intended to be only used together with Amber biomolecular forcefields (if your system also includes protein/DNA).


-------------------------------
write_nonbonded_FF_for_ligand
-------------------------------

.. code-block:: python

  def write_nonbonded_FF_for_ligand(fragment=None, charge=None, mult=None, coulomb14scale=1.0, lj14scale=1.0, 
    ff_type="CHARMM", charge_model="CM5", theory=None, LJ_model="UFF", resname="LIG", numcores=1):


ASH features a function (**write_nonbonded_FF_for_ligand**) that allows one to quickly create an OpenMM-style XML forcefield file for any ligand/molecule
with only nonbonded parameters specified which can be sufficient for QM/MM simulations or classical simulations where the ligand/molecule is rigid (all bonds constrained).

One can choose to derive the atom charges from either an xTB-calculation (using the xTB interface) or a DFT-calculation (ORCA interface).
The charge_model options are: CM5 charges or DDEC3/DDEC6 charges (requires DDEC3/DDEC6).
The Lennard-Jones parameters can either come from UFF (very crude: element-specific LJ parameters) or via DDEC3/DDEC6 population analysis.


.. warning:: It is up to you the user to make sure that the nonbonded parameters from this procedure are sensible and compatible with other molecules present in your system (described by another forcefield).
  You may have to change the parameters manually 

*Example:*

.. code-block:: text

    from ash import *

    frag=Fragment(xyzfile="ligand.xyz")

    #Script to get nonbonded model parameters for a ligand
    orcatheory=ORCATheory(orcasimpleinput="!r2scan ZORA ZORA-def2-TZVP tightscf CPCM", numcores=8)

    write_nonbonded_FF_for_ligand(fragment=frag, resname="MCMtest", charge=0, mult=1,
        coulomb14scale=1.0, lj14scale=1.0, charge_model="CM5_ORCA", theory=orcatheory, LJ_model="UFF", ff_type="CHARMM")


**Options:**

- charge_model: Options are 'CM5', 'xTB', 'DDEC3', 'DDEC6'
- LJ_model: Options are 'UFF', 'DDEC3', 'DDEC6'
- The ff_type keyword (options: 'CHARMM', 'AMBER', 'None'), writes the forcefield file so that it is compatible with the CHARMM, Amber biomolecular forcefields. Choose 'None' if not needed.
- coulomb14scale and lj14scale parameters can be changed, depending on what other forcefield this ligand-forcefield will be combined with  (OpenMM requires compatibility)

**NOTES**

- Parameters will be derived for each atom in the XYZ-file. Symmetry is currently not incorporated and this means that very 
  similar atoms in the structure will have their own charge/LJ parameters. Since this is not always desired, the user
  should take care to combine and symmetrize the parameters in the XML-file manually.
- For a ligand bound to the protein, special care must be taken. Charges are best derived from a ligand structure with all metal ions
  coordinated (e.g. including an amino acid side chain) but then the calculation will contain those extra atoms.
  This requires manual tweaking of the final charges (make sure that the sum of atom charges add up to the correct total charge).
- DDEC3/DDEC6: Both atom charges and LJ parameters can be determined from a DFT-calculation and a DDEC3/DDEC6 population analysis using the Chargemodel. This options has not been well tested and requires external programs (Chargemol and mol2aim)


-------------------------------
small_molecule_parameterizer
-------------------------------

.. code-block:: python

  def small_molecule_parameterizer(charge=None, xyzfile=None, pdbfile=None, molfile=None, sdffile=None, 
                                  smiles_string=None, resname="LIG", forcefield_option='GAFF', 
                                  gaffversion='gaff-2.11', openff_file="openff-2.0.0.offxml",
                                  expected_coul14=0.8333333333333334, expected_lj14=0.5, allow_undefined_stereo=None):

**small_molecule_parameterizer** allows you to quickly create an OpenMM XML forcefield file with bonded and nonbonded parameters for your molecule.
You can choose between two general forcefields: `GAFF <https://ambermd.org/antechamber/gaff.html>`_  or `OpenFF <https://openforcefield.org>`_. 
Different GAFF and OpenFF versions are also available. The limitation is that creating the small-molecule forcefield from these general forcefields can only be done for "organic" chemical elements (H,C,N,O,S,P,F,Cl,Br,I; also ions such as 
Li+, Na+, K+, Rb+, F-, Cl-, Br-, and I-).
These small-molecule forcefields are intended to be only used together with Amber biomolecular forcefields (if your system also includes protein/DNA).

The program depends on a few Python libraries that have to be installed when prompted.
It should be enough to install `openmmforcefields <https://github.com/openmm/openmmforcefields>`_ as it will automatically install: `openff-toolkit <https://github.com/openforcefield/openff-toolkit>`_ , `RDKit <https://github.com/rdkit/rdkit>`_, `parmed <https://github.com/ParmEd/ParmEd>`_.
ASH will tell you which libraries are missing and how to install them when you try to use the function.
Specifically we use the OpenFF toolkit to create a Molecule object (see `OpenFF Molecule <https://open-forcefield-toolkit.readthedocs.io/en/0.9.2/api/generated/openff.toolkit.topology.Molecule.html>`_ and `Molecule Cookbook <https://docs.openforcefield.org/projects/toolkit/en/stable/users/molecule_cookbook.html>`_) .

**small_molecule_parameterizer** is very easy to use most of the time.
You simply need to provide molecular structure information in the form of either an XYZ-file, PDB-file, a `SMILES string <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system>`_ , MDL Mol-file or SDF-file.
Additionally the total charge of the molecule needs to be specified.

There are cases where parsing the molecular information from a coordinate-file fails and you may have to input a SMILES-string directly.
See `Molecule Cookbook <https://docs.openforcefield.org/projects/toolkit/en/stable/users/molecule_cookbook.html>`_  and `SMILES tutorial: <http://hjkgrp.mit.edu/tutorials/2013-10-29-geometries-strings-smiles-and-openbabel>`_ . 


*Example using OpenFF and XYZ-file*

.. code-block:: python

  from ash import *
  #Creating forcefield for nitrate using OpenFF. Here providing xyz-file as input
  small_molecule_parameterizer(forcefield_option="OpenFF", xyzfile="no3.xyz", charge=-1)

The output is an XML-file that can then be used as input to **OpenMMTheory**, **OpenMM_Modeller** or **solvate_small_molecule** functions (see below).
Additionally a PDB-file is written out for convenience (matches information in XML-file).

*Example using GAFF and SMILES string*

.. code-block:: python

  from ash import *
  #Creating forcefield for nitrate using GAFF. Here providing a SMILES string as input
  small_molecule_parameterizer(forcefield_option="GAFF", smiles_string="[N+](=O)([O-])[O-]")
  #Note: no PDB-file will be created in this case.


.. warning:: The XML-file created by this function will contain bonded parameters and it is thus important that the topology of the molecule is available when using the XML-file
  together with OpenMM. Otherwise, the pairing of molecule and small-molecule forcefield in the XML-file will not work. As OpenMM will typically get the topology from a PDB-file you must ensure 
  to have a PDB-file that contains CONECT lines at the bottom of the PDB-file that describes the connectivity of the small molecule. A PDB-file with connectivity is automatically created if you read in an XYZ-file
  to small_molecule_parameterizer above. You can also use the  **xyz_to_pdb_with_connectivity** function.


The following error can sometimes occur: 

.. code-block:: text

  ValueError: Final molecular charge does not match input; could not find valid bond ordering

This means that RDKit failed to understand the bonds in the molecule. Often this occurs if the charge of the molecule is wrong.


See also :doc:`protein_ligand_binding` for a demonstration on using the **small_molecule_parameterizer** for setting up a protein-ligand complex.