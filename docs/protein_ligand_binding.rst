Modelling protein-ligand binding in ASH: Tutorial
====================================================

This is a tutorial on how to set up a model for studying organic ligand binding to a protein using ASH, 
starting with  a classical approach and ending with a QM/MM approach.
All the steps about how to set up the system, parameterize the ligand, solvate the system and run MD simulations
are performed within the ASH environment with the help of libraries such as OpenMM, OpenFF, OpenBabel, ParmEd, MDTraj, etc.
Later we will also show how to run metadynamics simulations to study the free-energy surface of the ligand binding.

The instructions below will be given in a general way with specific example shown for the trypsin-benzamidine protein-ligand complex with PDB-ID: 2OXS

######################################################
**1. Preparing initial files**
######################################################

**Protein PDB-file**

To get started we need a PDB file of the protein. This can be an initial crystal structure from the PDB database or a previous model. Hydrogens do not need to be present in the PDB file, as they will be added during the setup if missing.
You can use a PDB-file of a protein-ligand complex but you will have to prepare a version of the PDB-file that only contains the protein (no ligand).
We will use here the trypsin-benzamidine complex with PDB-ID: 2OXS. In addition to the protein, it contains benzamidine (BEN residue), a sulfate ion (SO4 residue) and a calcium ion (CA).
We will not model the sulfate ion or the calcium ion and will delete them. Additionally we split the file so that we have files for the protein separately and the ligand separately.

We will here create files with the following naming convention:

- xray_2OXS_full.pdb #This is the original PDB-file
- xray_2OXS_protein.pdb #This file contains only protein. Ligand and ion (SO4) has been manually removed
- xray_2OXS_ligand.pdb #This file contains only the ligand at it's original coordinates in X-ray structure.

**Ligand PDB-file**

The ligand in an X-ray structure is unlikely to contain H-atoms but we need them in order to create a reliable forcefield that can be used in MD simulations.
We need to create a separate PDB-file for the ligand that contains all the hydrogen atoms and has a sensible internal geometry.
This PDB-file can e.g. be created from an XYZ-file that should have an acceptable initial geometry (ideally DFT-optimized)
or from some other file format such as a MDL Mol-file or MDL SDF file. We could in principle use the original coordinates from the X-ray structure for the ligand, but it is not necessary.
This PDB-file needs to contain all the hydrogen atoms and needs to have a sensible internal geometry.
The file will be referered to as ligand.pdb here.

Note that since we need the PDB-file of the ligand to contain connectivity information (CONECT statements) we need to make sure that the PDB-file is created correctly.
Here we will use OpenBabel to create the PDB-file from either an XYZ-file, Mol-file, SDF-file or another PDB-file.

.. code-block:: python

    from ash import *

    #Example: Create PDB-file from an XYZ-file. Creates file: ligand.pdb (requires OpenBabel).
    pdbfile = xyz_to_pdb_with_connectivity("ligand.xyz")

    #Example: Create PDB-file with correct connectivity. Creates file: ligand__withcon.pdb (requires OpenBabel).
    pdbfile = writepdb_with_connectivity("ligand.pdb")

    #Other options:
    #Example: Create PDB-file from an MDL Mol-file (requires OpenBabel). Creates file: ligand.pdb
    pdbfile = mol_to_pdb("ligand.mol")
    #Example: Create PDB-file from an MDL SDF-file (requires OpenBabel). Creates file: ligand.pdb
    pdbfile = sdf_to_pdb("ligand.sdf")


######################################################
**2. Preparing ligand forcefield**
######################################################

We then need to think about the forcefield. Various protein/nucleic-acid forcefields are available in ASH (CHARMM, Amber etc.) and can be used automatically.
However, the forcefield for the ligand is the main issue as it is unlikely present in the biomolecular forcefields.
We also need to consider the compatibility between the forcefield for the ligand and the forcefield for the protein.

Here we choose to use the Amber14 forcefield for the protein and the GAFF (Generalized Amber force field) forcefield for the ligand as this can be conveniently set up using ASH.
Another option is to use one of the OpenFF forcefields for the ligands (also compatible with Amber14).

ASH features a convenient function : **small_molecule_parameterizer** that can automatically generate the forcefield for the ligand
by determining the topology of the input ligand and matching it to general parameters available for either GAFF or OpenFF.
This is made possible by functionality available in the **openmm-forcefields** package (https://github.com/openmm/openmmforcefields) 
which needs to be installed when prompted.

The **small_molecule_parameterizer** function requires in principle only the input PDB-file for the ligand.
The PDB-file will be automatically converted into a SMILES string that is then used to generate the topology and suitable parameters
are found in the GAFF (or OpenFF) forcefield. If for some reason this does not work it is also possible to provide other inputfiles to
**small_molecule_parameterizer** such as : xyzfile, molfile, sdffile. One can also provide a SMILES string (smiles_string keyword).


.. code-block:: python

    from ash import *
    #Get Amber compatible forcefield for small molecule using GAFF or OpenFF
     small_molecule_parameterizer(pdbfile="ligand.pdb", forcefield_option='GAFF', output_xmlfile="ligand.xml")


The function returns an OpenMM forcefield object (that assumes Amber14 for protein and solvent and GAFF for the ligand)
but also writes out an XML-file with the forcefield parameters for the ligand (ligand.xml). 
It is usually best to use the ligand.xml file directly.


######################################################
**3. Merge and align protein and ligand**
######################################################

Now that we have a ligand.pdb file (containing the ligand with H-atoms and correct connectivity) and a ligand.xml file (containing the forcefield parameters for the ligand)
we could in principle proceed to set up the system. However, first we need to merge the protein and ligand into one PDB-file.
If we don't care about the ligand being in a specific position w.r.t. the protein, we could simply visualize xray_2OXS_protein.pdb and ligand.pdb in e.g. VMD and make sure that protein and ligand do not clash. Otherwise modify the coordinates of the ligand in the ligand.pdb file.
This would be fine if want to initially study the unbound form of the system or possible predict binding by MD later.

However, if we want start a simulation with the ligand in the original binding site according to the X-ray structure then we have to make sure that the new hydrogenated-ligand we created is properly aligned in the protein.
This would require either modifying the coordinates of the ligand in the ligand.pdb file using a suitable visualization program (e.g. VMD), perform docking,  or alternatively we could superimpose the new hydrogenated ligand onto the original ligand-position in the X-ray structure.
Here we will show how to do the latter using ASH using the **flexible_align_pdb** function.

**Align the ligand onto the desired previous position**

.. code-block:: python

    #a. Read hydrogenated ligand PDB-file into ASH
    new_ligand_pdb="ligand.pdb"
    newligand = Fragment(pdbfile=new_ligand_pdb)
    print("New ligand coords:", newligand.print_coords())

    #b. Read ligand from a file containing only the ligand ATOM/HETATM lines from original PDB-structure (e.g. an X-ray structure with a bound-ligand)
    old_ligand_pdb="xray_2OXS_ligand.pdb" #This file should only contain the ligand. Probably missing H-atoms.
    oldligand = Fragment(pdbfile=old_ligand_pdb)
    print("Old ligand:", oldligand.print_coords())

    #c. Define the atoms in common in new and old ligand (at least carbon skeleton, all nonH-atoms should work)
    #Here defining a list of lists that contain the atom indices in new_ligand (system A) and old_ligand (systemB)
    subsetA=newligand.get_nonH_atomindices() #Getting atom indices of non-H atoms
    subsetB=oldligand.get_nonH_atomindices() #Getting atom indices of non-H atoms
    subset=[subsetA,subsetB] #Combining lists into a list-of-lists

    #d. Align new ligand (with H-atoms and matching XML-file) so that it matches (as well as possible) the position of the old-ligand atoms
    #Note: subset needs to be properly chosen. Reordering is usuaully necessary for alignment (because atom order may differ)
    newligand_aligned = flexible_align_pdb(new_ligand_pdb, old_ligand_pdb, subset=subset, reordering=True, reorder_method='brute')


**Merged protein-ligand PDB-file**

Now that we have the ligand.pdb file, oriented and aligned the way we want, we can merge protein and ligand back together in to a single PDB-file.
We can use the **merge_pdb_files** function in ASH to do this. This function is convenient as it will preserve and update the CONECT statements of the ligands which is important for the OpenMM_Modeller step later.

.. code-block:: python

    from ash import *
    protein_pdbfile="xray_2OXS_protein.pdb"
    ligand_pdbfile="ligand_aligned.pdb" #This is the aligned ligand PDB-file (i.e. having the geometry we want). Atom-order needs to match information in ligand.xml
    merged_pdbfile = merge_pdb_files(protein_pdbfile,ligand_pdbfile, outputname="2OXS_protein_ligand_merged.pdb")



######################################################
**4. Prepare system using OpenMM_Modeller**
######################################################

Now we should have a merged PDB-file (containing both protein and ligand) and a forcefield for the ligand (ligand.xml).
We can now proceed to use the **OpenMM_Modeller** function to set up the system. We use the merged protein-ligand PDB-file to define the system geometry and topology, 
we specify an Amber14 forcefield for the protein (needs to be compatible with the ligand FF), TIP3P-FB forcefield for water (compatible with Amber14) and the ligand forcefield (GAFF or OpenFF) for the ligand via the 
ligand.xml file previously created.

See :doc:`OpenMM-interface` for more information on using **OpenMM_Modeller**.

.. code-block:: python

    from ash import *

    merged_pdbfile="merged.pdb"
    #Setup system using OpenMM_Modeller using merged PDB-file
    OpenMM_Modeller(pdbfile=merged_pdbfile, forcefield="Amber14",
        extraxmlfile="ligand.xml", residue_variants={}, watermodel="tip3p-fb", pH=7.0, solvent_padding=10.0, ionicstrength=0.1)

**OpenMM_Modeller** will apply the Amber14 protein forcefield to the protein and the GAFF/OpenFF forcefield to the ligand.
Note that one must make sure that the merged PDB-file of the protein and ligand contains the correct connectivity information for the ligand (CONECT lines).
Additionally one must make sure that any residues in the protein are correctly treated (with respect to protonation states, disulfide bridges, metal ions etc.).

If the **OpenMM_Modeller** function is successful, a final PDB-file, "finalsystem.pdb" will be created that contains the solvated protein-ligand system with
protein and ligand oriented according to the initial coordinates of "merged.pdb". The coordinates in the input "merged.pdb" file 
can contain the system in either bound or unbound form and can be modified before running **OpenMM_Modeller**. 
Note that due to the present of solvent it is trickier to change the ligand position of the solvated system after the **OpenMM_Modeller** step
(would require running a biased MD simulation).

.. warning:: Make sure that the ligand geometry in the merged PDB-file matches the information in the ligand.xml file. Otherwise the ligand will not be recognized.


######################################################
**5. STEPS 1-4 COMBINED**
######################################################

Here we show a script that combines the steps 1-4 into a single ASH script that could in principle be used to conveniently perform all the steps in one go.

.. code-block:: python

    from ash import *

    #############################################################
    #1. Parameterize ligand using a hydrogenated XYZ-structure
    #############################################################
    #Here choosing GAFF
    small_molecule_parameterizer(xyzfile="ligand.xyz",forcefield_option="GAFF",
        allow_undefined_stereo=True, resname="BEN")
    #Note: small_molecule_parameterizer creates a PDB-file: ligand.pdb (with conect lines)

    #############################################################
    #2. Orientation of new hydrogenated ligand (with a matching
    #FF XML file) into protein-ligand complex
    #############################################################
    #a. Read ligand PDB-file into ASH
    new_ligand_pdb="ligand.pdb"
    newligand = Fragment(pdbfile=new_ligand_pdb)
    print("New ligand coords:", newligand.print_coords())

    #b. Read ligand from a file containing ligand ATOM/HETATM lines from original PDB-structure (e.g. an X-ray structure with a bound-ligand)
    old_ligand_pdb="xray_2OXS_ligand.pdb" #This file should only contain the ligand
    oldligand = Fragment(pdbfile=old_ligand_pdb)
    print("Old ligand:", oldligand.print_coords())

    #c. Define the atoms in common in new and old ligand (at least carbon skeleton, all nonH-atoms should work)
    #Here defining a list of lists that contain the atom indices in new_ligand (system A) and old_ligand (systemB)
    subsetA=newligand.get_nonH_atomindices() #Getting atom indices of non-H atoms
    subsetB=oldligand.get_nonH_atomindices() #Getting atom indices of non-H atoms
    subset=[subsetA,subsetB] #Combining lists into a list-of-lists

    #d. Align new ligand (with H-atoms and matching XML-file) so that it matches (as well as possible) the position of the old-ligand atoms
    #Note: subset needs to be properly chosen. Reordering is usuaully necessary for alignment (because atom order may differ)
    newligand_aligned = flexible_align_pdb(new_ligand_pdb, old_ligand_pdb, subset=subset, reordering=True, reorder_method='brute')

    #############################################################
    #3. Merging protein and new aligned ligand
    #############################################################
    protein_pdbfile="xray_2OXS_protein.pdb"
    ligand_pdbfile="ligand_aligned.pdb"
    merged_pdbfile = merge_pdb_files(protein_pdbfile,ligand_pdbfile, outputname="2OXS_protein_ligand_merged.pdb")

    #############################################################
    #4. Finally  using OpenMM_Modeller to setup system
    #############################################################
    #The inputfiles required
    pdbfile="2OXS_protein_ligand_merged.pdb" #A merged protein-ligand complex PDB-file (needs to contain a ligand with all hydrogens)
    ligand_xmlfile="ligand.xml" #An XML-file containing the FF for the ligand

    #Calling OpenMM_Modeller
    openmmobject, ashfragment = OpenMM_Modeller(pdbfile=pdbfile, forcefield='Amber14', watermodel="tip3pfb",pH=7.0,
        solvent_padding=10.0, ionicstrength=0.1, extraxmlfile=ligand_xmlfile)



######################################################
**6. Run initial preparatory MD simulations**
######################################################

Before we can start running production MD simulations to explore protein-ligand binding scenarios or even free-energy simulations we must 
first run some initial preparatory MD simulations to equilibrate the system and remove any clashes between the protein and ligand and make sure the solvent is properly equilibrated.

The following script can be used to conveniently warm up the system (**Gentle_warm_up_MD** function) using a series of MD simulations 
with increasing temperature and time step before switching to **OpenMM_box_equilibration** which performs an NPT simulation until the 
density and volume of the system has converged.


.. code-block:: python

    from ash import *

    #Defining fragment containing coordinates
    pdbfile="finalsystem.pdb"
    fragment=Fragment(pdbfile=pdbfile)

    #Creating an OpenMMTheory object using XML-files and PDB-file (only used to define topology)
    omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "amber14/tip3pfb.xml", "ligand.xml"], 
                pdbfile=pdbfile, periodic=True,
                autoconstraints='HBonds', rigidwater=True)

    #Gentle warmup MD (3 MD simulations: 10/50/200 steps with timesteps 0.5/1/4 fs at 1 K/10K/300K)
    Gentle_warm_up_MD(fragment=fragment, theory=omm, time_steps=[0.0005,0.001,0.004], 
                steps=[10,50,200], temperatures=[1,10,300])
    
    #Run NPT simulation until density and volume converges
    OpenMM_box_equilibration(fragment=fragment, theory=omm, datafilename="nptsim.csv", numsteps_per_NPT=10000,
                      temperature=300, timestep=0.001, traj_frequency=100, trajfilename='equilbox_NPT', 
                      trajectory_file_option='DCD', coupling_frequency=1)

It is of course also possible to split this script up into 2 scripts. Just make sure to redfine the fragment object so that it reads a PDB-file that contains updated coordinates.


Inside the scripts directory of the main ASH source-code directory there is a script called **plot_md_data.py** 
that can be used to conveniently visualize the convergence of the density and volume data from the nptsim.csv file (created by **OpenMM_box_equilibration**)

.. code-block:: text

    #Plot density and volume from nptsim.csv via MatplotLib
    python3 plot_md_data.py nptsim.csv


######################################################
**7. Run long time-scale NVT simulation**
######################################################

Once the system has been properly equilibrated we can start running longer time-scale simulations to explore protein-ligand binding scenarios.
Here we will run a 1 ns NVT simulation using the LangevinMiddleIntegrator integrator.

.. note:: OpenMM MD simulations in general run much faster using a GPU than on the CPU. Use platform='CUDA' or platform='OpenCL' to run on the GPU.
    Using a modern graphics card, 1000 ns simulations should be achievable on a desktop in 1-3 days.

.. code-block:: python

    from ash import *

    #Defining fragment containing coordinates
    pdbfile="equilbox_NPT.pdb"
    fragment=Fragment(pdbfile=pdbfile)

    #Creating an OpenMMTheory object using XML-files and PDB-file (only used to define topology)
    omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "amber14/tip3pfb.xml", "ligand.xml"], 
                pdbfile=pdbfile, periodic=True,
                autoconstraints='HBonds', rigidwater=True)

    #Run a NVT MD simulation (NPT can also be performed if you add a barostat)
    OpenMM_MD(fragment=fragment, theory=omm, timestep=0.001, simulation_time=1000, traj_frequency=10, 
        temperature=30, platform='OpenCL', integrator='LangevinMiddleIntegrator', coupling_frequency=1, 
        trajfilename='NVT-MD',trajectory_file_option='DCD')

    #Re-image trajectory so that protein is in middle
    MDtraj_imagetraj("NVT-MD.dcd", "NVT-MD.pdb", format='DCD')


The resulting trajectory can be visualized using e.g. VMD. 
It is then best to use the "imaged" versions (requires **mdtraj**) of the trajectory file (NVT-MD_imaged.dcd) where the 
protein is in the middle of the box.

The usefulness of the unbiased MD trajectory depends on whether any kind of binding of the ligand to a protein pocket can be observed.
It is likely that a few hundred ns of unbiased MD simulations are required to even see any spontaneous binding event.


#########################################################
**8. Funnel metadynamics of the protein-ligand system**
#########################################################

In order to a realistically explore protein-ligand binding scenarios we need to use enhanced sampling methods.
Metadynamics is a general free-energy simulation method that is in principle well suited to study protein-ligand binding
as we could sample the free-energy surface of the bound vs. unbound conformation. Metadynamics use a history-dependent biasing potential
that is built-up using Gaussians during the simulation, preventing the simulation from visiting previous parts of the free-energy surface.
Metadynamics require the definition of one or more collective variables (CVs) that act as "reaction coordinates" for the biasing potential.

A metadynamics simulation for a binding reaction such as here, however, creates a problem as the ligand encounters 
the "unbound" part of the free energy surface (when the ligand is far away from the protein binding site).
The simulation can not realistically converge as the ligand will encounter a practically infinite amount of conformations 
outside the protein binding site.

To combat this problem we will use funnel metadynamics (https://www.pnas.org/doi/10.1073/pnas.1303186110) 
which adds a restraing potential with a funnel shape that prevents the ligand from escaping too far away from the protein binding site.

**THIS IS NOT YET COMPLETE**


#########################################################
**9. QM/MM  of the protein-ligand system**
#########################################################

**THIS IS NOT YET COMPLETE**