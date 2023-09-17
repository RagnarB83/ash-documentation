Modelling protein-ligand binding in ASH: Tutorial
====================================================

This is a tutorial on how to set up a model for studying organic ligand binding to a protein using ASH.
All the steps about how to set up the system, parameterize the ligand, solvate the system and run MD simulations
are performed within the ASH environment with the help of libraries such as OpenMM, OpenBabel, OpenFF, ParmEd, MDTraj, etc.
Later we will also show how to run metadynamics simulations to study the free-energy surface of the ligand binding.

Test system:

######################################################
**1. Preparing initial files**
######################################################

**Protein PDB-file**

To get started we need a PDB file of the protein. This can be an initial crystal structure from the PDB database or a previous model. Hydrogens do not need to be present in the PDB file, as they will be added during the setup if missing.

**Ligand PDB-file**

We also need a PDB file of the ligand. This PDB-file can be created from an XYZ-file that should have a decent-looking geometry (ideally DFT-optimized)
or from some other file format such as a MDL Mol-file, SDF file.
This PDB-file needs to contain all the hydrogen atoms and needs to have The file will be referered to as ligand.pdb here.

.. code-block:: python

    from ash import *

    #Example: Create PDB-file from an XYZ-file. Creates file: ligand.pdb
    pdbfile = ligand_pdbfile = xyz_to_pdb_with_connectivity("ligand.xyz")

    #Example: Create PDB-file with correct connectivity. Creates file: ligand__withcon.pdb
    pdbfile = writepdb_with_connectivity("ligand.pdb")

    #Other options:
    #Example: Create PDB-file from an MDL Mol-file (requires OpenBabel). Creates file: ligand.pdb
    pdbfile = mol_to_pdb("ligand.mol")
    #Example: Create PDB-file from an MDL SDF-file (requires OpenBabel). Creates file: ligand.pdb
    pdbfile = sdf_to_pdb("ligand.sdf")


**Merged protein-ligand PDB-file**

We also need a combined PDB-file that contains both the protein and the ligand combined. Such a PDB-file could be prepared by a GUI visualization program (VMD, Chimera). 
The ligand coordinates can not clash with the protein and can be either away from protein or in a prospective binding site. 
Docking could also be used to prepare such a system. It is important that the connectivity of the ligand is correct in the PDB-file by inclusion of CONECT statements.

If the coordinates of protein and ligand in their respective files are already compatible (not clashing) it is possible to use the 
ASH **merge_pdb_files** function to combine the two PDB-files into one. This function will preserve and update the CONECT statements of the ligands.

.. code-block:: python

    from ash import *
    protein_pdbfile="protein.pdb"
    ligand_pdbfile="ligand.pdb"
    merged_pdbfile = merge_pdb_files(protein_pdbfile,ligand_pdbfile, outputname="merged.pdb")

######################################################
**2. Preparing ligand forcefield**
######################################################

We then need to think about the forcefield. Various protein/nucleic-acid forcefields are available in ASH (CHARMM, Amber etc.) and can be used automatically.
However, the forcefield for the ligand is the main issue as it is unlikely present in the biomolecular forcefields.
We also need to consider the compatibility between the forcefield for the ligand and the forcefield for the protein.

Here we choose to use the Amber14 forcefield for the protein and the GAFF (Generalized Amber force field) forcefield for the ligand as this can be conveniently set up using ASH.
Another option is to use one of the OpenFF forcefields for the ligands (also compatible with Amber14)

ASH features a convenient function : **small_molecule_parameterizor** that can automatically generate the forcefield for the ligand
by determining the topology of the input ligand and matching it to general parameters available for either GAFF or OpenFF.
This is made possible by functionality available in the **openmm-forcefields** package (https://github.com/openmm/openmmforcefields) 
which needs to be installed when prompted.

The **small_molecule_parameterizor** function requires in principle only the input PDB-file for the ligand.
The PDB-file will be automatically converted into a SMILES string that is then used to generate the topology and suitable parameters
are found in the GAFF (or OpenFF) forcefield. If reason does not work it is also possible to provide other inputfiles to
small_molecule_parameterizor such as : xyzfile, molfile, sdffile. One can also provide a SMILES string (smiles_string keyword).


.. code-block:: python

    from ash import *
    #Get Amber compatible forcefield for small molecule using GAFF or OpenFF
     small_molecule_parameterizor(pdbfile="ligand.pdb", forcefield_option='GAFF', output_xmlfile="ligand.xml")


.. warning:: Make sure that the PDB-file, MOL-file or SDF-file atom ordering and connectivity matches the information present in the PDB-file that will be used to 
    to setup the system later. Otherwise the ligand will not be recognized.

The function returns an OpenMM forcefield object (that assumes Amber14 for protein and solvent and GAFF for the ligand)
and also writes out an XML-file with the forcefield parameters for the ligand (ligand.xml). 
It is usually best to use the ligand.xml file.


######################################################
**3. Prepare system using OpenMM_Modeller**
######################################################

Now we should have a merged PDB-file (containing both protein and ligand) and a forcefield for the ligand (ligand.xml).
We can now proceed to use the **OpenMM_Modeller** function to set up the system. We use the merged protein-ligand PDB-file to define the system geometry and topology, 
we specify an Amber14 forcefield for the protein, TIP3P forcefield for water (compatible with Amber14) and the ligand forcefield (GAFF or OpenFF) for the ligand via the 
ligand.xml file previously created.

See :doc:`OpenMM-interface` for detail on using OpenMM_Modeller.

.. code-block:: python

    from ash import *

    merged_pdbfile="merged.pdb"
    #Setup system using OpenMM_Modeller using merged PDB-file
    OpenMM_Modeller(pdbfile=merged_pdbfile, forcefield="Amber14",
        extraxmlfile="ligand.xml", residue_variants={}, watermodel="tip3p", pH=7.0, solvent_padding=10.0, ionicstrength=0.1)

OpenMM_Modeller will apply the Amber14 protein forcefield to the protein and the GAFF/OpenFF forcefield to the ligand.
Note that one must make sure that the merged PDB-file of the protein and ligand contains the correct connectivity information for the ligand (CONECT lines).
Additionally one must make sure that any residues in the protein are correctly treated (with respect to protonation states, disulfide bridges, metal ions etc.).

If the OpenMM_Modeller function is successful a final PDB-file, "finalsystem.pdb" will be created that contains the solvated protein-ligand system with
protein and ligand oriented according to the coordinates of "merged.pdb". The coordinates in the input "merged.pdb" file 
can contain the system in either bound or unbound form.

######################################################
**4. Run initial preparatory MD simulations**
######################################################

Before we can start running production MD simulations to explore protein-ligand binding scenarios or even free-energy simulations we must 
first run some initial preparatory MD simulations to equilibrate the system and remove any clashes between the protein and ligand and make sure the solvent is properly equilibrated.

The following script can be used to conveniently warm up the system (**Gentle_warm_up_MD** function) using a series of MD simulations 
with increasing temperature and time step before switching to OpenMM_box_equilibration which performs an NPT simulation until the 
density and volume of the system has converged.


.. code-block:: python

    from ash import *

    #Defining fragment containing coordinates
    pdbfile="finalsystem.pdb"
    fragment=Fragment(pdbfile=pdbfile)

    #Creating an OpenMMTheory object using XML-files and PDB-file (only used to define topology)
    omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "amber14/tip3pfb.xml", "gaff_ligand.xml"], 
                pdbfile=pdbfile, periodic=True,
                autoconstraints='HBonds', rigidwater=True)

    #Gentle warmup MD (3 MD simulations: 10/50/200 steps with timesteps 0.5/1/4 fs at 1 K/10K/300K)
    Gentle_warm_up_MD(fragment=fragment, theory=omm, time_steps=[0.0005,0.001,0.004], 
                steps=[10,50,200], temperatures=[1,10,300])
    
    #Run NPT simulation until density and volume converges
    OpenMM_box_equilibration(fragment=fragment, theory=omm, datafilename="nptsim.csv", numsteps_per_NPT=10000,
                      temperature=300, timestep=0.001, traj_frequency=100, trajfilename='equilbox_NPT', trajectory_file_option='DCD', coupling_frequency=1)

It is of course also possible to split this script up into 2 scripts. Just make sure to redfine the fragment object so that it reads a PDB-file that contains updated coordinates.


Inside the scripts directory of the main ASH source-code directory there is a script called **plot_md_data.py** 
that can be used to conveniently visualize the convergence of the density and volume data from the nptsim.csv file (created by **OpenMM_box_equilibration**)

.. code-block:: text

    #Plot density and volume from nptsim.csv via MatplotLib
    python3 plot_md_data.py nptsim.csv


######################################################
**5. Run long time-scale NVT simulation**
######################################################

One the system has been properly equilibrated we can start running longer time-scale simulations to explore protein-ligand binding scenarios.
Here we will run a 1 ns NVT simulation using the LangevinMiddleIntegrator integrator.

.. note:: OpenMM MD simulations in general run much faster using a GPU than on the CPU. Use platform='CUDA' or platform='OpenCL' to run on the GPU.

.. code-block:: python

    from ash import *

    #Defining fragment containing coordinates
    pdbfile="equilbox_NPT.pdb"
    fragment=Fragment(pdbfile=pdbfile)

    #Creating an OpenMMTheory object using XML-files and PDB-file (only used to define topology)
    omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "amber14/tip3pfb.xml", "gaff_ligand.xml"], 
                pdbfile=pdbfile, periodic=True,
                autoconstraints='HBonds', rigidwater=True)

    #Run a NVT MD simulation (NPT can also be performed if you add a barostat)
    OpenMM_MD(fragment=fragment, theory=omm, timestep=0.001, simulation_time=1000, traj_frequency=10, temperature=30, platform='OpenCL',
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajfilename='NVT-MD',trajectory_file_option='DCD')

    #Re-image trajectory so that protein is in middle
    MDtraj_imagetraj("NVT-MD.dcd", "NVT-MD.pdb", format='DCD')


The resulting trajectory can be visualized using e.g. VMD. 
It is then best to use the "imaged" versions (requires **mdtraj**) of the trajectory file (NVT-MD_imaged.dcd) where the 
protein is in the middle of the box.

The usefulness of the unbiased MD trajectory depends on whether any kind of binding of the ligand to a protein pocket can be observed.



#########################################################
**6. Funnel metadynamics of the protein-ligand system**
#########################################################

In order to a realistically explore protein-ligand binding scenarios we need to use enhanced sampling methods.
Metadynamics is a general free-energy simulation method that is in principle well suited to study protein-ligand binding
as we could sample the free-energy surface of the bound vs. unbound conformation.

The trouble is that when a metadynamics simulation encounters the "unbound" part of the free energy surface
(when the ligand is far away from the protein binding site) the simulation can not realistically converge as the ligand
will encounter a practically infinite amount of conformations outside the protein binding site.

To combat this problem we turn to funnel metadynamics (https://www.pnas.org/doi/10.1073/pnas.1303186110) 
which adds a restraing potential with a funnel shape that prevents the ligand from escaping too far away from the protein binding site.

**THIS IS NOT YET COMPLETE**


#########################################################
**7. QM/MM  of the protein-ligand system**
#########################################################

**THIS IS NOT YET COMPLETE**