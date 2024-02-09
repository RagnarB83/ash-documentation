Modelling protein-ligand binding in ASH: Tutorial
====================================================

This is a tutorial on how to set up a model for studying organic ligand binding to a protein using ASH, 
both via a classical approach, a QM/MM approach and also via optimizations, MD simulations and free-energy simulations.

All the steps about how to set up the system, parameterize the ligand, solvate the system and run MD simulations
are performed within the ASH environment with the help of libraries such as OpenMM, OpenFF, RDKit, OpenBabel, ParmEd, MDTraj, etc.
Later we will also show how to run metadynamics simulations to study the free-energy surface of the ligand binding.

The tutorial below will be discussed in a general way, with, however, a specific example shown: 
the trypsin-benzamidine protein-ligand complex with PDB-ID:  `2OXS <https://www.rcsb.org/structure/2OXS>`_.

Note that in order to use the small_molecule_parameterizer function in ASH you need to install the `openmm-forcefields package <https://github.com/openmm/openmmforcefields>`_ like this:
.. code-block:: python

    mamba install openmmforcefields
    #Or: conda install --yes -c conda-forge openmmforcefields

This will additionally install the `OpenFF Toolkit <http://github.com/openforcefield/openff-toolkit>`_ and `RDKit <https://github.com/rdkit/rdkit>`_ which are required for creating the ligand forcefield.

######################################################
**1. Preparing initial files**
######################################################

**Protein PDB-file**

To get started we need a PDB file of the protein. This can be an initial crystal structure from the PDB database or a previous model. Hydrogens do not need to be present in the PDB file, as they will be added during the setup if missing.
You can use a PDB-file of a protein-ligand complex but you will have to prepare a version of the PDB-file that only contains the protein (no ligand).
We will use here the trypsin-benzamidine complex with PDB-ID: 2OXS. In addition to the protein, it contains benzamidine (BEN residue), a sulfate ion (SO4 residue) and a calcium ion (CA).
We will not model the sulfate ion or the calcium ion (crystallization species) and will delete them. Additionally we split the file so that we have files for the protein separately and the ligand separately.
This is best to do manually (note that you only need to preserve ATOM/HETATM lines, everything else, including the long header can usually be deleted)

We will here create files with the following naming convention:

- 2oxs_full.pdb #This is the original PDB-file
- 2oxs-protein.pdb #This file contains only ATOM lines for the protein. Ligand and ion (SO4) has been manually removed
- 2oxs-ligand.pdb #This file contains only ATOM (or HETATM) lines for the ligand at it's original coordinates in X-ray structure.

The 2oxs-ligand.pdb should look like below:

.. toggle::

    .. code-block:: text

        HETATM 1649  C1  BEN A 801      -1.824  14.356 -17.123  1.00 11.46           C
        HETATM 1650  C2  BEN A 801      -2.035  15.726 -17.073  1.00 14.31           C
        HETATM 1651  C3  BEN A 801      -1.737  16.436 -15.918  1.00 15.80           C
        HETATM 1652  C4  BEN A 801      -1.246  15.789 -14.795  1.00 13.14           C
        HETATM 1653  C5  BEN A 801      -1.020  14.434 -14.843  1.00 14.33           C
        HETATM 1654  C6  BEN A 801      -1.308  13.728 -15.998  1.00 13.03           C
        HETATM 1655  C   BEN A 801      -2.138  13.578 -18.351  1.00 12.13           C
        HETATM 1656  N1  BEN A 801      -2.817  14.063 -19.305  1.00 12.82           N
        HETATM 1657  N2  BEN A 801      -1.742  12.315 -18.449  1.00 11.04           N

**Ligand coordinate-file**

The ligand in an X-ray structure is unlikely to contain H-atoms but we need them in order to create a reliable forcefield that can be used in MD simulations.
We need to create a separate coordinate file for the ligand that contains all the hydrogen atoms and has a sensible internal geometry.
This can be done in a visualization program such as Chemcraft or Avogadro.
Here we recommend simply creating XYZ-file that has an acceptable initial geometry (ideally DFT-optimized from e.g. ORCA) to be used for the next step.

Since the benzamidine ligand likely exists as a protonated bound cation in the protein-ligand complex we will model it as such: 2 H-atoms on each nitrogen and a total charge of +1.

.. code-block:: python
    
    from ash import *
    frag = Fragment(xyzfile="BEN_cation_initial.xyz")
    #Simple r2SCAN-3c CPCM theory
    theory = ORCATheory(orcasimpleinput="! r2SCAN-3c CPCM tightscf")
    Optimizer(theory=theory, fragment=frag)

where BEN_cation_initial.xyz is:

.. toggle::

    .. code-block:: text

        18
        initial geometry drawn in Chemcraft
        C      -1.762630      1.524696      1.900268
        N      -1.257078      2.198475      2.972357
        N      -3.077175      1.655179      1.654487
        H      -1.882798      2.685338      3.601707
        H      -3.489718      0.993381      1.009265
        H      -3.521783      2.559541      1.569628
        C      -0.905586      0.690315      1.077732
        C       0.307642      0.191579      1.588217
        C      -1.274647      0.390052     -0.243916
        C       1.138590     -0.567317      0.782632
        C      -0.440260     -0.376761     -1.040084
        C       0.765119     -0.853149     -0.529881
        H       0.597455      0.343481      2.624281
        H      -2.184746      0.790936     -0.680460
        H       2.070072     -0.958332      1.177935
        H      -0.721899     -0.594827     -2.064652
        H       1.417349     -1.455059     -1.155223
        H      -0.266968      2.319809      3.134665

and after optimization:

.. toggle::

    .. code-block:: text

        18
        Coordinates from ORCA-job BEN-cation-opt
        C   -1.76637201220656      1.57665753553453      1.86821779542988
        N   -1.23141932878006      2.36306049909141      2.78425220431610
        N   -3.07242649498273      1.53850181376330      1.67620767483102
        H   -1.80136345579512      2.91828994571365      3.40837328515338
        H   -3.47747372898837      0.84294028894981      1.06844067153159
        H   -3.70292643677967      2.14298185221503      2.18540026161394
        C   -0.88974908575088      0.72549192815858      1.04457039055598
        C   0.23985943059254      0.12816154736947      1.61600703579234
        C   -1.18899041730278      0.51501929822992     -0.30667188687504
        C   1.05984922036702     -0.67813361570197      0.83700910420152
        C   -0.35798964752115     -0.28690904409586     -1.07832625926310
        C   0.76398685496581     -0.88514602803884     -0.50863885456430
        H   0.46015869768885      0.26462924366790      2.67082627889920
        H   -2.04721387592824      1.00011436526987     -0.76256088903936
        H   1.92829015953726     -1.15180624038479      1.28439747883171
        H   -0.58399911032055     -0.43921493182412     -2.12926136220035
        H   1.40931983776603     -1.51423945878798     -1.11481276840237
        H   -0.23060160656140      2.47693800087007      2.83552783918784


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

The **small_molecule_parameterizer** function requires in principle only an XYZ-file of the ligand (containing all H-atoms), the desired total charge and the forcefield option (GAFF or OpenFF)
The XYZ-coordinates will be fed to RDKit (installed ) which will generate the correct connectivity and bond orders that can then be passed onto the OpenFF toolkit.
OpenFF toolkit can next generate a forcefield for the ligand (either GAFF or OpenFF).
**small_molecule_parameterizer** can read as input : xyzfile, molfile, sdffile. One can also provide a SMILES string (smiles_string keyword).
Generally we recommend an XYZ-file.


.. code-block:: python

    from ash import *
    #Create an Amber-compatible forcefield for a small molecule using GAFF or OpenFF
     small_molecule_parameterizer(xyzfile="BEN-cation-opt.xyz", forcefield_option='GAFF', resname="BEN", charge=1)
    #This will create a BEN.pdb file and a gaff_BEN.xml file

The function writes out an XML-file with the forcefield parameters for the ligand (here BEN.xml) and also writes out a compatible PDB-file (here BEN.pdb).
Do note that the atom ordering may have changed compared to the input XYZ-file. This PDB-file will contain CONECT lines for the ligand (necessary for OpenMM to recognize the ligand).

.. note:: If you don't wish to use small_molecule_parameterizer (or if it fails; contact us if that is the case) you could prepare an OpenMM XML-file for the ligand in some other way. Make sure that the PDB-file atom ordering and names match the XML-file.
    

######################################################
**3. Merge and align protein and ligand**
######################################################

We now have a PDB-file for the ligand (BEN.pdb) that contains the ligand with H-atoms and correct connectivity and an OpenMM XML file (gaff_BEN.xml, containing the forcefield parameters for the ligand).
We could in principle proceed to set up the system. However, first we need to merge the protein and ligand into one PDB-file (as OpenMM_Modeller expects a single PDB-file) and we need to make sure that the ligand is properly aligned in the protein.
If we don't care about the ligand being in a specific position w.r.t. the protein, we could simply visualize 2oxs_protein.pdb and the ligand PDB-file in e.g. VMD, to make sure that protein and ligand do not clash and are reasonably close.
Otherwise modify the coordinates of the ligand in the ligand PDB file. This would be fine if want to initially study the unbound form of the system or possible predict binding by MD later.

However, if we want to start a simulation with the ligand in the original binding site according to the X-ray structure then we have to make sure that the new hydrogenated-ligand we created is properly aligned in the protein.
This would require either modifying the coordinates of the ligand in the ligand.pdb file using a suitable visualization program (e.g. VMD), perform docking,  or alternatively we could superimpose the new hydrogenated ligand onto the original ligand-position in the X-ray structure.
Here we will show how to do the latter using ASH using the **flexible_align_pdb** function in ASH.

**Align the ligand onto the desired previous position**

.. code-block:: python

    #a. Read hydrogenated ligand PDB-file into ASH
    new_ligand_pdb="BEN.pdb"
    newligand = Fragment(pdbfile=new_ligand_pdb)
    print("New ligand coords:", newligand.print_coords())

    #b. Read ligand from a file containing only the ligand ATOM/HETATM lines from original PDB-structure (e.g. an X-ray structure with a bound-ligand)
    old_ligand_pdb="2oxs_ligand.pdb" #This file should only contain the ligand. Probably missing H-atoms.
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

The **flexible_align_pdb** function creates a new PDB-file called BEN_aligned.pdb that contains the ligand in the same position as the old ligand. Unlike before, the new ligand contains all H-atoms and has a corresponding forcefield XML-file (same atomordering).


**Merged protein-ligand PDB-file**

Now that we have the ligand PDB-file, oriented and aligned the way we want, we can merge protein and ligand back together into a single PDB-file.
We can use the **merge_pdb_files** function in ASH to do this. This function is convenient as it will preserve and update the CONECT statements of the ligands which is important for the **OpenMM_Modeller** step later.

.. code-block:: python

    from ash import *
    protein_pdbfile="2oxs_protein.pdb"
    ligand_pdbfile="BEN_aligned.pdb" #This is the aligned ligand PDB-file (i.e. having the geometry we want). Atom-order needs to match information in ligand.xml
    merged_pdbfile = merge_pdb_files(protein_pdbfile,ligand_pdbfile, outputname="merged.pdb")



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
        extraxmlfile="gaff_BEN.xml", residue_variants={}, watermodel="tip3p-fb", pH=7.0, solvent_padding=10.0, ionicstrength=0.1)

**OpenMM_Modeller** will apply the Amber14 protein forcefield to the protein and the GAFF/OpenFF forcefield to the ligand.
Note that one must make sure that the merged PDB-file of the protein and ligand contains the correct connectivity information for the ligand (CONECT lines).
Additionally one must make sure that any residues in the protein are correctly treated (with respect to protonation states, disulfide bridges, metal ions etc.). 
We will not go into this aspect in this tutorial but we emphasize that this is a vital step in setting up any biomolecular system correctly.

If the **OpenMM_Modeller** function is successful, a final PDB-file, "finalsystem.pdb" will be created that contains the solvated protein-ligand system with
protein and ligand oriented according to the initial coordinates of "merged.pdb". The coordinates in the input "merged.pdb" file 
can contain the system in either bound or unbound form and can be modified before running **OpenMM_Modeller**. 
Note that due to the present of the solvent, it is trickier to change the ligand position of the solvated system after the **OpenMM_Modeller** step
(would require running a biased MD simulation).

.. warning:: Make sure that the ligand geometry in the merged PDB-file matches the information in the ligand.xml file. Otherwise the ligand will not be recognized by OpenMM.


######################################################
**5. STEPS 1-4 COMBINED**
######################################################

Here we show a script that combines the steps 1-4 into a single ASH script that could in principle be used to conveniently perform all the steps in one go.

.. code-block:: python

    from ash import *

    original_protein_pdbfile="2oxs-protein.pdb" #This file should only contain the protein
    original_ligand_pdbfile="2oxs-ligand.pdb" #This file should only contain the ligand
    #############################################################
    #1. Parameterize ligand using a hydrogenated XYZ-structure
    #############################################################
    residue_name="BEN" #A 3-letter name for ligand-residue (used to name files as well)
    #Here choosing GAFF
    small_molecule_parameterizer(xyzfile="BEN-cation-opt.xyz",forcefield_option="GAFF", resname=residue_name, charge=1)
    #Note: small_molecule_parameterizer creates a PDB-file: BEN.pdb (with conect lines)

    #############################################################
    #2. Orientation of new hydrogenated ligand (with a matching
    #FF XML file) into protein-ligand complex
    #############################################################
    #a. Read ligand PDB-file into ASH
    new_ligand_pdb=f"{residue_name}.pdb"
    newligand = Fragment(pdbfile=new_ligand_pdb)
    print("New ligand coords:", newligand.print_coords())

    #b. Read ligand from a file containing ligand ATOM/HETATM lines from original PDB-structure (e.g. an X-ray structure with a bound-ligand)
    old_ligand_pdb=original_ligand_pdbfile #This file should only contain the ligand
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
    protein_pdbfile=original_protein_pdbfile
    ligand_pdbfile=f"{residue_name}_aligned.pdb"
    merged_pdbfile = merge_pdb_files(protein_pdbfile,ligand_pdbfile, outputname="merged.pdb")

    #############################################################
    #4. Finally  using OpenMM_Modeller to setup system
    #############################################################
    #The inputfiles required
    pdbfile="merged.pdb" #A merged protein-ligand complex PDB-file (needs to contain a ligand with all hydrogens)
    ligand_xmlfile=f"gaff_{residue_name}.xml" #An XML-file containing the FF for the ligand

    #Calling OpenMM_Modeller
    openmmobject, ashfragment = OpenMM_Modeller(pdbfile=pdbfile, forcefield='Amber14', watermodel="TIP3P",pH=7.0,
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
    omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "amber14/tip3pfb.xml", "gaff_ligand.xml"], 
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
    omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "amber14/tip3pfb.xml", "gaff_ligand.xml"], 
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