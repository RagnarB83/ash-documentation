Explicit solvation (small molecule)
======================================

ASH allows you to easily create explicit solvation models for small molecules that can then be either
used either for classical MD simulations or QM/MM MD simulations.


################################################################################################
Creating a small-molecule forcefield
################################################################################################

An important step is to acquire a forcefield for the small molecule of interest.
Ideally we want a forcefield in the OpenMM XML format that can be easily used to define an OpenMMTheory object.

There are 2 approaches possible, a nonbonded-only forcefield and full-fledged forcefield.

**write_nonbonded_FF_for_ligand**

This function can create a simple nonbonded forcefield. This is useful for inorganic molecules (e.g. metal complexes) where a full
forcefield is hard to create. Typically the easiest option if the final goal is to perform QM/MM MD simulations. 
Note that for a classical MD simulation with a nonbonded forcefield, the internal degrees of freedom of the molecule should be constrained (either by freezing atoms or constraining internal coordinates).

**small_molecule_parameterizer**

This function is capable of creating/deriving a full forcefield of an organic molecule with both bonded and nonbonded parameters. 
Only available for molecules with common light maingroup elements. Options are the GAFF or OpenFF forcefields.

See :doc:`OpenMM-interface` for more details on how to use these functions.

################################################################################################
Example 1. Modelling of an organic molecule in explicit water with a ligand forcefield
################################################################################################

.. note:: Do note that when running OpenMM simulations you typically want to run simulations on the GPU if possible (much faster than the default CPU).
    You then have to set the platform keyword equal to either platform='CUDA' or platform='OpenCL' when defining the OpenMM_Theory object (default is 'CPU').


1. Create a forcefield file

*create_ff.py:*

.. code-block:: python
        
    from ash import *

    #Parameterize small molecule using OpenFF
    small_molecule_parameterizer(xyzfile="isobutyraldehyde.xyz", forcefield_option="OpenFF", charge=0)

The function will create a file, called openff_LIG.xml.

2. Solvate the system

*solvate.py:*

.. code-block:: python
        
    from ash import *

    mol = Fragment(xyzfile="isobutyraldehyde.xyz", charge=0, mult=1)
    #Solvate molecule in a 30x30x30 Å TIP3P water box.
    solvate_small_molecule(fragment=mol, xmlfile="openff_LIG.xml", watermodel='tip3p', solvent_boxdims=[30,30,30])

This creates files: system_aftersolvent.xyz and system_aftersolvent.pdb

3. Perform initial classical MD simulation to warmup and equilibrate 

*run_mm_md.py:*

.. code-block:: python

    from ash import *

    #Read in coordinates of the full system
    pdbfile="system_aftersolvent.pdb"
    fragment = Fragment(pdbfile=pdbfile)
    #Create an OpenMMTheory object based on PDB-file and XML-files for water and small-molecule
    omm =OpenMMTheory(xmlfiles=["openff_LIG.xml", "amber/tip3p_standard.xml"],
                pdbfile=pdbfile, periodic=True, rigidwater=True, autoconstraints='HBonds')

    #Gently warms up the system
    Gentle_warm_up_MD(fragment=fragment, theory=omm)

    #Equilibrates the system via a multi-step NPT simulation. 
    #This changes the box size of the system until volume and density have converged.
    #Note that thresholds for volume and density may have to be adjusted
    OpenMM_box_equilibration(fragment=fragment, theory=omm, datafilename="nptsim.csv", timestep=0.004,
                                numsteps_per_NPT=10000,max_NPT_cycles=10,traj_frequency=100,
                                volume_threshold=1.0, density_threshold=0.01, temperature=300)

You can use the script **plot_md_data.py** (present in the scripts directory of the main ASH directory) to conveniently plot the data (requires Matplotlib).

.. code-block:: python

    #Plot density, volume, temperature from nptsim.csv via MatplotLib
    python3 plot_md_data.py nptsim.csv

Re-imaging of the trajectory is often desirable for visualization purposes. This centers the system on the molecule.
This is also required for QM/MM MD simulation (where the molecule must be in the middle)

.. code-block:: python

    #Re-image trajectory so that protein is in middle
    MDtraj_imagetraj("equilibration_NPT.dcd", "equilibration_NPT.pdb", format='DCD')
    
    #Sometimes the procedure fails unless you specify that solute_anchor=True
    MDtraj_imagetraj("warmup_MD_cycle2.pdb","warmup_MD_cycle2.pdb", solute_anchor=True)


4. Run QM/MM MD simulation

*run_qmmm_md.py:*

.. code-block:: python

    from ash import *

    #Read in coordinates of the full system
    #Note that for QM/MM you must use a box where the molecule is centered. Re-image the file using MDtraj_imagetraj if necessary
    pdbfile="equilibration_NPT_imaged.pdb"
    fragment = Fragment(pdbfile=pdbfile)
    #Create an OpenMMTheory object based on PDB-file and XML-files for water and small-molecule
    omm =OpenMMTheory(xmlfiles=["openff_LIG.xml", "amber/tip3p_standard.xml"],
                pdbfile=pdbfile, periodic=True, rigidwater=True, autoconstraints='HBonds')
    #Create a QM/MM object
    qm = xTBTheory(xtbmethod='GFN2')
    #Defining QM-atoms to be the solute.  Note that the atom indices are 0-based
    qmatomlist = list(range(0,13))
    #QM/MM from QM and MM objects. Setting QM-region charge and multiplicity
    qm_mm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=fragment, qmatoms=qmatomlist,
            qm_charge=0, qm_mult=1)

    #Run a NVT MD simulation (NPT could also be performed if you add a barostat)
    #Note: timesteps for QM/MM must be much smaller than in MM
    OpenMM_MD(fragment=fragment, theory=qm_mm, timestep=0.001, simulation_time=10, traj_frequency=10, 
        temperature=300, integrator='LangevinMiddleIntegrator', coupling_frequency=1, 
        trajfilename='QM_MM_NVT-MD',trajectory_file_option='DCD')
    


#########################################################################################################
Example 2. Modelling of an inorganic molecule in explicit water using a simple non-bonded forcefield
#########################################################################################################

For inorganic molecules (e.g. metal complexes) it is trickier to parameterize a full forcefield and would 
have to be performed to some extent manually. However, it is possible to use a simple non-bonded forcefield,
i.e. without any bonded parameters.  This then requires classical simulations to be performed with internal degrees of freedom frozen 
(bonds, angles, dihedrals) while QM/MM MD simulations can be performed as normal.

Here we use the **write_nonbonded_FF_for_ligand** function to define the nonbonded parameters (charges and LJ parameters) and 
create an OpenMM XML-file. The molecule will be an FeCl4- complex (S=5/2).


1. Create a forcefield file

*create_ff.py:*

.. code-block:: python
        
    from ash import *

    #Create a nonbonded FF for molecule
    frag = Fragment(xyzfile="fecl4.xyz", charge=-1, mult=6)
    #Defining QM-theory to be used for charge calculation
    orca_theory = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf")
    #
    write_nonbonded_FF_for_ligand(fragment=frag, resname="LIG", theory=orca_theory,
            coulomb14scale=1.0, lj14scale=1.0, charge_model="CM5_ORCA")

The function will create a file, here called: "LIG.xml". By default it uses the ff_type to be "AMBER". This means the XML-file will
use a form of the nonbonded potential that is compatible with Amber-style forcefield. This is recommended im general but can be changed to "CHARMM" or "None" if required.


2. Solvate the system

We can solvate the system as we did before

*solvate.py:*

.. code-block:: python
        
    from ash import *

    mol = Fragment(xyzfile="fecl4.xyz", charge=-1, mult=6)
    #Solvate molecule in a 30x30x30 Å TIP3P water box.
    solvate_small_molecule(fragment=mol, xmlfile="LIG.xml", watermodel='tip3p', solvent_boxdims=[30,30,30])

This creates files: system_aftersolvent.xyz and system_aftersolvent.pdb

3. Perform initial classical MD simulation to warmup and equilibrate.

Here we run some initial classical MD. Unlike before, however, we have to constrain the internal degrees of freedom of the ligand
as there are no bonded parameters available. The simplest way is to add constraints for all the Fe-Cl bonds.
Additional angle constraints or dihedral constraints may also be required for other molecules (or )

*run_mm_md.py:*

.. code-block:: python

    from ash import *

    #Read in coordinates of the full system
    pdbfile="system_aftersolvent.pdb"
    fragment = Fragment(pdbfile=pdbfile)

    #Constrain the bonds of the ligand.
    #Note: additional angle and dihedral constraints may also be appropriate
    bondconstraints = [[0,1],[0,2],[0,3],[0,4]]

    #Create an OpenMMTheory object based on PDB-file and XML-files for water and small-molecule
    omm =OpenMMTheory(xmlfiles=["LIG.xml", "charmm36/water.xml"],
                pdbfile=pdbfile, periodic=True, rigidwater=True, autoconstraints='HBonds',
                constraints=bondconstraints)

    #Gently warms up the system
    Gentle_warm_up_MD(fragment=fragment, theory=omm, use_mdtraj=False)

    #Equilibrates the system via a multi-step NPT simulation.
    #This changes the box size of the system until volume and density have converged.
    #Note that thresholds for volume and density may have to be adjusted
    OpenMM_box_equilibration(fragment=fragment, theory=omm, datafilename="nptsim.csv", timestep=0.004,
                                numsteps_per_NPT=10000,max_NPT_cycles=10,traj_frequency=100,
                                volume_threshold=1.0, density_threshold=0.01, temperature=300)


You can use the script **plot_md_data.py** (present in the scripts directory of the main ASH directory) to conveniently plot the data (requires Matplotlib).

.. code-block:: python

    #Plot density, volume, temperature from nptsim.csv via MatplotLib
    python3 plot_md_data.py nptsim.csv

Re-imaging of the trajectory is often desirable for visualization purposes. This centers the system on the molecule.
This is also required for QM/MM MD simulation (where the molecule must be in the middle)

.. code-block:: python

    #Re-image trajectory so that protein is in middle
    MDtraj_imagetraj("equilibration_NPT.dcd", "equilibration_NPT.pdb", format='DCD')
    #Sometimes the procedure fails unless you specify that solute_anchor=True
    MDtraj_imagetraj("warmup_MD_cycle2.pdb","warmup_MD_cycle2.pdb", solute_anchor=True)

This is required for the QM/MM MD.

4. Run QM/MM MD simulation

*run_qmmm_md.py:*

.. code-block:: python

    from ash import *

    #Read in coordinates of the full system
    #Note that for QM/MM you must use a box where the molecule is centered. Re-image the file using MDtraj_imagetraj if necessary
    pdbfile="equilibration_NPT_imaged.pdb"
    fragment = Fragment(pdbfile=pdbfile)
    #No constraints necessary anymore as the solute will be in the QM-region
    
    #Create an OpenMMTheory object based on PDB-file and XML-files for water and small-molecule
    omm =OpenMMTheory(xmlfiles=["LIG.xml", "charmm36/water.xml"],
                pdbfile=pdbfile, periodic=True, rigidwater=True, autoconstraints='HBonds')
    #Create a QM/MM object
    qm = xTBTheory(xtbmethod='GFN2')
    #Defining QM-atoms to be the solute.  Note that the atom indices are 0-based
    qmatomlist = list(range(0,5))
    #QM/MM from QM and MM objects. Setting QM-region charge and multiplicity
    qm_mm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=fragment, qmatoms=qmatomlist,
            qm_charge=0, qm_mult=1)

    #Run a NVT MD simulation (NPT could also be performed if you add a barostat)
    #Note: timesteps for QM/MM must be much smaller than in MM
    OpenMM_MD(fragment=fragment, theory=qm_mm, timestep=0.001, simulation_time=10, traj_frequency=10, 
        temperature=300, integrator='LangevinMiddleIntegrator', coupling_frequency=1, 
        trajfilename='QM_MM_NVT-MD',trajectory_file_option='DCD')
    


#########################################################################################################
Example 3. Setting up an explicit non-aqueous solution system (with an inorganic solute)
#########################################################################################################

It is also possible to easily setup a non-aqueous solution system thanks to an interface to Packmol program.
Here we will use acetonitrile as an organic solvent and a metal complex as a solute. 
We will utilize a full forcefield for the solvent but a nonbonded forcefield for the solute.
This type of workflow is e.g. suitable for modelling transition metal catalysis in explicit solution.

1. Acquire the forcefield for the solvent (acetonitrile). 

This can be done as in Example 1.

*create_solvent_ff.py:*

.. code-block:: python

    #Parameterize small molecule using OpenFF
    small_molecule_parameterizer(forcefield_option="OpenFF", xyzfile="acetonitrile.xyz", 
                charge=0,resname="ACN")

The function will create an XML file, called openff_ACN.xml and a PDB-file called ACN.pdb

2. Creating a solvation box with a non-aqueous solvent via Packmol interface

See :doc:`helper_programs` for details on the Packmol interface.

Starting from a PDB-file containing a single acetonitrile molecule, we can create a cubic box of acetonitrile molecules with a specified density.
Here we specify the minimum and maximum coordinates of the box to be [0,0,0] and [50,50,50] Å.

.. code-block:: python

    #Create a 50 Å cubic box of acetonitrile molecules corresponding to a density of 0.786 g/ml
    packmol_solvate(inputfiles=["ACN.pdb"], density=0.786,
        min_coordinates=[0,0,0], max_coordinates=[50,50,50])

    #NPT equilibration. Will give optimal box dimensions
    pdbfile="final_withcon.pdb"
    fragment = Fragment(pdbfile=pdbfile)
    #Note: using slightly larger box dimensions (55 instead of 50) to avoid initial periodicity problems at boundary
    omm = OpenMMTheory(xmlfiles=["openff_ACN.xml"],pdbfile=pdbfile, platform="OpenCL",
                periodic=True, autoconstraints='HBonds', periodic_cell_dimensions=[55.0,55.0,55.0,90.0,90.0,90.0])
    #NPT equilibration. Note: platform='OpenCL' (or CUDA if NVIDIA GPU) runs OpenMM on GPU, should run quite fast even on laptop. 
    #Use platform='CPU' if no GPU available
    OpenMM_box_equilibration(fragment=fragment, theory=omm, datafilename="nptsim.csv", timestep=0.001, 
                                    numsteps_per_NPT=10000,max_NPT_cycles=10,traj_frequency=100,
                                    volume_threshold=1.0, density_threshold=0.01, temperature=300)

**packmol_solvate** will create a PDB-file called "final_withcon.pdb" containing coordinates for a 50 Å cubic box of acetonitrile molecules.
Note that connectivity lines are written at the end of the file which are necessary for OpenMM to recognize the topology.
Here we also immediately launch a classical NPT-equilibration simulation both to check that the solvent XML-file works and to get optimal box dimensions.
**OpenMM_box_equilibration** will create file equilibration_NPT_lastframe.pdb containing last frame of NPT-workflow and will also contain the equilibrated PBC box dimensions.

.. note:: It is important to read in the same single-molecule solvent PDB-file (here ACN.pdb) that you got from the **small_molecule_parameterizer** function,
    otherwise OpenMMTheory will probably not be able to recognize the topology.

3. Create a nonbonded forcefield for the solute

This can be done as in Example 2 above.

*create_solute_ff.py:*

.. code-block:: python
        
    from ash import *

    #Create a nonbonded FF for molecule
    frag = Fragment(xyzfile="fecl4.xyz", charge=-1, mult=6)
    #Defining QM-theory to be used for charge calculation
    orca_theory = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf")
    #
    write_nonbonded_FF_for_ligand(fragment=frag, resname="TMC", theory=orca_theory,
            coulomb14scale=1.0, lj14scale=1.0, charge_model="CM5_ORCA")
    #Write PDB-file (TMC.pdb) for later use, skipping connectivity lines
    frag.write_pdbfile_openmm("TMC", skip_connectivity=True)

This will create an XML-file of the solute called TMC.xml and a basic PDB-file called TMC.pdb.

.. note:: The PDB-file created here for the metal complex should not contain any connectivity lines (CONECT records), which is why we select the \
    skip_connectivity=True option. This is because we have a nonbonded FF and no bonded parameters in the XML-file and presence of connectivity lines would cause OpenMM to expect bond terms present. 

4. Inserting the solute into the box of solvent

See documentation of *insert_solute_into_solvent* at :doc:`OpenMM-interface` for more details.

Finally we want to combine the solute and solvent into a single PDB-file of the solution with the solute in the center of the box.
For this we need PDB-files for the single solute-molecule and solvent-box.

.. code-block:: python

    from ash import *

    #Input files
    solute_pdbfile="TMC.pdb" #should have no CONECT lines
    solvent_pdbfile="equilibration_NPT_lastframe.pdb" # using NPT-equilibrated solvent-box
    solute_xmlfile="TMC.xml"
    solvent_xmlfile="openff_ACN.xml"

    #Inserting solute into solvent-box and get new solution fragment and file solution.pdb
    solution = insert_solute_into_solvent(solvent_pdb=solvent_pdbfile, solute_pdb=solute_pdbfile,
                write_pdb=True)

    #Test to see if OpenMMTheory object can be defined from XML-files and final system file (solution.pdb)
    omm = OpenMMTheory(xmlfiles=[solute_xmlfile,solvent_xmlfile],pdbfile="solution.pdb", periodic=True)

**insert_solute_into_solvent** will create a PDB-file called "solution.pdb" containing the set-up system.
Note that the solution.pdb file will serve as both a topology for OpenMMTheory and as starting coordinates for an MD simulation.

We can now in principle run a QM/MM MD simulation of the system.

*run_qmmm_md.py:*

.. code-block:: python

    from ash import *

    #Read in coordinates of the full system
    pdbfile="solution.pdb"
    fragment = Fragment(pdbfile=pdbfile)

    solute_xmlfile="TMC.xml"
    solvent_xmlfile="openff_ACN.xml"

    omm = OpenMMTheory(xmlfiles=[solute_xmlfile,solvent_xmlfile],pdbfile=pdbfile, platform='OpenCL',
                periodic=True, autoconstraints='HBonds')
    #Create a QM/MM object
    qm = xTBTheory(xtbmethod='GFN2')
    #Defining QM-atoms to be the solute.  Note that the atom indices are 0-based
    qmatomlist = list(range(0,5))
    #QM/MM from QM and MM objects. Setting QM-region charge and multiplicity
    qm_mm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=fragment, qmatoms=qmatomlist,
            qm_charge=-1, qm_mult=6)

    #Run a NVT QM/MM MD simulation
    OpenMM_MD(fragment=fragment, theory=qm_mm, timestep=0.001, simulation_time=10, traj_frequency=1,
        temperature=30, integrator='LangevinMiddleIntegrator', coupling_frequency=1,
        trajfilename='QM_MM_NVT-MD',trajectory_file_option='DCD')
    
.. note:: Rather than starting QM/MM MD directly like above, it may also be a good idea to equilibrate the system after insertion by e.g.
    a classical NVT simulation of the solution system with solute frozen.


#########################################################################################################
Running QM/MM MD simulations with periodic boundary conditions
#########################################################################################################

Do note that in the QM/MM MD simulation examples above, we typically run the MM-part in a periodic fashion
while the QM-part runs without periodic boundary conditions (since the QM-code xtb does not support PBC).
This is an approximation that is reasonably accurate as long as the QM-system remains in the center of the box (reasonably homogeneous QM-MM interactions).
As the molecule may start to migrate to the box boundaries during a longer simulation it is important to deal with this issue in real-world simulations.
The options are:

1. Run non-periodic simulations by setting periodicity=False in OpenMMTheory.

This may also require freezing parts of the box boundaries to avoid solvent molecules from evaporating or adding a boundary force.
As the solute may migrate towards the frozen boundaries, resulting in artificial dynamics, a centering potential may also be required for the solute.
Non-periodic MM may result in somewhat slower MM-part execution (which may not matter if the QM-code execution is the bottleneck).

2. Stick with periodic MM-part but force the molecule to be in the center of the box.

To make sure of this one can add a flat-bottom centering potential like below. 
See section "Adding a flat-bottom centering potential" in :doc:`module_dynamics` for more information.

.. code-block:: python

  # Here assuming a qm_mm theory, a solute+solvent system have already been created and that solute has indices 0-5.
  # Force constant : 10 kcal/mol/Å^2 Distance from where force acts: 10.0 Å
  MolecularDynamics(theory=qm_mm, fragment=solution, simulation_time=50, timestep=0.001, traj_frequency=100,
    add_centerforce=True, centerforce_atoms=[0,1,2,3,4,5], centerforce_constant=10, centerforce_distance=10.0)



3. Use a QM-code supporting periodic boundary conditions.

Currently CP2K is the primary QM-code in ASH supporting periodic boundary conditions.



#########################################################################################################
Issues
#########################################################################################################

If you get an error like this from OpenMM:

.. code-block:: text

    ValueError: Found multiple NonbondedForce tags with different 1-4 scales

This indicates that there is an incompatibility between the small-molecule XML-file and the water-forcefield XML-file.
Most likely you have selected the wrong XML-file for your solvent in OpenMMTheory. For GAFF and OpenFF you typically want to select the Amber water XML-file:
"amber/tip3p_standard.xml" . For nonbondedFF-only XML-files in CHARMM-style you typically want "charmm36/water.xml".
These files are available globally (if OpenMM is installed) and can be inspected :
dirname $(which python3)
/Users/rb269145/miniconda3/envs/ASH_openmm/lib/python3.11/site-packages/openmm/app/data/

For CHARMM and normal OpenMM XML-files the NonbondedForce line should look like this:

.. code-block:: text

    <NonbondedForce coulomb14scale="1.0" lj14scale="1.0">

where coulomb14scale and lj14scale are set to 1.0

For Amber, GAFF and OpenFF  XML-files the NonbondedForce line should look like this:

.. code-block:: text

  <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">

Additionally, CHARMM XML files contain in addition to NonBondedForce an extra block:

.. code-block:: text

    <LennardJonesForce lj14scale="1.0">