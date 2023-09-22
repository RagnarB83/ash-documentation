Explicit solvation (small molecule)
======================================

ASH allows you to easily create explicit solvation models for small molecules that can then be either
used either for classical MD simulations or QM/MM MD simulations.

Current limitations:

- Only water solvent possible for now

################################################################################################
Creating a small-molecule forcefield
################################################################################################

The initial step will be to acquire a forcefield for the small molecule of interest.
ASH features 2 main options:

**write_nonbonded_FF_for_ligand**

Creates a simple nonbonded forcefield. Useful for inorganic molecules (e.g. metal complexes) where a full
forcefield is hard to create. Good option for QM/MM MD simulations. Classical MD simulations will require the internal
degrees of freedom of the molecule to be constrained.

**small_molecule_parameterizor**

Creates a full forcefield of the molecule with bonded and nonbonded parameters. Works when molecule is organic or drug-like.
Can either use a GAFF or OpenFF forcefield.


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
    small_molecule_parameterizor(xyzfile="isobutyraldehyde.xyz", forcefield_option="OpenFF")

The function will create a file, called openff_ligand.xml.

2. Solvate the system

*solvate.py:*

.. code-block:: python
        
    from ash import *

    mol = Fragment(xyzfile="isobutyraldehyde.xyz", charge=0, mult=1)
    #Solvate molecule in a 30x30x30 Å TIP3P water box.
    solvate_small_molecule(fragment=mol, xmlfile="openff_ligand.xml", watermodel='tip3p', solvent_boxdims=[30,30,30])

This creates files: system_aftersolvent.xyz and system_aftersolvent.pdb

3. Perform initial classical MD simulation to warmup and equilibrate 

*run_mm_md.py:*

.. code-block:: python

    from ash import *

    #Read in coordinates of the full system
    pdbfile="system_aftersolvent.pdb"
    fragment = Fragment(pdbfile=pdbfile)
    #Create an OpenMMTheory object based on PDB-file and XML-files for water and small-molecule
    omm =OpenMMTheory(xmlfiles=["openff_ligand.xml", "amber/tip3p_standard.xml"],
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
    omm =OpenMMTheory(xmlfiles=["openff_ligand.xml", "amber/tip3p_standard.xml"],
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
        temperature=30, platform='OpenCL', integrator='LangevinMiddleIntegrator', coupling_frequency=1, 
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

The function will create a file, here called: "LIG.xml". By default it uses the ff_type to be "CHARMM" which means that the XML-file will
use a form of the nonbonded potential that is compatible with CHARMM protein forcefield in general. This can be changed to "AMBER" or "None" if required.


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
        temperature=30, platform='OpenCL', integrator='LangevinMiddleIntegrator', coupling_frequency=1, 
        trajfilename='QM_MM_NVT-MD',trajectory_file_option='DCD')
    







#########################################################################################################
Issues
#########################################################################################################

If you get an error like this from OpenMM:

.. code-block:: text

    ValueError: Found multiple NonbondedForce tags with different 1-4 scales

This indicates that there is an incompatibility between the small-molecule XML-file and the water-forcefield XML-file.
Most likely you have selected the wrong XML-file for your solvent in OpenMMTheory. For GAFF and OpenFF you typically want to select the Amber water XML-file:
amber/tip3p_standard.xml . For nonbondedFF-only XML-files in CHARMM-style you typically want "charmm36/water.xml".

For CHARMM and normal OpenMM XML-files the NonbondedForce line should look like this:

.. code-block:: text

    <NonbondedForce coulomb14scale="1.0" lj14scale="1.0">

where coulomb14scale and lj14scale are set to 1.0

For Amber, GAFF and OpenFF  XML-files the NonbondedForce line should look like this:

  <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">

Additionally, CHARMM XML files contain in addition to NonBondedForce an extra block:

.. code-block:: text

    <LennardJonesForce lj14scale="1.0">