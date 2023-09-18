Explicit solvation (small molecule)
======================================

ASH allows you to easily create explicit solvation models for small molecules that can then be either
used either for classical MD simulations or QM/MM MD simulations.

Current limitations:

- Only water solvent possible for now

################################################################################################
Example 1. Modelling of an organic molecule in explicit water with a ligand forcefield
################################################################################################

ASH allows one to easily parameterize a full forcefield for an organic molecule using either GAFF or OpenFF
models via **small_molecule_parameterizor**

.. code-block:: python
        
    from ash import *

    #Parameterize small molecule using OpenFF
    small_molecule_parameterizor(xyzfile="isobutyraldehyde.xyz", forcefield_option="OpenFF")

The function will created a file, called openff_ligand.xml.


#########################################################################################################
Example 2. Modelling of an inorganic molecule in explicit water using a simple non-bonded forcefield
#########################################################################################################

For inorganic molecules (e.g. metal complexes) it is trickier to parameterize a full forcefield and would
have to be performed to some extent manually. However, it is possible to use a simple non-bonded forcefield,
i.e. without bonded parameters. 
This then requires classical simulations to be performed with internal degrees of freedom frozen 
(bonds, angles, dihedrals) while QM/MM MD simulations can be performed as normal.




.. code-block:: python

    from ash import *

    numcores=4

    #Read in coordinates for small molecule
    mol=Fragment(xyzfile="3fgaba.xyz")
    mol.charge=0;mol.mult=1

    #Solvate molecule (70x70x70 Å TIP3P water box), similar to OpenMM_Modeller
    forcefield, topology, ashfragment = solvate_small_molecule(fragment=mol, charge=mol.charge, mult=mol.mult, watermodel='tip3p', 
        solvent_boxdims=[70,70,70], nonbonded_pars="CM5_UFF", numcores=numcores)

    #Creating new OpenMM object from forcefield, topology and and fragment
    soluteatoms=[i for i in range(0,mol.numatoms)]
    openmmobject =OpenMMTheory(platform='OpenCL', numcores=numcores, forcefield=forcefield, topoforce=True, topology=topology, periodic=True, frozen_atoms=soluteatoms, 
        autoconstraints='HBonds', rigidwater=True)


    #MM minimization for 100 steps
    OpenMM_Opt(fragment=ashfragment, theory=openmmobject, maxiter=100, tolerance=1, enforcePeriodicBox=True)

    #Classical MD simulation for 100 ps
    OpenMM_MD(fragment=ashfragment, theory=openmmobject, timestep=0.001, simulation_time=100, traj_frequency=10, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD', enforcePeriodicBox=True)

    #MDAnalysis transforming of trajectory
    lastframe_elems, lastframe_coords = MDAnalysis_transform("final_MDfrag_laststep.pdb","trajectory.dcd", solute_indices=soluteatoms, 
        trajoutputformat='PDB')

    #new ASH fragment after the classical prep
    system_after_classical_prep=Fragment(elems=lastframe_elems, coords=lastframe_coords)
    system_after_classical_prep.write_xyzfile(xyzfilename="Finalframe_system.xyz")

