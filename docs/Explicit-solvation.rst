Explicit solvation (small molecule)
======================================

ASH allows you to easily create explicit solvation models for small molecules that can then be combined with QM/MM to perform QM/MM molecular dynamics, QM/MM metadynamics as well as extracting of snapshots for single-point energy/property calculations.

Current limitations:

- Only water solvent possible for now

################################################################################################
Example 1. Solvation, MM minimization and classical MM dynamics for 100 ps.
################################################################################################

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


################################################################################################
Example 2. Solvation, MM minimization, MM MD and QM/MM MD
################################################################################################

Not ready yet