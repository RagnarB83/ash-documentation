======================================
OpenMM interface
======================================

OpenMM is an open-source molecular mechanics library written in C++. It comes with a handy Python interface that was easily adapted for use with ASH. It has been designed for both CPU and GPU codes.



######################################
The OpenMMTheory theory level
######################################

OpenMM can run either on the CPU or on the GPU platform by specifying platform as either: 'CPU', 'OpenCL' or 'CUDA' depending on whether there is a GPU available on the computer. For platform='CPU' the numcores can additionally be specified in order to control the number of CPU cores. Alternatively the shell environment variable OPENMM_CPU_THREADS can be set to control the number of CPU cores that OpenMM will use. If neither numcores keyword is provided or OPENMM_CPU_THREADS variable set, OpenMM will use the number of physical cores present.

The OpenMMTheory class:

.. code-block:: python

    class OpenMMTheory:
        def __init__(self, printlevel=2, platform='CPU', numcores=None, 
                     CHARMMfiles=False, psffile=None, charmmtopfile=None, charmmprmfile=None,
                     GROMACSfiles=False, gromacstopfile=None, grofile=None, gromacstopdir=None,
                     Amberfiles=False, amberprmtopfile=None,
                     xmlfile=None, pdbfile=None, use_parmed=False,
                     do_energy_decomposition=False,
                     periodic=False, charmm_periodic_cell_dimensions=None, customnonbondedforce=False,
                     periodic_nonbonded_cutoff=12, dispersion_correction=True, 
                     switching_function_distance=10,
                     ewalderrortolerance=1e-5, PMEparameters=None,
                     delete_QM1_MM1_bonded=False, applyconstraints=True,
                     autoconstraints=None, hydrogenmass=None, rigidwater=True):




It is possible to read in multiple types of forcefield files: AmberFiles, CHARMMFiles, GROMACSFiles or an OpenMM XML forcefieldfile.
Note: In rare cases OpenMM fails to read in Amber/CHARMM/GROMACS files correctly. In those cases the Parmed library may be more successful (use_parmed=True). Requires ParMed (pip install parmed).

Example creation of an OpenMMtheory object with CHARMM-files:

.. code-block:: python

    forcefielddir="/path/to/dir"
    topfile=forcefielddir+"top_all36_prot.rtf"
    parfile=forcefielddir+"par_all36_prot.prm"
    psffile=forcefielddir+"new-XPLOR-psffile.psf"
    openmmobject = OpenMMTheory(CHARMMfiles=True, psffile=psffile, charmmtopfile=topfile,
                               charmmprmfile=parfile)

Example creation of an OpenMMtheory object with GROMACS-files:

.. code-block:: python

    openmmobject = OpenMMTheory(GROMACSfiles=True, gromacstopdir="/path/to/gromacstopdir",
                    gromacstopfile="gromacstopfile.top", grofile="grofile.gro")

Example creation of an OpenMMtheory object with AMBER files:

.. code-block:: python

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile="/path/to/amberprmtopfile")

Example creation of an OpenMMtheory object with OpenMM XML file:

.. code-block:: python

    openmmobject = OpenMMTheory(xmlfile="example.xml")


An openmmtheory object can then be used to create a QM/MM theory object. See :doc:`module_QM-MM` page.

Periodic boundary conditions:

- If periodic boundary conditions are chosen (periodic=True) then the PBC box parameters are automatically found in the Amber PRMTOP file or the GROMACS Grofile or in the case of CHARMM-files they need to be provided: charmm_periodic_cell_dimensions
- PME parameters can be modified: PMEparameters=[alpha_separation,numgridpoints_X,numgridpoints_Y,numgridpoints_Z] 
- The ewalderrortolerance can be modified (default: 1e-5)
- The periodic nonbonded cutoff can be modified. Default: 12 Å
- Long-range dispersion correction can be turned on or off.
- The switching function distance can be changed. Default: 10 Å. Used for CHARMM and XML files.

######################################
Molecular Dynamics via OpenMM
######################################

It is possible to run MM molecular dynamics of system using the OpenMMTheory object created.
This is accomplished directly via the MD algorithms present in the OpenMM library.
The OpenMM_MD function takes as argument an ASH fragment, an OpenMMTheory object and then the user can select an integrator of choice, simulation temperature, simulation length, timestep, optional additional thermostat, barostat etc.

Most general options available in OpenMM are available in this interface. 
See OpenMM documentation page: http://docs.openmm.org/latest/userguide/application.html#integrators  for details about the integrators, thermostats, barostats etc.

- Available Integrators: Langevin, LangevinMiddleIntegrator, NoseHooverIntegrator, VerletIntegrator, VariableLangevinIntegrator, VariableVerletIntegrator
- Available Barostat: MonteCarloBarostat
- Optional additional thermostat: Anderson

.. code-block:: python

    def OpenMM_MD(fragment=None, openmmobject=None, timestep=0.001, simulation_steps=None, 
    simulation_time=None, traj_frequency=1000, temperature=300, integrator=None, barostat=None, 
    trajectory_file_option='PDB', coupling_frequency=None, anderson_thermostat=False, 
    enforcePeriodicBox=False, frozen_atoms=None):

Options:

- fragment: ASH Fragment object.
- openmmobject: ASH OpenMMTheory object. 
- timestep: float (default: 0.001 ps). Size of timestep in picoseconds.
- simulation_steps: integer. Number of steps to take. (Use either simulation_steps or simulation_time)
- simulation_time: integer. Length of simulation time in ps. (Use either simulation_steps or simulation_time)
- temperature: integer (default:300). Temperature in Kelvin
- integrator: string (regular integrator or integrator+thermostat, e.g. 'LangevinMiddleIntegrator')
- barostat: string (e.g. 'MonteCarloBarostat'). Whether to add barostat to simulation for NPT simulations.
- coupling_frequency: frequency (ps^-1) to update thermostat/integrator. Applies to Nose-Hoover/Langevin.
- anderson_thermostat: Boolean (default: False)
- trajectory_file_option: 'PDB' or 'DCD'. Creates an ASCII PDB-trajectory or a compressed DCD trajectory.
- traj_frequency: integer (default: 1000). How often to write coordinates to trajectory file (every nth step)
- enforcePeriodicBox: Boolean (default: False). Option to fix PBC-image artifacts in trajectory.
- frozen_atoms: list (default: None). What atom indices to freeze in simulation (masses = zero). Note: ASH counts from zero.
- constraints: list of lists (default: None). [[atom_i,atom_j,distance]] Each list defines an atom-pair that is constrained with optionally the bond distance specified.   Example: constraints=[[827,830], [830,833, 1.010]].  Only bond constraints available for now.
- restraints: list of lists (default: None). [[atom_i,atom_j,distance,force_constant]] Example: restraints=[[830,833, 1.010, 1000.0]]  where 830,833 are atom indices, 1.010 is the distance in Å and 1000.0 is the force-constant in kcal/mol \*Å^-2. Only bond restraints available for now.
- autoconstraints: string (default: None) Options: 'HBonds' (X-H bonds constrained), 'AllBonds' (all bonds constrained), 'HAngles' (all bonds and H-X-H and H-O-X angles constrained) or None (default).  Only affects OpenMM_MD runs.
- hydrogenmass: integer (default: None) Mass of hydrogens (e.g. set to 2 for deuterium). Only affects OpenMM_MD runs.


Example:

.. code-block:: python

    from ash import *

    #Forcefield parameters
    forcefielddir="/home/bjornsson/ASH-DEV_GIT/testsuite/OpenMM-files-for-tests/dhfr/charmm/"
    psffile=forcefielddir+"step3_pbcsetup.psf"
    topfile=forcefielddir+"top_all36_prot.rtf"
    prmfile=forcefielddir+"par_all36_prot.prm"

    #Defining fragment
    xyzfile=forcefielddir+"file.xyz"
    frag = Fragment(xyzfile=xyzfile, conncalc=False)

    #Defining OpenMM theory object: CHARMM forcefield with periodic boundary conditions
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=prmfile, periodic=True, charmm_periodic_cell_dimensions=[80, 80, 80, 90, 90, 90],
        dispersion_correction=False, periodic_nonbonded_cutoff=12, switching_function_distance=10,
        PMEparameters=[1.0/0.34, 90, 90, 90])

    #Launching a molecular dynamics simulation
    OpenMM_MD(fragment=frag, openmmobject=openmmobject, timestep=0.001, simulation_steps=20, traj_frequency=1, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')


**General constraints or H-mass modification:**

- In order to allow shorter timesteps in MD simulations it is common to utilize some general constraints in biomolecular simulations, e.g. all X-H bonds, all bonds or even all-bond and some angles. This can be accomplished  via the autoconstraints option (NOTE: an option to OpenMMTheory rather than OpenMM_MD). autoconstraints can be set to: 'HBonds' (X-H bonds constrained), 'AllBonds' (all bonds constrained), 'HAngles' (all bonds and H-X-H and H-O-X angles constrained) or None (default)
- An alternative (or addition) is to change the masses of the hydrogen atoms (fastest-moving atoms). This is also an option to OpenMMTheory. hydrogenmass keyword takes an integer and can e.g. be 2 (mass of deuterium) or heavier. hydrogenmass=None is default (no changes to hydrogen masses).



General X-H constraints and deuterium-mass example:

.. code-block:: python

    from ash import *

    #Forcefield parameters
    forcefielddir="/home/bjornsson/ASH-DEV_GIT/testsuite/OpenMM-files-for-tests/dhfr/charmm/"
    psffile=forcefielddir+"step3_pbcsetup.psf"
    topfile=forcefielddir+"top_all36_prot.rtf"
    prmfile=forcefielddir+"par_all36_prot.prm"

    #Defining fragment
    xyzfile=forcefielddir+"file.xyz"
    frag = Fragment(xyzfile=xyzfile, conncalc=False)

    #Defining OpenMM theory object: CHARMM forcefield with periodic boundary conditions
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=prmfile, periodic=True, charmm_periodic_cell_dimensions=[80, 80, 80, 90, 90, 90], autoconstraints='HBonds', hydrogenmass=2)

    #Launching a molecular dynamics simulation
    OpenMM_MD(fragment=frag, openmmobject=openmmobject, timestep=0.001, simulation_steps=20, traj_frequency=1, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')



Dealing with PBC image problems in trajectory. See https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#how-do-periodic-boundary-conditions-work
To obtain a more pleasing visualization of the trajectory you can "reimage" the trajectory afterwards using the program mdtraj (requires installation of mdtraj: pip install mdtraj)

Example:

.. code-block:: python

    from ash import *
    #Provide trajectory file, PDB topology file and final format of trajectory
    MDtraj_imagetraj("output_traj.dcd", "final_MDfrag_laststep.pdb", format='DCD')
    
    #If periodic box info is missing from trajectory file (can happen with CHARMM files):
    MDtraj_imagetraj("out", pdbtopology, format='DCD', unitcell_lengths=[100.0,100.0,100.0], unitcell_angles=[90.0,90.0,90.0])




######################################
Simple minimization via OpenMM
######################################


Example:

.. code-block:: python

    from ash import *

    #Forcefield parameters
    forcefielddir="/home/bjornsson/ASH-DEV_GIT/testsuite/OpenMM-files-for-tests/dhfr/charmm/"
    psffile=forcefielddir+"step3_pbcsetup.psf"
    topfile=forcefielddir+"top_all36_prot.rtf"
    prmfile=forcefielddir+"par_all36_prot.prm"

    #Defining fragment
    xyzfile=forcefielddir+"file.xyz"
    frag = Fragment(xyzfile=xyzfile, conncalc=False)

    #Defining OpenMM theory object: CHARMM forcefield with periodic boundary conditions
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=prmfile, periodic=True, charmm_periodic_cell_dimensions=[80, 80, 80, 90, 90, 90])

    #Launching a minimization
    OpenMM_Opt(fragment=frag, openmmobject=openmmobject, maxiter=1000, tolerance=1)
    #After minimization, the ASH fragment is updated, a PDB-file is written out: frag-minimized.pdb
    #Alternative XYZ write-out:
    frag.write_xyzfile(xyzfilename="frag_afteropt.xyz")


If you want to do a simple minimization of only the H-atoms of your system (e.g. your protein with newly added H-atoms),
you can do this by freezing all non-H atoms. An ASH fragment can conveniently give you lists of atom indices by the built-in functions:

- fragment.get_atomindices_for_element('C')   #List of atom-indices for carbon atoms in the system
- fragment.get_atomindices_except_element('H')   #List of atom-indices for all atoms except the chosen element (here H).

Note: all constraints in the OpenMM object needs to be turned off for (autoconstraints=None, rigidwater=False) for this many frozen atoms (frozen atoms can not have constraints).

.. code-block:: python

    from ash import *

    numcores = 4

    pdbfile = "ash_inp.pdb"
    prmtopfile = "prmtop"

    frag = Fragment(pdbfile=pdbfile)

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile=prmtopfile, periodic=True,
            platform='CPU', autoconstraints=None, rigidwater=False)

    #List of all non-H atoms
    allnonHatoms=frag.get_atomindices_except_element('H')

    OpenMM_MD(fragment=frag, openmmobject=openmmobject, timestep=0.001, simulation_steps=100,
            traj_frequency=1, temperature=300, integrator="LangevinIntegrator",
            coupling_frequency=1, trajectory_file_option="PDB", frozen_atoms=allnonHatoms,)


######################################
System setup via OpenMM: Modeller
######################################

OpenMM features a convenient PDBfixer program (https://github.com/openmm/pdbfixer) and a Modeller tool (http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.modeller.Modeller.html)
that is capable of setting up a new biomolecular system from scratch. See also: http://docs.openmm.org/7.2.0/userguide/application.html#model-building-and-editing
ASH features a highly convenient interface to these programs and allows near-automatic system-setup for favorable systems.

.. code-block:: python

    def OpenMM_Modeller(pdbfile=None, forcefield=None, xmlfile=None, waterxmlfile=None, watermodel=None, pH=7.0, 
                    solvent_padding=10.0, solvent_boxdims=None, extraxmlfile=None, residue_variants=None,
                    ionicstrength=0.1, iontype='K+'):



Lysozyme example (simple, no modifications required):

.. code-block:: python

    from ash import *

    #Original raw PDB-file (no hydrogens, nosolvent)
    #Download from https://www.rcsb.org/structure/1AKI
    pdbfile="1aki.pdb"


    #Defining residues with special user-wanted protonation states
    #Example: residue_variants={0:'LYN', 17:'CYX', 18:'ASH', 19:'HIE' } 
    #residue 0 neutral LYS, residue 17, deprotonated CYS, residue 18 protonated ASP, residue 19 epsilon-protonated HIS.
    residue_variants={}

    #Setting up new system, adding hydrogens, solvent, ions and defining forcefield, topology
    forcefield, topology, ashfragment = OpenMM_Modeller(pdbfile=pdbfile, forcefield='CHARMM36', watermodel="tip3p", pH=7.0, 
        solvent_padding=10.0, ionicstrength=0.1, iontype="Na+", residue_variants=residue_variants)

    #Creating new OpenMM object from forcefield, topology and and fragment
    openmmobject =OpenMMTheory(platform='CPU', numcores=numcores, Modeller=True, forcefield=forcefield, topology=topology,
                     pdbfile=None, do_energy_decomposition=True, periodic=True,
                     autoconstraints='HBonds', rigidwater=True)

    #MM minimization for 100 steps
    OpenMM_Opt(fragment=ashfragment, openmmobject=openmmobject, maxiter=100, tolerance=1)

    #Classical MD simulation for 10 ps
    OpenMM_MD(fragment=ashfragment, openmmobject=openmmobject, timestep=0.001, simulation_time=10, traj_frequency=100, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')



If the protein contains nonstandard residues (e.g. metalcofactors) that are not present in a typical protein forcefield (OpenMM_Modeller will exit with errors),
then these need to be provided using the extraxmlfile option.

.. code-block:: python

    forcefield, topology, ashfragment = OpenMM_Modeller(pdbfile=pdbfile, forcefield='CHARMM36', watermodel="tip3p", pH=7.0, 
        solvent_padding=10.0, ionicstrength=0.1, iontype="Na+", residue_variants=residue_variants, extraxmlfile="cofactor.xml")


The cofactor.xml file needs to define a forcefield (a nonbonded one at least) for the residue. 
Here defining a dummy molybdenum ion:

.. code-block:: 

    <ForceField>
    <AtomTypes>
    <Type name="MOX" class="Mo" element="Mo" mass="99.0"/>
    </AtomTypes>
    <Residues>
    <Residue name="FEM">
    <Atom name="MOD" type="MOX"/>
    </Residue>
    </Residues>
    <NonbondedForce coulomb14scale="1.0" lj14scale="1.0">
    <Atom type="MOX" charge="3" sigma="0.375" epsilon="0.439"/>
    </NonbondedForce>
    </ForceField>


Advanced example (additional forcefield parameters required):

.. code-block:: python

    from ash import *


    #TODO




