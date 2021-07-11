======================================
OpenMM interface
======================================

OpenMM is an open-source molecular mechanics library written in C++. It comes with a handy Python interface that was easily adapted for use with ASH. It has been designed for both CPU and GPU godes.



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
                     delete_QM1_MM1_bonded=False, applyconstraints=False):




It is possible to read in multiple types of forcefield files: AmberFiles, CHARMMFiles, GROMACSFiles or OpenMMXMLFile.
Note: In rare cases OpenMM fails to read in the files correctly. In those cases the Parmed library may be more successful (use_parmed=True). Requires ParMed (pip install parmed).

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

    openmmobject = OpenMMTheory(xmlfile="exampl.xml")


An openmmtheory object can then be used to create a QM/MM theory object. See :doc:`module_QM-MM` page.
If platform is 'CPU' (default) then the number of cores can be controlled by the numcores keyword (highest preference) or alternatively by the shell variable OPENMM_CPU_THREADS. If neither option is set then OpenMM will use the number of physical cores present (not always desired).


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

    def OpenMM_MD(fragment=None, openmmobject=None, timestep=0.001, simulation_steps=None, simulation_time=None, traj_frequency=1000, temperature=300, integrator=None, barostat=None, trajectory_file_option='PDB', coupling_frequency=None, anderson_thermostat=False, enforcePeriodicBox=False, frozen_atoms=None):

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
        charmmprmfile=prmfile, periodic=True, periodic_cell_dimensions=[80, 80, 80, 90, 90, 90],
        dispersion_correction=False, periodic_nonbonded_cutoff=12, switching_function_distance=10,
        PMEparameters=[1.0/0.34, 90, 90, 90])

    #Launching a molecular dynamics simulation
    OpenMM_MD(fragment=frag, openmmobject=openmmobject, timestep=0.001, simulation_steps=20, traj_frequency=1, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')







