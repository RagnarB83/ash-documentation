======================================
OpenMM interface
======================================



######################################
The OpenMMTheory theory level
######################################

The OpenMMTheory class:

.. code-block:: python

    class OpenMMTheory:
        def __init__(self, pdbfile=None, platform='CPU', active_atoms=None, frozen_atoms=None,
                     CHARMMfiles=False, psffile=None, charmmtopfile=None, charmmprmfile=None,
                     GROMACSfiles=False, gromacstopfile=None, grofile=None, gromacstopdir=None,
                     Amberfiles=False, amberprmtopfile=None, printlevel=2, do_energy_composition=False,
                     xmlfile=None, periodic=False, periodic_cell_dimensions=None, customnonbondedforce=False,
                     delete_QM1_MM1_bonded=False, watermodel=None, use_parmed=False, periodic_nonbonded_cutoff=12,
                     dispersion_correction=True, switching_function=False, switching_function_distance=10,
                     ewalderrortolerance=1e-5, applyconstraints=False, PMEparameters=None):


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


######################################
Molecular Dynamics via OpenMM
######################################

It is possible to run MM molecular dynamics of system using the OpenMMTheory object created.
This is accomplished directly via the MD algorithms present in the OpenMM library.
The OpenMM_MD function takes as argument an ASH fragment, an OpenMMTheory object and then the user can select an integrator of choice, simulation temperature, simulation length, timestep, optional additional thermostat, barostat etc.

Most general options available in OpenMM are available in this interface. 
See OpenMM documentation page: http://docs.openmm.org/latest/userguide/application.html#integrators  for details about the integrators, thermostats, barostats etc.

- Available Integrators: LangevinMiddleIntegrator, NoseHooverIntegrator, VerletIntegrator, VariableLangevinIntegrator, VariableVerletIntegrator
- Available Barostat: MonteCarloBarostat
- Optional additional thermostat: Anderson

.. code-block:: python

    def OpenMM_MD(fragment=None, openmmobject=None, timestep=0.001, simulation_steps=None, simulation_time=None,
     temperature=300, integrator=None, coupling_frequency=None, barostat=None, anderson_thermostat=None, 
     trajectory_file_option='PDB', traj_frequency=1000):

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







