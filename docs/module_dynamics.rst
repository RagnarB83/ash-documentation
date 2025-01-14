
Molecular dynamics
======================================

Molecular dynamics in ASH is performed with the help of MD algorithms from the OpenMM library.
It requires the OpenMM library to be installed.

This approach not only allows MM dynamics, but dynamics at any ASH level of theory (QM and QM/MM etc.). 
For non-classical MD, ASH uses an OpenMM CustomExternalForce object to communicate forces from the ASH theory to a dummy OpenMM system.

Workflows to perform enhanced or biased sampling molecular dynamics is also available. See :doc:`Biased-sampling`.

######################################################
MolecularDynamics (via OpenMM library routines)
######################################################

Dynamics are performed by calling the **MolecularDynamics** function (alias for the **OpenMM_MD** function)

- Available Integrators: Langevin, LangevinMiddleIntegrator, NoseHooverIntegrator, VerletIntegrator, VariableLangevinIntegrator, VariableVerletIntegrator
- Available Barostat: MonteCarloBarostat
- Optional additional thermostat: Anderson

See `OpenMM documentation page <http://docs.openmm.org/latest/userguide/application.html#integrators>`_  for details about the integrators, thermostats, barostats etc.

The **OpenMM_MD** function (also called **MolecularDynamics**) takes as argument an ASH fragment, a theory object and the user then selects an integrator of choice, 
simulation temperature, simulation length, timestep, optional additional thermostat, barostat and possible other options.
The theory level can be **OpenMMTheory**, **QMMMTheory** or even a **QMTheory**.

.. code-block:: python

    #MolecularDynamics is alias for OpenMM_MD
    def OpenMM_MD(fragment=None, theory=None, timestep=0.004, simulation_steps=None, simulation_time=None,
                traj_frequency=1000, temperature=300, integrator='LangevinMiddleIntegrator',
                barostat=None, pressure=1, trajectory_file_option='DCD', trajfilename='trajectory',
                coupling_frequency=1, charge=None, mult=None,
                anderson_thermostat=False,
                enforcePeriodicBox=True, dummyatomrestraint=False, center_on_atoms=None, solute_indices=None,
                datafilename=None, dummy_MM=False, plumed_object=None, add_center_force=False,
                center_force_atoms=None, centerforce_constant=1.0, barostat_frequency=25, specialbox=False):


**OpenMM_MD** options:

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``fragment``
     - ASH Fragment
     - None
     - The ASH fragment.
   * - ``theory``
     - ASH Theory
     - None
     - The ASH Theory object.
   * - ``timestep``
     - float
     - 0.004
     - | The timestep . Default: 0.004 ps (suitable for LangevinMiddleIntegrator 
       | dynamics with frozen X-H bonds)
   * - ``simulation_steps``
     - integer
     - None
     - Number of simulation steps to take. Alternative to simulation_time below
   * - ``simulation_time``
     - float
     - None
     - Length of simulation in picoseconds. Alternative to simulation_time above.
   * - ``temperature``
     - integer
     - 300
     - The temperature in Kelvin.
   * - ``integrator``
     - string
     - LangevinMiddleIntegrator
     - | The integrator to use. Options: 'Langevin', 'LangevinMiddleIntegrator', 
       | 'NoseHooverIntegrator', 'VerletIntegrator', 'VariableLangevinIntegrator',
       | 'VariableVerletIntegrator'
   * - ``coupling_frequency``
     - integer
     - 1
     - The coupling frequency of thermostat (in ps^-1 for Nosé-Hoover and Langevin-type)
   * - ``barostat``
     - string
     - None
     - Barostat to use for NPT simulations. Options: 'MonteCarloBarostat'
   * - ``barostat_frequency``
     - int
     - 25
     - Frequency of barostat update.
   * - ``pressure``
     - int
     - 1
     - Pressure to enforce by barostat
   * - ``anderson_thermostat``
     - Boolean
     - False
     - Whether to use Andersen Thermostat
   * - ``trajectory_file_option``
     - string
     - 'DCD'
     - | Type of trajectory file. Options: 'DCD' (compressed), 'PDB', 'NetCDFReporter' 
       | (compressed), 'HDF5Reporter' (compressed). Applies only to pure MM simulations.
   * - ``traj_frequency``
     - integer
     - 1000
     - Frequency of writing trajectory (every Xth timestep).
   * - ``trajfilename``
     - string
     - None
     - Name of trajectory file (without suffix).
   * - ``enforcePeriodicBox``
     - Boolean
     - True
     - Enforce PBC image during simulation. Fixes PBC-image artifacts in trajectory.
   * - ``center_on_atoms``
     - list
     - None
     - Expert options: Center system on these atoms.
   * - ``dummyatomrestraint``
     - Boolean
     - False
     - Expert options: Dummy atom restraints.
   * - ``solute_indices``
     - list
     - None
     - Expert options: solute_indices
   * - ``add_center_force``
     - Boolean
     - False
     - Whether to add a spherical force that pushes atoms to the center.
   * - ``center_force_atoms``
     - list
     - None
     - List of atom indices that the center force acts on.
   * - ``centerforce_constant``
     - float
     - None
     - Value of the spherical center force in kcal/mol/Ang^2.
   * - ``specialbox``
     - Boolean
     - False
     - Expert option: Special box for QM/MM.
   * - ``plumed_object``
     - ASH-Plumed object
     - None
     - Expert option: Plumed object for biased dynamics.


######################################################
Examples: how to use
######################################################

The simplest way to run MD is to simply call the **OpenMM_MD** function, also named **MolecularDynamics**.

.. code-block:: python
  
  #Fragment and Theory object creation not shown here
  #OpenMM_MD(fragment=frag, theory=theory, timestep=0.001, simulation_time=2)
  MolecularDynamics(fragment=frag, theory=theory, timestep=0.001, simulation_time=2)


**Defing an object from class instead of using the function**

The **OpenMM_MD** (MolecularDynamics) function is actually a wrapper around a class named **OpenMM_MDclass**.
Sometimes more flexiblity can be achieved by using the class instead of calling the function and this can be useful when developing
more complex MD-based workflows in ASH.

Below we see create an MD-object and run it using the run method instead of calling **MolecularDynamics**.
The 2 approaches are equivalent, the simple functiona-call approach is simpler and typically recommended, the latter is more explicit and can be useful (see use-case later in restraints section).

.. code-block:: python
  
  #Fragment and Theory object creation not shown here
  #Alternative: Creating object from class and running
  md_object = OpenMM_MDclass(fragment=frag, theory=theory, timestep=0.001)
  md_object.mdobj.run(simulation_time=2.0)

**Pure MM example:**

For pure classical MM-based MD the dynamics are run by calling directly the MD routines in the OpenMM library on an OpenMM system.
This approach results in no data transfer between the C++/CUDA/OpenCL layer (of OpenMM) and Python layer (of OpenMM and ASH) while running,
and thus is just as fast as using OpenMM directly.
Running OpenMM via the GPU code (if a GPU is available) will allow particularly fast MM dynamics.
See :doc:`OpenMM-interface` for additional details on running classical MD simulations on OpenMMTheory objects.

.. code-block:: python

	from ash import *

	#Fragment from XYZ-file
	frag=Fragment(xyzfile="frag.xyz")
	#Defining frozen region. Taking difference of all-atom region and active region
	actatoms=[14,15,16]
	frozen_atoms=listdiff(frag.allatoms,actatoms)
	#Defining OpenMM object 
	openmmobject = OpenMMTheory(cluster_fragment=frag, ASH_FF_file="Cluster_forcefield.ff", frozen_atoms=frozen_atoms)
	#Calling OpenMM_MD function
	MolecularDynamics(fragment=frag, theory=openmmobject, timestep=0.001, simulation_time=2, traj_frequency=10, temperature=300,
	    integrator='LangevinIntegrator', coupling_frequency=1)


**QM/MM example:**

For a QM/MM system that utilizes OpenMMTheory as mm_theory and any QM-theory as qm_theory, it is also possible to use **MolecularDynamics** to do QM/MM dynamics. 
In this case the QM+PC gradient is used to update the forces of the OpenMM system (as a CustomExternalForce).
While some data transfer between the OpenMM C++/CUDA/OpenCL and ASH Python layers is necessary, the performance is still very good
and should be negligible compared to the speed of the QM-step, which should always be the bottleneck.
See :doc:`OpenMM-interface` for details.

.. code-block:: python

	from ash import *

	#Fragment
	frag=Fragment(xyzfile="frag.xyz")
	#Defining frozen region. Taking difference of all-atom region and active region
	actatoms=[14,15,16]
	frozen_atoms=listdiff(frag.allatoms,actatoms)

	xtbtheory = xTBTheory(runmode='inputfile', xtbmethod='GFN2', numcores=numcores)
	openmmobject = OpenMMTheory(cluster_fragment=frag, ASH_FF_file="Cluster_forcefield.ff", frozen_atoms=frozen_atoms)
	QMMMTheory = QMMMTheory(fragment=frag, qm_theory=xtbtheory, mm_theory=openmmobject,
		qmatoms=qm_region, embedding='Elstat', numcores=numcores)

	MolecularDynamics(fragment=frag, theory=QMMMTheory, timestep=0.001, simulation_time=2, traj_frequency=10, temperature=300,
	    integrator='LangevinIntegrator', coupling_frequency=1, charge=0, mult=1)

**QM example:**

It is even possible to use the dynamics routines of the OpenMM library to drive an MD simulation at the QM-level. 
This is possible by automatic creation of a dummy OpenMM system, containing no MM force except CustomExternalForce,
which is then updated directly by the gradient resulting from running the ASH-Theory object.


.. code-block:: python

	from ash import *
	
	numcores=12
	#Simple n-butane system
	butane=Fragment(xyzfile="butane.xyz", charge=0, mult=1)

	# Creating xTBTheory object (Note: runmode='library' runs faster) that is parallelized. Using GFN1-xTB.
	xtbcalc = xTBTheory(xtbmethod='GFN1', runmode='library', numcores=numcores)
	
	#Running NVE dynamics (initial temp=300 K) on butane using xTBTheory.
	# 0.001 ps timestep, 2 ps , writing every 10th step to trajectory. A velocity Verlet algorithm is used.
	MolecularDynamics(fragment=butane, theory=xtbcalc, timestep=0.001, simulation_time=2, traj_frequency=10, temperature=300,
		integrator='LangevinIntegrator', coupling_frequency=1)

######################################################
Running MD in different ensembles
######################################################

Molecular dynamics can be performed in either a NVE, NVT or NPT ensemble. 
Which ensemble is simulated, depends on whether a thermostat-integrator and/or barostat is present.
The thermostat is specified by the *integrator* keyword (i.e. one should choose an integrator which is also a thermostat).

**NVT ensemble (default)**

If only a integrator-thermostat is present, a NVT molecular dynamics simulation is run. 
This is the default.
A `LangevinMiddleIntegrator <http://docs.openmm.org/latest/api-python/generated/openmm.openmm.LangevinMiddleIntegrator.html?highlight=langevinmiddle#openmm.openmm.LangevinMiddleIntegrator>`_  
integrator-thermostat is automatically selected, with a default temperature of 300 K for the heatbath.
A friction coefficient (coupling frequency with heat bath) of 1 inverse picosecond is the default and can be changed by the  **coupling_frequency** keyword.

Specifying a different thermostat, temperature and friction can be done like this:

.. code-block:: python

  #NVT ensemble. Changing the the integrator to Nose-Hoover here
  MolecularDynamics(fragment=butane, theory=xtbcalc, timestep=0.001, simulation_time=2, 
        integrator='NoseHooverIntegrator', temperature=300, coupling_frequency=1)

**NPT ensemble**

Running a NPT simulation requires having both a thermostat (see above) and a barostat present. 
The barostat is added by specifying the barostat keyword in the **MolecularDynamics** function which adds a
`MonteCarloBarostat <http://docs.openmm.org/latest/api-python/generated/openmm.openmm.MonteCarloBarostat.html?highlight=montecarlo#openmm.openmm.MonteCarloBarostat>`_ to the OpenMM object.
A pressure of 1 bar is the default and the default barostat coupling frequency is every 25 steps.

.. code-block:: python

  #NPT ensemble with a LangevinMiddleIntegrator integrator-thermostat and a Monte Carlo barostat
  MolecularDynamics(fragment=butane, theory=xtbcalc, timestep=0.001, simulation_time=2, 
        integrator='LangevinMiddleIntegrator', temperature=300, coupling_frequency=1,
        barostat=None, pressure=1, barostat_frequency=25)

NPT simulations are often performed for the purpose of equilibrating the volume of a system under applied pressure and temperature.
The periodic box vectors (and hence volume of box) is allowed to fluctuate during the simulation.
For classical MD simulations (i.e. theory is OpenMMTheory), ASH features a workflow called **OpenMM_box_equilibration** (see :doc:`OpenMM-interface`)
for conveniently equilibrating a system until convergence is reached in volume of the box and density.

**NVE ensemble**

If no integrator-thermostat or barostat is present, a NVE simulation will be run.
Below we specify no barostat and we additionally change the integrator keyword to be "VerletIntegrator" 
which corresponds to a leap-frog Verlet algorithm (instead of the default LangevinMiddleIntegrator which is an integrator+thermostat).

.. code-block:: python

  #NVE ensemble. Changing the the integrator to VerletIntegrator here
  MolecularDynamics(fragment=butane, theory=xtbcalc, timestep=0.001, simulation_time=2, 
        integrator='VerletIntegrator')

######################################################
Using general constraints in MD simulations
######################################################

For classical MD simulations as well as QM/MM MD simulation it is common to use bond constraints during the simulation.
Common waterforcefields (e.g. TIP3P) are typically designed to be completely rigid and it is common in biomolecular simulations to contrains all X-H bonds (where X is a heavy atom like e.g. C).

Constraints during MD can be implemented in a few ways:

**Automatic XH constraints with rigidwater in OpenMMTHeory**

By using the autoconstraints keyword (options: 'HBonds', 'AllBonds', 'HAngles') in **OpenMMTheory** one can constrain the XH-bonds ('HBonds'), all bonds ('AllBonds') or all-bonds + all angles ('HAngles').
Furthermore the rigidwater keyword (True or False) sets the constraints for water molecules if present in the system. The autoconstraints keyword should work to add constraints for many forcefields but may fail to
work for some manually defined forcefields.

.. code-block:: python

  #OpenMMTheory object with constraints enabled.
  openmmobject = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "specialresidue.xml"],
                    autoconstraints='HBonds', rigidwater=True)

**Defining constraints manually or semi-automatically**

OpenMMTheory also allows bond-constraints to be added manually by providing a list of lists of pairs of atomindices.

.. code-block:: python

  #Manual definition of the constraints list-of-lists: here constraining bond between atoms 17 & 18 as well as 235 & 236
  con_list = [[17,18], [235,236]]
  #OpenMMTheory object with constraints enabled.
  openmmobject = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "specialresidue.xml"],
                    bondconstraints=con_list)

To avoid manually defining a long list of constraints for a large system it is possible to use the **define_XH_constraints** function
to conveniently generate this list of constraints for either the full system or a subset of it (e.g. an active-region). 

.. code-block:: python

  frag = Fragment(pdbfile="system.pdb")
  # Automatically define a list of lists of all X-H constraints (constraints for all bonds involving H atoms) for all system atoms (actatoms=frag.allatoms)
  con_list = define_XH_constraints(frag, actatoms=frag.allatoms)
  # Same but with the option of avoiding constraints for a group of atoms (e.g. a QM-region)
  qmatoms= [17,18,19,20]
  con_list = define_XH_constraints(frag, actatoms=frag.allatoms, excludeatoms=qmatoms)
  #OpenMMTheory object with constraints enabled.
  openmmobject = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "specialresidue.xml"],
                    bondconstraints=con_list)

######################################################
Adding custom forces to MD simulation
######################################################

Sometimes it is desirable to add additional forces or potentials to and MD simulation for the purposes of avoiding or enhancing certain behaviour.
For example in enhanced sampling methods like metadynamics or umbrella sampling we add biasing potentials to steer the system away from certain regions.
See :doc:`Biased-sampling` for more info on how to perform metadynamics or umbrella sampling in ASH.

ASH allows a few different ways of adding additional forces to a system.

.. note::
  If you want to add a force in QM-based MD simulation where no **OpenMMTheory** object has been previously defined, 
  you may first have to create the MD-simulation object from the **OpenMM_MDclass** . You can then access the OpenMMTheory object 
  as the *openmmobject* attribute of the MD-simulation object.


**Adding an internal-coordinate restraint potential**

Internal-coordinate restraint potentials are e.g. used in umbrella sampling simulations.

It is easy to add internal-coordinate restraints to an MD simulation and can be accomplished in 2 ways:
i ) by providing list of restraints to the **MolecularDynamics** function or ii) by calling the 
**add_custom_bond_force**, **add_custom_angle_force** or **add_custom_torsion_force** methods of the OpenMMTheory object.


The recommended way is to use the *restraints* keyword in the **MolecularDynamics** function.
This option allows one to provide a list of restraints (as a list of lists) where each restraint is defined by the following syntax:

- bondrestraint: [atom1, atom2, equilibrium_distance, forceconstant]
- anglerestraint: [atom1, atom2, atom3, equilibrium_angle, forceconstant]
- torsionrestraint: [atom1, atom2, atom3, atom4, equilibrium_torsion, forceconstant]

.. warning::
  Angles and torsions should be in radian units while bonds/distances are in Ångström.

Example below defines a bond-restraint between atom 0 and atom 3 in butane with an equilibrium value of 2.5 and 
a force constant of 100.0 kcal/mol/Å^2.

.. code-block:: python

  from ash import *
  frag = Fragment(databasefile="butane.xyz")
  xtbcalc = xTBTheory()
  MolecularDynamics(fragment=frag, theory=xtbcalc, timestep=0.001, simulation_time=2, traj_frequency=10,
                      restraints=[[0,3,2.5,100.0]])

Alternatively, one can also add restraints by calling the respective method of the **OpenMMTheory** object.


.. code-block:: python

  # Bond,angle and torsion restraint methods inside OpenMMTheory
  def add_custom_bond_force(self,i,j,value,forceconstant):
  def add_custom_angle_force(self,i,j,k,value,forceconstant):
  def add_custom_torsion_force(self,i,j,k,l,value,forceconstant):

If an OpenMMTheory object is already defined then we can add a restraint potential acting on an internal coordinate of the system.

.. code-block:: python

  #Assuming an OpenMMTheory object called omm already defined

  #Adding a bond/distance restraint between atom 5 and 11
  #equilibrium distance of 3.0 and force-constant of 5.0 kcal/mol/Å^2
  omm.add_custom_bond_force(5,11,3.0,5.0)
  #Adding an angle restraint defined by atoms 5,11,12
  #equilibrium angle of 3.0 radians and force-constant of 4.0 kcal/mol/Å^2
  omm. add_custom_angle_force(5,11,12,100.0,4.0)
  #Adding a torsion restraint between atoms 5,11,12,13
  #equilibrium torsion of 3.3 radians and a force-constant of 3.0 kcal/mol/Å^2
  omm.add_custom_torsion_force(5,11,12,13,3.0,3.0)


Below is an example for a QMTheory where we first define the OpenMM_MD object (instead of calling MolecularDynamics/OpenMM_MD ) 
and then access the OpenMMTheory object inside  (necessary since we didn't have an OpenMMTheory object before).
We can then use the **add_custom_bond_force** method inside to directly define the restraint.

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="butane.xyz")
  xtbcalc = xTBTheory()

  #Since QM-theory we have to first define the OpenMM_MD object
  mdobj = OpenMM_MDclass(fragment=frag, theory=xtbcalc, timestep=0.001, traj_frequency=10)
  #Then call the method inside the internal OpenMMTheory object
  mdobj.openmmobject.add_custom_bond_force(0,3,2.5,100.0)
  #Then run the simulation
  mdobj.run(simulation_time=2.0)

**Adding a funnel restraint**

Funnel restraints are e.g. used in funnel metadynamics.

NOT READY YET!

**Adding a custom-external force (EXPERT)**

It is possible to add a user-specific custom force to an OpenMMTheory object by using the :
*add_custom_external_force* and *update_custom_external_force* methods. This functionality is used to run MD simulations using non-MM theories in ASH (QM and QM/MM).

This feature can be used to create a manual MD-simulation like below:

.. code-block:: python

  from ash import *
  import numpy as np

  frag = Fragment(databasefile="h2o.xyz")
  #Create an OpenMMTheory object. Here doing one without any forces (dummysystem option)
  omm = OpenMMTheory(fragment=frag, dummysystem=True)
  #Add a custom force for each particle of the system. Initial forces at value zero.
  customforce = omm.add_custom_external_force() #Updates omm.system, force also returned
  #Create simulation 
  simulation = omm.create_simulation() #Create a simulation object

  #Initial positions (Numpy array in Å)
  omm.set_positions(frag.coords,simulation)

  # Run simulation step by step
  for i in range(100):

    #Get system-state (current coordinates and energy) in each step
    current_state=simulation.context.getState(getPositions=True, getEnergy=True)
    current_coords = np.array(current_state.getPositions(asNumpy=True))*10 # coords in Å

    #Create gradient using some user-function
    get_gradient_from_coordinates = np.random.random #Here just generating a random array
    full_system_gradient = get_gradient_from_coordinates((frag.numatoms,3))

    #Update the customforce using full_system_gradient, should be a Numpy array in Eh/Bohr
    omm.update_custom_external_force(customforce,full_system_gradient,simulation)

    #Take a simulation step
    simulation.step(1)



**Adding a flat-bottom centering potential**

For simulations involving e.g. a small molecule in a solvation box (also a ligand in a protein-ligand complex), 
it can sometimes be desirable to restrain the motion of the solute/ligand in order to prevent it from either 
exploring conformational space far away from the site-of-interest or to prevent it from going too far away from the center of the box,
the latter creating problems for QM/MM simulations where the QM-code is not periodic.

We can use the **add_centerforce** option to the **MolecularDynamics** function, to add a restraining force.
This option adds a harmonic restraining force acting only on selected atoms (e.g. solute) and features a "flat-bottom" shape, 
meaning that the force only contributes when the solute is at a chosen distance from a center of interest. 
The restraining force thus has the effect of giving a gentle or hard (controlled by the force-constant) 
kick to the solute when it travels too far from the center.

The **add_centerforce** option is easy to use. 
By default it restrains the molecule/group (defined by *centerforce_atoms* keyword), by a forceconstant of 1.0 kcal/mol/Å^2, 
when the molecule/group is more than 10 Å away from the center (defined as the geometric center of the system).
Do note that the molecule/group should usually already be present in the center in this case (or the force starts immediately acting on it).

.. code-block:: python

  # Here assuming a qm_mm theory, a solute+solvent system have already been created and that solute has indices 0-5.
  # Force constant : 100 kcal/mol/Å^2 Distance from where force acts: 10.0 Å
  MolecularDynamics(theory=qm_mm, fragment=solution, simulation_time=50, timestep=0.001, traj_frequency=100,
    add_centerforce=True, centerforce_atoms=[0,1,2,3,4,5], centerforce_constant=100, centerforce_distance=10.0)

In certain cases it can also be desirable to restrain the molecule to a specific XYZ position within the box:

.. code-block:: python

  #Here restraining the molecule to be at position [11.0,11.0,11.0] in the box.
  MolecularDynamics(theory=qm_mm, fragment=solution, simulation_time=50, timestep=0.001, traj_frequency=100,
    add_centerforce=True, centerforce_atoms=[0,1,2,3,4,5], centerforce_constant=1, centerforce_distance=10.0)
    centerforce_center=[11.0,11.0,11.0],


Choosing the value of the *centerforce_constant* as well as the *centerforce_distance*, may require some experimentation.


**Adding any custom OpenMM force to an OpenMMTheory object**

If one defines some kind of special custom force, either something built-in from the OpenMM library (i.e. defined by the OpenMM Python API instead of ASH), 
a force from `OpenMMTools Forces <https://openmmtools.readthedocs.io/en/stable/forces.html>`_ library or perhaps a force from an OpenMM plugin,
you can add it to the OpenMMTheory object and it will be included in the MD simulation.

Useful reading:

- `OpenMM Standard Forces <http://docs.openmm.org/development/userguide/theory/02_standard_forces.html>`_
- `OpenMM Custom Forces <http://docs.openmm.org/development/userguide/theory/03_custom_forces.html>`_
- `OpenMM Example Plugin <https://github.com/openmm/openmmexampleplugin>`_ 


.. code-block:: python

  from ash import *
  import openmmtools # needs to be installed: mamba install openmmtools

  frag = Fragment(databasefile="h2o.xyz")
  # Defining OpenMM object
  omm = OpenMMTheory(fragment=frag, dummysystem=True)
  # Defining a force directly. Here a reaction-field from OpenMMTools
  customforce = openmmtools.forces.UnshiftedReactionFieldForce(
          reaction_field_dielectric=78.3)

  # Add the force to OpenMMTheory object. This adds a CustomNonbondedForce to omm.system
  omm.add_force(customforce)
  # Get defined forces in OpenMMTheory system object
  print("omm system forces:", omm.system.getForces())
  #Remove force by index (if we know it)
  #omm.remove_force(1)
  #Remove force by name (if we know it)
  omm.remove_force_by_name("CustomNonbondedForce")




######################################################
mdtraj interface
######################################################

Postprocessing of MD trajectories is often necessary.
As this can be a computationally intensive task, ASH contains an interface to
the `MDtraj <https://www.mdtraj.org>`_ library that is capable of various trajectory analysis.  Requires installation of mdtraj: pip install mdtraj.

*Analyze internal coordinates of trajectory*

Often one wants to inspect how a distance, angle or torsion varies during an MD trajectory.
This can be conveniently done using the function **MDtraj_coord_analyze**

.. code-block:: python

  def MDtraj_coord_analyze(trajectory, pdbtopology=None, periodic=True, indices=None):

You need to simply provide the trajectory, pdbtopology and give the atom indices as a list.
If you provide 2 atom indices the function will grab a distance, if 3 then an angle, if 4 then a dihedral angle.

Example:

.. code-block:: python

  #Get distance between atoms 50 and 67
  MDtraj_coord_analyze("trajectory.dcd", pdbtopology="trajectory.pdb", indices=[50,67])
  #Get angle between atoms 4,7 and 10
  MDtraj_coord_analyze("trajectory.dcd", pdbtopology="trajectory.pdb", indices=[4,7,10])
    #Get dihedral angle between atoms 10,11,12,13
  MDtraj_coord_analyze("trajectory.dcd", pdbtopology="trajectory.pdb", indices=[10,11,12,13])

*Slice trajectory*

To obtain specific frames from a trajectory you use the **MDtraj_slice** function.
This will create a sliced trajectory file (format can be 'PDB' or 'DCD') containing only those frames.

.. code-block:: python

  #This will grab all frames between steps 50 and 100
  MDtraj_slice(trajectory, pdbtopology, format='PDB', frames=[50,100])

*Re-imaging trajectory*

Periodic MD trajectories from OpenMM sometimes contain the molecule split between periodic images rather than showing a whole molecule in each periodic box.
This is just a visualization artifact but to obtain a more pleasing visualization of the trajectory you can "reimage" the trajectory as shown below.

Example:

.. code-block:: python

  from ash import *
  #Provide trajectory file, PDB topology file and final format of trajectory
  MDtraj_imagetraj("output_traj.dcd", "final_MDfrag_laststep.pdb", format='DCD')
  
  #If periodic box info is missing from trajectory file (can happen with CHARMM files):
  MDtraj_imagetraj("output_traj.dcd", "final_MDfrag_laststep.pdb", format='DCD', 
    unitcell_lengths=[100.0,100.0,100.0], unitcell_angles=[90.0,90.0,90.0])


*Calculating root-mean-square fluctations (RMSF) in trajectory*

Calculates the RMSF and prints out the atoms that move the most during the trajectory.
See `mdtraj RMSF page <https://mdtraj.org/1.9.4/api/generated/mdtraj.rmsf.html>`_  for more info.

.. code-block:: python
  
  indices = MDtraj_RMSF(trajectory, pdbtopology, print_largest_values=True, threshold=0.005, largest_values=10)
  #Returns atom indices with the largest RMSF values

######################################################
Dynamics via ASE (not recommended)
######################################################

The ASE approach uses dynamics routines of the ASE library. 
This also allows molecular dynamics to be performed via any theory level in ASH: QM, MM or QM/MM theory. 
This requires the `ASE <https://wiki.fysik.dtu.dk/ase/>`_  library to be installed (simple Python pip installation).
Generally, this is an old option that is no longer recommended as the OpenMM interface is faster and more versatile.

.. code-block:: python

	def Dynamics_ASE(fragment=None, PBC=False, theory=None, temperature=300, timestep=None,
          simulation_steps=None, simulation_time=None, thermostat=None, ttime_nosehoover=5,
          barostat=None, trajectoryname="Trajectory_ASE", traj_frequency=1, coupling_freq=0.002,
          frozen_atoms=None, frozen_bonds=None,frozen_angles=None, frozen_dihedrals=None,
          plumed_object=None, multiple_walkers=False, numwalkers=None,
          gpaw_rattle_constraints=False, charge=None, mult=None,
          safires=False, safires_solute=None, safires_inner_region=None, safires_solvent_atomsnum=3):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``fragment``
     - ASH Fragment
     - None
     - ASH Fragment object.
   * - ``theory``
     - ASHTheory object
     - None
     - ASH Theory object.
   * - ``PBC``
     - Boolean
     - False
     - Whether periodic boundary conditions are used or not during simulation. UNTESTED.
   * - ``temperature``
     - integer
     - 300
     - The chosen temperature to start simulation with and maintain (if using a thermostat).
   * - ``thermostat``
     - string
     - None
     - The thermostat to use. Options: 'Langevin', 'Andersen', 'NoseHoover', 'Berendsen'.
   * - ``timestep``
     - float
     - None
     - The timestep in ps.
   * - ``simulation_steps``
     - integer
     - None
     - The number of simulation steps to carry out.
   * - ``simulation_time``
     - float
     - None
     - Alternative to simulation_steps: the simulation time.
   * - ``barostat``
     - string
     - None
     - Name of barostat. CURRENTLY INACTIVE OPTION.
   * - ``trajectoryname``
     - string
     - Trajectory_ASE
     - Name of trajectoryfile created by ASE.
   * - ``traj_frequency``
     - integer
     - 1
     - Interval between simulation trajectory updates.
   * - ``coupling_freq``
     - float
     - 0.002
     - | coupling_freq determines friction-coefficient if thermostat is Langevin, but 
       | collision probability if thermostat is Andersen.
   * - ``frozen_atoms``
     - list of integers
     - None
     - Freeze the atoms defined according to the list of atom indices.
   * - ``frozen_bonds``
     - list of integers
     - None
     - Freeze the bond/distance defined by the list of atom indices.
   * - ``frozen_angles``
     - list of integers
     - None
     - Freeze the angle defined by the list of atom indices.
   * - ``frozen_dihedrals``
     - list of integers
     - None
     - Freeze the dihedral defined by the list of atom indices.
   * - ``plumed_object``
     - <ASH Plumed object>
     - None
     - Optional ASH Plumed object for running metadynamics via the Plumed library. See later.
   * - ``multiple_walkers``
     - Boolean
     - False
     - Whether to use multiple_walkers during Plumed metadynamics simulation.
   * - ``numwalkers``
     - integer
     - None
     - How manu walkers to use with multiple_walkers option.
   * - ``ttime_nosehoover``
     - float
     - 5
     - Coupling time in fs for Nose-Hoover thermostat. CURRENTLY INACTIVE:
   * - ``safires``
     - X
     - None
     - Whether to turn SAFIRES option on.
   * - ``safires_solute``
     - list of integers
     - None
     - SAFIRES: The definition of the SAFIRES solute region.
   * - ``safires_inner_region``
     - list of integers
     - None
     - SAFIRES: The definition of the SAFIRES inner region.
   * - ``safires_solvent_atomsnum``
     - integer
     - 3
     - SAFIRES: How many atoms are in solvent molecule.
   * - ``gpaw_rattle_constraints``
     - Boolean
     - False
     - | Whether to define constraints via GPAW instead of ASE (faster).
       | Requires GPAW installation.
   * - ``charge``
     - integer
     - None
     - Optional charge. Will override charge attribute of ASH Fragment.
   * - ``mult``
     - integer
     - None
     - Optional spin multiplicity. Will override mult attribute of ASH Fragment.
  


In order to use the **Dynamics_ASE** function, `ASE <https://wiki.fysik.dtu.dk/ase/>`_ must have been installed before to the same Python environment 
that ASH uses (easiest done via: pip install ase).

**Examples:**

*Simple NVE example:*

.. code-block:: python

	from ash import *
	
	numcores=12
	#Simple n-butane system
	butane=Fragment(xyzfile="butane.xyz", charge=0, mult=1)

	# Creating xTBTheory object (Note: runmode='library' runs faster) that is parallelized. Using GFN1-xTB.
	xtbcalc = xTBTheory(xtbmethod='GFN1', runmode='library', numcores=numcores)
	
	#Running NVE dynamics (initial temp=300 K) on butane using xTBTheory.
	# 0.001 ps timestep, 100000 steps, writing every 10th step to trajectory. A velocity Verlet algorithm is used.
	Dynamics_ASE(fragment=butane, theory=xtbcalc, temperature=300.0, timestep=0.001, simulation_steps=100000, traj_frequency=10)


*Simple NVT (Langevin thermostat) example:*

.. code-block:: python

	#Running NVT dynamics with a Langevin thermostat on butane using xTBTheory
	# 0.001 ps timestep, 100000 steps, writing every 10th step to trajectory.
	Dynamics_ASE(fragment=butane, theory=xtbcalc, thermostat='Langevin', coupling_freq=0.002, 
		temperature=300.0, timestep=0.001, simulation_steps=100000, traj_frequency=10)

