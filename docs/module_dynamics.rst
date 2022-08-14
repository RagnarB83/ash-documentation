
Molecular dynamics
======================================

Molecular dynamics in ASH is currently available via 2 different approaches: 1) dynamics routines via the OpenMM library or 2) the ASE library.
Both approaches support all available ASH Theory levels.

- The OpenMM approach uses an interface to the MM dynamics routines of OpenMM library. It is particularly recommended for running MM and QM/MM simulations when system size is larger. This is by far the best option when a significant part of time is spent on calculating the MM energy+gradient. This approach also allows QM dynamics. It requires the OpenMM library to be installed.
- The ASE approach uses dynamics routines of the ASE library. This also allows molecular dynamics to be performed via any theory level in ASH: QM, MM or QM/MM theory. This requires the `ASE <https://wiki.fysik.dtu.dk/ase/>`_  library to be installed (simple Python pip installation). 

In the future, ASH may feature it's own native dynamics.

Support for biased sampling e.g. metadynamics is also available. See :doc:`Biased-sampling`.

######################################################
Dynamics via OpenMM_MD
######################################################

Dynamics via **OpenMM_MD** is becoming the recommended way of running molecular dynamics in ASH via any Hamiltonian: QM, MM or QM/MM.
See :doc:`OpenMM-interface` for details on the **OpenMM_MD** function.

**Pure MM example:**

For pure classical forcefield-based MD it is strongly recommended to run the dynamics via the OpenMM library via the ASH **OpenMM_MD** function. 
The reason is that this results in essentially no data transfer between the C++ layer (of OpenMM) and Python layer (of OpenMM and ASH) while this is the case if running via **Dynamics_ASE**. 
Running OpenMM via the GPU code (if a GPU is available) will allow particularly fast MM dynamics.

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
	OpenMM_MD(fragment=frag, theory=openmmobject, timestep=0.001, simulation_time=2, traj_frequency=10, temperature=300,
	    integrator='LangevinIntegrator', coupling_frequency=1)



**QM/MM example:**

For a QM/MM system that utilizes OpenMMTheory as mm_theory and any QM-theory as qm_theory, it is also possible to use **OpenMM_MD** to do QM/MM dynamics. In this case the QM+PC gradient is used to update the forces of the OpenMM system (as a CustomExternalForce).
This is beneficial if a considerable amount of time of the QM/MM energy+gradient is spent on calculating the MM energy+gradient and then there the reduced data transfer (and unnecessary data conversion) between the Python and C++ layers results in faster MM energy+gradient steps. This is only the case if the QM-theory is really cheap (i.e. a semi-empirical method like xTB or AM1, PM3), otherwise the QM energy+gradient will dominate the total cost. See :doc:`OpenMM-interface` for details.

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

	OpenMM_MD(fragment=frag, theory=QMMMTheory, timestep=0.001, simulation_time=2, traj_frequency=10, temperature=300,
	    integrator='LangevinIntegrator', coupling_frequency=1, charge=0, mult=1)

**QM example:**

It is even possible to use the dynamics routines of the OpenMM library to drive an MD simulation at the QM-level. This is possible by setting up a dummy MM system and reading in the QM-theory forces (via ASH) as a custom external force to the OpenMM theory.


.. code-block:: python

	from ash import *
	
	numcores=12
	#Simple n-butane system
	butane=Fragment(xyzfile="butane.xyz", charge=0, mult=1)

	# Creating xTBTheory object (Note: runmode='library' runs faster) that is parallelized. Using GFN1-xTB.
	xtbcalc = xTBTheory(xtbmethod='GFN1', runmode='library', numcores=numcores)
	
	#Running NVE dynamics (initial temp=300 K) on butane using xTBTheory.
	# 0.001 ps timestep, 2 ps , writing every 10th step to trajectory. A velocity Verlet algorithm is used.
	OpenMM_MD(fragment=butane, theory=xtbcalc, timestep=0.001, simulation_time=2, traj_frequency=10, temperature=300,
		integrator='LangevinIntegrator', coupling_frequency=1)



######################################################
Dynamics via ASE
######################################################

The Dynamics_ASE function allows NVE and NVT based molecular dynamics in ASH using any available QM, MM or QM/MM theory.

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

