=================================================
Molecular dynamics
=================================================

Molecular dynamics in ASH is currently available via 2 different approaches.

- The more general approach make use of an interface to the dynamics routines of the ASE library. This allows molecular dynamics to be performed via any theory level in ASH: QM, MM or QM/MM theory. This requires the `ASE <https://wiki.fysik.dtu.dk/ase/>`_  library to be installed (simple Python pip installation). In the future, ASH will feature it's own native dynamics.

- The more specialized approach uses an interface to the MM dynamics routines of OpenMM library. This is only available for theory levels: OpenMMTheory (i.e. pure classical simulations) and QMMMTheory (with mm_theory=OpenMMTheoryobject). This approach makes the most sense when a considerable amount of time is spent on calculating the MM energy+gradient via OpenMM (faster than the Dynamics_ASE approach)


######################################################
Dynamics via ASE
######################################################

The Dynamics_ASE function allows NVE and NVT based molecular dynamics in ASH using any available QM, MM or QM/MM theory.

.. code-block:: python

	def Dynamics_ASE(fragment=None, theory=None, temperature=300, timestep=None, thermostat=None, simulation_steps=None, 
				simulation_time=None, barostat=None, trajectoryname="Trajectory_ASE", traj_frequency=1, coupling_freq=0.002, 
				frozen_atoms=None, frozen_bonds=None, frozen_angles=None, frozen_dihedrals=None, plumed_object=None):

In order to use the Dynamics_ASE function, `ASE <https://wiki.fysik.dtu.dk/ase/>`_ must have been installed before to the same Python that ASH uses (easiest done via: pip install ase).

Simple NVE example:

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


Simple NVT (Langevin thermostat) example:

.. code-block:: python

	#Running NVT dynamics with a Langevin thermostat on butane using xTBTheory
	# 0.001 ps timestep, 100000 steps, writing every 10th step to trajectory.
	Dynamics_ASE(fragment=butane, theory=xtbcalc, thermostat='Langevin', coupling_freq=0.002, 
		temperature=300.0, timestep=0.001, simulation_steps=100000, traj_frequency=10)


Thermostat options: 'Langevin', 'Andersen', 'NoseHoover', 'Berendsen'.

coupling_freq determines friction-coefficient in Langevin, collions probability in Andersen.

It is possible to freeze atoms using frozen_atoms= option. Provide list of atom indices.

Bonds, angles and dihedrals can be frozen using frozen_bonds=, frozen_angles= and frozen_dihedrals= options.


######################################################
Dynamics via OpenMM_MD
######################################################

For pure classical forcefield-based MD it is recommended to run the dynamics via the OpenMM library instead via the ASH **OpenMM_MD** function. The reason is that this results in essentially no data transfer between the C++ layer (of OpenMM) and Python layer (of OpenMM and ASH) while this is the case if running via Dynamics_ASE. Dynamics via OpenMM_MD requires the system to have been set up using OpenMMTheory and utilizes the OpenMM library for energy, forces and dynamics. See :doc:`OpenMM-interface` for details.

**Pure MM example:**

.. code-block:: python

	from ash import *

	#Fragment
	frag=Fragment(xyzfile="frag.xyz")
	#Defining frozen region. Taking difference of all-atom region and active region
	actatoms=[14,15,16]
	frozen_atoms=listdiff(frag.allatoms,actatoms)

	openmmobject = OpenMMTheory(cluster_fragment=frag, ASH_FF_file="Cluster_forcefield.ff", frozen_atoms=frozen_atoms)

	OpenMM_MD(fragment=frag, theory=openmmobject, timestep=0.001, simulation_time=2, traj_frequency=10, temperature=300,
	    integrator='LangevinIntegrator', coupling_frequency=1)


For a QM/MM system that utilizes OpenMMTheory as mm_theory and any QM-theory as qm_theory, it is also possible to use OpenMM_MD to do QM/MM dynamics. In this case the QM+PC gradient is used to update the forces of the OpenMM system (as a CustomExternalForce)
This is beneficial if a considerable amount of time of the QM/MM energy+gradient is spent on calculating the MM energy+gradient and then there the reduced data transfer (and unnecessary data conversion) between the Python and C++ layers results in faster MM energy+gradient steps. This is only the case if the QM-theory is really cheap (i.e. a semi-empirical method like xTB or AM1, PM3), otherwise the QM energy+gradient will dominate the total cost. See :doc:`OpenMM-interface` for details.

**QM/MM example:**

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



######################################################
Metadynamics via ASE and Plumed
######################################################

Via an interface to the `Plumed <https://www.plumed.org>`_ library it is possible to perform metadynamics in ASH. Any theory level in ASH is supported (including QM/MM theories).

Requirements:

- `ASE <https://wiki.fysik.dtu.dk/ase/>`_ library (see above)
- `Plumed <https://www.plumed.org>`_ installation (requires compilation). Alternatively it might be possible to install via `conda-forge <https://anaconda.org/conda-forge/plumed>`_ (untested)
- Plumed Python wrappers (pip install plumed)


1D metadynamics example (torsion):

.. code-block:: python

	from ash import *

	numcores=12

	#Simple n-butane system
	frag=Fragment(xyzfile="butane.xyz", charge=0, mult=1)
	# Creating xTBTheory object (Note: runmode='library' runs faster) that is parallelized 
	xtbcalc = xTBTheory(xtbmethod='GFN1', runmode='library', numcores=numcores)

	#Create ASH-Plumed object. Points to Plumed kernel and defines collective variables etc.
	plumed_object = plumed_ASH(path_to_plumed_kernel="/home/bjornsson/plumed-install-serial/lib/libplumedKernel.so", 
					bias_type="MTD", fragment=frag, CV1_type="TORSION", CV1_indices=[0,3,7,10],
	                temperature=298.15, hills_file="HILLS", colvar_file="COLVAR", height=0.012, 
	                sigma1=0.35, biasfactor=6.0, timestep=0.001, stride_num=1, pace_num=1)

	#Call ASH-ASE dynamics with plumed_object. Here running 100K steps with 1 fs timstep, writing trajectory every 10th step.
	Dynamics_ASE(fragment=frag, theory=xtbcalc, timestep=0.001, simulation_steps=100000, traj_frequency=10, plumed_object=plumed_object)

	#Analyze the results of the metadynamics
	MTD_analyze(path_to_plumed="/home/bjornsson/plumed-install-serial", Plot_To_Screen=False, 
		colvar_type="Torsion", temperature=298.15, CV1atoms=[0,3,7,10])


2D metadynamics example (torsion,distance):

.. code-block:: python

	from ash import *

	numcores=12

	#Simple n-butane system
	frag=Fragment(xyzfile="butane.xyz", charge=0, mult=1)
	# Creating xTBTheory object (Note: runmode='library' runs faster) that is parallelized 
	xtbcalc = xTBTheory(xtbmethod='GFN1', runmode='library', numcores=numcores)

	#Create ASH-Plumed object. Points to Plumed kernel and defines collective variables etc.
	plumed_object = plumed_ASH(path_to_plumed_kernel="/home/bjornsson/plumed-install-serial/lib/libplumedKernel.so", 
					bias_type="MTD", fragment=frag, CV1_type="TORSION", CV1_indices=[0,3,7,10], CV2_type="DISTANCE", CV2_indices=[1,2],
	                temperature=298.15, hills_file="HILLS", colvar_file="COLVAR", height=0.012, 
	                sigma1=0.35, sigma2=0.5, biasfactor=6.0, timestep=0.001, stride_num=1, pace_num=1)

	#Call ASH-ASE dynamics with plumed_object. Here running 100K steps with 1 fs timstep, writing trajectory every 10th step.
	Dynamics_ASE(fragment=frag, theory=xtbcalc, timestep=0.001, simulation_steps=100000, traj_frequency=10, plumed_object=plumed_object)

	#Analyze the results of the metadynamics
	MTD_analyze(path_to_plumed="/home/bjornsson/plumed-install-serial", Plot_To_Screen=False, 
		colvar_type="Torsion", temperature=298.15, CV1atoms=[0,3,7,10])


ASH Plumed class keywords:

- path_to_plumed_kernel (string). Should give full path to the libplumedKernel.so file in Plumed installation.
- bias_type (string). Current options: "MTD" (for metadynamics job)   (more to come...)
- fragment (ASH fragment). The ASH fragment for the system.
- CV1_type/CV2_type (string). Type of collective variable 1 (Plumed keyword). Options: TORSION, DISTANCE, ANGLE, RMSD (and more in principle)
- CV1_indices/CV2_indices (list). List of atom indices that defines the chosen torsion, distance, angle (note: use 0-based indexing)
- temperature (float). The temperature provided to Plumed (in Kelvin). Used in well-tempered MTD
- hills_file (string). Name of HILLS-file (default HILLS).
- colvar_file (string). Name of COLVAR-file (default COLVAR). 
- height (float). The height of the Gaussian in energy-unit eV. Default: 0.01243 eV (= 1.2 kJ/mol)
- sigma1/sigma2 (float). The width of the Gaussian in CV units for each CV defined. Depends on CV-type. Example: sigma1=0.35 radians(torsion), sigma=0.5 Å (distance).
- biasfactor (float). Parameter used in well-tempered metadynamics. Default: 6.0
- timestep (float). The timestep (in ps) provided to Plumed.
- stride_num (int). Frequency of writing to COLVAR file. Default: 10
- pace_num (int). Frequency of writing to HILLS file. Default: 500
- numwalkers (int). Number of walkers used for multiple walker metadynamics. CURRENTLY INACTIVE


About the ASH-Plumed interface:

- Well-tempered metadynamics is always specified in the current interface (regular metadynamics is a largely obsolete method).
- Only 1D and 2D metadynamics currently possible.
- ASH uses the same units for distance (Å), energy (eV) and time (ps) as the dynamics program (currently ASE). Radians are used for torsions. This is different from the default Plumed units (nm for distances and kJ/mol for energy). Keep this in mind when defining sigma (width of Gaussian in CV-unit) and height (of Gaussian in energy-unit).



.. note:: Not yet available: multiple-walker metadynamics


######################################################
Metadynamics via OpenMM_MD and Plumed
######################################################

to be documented

