=================================================
Molecular dynamics
=================================================

Molecular dynamics in ASH is currently available via 2 different approaches: 1) dynamics routines via the OpenMM library or 2) the ASE library.
Both approaches support all available ASH Theory levels.

- The OpenMM approach uses an interface to the MM dynamics routines of OpenMM library. It is particularly recommended for running MM and QM/MM simulations when system size is larger. This is by far the best option when a significant part of time is spent on calculating the MM energy+gradient. This approach also allows QM dynamics. It requires the OpenMM library to be installed.
- The ASE approach uses dynamics routines of the ASE library. This also allows molecular dynamics to be performed via any theory level in ASH: QM, MM or QM/MM theory. This requires the `ASE <https://wiki.fysik.dtu.dk/ase/>`_  library to be installed (simple Python pip installation). 

In the future, ASH may feature it's own native dynamics.



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
Metadynamics via OpenMM_MD and Plumed
######################################################

to be documented



######################################################
Dynamics via ASE
######################################################

The Dynamics_ASE function allows NVE and NVT based molecular dynamics in ASH using any available QM, MM or QM/MM theory.

.. code-block:: python

	def Dynamics_ASE(fragment=None, PBC=False, theory=None, temperature=300, timestep=None, thermostat=None, simulation_steps=None, simulation_time=None,
					barostat=None, trajectoryname="Trajectory_ASE", traj_frequency=1, coupling_freq=0.002, frozen_atoms=None, frozen_bonds=None,
					frozen_angles=None, frozen_dihedrals=None, plumed_object=None, multiple_walkers=False, numwalkers=None,
					ttime_nosehoover=5, safires=False, safires_solute=None, safires_inner_region=None, safires_solvent_atomsnum=3,
					gpaw_rattle_constraints=False, charge=None, mult=None):

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
     - coupling_freq determines friction-coefficient if thermostat is Langevin, but collision probability if thermostat is Andersen.
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
     - Whether to define constraints via GPAW instead of ASE (faster). Requires GPAW installation.
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



######################################################
Metadynamics via ASE and Plumed
######################################################


Via an interface to the `Plumed <https://www.plumed.org>`_ biased sampling library it is possible to perform metadynamics in ASH with
the help of the ASE dynamics routines. Any theory level in ASH is supported (including QM/MM theories).

About the ASH-Plumed interface:

- Only tested for metadynamics. Only 1D and 2D metadynamics currently possible.
- Well-tempered metadynamics is always specified in the current interface (regular metadynamics is a largely obsolete method).
- ASH uses the same units for distance (Å), energy (eV) and time (ps) as the dynamics program (currently ASE). Radians are used for torsions. This is different from the default Plumed units (nm for distances and kJ/mol for energy). Keep this in mind when defining sigma (width of Gaussian in CV-unit) and height (of Gaussian in energy-unit).

.. note:: Not yet available: multiple-walker metadynamics

**Requirements:**

- `ASE <https://wiki.fysik.dtu.dk/ase/>`_ library (see above)
- `Plumed <https://www.plumed.org>`_ installation (requires compilation). Alternatively it might be possible to install via `conda-forge <https://anaconda.org/conda-forge/plumed>`_ (untested)
- Plumed Python wrappers (pip install plumed)


.. code-block:: python

  class plumed_ASH():
      def __init__(self, path_to_plumed_kernel=None, bias_type="MTD", fragment=None, CV1_type=None, CV1_indices=None,
                  CV2_type=None, CV2_indices=None, temperature=300.0, hills_file="HILLS", colvar_file="COLVAR", 
                  height=0.01243, sigma1=None, sigma2=None, biasfactor=6.0, timestep=None,
                  stride_num=10, pace_num=500, dynamics_program=None, numwalkers=None, debug=False):

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
   * - ``path_to_plumed_kernel``
     - string
     - None
     - Should give full path to the libplumedKernel.so file in Plumed installation.
   * - ``bias_type``
     - string
     - 'MTD
     - Current options: "MTD" (for metadynamics job)   (more to come...)
   * - ``CV1_type``
     - string
     - None
     - Type of collective variable 1 (Plumed keyword). Options: TORSION, DISTANCE, ANGLE, RMSD (and more in principle)
   * - ``CV1_indices``
     - list of integers
     - None
     - List of atom indices that defines the chosen torsion, distance, angle (note: use 0-based indexing)
   * - ``CV2_type``
     - string
     - None
     - Type of collective variable 2 (Plumed keyword). Options: TORSION, DISTANCE, ANGLE, RMSD (and more in principle)
   * - ``CV2_indices``
     - list of integers
     - None
     - List of CV2 atom indices that defines the chosen torsion, distance, angle (note: use 0-based indexing)
   * - ``temperature``
     - float
     - 300.0
     - The temperature provided to Plumed (in Kelvin). Used in well-tempered MTD.
   * - ``hills_file``
     - string
     - 'HILLS'
     - Name of HILLS-file (default HILLS).
   * - ``colvar_file``
     - string
     - 'COLVAR'
     - Name of COLVAR-file (default COLVAR). 
   * - ``height``
     - float
     - 0.01243
     -  The height of the Gaussian in energy-unit eV. Default: 0.01243 eV (= 1.2 kJ/mol)
   * - ``sigma1/sigma2``
     - X
     - None
     - he width of the Gaussian in CV units for each CV defined. Depends on CV-type. Example: sigma1=0.35 radians(torsion), sigma1=0.5 Å (distance).
   * - ``biasfactor``
     - float
     - 6.0
     - Parameter used in well-tempered metadynamics. Default: 6.0
   * - ``timestep``
     - float
     - None
     - The timestep (in ps) provided to Plumed.
   * - ``stride_num``
     - integer
     - 10
     - Frequency of writing to COLVAR file.
   * - ``pace_num``
     - integer
     - 500
     -  Frequency of writing to HILLS file.
   * - ``numwalkers``
     - integer
     - None
     - Number of walkers used for multiple walker metadynamics. CURRENTLY INACTIVE
   * - ``dynamics_program``
     - string
     - None
     - Name of dynamics program used.
   * - ``debug``
     - Boolean
     - False
     - Debug mode.




Function to analyze the results of the metadynamics

.. code-block:: python

  def MTD_analyze(plumed_ash_object=None, path_to_plumed=None, Plot_To_Screen=False, CV1_type=None, CV2_type=None, temperature=None,
                  CV1_indices=None, CV2_indices=None):


.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``plumed_ash_object``
     - plumed_ASH
     - None
     - An object of class plumed_ASH.
   * - ``path_to_plumed``
     - string
     - None
     - Path to Plumed directory (containing lib dir etc.)
   * - ``Plot_To_Screen``
     - Boolean
     - False
     - Whether to plot graph to screen or not.
   * - ``CV1_type``
     - string
     - None
     - Type of CV1.
   * - ``CV2_type``
     - string
     - None
     - Type of CV1.
   * - ``temperature``
     - float
     - None
     - Temperature in Kelvin
   * - ``CV1_indices``
     - list of integers
     - None
     - List of integers defining CV1
   * - ``CV2_indices``
     - list of integers
     - None
     - List of integers defining CV2


**Examples:**

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
		CV1_type="Torsion", temperature=298.15, CV1atoms=[0,3,7,10])


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
		CV1_type="Torsion", temperature=298.15, CV1atoms=[0,3,7,10])


