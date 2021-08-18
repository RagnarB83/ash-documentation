=================================================
Dynamics module
=================================================

Molecular dynamics in ASH is currently available via an interface to the dynamics routines of the ASE library. This requires the `ASE <https://wiki.fysik.dtu.dk/ase/>`_  library to be installed (simple Python pip installation). In the future, ASH will feature it's own native dynamics.


.. note:: For pure classical forcefield-based MD it is recommended to instead run the dynamics via OpenMM (faster). This requires the system to have been set up using OpenMMTheory and utilizes the OpenMM library for energy, forces and dynamics. See :doc:`OpenMM-interface`


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
	butane=Fragment(xyzfile="butane.xyz")

	# Creating xTBTheory object (Note: runmode='library' runs faster) that is parallelized. Using GFN1-xTB.
	xtbcalc = xTBTheory(charge=0, mult=1, xtbmethod='GFN1', runmode='library', numcores=numcores)
	
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
Metadynamics via ASE and Plumed
######################################################

Via an interface to the `Plumed <https://www.plumed.org>`_ library it is possible to perform metadynamics in ASH. Any theory level in ASH is supported (including QM/MM theories).

Requirements:

- `ASE <https://wiki.fysik.dtu.dk/ase/>`_ library (see above)
- `Plumed <https://www.plumed.org>`_ installation (requires compilation). Alternatively it might be possible to install via `conda-forge <https://anaconda.org/conda-forge/plumed>`_ (untested)
- Plumed Python wrappers (pip install plumed)


.. code-block:: python

	from ash import *

	numcores=12

	#Simple n-butane system
	frag=Fragment(xyzfile="butane.xyz")
	# Creating xTBTheory object (Note: runmode='library' runs faster) that is parallelized 
	xtbcalc = xTBTheory(charge=0, mult=1, xtbmethod='GFN1', runmode='library', numcores=numcores)

	#Create ASH-Plumed object. Points to Plumed kernel and defines collective variables etc.
	plumed_object = plumed_ASH(path_to_plumed_kernel="/home/bjornsson/plumed-install-serial/lib/libplumedKernel.so", 
					bias_type="1D_MTD", fragment=frag, colvar_type="torsion", colvar_indices=[0,3,7,10],
	                temperature=298.15, hills_file="HILLS", colvar_file="COLVAR", height=0.012437126761597656, 
	                sigma=0.35, biasfactor=6.0, timestep=0.001, stride_num=1, pace_num=1)

	#Call ASH-ASE dynamics with plumed_object. Here running 100K steps with 1 fs timstep, writing trajectory every 10th step.
	Dynamics_ASE(fragment=frag, theory=xtbcalc, timestep=0.001, simulation_steps=100000, traj_frequency=10, plumed_object=plumed_object)

	#Analyze the results of the metadynamics
	MTD_analyze(path_to_plumed="/home/bjornsson/plumed-install-serial", Plot_To_Screen=False, 
		colvar_type="Torsion", temperature=298.15, CV1atoms=[0,3,7,10])


.. note:: Not yet available: multiple-walker metadynamics

