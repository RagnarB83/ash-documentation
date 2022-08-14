Biased sampling MD
======================================

EXPERIMENTAL: THIS IS NOT YET READY UNFORTUNATELY

Biased sampling e.g. metadynamics is possible via an interface to the `Plumed <https://www.plumed.org>`_ biased sampling library.

Such simulations are performed using the dynamics functions documend in :doc:`module_dynamics`
and it works by providing an object of class plumed_ASH to functions: **Dynamics_ASE** or **OpenMM_MD** using the plumed_object keyword.
Any theory level in ASH can be used(including QM/MM theories).

**Current status:**

- Technically working for Dynamics_ASE but requires testing and confirm units.
- Available for OpenMM_MD but untested and probably not working.
- Only tested for metadynamics. Only 1D and 2D metadynamics currently possible.
- Well-tempered metadynamics is always specified in the current interface (regular metadynamics is a largely obsolete method).
- ASH uses the same units for distance (Å), energy (eV) and time (ps) as the dynamics program (currently ASE). Radians are used for torsions. This is different from the default Plumed units (nm for distances and kJ/mol for energy). Keep this in mind when defining sigma (width of Gaussian in CV-unit) and height (of Gaussian in energy-unit).

.. note:: Not yet available: multiple-walker metadynamics

**Requirements:**

- `Plumed <https://www.plumed.org>`_ installation (requires compilation). Alternatively it might be possible to install via `conda-forge <https://anaconda.org/conda-forge/plumed>`_ (untested)
- Plumed Python wrappers (pip install plumed)


######################################################
Plumed_ASH class
######################################################

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
     - | Type of collective variable 1 (Plumed keyword). 
       | Options: TORSION, DISTANCE, ANGLE, RMSD (and more in principle)
   * - ``CV1_indices``
     - list of integers
     - None
     - List of atom indices that defines the chosen torsion, distance, angle.
   * - ``CV2_type``
     - string
     - None
     - | Type of collective variable 2 (Plumed keyword). 
       | Options: TORSION, DISTANCE, ANGLE, RMSD (and more in principle)
   * - ``CV2_indices``
     - list of integers
     - None
     - List of CV2 atom indices that defines the chosen torsion, distance, angle.
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
     - float
     - None
     - | The width of the Gaussian in CV units for each CV defined. Depends on CV-type. 
       | Example: sigma1=0.35 radians(torsion), sigma1=0.5 Å (distance).
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


######################################################
MTD_analyze: Analyze the results
######################################################

MTD_analyze is a function to analyze the results of the metadynamics simulation.

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
     - Temperature in Kelvin.
   * - ``CV1_indices``
     - list of integers
     - None
     - List of integers defining CV1.
   * - ``CV2_indices``
     - list of integers
     - None
     - List of integers defining CV2.


######################################################
Examples:
######################################################


*1D metadynamics example (torsion):*

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


*2D metadynamics example (torsion,distance):*

.. code-block:: python

	from ash import *

	numcores=12

	#Simple n-butane system
	frag=Fragment(xyzfile="butane.xyz", charge=0, mult=1)
	# Creating xTBTheory object (Note: runmode='library' runs faster) that is parallelized 
	xtbcalc = xTBTheory(xtbmethod='GFN1', runmode='library', numcores=numcores)

	#Create ASH-Plumed object. Points to Plumed kernel and defines collective variables etc.
	plumed_object = plumed_ASH(path_to_plumed_kernel="/home/bjornsson/plumed-install-serial/lib/libplumedKernel.so", 
          bias_type="MTD", fragment=frag, hills_file="HILLS", colvar_file="COLVAR", 
          CV1_type="TORSION", CV1_indices=[0,3,7,10], CV2_type="DISTANCE", CV2_indices=[1,2],
          sigma1=0.35, sigma2=0.5, biasfactor=6.0, temperature=298.15, height=0.012,
          timestep=0.001, stride_num=1, pace_num=1)

	#Call ASH-ASE dynamics with plumed_object. Here running 100K steps with 1 fs timstep, writing trajectory every 10th step.
	Dynamics_ASE(fragment=frag, theory=xtbcalc, timestep=0.001, simulation_steps=100000, traj_frequency=10, plumed_object=plumed_object)

	#Analyze the results of the metadynamics
	MTD_analyze(path_to_plumed="/home/bjornsson/plumed-install-serial", Plot_To_Screen=False, 
		CV1_type="Torsion", temperature=298.15, CV1atoms=[0,3,7,10])


