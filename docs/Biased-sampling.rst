Biased sampling MD & Free energy simulations
===============================================

Biased sampling or free-energy simulations are possible in ASH via the OpenMM molecular dynamics routines (see :doc:`module_dynamics` and :doc:`OpenMM-interface`).
Umbrella sampling could in principle be performed by adding restraint potentials to the OpenMMTheory object before running MD. Convenient workflows in ASH for umbrella sampling are currently missing, however.

Metadynamics has become a very popular biased sampling / free energy simulation method due to its ease-of-use and handy parallelization strategy.
It is possible to perform metadynamics simulations in ASH using essentially any theory-level, including an OpenMMTheory level, any type of QMTheory level and a QMMMTheory level.
The OpenMM routines are used for the simulations, and in the case of QM and QM/MM Theories, the energy and forces are passed onto OpenMM in each timestep.

The **OpenMM_metadynamics** function is used to start a metadynamics simulation from an ASH Fragment and ASH Theory level and some collective variable information.
The function sets up the necessary collective-variable bias potentials before launching an MD simulation using the **OpenMM_MD** class.
Two different metadynamics modes are available using this function:

1. OpenMM native metadynamics (uses the `OpenMM implementation <http://docs.openmm.org/development/api-python/generated/openmm.app.metadynamics.Metadynamics.html>`_)
2. PLUMED-OpenMM metadynamics (requires installation of the `OpenMM PLUMED plugin <https://github.com/openmm/openmm-plumed>`_ and `Plumed <https://www.plumed.org>`_ , the biased sampling library.

The first option is generally recommended as it is faster (no OpenMM-Plumed communication per timestep required) and the most useful collective variable options are available.
The PLUMED-OpenMM option may offer more flexiblity with respect to collective variable options but has been tested less.

######################################################
OpenMM_metadynamics
######################################################

See also :doc:`mtd_tutorial` for working examples on how to perform metadynamics in ASH.


**OpenMM_metadynamics**  has all the same options as the **OpenMM_MD** function (see :doc:`module_dynamics`) but has in addition some special metadynamics keyword options.

.. code-block:: python

  def OpenMM_metadynamics(fragment=None, theory=None, timestep=0.004, simulation_steps=None, simulation_time=None,
                traj_frequency=1000, temperature=300, integrator='LangevinMiddleIntegrator',
                barostat=None, pressure=1, trajectory_file_option='DCD', trajfilename='trajectory',
                coupling_frequency=1, charge=None, mult=None, platform='CPU',
                anderson_thermostat=False, restraints=None, 
                enforcePeriodicBox=True, dummyatomrestraint=False, center_on_atoms=None, solute_indices=None,
                datafilename=None, dummy_MM=False, plumed_object=None, add_center_force=False,
                center_force_atoms=None, centerforce_constant=1.0, barostat_frequency=25, specialbox=False,
                CV1_atoms=None, CV2_atoms=None, CV1_type=None, CV2_type=None, biasfactor=6, 
                height=1, flatbottom_restraint_CV1=None, flatbottom_restraint_CV2=None,
                CV1_biaswidth=0.5, CV2_biaswidth=0.5, CV1_range=None, CV2_range=None,
                frequency=1, savefrequency=10,
                biasdir='.', use_plumed=False, plumed_input_string=None,):


.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``CV1_atoms``
     - list
     - None
     - List of atoms defining the 1st CV.
   * - ``CV2_atoms``
     - list
     - None
     - List of atoms defining the 2nd CV (if desired)
   * - ``CV1_type``
     - string
     - None
     - Name of CV-type for CV1. Options: 'bond' (or 'distance'), 'angle', 'dihedral' (or 'torsion') or 'rmsd'.
   * - ``CV2_type``
     - string
     - None
     - Name of CV-type for CV2. Options: 'bond' (or 'distance'), 'angle', 'dihedral' (or 'torsion') or 'rmsd'.
   * - ``CV1_biaswidth``
     - float
     - None
     - Biaswidth for CV1 (sigma) in CV units. CV-unit is Angstrom for 'bond' and 'rmsd' but radians for 'angle' and 'dihedral'. 
   * - ``CV2_biaswidth``
     - string
     - None
     - Biaswidth for CV2 (sigma) in CV units. CV-unit is Angstrom for 'bond' and 'rmsd' but radians for 'angle' and 'dihedral'.
   * - ``biasfactor``
     - int
     - 6
     - The biasfactor in well-tempered metadynamics. Default is 6.
   * - ``height``
     - float
     - 1
     - The height of the Gaussian hill added in kJ/mol. Default is 1 kJ/mol.
   * - ``CV1_range``
     - list
     - None
     - Optional setting of the min-max range that CV1 can take (in CV-units: Angstrom for 'bond' and 'rmsd' and radians for 'angle' and 'dihedral')
   * - ``CV2_range``
     - list
     - None
     - Optional setting of the min-max range that CV2 can take (in CV-units: Angstrom for 'bond' and 'rmsd' and radians for 'angle' and 'dihedral')
   * - ``frequency``
     - integer
     - 1
     - The interval in time steps at which Gaussians will be added to the bias potential
   * - ``savefrequency``
     - integer
     - 10
     - The interval in time steps at which to write out the current biases to disk. At the same time it writes biases, it also checks for updated biases written by other processes and loads them in. This must be a multiple of frequency.
   * - ``biasdirectory``
     - string
     - '.'
     - The name or path of the biasdirectory where biases are written (and read) during simulation. Can be a local directory or global directory.
   * - ``flatbottom_restraint_CV1``
     - list
     - None
     - List of parameters (max value in Ang unit and force constant in kcal/mol/Ang^2) for an optional flatbottom restraint (only for bond and rmsd) for CV1 that prevents the simulation from straying too far from a max value.
   * - ``flatbottom_restraint_CV2``
     - list
     - None
     - List of parameters (max value in Ang unit and force constant in kcal/mol/Ang^2) for an optional flatbottom restraint (only for bond and rmsd) for CV2 that prevents the simulation from straying too far from a max value.
   * - ``use_plumed``
     - Boolean
     - False
     - Whether to use the OpenMM-Plumed interface (requires installation of plugin.)
   * - ``plumed_input_string``
     - string
     - None
     - Optional multi-line string containing the PLUMED input syntax.

######################################################
metadynamics_plot_data
######################################################

If using the native metadynamics implementation inside OpenMM then it is convenient to use the **metadynamics_plot_data** 
function to analyze the bias-files (after or during simulation), get the free-energy surface and plot the data. 
Plotting requires a Matplotlib installation.
The function reads the metadynamics simulation parameter files from the 'ASH_MTD_parameters.txt' that should be present in the biasdirectory.

.. code-block:: python

  def metadynamics_plot_data(biasdir=None, dpi=200, imageformat='png', plot_xlim=None, plot_ylim=None ):

To use, just write a simple input-script, call metadynamics_plot_data, giving the location of the biasdirectory:

.. code-block:: python

  from ash import *
  metadynamics_plot_data(biasdir='/path/to/biasdirectory')

The function will automatically detect whether the simulation used 1 or 2 CVs and will convert data to suitable units.
The script will write the following files that can be used on their own:

- CVn_coord_values.txt # File(s) containing the CV1/CV2 values on the original grid created during the simulation setup (controlled by CV1_range keyword)
- MTD_free_energy.txt # File containing a numpy array of the free-energy per gridpoint in kcal/mol.
- MTD_free_energy_rel.txt # File containing a numpy array of the relative free-energy per gridpoint in kcal/mol.
- MTD_CV1.png/MTD_CV1_CV2_.png # Image containing the final plot (requires Matplotlib)


######################################################
Restraining CVs
######################################################

For CVs such as 'bond'/'distance' or 'rmsd' it is possible for the metadynamics simulations to wander too far off the region of interest
and furthermore if the simulation involves a reaction where a substrate/product is dissociated from a fragment this may cause problems in sampling the region of primary importance.
For dealing with such scenarios it is possible to add a flatbottom restraining potential that pushes the system away from such bad regions by a harmonic resraint.

To use, you simply add a flatbottom_restraint_CV1 or flatbottom_restraint_CV2 keyword to the OpenMM_metadynamics function, specifying the max value that the CV
can take and the force-constant of the restraint.
Example:

.. code-block:: python


  OpenMM_metadynamics(..., CV1_atoms=[0,10], CV1_type='distance', CV1_biaswidth=0.5, flatbottom_restraint_CV1=[5.0, 7.0])

This example adds a 'distance' CV between atom 0 and atom 10, with a biaswidth of 0.5 Angstrom and a restraint has been added
so that if the CV1 takes a value above 5.0 Angstrom, it will feel a restraining potential of 7.0 kcal/mol/Angstrom^2 that will push it back.

This type of restraint is currently only possible for 'bond'/'distance' and 'rmsd' restraints.

######################################################
Parallelization of metadynamics
######################################################

One can of course control the number of CPU-cores in the ASHTheory level as usual which will affect how long each timestep will take.
For an MM simulation, it is best to run OpenMM on the GPU instead of CPU (platform='CUDA' or 'OpenCL').

However, it is even better to parallelize a metadynamics simulation via the multiwalker strategy:

By launching equivalent metadynamics simulation jobs (either simultaneously or at different times) but choosing a common biasdirectory 
that the different simulations will both read and write biasfiles from and to, one can extensive speed-up the exploration of the free-energy surface
and aid convergence.
Make sure to select a global biasdirectory that is available to all computing nodes  and then launch as many ASH-metadynamics jobs 
(e.g. using the **subash** submission script, see :doc:`basics`) as desired. 
Each "walker" simulation will write its bias-files to the common biasdirectory (according to the savefrequency keyword) and during each write-step
it will also read all bias-files and update the bias-potential. This will then influence the trajectory of each simulation and speed-up the build-up of the bias-potential.


######################################################
Examples:
######################################################

See also :doc:`mtd_tutorial` for working examples.


*1-CV metadynamics using the OpenMM native implementation:*

.. code-block:: python

  from ash import *

  #Name of biasdirectory (must exist)
  biasdir="./biasdirectory"

  #Creation of the ASH fragment
  frag = Fragment(databasefile="butane.xyz", charge=0, mult=1)

  #Create theory level. Here xTB using the in-memory library approach (no disk-based input or output)
  theory = xTBTheory(runmode='library')

  #Call the OpenMM metadynamics for 10K steps (each step being 0.001 ps = 1 fs)
  OpenMM_metadynamics(fragment=frag, theory=theory, timestep=0.001,
                simulation_steps=10000,
                temperature=300, integrator='LangevinMiddleIntegrator',
                coupling_frequency=1, traj_frequency=1, 
                CV1_atoms=[0,1,2,3], CV1_type='dihedral', CV1_biaswidth=0.5,
                biasfactor=6, height=1,
                frequency=1, savefrequency=1, biasdir=biasdir)


######################################################
MTD_analyze (for Plumed run): Analyze the results
######################################################

For metadynamics simulations utilizing the Plumed plugin, where the metadynamics results are available in the form of HILLS and COLVAR files it is possible
to use the **MTD_analyze** function to analyze the results and plot the data.

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