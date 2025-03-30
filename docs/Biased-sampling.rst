Biased sampling MD & Free energy simulations
===============================================

Biased/enhanced sampling or free-energy simulations are possible in ASH via the OpenMM molecular dynamics routines (see :doc:`module_dynamics` and :doc:`OpenMM-interface`).

Metadynamics has become a very popular biased sampling / free energy simulation method due to its ease-of-use and handy parallelization strategy.
It is possible to perform metadynamics simulations in ASH using essentially any theory-level, including an OpenMMTheory level, any type of QMTheory level and a QMMMTheory level.
The OpenMM native metadynamics routines are used for the simulations, and in the case of QM and QM/MM Theories, the energy and forces are passed onto OpenMM in each timestep.

The **OpenMM_metadynamics** function can be used to start a metadynamics simulation from an ASH Fragment and ASH Theory level and some collective variable information, using the metadynamics functionality inside the OpenMM library.
See `OpenMM implementation <http://docs.openmm.org/development/api-python/generated/openmm.app.metadynamics.Metadynamics.html>`_ .
The function sets up the necessary collective-variable bias potentials before launching an MD simulation using the **OpenMM_MD** class.

Umbrella sampling simulations in ASH can be performed by adding restraint potentials to the OpenMMTheory object before running MD. 

Additionally ASH features an interface (via OpenMM) to the popular enhanced-sampling library PLUMED that has extensive support for various algorithms, free-energy methods, analyzis tools etc.
Plumed-based simulations use the  **OpenMM_MD_plumed** function which activates the interface between OpenMM and the powerful `Plumed library <https://www.plumed.org>`_ .
This option requires the installation of the `OpenMM PLUMED plugin <https://github.com/openmm/openmm-plumed>`_ and also `Plumed <https://www.plumed.org>`_ .

The first option is generally recommended when possible, as it is faster (no OpenMM-Plumed communication per timestep required) and does not require additional software.
The most useful collective variable options are available (distances, angles, torsion, RMSD, coordination-number).
In addition it is possible to define your own collective variables (see later on this page) by the useful Custom-forces options inside OpenMM.
However, the PLUMED-OpenMM option offers even more flexiblity as in principle all collective variables inside the Plumed library can be used 
(as long as they depend only on geometric information). Additionally, other enhanced sampling methods inside Plumed can in principle also be used.

######################################################
OpenMM_metadynamics 
######################################################

This function utilizes the metadynamics options that are built into the OpenMM library.
See also :doc:`mtd_tutorial` for working examples on how to perform metadynamics in ASH.

**OpenMM_metadynamics**  has all the same options as the **OpenMM_MD** function (see :doc:`module_dynamics`) but has in addition some special metadynamics keyword options.

.. code-block:: python

  def OpenMM_metadynamics(fragment=None, theory=None, timestep=0.004, simulation_steps=None, simulation_time=None,
                traj_frequency=1000, temperature=300, integrator='LangevinMiddleIntegrator',
                barostat=None, pressure=1, trajectory_file_option='DCD', trajfilename='trajectory',
                coupling_frequency=1, charge=None, mult=None, platform='CPU',
                anderson_thermostat=False, restraints=None, 
                enforcePeriodicBox=True, dummyatomrestraint=False, center_on_atoms=None, solute_indices=None,
                datafilename=None, dummy_MM=False, add_center_force=False,
                center_force_atoms=None, centerforce_constant=1.0, barostat_frequency=25, specialbox=False,
                CV1_atoms=None, CV2_atoms=None, CV1_type=None, CV2_type=None, biasfactor=6, 
                height=1, flatbottom_restraint_CV1=None, flatbottom_restraint_CV2=None,
                CV1_biaswidth=0.5, CV2_biaswidth=0.5, CV1_range=None, CV2_range=None,
                frequency=1, savefrequency=10,
                biasdir='.'):


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
     - Name of CV-type for CV1. Options: 'bond' (or 'distance'), 'angle', 'dihedral' (or 'torsion'), 'rmsd' or 'custom'.
   * - ``CV2_type``
     - string
     - None
     - Name of CV-type for CV2. Options: 'bond' (or 'distance'), 'angle', 'dihedral' (or 'torsion'), 'rmsd' or 'custom'.
   * - ``CV1_biaswidth``
     - float
     - None
     - Biaswidth for CV1 (sigma) in CV units. CV-unit is Angstrom for 'bond' and 'rmsd' but radians for 'angle' and 'dihedral'. 
   * - ``CV2_biaswidth``
     - string
     - None
     - Biaswidth for CV2 (sigma) in CV units. CV-unit is Angstrom for 'bond' and 'rmsd' but radians for 'angle' and 'dihedral'.
   * - ``reference_xyzfile``
     - string
     - None
     - The reference structure (path to XYZ-file) to use when CV is 'rmsd'. If not set then initial coordinates are used as reference.
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

-----------------
Defining CVs
-----------------

The *CV1_type* and *CV2_type* keyword should specify the type of CV as a string.
The built-in options are: 'bond' (a.k.a. 'distance'), 'angle', 'dihedral' (a.k.a. 'torsion'), 'rmsd' and 'cn'. And 'custom' allows specifying your own CV.

Depending on the CV-type you then must specify the atoms that define that geometric variable 
The biaswidth for the CV must also be chosen and it is important to note that an appropriate value depends very much on the CV-type.
The unit for the biaswidth is Angstrom for 'bond' and 'rmsd' CVs, but radians for 'angle' and 'dihedral'. 
The *CV1_range* and *CV2_range* keywords define the min-to-max range that the CV can take (same unit as biaswidth).

Distance, angle and torsion CV definitions:

.. code-block:: python

  #Defining a distance CV between atom 0 and 1. Biaswidth 0.01 Ang
  CV1_atoms=[0,1], CV1_type='distance', CV1_biaswidth=0.01,
  #Defining an angle CV between atom 0, 1 and 2. Biaswidth 0.5 radians
  CV1_atoms=[0,1,2], CV1_type='angle', CV1_biaswidth=0.5,
  #Defining a torsion CV between atom 0 and 1. Biaswidth 0.5 radians
  CV1_atoms=[0,1,2,3], CV1_type='torsion', CV1_biaswidth=0.5,

The 'rmsd' option uses the RMS-difference between the current structure and a reference structure as a CV.
One should point *reference_xyzfile* keyword to an XYZ-file in this case.
If *reference_xyzfile* is not set, then the starting structure is used as reference structure instead.

.. code-block:: python

  #Defining an RMSD CV, using atoms 0,2,5,7 . Biaswidth 0.01 Ang
  CV1_atoms=[0,2,5,7], CV1_type='rmsd', CV1_biaswidth=0.01, reference_xyzfile="ref.xyz"


The 'cn' option defines a coordination number CV that can be highly useful, e.g. for protonation reactions.
In example below, we define a CN-CV for e.g. a nitrogen group that can be either R-NH3+ or R-NH2.
We define bonds between the nitrogen (index 6 below) and all acidic hydrogens (13,14,15).

.. code-block:: python

  #Defining a CN-CV, defining bonds between atoms (6,13),(6,14),(6,15). Biaswidth 0.01 (unitless)
  CV1_atoms=[[6,13],[6,14],[6,15]], CV1_type='cn', CV1_biaswidth=0.01, CV1_parameters=[2.0],

If all H-atoms are close to N, the CN will be close to 3 while if one of the hydrogens transfer to another site, the CN will be closer to 2.
A threshold value (r0 in equation below) needs to also be passed (roughly defines when the distance is no longer a bond). Here we do this using the CV1_parameters list.
The mathematical expression used to define the coordination number: 

.. math::

   N(r) = \Sigma_{i \in g_{1}} \Sigma_{j \in g_{2}} S ( \frac{|| r_{j} - r_{i} ||} {r_{0}})  

where S is a simplification of the common step function used for CNs:

.. math::

   S = \frac{1-x^{6}} {1-x^{12}} = \frac{1}{1+x^{6}}


Finally, the 'custom' option allows the user to define their own CV by utilizing the powerful Custom-force options inside the OpenMM Library.
First an OpenMM Force object must be created (using an appropriate OpenMM CustomForce and energy expression) and appropriate atoms or groups added to the Force.
See `OpenMM Custom Forces <http://docs.openmm.org/latest/userguide/theory/03_custom_forces.html>`_ for more information.
The OpenMM BiasVariable should then be defined: see `BiasVariable documentation <http://docs.openmm.org/latest/api-python/generated/openmm.app.metadynamics.BiasVariable.html?highlight=biasvariable>`_
Both the CV-force created and the biasvariable-object should then be passed as arguments to *user_cvforce1* and  *user_biasvar1* (or *user_cvforce2* and  *user_biasvar2*).

Example below shows how to define the coordination number CV as a custom-force instead of the built-in option.
This is just one example, the OpenMM CustomForces are flexible enough to allow definitions of highly complex CVs.

.. code-block:: python

  import openmm
  #Defining custom cvforce
  energy_expression="1/(1+x^6) ; x=r/threshold"
  cvforce = openmm.CustomBondForce(energy_expression)
  #Threshold that defines when a bond is present
  cvforce.addGlobalParameter("threshold", 2.0*openmm.unit.angstrom)

  #Adding the atoms that define each bonds
  cvforce.addBond(6,13)
  cvforce.addBond(6,14)
  cvforce.addBond(6,15)

  #Creating Biasvariable: forceobj, minval, maxval, biaswidth
  biasvar=openmm.app.BiasVariable(cvforce, 00, 5, 0.05, periodic=False)

  # Choosing custom CV option and and specifying the CV-force and Biasvariable
  OpenMM_metadynamics(...,CV1_type='custom', user_cvforce1=cvforce, user_biasvar1=biasvar)


-----------------
Restraining CVs
-----------------

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
OpenMM_MD_plumed 
######################################################


**OpenMM_MD_plumed** works similarly to **OpenMM_metadynamics** above but instead of defining collective variables using keyword arguments (i.e. CV1_atoms and CV1_type) we instead define a multi-line Plumed-string (*plumed_input_string* keyword) that will be passed to the Plumed libary.
This string should define all collective variable information by Plumed-syntax as well as defining specific Plumed metadynamics options etc.
Since the interface is quite general and Plumed handles all CV-setup and passes the biasing potential back to OpenMM, we can in principle define almost anything inside this Plumed-string.
This means that in principle some other enhanced sampling (not just metadynamics) methods available inside the Plumed library can also be used (not tested).

.. code-block:: python

  def OpenMM_MD_plumed(fragment=None, theory=None, timestep=0.001, simulation_steps=None, simulation_time=None,
                traj_frequency=1000, temperature=300, integrator='LangevinMiddleIntegrator',
                barostat=None, pressure=1, trajectory_file_option='DCD', trajfilename='trajectory',
                coupling_frequency=1, charge=None, mult=None, platform='CPU', hydrogenmass=1.5, constraints=None,
                anderson_thermostat=False, restraints=None, 
                enforcePeriodicBox=True, dummyatomrestraint=False, center_on_atoms=None, solute_indices=None,
                datafilename=None, dummy_MM=False, add_centerforce=False,
                centerforce_atoms=None, centerforce_distance=10.0, centerforce_constant=1.0, centerforce_center=None,
                barostat_frequency=25, specialbox=False,
                plumed_input_string=None, printlevel=2, numcores=1):

The *plumed_input_string* variable should define a multiline-string, this will be written to disk as a *plumedinput.in* file and the information is passed to the Plumed library.
The information in this string should define the desired collective variables as well as define all metadynamics options. Note also that the update frequency of the bias-potential and print-out (PACE and STRIDE options) should be defined here.
Information about coordinates should generally not be present here (an exception is when defining a reference geometry such as when using the RMSD CV).
In short: everything to do with the CVs and bias-potential needs to be defined in the Plumed input-string. Bias-widths and heights etc. will of course use Plumed units so the Plumed documentation should be used for help in defining things correctly.

During each MD simulation step information about coordinates and velocities is passed to Plumed that in turn passes information about the bias potential back.
Plumed additionally will write the information about the bias-potential and current CV-values to disk using it's own syntax (typically the files are called HILLS and COLVAR). These files can be used for deriving and plotting the free-energy surface, see later.


.. warning:: Be aware that unlike ASH and OpenMM (where atom indices are counted from zero), PLUMED atom indices start from 1. This needs to be taken into account when defining CVs by atom indices. 




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
                frequency=1, savefrequency=10, biasdir=biasdir)


*1-CV metadynamics using the OpenMM-Plumed implementation:*

Note that here we perform an identical MTD simulation on butane with the same theory-level but using the Plumed interface instead.

.. code-block:: python

  from ash import *

  #Creation of the ASH fragment
  frag = Fragment(databasefile="butane.xyz", charge=0, mult=1)

  #Create theory level. Here xTB using the in-memory library approach (no disk-based input or output)
  theory = xTBTheory(runmode='library')

  # Warning: atom indices inside plumed-string start from 1 (not 0 like in ASH)
  plumedstring="""
  # set up two variables for Phi and Psi dihedral angles
  phi: TORSION ATOMS=1,2,3,4

  # Activate metadynamics in phi
  # depositing a Gaussian every 500 time steps,
  # with height equal to 1.2 kJ/mol,
  # and width 0.35 rad
  #
  metad: METAD ARG=phi PACE=500 HEIGHT=1.2 SIGMA=0.35 FILE=HILLS
  # monitor the two variables and the metadynamics bias potential
  PRINT STRIDE=10 ARG=phi,metad.bias FILE=COLVAR
  """

  #Call the OpenMM_MD_plumed function
  OpenMM_MD_plumed(fragment=frag, theory=theory, timestep=0.001,
                simulation_steps=10000,
                temperature=300, integrator='LangevinMiddleIntegrator',
                coupling_frequency=1, traj_frequency=1,
                plumed_input_string=plumedstring)

###########################################################
Parallelization: Running multiple walker MTD simulations
###########################################################

In order to parallelize metadynamics simulations, one can of course control the number of CPU-cores in the ASHTheory level 
as usual which will affect how long each timestep will take (note that MM simulations, running on the GPU (platform='CUDA' or 'OpenCL') is much preferable to the CPU).
However, the multiple walker strategy works much better than the Theory parallelization.

As shown in the :doc:`mtd_tutorial` tutorial, running multiple walker metadynamics simulations quickly reduces the sampling error and 
allows faster convergence of the free energy surface.
This can be accomplished in a very simple way, one simply has to launch multiple simulations at the same time, while making sure
that each simulation uses the same shared biasdirectory, allowing information about the continuously built-up biasing potential
to be shared among each simulation.

In the example below, we use the same input as above but change the biasdirectory to a global directory accessible by all nodes.

.. code-block:: python

  from ash import *

  #Name of a globally shared biasdirectory (must be accessible by all nodes)
  # Best to use something that exists in e.g. your /home or /work directory if you are on a cluster
  biasdir="/home/username/metadynamics_simulations/biasdirectory"

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
                frequency=1, savefrequency=10, biasdir=biasdir)

We could then go ahead and submit this script to a cluster multiple times.
In order for each job to write to a unique outputfile, it might be best to create copies of the ASH inputscript, 
and then submit like this (here using the **subash** (see :doc:`basics` ) submission script):

.. code-block:: shell
  
  subash mtd_sim1.py -p 1 # Here requesting 1 CPU core

Each job will probably end up on a different node, writing most temporary files to its own local scratch but will 
read and write bias-potential information on the shared biasdirectory. Pay attention to the *savefrequency* variable of **OpenMM_metadynamics** as it controls
how often the bias is read and written to disk. The more often, the more up-to-date the bias-potential will be but this may read to excessive read/write operations that will 
slow down the simulation and may lead to excessive network traffic on the cluster (especially if you are running MM metadynamics).

.. note:: The multiple-walker approach should also work for OpenMM_MD_plumed jobs but requires more setup. 
  The biasdirectory should be set by PLUMED keyword WALKERS_DIR, the number of walkers by WALKERS_N etc. See Plumed documentation.

The advantage of the approach above is that you can submit multiple walker-jobs, 
perhaps using different CPU cores for each simulation (to speed up the theory energy+gradient step), 
depending on the resources that are available. You can submit jobs whenever you want, even after all other jobs have finished and continue a previous metadynamics simulation.

Another scenario might come up where you might want to submit to a single computing node that has e.g. 24 cores and you wish to run 24 walkers on that node automatically.
Here we assume each walker will run with 1 CPU core.

The **subash** script (see :doc:`basics` ) has an option to automatically submit multiple-walker ASH calculations on a single node, using a single job submissions.

.. code-block:: shell
  
  subash mtd_sim.py -mw -p 24 # mw will launch multiple walkers on a single node, -p 24 will request 24 CPU cores


If you inspect the `subash script <https://github.com/RagnarB83/ash/blob/master/scripts/subash.sh>`_  (search for multiwalker) 
you can see the logic of what will be done on the node.
Briefly: upon submission of a job to the queuing system (subash assumes SLURM), before launching ASH, a separate directory will be 
created for each walker (here 24 in total), named walkersim1, walkersim2, etc...
All files originally copied to local scratch will be copied into each directory and then an ASH calculation will start inside each directory 
simultaneously (24 in total in this example). Note that only a single job will be submitted to the queuing system, however, but 24 Python processes will be running on the node.
Once each simulation has finished, the job ends and all contents are copied back to the submission directory.
This is a highly convenient way of launching multiple walkers on the same single node.

.. note:: Running as many multiple walkers as possible should generally be preferable to speeding up the energy+gradient step of the Theory level.


######################################################
Plotting the results
######################################################


------------------------
metadynamics_plot_data
------------------------
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
- MTD_CV1.png/MTD_CV1_CV2.png # Image containing the final plot (requires Matplotlib)


------------------------------------------------------------
plumed_MTD_analyze (for Plumed run): Analyze the results
------------------------------------------------------------

For metadynamics simulations utilizing the Plumed plugin, where the metadynamics results are available in the form of HILLS and COLVAR files it is possible
to use the **plumed_MTD_analyze** function to analyze the results and plot the data.

.. code-block:: python

  def plumed_MTD_analyze(plumed_ash_object=None, path_to_plumed=None, Plot_To_Screen=False, CV1_type=None, CV2_type=None, temperature=None,
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
     - Type of CV1. Options: 'bond' (or 'distance'), 'angle', 'dihedral' (or 'torsion') or 'rmsd'.
   * - ``CV2_type``
     - string
     - None
     - Type of CV2. Options: 'bond' (or 'distance'), 'angle', 'dihedral' (or 'torsion') or 'rmsd'.
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
Umbrella sampling in ASH
######################################################

ASH can also be used for umbrella sampling simulations by adding a restraint potential to the system
before running MD. See section "Adding custom forces to MD simulation" in  :doc:`module_dynamics` for more information.


A workflow and tutorial for running umbrella sampling is currently missing but a basic example is shown below.
The post-processing could then in principle be performed by a method such as MBAR or WHAM using a program such as `FastMBAR <https://fastmbar.readthedocs.io>`_

.. code-block:: shell

  from ash import *
  import os
  import math

  # Example script for performing a basic umbrella sampling simulation using ASH
  # Note:
  # System: Butane torsion using GFN1-XTB

  ####################################################################
  # Creating the ASH fragment
  frag = Fragment(databasefile="butane.xyz", charge=0, mult=1)
  # Defining the xTB theory (GFN1-xTB)
  theory = xTBTheory(runmode='library')
  ####################################################################

  # US restraint potential settings
  RC_atoms=[0,1,2,3]
  RC_FC=2000 #Unit?
  traj_frequency=10 #Frames saved to trajectory and used in US
  filename_prefix="US_window" # Used for created files
  M = 20 # M centers of harmonic biasing potentials
  theta0 = np.linspace(-math.pi, math.pi, M, endpoint = False) # array of values

  # Save US settings to parameterfile
  import json
  json.dump({'M':M, 'RC_atoms':RC_atoms,'RC_FC':RC_FC, 'traj_frequency':traj_frequency,
            'theta0':list(theta0), 'filename_prefix':filename_prefix}, open(f"ASH_US_parameters.txt",'w'))

  # Loop over windows and run biased simulation in each
  # Note: More efficient to run these as independent simulations in parallel
  for ind,RC_val in enumerate(theta0):
      print("="*50)
      print(f"NEW UMBRELLA WINDOW. Value: {RC_val}")
      print("="*50)
      # Setting restraint potential as a list: [atom_indices, value, force constant]
      restraint=RC_atoms+[RC_val]+[RC_FC] # Combining into 1 list

      # Calling OpenMM_MD with a restraint potential
      OpenMM_MD(fragment=frag, theory=theory,
              timestep=0.001, simulation_time=1, traj_frequency=traj_frequency,
              temperature=300, restraints=[restraint])

      os.rename("trajectory.dcd", f"{filename_prefix}_{ind}.dcd")