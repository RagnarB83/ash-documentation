Job Types
==========================

Jobs in ASH are performed by calling functions that take some input (usually an ASH **Fragment** and an ASH **Theory** object).

The primary job-types available in ASH:

* Single-point energy/property job
* Geometry optimization
* Numerical frequencies
* Analytic frequencies
* Nudged elastic band optimization (or saddlepoint search)
* Molecular dynamics
* Surface scan

These job-types are simple Python functions that are called with typically an ASH Theory and ASH Fragment as input.

Additionally ASH features various workflows that perform some combination of jobtypes above.
The job-types can be used with any theory object available (with some exceptions), e.g. one of the QMTheories in :doc:`QM-interfaces` or using
a QM/MM Theory object from :doc:`module_QM-MM`

################################
Output object of ASH Job-types
################################

Almost all the main job types in ASH now return the same object. 
This differs from previous versions where either an energy, geometry, fragment or dictionary might have been returned.

The return object is a `dataclass <https://realpython.com/python-data-classes/>`_  object of class ASH_Results. 
This object has various relevant attributes (energy, geometry, Hessian etc.) defined
with are by default set to None, while the job functions will set the relevant attributes depending on what was calculated.


.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``label``
     - string
     - None
     - A label that usually indicates what job-function created the object.
   * - ``energy``
     - float
     - None
     - The final energy. Set by Singlepoint, Optimizer, NEBTS.
   * - ``gradient``
     - Numpy array
     - None
     - A calculated gradient. Set by Singlepoint with Grad=True.
   * - ``reaction_energy``
     - float
     - None
     - A reaction energy in a chosen unit. Set by Singlepoint_reaction.
   * - ``energy_contributions``
     - dict
     - None
     - A dictionary of various energy contributions. Set by Singlepoint_reaction if a multistep theory like ORCA_CC_CBS_Theory was used.
   * - ``energies``
     - list of floats
     - None
     - A list of energies calculated. Set by Singlepoint_fragments, Singlepoint_theories, Singlepoint_fragments_and_theories and Singlepoint_reaction.
   * - ``reaction_energies``
     - list of floats
     - None
     - A list of reaction energies calculated. Set by Singlepoint_reaction.
   * - ``gradients``
     - List of numpy arrays.
     - None
     - A list of gradients calculated. Set by Singlepoint_parallel.
   * - ``energies_dict``
     - dict
     - None
     - Dictionary containing multiple energies for multiple single-point energy calculations. Set by Singlepoint_parallel.
   * - ``gradients_dict``
     - dict
     - None
     - Dictionary containing multiple gradients for multiple single-point energy calculations. Set by Singlepoint_parallel.
   * - ``geometry``
     - Numpy array
     - None
     - The final geometry from job. Set by Optimizer.
   * - ``initial_geometry``
     - Numpy array
     - None
     - The initial geometry from job. Set by Optimizer.
   * - ``hessian``
     - Numpy array
     - None
     - Hessian matrix. Set by NumFreq.
   * - ``frequencies``
     - List
     - None
     - List of vibrational frequencies in cm**-1. Set by NumFreq and AnFreq.
   * - ``vib_eigenvectors``
     - Numpy array
     - None
     - Eigenvectors from a mass-weighed Hessian diagonalization. Set by NumFreq.
   * - ``normal_modes``
     - Numpy array
     - None
     - Normal modes (unweighted eigenvectors). Set by NumFreq.
   * - ``thermochemistry``
     - dict
     - None
     - A dictionary containing various thermochemistry components ('ZPVE','Gcorr' etc.). Set by NumFreq and AnFreq.
   * - ``surfacepoints``
     - dict
     - None
     - Dictionary of energies of surfacepoints. Set by calc_surface and calc_surface_fromXYZ.
   * - ``saddlepoint_fragment``
     - ASH Fragment
     - None
     - An ASH fragment for the saddlepoint found. Set by NEB and NEBTS.
   * - ``MEP_energies_dict``
     - dict
     - None
     - Dictionary of total energies for each image. Set by NEB and NEBTS.
   * - ``barrier_energy``
     - float
     - None
     - The barrier height in kcal/mol (reactant->SP). Set by NEBTS.


If you are unsure what the Results object contains you can simply print it out and inspect

.. code-block:: python

  from ash import *

  #Geometry optimization of HF using ORCA
  HF_frag=Fragment(databasefile="hf.xyz", charge=0, mult=1)
  ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf")
  result = Optimizer(fragment=HF_frag, theory=ORCAcalc, coordsystem='tric')

  print(result) #Print whole object
  print(result.print_defined()) #Print only defined attributes (i.e. not None)

.. code-block:: text

  ASH_Results(label='Optimizer', energy=-100.354499689476, gradient=None, reaction_energy=None, energy_contributions=None, 
    energies=None, reaction_energies=None, gradients=None, energies_dict=None, gradients_dict=None, 
    geometry=array([[ 2.41430696e-10,  2.39715804e-09, -7.96015314e-03],[-2.41430696e-10, -2.39715804e-09,  9.24960153e-01]]),
    initial_geometry=None, charge=None, mult=None, hessian=None, frequencies=None, normal_modes=None, 
    vib_eigenvectors=None, thermochemistry=None, surfacepoints=None, reactant_geometry=None, product_geometry=None, 
    saddlepoint_geometry=None, saddlepoint_fragment=None, MEP_energies_dict=None, barrier_energy=None)

  Printing defined attributes of ASH_Results dataclass
  label: Optimizer
  energy: -100.354499689474
  geometry: [[ 6.66989696e-10  3.31567055e-09 -7.96014574e-03]
  [-6.66989696e-10 -3.31567055e-09  9.24960145e-01]]
  

**Using ASH_Results object:**

The ASH_results object can also be used for your own workflows

.. code-block:: python

  result = ASH_Results(label="MyLabel", energy=-100.354499689476)

**Writing ASH_Results object to disk:**

You can write the ASH_Results object to disk using a JSON format. Note that ASH usually does this automatically for all major jobtypes

.. code-block:: python

  result = ASH_Results(label="MyWorkflow", energy=-100.354499689476)
  result.write_to_disk(filename="ASH_myworkflow.result")

**Reading ASH_Results file :**

You can also read a previously written ASH_Results file back into a ASH_Results dataclass object.
This is a convenient option for reusing results from previous calculations.

.. code-block:: python

  result = ASH_Results(label="MyWorkflow", energy=-100.354499689476)
  result.write_to_disk(filename="ASH_myworkflow.result")
  read_result = read_results_from_file(filename="ASH_myworkflow.result")


###########################
Single-point calculation
###########################

The most basic jobtype. See :doc:`singlepoint`
In addition to the basic **Singlepoint** jobtype, there are also specialized functions: **Singlepoint_fragments**, **Singlepoint_theories**, 
**Singlepoint_fragments_and_theories**, **Singlepoint_reaction** and **Singlepoint_parallel** that are used to run single-point calculations
on multiple fragments or with multiple theories.

**Example:**

.. code-block:: python

    from ash import *

    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1) #Fragment object creation
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf") #ORCATheory object creation

    result = Singlepoint(fragment=HF_frag, theory=ORCAcalc)

###########################
Geometry optimization
###########################

Geometry optimizations are easily performed in ASH due to availability of the flexible optimizer: geomeTRIC (https://github.com/leeping/geomeTRIC): 

| See  :doc:`Geometry-optimization` documentation for all options.

The geomeTRIC **Optimizer** allows efficient optimization in multiple coordinate systems: TRIC, HDLC, DLC, Cartesian, redundant internals. Supports constraints as well as frozen atoms natively. 
Allows an active-region definition which enables efficient QM/MM optimizations of a part of large systems (where most atoms are frozen).
ASH features a full-featured interface to geomeTRIC that allows flexible constraint input, QM/MM optimizations, relaxed and unrelaxed 1D/2D surface scans and more.

**Example:**

.. code-block:: python

    from ash import *

    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1) #Fragment object creation
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf") #ORCATheory object creation

    #Geometry optimization of the ORCA using geomeTRIC optimizer
    Optimizer(fragment=HF_frag, theory=ORCAcalc, coordsystem='tric')
    #Optimizer and Opt are aliases for the geomeTRICOptimizer function name.

See :doc:`Geometry-optimization` for all features.

Other optimizers:

- An internal optimizer is available (called **SimpleOpt**) that can optimize the system in Cartesian coordinates only using the LBFGS algorithm. While frozen atoms are supported, no other constraints are supported.

################################
Numerical frequencies (Hessian)
################################

Numerical frequencies can be performed with ASH using any QM, MM or QM/MM theory object. Parallelization is available.
See :doc:`module_freq` documentation for all options.

**Example:**

.. code-block:: python

    from ash import *

    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1) #Fragment object creation
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf") #ORCATheory object creation

    NumFreq(fragment=HF_frag, theory=ORCAcalc)

################################
Analytical frequencies (Hessian)
################################
Analytical frequencies can be requested in some cases if supported by the theory-level interface as well as the Hamiltonian inside program.
See :doc:`module_freq`


##################################
Nudged Elastic Band Calculations
##################################

Through an interface to the external code Knarr, nudged elastic band (NEB) calculations are possible.
This enables one to calculate minimum energy paths and locate saddlepoints ("transition states") using any QM, MM or QM/MM Theory in ASH.

See :doc:`neb` for documentation

**Example:**

.. code-block:: python

    from ash import *

    Reactant=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    Product=Fragment(xyzfile="prod.xyz", charge=0, mult=1)

    #Calculator object without frag
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library')

    NEB(reactant=Reactant, product=Product, theory=xtbcalc, images=10, CI=True)


###########################
Surface scans
###########################
Potential Energy Surfaces can be conveniently scanned in ASH using the **calc_surface function** that uses the **geometric** optimization library.
Both unrelaxed and relaxed scans be calculated, using either 1 and 2 reaction coordinates.

See :doc:`surfacescan`



###########################
Saddle-point optimization
###########################

Saddle-points searches can be be performed in ASH via a double-ended strategy (requiring reactant and product starting points) and a single-ended strategy (requiring only a single geometry).
The double-ended strategy involves use of the climbing image NEB method which also results in a minimum energy path between reactant and product.
See :doc:`neb` for documentation.

An eigenvector-following algorithm is also available via the geomeTRIC library (TSOpt=True option). This option is only feasible when a good guess for the 
saddlepoint geometry is available, e.g. from a surface scan, previous NEB/NEB-CI job etc. It furthermore requires a good initial approximation to the Hessian (default: exact Hessian in first step).
See :doc:`Geometry-optimization` for all features.

**Example:**

.. code-block:: python

    from ash import *

    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1) #Fragment object creation
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf") #ORCATheory object creation

    #TSOpt=True enables saddlepoint optimization in geomeTRIC. Note: Exact Hessian is calculated in the first step by default.
    Optimizer(fragment=HF_frag, theory=ORCAcalc, coordsystem='tric', TSOpt=True)


-----------------------------------
**NEB-TS**
-----------------------------------

A combination of the double-ended NEB strategy and a single-ended eigenvector-following method is also available in ASH in the form of the NEB-TS method.
This is probably one of the most efficient and accurate method for finding a saddlepoint as discussed in the article:

V. Ásgeirsson, B. Birgisson, R. Bjornsson, U. Becker, F. Neese, C: Riplinger,  H. Jónsson, J. Chem. Theory Comput. 2021,17, 4929–4945.
DOI: 10.1021/acs.jctc.1c00462

See :doc:`neb` for documentation on the NEB-TS function.

**Example:**

.. code-block:: python

    from ash import *

    Reactant=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    Product=Fragment(xyzfile="prod.xyz", charge=0, mult=1)
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf") #ORCATheory object creation

    #NEB-TS combines a CI-NEB job (note: looser thresholds than default CI-NEB) and a Optimizer(TSOpt=True) job.
    result = NEBTS(reactant=Reactant, product=Product, theory=ORCAcalc, images=12, printlevel=0)


###########################
Molecular Dynamics
###########################

It is possible to perform molecular dynamics in ASH using the interface to OpenMM that can be used for MM, QM and QM/MM simulations.
 
See :doc:`module_dynamics` and also  :doc:`Biased-sampling`



