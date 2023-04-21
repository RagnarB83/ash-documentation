Parallelization in ASH
======================================

ASH can utilize parallelization in a few different ways: either via independent parallelization of the external QM program or MM program or via Python multiprocessing (running independent jobs in parallel).
For example if you create an **ORCATheory** object with numcores=X option, when ASH tells ORCA to run a calculation, ORCA will launch in parallel mode (as the ORCA inputfile created by ASH will contains parallelization information)
and will run its calculations in parallel using X cores using OpenMPI parallelization. 
This requires, however, OpenMPI to be set up (define PATH and LD_LIBRARY_PATH) correctly in the environment that ASH runs in (typically the jobscript used to submit calculations to the queuing system).
OpenMM and some QM programs may utilize simpler threads-based parallelization. The number of threads launched can usually be controlled by ASH via. numcores=X option. 

.. note:: for OpenMM the number of threads needs to be set up outside ASH (jobscript)

Some parts of ASH are parallelized by the Python multiprocessing library. This allows many independent calculations to be run simultaneously via the Pool feature of the multiprocessing library.
ASH job-functions that use multiprocessing parallelization are currently: **NumFreq**, **NEB**, **calc_surface**, **calc_surface_fromXYZ**.
These job-functions utilize the common **Job_parallel** function to parallelize over either multiple fragments or theories.
The **Job_parallel** function can also be run on its own in your own script.

######################################
Job_parallel
######################################

An alias for Job_parallel is Singlepoint_parallel.

.. code-block:: python

	def Job_parallel(fragments=None, fragmentfiles=None, theories=None, numcores=None, mofilesdir=None, 
							allow_theory_parallelization=False, Grad=False, printlevel=2, copytheory=False,
							version='multiprocessing', Opt=False, optimizer=None):


The **Job_parallel** function allows one to run many independent jobs in complete parallelization via the Python multiprocessing library. 
Typically the QM program parallelization is turned off in this case as it is more efficient to run as many calculations simultaneously as possible with each calculation utilizing a single core.
By default **Job_parallel** will run single-point energy calculations only but it is also possible to run Energy+Gradient calculations (use Grad=True)
or even geometry optimizations (use Opt=True).

Note that due to Job_parallel running multiple jobs simultaneously the printing during the job can be quite erratic.
This is normal.

**Example:**
if you have 120 single-point jobs to do (with roughly equivalent cost) and 24 cores available, it scales perfectly to occupy all CPU cores by 24 jobs at once (each job utilizing 1 core) and thus run through the list of 120 jobs in 5 batches.
This would be faster than running each job 1-by-1 utilizing QM-program parallelization (using 24 cores) as the QM-program parallelization will simply not scale as well (due to intrinsic parallelization limitations of the QM algorithms).
**Job_parallel** allows you to conveniently launch such parallelization jobs. The function distinguishes betwen 2 types of jobs: multiple fragments vs. multiple theories

- **multiple fragments with 1 theory**


For the more common case of multiple fragments, you may have a directory of e.g. 120 XYZ-files of different molecules and you want to run a singlepoint-energy job for each one. 
To do this using **Job_parallel** you would just need to create a list of ASH fragments for each XYZ-file and then pass the list and an ASH Theory level object to **Job_parallel**.
This can be easily accomplished via the script below where we make use of the convenient function **read_xyzfiles** to get an automatic list of all ASH fragment from a collection of XYZ-files.

.. code-block:: python

	from ash import *

	#Directory of XYZ files. Can be full path or relative path.
	xyzdir = '/path/to/xyz_files'

	#Creating list of ASH fragments from XYZ files. Using filename as label. 
	#NOTE: Using readchargemult=True, charge and mult will be read from the header of each XYZ-file.
	fragments = read_xyzfiles(xyzdir,readchargemult=True, label_from_filename=True)

	#Theory object
	ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=1)

	#Calling the Job_parallel function and providing list of fragments and theory:
	results = Job_parallel(fragments=fragments, theories=[ORCAcalc], numcores=4)

Each XYZ-file must have charge/mult information in the header of each file like this:

.. code-block:: text

    3
    0 1
    O       -1.377626260      0.000000000     -1.740199718
    H       -1.377626260      0.759337000     -1.144156718
    H       -1.377626260     -0.759337000     -1.144156718


- **multiple theories for 1 fragment**

For multiple theories you instead have to create a list of multiple **Theory** objects to be run on the single fragment.

.. code-block:: python

	from ash import *

	#Fragment for HBr
	hbr = Fragment(xyzfile="hbr.xyz", charge=0, mult=1)
	list_of_fragments=[hbr]

	#Create list of ORCATheory objects via for-loop
	list_of_orcaobjects=[]
	for functional in ['B3LYP', 'BP86', 'PBE0', 'M06', 'M06-2X', 'r2SCAN', 'SCAN', 'TPSS', 'PBE', 'PWLDA']:
	    ORCAcalc = ORCATheory(orcasimpleinput="! def2-SVP def2/J "+functional, orcablocks="", label=functional)
	    list_of_orcaobjects.append(ORCAcalc)

	#Calling the Job_parallel function 
	results = Job_parallel(fragments=list_of_fragments, theories=list_of_orcaobjects, numcores=4)

- **multiple theories for multiple fragments**

This option is currently not available for **Job_parallel**.


**Enabling QM-code parallelization**

There is also an option that allows both Python multiprocessing parallelization and the QMTheory parallelization to be active in a **Job_parallel** job. This option is turned off by default but can be enabled by the
*allow_theory_parallelization=True* keyword argument. However, care needs to be taken to make sure that the number of used CPU cores by ASH does not exceed the number of available CPU cores to the job (e.g. that requested by the queuing system). 

.. code-block:: python

	from ash import *

	#Defining some useful variables
	numcores = 8 #Total number of cores used by ASH. Should be equal to poolcores*QMcores. If using the subash script then this line is grepped.
	poolcores = 4 #The cores used by Job_parallel to run that many simultaneous jobs
	QMcores = 2 #How many cores are available to the external QM-code 

	xyzfiles_dir="/path/to/xyzfiles"

	#Creating list of ASH fragments from XYZ files. Using filename as label. 
	#Using readchargemult=True, charge and mult will be read from header of XYZ-file.
	fragments = read_xyzfiles(xyzfiles_dir,readchargemult=True, label_from_filename=True)

	orcacalc=ORCATheory(orcasimpleinput="! HF def2-SVP", numcores=QMcores)
	results = Job_parallel(theories=[orcacalc], fragments=fragments, numcores=poolcores, allow_theory_parallelization=True)


**Getting the results **

**Job_parallel**  like most other job-functions in ASH returns an ASH Results object (see :doc:`job-types` ).

The object always contains energies_dict which is a dictionaries (a label as key)
containing the energy of each job (e.g. fragment). If using Grad=True option there is also as gradients_dict.

The object also contains a dictionary (worker_dirnames) with the name of directories used for each job.

.. code-block:: python

	from ash import *

	results = Job_parallel(theories=[orcacalc], fragments=fragments, numcores=poolcores, allow_theory_parallelization=True)

	print(results)
	print(results.energies_dict)
	print(results.worker_dirnames)



**Running geometry optimizations in parallel **

Geometry optimizations (using geomeTRIC library) can also be run in parallel.
By setting Opt=True the function will run a geometry optimization instead of a Single-point calculation.
Default geometry optimization settings are employed. If you want to modify geometry optimization settings, however,
you need to provide an optimizer object. This could for example be used to set constraints for the optimization.

The optimized geometries are available as XYZ-files inside each worker dictory of the calculation.

.. code-block:: python

	from ash import *

	#Directory of XYZ files. Can be full path or relative path.
	xyzdir = '/path/to/xyz_files'

	#Creating list of ASH fragments from XYZ files. Using filename as label. 
	#NOTE: Using readchargemult=True, charge and mult will be read from the header of each XYZ-file.
	fragments = read_xyzfiles(xyzdir,readchargemult=True, label_from_filename=True)

	#Theory object
	ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=1)

	#Calling the Job_parallel function and providing list of fragments and theory:
	optimizer = GeomeTRICOptimizerClass(coordsystem='tric') #Creating an optimizer object
	results = Job_parallel(fragments=fragments, theories=[ORCAcalc], numcores=4, Opt=True, optimizer=optimizer)

	print("results object:", results)
	print("Energy dictionary:", result.energies_dict)

	#Getting fragments with optimized geometries from each worker-directory
	opt_fragments=[]
	for label,dirname in results.worker_dirnames.items():
		frag = Fragment(xyzfile=f"{dirname}/Fragment-optimized.xyz", label=label)
		opt_fragments.append(frag)
	print("opt_fragments:", opt_fragments)
