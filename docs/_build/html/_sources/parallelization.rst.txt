======================================
Parallelization in ASH
======================================

ASH can utilize parallelization in a few different ways: either via independent parallelization of the external QM program or MM program or via Python multiprocessing (running independent jobs in parallel).
For example if you create an ORCATheory object with numcores=X option, when ASH tells ORCA to run a calculation, ORCA will launch in parallel mode (as the ORCA inputfile created contains parallelization information)
and will run its calculations in parallel using OpenMPI parallelization. 
This requires, however, OpenMPI to be set up (define PATH and LD_LIBRARY_PATH) correctly in the environment that ASH runs in (typically the jobscript used to submit calculations to the queuing system).
OpenMM and some QM programs may utilize simpler threads-based parallelization. ASH can sometimes control the number of threads launched (also via. numcores=X option). 

note: for OpenMM the number of threads needs to be set up outside ASH (jobscript)


Some parts of ASH are parallelized by the Python multiprocessing library. This allows many independent calculations to be run simultaneously via the Pool feature of the multiprocessing library.


ASH modules that use multiprocessing parallelization are currently: NumFreq and Singlepoint_parallel (see below).


######################################
Singlepoint_parallel
######################################


.. code-block:: python

	def Singlepoint_parallel(fragments=None, fragmentfiles=None, theories=None, numcores=None, mofilesdir=None):


The Singlepoint_parallel function allows one to run many independent single-point energy jobs in complete parallelization via the Python multiprocessing library. 
Typically the QM program parallelization is turned off in this case as it is more efficient to run as many calculations simultaneously as possible with each calculation utilizing a single core.
Example: if you have 120 single-point jobs to do (with roughly equivalent cost) and 24 cores available, it scales perfectly to occupy all CPU cores by 24 jobs at once (each utilizing 1 core) and thus run through the list of 120 jobs in 5 batches.
This would be faster than running each job 1-by-1 utilizing QM-program parallelization (using 24 cores) as the QM-program parallelization will simply not scale as well.

Singlepoint_parallel allows you to conveniently launch such parallelization jobs. The function distinguishes betwen 2 job-types: multiple fragments vs. multiple theories

- **multiple fragments with 1 theory**


For the more common case of multiple fragments, you may have a directory of e.g. 120 XYZ-files of different molecules and you want to run a singlepoint-energy job for each one. 
To do this using Singlepoint_parallel you would just need to create a list of ASH fragments for each XYZ-file and then pass the list and an ASH Theory level object to Singlepoint_parallel.

.. code-block:: python

	from ash import *

	#Directory of XYZ files. Can be full path or relative path.
	xyzdir = '/path/to/xyz_files'

	list_of_molecules=[]
	#Creating list of ASH fragments from XYZ files. Using filename as label
	for file in glob.glob(xyzdir+'/*.xyz'):
	    print("XYZ-file:", file)
	    basename=os.path.basename(file)
	    label=os.path.splitext(basename)[0]
	    molecule=Fragment(xyzfile=file,label=label)
	    #molecule=Fragment(xyzfile=file,label=label, readchargemult=true)
	    list_of_molecules.append(molecule)

	#Theory object
	ORCAcalc = ORCATheory(charge=0, mult=1, orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", numcores=1)

	#Calling the Singlepoint_parallel function and providing list of fragments and theory:
	results = Singlepoint_parallel(fragments=list_of_molecules, theories=[ORCAcalc], numcores=4)

NOTE: Here assuming that each molecule has charge=0 and mult=1. If this is not the case, one can instead read in the charge and multiplicity from the XYZ-title line (line 2) for each molecule like this in the for-loop above:

molecule=Fragment(xyzfile=file,label=label, readchargemult=true).

Each XYZ_file should then contain charge and multiplicity in line number 2 (space separated).



- **multiple theories for 1 fragment**

For multiple theories you instead have to create a list of multiple Theory objects to be run on the single fragment.

.. code-block:: python

	from ash import *

	#Fragment for HBr
	hbr = Fragment(xyzfile="hbr.xyz")
	list_of_fragments=[hbr]

	#Create list of ORCATheory objects
	list_of_orcaobjects=[]
	for functional in ['B3LYP', 'BP86', 'PBE0', 'M06', 'M06-2X', 'r2SCAN', 'SCAN', 'TPSS', 'PBE', 'PWLDA']:
	    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput="! def2-SVP def2/J "+functional, orcablocks="", label=functional)
	    list_of_orcaobjects.append(ORCAcalc)

	#Calling the Singlepoint_parallel function 
	results = Singlepoint_parallel(fragments=list_of_fragments, theories=list_of_orcaobjects, numcores=4)

- **multiple theories for multiple fragments**

TODO: This is currently not available


