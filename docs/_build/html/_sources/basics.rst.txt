==========================
Basic usage
==========================

#####################
Input structure
#####################
You create a Python3 script (e.g. called system.py) and import the Ash functionality:

.. code-block:: python

    from ash import *


Ash functionality can only be imported if the Ash source dir is in the PYTHONPATH. Make sure you have set in the shell:

.. code-block:: shell

    export PYTHONPATH=/path/to/ash_dir:$PYTHONPATH

For convenience you may want to initalize standard global settings (connectivity etc.):

.. code-block:: python

    settings_ash.init()

The global settings are stored in your *ash-dir/settings_ash.py* and can be modified.

From then on you have the freedom of writing a Python script in whatever way you prefer but taking the advantage
of Ash functionality. Typically you would first create one (or more) molecule fragments, then define a theory
object and then call a specific job-module (an optimizer, numerical-frequencies, MD).
See  :doc:`coordinate-input` for various ways of dealing with coordinates and fragments.

#####################
Example script
#####################

Here is a basic Ash Python script, e.g. named: ashtest.py

.. code-block:: python

    from ash import *
    settings_ash.init()

    #Create fragment
    Ironhexacyanide = Fragment(xyzfile="fecn6.xyz")

    #Defining ORCA-related variables
    orcadir='/opt/orca_4.2.1'
    orcasimpleinput="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    orcablocks="%scf maxiter 200 end"

    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1,
                                orcasimpleinput=orcasimpleinput, orcablocks=orcablocks)

    #Basic Cartesian optimization with KNARR-LBFGS
    Optimizer(fragment=Ironhexacyanide, theory=ORCAcalc, optimizer='KNARR-LBFGS')


The script above loads Ash, creates a new fragment from an XYZ file (see :doc:`coordinate-input` for other ways),
defines variables related to the ORCA-interface , creates an ORCA-theory object
(see :doc:`QM-interfaces`), defines an Optimizer object and finally runs a geometry
optimization  (see :doc:`job-types` for other options).

########################
Running script directly
########################

For a simple job we can just run the script directly

.. code-block:: shell

    python3 ashtest.py

The output will be written to standard output (i.e. your shell). You can redirect the output to a file.

.. code-block:: shell

    python3 ashtest.py >& ashtest.out


#####################
Submitting job
#####################

For a more complicated job we would probably want to create a job-script that would handle various environmental variables,
dealing with local scratch, copy files back when done etc.
Here is an example SLURM jobscript. Remember to go through all the lines and change the various things like the path to
local scratch, set the correct PATH variables, load modules etc.

.. code-block:: shell

    #!/bin/bash

    #SBATCH -N 1
    #SBATCH --tasks-per-node=12
    #SBATCH --time=8760:00:00
    #SBATCH -p compute
    #SBATCH --mem-per-cpu=3000
    #SBATCH --job-name=ASHJOB
    #SBATCH --output=%x.o%j
    #SBATCH --error=%x.o%j

    #Use like this:
    #sbatch -J inputfile.py jobscript.sh

    export job=$SLURM_JOB_NAME
    export job=$(echo ${job%%.*})

    #Outputname
    outputname="$job.out"

    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export OMP_STACKSIZE=1G
    export OMP_MAX_ACTIVE_LEVELS=1

    #Create scratch directory on local scratch
    path_to_scratch=/scratch
    if [ ! -d $path_to_scratch/$USER ]
    then
      mkdir -p $path_to_scratch/$USER
    fi
    tdir=$(mktemp -d $path_to_scratch/$USER/ashjob__$SLURM_JOB_ID-XXXX)
    chmod +xr $tdir


    #Copy all relevant inputfiles for ASH: python scripts, CIF-files, XYZ files etc.
    cp $SLURM_SUBMIT_DIR/*.py $tdir/
    cp $SLURM_SUBMIT_DIR/*.cif $tdir/
    cp $SLURM_SUBMIT_DIR/*.xyz $tdir/
    cp $SLURM_SUBMIT_DIR/*.xtl $tdir/
    cp $SLURM_SUBMIT_DIR/*.ff $tdir/
    cp $SLURM_SUBMIT_DIR/*.ygg $tdir/
    cp $SLURM_SUBMIT_DIR/*.pdb $tdir/
    cp $SLURM_SUBMIT_DIR/*.hess $tdir/
    cp $SLURM_SUBMIT_DIR/*.info $tdir/
    cp $SLURM_SUBMIT_DIR/Centralmainfrag $tdir/

    # cd to scratch
    cd $tdir
    echo "tdir is $tdir"

    # Copy job and node info to beginning of outputfile
    echo "Starting job in scratch dir: $tdir" > $SLURM_SUBMIT_DIR/$outputname
    echo "Job execution start: $(date)" >> $SLURM_SUBMIT_DIR/$outputname
    echo "Shared library path: $LD_LIBRARY_PATH" >> $SLURM_SUBMIT_DIR/$outputname
    echo "Slurm Job ID is: ${SLURM_JOB_ID}" >> $SLURM_SUBMIT_DIR/$outputname
    echo "Slurm Job name is: ${SLURM_JOB_NAME}" >> $SLURM_SUBMIT_DIR/$outputname
    echo $SLURM_NODELIST >> $SLURM_SUBMIT_DIR/$outputname

    #ASH environment

    #Load necessary modules.
    #If using modules for Python/OpenMPI/ORCA etc. then that all should be loaded here.

    # Load or set Python environment here:
    # e.g. module load python37  or:
    export PATH=/path/to/python/bin:$PATH
    # If using Conda, activate desired Conda environment.
    # May have to add conda bin directory to $PATH first.
    #conda activate ashpy37



    #Add path to Julia
    export PATH=/path/to/julia/bin:$PATH

    #Put ASH in PYTHONPATH and LD_LIBRARY_PATH
    export PYTHONPATH=/path/to/ash:$PYTHONpath
    export LD_LIBRARY_PATH=/path/to/ash:/path/to/ash/lib:$LD_LIBRARY_PATH

    #Print out environment variables for debuggin.
    echo "PATH is $PATH"
    echo "PYTHONPATH is $PYTHONPATH"
    echo "LD_LIBRARY_PATH is $LD_LIBRARY_PATH"
    echo ""
    echo "Running Ash  job"

    #Put ORCA in PATH and LD_LIBRARY_PATH
    export PATH=/path/to/orca:$PATH
    export LD_LIBRARY_PATH=/path/to/orca:$LD_LIBRARY_PATH

    #OpenMPI path for ORCA
    export PATH=/opt/openmpi-2.1.5/bin:$PATH
    export LD_LIBRARY_PATH=/opt/openmpi-2.1.5/lib:$LD_LIBRARY_PATH


    #Start Ash job from scratch dir.  Output file is written directly to submit directory
    export PYTHONUNBUFFERED=1
    python-jl $job.py >>& $SLURM_SUBMIT_DIR/$outputname

    # Ash has finished. Now copy important stuff back.
    outputdir=$SLURM_SUBMIT_DIR/${job}_${SLURM_JOB_ID}
    cp -r $tdir $outputdir
    #mkdir $outputdir
    #cp -r $tdir/*xyz $outputdir
    #cp -r $tdir/*txt $outputdir
    #cp -r $tdir/*xtl $outputdir
    #cp -r $tdir/*charges $outputdir
    #cp -r $tdir/orca*inp $outputdir
    #cp -r $tdir/orca*out $outputdir
    #cp -r $tdir/*.ygg $outputdir
    #cp -r $tdir/*.ff $outputdir
    #cp -r $tdir/*.info $outputdir

    # Removing scratch folder
    rm -rf $tdir





