==========================
Basic usage
==========================

#####################
Input structure
#####################
You create a Python3 script (e.g. called ashtest.py) and import the Ash functionality:

.. code-block:: python

    from ash import *   # This will import the most important ASH functionality into your namespace
    #or
    import ash   # If you use this option you will have to add the "ash." prefix in front of ASH functions/classes.


Ash functionality can only be imported if the Ash source dir is in the PYTHONPATH.
Make sure you have already set in the shell (part of Setup):

.. code-block:: shell

    export PYTHONPATH=/path/to/ash_dir:$PYTHONPATH


Global settings are stored in your *ash-dir/settings_ash.py* and can in principle be modified. However, it is better to instead create a settings file, **ash_user_settings.ini** for your user in your home-directory that should look like below.
Here you can set whether to use ANSI colors in output, whether to print inputfile and logo, timings etc.

.. code-block:: shell

    [Settings]
    scale = 1.0
    tol = 0.2
    use_ANSI_color = True
    print_input = True
    print_logo = True
    load_julia = True
    debugflag = False
    print_exit_footer = True
    print_full_timings = True
    nonbondedMM_code = "julia"
    connectivity_code = "julia"
    orcadir = '/path/to/orcadir'

In addition to options above it is also possible to specify the paths to various external codes.
If these paths are set in the settings file, one can avoid setting them in the inputfiles.

.. code-block:: shell

    [Settings]
    orcadir = '/path/to/orcadir'
    daltondir = '/path/to/daltondir'
    xtbdir = '/path/to/xtbdir'
    psi4dir = '/path/to/psi4dir'
    cfourdir = '/path/to/cfourdir'
    crestdir = '/path/to/crestdir'


You then have the freedom of writing a Python script in whatever way you prefer but taking the advantage
of ASH functionality. Typically you would first create one (or more) molecule fragments, then define a theory
object and then call a specific job-module (singlepoint, an optimizer, numerical-frequencies etc.).
See  :doc:`coordinate-input` for various ways of dealing with coordinates and fragments.

#####################
Example script
#####################

Here is a basic Ash Python script, e.g. named: ashtest.py

.. code-block:: python

    from ash import *

    #Setting numcores. Used by ORCA.
    numcores=4
    #Create fragment
    fragstring="""
    H 0.0 0.0 0.0
    F 0.0 0.0 1.0
    """
    molecule = Fragment(coordsstring=fragstring)

    #Defining ORCA-related variables
    orcadir='/opt/orca_4.2.1'
    orcasimpleinput="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    orcablocks="%scf maxiter 200 end"

    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1,
                                orcasimpleinput=orcasimpleinput, orcablocks=orcablocks, nprocs=numcores)

    #Geometry Optimization using geomeTRIC
    geomeTRICOptimizer(fragment=molecule, theory=ORCAcalc, coordsystem='tric')


The script above loads Ash, creates a new fragment from a string (see :doc:`coordinate-input` for other ways),
defines variables related to the ORCA-interface , creates an ORCA-theory object
(see :doc:`QM-interfaces`), and runs a geometry optimization using the SimpleOpt optimizer function  (see :doc:`job-types` for other better options).

########################
Running script directly
########################

For a simple job we can just run the script directly

.. code-block:: shell

    python3 ashtest.py
    #or (for full Python-Julia support)
    python3_ash ashtest.py

The output will be written to standard output (i.e. your shell). You can redirect the output to a file.

.. code-block:: shell

    python3 ashtest.py >& ashtest.out

#####################################################
Interactive ASH in a REPL or iPython environment
#####################################################
It is also possible to run ASH within a read-eval-print-loop environment such as iPython.
This allows for interactive use of ASH. See video below for an example.

If ASH has been set up correctly (PYTHONPATH etc.) then ASH within iPython should be straightforward.
Make sure to use the iPython that uses the same Python as ASH.

.. raw:: html

    <div align=center>
   <script id="asciicast-MUrhNGhDx9mAjdqomBppIGWsI" src="https://asciinema.org/a/MUrhNGhDx9mAjdqomBppIGWsI.js" async></script>
    </div>

#####################
Submitting job
#####################

For a more complicated job we would probably want to create a job-script that would handle various environmental variables,
dealing with local scratch, copy files back when done etc.
Here is an example SLURM jobscript. Remember to go through all the lines and change the various things like the path to
local scratch, set the correct PATH variables, load modules etc.

Use like this:

.. code-block:: shell

    sbatch -J ashtest.py jobscript.sh


where jobscript.sh is:

.. code-block:: shell

    #!/bin/zsh

    #SBATCH -N 1
    #SBATCH --tasks-per-node=1
    #SBATCH --time=8760:00:00
    #SBATCH -p compute
    #SBATCH --mem-per-cpu=3000

    #Use like this:
    #sbatch -J inputfile.py jobscript.sh

    export job=$SLURM_JOB_NAME
    export job=$(echo ${job%%.*})
    outputname="$job.out"

    #Controlling threading
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

    #Python and ASH environment

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
    export PYTHONPATH=/path/to/ash:$PYTHONPATH
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
    python3_ash $job.py >> $SLURM_SUBMIT_DIR/$outputname 2>&1

    # Ash has finished. Now copy important stuff back.
    outputdir=$SLURM_SUBMIT_DIR/${job}_${SLURM_JOB_ID}
    cp -r $tdir $outputdir

    # Removing scratch folder
    rm -rf $tdir

For even more convenient job-submissions one can utilize a **subash** wrapper script that copies the jobscript.sh file (above)
to the current directory, modifies the number of cores requested and then submits.
The number of cores can be provided in the command-line (should match the number of cores requested in the ASH Python script, e.g. as in ashtest.py above)
or alternatively it can read the numcores variable in ashtest.py (if present). For the latter: make sure to have a line containing:
"numcores=X"
in the Python script (as in ashtest.py above).
Make sure to change path_to_jobscript variable in line 5.

.. code-block:: shell

    subash ashtest.py
    # or:
    subash ashtest.py -p 8  #for requesting an 8-core job.


.. code-block:: shell

    #!/bin/zsh
    #subash
    #Wrapper script for ASH job-script

    path_to_jobscript=/home/bjornsson/jobscripts/job-ash.sh

    green=`tput setaf 2`
    yellow=`tput setaf 3`
    normal=`tput sgr0`
    cyan=`tput setaf 6`
    if [[ "$1" == "" ]]
    then
      echo "${green}subash${normal}"
      echo "${yellow}Usage: subash input.py      Dir should contain .py Python script.${normal}"
      echo "${yellow}Or: subash input.py -p 8      Submit with 8 cores.${normal}"
      exit
    fi

    export file=$1


    if [[ "$2" == "-p" ]]
    then
      export NPROC=$3
    else
      #Grabbing numcores from input-file.py if not using -p flag
      echo "No -p N provided. Grabbing cores from Python script (searches for line beginning with numcores= )"
      var=$(grep '^numcores' $file)
      export NPROC=$(echo $var | awk -F'=' '{print $2}')
      #export NPROC=$(grep -m 1 numcores $file | awk -F'=' '{print $2}')
      if ((${#NPROC} == 0))
      then
        echo "No numcores variable in Python script found. Exiting..."
        exit
      fi
    fi

    #Copying job-script to dir:
    cp $path_to_jobscript .
    #Note: jobscript should have tasks-per-node set to 1 for the sed substitution to work
    sed -i "s/#SBATCH --tasks-per-node=1/#SBATCH --tasks-per-node=$NPROC/g" job-ash.sh

    #Submit job.
    sbatch -J $file job-ash.sh
    echo "${cyan}ASH job submitted using $NPROC cores using file $file.$mult ${normal}"




