
Basic usage
======================================

This page shows how to use ASH once you have correctly set up ASH as a Python library.
See :doc:`setup` for information on installing ASH.

#####################
Input structure
#####################

You create a Python3 script (e.g. called ashtest.py) and import the ASH library:

.. code-block:: python

    from ash import *   # This will import the most important ASH functionality into your namespace
    #or
    import ash   # If you use this option you will have to add the "ash." prefix in front of ASH functions/classes.


ASH functionality can only be imported if the ASH source dir is in the PYTHONPATH.
Make sure you have already set in the shell (part of Setup):

.. code-block:: shell

    export PYTHONPATH=/path/to/ash_dir:$PYTHONPATH


You then have the freedom of writing a Python script in whatever way you prefer but taking the advantage
of ASH functionality. Typically you would first create one (or more) molecule fragments, then define a theory
object and then call a specific job-module (singlepoint, an optimizer, numerical-frequencies etc.).
See  :doc:`basic-examples` for examples.

#####################
Example script
#####################

Here is a basic ASH Python script, e.g. named: ashtest.py

.. code-block:: python

    from ash import *

    #Defining a numcores variable with the number of available CPU cores. Will be given to ORCA object.
    numcores=4
    #Create fragment from a multi-line Python string.
    fragstring="""
    H 0.0 0.0 0.0
    F 0.0 0.0 1.0
    """
    molecule = Fragment(coordsstring=fragstring, charge=0, mult=1)

    #Defining ORCA-related variables
    orcasimpleinput="! BP86 def2-SVP tightscf"
    orcablocks="%scf maxiter 200 end"

    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput, orcablocks=orcablocks, numcores=numcores)

    #Geometry Optimization using geomeTRIC
    Optimizer(fragment=molecule, theory=ORCAcalc, coordsystem='tric')


The script above loads ASH, creates a new fragment from a string (see :doc:`coordinate-input` for other ways),
defines variables related to the ORCA-interface , creates an ORCA-theory object
(see :doc:`QM-interfaces` and :doc:`ORCA-interface` ), and runs a geometry optimization using the geomeTRICOptimizer optimizer function  (see :doc:`job-types` and :doc:`Geometry-optimization` ).


######################################
Running script directly in the shell
######################################

For a very simple short job we can just run the script directly

.. code-block:: shell

    python3 ashtest.py

The output will be written to standard output (i.e. your shell). 
To save the output it is better to redirect the output to a file.

.. code-block:: shell

    python3 ashtest.py >& ashtest.out

For a really long job you would typically submit a jobscript to a queuing system instead.

#####################################################
Interactive ASH in a REPL or iPython environment
#####################################################
It is also possible to run ASH within a read-eval-print-loop environment such as iPython.
This allows for interactive use of ASH. See video below for an example.

If ASH has been set up correctly (PYTHONPATH etc.) and iPython is available (pip install ipython), then ASH within iPython should be straightforward.
Make sure to use the iPython that uses the same Python as ASH.

.. raw:: html

    <div align=center>
   <script id="asciicast-MUrhNGhDx9mAjdqomBppIGWsI" src="https://asciinema.org/a/MUrhNGhDx9mAjdqomBppIGWsI.js" async></script>
    </div>


############################
ASH with Julia support
############################

Some ASH functionality (primarily the molecular crystal QM/MM code) requires a working Python-Julia interface as some of the Python routines are too slow.
ASH can use both `PythonCall <https://cjdoris.github.io/PythonCall.jl/stable/pycall/>`_ and `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_
The PythonCall/juliacall library is recommended.

#####################
ASH settings
#####################

Global settings are stored in  */path/to/ash/ash/settings_ash.py* and can in principle be modified. However, it is better to instead create a settings file called **ash_user_settings.ini** for your user in your home-directory that should look like below.
Here you can set whether to use ANSI colors in output, whether to print inputfile and logo, timings etc.
ASH will attempt to read this file on startup.

.. code-block:: text

    [Settings]
    scale = 1.0
    tol = 0.2
    use_ANSI_color = True
    print_input = True
    print_logo = True
    load_julia = True
    julia_library = pythoncall #pyjulia is another option
    debugflag = False
    print_exit_footer = True
    print_full_timings = True
    nonbondedMM_code = julia
    connectivity_code = julia

.. warning:: The file ~/ash_user_settings.ini should not contain ' or "" symbols when defining strings.

In addition to options above it is also possible to specify the paths to various external codes.
If these paths are set in the settings file, one can avoid defining them in the inputfiles.


.. code-block:: text

    [Settings]
    orcadir = /path/to/orcadir
    daltondir = /path/to/daltondir
    xtbdir = /path/to/xtbdir
    psi4dir = /path/to/psi4dir
    cfourdir = /path/to/cfourdir
    crestdir = /path/to/crestdir




###############################
Use of colors in ASH output
###############################

ASH can display ANSI colors in output if  use_ANSI_color = True   is used in the settings file (see above). 
This makes the output more readable.

Note, however, that colors will only display properly if using a text reader that supports it:

| - less may require the -R flag: less -R outputfile. Or setting: export LESS=-R
| - vim and emacs require plugins


#####################
Submitting job
#####################

For a more complicated job to be submitted to the queuing system on a cluster, one would like to have a submission
script that also configures the environment, copies inputfiles to local scratch and results from scratch to back etc.

The **subash**  script below can be used for this purpose.
The number of cores can be provided in the command-line (should match the number of cores requested for the Theory level).
Alternatively this script can read the numcores variable in the inputfile.py (if present). 
Make sure to have a line containing: "numcores=X" in the Python script (as in the ashtest.py example above).

Using **subash**:

.. code-block:: text

    subash input.py   (where input.py is your ASH script)
    subash input.py -p 8  (where -p 8 indicates 8 cores to request from the queuing system)

**subash**:
This script is written for the SLURM queuing system.
To use, one needs to change the name of the queue, default walltime, name of local scratch location on each node in the first lines of the script below.
Also one needs provide the path to an environment-file that configures the ASH environment (PYTHONPATH, PATH, LD_LIBRARY_PATH etc.).
If you set up ASH using the Conda instructions `Conda <https://ash.readthedocs.io/en/latest/setup.html#b1-semi-automatic-miniconda-setup-easiest>`_
you should have a file in the ASH directory: /path/to/ash/set_environment_ash.sh

.. code-block:: text

    #!/bin/bash
    #subash: Submission script for ASH
    #Usage: subash ash_script.py

    #######################################
    #CLUSTER SETTINGS (to be modified by user)

    #Default name of cluster queue to submit to. Can be verridden by -q or --queue option
    queue=default1

    #Default walltime. Overridden by -w or --walltime
    walltime=8760

    # Path to local scratch on nodes (Script will create temporary user and job-directory here)
    SCRATCHLOCATION=/nobackup

    #Path to bash file containing PATH and LD_LIBRARY_PATH definitions
    #Should define environment for ASH and external QM programs
    ENVIRONMENTFILE=/home/rb269145/ASH/NEW/ash/set_environment_ash.sh

    #Default ASH branch to use
    ashbranch="new"

    #######################################
    # End of user modifications (hopefully)
    #######################################


    #Colors
    green=`tput setaf 2`
    yellow=`tput setaf 3`
    normal=`tput sgr0`
    cyan=`tput setaf 6`

    print_usage () {
    echo "${green}subash${normal}"
    echo "${yellow}Usage: subash input.py      Dir should contain .py Python script.${normal}"
    echo "${yellow}Or: subash input.py -p 8      Submit with 8 cores.${normal}"
    echo "${yellow}Or: subash input.py -master      Submit using ASH masterbranch.${normal}"
    echo "${yellow}Or: subash input.py -new      Submit using ASH newbranch.${normal}"
    exit
    }

    arguments=$@
    argument_first=$1
    file=$argument_first
    argumentnum=$#
    #echo "Arguments provided : $arguments"

    #If positional argument not .py then exit
    if [[ $argument_first != *".py"* ]]; then
    echo "No .py file provided. Exiting..."
    echo
    print_usage
    fi

    #Go through arguments
    while [[ $# -gt 0 ]]
    do
    key="$1"
    case $key in
        -branch|--branch)
        ashbranch="$2"
        shift # past argument
        ;;
        -p|--procs|--cores|--numcores) #Number of cores
        numcores="$2"
        shift # past argument
        shift # past value
        ;;
        -q|--queue) #Name of queue
        queue="$2"
        shift # past argument
        shift # past value
        ;;
        -w|--walltime) #Name of queue
        walltime="$2"
        shift # past argument
        shift # past value
        ;;
        --default)
        DEFAULT=YES
        ;;
        *)    #
        shift # past argument
    esac
    done

    #Now checking if numcores are defined
    if [[ $numcores == "" ]]
    then
    #Grabbing numcores from input-file.py if not using -p flag
    echo "Numcores not provided (-p option). Trying to grab cores from Python script."
    var=$(grep '^numcores' $file)
    NPROC=$(echo $var | awk -F'=' '{print $NF}')
    numcores=$(echo $NPROC | sed -e 's/^[[:space:]]*//')
    if ((${#numcores} == 0))
    then
        echo "No numcores variable defined Python script. Exiting..."
        exit
    fi
    fi

    ######################
    #Job-script creation
    ######################
    rm -rf ash.job
    cat <<EOT >> ash.job
    #!/bin/bash

    #SBATCH -N 1
    #SBATCH --tasks-per-node=$numcores
    #SBATCH --time=$walltime:00:00
    #SBATCH -p $queue
    #SBATCH --output=%x.o%j
    #SBATCH --error=%x.o%j

    export job=\$SLURM_JOB_NAME
    export job=\${job%%.*}

    #Outputname
    outputname="\$job.out"

    #NUM_CORES
    NUM_CORES=\$((SLURM_JOB_NUM_NODES*SLURM_CPUS_ON_NODE))

    #For OpenMM
    export OPENMM_CPU_THREADS=\$NUM_CORES

    # Usage:
    #ulimit -u unlimited
    #limit stacksize unlimited

    #Necessary?
    #setopt EXTENDED_GLOB
    #setopt NULL_GLOB
    export MKL_NUM_THREADS=\$NUM_CORES
    export OMP_NUM_THREADS=\$NUM_CORES
    export OMP_STACKSIZE=1G
    export OMP_MAX_ACTIVE_LEVELS=1

    echo "OPENMM_CPU_THREADS: \$OPENMM_CPU_THREADS"
    echo "MKL_NUM_THREADS: \$MKL_NUM_THREADS"
    echo "OMP_NUM_THREADS: \$OMP_NUM_THREADS"


    #Create scratch
    scratchlocation=$SCRATCHLOCATION
    echo "scratchlocation: \$scratchlocation"
    #Checking if scratch drive exists
    if [ ! -d \$scratchlocation ]
    then
    echo "Problem with scratch directory location: \$scratchlocation"
    echo "Is scratchlocation in subash script set correctly ?"
    echo "Exiting"
    exit
    fi

    #Creating user-directory on scratch if not available
    if [ ! -d \$scratchlocation/\$USER ]
    then
    mkdir -p \$scratchlocation/\$USER
    fi
    #Creating temporary dir on scratch
    tdir=\$(mktemp -d \$scratchlocation/\$USER/ashjob__\$SLURM_JOB_ID-XXXX)
    echo "Creating temporary tdir : \$tdir"

    #Checking if directory exists
    if [ -z \$tdir ]
    then
    echo "tdir variable empty: \$tdir"
    echo "Problem creating temporary dir: \$scratchlocation/\$USER/ashjob__\$SLURM_JOB_ID-XXXX"
    echo "Is scratch-disk (\$scratchlocation) writeable on node: \$SLURM_JOB_NODELIST  ?"
    echo "Exiting"
    exit
    fi

    #Checking if tdir exists
    if [ ! -d \$tdir ]
    then
    echo "Problem creating temporary dir: \$scratchlocation/\$USER/ashjob__\$SLURM_JOB_ID-XXXX"
    echo "Is scratch-disk (\$scratchlocation) writeable on node: \$SLURM_JOB_NODELIST  ?"
    echo "Exiting"
    exit
    fi
    chmod +xr \$tdir
    echo "tdir: \$tdir"

    cp \$SLURM_SUBMIT_DIR/*.py \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.cif \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.xyz \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.c \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.gbw \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.xtl \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.ff \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.ygg \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.pdb \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.info \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/qmatoms \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/hessatoms \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/act* \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.xml \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.txt \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.rtf \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.prm \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.gro \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.psf \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.rst7 \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.top \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.itp \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*prmtop \$tdir/ 2>/dev/null

    # cd to scratch
    echo "Entering scratchdir: \$tdir"
    cd \$tdir
    header=\$(df -h | grep Filesy)
    scratchsize=\$(df -h | grep \$scratchlocation)

    # Copy job and node info to beginning of outputfile
    echo "Starting job in scratch dir: \$tdir" > \$SLURM_SUBMIT_DIR/\$outputname
    echo "Job execution start: \$(date)" >> \$SLURM_SUBMIT_DIR/\$outputname
    echo "Shared library path: \$LD_LIBRARY_PATH" >> \$SLURM_SUBMIT_DIR/\$outputname
    echo "Slurm Job ID is: \${SLURM_JOB_ID}" >> \$SLURM_SUBMIT_DIR/\$outputname
    echo "Slurm Job name is: \${SLURM_JOB_NAME}" >> \$SLURM_SUBMIT_DIR/\$outputname
    echo "Nodes: \$SLURM_JOB_NODELIST" >> \$SLURM_SUBMIT_DIR/\$outputname
    echo "Scratch size before job:" >> \$SLURM_SUBMIT_DIR/\$outputname
    echo "\$header" >> \$SLURM_SUBMIT_DIR/\$outputname
    echo "\$scratchsize" >> \$SLURM_SUBMIT_DIR/\$outputname

    #ASH environment
    #This activates the correct Python, Julia, ASH environment
    source $ENVIRONMENTFILE

    echo "PATH is \$PATH"
    echo "LD_LIBRARY_PATH is \$LD_LIBRARY_PATH"
    export OMPI_MCA_btl=vader,self
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
    echo "Running ASH job"

    #Start ASH job from scratch dir.  Output file is written directly to submit directory
    export PYTHONUNBUFFERED=1
    python3 \$job.py >> \$SLURM_SUBMIT_DIR/\$outputname 2>&1

    header=\$(df -h | grep Filesy)
    echo "header: \$header"
    scratchsize=\$(df -h | grep \$scratchlocation)
    echo "Scratch size after job: \$scratchsize"

    # Ash has finished. Now copy important stuff back.
    outputdir=\$SLURM_SUBMIT_DIR/\${job}_\${SLURM_JOB_ID}
    cp -r \$tdir \$outputdir

    # Removing scratch folder
    rm -rf \$tdir

    EOT
    ######################
    #Submit job.
    sbatch -J $file ash.job
    echo "${cyan}ASH job: $file submitted using $numcores cores.${normal}"
    echo "Queue: $queue and walltime: $walltime"
