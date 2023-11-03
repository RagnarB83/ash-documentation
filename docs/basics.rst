
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
    #Create fragment from a multi-line Python string (Cartesian coordinates in Angstrom)
    fragstring="""
    H 0.0 0.0 0.0
    F 0.0 0.0 1.0
    """
    molecule = Fragment(coordsstring=fragstring, charge=0, mult=1)

    #Defining string containing ORCA Simple-input line
    orcasimpleinput="! BP86 def2-SVP tightscf"
    #Multi-line string containing ORCA block-input:
    orcablocks="""
    %scf maxiter 200
    end
    """"
    #Defining ORCATheory object
    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput, orcablocks=orcablocks, numcores=numcores)

    #Geometry Optimization using the Optimizer (requires geomeTRIC library)
    Optimizer(fragment=molecule, theory=ORCAcalc, coordsystem='tric')


The script above loads ASH, creates a new fragment from a string (see :doc:`coordinate-input` for other ways),
defines variables related to the ORCA-interface , creates an ORCA-theory object
(see :doc:`QM-interfaces` and :doc:`ORCA-interface` ), and runs a geometry optimization using the Optimizer optimizer function  (see :doc:`job-types` and :doc:`Geometry-optimization` ).


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
See e.g. **subash** script at the bottom of this page.

#####################################################
Interactive ASH in a REPL or iPython environment
#####################################################
It is also possible to run ASH within a read-eval-print-loop environment such as iPython.
This allows for interactive use of ASH. See video below for an example.

If ASH has been set up correctly (setting PYTHONPATH etc.) and iPython is available (pip install ipython), then ASH within iPython should be straightforward.
Make sure to use the iPython that uses the same Python environment as ASH.

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

Global settings are stored in  */path/to/ash/ash/settings_ash.py* and can in principle be modified. 
However, it is better to instead create a settings file called **ash_user_settings.ini** for your user 
in your home-directory that should look like below.
Here you can set whether to use ANSI colors in output, whether to print inputfile and logo at the top, timings at the bottom etc.
ASH will attempt to read this file on startup.

.. code-block:: text

    [Settings]
    scale = 1.0
    tol = 0.2
    use_ANSI_color = True
    print_input = True
    print_logo = True
    load_julia = True
    julia_library = pythoncall
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

| - less may require the -R flag: less -R outputfile. Or use the global setting: export LESS=-R
| - vim and emacs require plugins


#####################
Submitting job
#####################

For a more complicated job to be submitted to the queuing system on a cluster, 
it is best to have a submission script that configures the environment, 
copies inputfiles to local scratch and copies results back from scratch to submission directory.

The **subash**  script below can be used for this purpose (with appropriate modifications for your cluster)
The number of cores can be provided in the command-line (should match the number of cores requested for the Theory level).
Alternatively the **subash** script will automatically look for and read a numcores variable in the ASH inputfile (if present). 
Make sure to have a line containing: "numcores=X" in the Python script (as in the ashtest.py example above).

Using **subash**:

.. code-block:: text

    subash input.py   (where input.py is your ASH script)
    subash input.py -p 8  (where -p 8 indicates 8 cores to request from the queuing system)

**subash**:
This script is written for the SLURM queuing system (could be adapted easily to others).
To use, one needs to change the name of the queue, default walltime, name of local scratch location on each node etc. 
in the first lines of the script below.
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
    queue=LocalQ

    #Default walltime. Overridden by -w or --walltime
    walltime=900

    #Default memory per CPU
    memory_per_cpu=10

    #Default number of threads (OMP_NUM_THREADS, MKL_NUM_THREADS, OPENMM_CPU_THREADS)
    threads=1  #If threads=auto then threads are set to SLURM CPU cores

    #Default number of GPU slots
    number_of_gpus=0
    gpu_memory=1
    #-g gpu_mem:4G

    # Path to local scratch on nodes (Script will create temporary user and job-directory here)
    #8TB HDD scratch: /data-hdd/SCRATCH
    #1 TB SSD scratch. /data-nvme/SCRATCH
    SCRATCHLOCATION=/data-hdd/SCRATCH

    #Path to bash file containing PATH and LD_LIBRARY_PATH definitions
    #Should define environment for ASH and external QM programs
    ENVIRONMENTFILE=/homelocal/rb269145/scripts/set_environment_ash.sh

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
    echo "${yellow}Or: subash input.py -g 1:1GB      Submit with 1 GPU core and request 1 GB of memory.${normal}"
    echo "${yellow}Or: subash input.py -master      Submit using ASH masterbranch.${normal}"
    echo "${yellow}Or: subash input.py -new      Submit using ASH newbranch.${normal}"
    echo "${yellow}Or: subash input.py -m X      Memory setting per CPU (in GB): .${normal}"
    echo "${yellow}Or: subash input.py -t T      Submit using T threads (for OpenMP/OpenMM/MKL).${normal}"
    echo "${yellow}Or: subash input.py -s /path/to/scratchdir      Submit using specific scratchdir .${normal}"
    echo "${yellow}Or: subash input.py -mw            Submit multi-Python job (multiple walkers) .${normal}"
    echo "${yellow}Or: subash input.py -n aar154      Submit to specific node (-n, --node).${normal}"
    echo "${yellow}Or: subash input.py -q queuname    Submit to specific queue: .${normal}"
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

    #multiwalker default false
    multiwalker=false

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
        -t|--threads) #Number of threads
        threads="$2"
        shift # past argument
        shift # past value
        ;;
        -m|--mempercpu) #Memory per core (GB)
        memory_per_cpu="$2"
        shift # past argument
        shift # past value
        ;;
        -g|--gpu) #Number of GPUcores and memory
        gpustuff="$2"
        gpuoptions=(${gpustuff//:/ })
        number_of_gpus=${gpuoptions[0]}
        gpu_memory=${gpuoptions[1]}
        shift # past argument
        shift # past value
        ;;
        -mw|--multiwalker) #Multiwalker
        multiwalker=true
        shift # past argument
        #shift # past value
        ;;
        -w|--walltime) #Walltime
        walltime="$2"
        shift # past argument
        shift # past value
        ;;
        -n|--node) #Name of node
        specificnode="$2"
        shift # past argument
        shift # past value
        ;;
        -s|--scratchdir) #Name of scratchdir
        SCRATCHLOCATION="$2"
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

    #Memory setting
    echo "Memory setting per CPU: $memory_per_cpu GB"
    slurm_mem_line="#SBATCH --mem-per-cpu=${memory_per_cpu}G"

    #THREADS-setting. Applies to programs using multithreading that needs to be controlled by:
    #OMP_NUM_THREADS, MKL_NUM_THREADS or OPENMM_CPU_THREADS variables
    #If threads is auto
    if [[ $threads == auto ]]
    then
    echo "Threads-setting is auto. Setting threads equal to numcores: $numcores"
    threads=$numcores
    else
    echo "Threads set to $threads"
    fi

    #Possible request of GPUs
    if [[ $number_of_gpus != 0 ]]
    then
    echo "GPUs requested: $number_of_gpus"
    echo "GPU memory: $gpu_memory"
    slurm_gpu_line1="#SBATCH --gres=gpu:$number_of_gpus"
    slurm_gpu_line2="#SBATCH --mem-per-gpu $gpu_memory"
    else
    slurm_gpu_line1=""
    slurm_gpu_line2=""
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
    $slurm_gpu_line1
    $slurm_gpu_line2
    #SBATCH --output=%x.o%j
    #SBATCH --error=%x.o%j
    $slurm_mem_line

    export job=\$SLURM_JOB_NAME
    export job=\${job%%.*}

    #Outputname
    outputname="\$job.out"

    #Multiwalker option
    multiwalker=$multiwalker


    #NUM_CORES
    NUM_CORES=\$((SLURM_JOB_NUM_NODES*SLURM_CPUS_ON_NODE))

    #For OpenMM we set this variable
    export OPENMM_CPU_THREADS=$threads

    #Setting MKL_NUM_THREADS and OMP_NUM_THREADS to 1 by default.
    #Note: pyscf threading will not work without this.
    #OpenMM uses variable above instead.
    #Wfoverlap should be set by subprocess call. Should be possible to control all subprocess calls by embedding export OMP_NUM_THREADS into subprocess call
    #Python libraries are trickier
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export OMP_STACKSIZE=1G
    export OMP_MAX_ACTIVE_LEVELS=1

    echo "OPENMM_CPU_THREADS: \$OPENMM_CPU_THREADS"
    echo "MKL_NUM_THREADS: \$MKL_NUM_THREADS"
    echo "OMP_NUM_THREADS: \$OMP_NUM_THREADS"


    # Usage:
    #ulimit -u unlimited
    #limit stacksize unlimited

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
    cp \$SLURM_SUBMIT_DIR/*nat \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.chk \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.xtl \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.ff \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.ygg \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.pdb \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/*.info \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/POTENTIAL \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/BASIS_MOLOPT \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/qmatoms \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/hessatoms \$tdir/ 2>/dev/null
    cp \$SLURM_SUBMIT_DIR/Hessian* \$tdir/ 2>/dev/null
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

    echo "Node(s): \$SLURM_JOB_NODELIST"
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


    # Multiple walker ASH run (intended for multiwalker metadynamics primarily

    if [ "\$multiwalker" = true ]
    then
    echo "Multiwalker True! NUM_CORES: \$NUM_CORES"
    #Creating multiple subdir walkersim$i HILLS files are stored in $tdir
    for (( i=0; i<\$NUM_CORES; i++ ))
    do
        echo "Creating dir: walkersim\$i"  >> \$SLURM_SUBMIT_DIR/\$outputname
        mkdir walkersim\$i
        echo "Copying files to dir: walkersim\$i"  >> \$SLURM_SUBMIT_DIR/\$outputname
        cp * walkersim\$i/
        cd walkersim\$i
        echo "Entering dir: walkersim\$i"  >> \$SLURM_SUBMIT_DIR/\$outputname
        echo "Process launched : \$i"  >> \$SLURM_SUBMIT_DIR/\$outputname
        sleep 2
        python3 \$job.py >> \$SLURM_SUBMIT_DIR/\${job}_walker\${i}.out 2>&1 &
        declare P\$i=\$!
        cd ..
    done
    wait
    # \$P1 \$P2  #Does not matter?

    else
    #Regular job
    python3 \$job.py >> \$SLURM_SUBMIT_DIR/\$outputname 2>&1

    fi

    #Making sure to delete potentially massive  files before copying back
    rm -rf core.*  #If  segfaults
    rm -rf orca.*tmp* #ORCA tmp files from e.g. MDCI

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
    if [[  -z "$specificnode" ]]; then
    sbatch -J $file ash.job
    else
    #Submit to a specific node
    echo "Submitting to specific node: $specificnode"
    sbatch -J $file -w $specificnode ash.job
    fi
    echo "${cyan}ASH job: $file submitted using $numcores cores.${normal}"
    echo "Queue: $queue and walltime: $walltime"

    #Multiwalker
    if [[ "$multiwalker" == true ]]
    then
    echo "Multiwalker option chosen. ASH will create multiple dirs on scratch and submit $numcores jobs"
    echo "Make sure to adjust numcores inside ASH script!"
    fi
