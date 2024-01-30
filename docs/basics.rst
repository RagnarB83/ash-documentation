
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


ASH can only be imported if the ASH source dir has either been installed (or manually set in the PYTHONPATH). See :doc:`setup`.

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
See discussion about the **`subash <https://github.com/RagnarB83/ash/raw/master/scripts/subash.sh>`_** submission script
at the bottom of this page.

#####################################################
Interactive ASH in a REPL or iPython environment
#####################################################

It is also possible to run ASH within a read-eval-print-loop environment such as iPython.
This allows for interactive use of ASH. See video below for an example.

If ASH has been set up correctly and iPython is available (pip install ipython), then ASH within iPython should be straightforward.
Make sure to use the iPython that uses the same Python environment as ASH.

.. raw:: html

    <div align=center>
   <script id="asciicast-MUrhNGhDx9mAjdqomBppIGWsI" src="https://asciinema.org/a/MUrhNGhDx9mAjdqomBppIGWsI.js" async></script>
    </div>

#####################################################
Interactive ASH in a Google Colab notebook
#####################################################


Try out the ASH basics in this Google Colab notebook: 

`ASH in Google Colab <https://colab.research.google.com/drive/11-FG7eTElCvcMNZiTIEXcdWjcR4YWRS-#scrollTo=ViPg1cGuck_a>`_

Requires a Google account.


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


##################################
Submitting ASH jobs on a cluster
##################################

Once you start using ASH in general you would of course want to submit jobs to a queuing system on a cluster. 
Then it is best to have a submission script that sets up the ASH-Python environment, 
copies inputfiles to local scratch on the computing node, runs the ASH-Python script and copies results back from scratch to submission directory.

A submission script,  `subash <https://github.com/RagnarB83/ash/blob/master/scripts/subash.sh>`_ ,  
has been created for this purpose which can be found in the scripts directory of the ASH source code.
It assumes the Slurm queuing system but should be fairly easily be adaptable to other systems.

Simply `download subash <https://github.com/RagnarB83/ash/blob/master/scripts/subash.sh>`_  or copy/paste the script into a file called e.g. subash.sh or subash:

.. code-block:: text

    #1. Download subash script from repository
    wget https://github.com/RagnarB83/ash/raw/master/scripts/subash.sh
    #2. Copy or move it to a suitable location on your cluster
    #3. Make it executable
    chmod +x subash
    #4. Make some changes to the script for your cluster (see below)

You will have to make the following changes to the **subash** script (lines 6-35):

- Change the name of the default queue to submit to (depends on cluster)
- Change the default walltime (depends on cluster and queue)
- Change the name of the local/global scratch directory on each node (depends on cluster and maybe queue)
- Change the path to the ASH environment file (set_environment_ash.sh). See more info below.

Using **subash**:

.. code-block:: shell

    subash input.py   # (where input.py is your ASH script). subash will figure out cores if numcores variable define in input.py
    subash input.py -p 8  # (where -p 8 indicates 8 cores to request from the queuing system)

When you use subash it will perform some checks such as whether the inputfile exists, whether you have provided CPU-core information etc.
It will create a Slurm submission file in the current directory (called ash.job) and then submit the job to the queuing system.
The number of CPU cores to use when submitting, can be provided as a flag to subash (should then match the number of cores requested for the Theory level in the inputfile).

However, our preferred way is to have the **subash** script automatically find out suitable number of CPU cores from the inputfile.
This will work if you have a line containing: "numcores=X" in the Python script (as in the ashtest.py example above and throughout the documentation) and then make sure 
that this variable is used to define the number of cores in the Theory object (see e.g. how the ORCAcalc object is defined in the example script at the top of this page). 
Many ASH examples in the documentation use this approach.

**subash** options:

.. code-block:: text

    subash
    No .py file provided. Exiting...

    subash
    Usage: subash input.py      Dir should contain .py Python script.
    Or: subash input.py -p 8      Submit with 8 cores.
    Or: subash input.py -g 1:1GB      Submit with 1 GPU core and request 1 GB of memory.
    Or: subash input.py -m X      Memory setting per CPU (in GB): .
    Or: subash input.py -t T      Submit using T threads (for OpenMP/OpenMM/MKL).
    Or: subash input.py -s /path/to/scratchdir      Submit using specific scratchdir .
    Or: subash input.py -mw            Submit multi-Python job (multiple walkers) .
    Or: subash input.py -n node54      Submit to specific node (-n, --node).
    Or: subash input.py -q queuename    Submit to specific queue.


**Configuring the ASH environmentfile:**

If you followed the installation instructions  for ASH (:doc:`setup`) you should have created an ASH environment file : ~/set_environment_ash.sh .
This file can be moved anywhere you like but make sure you modify the ENVIRONMENTFILE variable in the **subash** script to point to it.
It will be sourced (source $ENVIRONMENTFILE) right before python3 is launched (python3 input.py) in the job-submission file.

If you want to make sure ASH can automatically find external QM programs e.g. ORCA, MRCC, CFour, CP2K etc. as well as OpenMPI (used e.g. by ORCA) then it is best to add 
the necessary *PATH* and *LD_LIBRARY_PATH* definitions (or alternatively module load commands) to this ASH environment file.

Below is an example used on our cluster:

.. code-block:: text

    #!/bin/bash

    #ASH and Python environment
    ASHPATH=/data-hdd/PROGRAMS/ASH-PROGRAM/ash/ash
    python3path=/homelocal/rb269145/miniforge3/envs/ASHv1/bin

    #PYTHONPATH for finding ASH generally not recommended. commented out
    #export PYTHONPATH=$ASHPATH:\$ASHPATH/lib:$PYTHONPATH
    export PATH=$python3path:$PATH
    export LD_LIBRARY_PATH=$ASHPATH/lib:$LD_LIBRARY_PATH

    ulimit -s unlimited

    #OTHER PROGRAMS

    #MKL
    source /homelocal/rb269145/PROGRAMS/intel/oneapi/mkl/latest/env/vars.sh

    #Multiwfn
    Multiwfndir=/homelocal/rb269145/PROGRAMS/Multiwfn_3.8dev_bin_Linux_noGUI_jul2023

    #MRCC
    MRCCDIR=/data-hdd/PROGRAMS/MRCC-binaries

    #CFour
    CFOURDIR=/data-hdd/PROGRAMS/CFOUR/cvfour-v21dev-install-mkl-serial/bin

    #ORCA
    source /homelocal/rb269145/scripts/set_environment_orca5.sh

    #CP2K
    CP2KPATH=/homelocal/rb269145/PROGRAMS/CP2K/CP2K-2023.2-bin

    #Dice
    source /homelocal/rb269145/scripts/set_environment_dice_pyscf.sh

    #OPENMPI
    OPENMPIDIR=/data-hdd/PROGRAMS/openmpi-411-install

    #PATH and LD_LIBRARY_PATH modifications for program paths defined above
    export PATH=$OPENMPIDIR/bin:$MRCCDIR:$ORCADIR:$Multiwfndir:$CP2KPATH:$CFOURDIR:$PATH
    export LD_LIBRARY_PATH=$OPENMPIDIR/lib:$ORCADIR:$LD_LIBRARY_PATH

