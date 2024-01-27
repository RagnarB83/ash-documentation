.. raw:: html

     <style> .red {color:#aa0060; font-weight:bold; font-size:16px} </style>

.. role:: red

Setup
======================================

To install/setup ASH you need to download the code from the `Github <https://github.com/RagnarB83/ash>`_ repository or alternatively install via pip (see later).

ASH is 99% Python with 1 % Julia.
A Python3 distribution (version >= 3.7 or higher) is required and you need to be able to install Python packages via package managers such as mamba/conda or pip.

It is recommended to use a Miniforge/Minconda package manager to install Python and the required packages (OpenMM in particular)
Some functionality (primarily the molecular crystal QM/MM part) require a Julia installation (as the Python routines will be too slow).
Future versions may make the Julia interface a requirement.

Strict dependencies:

* `Python version 3.7 <https://www.python.org>`_ or higher
* `Numpy <https://numpy.org>`_ library.


Strongly recommended (necessary for some ASH functionality):

* `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ (Python package via pip). Required for geometry optimizations.
* `OpenMM <http://openmm.org>`_ version 7.6 or later. Required for most MM and MD functionality in ASH.

For Molecular crystal QM/MM functionality:
* `Julia 1.7 <https://julialang.org/downloads>`_ installation for fast routines in MolCrys QM/MM
* Python-Julia library: `PythonCall <https://cjdoris.github.io/PythonCall.jl/stable/pycall/>`_ or `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_

Recommended external QM codes (many ASH examples will use these) for basic semi-empirical and DFT calculations:

* `xTB program <https://xtb-docs.readthedocs.io/en/latest/>`_ The semi-empirical tightbinding DFT code by Grimme and coworkers.
* `ORCA program <https://orcaforum.kofo.mpg.de>`_ The popular free-for-academic-use HF,DFT,WFT program by Neese and coworkers.
* `PySCF program <https://pyscf.org>`_  The powerful open-source Python-based electronic structure library.

Useful libraries for specific functionality:

* `MDtraj <https://www.mdtraj.org>`_ MD trajectory analysis
* `Matplotlib <https://matplotlib.org>`_ library. Used to plot graphs/surfaces.
* `Scipy <https://www.scipy.org>`_ library. Used for interpolation routines when plotting surfaces and in molecular crystal QM/MM.
* `Plumed <https://www.plumed.org>`_ Plumed library
* `Parmed <https://parmed.github.io/ParmEd/html/index.html>`_ May be used by OpenMM interface.
* `MDAnalysis <https://www.mdanalysis.org>`_ MD trajectory analysis
* `ASE <https://wiki.fysik.dtu.dk/ase/>`_ Atomic Simulation Environment


##################################################
A. Python environment setup
##################################################

Installing `Miniforge <https://github.com/conda-forge/miniforge>`_ is recommended to handle Python and package installations (openMM in particular).
This installs both the `mamba <https://github.com/mamba-org/mamba>`_ and conda package manager and sets it up for use with `conda-forge <https://conda-forge.org>`_ 
collection of repositories.
Another option is: `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_  or some conda setup on a cluster.


Here is a simple setup for a Linux/Unix system that you can copy-paste into your shell:

.. code-block:: shell

    #Download Miniforge (mamba,conda)
    #Sites: https://github.com/conda-forge/miniforge and https://github.com/mamba-org/mamba
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    #Install without prompt (this accepts the license). Remove -b if you want to control where it installs miniforge3
    sh Miniforge3-Linux-x86_64.sh -b #creates a ~/miniforge3 directory
    #Initialize mamba/conda for your shell (if desired). 
    ~/miniforge3/bin/mamba init #This will add code to your shell-config file (.bashrc or equivalent)
    #Log out and log in again to activate the new shell-config
    #source-ing the shell-config file should also work:
    source ~/.bashrc
    #Create new environment for ASH (recommended):
    mamba create --name ASH #This creates a new environment called ASH
    mamba activate ASH #This activates ASH.
    #Install python and numpy in the ASH environment
    #Note: python 3.12 problematic at the moment
    mamba config --add channels conda-forge #Add conda-forge channel
    mamba install -c conda-forge python=3.11 numpy -y  #Installs approx. 50 MB of packages

.. note:: Note that mamba and conda are both installed by Miniforge. Mamba is a faster version of conda and is usually better behaved. They are pretty much interchangeable.

.. note:: If you do use a previous conda installation it is recommended to add the conda-forge channel like this :  conda config --add channels conda-forge

***************************************************************
B. The lazy/impatient way to set up ASH (easy but incomplete)
***************************************************************

If you are impatient and want to get ASH going immediately without all features enabled. 
Make sure you have a suitable Python interpreter available, ideally in a conda/mamba environment (see above).
For OpenMM functionality, you need to install OpenMM via conda/mamba. See below.

Option 1 : Install ASH via pip (recommended):
This will also add the Numpy and geometric dependency.

.. code-block:: shell

    #Install ASH using pip (default main branch). Approx. 390 MB
    python -m pip install git+https://github.com/RagnarB83/ash.git
    #Install the NEW (development) branch of ASH. Approx. 390 MB
    python -m pip install git+https://github.com/RagnarB83/ash.git@NEW

Option 2 (mostly if you want to help develop ASH): Download ASH from Github and set PYTHONPATH.

.. code-block:: shell

    #Download ASH from Github
    git clone https://github.com/RagnarB83/ash.git 
    #Do next: git checkout NEW if you want the development branch
    #Set PYTHONPATH to the ASH directory
    export PYTHONPATH=/path/to/ash:$PYTHONPATH   (where /path/to/ash is the directory containing README.md)


Test ASH immediately by launching: **python3**  (same python as used above!) and then do: 

.. code-block:: python

    from ash import *
    create_ash_env_file()  #This creates a file: set_environment_ash.sh

You can then do the following to activate the ASH environment for future shell sessions:

.. code-block:: shell

    source ~/set_environment_ash.sh 

.. note:: ASH will complain when you try to use features that require additional installations (e.g. OpenMM, julia, etc). You then have to install them via conda/mamba or pip. 
    Note that OpenMM requires a conda/mamba environment. See below.


See :doc:`basics` for information on how to use ASH, including how to submit ASH jobs to a cluster (e.g. using the **subash** submission script).

*****************************************************
C. Semi-Automatic Miniconda setup (recommended)
*****************************************************

This is the recommended way for a fully functioning ASH. 
Required if you intend to do MM or QM/MM using the OpenMM package (as OpenMM has to be installed via conda/mamba).

1. Install Miniforge or Miniconda (see section A above).  Install it in a location where your user has access (e.g. your home-directory)
2. Create new environment (recommended): **mamba create --name ASH** (you can also use conda)
3. Load environment: **mamba activate ASH** #IMPORTANT
4. python -m pip install git+https://github.com/RagnarB83/ash.git #This installs ASH in your environment
5. Install some of the desired packages listed in: `ASH-packages.sh <https://github.com/RagnarB83/ash/blob/master/ASH-packages.sh>`_ (inside ASH source code directory) via conda or pip.
   You can always come back to this step (just remember to do **mamba activate ASH** first).

Test ASH immediately (with **mamba activate ASH**  activated) by launching in the same shell session: **python3**  and then do: 

.. code-block:: python

    from ash import * #If you get an error here then ASH is not installed correctly
    create_ash_env_file()  #This creates a file: set_environment_ash.sh

The *~/set_environment_ash.sh* file created by the **create_ash_env_file** function above is a convenient way to activate the ASH environment for future shell sessions.
It can be sourced in your shell environment startup file (e.g. *.bashrc*, *.bash_profile* or *.zshrc* ) and in your job-submission script. 
It sets the necessary PATHs for ASH and Python to work without having to load the conda/mamba environment each time. 
It is recommended to add PATH and LD_LIBRARY_PATH definitions for various external packages (e.g. ORCA) to this file.

See :doc:`basics` for information on how to use ASH, including how to submit ASH jobs to a cluster (e.g. using the **subash** submission script).

.. note:: If you want to add packages (using mamba/conda or pip) to your ASH environment (i.e. go back to step 5 above), always make sure you have activated the ASH environment first: **mamba activate ASH**. Otherwise the packages will be added to your base environment instead.
    Do **mamba info --envs** to see your environments and which one is active.

Only if molecular crystal QM/MM feature is needed:

- Optional: Make sure the Python-Julia interface works (only needed for MolCrys QM/MM functionality). PythonCall/JuliaCall is recommended. See section F for problems.


#########################################
D. Install External Programs
#########################################

See `ASH-packages.sh <https://github.com/RagnarB83/ash/blob/master/ASH-packages.sh>`_  in ASH source code directory!

**Step 1.** Install desired QM program(s):

.. warning:: Don't try to install everything all at once. Chances are you only need a select few of the QM-programs.

Examples:

* `ORCA <https://orcaforum.kofo.mpg.de>`_ is a recommended QM code (flexible interface in ASH). See installation instructions on the `ORCA Input Library <https://sites.google.com/site/orcainputlibrary/setting-up-orca>`_. The path to ORCA needs to be in PATH and LD_LIBRARY_PATH of your shell and later your jobscript.
* `pySCF <http://www.pyscf.org>`_ 
* `xTB <https://xtb-docs.readthedocs.io/en/latest/>`_ 
* `psi4 <https://psicode.org>`_


Some of these QM-programs are packages installable via either pip or conda/mamba:

.. code-block:: shell

    #pySCF
    python -m pip install pyscf       #PySCF QM program: http://www.pyscf.org
    #xtb: semi-empirical QM
    mamba install -c conda-forge xtb 
    #Psi4
    mamba install -c psi4 psi4 #Psi4 QM program: https://psicode.org


#########################################
E. Test ASH
#########################################

Example ASH script to try out with an external QM code (geometry optimization of H2O using ORCA):

.. code-block:: shell

    python3 first-ash-job.py


first-ash-job.py:

.. code-block:: python

    from ash import *

    #Create H2O fragment
    coords="""
    O       -1.377626260      0.000000000     -1.740199718
    H       -1.377626260      0.759337000     -1.144156718
    H       -1.377626260     -0.759337000     -1.144156718
    """
    H2O=Fragment(coordsstring=coords, charge=0, mult=1)

    #Defining ORCA-related variables
    orcasimpleinput="! BP86 def2-SVP def2/J tightscf"
    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput)

    #Geometry optimization
    geomeTRICOptimizer(fragment=H2O, theory=ORCAcalc, coordsystem='tric')

This will only work if ORCA is available in the shell session. It is usually best to add PATH and LD_LIBRARY_PATH definitions for ORCA to your *~/set_environment_ash.sh* file.



#########################################
F. Installation problems
#########################################

**ASH library not found by Python interpreter**

Error message:

.. code-block:: text

    ModuleNotFoundError: No module named 'ash'

This means that you have not correctly told your Python environment where ASH exists. If you downloaded or cloned the code you need to either do:

.. code-block:: shell

    #Option 1: Set PYTHONPATH
    export PYTHONPATH=/path/to/ash:$PYTHONPATH 

    #Option 2: Locally install using pip
    cd /path/to/ash #Where the README.md file is located
    python -m pip install .

However, it is usually better to install directly from the repository:

.. code-block:: shell

    python -m pip install git+https://github.com/RagnarB83/ash.git


**Module numpy not found**

Error message:

.. code-block:: text

    ModuleNotFoundError: No module named 'numpy'

Your Python environment requires the numpy library to be installed. Install either via mamba/conda or pip.
Make sure that you have activated your ASH environment  (**mamba activate ASH** or **conda activate ash**).


**OpenMM or QM/MM or MD is not working in ASH**

For general MM, QM/MM and MD functionality in ASH,  the `OpenMM program <http://openmm.org>`_ must be available.
It can be installed using mamba/conda.

.. code-block:: shell

    mamba install -c conda-forge openmm
    #or :
    conda install -c conda-forge openmm

**Julia-Python interface not working**

ASH requires a Python-Julia library in order to enable communication between Python and Julia.
The recommended option  is: `PythonCall/julicall <https://cjdoris.github.io/PythonCall.jl/stable/pycall/>`_ 

It is best to have PythonCall handle the Julia installation.

.. code-block:: shell

    python -m pip install juliacall
    
Once juliacall is installed, check that it is working correctly by: 

1. Launch python3 interactive session : 

.. code-block:: shell

    python3 # in shell

2. Run in python3 session: 

.. code-block:: python3

    import juliacall   #This will try to import the PythonCall/Juliacall interface, will check for Julia availability etc. 
    #This may take a while. Once done:
    juliacall.Main.sin(34.5) #This will call the Julia sin function.

If no errors then things should be good to go for ASH.