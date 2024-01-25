.. raw:: html

     <style> .red {color:#aa0060; font-weight:bold; font-size:16px} </style>

.. role:: red

Setup
======================================

To install/setup ASH you need to download the code from the `Github <https://github.com/RagnarB83/ash>`_ repository or alternatively install via pip (see later).

ASH is 99% Python with 1 % Julia.
A Python3 distribution (version >3.6 or higher) is required and you need to be able to install Python packages via package managers such as mamba/conda or pip.

It is recommended to use a Miniforge/Minconda package manager to install Python and the required packages (OpenMM in particular)
Some functionality (primarily the molecular crystal QM/MM part) require a Julia installation (as the Python routines will be too slow).
Future versions may make the Julia interface a requirement.

Strict dependencies:

* `Python version 3.6 <https://www.python.org>`_ or higher
* `Numpy <https://numpy.org>`_ library.


Strongly recommended (necessary for some ASH functionality):

* `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ (Python package via pip). Required for geometry optimizations.
* `OpenMM <http://openmm.org>`_ version 7.6 or later. Required for most MM and MD functionality in ASH.

For Molecular crystal QM/MM functionality:
* `Julia 1.7 <https://julialang.org/downloads>`_ installation for fast routines in MolCrys QM/MM
* Python-Julia library: `PythonCall <https://cjdoris.github.io/PythonCall.jl/stable/pycall/>`_ or `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_

Recommended external QM codes (most ASH examples will use these):

* `xTB <https://xtb-docs.readthedocs.io/en/latest/>`_ The semi-empirical tightbinding DFT code by Grimme and coworkers.
* `ORCA <https://orcaforum.kofo.mpg.de>`_ The popular free-for-academic-use HF,DFT,WFT program by Neese and coworkers.


Useful libraries for specific functionality:

* `Matplotlib <https://matplotlib.org>`_ library. Used to plot graphs/surfaces.
* `Scipy <https://www.scipy.org>`_ library. Used for interpolation routines when plotting surfaces and in molecular crystal QM/MM.
* `Plumed <https://www.plumed.org>`_ Plumed library
* `Parmed <https://parmed.github.io/ParmEd/html/index.html>`_ May be used by OpenMM interface.
* `MDAnalysis <https://www.mdanalysis.org>`_ MD trajectory analysis
* `MDtraj <https://www.mdtraj.org>`_ MD trajectory analysis
* `ASE <https://wiki.fysik.dtu.dk/ase/>`_ Atomic Simulation Environment


##################################################
A. Python environment setup
##################################################

Installing `Miniforge <https://github.com/conda-forge/miniforge>`_ is recommended to handle Python and package installations (openMM in particular).
This installs both the `mamba <https://github.com/mamba-org/mamba>`_ and conda package manager and sets it up for use with `conda-forge <https://conda-forge.org>`_ 
collection of repositories.
Another option is: `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ .


Here is a simple workflow for a Linux/Unix system:

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
    mamba install python=3.11 numpy  #Installs approx. 50 MB of packages

.. note:: Note that mamba and conda are both installed by Miniforge. Mamba is a faster version of conda and is usually better behaved. They are pretty much interchangeable.

***************************************************************
B0. The lazy/impatient way to set up ASH (easy but incomplete)
***************************************************************

If you are impatient and want to get ASH going immediately without all features enabled. 
Make sure you have a suitable Python interpreter available, ideally in a conda/mamba environment (see above).
For OpenMM functionality, you need to install OpenMM via conda/mamba. See below.

Option1 : Install ASH via pip (recommended):
This will also add the Numpy and geometric dependency.

.. code-block:: shell

    #Install ASH using pip (default main branch). Approx. X MB of packages
    pip install git+https://github.com/RagnarB83/ash.git
    #Install the NEW (development) branch of ASH. Approx. X MB of packages
    pip install git+https://github.com/RagnarB83/ash.git@NEW

Option 2 (mostly if you want to help develop ASH): Download ASH from Github and set PYTHONPATH.

.. code-block:: shell

    #Download ASH from Github
    git clone https://github.com/RagnarB83/ash.git 
    #git checkout NEW if you want the development branch
    #Set PYTHONPATH to the ASH directory
    export PYTHONPATH=/path/to/ash:$PYTHONPATH   (where /path/to/ash is the directory containing README.md)


Test ASH immediately by launching: **python3**  and then do: 

.. code-block:: python

    from ash import *
    create_ash_env_file()  #This creates a file: set_environment_ash.sh

.. note:: ASH will complain when you try to use features that require additional installations (e.g. OpenMM, julia, etc). You then have to install them via conda/mamba or pip. 
    Note that OpenMM requires a conda/mamba environment. See below.


*****************************************************
B1. Semi-Automatic Miniconda setup (recommended)
*****************************************************

This is the recommended way for a fully functioning ASH. 
Required if you intend to do MM or QM/MM using the OpenMM package (as OpenMM has to be installed via conda/miniconda).

1. Install Miniforge or Miniconda (see A above).  Install it in a location where your user has access (e.g. your home-directory)
2. Create new environment (recommended): **mamba create --name ASH**
3. Load environment: **mamba activate ASH** #IMPORTANT
4. pip install git+https://github.com/RagnarB83/ash.git #This installs ASH in your environment
5. Install some of the desired packages listed in: `ASH-packages.sh <https://github.com/RagnarB83/ash/blob/master/ASH-packages.sh>`_ (inside ASH source code directory) via conda or pip. 

Test ASH immediately by launching in the same shell session: **python3**  and then do: 

.. code-block:: python

    from ash import * #If you get an error here then ASH is not installed correctly
    create_ash_env_file()  #This creates a file: set_environment_ash.sh

The set_environment_ash.sh file is a convenient way to activate the ASH environment in general. 
It can be sourced in your shell environment startup file (.bashrc, .bash_profile, .zshrc etc.) or in your jobscript. It sets the necessary PATHs for ASH to work
without having to load the conda/mamba environment.

If molecular crystal QM/MM feature is needed:

- Optional: Make sure the Python-Julia interface works (only needed for MolCrys QM/MM functionality). PythonCall/JuliaCall is recommended. See Section B3: Step 5a and 5b below for details.

*****************************************************
B2. Semi-Automatic non-Mamba/Conda setup
*****************************************************

This option is not recommended if you intend to use ASH and OpenMM for MM, MD and QM/MM functionality. 
This is because OpenMM is only easily installable via mamba/conda. See section B1 above.

This uses the nonconda_install_ash.sh script inside the ASH directory.
The script downloads and installs Python packages (numpy, geometric) as well as Julia and packages and creates a convenient script for setting up the ASH environment. It requires a working Python3 installation.

**Step 1.** Make sure the desired python3 is in your environment ('which python3' in the shell) or set path_to_python3_dir in the ./nonconda_install_ash.sh script to the Python3 installation you want to use. Script has a few possible settings in the beginning.
Note: You need to be able to install packages to this installation via pip 

**Step 2.** cd to ASH directory, make nonconda_install_ash.sh executable and run: 

- chmod +x ./nonconda_install_ash.sh
- ./nonconda_install_ash.sh

**Step 3.** If installation is successful:

- source ./set_environment_ash.sh    to activate ASH environment.


*****************************************************
B3. Manual
*****************************************************

(Use only if optons B0, B1 or B2 do not work)

**Step 1.** 

Make sure ASH has been downloaded and moved to some location where it will stay.
The location of the ASH directory will be referred to as /path/to/ash below (substitute /path/to/ash for the actual location on your machine).

**Step 2.** 

Check if a suitable Python3 installation is available (globally available or maybe via a module on your cluster). It needs to be relatively new (version 3.6 and above) contain Numpy and you will need to be able to install Python packages to it using the package manager pip. 

.. code-block:: shell

    #Check where python3 is:
    which python3
    #Check Python3 version
    python3 --version
    #Check that pip/pip3 is available (sometimes pip3 should be called instead of pip)
    which pip3  #Make sure the pip path is the same as python3 path)
    #Check that numpy is available inside the Python3 installation
    pip3 list | grep numpy


If you have a suitable Python3 with numpy then make sure it is loaded in your environment when using ASH.
A miniconda/miniforge distribution can be used. Make sure the conda/mamba environment is loaded.

If you don't already have a suitable Python3 distribution, go to Step 2b.


*Option 2: Miniforge/Miniconda Python3 setup*

Download `Miniforge <https://github.com/conda-forge/miniforge>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ and install in e.g. your user directory.
Follow Miniforge/Miniconda installation instructions. Install numpy unless already installed. It could be a good idea to create your own conda environment for ASH but this is not strictly necessary.


**Step 3.** To make ASH available to Python3, set the environment variables:

.. code-block:: shell
    
    export ASHPATH=/path/to/ash  # Change /path/to/ash to the actual ASH directory location on your machine
    export PYTHONPATH=$ASHPATH:$ASHPATH/lib:$PYTHONPATH
    export PATH=$ASHPATH:$PATH
    export LD_LIBRARY_PATH=$ASHPATH/lib:$LD_LIBRARY_PATH

where */path/to/ash* is the dir that contains README.md.
Put these environment definitions in your shell environment startup file e.g. .bashrc, .bash_profile or .zshrc.
This step will be necessary for each user on the cluster.

**Step 4.** Install the recommended Python packages via pip/pip3:

.. code-block:: shell

    pip3 install geometric   (geomeTRIC optimizer)

This requires you to be able ot install packages to your Python installation. It may also be possible to install Python packages
locally to your user's home directory by the "--user" option:  pip3 install geometric --user


**Step 5a.** Install Julia from the `Julia official site <https://julialang.org/downloads>`_.

Julia is necessary for some fast QM/MM functionality inside ASH (e.g. Molcrys). This step can be skipped if you won't be using the molecular crystal QM/MM feature.

 i) Download appropriate binaries from the official Julia website. Version 1.7 or higher. Extract archive.
 ii) Add Julia binaries to path: e.g. export PATH=/path/to/julia-1.7.1/bin:$PATH . Put this PATH definition in your shell startup file.
 iii) Run Julia using the ASH sourcefile julia-packages-setup.jl (inside ASH source directory) as input to download and install the  required Julia packages. Currently: PyCall, Hungarian, Distances

.. code-block:: shell

    julia julia-packages-setup.jl  #This launches the julia interpreter and requests installation of required Julia packages for ASH.

This will download and install required Julia packages.

.. note:: To avoid having to setup the Julia packages for each user on a computing cluster, one can specify a global Julia package-store-location: export JULIA_DEPOT_PATH=/path/to/julia-packages-dir  before running :  julia julia-packages-setup.jl


If there is an error like this: ERROR: SystemError: opening file "/path/to/.julia/registries/General/Registry.toml": No such file or directory
Then execute in shell: rm -rf ~/.julia/registries/General

**Step 5b.** Install Julia-Python interface

ASH requires a Python-Julia library in order to enable communication between Python and Julia.
The options are: `PythonCall <https://cjdoris.github.io/PythonCall.jl/stable/pycall/>`_ and `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_
ASH currently supports both but the newer PythonCall is currently recommended due to PyJulia currently requiring to call ASH with a modified Python interpreter (python-jl) due to static libpython issues.

:red:`Important:` Make sure the correct Python environment is active before proceeding. Check that the pip or pip3 executable is available and corresponds to the Python you want:

:red:`Important:` Make sure the Julia executable is in your PATH already.

.. code-block:: shell

    which pip
    which pip3

Then install using pip/pip3:

**PythonCall/JuliaCall option:**

.. code-block:: shell

    pip3 install juliacall

Once juliacall is installed, check that it is working correctly by: 

1. Launch python3 interactive session : 

.. code-block:: shell

    python3 # in shell

2. Run in python3 session: 

.. code-block:: python3

    import juliacall   #This will try to import the PythonCall/Juliacall interface, will check for Julia availability etc. 
    juliacall.Main.sin(34.5) #This will call the Julia sin function.

If no errors then things should be good to go for ASH.

* Make sure the correct Python3 environment is active. Otherwise ASH will not work.

* The regular Python3 executable, *python3*  can also be used to run ASH scripts and is recommended if you don't require ASH to launch Julia routines (molcrystal-QM/MM primarily). There may be warnings about the Python-Julia-interface not working. These warnings can be ignored . For large systems or when using QM/MM-Molcrys, this is not a good option, however, as very slow Python routines will be used for time-consuming steps.


#########################################
C. Install External Programs
#########################################

See also ASH-packages.sh in ASH source code directory!

**Step 1.** Install desired QM program(s):

* `ORCA <https://orcaforum.kofo.mpg.de>`_ is a recommended QM code (flexible interface in ASH). See installation instructions on the `ORCA Input Library <https://sites.google.com/site/orcainputlibrary/setting-up-orca>`_. The path to ORCA needs to be in PATH and LD_LIBRARY_PATH of your shell and later your jobscript.
* `xTB <https://xtb-docs.readthedocs.io>`_ needs to be in PATH and later your jobscript.


Optional Python packages to install via pip (depends on whether you will use the interfaces to PySCF and PyFrame):

* `PySCF <http://www.pyscf.org/>`_
* `PyFrame <https://gitlab.com/FraME-projects/PyFraME>`_:


.. code-block:: shell

    pip3 install pyscf       #PySCF QM program
    pip3 install pyframe     #polarizable embedding helper tool


**Step 2.** Optional: Install OpenMM

For general MM, QM/MM and MD functionality in ASH,  the `OpenMM program <http://openmm.org>`_ must be available.
It can be installed using mamba/conda.

.. code-block:: shell

    mamba install -c conda-forge openmm

#########################################
D. Test ASH
#########################################

Test if things work in general:
python3 /path/to/ash/ash/test_ash.py   #This runs a basic test job using the regular Python interpreter


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



#########################################
E. Installation problems
#########################################

**ASH library not found by Python interpreter**

Error message:

.. code-block:: text

    ModuleNotFoundError: No module named 'ash'

This means that you have not correctly told your Python environment where ASH exists. If you downloaded or cloned the code you need to either do:

.. code-block:: shell

    #Option 1: Set PYTHONPATH
    export PYTHONPATH=/path/to/ash:$PYTHONPATH 

    #Locally install using pip
    cd /path/to/ash #Where the README.md file is located
    pip install .

However, it is usually better to install directly from the repository:

.. code-block:: shell

    pip install git+https://github.com/RagnarB83/ash.git


**Module numpy not found**

Error message:

.. code-block:: text

    ModuleNotFoundError: No module named 'numpy'

Your Python environment requires the numpy library to be installed. Install either via mamba/conda or pip.
