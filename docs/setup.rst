.. raw:: html

     <style> .red {color:#aa0060; font-weight:bold; font-size:16px} </style>

.. role:: red

Setup
======================================
To install/setup ASH you need to download the code from the `Github <https://github.com/RagnarB83/ash>`_ repository.

ASH is 99% Python with 1 % Julia.
A Python3 distribution (version >3.6 or higher) is required and you need to be able to install Python packages via conda or pip.

Scientific Python distributions like Anaconda/miniconda can be convenient: https://www.anaconda.com/distribution/
Anaconda comes with Numpy, SciPy, Matplotlib.

Some functionality (primarily the molecular crystal QM/MM part) require a Julia installation (as the Python routines will be too slow).
Future versions may make the Julia interface a requirement.

Strict dependencies:

* `Python version 3.6 <https://www.python.org>`_ or higher
* `Numpy <https://numpy.org>`_ library.


Strongly recommended (necessary for some ASH functionality):

* `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ (Python package via pip). Required for geometry optimizations.
* `OpenMM <http://openmm.org>`_ version 7.6 or later. Required for most MM and MD functionality in ASH.
* `Julia 1.7 <https://julialang.org/downloads>`_ installation for fast routines in MolCrys QM/MM
* Python-Julia library: `PythonCall <https://cjdoris.github.io/PythonCall.jl/stable/pycall/>`_ or `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_

Recommended external codes (most ASH examples will use these):

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


##############################################
A. Download ASH
##############################################
Clone or download the archive from `Github <https://github.com/RagnarB83/ash>`_ containing ASH and put the directory (named ash) in some appropriate location.
The ash directory contains the Python source code files, e.g. ash_header.py and several directories.


##################################################
B. Installation and Configuration
##################################################

***************************************************************
B0. The lazy/impatient way to set up ASH (easy but incomplete)
***************************************************************

If you are impatient and want to get ASH going immediately without all features enabled.

1. Download the ASH repository from `Github <https://github.com/RagnarB83/ash>`_
2. Make sure there is a python3 interpreter on your system (with numpy)
3. Set: export PYTHONPATH=/path/to/ash:$PYTHONPATH   (where /path/to/ash is the directory containing the source-code directory (called "ash"))
4. Test ASH by launching: **python3**  and then do: from ash import *        (this will fail if PYTHONPATH is set incorrectly)

.. note:: ASH will complain when you try to use features that require additional installations (e.g. geometric, julia, OpenMM etc)

.. warning:: If your installation path is e.g. /home/ragnar/ASH-program/ash  (where the ash directory contains bunch of files including ash_header.py) then you should set PYTHONPATH like this: export PYTHONPATH=/home/ragnar/ASH-program:$PYTHONPATH   but NOT:  export PYTHONPATH=/home/ragnar/ASH-program/ash:$PYTHONPATH


*****************************************************
B1. Semi-Automatic Miniconda setup (easiest)
*****************************************************

This is the recommended way. Required if you intend to do MM or QM/MM using the OpenMM package (as OpenMM has to be installed via conda/miniconda).

1. Install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.  Install it in a location where your user has access (e.g. your home-directory)
2. Create new environment (not required but recommended): "conda create \-\-name ASHv2"
3. Load environment: **conda activate ASH**
4. Change directory to ASH location 
5. Install all desired packages listed in: ASH-packages.sh (inside ASH source code directory) via conda or pip (conda is preferred).
6. Make sure the chosen Python-Julia interface works. PythonCall/JuliaCall is recommended, PyJulia is another option. See Section B3: Step 5a and 5b below for details.
7. Optional: Run: **julia julia-packages-setup.jl** to install some required Julia packages. Note: Julia dependency only required for molecular-crystal QM/MM.
8. Run: **bash conda_setup_ash.sh** # This creates new files: set_environment_ash.sh and python3_ash
9. Run: **source set_environment_ash.sh**  (this sets necessary PATHs and should probably be put in each user's .bash_profile, job-submission script etc.)
10. Test ASH by launching: **python3**  and then do: from ash import *

.. note:: if PyJulia interface was specifically installed (not recommended) you must use **python-jl** for ASH to correctly call Julia routines.

*****************************************************
B2. Semi-Automatic non-Conda setup
*****************************************************

This option is not recommended if you intend to use ASH and OpenMM for MM, MD and QM/MM functionality. 
This is because OpenMM is only easily installable via conda. See section B1 above.

This uses the nonconda_install_ash.sh script inside the ASH directory.
The script downloads and installs Python packages (numpy, geometric,pyjulia) as well as Julia and packages and creates a convenient script for setting up the ASH environment. It requires a working Python3 installation.

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

(Use only if semi-automatic approach B1 or B2 does not work)

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
An Anaconda/miniconda distribution can be used. Make sure the conda environment is loaded.

If you don't already have a suitable Python3 distribution, go to Step 2b.


**Step 2b. Install Python if required** 

*Option 1: Python3 via system package manager*

.. note:: This option might be preferred if installing on a cluster for multiple users.

Linux: Install Python3 via a Linux package manager (e.g. Centos: yum -y install python3, Ubuntu: apt install python3).
Installing via a package manager is prefereable than compiling from source (see python.org for compile options).
Mac OS X: TODO
Windows: TODO

Install numpy via pip:

.. code-block:: shell

    pip3 install numpy


Make sure that the Python3 that you have installed is in your PATH environment while finishing the setup process and when using ASH:

.. code-block:: shell

    export PATH=/path/to/python3/bin:$PATH



*Option 2: Anaconda/Miniconda Python3 setup*

Download `Anaconda Python3 package <https://www.anaconda.com/products/individual>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ and install in e.g. your user directory.
Follow Anaconda/Miniconda installation instructions. Install numpy unless already installed. It could be a good idea to create your own conda environment for ASH but this is not strictly necessary.


**Step 3.** To make ASH available to Python3, set the environment variables:

.. code-block:: shell
    
    export ASHPATH=/path/to/ash  # Change /path/to/ash to the actual ASH directory location on your machine
    export PYTHONPATH=$ASHPATH:$ASHPATH/lib:$PYTHONPATH
    export PATH=$ASHPATH:$PATH
    export LD_LIBRARY_PATH=$ASHPATH/lib:$LD_LIBRARY_PATH

where */path/to/ash* is the dir that contains the "ash" source-code diretory .
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

**PyJulia option:**

.. code-block:: shell

    pip3 install julia



* The Pyjulia executable, *python-jl* (available after pip3 install julia) must generally be used if Julia routines are called by ASH. It is needed for the PyJulia interface to work properly.

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

Optional installation of the `Psi4 <http://www.psicode.org/>`_ QM code (if you intend to use it), best done via Conda:

.. code-block:: shell

    conda install psi4 psi4-rt -c psi4


**Step 2.** Optional: Install OpenMM

For general MM, QM/MM and MD functionality in ASH,  the `OpenMM program <http://openmm.org>`_ must be available.
It can be installed using conda.

.. code-block:: shell

    conda install -c conda-forge openmm



#########################################
D. Test ASH
#########################################

Test if things work in general:
python3 /path/to/ash/ash/test_ash.py   #This runs a basic test job using the regular Python interpreter
python-jl /path/to/ash/ash/test_ash.py   #Only required when PyJulia is used



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

**python-jl (PyJulia) problem**

If you get an error message when launching python-jl (only when PyJulia has been installed) that looks like the following:

.. code-block:: text

    File "/home/bjornsson/ash/python3_ash", line 9, in <module>
    sys.exit(main())
    File "/home/bjornsson/.local/lib/python3.8/site-packages/julia/python_jl.py", line 114, in main
    execprog([julia, "-e", script_jl, "--"] + unused_args)
    ...
    FileNotFoundError: [Errno 2] No such file or directory

This means that the Python-Julia interface is not completely working.
Check the following:

1. Is Julia accessible from the shell?, i.e. does typing *julia* in the shell, launch the Julia interpreter ? If not then the PATH to Julia bin dir needs to set: export PATH=/path/to/julia/bin:$PATH
2. Something went wrong in the installation of Julia or PyJulia. Go through these steps again.
3. Make sure you are using the same Python environment you used when you installed things.
4. Set up PyCall for each Julia user environment (this updates ~/.julia dir)

