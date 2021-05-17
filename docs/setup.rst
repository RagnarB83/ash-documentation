.. raw:: html

     <style> .red {color:#aa0060; font-weight:bold; font-size:16px} </style>

.. role:: red

Setup
======================================
Contact Ragnar if you want access to the code.

ASH is 99% Python with 1 % Julia.
A Python3 distribution (version >3.6 or higher) is required and as packages will have to be installed you need to be able to
install Python packages via pip.

Scientific Python distributions like Anacond/miniconda can be convenient: https://www.anaconda.com/distribution/
Anaconda comes with Numpy, SciPy, Matplotlib.

Treatment of large systems via QM/MM require a Julia installation (as the Python routines will be too slow).
Future versions will make Julia a requirement.

Strict dependencies:

* `Python version 3.6 <https://www.python.org>`_ >=
* `Numpy <https://numpy.org>`_ library.
* `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ (Python package via pip).

Strongly recommended (necessary for some parts):

* `Julia 1.6 <https://julialang.org/downloads>`_ installation for fast routines for large system treatment.
* `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_ installation (Python package via pip).

Useful:

* `Matplotlib <https://matplotlib.org>`_ library. Used to plot graphs/surfaces.
* `Scipy <https://www.scipy.org>`_ library. Used for interpolation routines when plotting surfaces.


##############################################
A. Download ASH
##############################################
Clone or download the archive containing ASH and put the directory (named ash) in your home directory or wherever you want it. The ash directory contains the Python source code files, named ash.py etc.


##################################################
B. Installation and Configuration
##################################################


*****************************************************
B1. Semi-Automatic
*****************************************************

(Experimental feature)

This uses the install_ash.sh script inside the ASH directory.
The script downloads and installs Python packages (numpy, geometric,pyjulia) as well as Julia and packages.
It requires a working Python3 installation.

**Step 1.** Set path_to_python3_dir in ./install_ash.sh script to the Python3 installation you want to use. Note: You need to be able to install packages to this installation via pip.

**Step 2.** cd to ASH directory and run: ./install_ash.sh

**Step 3.** Run: source ./set_environment_ash.sh

*****************************************************
B2. Manual
*****************************************************

**Step 1.** 

Clone or download an archive containing ASH and put the directory (named ash) in your home directory or wherever you want it. The ash directory contains the Python source code files, named ash.py etc.
The location of the ASH directory will be referred to as /path/to/ash below (substitute /path/to/ash for the actual loction on your machine).

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

If you don't already have a suitable Python3 distribution, go to Step 2b.


**Step 2b. Install Python if required** 

*Option 1: Python3 via system package manager*

.. note:: This option is preferred if installing on a cluster for multiple users.

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

where */path/to/ash* is the dir where all the ASH sourcefiles are (e.g. ash.py) .
Put these environment definitions in your shell environment startup file e.g. .bashrc, .bash_profile or .zshrc.
This step will be necessary for each user on the cluster.

**Step 4.** Install the recommended Python packages via pip/pip3:

.. code-block:: shell

    pip3 install geometric   (geomeTRIC optimizer)

**Step 5a.** Install Julia from the `Julia official site <https://julialang.org/downloads>`_.

Julia is necessary for some fast QM/MM functionality inside ASH (e.g. MolCrys). Step can be skipped if you won't be using QM/MM.

 i) Download appropriate binaries from the official Julia website. Version 1.6 or higher. Extract archive.
 ii) Add Julia binaries to path: e.g. export PATH=/path/to/julia-1.6.1/bin:$PATH . Put this PATH definition in your shell startup file.
 iii) Run Julia using the ASH sourcefile julia-packages-setup.jl (inside ASH source directory) as input to download and install the  required Julia packages. Currently: PyCall, Hungarian, Distances

.. code-block:: shell

    julia julia-packages-setup.jl  #This launches the julia interpreter and requests installation of required Julia packages for ASH.

This will download and install required Julia packages.

.. note:: To avoid having to setup the Julia packages for each user on a computing cluster, one can specify a global Julia package-store-location: export JULIA_DEPOT_PATH=/path/to/julia-packages-dir  before running :  julia julia-packages-setup.jl


If there is an error like this: ERROR: SystemError: opening file "/path/to/.julia/registries/General/Registry.toml": No such file or directory
Then execute in shell: rm -rf ~/.julia/registries/General

**Step 5b.** Install `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_


:red:`Important:` Make sure the correct Python environment is active before proceeding. Check that the pip or pip3 executable is available and corresponds to the Python you want:

.. code-block:: shell

    which pip
    which pip3

Then install using pip/pip3:

.. code-block:: shell

    pip3 install julia


**Step 6.** Activate python3_ash

Make the python3_ash executable (inside /path/to/ash): chmod +x /path/to/ash/python3_ash

* The ASH python3 executable, *python3_ash* should generally be used. It is required for the PyJulia interface to work
properly. The PyJulia interface is needed for treating large systems.

* Make sure the correct Python3 environment is active. Otherwise ASH will not work.

* The regular Python3 executable, *python3*  can also be used to run ASH scripts. There will, however, be a warning about the Python-Julia-interface not working. This warning can be ignored if fast Julia routines are not needed. For large systems or when using MolCrys, this is not a good option, however, as very slow Python routines will be used.


#########################################
C. Install External Programs
#########################################

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


**Step 2.** Optional: Install OpenMM (if needed)

For protein and explict solvation QM/MM in ASH, then the `OpenMM program <http://openmm.org>`_ is used as MM code.
It can be installed using conda.

.. code-block:: shell

    conda install -c omnia openmm



#########################################
D. Test ASH
#########################################

Test if things work in general:

python3_ash /path/to/ash/test_ash.py   #This runs a basic test job.



Example ASH script to try out with an external QM code (geometry optimization of H2O using ORCA):

.. code-block:: shell

    python3_ash first-ash-job.py


first-ash-job.py:

.. code-block:: python

    from ash import *

    #Create H2O fragment
    coords="""
    O       -1.377626260      0.000000000     -1.740199718
    H       -1.377626260      0.759337000     -1.144156718
    H       -1.377626260     -0.759337000     -1.144156718
    """
    H2Ofragment=Fragment(coordsstring=coords)
    #Defining ORCA-related variables
    orcadir='/opt/orca_4.2.1'
    orcasimpleinput="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    orcablocks="%scf maxiter 200 end"

    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1,
                                orcasimpleinput=orcasimpleinput, orcablocks=orcablocks)

    #Geometry optimization
    geomeTRICOptimizer(fragment=H2Ofragment, theory=ORCAcalc, coordsystem='tric')


If you get an error message when launching python3_ash that looks like the following:

.. code-block:: shell

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

