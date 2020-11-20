.. raw:: html

     <style> .red {color:#aa0060; font-weight:bold; font-size:16px} </style>

.. role:: red

Setup
======================================
Contact Ragnar if you want access to the code.

ASH is 99% Python with 1 % Julia.
A Python3 distribution (version >3.6 or higher) is required and as packages will have to be installed you need to be able to
install Python packages via pip.
We recommend Anaconda or miniconda (https://www.anaconda.com/distribution/) for a good scientific Python distribution.
Anaconda is inconvenient for installing Numpy, SciPy, Matplotlib.

Treatment of large systems via QM/MM require a Julia installation (as the Python routines will be too slow).
Future versions will make Julia a requirement.

Strict dependencies:

* `Python version 3.6 <https://www.python.org>`_ >=
* `Numpy <https://numpy.org>`_ library.


Strongly recommended:

* `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ (Python package via pip).
* `Julia 1.5.2 <https://julialang.org/downloads>`_ installation for fast routines for large system treatment.
* `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_ installation (Python package via pip).

Useful:

* `Scipy <https://www.scipy.org>`_ library. Used for interpolation routines when plotting surfaces.


###############################
Installation and Configuration
###############################
**Step 1.** Clone or download an archive containing ASH and put the directory (named ash) in your home directory.
The ash directory contains the Python source code files, named ash.py etc.

**Step 2.** Check if a suitable Python3 installation is available. It needs to contain numpy and you will need to be able to install
Python packages to it using pip. If you don't already have a suitable Python3 distribution, go to Step 2b.

**Step 2b.** Anaconda Python3 setup (recommended)

Download `Anaconda Python3 package <https://www.anaconda.com/products/individual>`_ and install in e.g. your user directory.
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ also works.
Follow Anaconda installation instructions.

**Step 2c.** Create a new conda Python3.7 virtual environment (here called ashpy37) that will be used for ASH:

.. code-block:: shell

    conda create -n ashpy37 python=3.7 numpy   # Alternatively you can use the default base environment

Select the environment:

.. code-block:: shell

    conda activate ashpy37 # or use base environment if preferred

Make sure this environment is active while you finish the installation process and use this same environment when running ASH.

**Step 3.** To make ASH available to Python, set the environment variables:

.. code-block:: shell

    export PYTHONPATH=/path/to/ash:/path/to/ash/lib:$PYTHONPATH
    export LD_LIBRARY_PATH=/path/to/ash/lib:$LD_LIBRARY_PATH

where */path/to/ash* is the dir where all the ASH sourcefiles are.
Put these environment definitions in your shell environment startup file e.g. .bashrc, .bash_profile or .zshrc.

**Step 4.** Install the recommended Python packages via pip:

.. code-block:: shell

    pip install geometric   (geomeTRIC optimizer)

**Step 5a.** Install Julia from the `Julia official site <https://julialang.org/downloads>`_.
Julia is necessary for fast QM/MM functionality inside ASH.

 i) Download appropriate binaries from the official Julia website. Version 1.4 or higher. Extract archive.
 ii) Add Julia binaries to path: e.g. export PATH=/path/to/julia-1.4.1/bin:$PATH . Put PATH definition to your shell startup file.
 iii) Run Julia using the ASH sourcefile julia-packages-setup.jl (inside ASH source directory) as input to download and install the  required Julia packages. Currently: PyCall, Hungarian, Distances

.. code-block:: shell

    julia julia-packages-setup.jl  #This launches the julia interpreter and requests installation of required Julia packages for ASH.

This will download and install required Julia packages.

If there is an error like this: ERROR: SystemError: opening file "/path/to/.julia/registries/General/Registry.toml": No such file or directory
Then execute in shell: rm -rf ~/.julia/registries/General

**Step 5b.** Install `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_


:red:`Important:` Make sure the correct Python environment (e.g. your ashpy37 conda environment) is active before proceeding.
Then install using pip:

.. code-block:: shell

    pip install julia

Activate PyJulia by opening up the python3 interpreter, import julia library and install:

.. code-block:: shell

    python3 #This launches the python3 interpreter

Inside the Python interpreter do:

.. code-block:: python


    import julia
    julia.install()

    #If this is successful then the python-jl binary (installed by PyJulia) should become available.


**Step 7a.** Install desired QM program(s):

* `ORCA <https://orcaforum.kofo.mpg.de>`_ is a recommended QM code (flexible interface in ASH). See installation instructions on the `ORCA Input Library <https://sites.google.com/site/orcainputlibrary/setting-up-orca>`_. The path to ORCA needs to be in PATH and LD_LIBRARY_PATH of your shell and later your jobscript.
* `xTB <https://xtb-docs.readthedocs.io>`_ needs to be in PATH and later your jobscript.


Optional Python packages to install via pip (depends on whether you will use the interfaces to PyBerny, PySCF and PyFrame):

* `PyBerny <https://jan.hermann.name/pyberny/index.html>`_
* `PySCF <http://www.pyscf.org/>`_
* `PyFrame <https://gitlab.com/FraME-projects/PyFraME>`_:


.. code-block:: shell

    pip install pyberny     #pyBerny geometry optimizer
    pip install pyscf       #PySCF QM program
    pip install pyframe     #polarizable embedding helper tool

Optional installation of the `Psi4 <http://www.psicode.org/>`_ QM code (if you intend to use it), best done via Conda:

.. code-block:: shell

    conda install psi4 psi4-rt -c psi4


**Step 7b.** Optional: Install OpenMM (if needed)

Note: Not yet documented...

For protein and explict solvation QM/MM in ASH, then the `OpenMM program <http://openmm.org>`_ is used as MM code.
It can be installed using conda.

.. code-block:: shell

    conda install -c omnia openmm


**Step 8.** Try it out.

* Make sure the correct Python3 environment is active (e.g. switch to the conda environment you created in Step2c).

* If not doing QM/MM: The regular Python3 executable, *python3*  can be used to run all ASH scripts.

* If doing QM/MM: The Python-Julia executable, *python-jl* should always be used (for fast treatment of large systems via Julia) to run scripts. The python-jl executable was installed in the same dir as the python3 executable (e.g. in the conda environment). python-jl can always be used.

Example ASH script to try out (geometry optimization of H2O using ORCA):

.. code-block:: shell

    python3 first-ash-job.py

or:

.. code-block:: shell

    python-jl first-ash-job.py


first-ash-job.py:

.. code-block:: python

    from ash import *
    settings_ash.init()

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

    #Basic Cartesian optimization with KNARR-LBFGS
    geomeTRICOptimizer(fragment=H2Ofragment, theory=ORCAcalc, coordsystem='tric')


If you get error message when launching python-jl or something similar:

.. code-block:: shell

    File "/path/to/envs/ashpy37/bin/python-jl", line 8, in <module>
    sys.exit(main())
    File "/path/to/miniconda3/envs/ashpy37/lib/python3.7/site-packages/julia/python_jl.py", line 114, in main
    execprog([julia, "-e", script_jl, "--"] + unused_args)
    FileNotFoundError: [Errno 2] No such file or directory

This means that the Python-Julia interface is not completely active yet.
Check the following:

1. Is Julia accessible from the shell?, i.e. does typing *julia* in the shell, launch the Julia interpreter ? If not then the PATH to Julia bin dir needs to set.
2. Something went wrong in the installation of Julia or PyJulia in Step 5a or 5b.
3. Make sure you are using the same Python-conda environment you used when you installed things.
4. Setup PyCall for each Julia user environment (updates ~/.julia)

