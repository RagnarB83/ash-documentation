Setup
======================================

ASH is 98% Python with 1 % Julia and 1% Fortran.
A Python3 distribution (version >3.6 or higher) is required and as packages will have to be installed you need to be able to
install Python packages via pip.
We recommend Anaconda (https://www.anaconda.com/distribution/) for a good scientific Python distribution.
Anaconda comes with Python3, Numpy, SciPy, Matplotlib.

Dependencies:

* `Python version 3.6 <https://www.python.org>`_ >=
* `Numpy <https://numpy.org>`_ library.
* `Julia 1.X <https://julialang.org/downloads>`_ installation. PyCall library also required.
* `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_ installation (Python package via pip).
* `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ (Python package via pip).


###############################
Installation and Configuration
###############################
**Step 1.** Clone or download an archive containing Ash and put the directory (named ash) in your home directory.
The ash directory contains the Python source code files, named ash.py etc.

**Step 2.** Check if a suitable Python3 distribution is available. Needs to contain numpy and you will need to be able to install
Python packages using pip. If you don't already have a suitable Python3 distribution one then go to Step 2b.

**Step 2b.** Anaconda Python3 setup (recommended)

Download `Anaconda Python3 package <https://www.anaconda.com/products/individual>`_ and install in e.g. your user directory.
Follow Anaconda installation instructions.

**Step 2c.** Create a new conda Python3.7 virtual environment (here called ashpy37) that will be used for Ash:

.. code-block:: shell

    conda create -n ashpy37 python=3.7 numpy   # Alternatively you can use the default base environment

Select the environment:

.. code-block:: shell

    conda activate ashpy37 # or use base environment if preferred

Make sure this environment is active while you finish the installation process and use this same environment when running Ash.

**Step 3.** To make ASH available to Python, set the environment variables:

.. code-block:: shell

    export PYTHONPATH=/path/to/ash:/path/to/ash/lib:$PYTHONPATH
    export LD_LIBRARY_PATH=/path/to/ash/lib:$LD_LIBRARY_PATH

where */path/to/ash* is the dir where all the ASH sourcefiles are.
Put these environment definitions in your shell environment startup file e.g. .bashrc, .bash_profile or .zshrc.

**Step 4.** Install the recommended Python packages via pip:

.. code-block:: shell

    pip install geometric   (geomeTRIC optimizer)

**Step 5a.** Install Julia from the `Julia official site <https://julialang.org/downloads>`_
Julia is necessary for fast QM/MM functionality inside ASH.

| i. Download appropriate binaries from the official Julia website. Extract archive.
| ii. Add Julia binaries to path: e.g. export PATH=/path/to/julia-1.4.1/bin:$PATH . Put PATH definition to your shell startup file.
| iii. Launch Julia to install PyCall:

.. code-block:: shell

    julia  #This launches the julia interpreter

Inside the Julia interpreter do:

.. code-block:: julia

        using Pkg  # activate Pkg manager
        Pkg.add("PyCall")  #Install PyCall library
        exit()


If there is an error like this: ERROR: SystemError: opening file "/path/to/.julia/registries/General/Registry.toml": No such file or directory
Then execute in shell: rm -rf ~/.julia/registries/General   (assuming Julia is installed in ~).

**Step 5b.** Install `PyJulia <https://pyjulia.readthedocs.io/en/latest/>`_

Install using pip:

.. code-block:: shell

    pip install julia

Activate PyJulia by opening up the python3 interpreter, import julia library and install:

.. code-block:: shell

    python3 #This launches the python3 interpreter

Inside the Python interpreter do:

.. code-block:: python


    import julia
    julia.install()

    #If this is successful then the python-jl binary (installed by PyJulia) should be available.

**Step 6.** Compile Fortran library. When inside ash dir, compile the LJCoulombv1 code using either gfortran or ifort:
The Fortran library is necessary for fast QM/MM functionality inside ASH.

.. code-block:: shell

    f2py -c -m LJCoulombv1 LJCoulombv1.f90 --fcompiler=gfortran
    #f2py -c -m LJCoulombv1 LJCoulombv1.f90 --fcompiler=intel

The f2py command (`Fortran to Python Interface <https://numpy.org/doc/stable/f2py/>`_) is available if Python3 and numpy has been installed correctly.
Rename the compiled library file (something like LJCoulombv1.cpython-36m-x86_64-linux-gnu.so) to LJCoulombv1.so
and move to lib dir: /path/to/ash/lib

**Step 7.** Make sure preferred QM packages are available:

* `ORCA <https://orcaforum.kofo.mpg.de>`_ is a recommended QM code (flexible interface in ASH). See installation instructions on the `ORCA Input Library <https://sites.google.com/site/orcainputlibrary/setting-up-orca>`_.
* The path to ORCA needs to be in PATH and LD_LIBRARY_PATH of your shell and later your jobscript.
* `xTB <https://xtb-docs.readthedocs.io>`_ needs to be in PATH and later your jobscript.


Optional Python packages to install via pip (depends on whether you will use the interfaces to PyBerny, PySCF and PyFrame):

* `PyBerny <https://jan.hermann.name/pyberny/index.html>`_
* `PySCF <http://www.pyscf.org/>`_
* `PyFrame <https://gitlab.com/FraME-projects/PyFraME>`_:


.. code-block:: shell

    pip install pyberny     #pyBerny geometry optimizer
    pip install pyscf       #PySCF QM program
    pip install pyframe     #polarizable embedding helper tool

Optional installation of `Psi4 <http://www.psicode.org/>`_ , best done via Conda:

.. code-block:: shell

    conda install psi4


**Step 8.** Try it out.

* If not doing QM/MM: The regular Python3 executable, *python3*  can be used to run all ASH scripts.

* If doing QM/MM: The Python-Julia executable, *python-jl* should always be used (for fast treatment of large systems via Julia).
The python-jl executable was installe

Example ASH script to try out (geometry optimization of H2O using ORCA):

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


If you get error message when launching python-jl:

.. code-block:: shell

    File "/path/to/envs/ashpy37/bin/python-jl", line 8, in <module>
    sys.exit(main())
    File "/path/to/miniconda3/envs/ashpy37/lib/python3.7/site-packages/julia/python_jl.py", line 114, in main
    execprog([julia, "-e", script_jl, "--"] + unused_args)
    FileNotFoundError: [Errno 2] No such file or directory

This means that the Python-Julia interface is not completely active yet.
Check the following:
1. Is Julia accessible from the shell?, i.e. does typing *julia* launch the Julia interpreter ? If not then the PATH to Julia needs to set.
2. Something went wrong in the installation of Julia or PyJulia in Step 5a or 5b.
3. Make sure you are using the same Python-conda environment you used when you installed things.