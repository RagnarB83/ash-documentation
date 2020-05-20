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

.. code-block:: shell

    f2py -c -m ljlib2 ljlib2.f90 --fcompiler=gfortran
    #f2py -c -m ljlib2 ljlib2.f90 --fcompiler=intel

The f2py command is available if Python3 and numpy has been installed correctly.
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