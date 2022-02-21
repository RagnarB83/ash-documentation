==========================
ASH Program Philosophy
==========================

ASH is a program in its infancy and is constantly under revisement as the programming skills and ideas of the main developer evolve.
It is inspired by multiple software projects. It is written in Python but is not necessarily very Pythonic.

While object-oriented (OO) programming is a highly useful programming paradigm and is used partially in ASH, the OO-use is actually fairly limited.
OO is used within the program to define simple classes such as : Fragment, Theory classes etc.
that define the molecular-system or theory level) and to organize the code in a useful way.
Meanwhile regular Python functions serve the purpose of jobs: i.e. we run jobs via ASH functions that take ASHTheory objects and ASH Fragment objects as input and return.

We believe that the minimal exposure of object orientation to the user, makes ASH a simpler and easier-to-use program, especially to those not used to much programming or even Python.


###########################
Some basic principles
###########################

- *Object-orientation should be mostly hidden from how a regular user interacts with ASH.* 
  
While the user learns to create molecular fragment objects e.g. :

.. code-block:: python

  H2O = Fragment(xyzfile="h2o.xyz")
  orcacalc = ORCATheory(orcasimpleinput="! PBE def2-SVP")

the user should then primarily use ASH via running simple Python functions that take those fragments and theory objects as input and return some simple variable: e.g.:

.. code-block:: python

  energy = Singlepoint(fragment=H2O, theory=orcacalc)
 
Methods inside classes should typically **never** be called by the user and instead there should be a function available for carrying out such functionality. 
This creates a simpler scripting environment for users who are not used to object orientation. If a class is needed for improved flexibility for a new type of job, a wrapper function should be provided for the user.

- *The list of dependencies should be kept to a minimum.*

ASH works via interfaces to various high-quality QM, MM programs, optimizer program etc. that have functionaly that can not be easily replicated. 
However, otherwise the philosophy is to keep dependencies to a minimum to avoid future problems and make the ASH installation process as simple as possible.
Unnecessary Python libraries outside the standard library are not used unless unavoidable. 

- *ASH should be easy to install*

This is work in progress. Currently the conda-forge route (for almost everything) seems to be the easiest for maintaining the Python, OpenMM and Julia dependencies.

- *Ideally nothing should require compilation*

By using Julia we can do just-in-time-compilation for routines that need to be written in a compiled language.
We should avoid Fortran/C code that requires compilation.