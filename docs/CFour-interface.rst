CFour interface
======================================

CFour is a powerful wavefunction theory program, in particular known for its availability of first and second order
derivatives of various coupled cluster methods including CCSD(T). 
This is one of very few programs that allow for geometry optimizations and other properties with high level WF methods.
Higher order coupled cluster methods are also available such as CCSDT, CCSDT(Q).

https://cfour.uni-mainz.de/cfour/

The CFour capabilities for energy calculations can be found here:
https://cfour.uni-mainz.de/cfour/index.php?n=Main.Single-pointEnergyCalculations
Methods for which analytic gradients are available (allowing geometry optimizations and first-order properties) are shown here:
https://cfour.uni-mainz.de/cfour/index.php?n=Main.GeometryOptimizations
Methods for which analytic second derivates are available:
https://cfour.uni-mainz.de/cfour/index.php?n=Main.HarmonicVibrationalFrequencies


ASH offers a convenient interface to CFour. CFour keywords should be provided as a Python dictionary
which allows for easy manipulation of the CFour inputfile, without having to write one manually.
The ASH interface also allows one to conveniently use CFour to perform geometry optimizations without having to write a Z-matrix (as required by CFour normally).
ASH can furthermore call on CFour to calculate analytic second derivates (i.e. the Hessian) using the **AnFreq** function.

The ASH-CFour interface even allows for QM/MM pointcharge embedding meaning that CFourTheory can be used to create a QMMMTheory object (see :doc:`module_QM-MM`).
QM/MM energies and gradients are available, which even allows for CCSD(T)/MM geometry optimizations. 

**CFourTheory class:**

.. code-block:: python

  class CFourTheory:
      def __init__(self, cfourdir=None, printlevel=2, cfouroptions=None, numcores=1,
                  filename='cfourjob', specialbasis=None, ash_basisfile=None,
                  parallelization='MKL'):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``cfourdir``
     - string
     - None
     - Path to CFour directory.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``cfouroptions``
     - dict
     - None
     - CFour keywords as a Python dictionary 
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores CFour will use
   * - ``parallelization``
     - string
     - 'MKL'
     - Type of parallelization used by CFour. Options: 'MKL', 'MPI'
   * - ``filename``
     - string
     - 'cfourjob'
     - Name of CFour inputfile
   * - ``specialbasis``
     - dictionary
     - None
     - Optional specialbasis option.
   * - ``ash_basisfile``
     - string
     - 'def2-SVP'
     - ASH-internal basis set file for CFour. Options: 'def2-SVP'

######################################################
Installation
######################################################

CFour requires a license, see: https://cfour.uni-mainz.de/cfour/index.php?n=Main.Download
Once downloaded it needs to be compiled by a suitable C and Fortran compiler.
See this page for installation instructions: https://cfour.uni-mainz.de/cfour/index.php?n=Main.Examples


*Example serial installation using GCC and Intel MKL BLAS:*

.. code-block:: text

  #GCC compiler
  export PATH=/path/to/gcc/bin:$PATH
  export LD_LIBRARY_PATH=/path/to/gcc/lib:$LD_LIBRARY_PATH
  #Intel MKL
  source /home/prog/intel/oneapi-2022.1.0/mkl/latest/env/vars.sh
  #Installdir
  export installdir=/home/ragnar/PROGRAMS/CFOUR/cvfour-v21dev-install-mkl-serial
  #Configure and make
  ./configure CC=gcc F77=gfortran --prefix=$installdir
  nohup CXXFLAGS='-std=c++11' make >>& make.log
  make install

*Example parallel installation using GCC and Intel MKL BLAS and OpenMPI:*

TODO

################################
Finding the CFour program
################################

ASH can find the CFour program in a few different ways.

- ASH will first check if the cfourdir argument has been set which should be a string that points to the directory
- If the cfourdir argument has not been provided ASH will next see if cfourdir has been provided in the ASH settings (~/ash_user_settings.ini file): See :doc:`basics`
- Otherise, ASH will next search the operating systems's PATH environment variable for an executable "xcfour" and if found, will set the cfourdir accordingly and use that CFour version.  This can be a convenient option if you make sure to define your shell environments carefully in your jobscript or shell-startup file. Be careful, however, if you have multiple versions of the program available.

################################
CC modules
################################

CFOUR can use different coupled cluster modules which have different strengths.
The CC module is specified by the CC_PROG keyword in the cfouroptions dictionary.
Consult the CFour manual about which CC program to use.
Generally:
CC_PROG='VCC' : Default in CFour,except for CCSDT(Q) and CCSDTQ.
CC_PROG='ECC' : Generally better performance than VCC for CCSD(T) and closed-shell CCSDT. Not recommended for ROHF gradients.
CC_PROG='NCC' : New program. Recommended for CCSDT(Q) and CCSDTQ. Only closed-shell.


######################################################
Parallelization
######################################################

CFour calculations can be parallelized using either MKL or MPI.
OpenMPI parallelization can only be used if CFour has been compiled for that purpose (see above).
If CFour has been compiled without MPI but using the MKL BLAS library then the only option is to use MKL parallelization.

The parallelization strategy is controlled by the parallelization keyword in the CFourTheory class.
It is by default set to 'MKL'. If CFour has been compiled with MPI then it can be set to 'MPI' to use MPI parallelization.
Both parallelization modes will use the number of cores specified by the numcores keyword in the CFourTheory class.

.. code-block:: python

  cfourcalc_mkl = CFourTheory(cfouroptions=cfouroptions, parallelization='MKL', numcores=4)
  cfourcalc_mpi = CFourTheory(cfouroptions=cfouroptions, parallelization='MPI', numcores=4)

######################################################
Examples
######################################################

**Single-point CCSD(T) calculation:**

.. code-block:: python

    from ash import *

    #Define fragment
    frag=Fragment(databasefile="hf.xyz", charge=0, mult=1)

    cfouroptions = {
    'CALC':'CCSD(T)',
    'BASIS':'PVTZ',
    'REF':'RHF',
    'FROZEN_CORE':'ON',
    'MEM_UNIT':'MB',
    'MEMORY':3100,
    'PROP':'FIRST_ORDER',
    'CC_PROG':'ECC',
    'SCF_CONV':10,
    'LINEQ_CONV':10,
    'CC_MAXCYC':300,
    'SYMMETRY':'OFF',
    'HFSTABILITY':'OFF'
    }

    cfourcalc = CFourTheory(cfouroptions=cfouroptions)

    #Simple Energy SP calc
    result = Singlepoint(theory=cfourcalc, fragment=frag)


**Geometry optimization at CCSD(T) and CCSDT levels of theory:**

CCSD(T)/cc-pVTZ:

.. code-block:: python

    from ash import *

    #Define fragment
    frag=Fragment(databasefile="hf.xyz", charge=0, mult=1)

    cfouroptions = {
    'CALC':'CCSD(T)',
    'BASIS':'PVTZ',
    'REF':'RHF',
    'FROZEN_CORE':'ON',
    'MEM_UNIT':'MB',
    'MEMORY':3100,
    'CC_PROG':'VCC',
    'SCF_CONV':10,
    'LINEQ_CONV':10,
    'CC_MAXCYC':300,
    'SYMMETRY':'OFF',
    'HFSTABILITY':'OFF'
    }
    cfourcalc = CFourTheory(cfouroptions=cfouroptions)

    #Geometry optimization
    result = Optimizer(theory=cfourcalc, fragment=frag)


CCSDT/cc-PVTZ:

.. code-block:: python

    from ash import *

    #Define fragment
    frag=Fragment(databasefile="hf.xyz", charge=0, mult=1)

    cfouroptions = {
    'CALC':'CCSDT',
    'BASIS':'PVTZ',
    'REF':'RHF',
    'FROZEN_CORE':'ON',
    'MEM_UNIT':'MB',
    'MEMORY':3100,
    'CC_PROG':'VCC',
    'SCF_CONV':10,
    'LINEQ_CONV':10,
    'CC_MAXCYC':300,
    'SYMMETRY':'OFF',
    'HFSTABILITY':'OFF'
    }
    cfourcalc = CFourTheory(cfouroptions=cfouroptions)

    #Geometry optimization
    result = Optimizer(theory=cfourcalc, fragment=frag)


**Harmonic vibrational frequencies at the CCSD(T) level of theory:**

CCSD(T)/cc-pVTZ:

.. code-block:: python

    from ash import *

    #Define fragment
    frag=Fragment(databasefile="hf.xyz", charge=0, mult=1)

    cfouroptions = {
    'CALC':'CCSD(T)',
    'BASIS':'PVTZ',
    'REF':'RHF',
    'FROZEN_CORE':'ON',
    'MEM_UNIT':'MB',
    'MEMORY':3100,
    'CC_PROG':'VCC',
    'SCF_CONV':10,
    'LINEQ_CONV':10,
    'CC_MAXCYC':300,
    'SYMMETRY':'OFF',
    'HFSTABILITY':'OFF'
    }
    cfourcalc = CFourTheory(cfouroptions=cfouroptions)

    #Geometry optimization
    result = Optimizer(theory=cfourcalc, fragment=frag)
    #Analytical Hessian calculation
    result = AnFreq(theory=cfourcalc, fragment=frag)


**CFour CCSD(T) density calculation and visualization:**

As the CFour program can calculate densities at all levels of theory for which analytic gradients are available
one can calculate and visualize densities associated with the CCSD, CCSD(T), CCSDT wavefunctions.
If one includes 'PROP':'FIRST_ORDER' in the cfouroptions dictionary input to CFourTheory, 
then the density will be calculated at the requested level of theory. 
This density can be used to define various electric properties at the CC level of theory (dipole, EFG etc.), population analysis
but the density can also be useful on its own.
Here we utilize the MOLDEN_NAT file that CFour creates, which contains the natural orbitals of the CC wavefunction
that define the correlated WF density.

You can then use the **multiwfn_run function** (See :doc:`Multiwfn_interface` for details) that creates the density in
realspace using Multiwfn. 
The function will create a Cube-file that can be visualized in VMD, Chemcraft or other programs.

.. code-block:: python

  from ash import *

  numcores=8

  #Define fragment
  frag=Fragment(databasefile="hf.xyz", charge=0, mult=1)

  #Define CFour options
  cfouroptions = {
  'CALC':'CCSD',
  'BASIS':'PVDZ',
  'REF':'RHF',
  'FROZEN_CORE':'ON',
  'MEM_UNIT':'MB',
  'MEMORY':3100,
  'PROP':'FIRST_ORDER',
  'CC_PROG':'ECC',
  'SCF_CONV':10,
  'LINEQ_CONV':10,
  'CC_MAXCYC':300,
  'SYMMETRY':'OFF',
  'HFSTABILITY':'OFF'
  }
  #Define CFourTheory object
  cfourcalc = CFourTheory(cfouroptions=cfouroptions,numcores=numcores)

  #Run CFour calculation
  result=Singlepoint(theory=cfourcalc,fragment=frag)

  #Files produced by CFOUR: MOLDEN (SCF WF) and MOLDEN_NAT (Natural Orbitals of the correlated WF)
  multiwfn_run("MOLDEN_NAT", option='density', grid=3, numcores=numcores)

  #The Cube-file created, MOLDEN_NAT_mwfn.cube, can next be visualized in e.g. VMD or Chemcraft

