Turbomole interface
======================================

`Turbomole <https://www.turbomole.org>`_  is a popular quantum chemistry code, particularly known for its fast implementation of many algorithms
and its speed of execution for both DFT and WFT implementations.

ASH features a simple interface to Turbomole that allows the easy use of Turbomole for basic DFT and MP2 calculation.
Not all Turbomole features are currently supported.

Electrostatic embedding QM/MM is supported. Analytical Hessian is also supported.

Periodic boundary conditions via the Turbomole riper module are now supported.
This allows periodic HF and DFT calculations with Turbomole to be carried out.
Periodic geometry optimizations can be carried out using the ASH optimizer options, see :doc:`Periodic-systems`.




**TurbomoleTheory class:**

.. code-block:: python
    
  class TurbomoleTheory:
      def __init__(self, TURBODIR=None, turbomoledir=None, filename='XXX', printlevel=2, label="Turbomole", uff=False,
                  numcores=1, parallelization='SMP', functional=None, dispersion=None, gridsize="m4", scfconv=7, symmetry="c1", rij=True,
                  basis=None, jbasis=None, scfiterlimit=50, maxcor=500, ricore=500, controlfile=None,skip_control_gen=False,
                  mp2=False, pointcharge_type=None, pc_gaussians=None,
                  periodic=False, periodic_cell_vectors=None, PBC_dimension=3,
                  periodic_cell_dimensions=None, kpoint_values=[1,1,1]):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``TURBODIR``
     - string
     - None
     - Path to main Turbomole directory.
   * - ``turbomoledir``
     - string
     - None
     - Alternative ways of specifying path to Turbomole, should point to the bin directory of Turbomole
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``numcores``
     - integer
     - 1
     - Number of cores to use for either SMP or MPI parallelization.
   * - ``functional``
     - string
     - None
     - Name of DFT-functional keyword that will be specified in control-file.
   * - ``dispersion``
     - string
     - None
     - Name of DFT-D dispersion correction to use. Options: 'D2', 'D3BJ', 'D3ZERO', 'D4'.
   * - ``gridsize``
     - string
     - "m4"
     - Grid-setting. See Turbomole manual.
   * - ``scfconv``
     - integer
     - 7
     - SCF-convergence setting. See Turbomole manual.
   * - ``symmetry``
     - string
     - 'c1'
     - Point-group symmetry label. Default C1.
   * - ``rij``
     - Boolean
     - True
     - Whether to use RI approximation for Coulomb-J integrals. Default True and recommended. Controls the Turbomole executable.
   * - ``basis``
     - string
     - None
     - Name of what internal basis set to use.
   * - ``jbasis``
     - string
     - None
     - Name of what internal auxiliary basis set to use for RI-J calculations.
   * - ``scfiterlimit``
     - integer
     - 50
     - Number of max SCF iterations.
   * - ``maxcor``
     - integer
     - 500
     - Core-memory in MB.
   * - ``ricore``
     - integer
     - 500
     - Additional memory setting for RI.
   * - ``controlfile``
     - string
     - None
     - Option to bypass control-file creation and point to a user-created controlfile.
   * - ``skip_control_gen``
     - Boolean
     - False
     - Another option to bypass control-file creation and use files available in current dir.
   * - ``mp2``
     - Boolean
     - False
     - Whethe to run MP2 after SCF or not.
   * - ``pointcharge_type``
     - string
     - None
     - For pointcharge-embedding (i.e. QM/MM), whether to use an alternative to a regular pointcharge. Option: 'gaussians', 'mxrank=Z' (where Z is multipole rank), 'pe' (polarizable embedding)
   * - ``pc_gaussians``
     - list
     - None
     - If pointcharge_type="gaussian", provide here a list of alphas for each MM-pointcharge.
   * - ``periodic``
     - Boolean
     - False
     - Enable periodic boundary conditions.
   * - ``periodic_cell_vectors``
     - 2d array
     - None
     - Cell vectors in Å as a Numpy 2D array.
   * - ``periodic_cell_dimensions``
     - list
     - None
     - Cell dimensions as a list of [a,b,c,alpha,beta,gamma] parameters in Å and °.
   * - ``kpoint_values``
     - list
     - [1,1,1]
     - Specify k-point grid by a list. [1,1,1] indicates gamma point in all three X,Y,Z directions.


################################################################################
Turbomole installation and using with ASH
################################################################################
Turbomole needs to be installed separately.
Once installed, and confirmed to be working on it's own it can be used with ASH in the following ways.

It is easiest to have ASH find Turbomole on it's own. 
This simply requires the directory containing the Turbomole executables (containing e.g. the ridft binary)
to be in the shell path. ASH will search for the ridft program and then set up the necessary internal paths.
For example:

.. code-block:: shell

    export TURBODIR=/Applications/TmoleX2024/TURBOMOLE
    export PATH=$TURBODIR/bin/i686-apple-darwin:$PATH

Alternatively you can provide the path of the main Turbomole directory to TURBODIR keyword argument.

.. code-block:: python

    theory = TurbomoleTheory(TURBODIR="/my_disk/my_name/TURBOMOLE")

or the path of the bin directory to turbomoledir keyword argument:

.. code-block:: python

    theory = TurbomoleTheory(turbomoledir="/my_disk/my_name/TURBOMOLE/bin")



################################################################################
Parallelization
################################################################################

The interface supports both the SMP and MPI parallelization of Turbomole.
This is controlled by the parallelization="SMP" or parallelization="MPI" options
and by providing the desired number of CPU-cores to use via the *numcores* keyword.
Note that for MPI-parallelization, Turbomole will use it's own MPI.

################################################################################
Turbomole interface examples
################################################################################

Due to the nature of Turbomole, and because the interface is relatively new, the use of all Turbomole features
can not be supported.
Because traditionally Turbomole relies on the creation of a control-file and other files via the define helper-program,
the ASH interface either bypasses this or allows the user to provide their own control-file.
There are hence a few different ways of using the interface.

**Simple basic DFT example via keywords**

The simplest way of using the interface for DFT-calculations is simply to specify the functional-name, basis and jbasis keywords.
The functional-name follows Turbomole-naming conventions (see Turbomol manual).
Hence for a B3LYP calculation use functional="b3-lyp" etc.
The *rij* Boolean keyword controls whether the **ridft** or **dscf** executable is run. *rij* is by default True and almost always recommended. 

.. code-block:: python

    from ash import *
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
    theory = TurbomoleTheory(functional="b3-lyp", basis="def2-SVP", jbasis="def2-SVP", rij=True, numcores=2)

    Singlepoint(theory=theory, fragment=frag)

ASH in this case automatically creates the control-file. The options in the control-file created by ASH can be modified somewhat by
keywords: *gridsize*, *scfconv*, *symmetry*, *rij*, *mp2*, *scfiterlimit*, *maxcor* and *ricore*.

**Simple basic dispersion-corrected DFT example **

Grimme-type dispersion corrections can be employed via the *dispersion* keyword.
Choose:

- dispersion="D2" to get the DFT-D2 correction
- dispersion="D3ZERO" to get the original DFT-D3(0) correction
- dispersion="D3BJ" to get the recommended DFT-D3 correction that uses Becke-Johnson damping.
- dispersion="D4" to get the modern DFT-D4 dispersion correction.

.. code-block:: python

    from ash import *
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
    theory = TurbomoleTheory(functional="b3-lyp", basis="def2-SVP", dispersion='D4', jbasis="def2-SVP", rij=True, numcores=2)

    Singlepoint(theory=theory, fragment=frag)


**DFT example via user-provided control-file**

To have more flexibility in setting up a Turbomole calculation, such as providing an alternative basis-set, ECPs, relativisic approximation,
setup a different MO-guess, modify other options etc. ASH also allows providing the path to an already created control-file by the user.
In this case, one simply provides the path to the control-file, neglecting almost all other keywords.
The *rij* Boolean keyword is still relevant since it controls whether the **ridft** or **dscf** executable is run.

.. code-block:: python

    from ash import *
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
    theory = TurbomoleTheory(controlfile="/path/to/controlfile", rij=True)

    Singlepoint(theory=theory, fragment=frag)

This option will only work if the controlfile has been correctly set up, either manually or by using the **define** program.
The control-file needs to also contain correct paths to basis-sets etc.

An alternative similar option, *skip_control_gen* is also available.
This option simply skips the creation of the control-file and other setup and immediately attempts to run Turbomole.
It assumes that the control-file and other necessary files are already present in the directory.
This might be the case due to a previous run or due to the user having already run the **define** program etc.

.. code-block:: python

    from ash import *
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
    theory = TurbomoleTheory(skip_control_gen=True)

    Singlepoint(theory=theory, fragment=frag)



**Electrostatic embedding QM/MM**

A TurbomoleTheory object can be combined into a QMMMTheory object using both mechanical and electrostatic embedding.
For regular pointcharge-embedding, things should work automatically.

Since Turbomole supports alternatives besides simple pointcharges such as Gaussians or multipoles (see *pointcharge_type* keyword) ,
it is in principle possible to go beyond regular pointcharges. However, ASH does not support these options as part of QMMMTheory for the moment.


################################################################################
UFF via Turbomole
################################################################################

Turbomole has an UFF (Universal Force Field) implementation that the ASH interface supports.
Simply set the *uff=True* keyword without specifying anything else and it can be used for some basic UFF jobs.

.. code-block:: python

    from ash import *
    frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
    theory = TurbomoleTheory(uff=True)

    Optimizer(theory=theory, fragment=frag)


################################################################################
PBC DFT calculations via Turbomole
################################################################################

Turbomole supports periodic boundary conditions for HF and DFT via the riper module.
To do PBC DFT calculations via the ASH interface see example below:

.. code-block:: python

  from ash import *

  numcores=8
  frag = Fragment(xyzfile="ammonia.xyz", charge=0, mult=1)

  cell_vectors = np.array([[5.01336,0.0,0.0],[0.0,5.01336,0.0],[0.0,0.0,5.01336]])

  theory = TurbomoleTheory(basis="def2-SVP", functional="PBE", dispersion="D3BJ", jbasis="def2-SVP", 
              periodic=True, periodic_cell_vectors=cell_vectors)

  Optimizer(theory=theory, fragment=frag, coordsystem="hdlc")