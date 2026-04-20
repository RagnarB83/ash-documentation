Sella Optimizer
======================================

ASH features an interface to the Sella saddlepoint optimization program,
which features a novel approach for optimizing saddlepoints.
This is accomplished via an in iterative partial Hessian diagonalization algorithm.
While this requires a good saddlepoint guess (e.g. from NEB, scan etc.) no Hessian input is required,
since an adaptive partial Hessian is determined during the job.


The Sella program can be be found on the Github repository: https://github.com/zadorlab/sella

The algorithm used by Sella is described in:
https://pubs.acs.org/doi/10.1021/acs.jctc.9b00869
https://pubs.acs.org/doi/10.1021/acs.jctc.2c00395

################################
Installation
################################

The Sella library is easily installed via pip:

.. code-block:: bash

    pip install sella

See the Sella documentation for more details: https://github.com/zadorlab/sella


############################################
SellaOptimizer function
############################################

The *SellaOptimizer* function is the ASH wrapper around the Sella library.
To use, Sella must first have been installed.

.. code-block:: python

  def SellaOptimizer(theory=None, fragment=None, charge=None, mult=None, printlevel=2, NumGrad=False,
                    convergence_gmax=1e-4, maxiter=150, result_write_to_disk=False,
                    constraints=None, actatoms=None, frozenatoms=None,
                    gamma=0.03, eta=1e-4):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``theory``
     - ASH THeory
     - None
     - An ASH Theory. This will be used to supply energy+gradient information during the optimization.
   * - ``fragment``
     - ASH Fragment
     - None
     - An ASH fragment.
   * - ``NumGrad``
     - Boolean
     - False
     - Whether to perform a numerical gradient optimization (uses NumGradClass)
   * - ``convergence_gmax``
     - float
     - 1e-4
     - Convergence criterion for the maximum gradient component in units of Hartree/Bohr.
   * - ``maxiter``
     - integer
     - 150
     - Maximum number of iterations.
   * - ``result_write_to_disk``
     - Boolean
     - False
     - Whether to write optimization results to disk.
   * - ``gamma``
     - float
     - 0.03
     - Gamma parameter determines the convergence criterion for the iterative diagonalization routine of Sella.
   * - ``eta``
     - float
     - 1e-4
     - Eta parameter is the size of the finite difference step in the iterative diagonalization routine of Sella.
   * - ``constraints``
     - list of ASH Constraints
     - None
     - A list of ASH Constraints to be applied during the optimization.
   * - ``actatoms``
     - list of atom indices
     - None
     - A list of atom indices that will be active during the optimization. Note: Sella will only see these atoms.
   * - ``frozenatoms``
     - list of atom indices
     - None
     - A list of atom indices that will be frozen during the optimization.


For more information on the algorithmic details of the Sella parameters, see the `Sella wiki <https://github.com/zadorlab/sella/wiki/Hyperparameters>`_ 

###################################################
Constraints and frozen atoms and active region
###################################################

The Sella interface supports bond, angle, dihedral constraints as well as frozen atoms and partial XYZ constraints.

Similar to other optimizers in ASH , a dictionary defining constraints (*constraints* keyword) should be provided.
Syntax to use for the constraints dictionary:

.. code-block:: python

    constraints={'bond':[[0,1]]} #This defines a bond/distance constraint between atoms 0 and 1
    constraints={'bond':[[0,1],[3,4]]} #This defines multiple bond constraints: between atoms 0 and 1 AND also between atoms 3 and 4
    constraints={'angle':[[98,99,100]]} #This defines a angle constraint between atoms 98,99 and 100
    constraints={'dihedral':[[98,99,100,101]]} #This defines a dihedral constraint between atoms 98,99,100 and 101.
    constraints={'bond':[[0,1],[3,4]], 'angle':[[98,99,100]]} #This defines 2 bond constraints and 1 angle constraint.
    constraints={'xyz':[5,6]} #This defines XYZ constraints for the indicated atoms (i.e. freeze atoms). Alternative to frozenatoms option.
    constraints={'x':[5,6]} #This defines a partial X Cartesian constraint for the indicated atoms.
    constraints={'xy':[5,6]} #This defines a partial XY Cartesian constraint for the indicated atoms.

To freeze atoms in space one can also specify a list of frozen atoms by the *frozenatoms* keyword.
This is equivalent to defining XYZ constraints.

If a large number of atoms needs to be frozen (e.g. most of a protein in a QM/MM calculation) 
one should instead define an Active Region via the *actatoms* keyword.
If this option is enabled, only the active atoms will be passed onto Sella but optimization trajectories and coordinate
files will still be created for the whole system as well as the active-region.


############################################
Examples
############################################

The use of the **SellaOptimizer** function is straightforward.
Typically one just need to provide an ASH Theory and Fragment object.

The main options to consider exploring are:

- the convergence criterion (*convergence_gmax*). Converging the saddlepoint more tightly than in a standard optimization can be a good idea. 
- the *gamma* parameter (default 0.03 in ASH) is the convergence criterion for the iterative eigensolver. A smaller value will lead to a more accurate Hessian at the cost of more iterative diagonalization steps, which may, however, lead to a more robust and quicker optimization.
- the *eta* parameter (default 1e-4 Angstrom) controls the step size for the iterative diagonalization.

During testing we found that converging the saddlepoint a bit more tightly (*convergence_gmax*=1e-4 compared to *convergence_gmax*=3e-4 like in geomeTRICOptimizer) is a bit more robust.
Additionally the default value for *gamma* in ASH is smaller (0.03) than the default values in Sella (0.4). 
Overall we observe robust optimizations and quicker convergence, but will obviously be system-dependent and the user should feel free to explore alternative values.

For information on the algorithmic details, see the `Sella wiki <https://github.com/zadorlab/sella/wiki/Hyperparameters>`_ 

.. code-block:: python

    from ash import *

    frag=Fragment(xyzfile="tsguess.xyz", charge=0, mult=1) #Fragment object creation
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf") #ORCATheory object creation

    SellaOptimizer(theory=ORCAcalc, fragment=frag, convergence_gmax=1e-4, gamma=0.03)




############################################
SellaOptimizer_class 
############################################

The *SellaOptimizer* function is a wrapper around the SellaOptimizer_class.
For more flexibility one can also use the SellaOptimizer_class.

.. code-block:: python

    class SellaoptimizerClass:
            def __init__(self,charge=None, mult=None, printlevel=2, constraints=None,
                        convergence_gmax=3e-4, maxiter=150, result_write_to_disk=False,
                        gamma=0.4):

Example:

.. code-block:: python

    # Create Sellaoptimizer object
    sellaobj = SellaOptimizer_class(convergence_gmax=3e-4, maxiter=150, gamma=0.4)

    # Run optimization with an input theory and fragment
    sellaobj.run(theory=ORCAcalc, fragment=frag)



























.. code-block:: python

class SellaoptimizerClass:
        def __init__(self,theory=None, charge=None, mult=None, printlevel=2, constraints=None,
                     convergence_gmax=3e-4, maxiter=150, result_write_to_disk=False,
                     gamma=0.4):



################################
Example
################################
