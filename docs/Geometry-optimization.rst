Geometry optimization
======================================

Geometry optimizations in ASH can be performed in a few different ways.
The oldest interface is to the `geomeTRIC optimizer library <https://github.com/leeping/geomeTRIC>`_ that allows several type of internal coordinates.
There is also an interface to the powerful `DL-FIND library <https://www.itheoc.uni-stuttgart.de/research/kaestner/research/dlfind/>`_ , 
which has very fast execution speed (being written in Fortran).


In addition there is a native Cartesian-coordinate optimizer, *Cart_optimizer*, 
that performs geometry optimizations exclusively using Cartesian coordinates that also has fast execution speed.

The geomeTRIC and DL-FIND optimizers support saddlepoint optimizations via PRFO algorithms (with various Hessian input options)
and additionally ASH supports the powerful Sella saddlepoint optimization program.

More recently, periodic cell optimizations have become possible for all 3 optimizers,
allowing simultaneous atom+cell optimizations for periodic systems.

######################################################
geomeTRIC library
######################################################

See :doc:`geometric_interface`

######################################################
DL-FIND library
######################################################

See :doc:`dlfind-interface`


######################################################
Sella
######################################################

An interface to the Sella saddlepoint optimization algorithm is also available in ASH.
The Sella algorithm is only intended for saddlepoint optimizations and is based on an iterative Hessian diagonalization approach,
that avoids computation of the exact Hessian.
See :doc:`Sella_interface`

######################################################
Cart_optimizer
######################################################

**Cart_optimizer** is a native Cartesian-coordinate based optimizer in ASH that is capable of treating both molecular systems and systems with periodic boundary conditions.
It minimizes atoms and cell vectors simultaneously by a BFGS algorithm by default. 
In addition to BFGS one can choose the following step algorithms: 'sd' (steepest descent), 'cg' (conjugate gradient), 'damped-md' (a damped MD algorithm), 'nesterov' (a Nesterov modified damped MD algorithm).

This optimizer is likely to be on par with other similar Cartesian optimization algorithms in e.g. periodic DFT programs.
It does suffer from the Cartesian coordinate representation that for many molecular systems that will make convergence slow
for molecular systems in terms of number of steps required. Because it avoids the overhead of internal coordinate transformations,
and does not have any blow-up problems associated with internal coordinates, it can offer a certain robustness and execution speed
for Theory levels where the energy+gradient evaluation is fast.

Convergence criteria are by default the same as in the geomeTRICOptimizer (convergence_grms':1e-4, 'convergence_gmax':3e-4).
Can be modified by passing a dictionary : conv_criteria = {'convergence_grms':1e-4, 'convergence_gmax':3e-4}

Soft constraints (i.e. restraints) have been implemented for bonds, angles, dihedrals.
Frozen Cartesian positions are also available; atoms can be frozen, 
either by specifying them via *frozen_atoms* keyword  or by defining them as XYZ constraints in the constraints dictionary.
Partial Cartesian constraints (e.g. 'X', 'Y', 'Z', 'XY', 'XZ' and 'YZ') are also possible.
The soft bond/angle/dihedral constraints are implemented as harmonic restraints and the force constant 
can be tuned by the user via the kf_bonds, kf_angles and kf_dihedrals keywords (units of Eh/Bohr^2 and Eh/rad^2). 
.. code-block:: python

  def Cart_optimizer(fragment=None, theory=None, rate=2.0, 
                                  scaling_rate_cell=1.0, maxiter=50, 
                                  step_algo="bfgs",
                                  max_step=0.25, momentum=0.5, constrain_method='soft',
                                  printlevel=2, conv_criteria=None, PBC_format_option="CIF",
                                  constraints=None, frozen_atoms=None, result_write_to_disk=True,
                                  kf_bonds=10.0, kf_angles=10.0, kf_dihedrals=10.0):

For more flexibility it is also possible to define a Cart_optimizer_class object instead.

.. code-block:: python

  class Cart_optimizer_class:

    def __init__(self,fragment=None, theory=None, rate=2.0, scaling_rate_cell=1.0, maxiter=50, step_algo="bfgs",
                                max_step=0.25, momentum=0.5, printlevel=2, conv_criteria=None, print_atoms_list=None,
                                PBC_format_option="CIF", constraints=None, constrain_method='soft',
                                frozen_atoms=None, result_write_to_disk=True,
                                kf_bonds=10.0, kf_angles=10.0, kf_dihedrals=10.0):


The object would be defined like this:

.. code-block:: python

  cartopt=Cart_optimizer_class(step_algo="bfgs", printlevel=1, frozen_atoms=[0,1], 
      constraints={'bond':[[0,1],[3,4]], 'angle':[[98,99,100]]}, conv_criteria={'convergence_grms':1e-4, 'convergence_gmax':3e-4})

The object could then e.g. be passed to *calc_surface* to perform a surface scan using the Cart_optimizer with the defined options.
Or it could be run directly by calling the *run* method of the Cart_optimizer_class object
and passing a Fragment and Theory object to it.

.. code-block:: python

  cartopt.run(fragment=frag, theory=theory)


######################################################
Constraints
######################################################

GeomeTRIC, DL-FIND and Cart_optimizer all support constraints and the syntax is generally the same:

.. code-block:: python

    constraints={'bond':[[0,1]]} #This defines a bond/distance constraint between atoms 0 and 1
    constraints={'bond':[[0,1],[3,4]]} #This defines multiple bond constraints: between atoms 0 and 1 AND also between atoms 3 and 4
    constraints={'angle':[[98,99,100]]} #This defines a angle constraint between atoms 98,99 and 100
    constraints={'dihedral':[[98,99,100,101]]} #This defines a dihedral constraint between atoms 98,99,100 and 101.
    constraints={'bond':[[0,1],[3,4]], 'angle':[[98,99,100]]} #This defines 2 bond constraints and 1 angle constraint.
    constraints={'xyz':[5,6]} #This defines XYZ constraints for the indicated atoms (i.e. freeze atoms). Alternative to frozenatoms option.
    constraints={'x':[5,6]} #This freezes the X-coordinate for the indicated atoms.
    constraints={'xy':[7,8]} #This freezes the X and Y coordinates for the indicated atoms.

Be aware that the specific algoriths used to enforce constraints differ. DL-FIND will employ hard in-place constraints, 
while geomeTRIC will use soft constraints (restraints) during the optimization and will only enforce exact hard constraints towards the end of the optimization.
Cart_optimizer currently features only soft constraints, that are more approximate but tunable by modifying the force constant of the restraints.

Additionally, the geomeTRIC interface allows to specify the target value of 
the constraint while DL-FIND can only freeze the geometry in-place. See individual documentation.

######################################################
Restraints
######################################################

The general use oh harmonic restraints can be useful alternatives to constraints in geometry optimizations, being less strict and are less likely to blow
up optimizations, especially if the force constant is tuned. They can be used to drive a molecule gently towards a specific structure.

ASH features a special theory, RestraintTheory, that can be used to apply restraints to a system in a general way, irrespective of the optimizer.

.. code-block:: python

  class RestraintTheory:
      def __init__(self, fragment=None, printlevel=None, numcores=1, label=None,
                  restraints=None, force_constant=10000.0):

The restraints should be provided as a list of dictionaries.
Units of restraint values: Å (bonds) or ° (angles and dihedrals)
Units of force constants: Eh/Bohr^2 (bonds) or Eh/rad^2 (angles and dihedrals).

.. code-block:: python

  # Here defining a single restraint as a list with a single dictionary.
  # A dihedral restraint between atoms 6,5,11,29 (zero-based indexing), centered at 0° with a force-constant of 5 Eh/rad^2
  restraints=[{'type':'dihedral', 'indices': [6,5,11,29], 'target': 0, 'force_constant':5.0}]

The RestraintTheory object can used on it's own as a Theory object (where the energy and gradient contributions come solely from defined resetraints),
or can be combined with other objects in a hybrid WrapTheory object (see :doc:`module_Hybrid_Theory`).
It should work with most regular jobs in ASH.
It could in principle be used to add restraints to a molecular dynamics simulation, though in practice, this can be accomplished faster via OpenMM (that drives MD in ASH).

It may be particularly useful for special geometry optimizations as shown below.

**Example: Driving a molecule to a particular coordinate value by optimizing the RestraintTheory**

.. code-block:: python

  # Define Restrainttheory that adds harmonic restraints to enforce a particular coordinates
  # Note: here setting an individual forceconstant for each of the 3 restraints.
  rest_theory = RestraintTheory(fragment=frag, restraints=[
  {'type':'dihedral', 'indices': [6,5,11,29], 'target': 0, 'force_constant':5.0},
  {'type':'angle', 'indices': [7,3,10], 'target': 120, 'force_constant':4.0},
  {'type':'bond', 'indices': [8,1], 'target': 1.5, 'force_constant':50.0}])

  #Optimizing using the RestraintTheory. Any Optimizer can be used, here geomeTRIC is used.
  Optimizer(theory=rest_theory, fragment=frag)


**Example: Driving a molecule to a particular coordinate value by optimizing a Hybrid theory**

While the above example will often work, the optimization will likely not succeed if the initial geometry is quite
far from the target geometry (where the restraints are satisfied) or if the coordinate changes required
are coupled to movements of other atoms (will cause atoms to be in strained positions that will blow up the Optimizer or Theory run.).

In this case, a better alternative is to build a hybrid theory where the restraints are simply added on top
of a regular Theory object. 
Below we build a WrapTheory object that combines an xTBTheory object with the RestraintTheory object.
The final energy expression (and gradient) that is minimized by the Optimizer is then the sum of the xTB-energy and the restraint-energy.

.. code-block:: python

  #1. Define Restrainttheory that adds harmonic restraints to enforce a particular coordinates
  # Note here using a global force_constant for the 3 dihedral restraints.
  rest_theory = RestraintTheory(fragment=frag, restraints=[
  {'type':'dihedral', 'indices': [6,5,11,29], 'target': 0},
  {'type':'dihedral', 'indices': [7,3,10,28], 'target': 120},
  {'type':'dihedral', 'indices': [8,1,9,27], 'target': 120}],  force_constant=5.0)

  # Creating a WrapTheory object that combines restraints with xTB
  xtb_theory = xTBTheory(xtbmethod="GFN2")
  wrap_theory = WrapTheory(theories=[xtb_theory,rest_theory])

  #Any Optimizer can be used, here DL-FIND is used.
  DLFIND_optimizer(theory=wrap_theory, fragment=frag, jobtype="opt")


######################################################
Deprecated: SimpleOpt
######################################################

SimpleOpt is an internal alternative to the geomeTRIC and DL-FIND external optimizers. 
It is rarely recommended except for very small systems where it can find some use.
It can only performs geometry optimizations using Cartesian coordinates,  via the following algorithms:

- steepest descent (optimizer="SD")
- LBFGS (via Knarr library, optimizer="KNARR-LBFGS")
- FIRE (via Knarr library, optimizer="KNARR-FIRE")

.. code-block:: python

  def SimpleOpt(fragment=None, theory=None, charge=None, mult=None, optimizer='KNARR-LBFGS', maxiter=50, 
                frozen_atoms=None, actatoms=None, RMSGtolerance=0.0001, MaxGtolerance=0.0003, FIRE_timestep=0.00009):


*Example: Geometry optimization of an H2O fragment at the xTB level, with default options.*

.. code-block:: python

    from ash import *

    frag=Fragment(databasefile="h2o.xyz",charge=0, mult=1)
    xtbcalc=xTBTheory(xtbmethod='GFN1')

    Optimizer(theory=xtbcalc, fragment=frag)


######################################################
Numerical Gradient Optimizations
######################################################

For some complicated quantum chemical methods as well as some hybrid methods, no analytical gradient might be available, 
typically preventing convenient geometry optimizations. However, a numerical gradient can in these cases be defined.
While a numerical gradient can sometimes be requested in programs like ORCA, ASH also allows calculation of the numerical gradient
for any level of theory.

In ASH a numerical-gradient can be requested by wrapping the Theory object by the *NumGradclass*.

.. code-block:: python
    
  # Numerical gradient class
  class NumGradclass:
      def __init__(self, theory, npoint=2, displacement=0.00264589,  runmode="serial", numcores=1, printlevel=2):

Once a NumGradclass object has been defined, one can use it to e.g. perform geometry optimizations or any other job-type where gradients are needed (MD, surface scan etc.).
Be aware of course that numerical gradients are by definition more noisy than analytical ones and require considerable effort. 
MD using numerical gradients is unlikely to work well due to noise.

The example below requests a numerical-gradient geometry optimization using a ORCA-DFT level (just as an example, the analytical gradient is of course preferable here):

.. code-block:: python

  from ash import *

  #Fragment
  frag = Fragment(databasefile="h2o.xyz")
  #Theory
  theory = ORCATheory(orcasimpleinput="! B3LYP def2-SVP tightscf")
  #Numgrad wrapper
  numgrd_theory = NumGradclass(theory=theory)
  #Optimization via geomeTRIC
  Optimizer(theory=numgrd_theory, fragment=frag)

A simpler alternative for numerical gradient optimization makes use of the NumGrad keyword in the *Optimizer* :

.. code-block:: python

  from ash import *

  #Fragment
  frag = Fragment(databasefile="h2o.xyz")
  #Theory
  theory = ORCATheory(orcasimpleinput="! B3LYP def2-SVP tightscf")
  #Optimization
  Optimizer(theory=theory, fragment=frag, NumGrad=True)

A good use-case for numerical-gradient optimizations would e.g. involve geometry optimization of a small molecule using a correlated wavefunction
where no analytic gradient is available.


######################################################
Minimum Energy Crossing Point Optimizations
######################################################

Minimum energy crossing point (MECP) optimizations are intented for finding the point where 2 potential energy surfaces cross each other.
The 2 energy surfaces might differ e.g. by spin multiplicity or an alternative SCF solution of the same multiplicity.
The gradient is defined following Harvey et al: Harvey, J. N.; Aschi, M.; Schwarz, H.; Koch, W. Theor. Chem. Acc., 1998, 99, 95.

ASH has a way of conveniently defining the MECP gradient by the *MECPGradclass*. 
One simply couples together 2 different theory objects and also specifies the charge and multiplicity of both energy surfaces.
In principle the 2 theory objects could even be interfaces to 2 different QM programs.

.. code-block:: python

  # MEPC-gradient class
  class MECPGradclass:
      def __init__(self, theory_1=None,theory_2=None, charge_1=None, charge_2=None, 
                  mult_1=None, mult_2=None, runmode="serial", numcores=1, printlevel=2):

Once an *MECPGradclass* object is defined, one can run a geometry optimization by using the *MECPGradclass* object as a Theory object.
In the limited testing done so far, the MECP-optimizations have been found to be more efficient using the Cartesian-based *SimpleOpt* optimizer.



*Example: MECP optimization of the quartet/sextet crossing of FeO+*

We define the system as a fragment, then define 2 identical theory levels and then specify the different spin multiplicity for each.

.. code-block:: python

  from ash import *

  #Define system
  frag = Fragment(diatomic="FeO", bondlength=1.67, charge=1, mult=6)

  # Define theory levels for both electronic states
  theory_1 = ORCATheory(orcasimpleinput="! B3LYP tzvp tightscf")
  theory_2 = ORCATheory(orcasimpleinput="! B3LYP tzvp tightscf")

  # Wrap the 2 theory levels into a MECPGradclass object
  mecpgrad = MECPGradclass(theory_1=theory_1, theory_2=theory_2, charge_1=1, charge_2=1, mult_1=6, mult_2=4)

  # Run using the basic optimizer
  SimpleOpt(fragment=frag, theory=mecpgrad, optimizer='KNARR-LBFGS', maxiter=50,
    RMSGtolerance=0.00001, MaxGtolerance=0.00003)

Starting from an FeO+ distance of 1.67 Angstrom (close to the sextet minimum) the MECP optimization converges to a distance of 1.99 Angstrom
which is where the sextet and quartet surfaces cross.

*Example: MECP optimization involving a non-ground-state SCF solution*

ASH allows some additional freedom in MECP optimization as the 2 electronic states are 
controlled by the 2 theory objects as well as specifying the multiplicity of them.
For example for ORCATheory one could control the specific SCF-state solution by the deltaSCF feature as ASH can turn that off and on for each theory object.

.. code-block:: python

  from ash import *

  frag = Fragment(diatomic="FeO", bondlength=1.67, charge=1, mult=6)

  theory_1 = ORCATheory(orcasimpleinput="! UKS B3LYP tzvp tightscf")
  theory_2 = ORCATheory(orcasimpleinput="! UKS B3LYP tzvp tightscf", deltaSCF=True,
            deltaSCF_PMOM=False, deltaSCF_confline="betaconf 0,1", deltaSCF_turn_off_automatically=True)

  mecpgrad = MECPGradclass(theory_1=theory_1, theory_2=theory_2, charge_1=1, charge_2=1, mult_1=6, mult_2=4)

  SimpleOpt(fragment=frag, theory=mecpgrad, optimizer='KNARR-LBFGS', maxiter=50,
    RMSGtolerance=0.00001, MaxGtolerance=0.00003)

As it not necessarily straightforward though to stay on the correct SCF solution throughout so that requires some experimentation.

######################################################
Optimizers for periodic boundary conditions
######################################################

The interfaces to geomeTRIC and DL-FIND as well as the native **Cart_optimizer** are all capable of optimizing systems with periodic boundary conditions with both atoms and cell vectors minimized.

See :doc:`Periodic-systems` for more information.

