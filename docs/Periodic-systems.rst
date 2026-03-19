Periodic boundary conditions in ASH
======================================

ASH was originally designed for treating molecular systems without any translational symmetry. 
In fact most interfaces to QM codes in ASH do not support periodic boundary conditions (PBC) while the main MM program (OpenMM) does.

There is now, however, increasing support for PBC-based QM codes in ASH and ASH is starting to support PBCs in a more general way.
Work is still in progress in this important area and is subject to interface changes.


**Current status**

- **Fragment objects:** It is important to realize that an ASH Fragment object currently has no knowledge of periodic boundary conditions  (regardless of whether PBC information is present in the file that was read-in).
- **Theory objects:** PBC information generally always needs to be provided to the supported Theory object to get PBCs.
- **Job functions:** Not all job-types in ASH may account for PBCs.

Only a few Theory interfaces in ASH currently support periodic boundary conditions:

    - OpenMMTheory: periodic MM Hamiltonian as well as PBC MD simulations and optimizations
    - CP2KTheory: periodic xTB, QM and QM/MM Hamiltonians
    - xTBTheory: periodic xTB 
    - tbliteTheory: periodic xTB
    - DFTBTheory: periodic DFTB and xTB methods
    - pySCFTheory: not yet supported but soon
    - MACETheory: periodic ML potentials
    - FairchemTheory: periodic UMA ML potentials
    - TorchTheory: periodic AimNet2 ML potentials

Hybrid theories currently do not support PBCs (QMMMTheory, ONIOMTHeory, WrapTheory).


######################################################
Periodic cell optimizations
######################################################

New in ASH is the ability to perform geometry optimizations of periodic systems in a general way (independent of what algorithm is implemented in the QM program) minimizing both atom positions and cell vectors.
This option will only work for Theory classes that support periodic boundary conditions natively (see supported codes at the top).
The theories listed above should all work.

Developer note: Specifically the Theory classes need to support *theory.update_cell* and  *theory.get_cell_gradient* methods and *theory.periodic_cell_vectors* attribues.

.. warning::
  Periodic geometry optimizations in ASH will only work for simple theories, not QM/MM.

-----------------------------------------------------------
*Periodic_optimizer_cart*: Native atoms+cell PBC optimizer
-----------------------------------------------------------

Periodic_optimizer_cart is a native PBC geometry optimizer in Cartesian coordinates.
It minimizes atoms and cell vectors simultaneously by a BFGS algorithm. 
It is also possible to choose a different algorithm for taking the step: 'sd' (steepest descent), 'cg' (conjugate gradient), 'damped-md' (a damped MD algorithm), 'nesterov' (a Nesterov modified damped MD algorithm).

This optimizer is likely to be on par with other similar Cartesian optimization algorithms in e.g. periodic DFT programs, 
but will ultimately suffer from the Cartesian coordinate representation for many systems that will make convergence slow.

Convergence criteria are by default the same as in the geomeTRICOptimizer (convergence_grms':1e-4, 'convergence_gmax':3e-4).
Can be modified by passing a dictionary : conv_criteria = {'convergence_grms':1e-4, 'convergence_gmax':3e-4}

.. code-block:: python

    def Periodic_optimizer_cart(fragment=None, theory=None, rate=2.0, 
                                scaling_rate_cell=1.0, maxiter=50, 
                                step_algo="bfgs",max_step=0.25, momentum=0.5, 
                                printlevel=2, conv_criteria=None):

-----------------------------------------------------------------------------------
*Modified geomeTRIC Optimizer* :  PBC geometry optimization in internal coordinates
-----------------------------------------------------------------------------------

The ASH interface to the geomeTRIC library now also contains a way of coaxing the geomeTRIC library to 
simultaneously minimize atom positions and cell vectors of a periodic system using geomeTRIC's internal coordinate system.

This option has not been rigorously tested but it is likely to work well for periodic systems that feature molecular 
units e.g. molecular crystals, i.e. systems where an internal coordinate representation (only HDLCs for now) will offer advantages.

.. code-block:: python

    from ash import *

    # Defining an ASH fragment by reading an XYZ-file 
    frag = Fragment(xyzfile="ammonia.xyz", charge=0, mult=1)

    # Defining the cell vectors 
    cell_vectors = np.array([[5.01336,0.0,0.0],[0.0,5.01336,0.0],[0.0,0.0,5.01336]])
    #periodic_cell_dimensions=[5.01336, 5.01336, 5.01336, 90.0, 90.0,90.0]

    #CP2K basis set and pseudopotential information
    basis_dict={'C':'DZVP-MOLOPT-SR-GTH','O':'DZVP-MOLOPT-SR-GTH','H':'DZVP-MOLOPT-SR-GTH', 'N':'DZVP-MOLOPT-SR-GTH'}
    potential_dict={'C':'GTH-PBE-q4','O':'GTH-PBE-q6','H':'GTH-PBE-q1', 'N':'GTH-PBE-q5'}

    #Periodic CP2KTheory definition with specified cell dimensions
    theory = CP2KTheory(cp2k_bin_name="cp2k.psmp",basis_dict=basis_dict,potential_dict=potential_dict,
                    basis_method='GPW', functional='PBE', ngrids=4, cutoff=600, numcores=numcores,
                    periodic=True,cell_vectors=cell_vectors, psolver='periodic', stress_tensor=True)

    # Calling the geomeTRIC optimizer. 
    # The optimizer will check for PBC support of the theory object and enable PBC optimization in HDLC (only coordsystem supported)
    Optimizer(theory=theory, fragment=frag, coordsystem="hdlc")


See  :doc:`Geometry-optimization` documentation on how to use the geomeTRICOptimizer in general.


----------------------------------
Alternating atoms+cell optimizer
-----------------------------------
Finally, a highly limited optimizer is also available that will alternate been optimizing the atom coordinates via geomeTRICOptimizer (with frozen cell vectors)
and then take cell-vector steps and go back and forth.
step_algo options: 'sd', 'cg', 'damped-md', 'nesterov'

This optimizer is rarely recommended.

.. code-block:: python

    def periodic_optimizer_alternating(fragment=None, theory=None, rate=0.5, maxiter=50, tol=1e-3, step_algo="sd",
                                    force_orthorhombic=True, max_step=0.25, momentum=0.5,
                                    atoms_tolsetting=None, atom_opt_maxiter=100):


######################################################
Surface scans under PBC
######################################################

Surface scans under PBCs can be performed using the **calc_surface** function.

See  :doc:`surfacescan` for more information. 


######################################################
Nudged elastic band calculations
######################################################

More information to come.

######################################################
Periodic MD simulations
######################################################

More information to come.