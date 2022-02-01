==========================
Job Types
==========================

The primary job-types available in ASH:

- Single-point energy/property jobs. Useful for doing electrostatically embedded QM/MM, running multiple energy/property calculations in parallel, creating advanced workflows etc.
- Geometry optimizations jobs. Typically performed using an interface to the powerful geomeTRIC library.
- Numerical frequencies. Can be performed for any Hamiltonian (QM, MM or QM/MM). Analytical frequencies available for some theories.
- Nudged elastic band calculations are available via an interface to the Knarr-NEB code.
- Molecular dynamics.
- Surface scans.


Additionally ASH features various workflows that perform some combination of multiple jobtypes.  
The job-types can be used with any theory object available (with some exceptions), e.g. one of the QMTheories in :doc:`QM-interfaces` or using
a QM/MM Theory object from :doc:`module_QM-MM`


###########################
Single-point calculation
###########################

See :doc:`singlepoint`

###########################
Geometry optimization
###########################

Geometry optimizations are easily performed in ASH due to availability of the flexible optimizer: geomeTRIC (https://github.com/leeping/geomeTRIC): :doc:`geomeTRIC-interface`
geomeTRIC allows efficient optimization in multiple coordinate systems: TRIC, HDLC, DLC, Cartesian, redundant internals. Supports constraints as well as frozen atoms natively. Furthermore, the "ActiveRegion" feature inside Ash allows definition of an active region that allows efficient QM/MM optimizations of large systems (where most atoms are frozen). Only the active region coordinates are passed to geomeTRIC.
ASH features a full-featured interface to geomeTRIC that allows flexible constraint input, QM/MM optimizations, relaxed and unrelaxed 1D/2D surface scans and more.

Example:

.. code-block:: python

    from ash import *

    HF_frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1)
    #ORCA
    orcasimpleinput="! BP86 def2-SVP  tightscf"
    ORCAcalc = ORCATheory(orcasimpleinput=orcasimpleinput)

    #Geometry optimization of the ORCA using geomeTRIC optimizer
    geomeTRICOptimizer(fragment=HF_frag, theory=ORCAcalc, coordsystem='tric')


See :doc:`geomeTRIC-interface` for all features.

Other optimizers:

- An internal optimizer is available (called "SimpleOpt") that can optimize the system in Cartesian coordinates only using the LBFGS algorithm. While frozen atoms are supported, no other constraints are supported.


################################
Analytical frequencies (Hessian)
################################
Analytical frequencies can be requested in some cases if supported by the theory-level interface as well as the Hamiltonian inside program.
See :doc:`module_freq`

################################
Numerical frequencies (Hessian)
################################

Numerical frequencies can be performed with Ash using any QM, MM or QM/MM theory object.
See :doc:`module_freq`


##################################
Nudged Elastic Band Calculations
##################################

Through an interface to an external code, nudged elastic band (NEB) calculations are possible.
Both regular NEB and CI-NEB calculations are possible.
See :doc:`knarr-interface` for documentation

Any QM or QM/MM Hamiltonian can be used.

.. code-block:: python

    from ash import *

    Reactant=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    Product=Fragment(xyzfile="prod.xyz", charge=0, mult=1)

    #Calculator object without frag
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library')

    NEB(reactant=Reactant, product=Product, theory=xtbcalc, images=10, CI=True)


###########################
Surface scans
###########################
Potential Energy Surfaces can be conveniently scanned in ASH using the **calc_surface function** that uses the **geometric** optimization library.
Both unrelaxed and relaxed scans be calculated, using either 1 and 2 reaction coordinates.

See :doc:`surfacescan`



###########################
Saddle-point optimization
###########################

Currently, saddle-points searches can only be performed in ASH using the NEB method.


###########################
Molecular Dynamics
###########################

It is possible to perform molecular dynamics in ASH in 2 ways.

1. Classical molecular dynamics using a forcefield via the OpenMM.
This is only available for a system that has been set up using OpenMMTheory and utilizes the OpenMM library for energy, forces and dynamics. See :doc:`OpenMM-interface`

2. General dynamics via an interface to ASE. 
   See :doc:`module_dynamics`



