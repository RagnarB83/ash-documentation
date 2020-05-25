==========================
Job Types
==========================

The job-types available in ASH:

- Single-point energy/property jobs in Ash (instead of using the QM code directly) are useful for the purpose of doing electrostatically embedded QM/MM, running multiple energy/property calculations in parallel etc.
- Geometry optimizations can be performed using a simple internal Optimizer or via more flexible external optimizers that can be easily installed.
- Numerical frequencies can be performed for any Hamiltonian (QM, MM or QM/MM).
- Nudged elastic band calculations are available via an interface to the Knarr-NEB code.
- Molecular dynamics (not ready).

The job-types can be used with any theory object available, e.g. one of the QMTheories in :doc:`QM-interfaces` or using
a QM/MM Theory object from :doc:`QM-MM`

###########################
Single-point calculation
###########################
A single-point calculation is the most basic job to perform.
After creating an Ash fragment, you create a Theory object, e.g. a QMTheory from: :doc:`QM-interfaces` an
MMTheory (see :doc:`MM-interfaces`) or a QM/MMTheory (see :doc:`QM-MM`).
The ORCATheory class is recommended as a QM code in general as this interface is more supported than others.
ORCA contains a large variety of DFT and WFT methods.
Below, ORCAobject is created from the ORCATheory class, passing various ORCA-specific variables to it
(location of ORCA dir and specifying how the inputfile should look).

For a single-point calculation only then simply passes the Theory object and the Fragment object to the **Singlepoint**
function.

.. code-block:: python

    from ash import *
    import sys
    settings_ash.init() #initialize

    HF_frag=Fragment(xyzfile="hf.xyz")
    #ORCA
    orcadir='/opt/orca_4.2.1'
    orcasimpleinput="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    orcablocks="%scf maxiter 200 end"
    ORCAobject = ORCATheory(orcadir=orcadir, charge=0, mult=1,
                        orcasimpleinput=orcasimpleinput, orcablocks=orcablocks, nprocs=4)

    #Simple Energy SP calc. Energy will be printed to output
    Singlepoint(theory=ORCAobject, fragment=HF_frag)


The **Singlepoint** function will run an ORCA calculation using the ORCAobject and the coordinates from the HF_frag fragment.
The energy will be printed to standard output by default

We can also run the calculation and store the energy as a new variable (to be used for anything):

.. code-block:: python

    #Simple Energy+Gradient SP calc
    # The function will return the energy that can be stored as a variable
    Energy = Singlepoint(theory=ORCAobject, fragment=HF_frag)
    print("Energy is", Energy)

It is also possible to request a gradient calculation :

.. code-block:: python

    #Simple Energy+Gradient SP calc
    Energy, Gradient = Singlepoint(theory=ORCAobject, fragment=HF_frag, Grad=True)
    print("Energy is", Energy)
    print("Gradient is:", Gradient)


By default, the files created by the Theory interface are not cleaned up. To have ORCA (in this example) clean up
temporary files (e.g. so they don't interfere with a future job), one can use the cleanup function:

.. code-block:: python

    #Clean up
    ORCAobject.cleanup()


The energy and gradient from the last Energy/Energy+Gradient run is also stored inside the Theory object and can be accessed:

.. code-block:: python

    print(ORCAobject.energy)
    print(ORCAobject.grad)




###########################
Geometry optimization
###########################
Geometry optimizations are easily performed in Ash due to availability of a few different optimization codes.

- An internal optimizer is available (called "Optimizer") that can optimize the system in Cartesian coordinates only using the LBFGS algorithm. While frozen atoms are supported, no other constraints are supported.
- An interface to the PyBerny optimization program (https://github.com/jhrmnn/pyberny) is available that allows efficient optimizations in redundant internal coordinates. No frozen atoms or constraints are available currently. PyBerny requires installation via pip.
- The **recommended** optimizer is geomeTRIC (https://github.com/leeping/geomeTRIC) for which there is full-featured Ash interface. geomeTRIC allows efficient optimization in multiple coordinate systems: TRIC, HDLC, DLC, Cartesian, redundant internals. Supports constraints as well as frozen atoms natively. Furthermore, the "ActiveRegion" feature inside Ash allows definition of an active region that allows efficient QM/MM optimizations of large systems (where most atoms are frozen). Only the active region coordinates are passed to geomeTRIC.


.. code-block:: python

    from ash import *
    import sys
    settings_ash.init() #initialize

    HF_frag=Fragment(xyzfile="hf.xyz")
    #ORCA
    orcadir='/opt/orca_4.2.1'
    orcasimpleinput="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    orcablocks="%scf maxiter 200 end"
    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1,
                        orcasimpleinput=orcasimpleinput, orcablocks=orcablocks)
    #Note: if fragment is passed to optimizer it is not necessary to pass it to the QMtheory (here ORCAcalc) object

    #Geometry optimization of the ORCA using geomeTRIC optimizer
    geomeTRICOptimizer(fragment=HF_frag, theory=ORCAcalc, coordsystem='tric')

    #PyBerny example:
    BernyOpt(ORCAcalc,HF_frag)

    # Internal Cartesian-LBFGS Optimizer:
    Optimizer(fragment=HF_frag, theory=ORCAcalc, optimizer='KNARR-LBFGS', frozen_atoms=[])


################################
Numerical frequencies (Hessian)
################################
Numerical frequencies can be performed with Ash using any QM, MM or QM/MM theory object.
Any method for which there is an analytical gradient (forces) available can be used (numerical 2nd derivative on top of numerical 1st derivative is not recommended).

Use the **NumFreq** function to request a numerical frequency job. The function requires a fragment object and a theory level at minimum.
The fragment object should typically contain a fragment with optimized coordinates at same level of theory (i.e. an already optimized minimum or saddlepoint).

*Type of Hessian*
Additionally you can select to do a 1-point Hessian or a 2-point Hessian by the *npoint* keyword (value of 1 or 2).
A 1-point Hessian makes a single displacement (+ direction) for each atom and each x,y and z-coordinate from the input geometry. This option is reasonably accurate and is the default.
A more accurate 2-point Hessian makes displacement in both + and - directions (for each x-, y- and z-coordinate of each atom), is twice as expensive (double the displacements)
but is more accurate.
The displacement step can be chosen if wanted. The default setting is: 0.0005 Å.

*Serial or parallel*
Two runmodes are available: 'serial' and 'parallel'. The 'serial' mode will run each displacement sequentially.
The Energy+Gradient step can still be run in parallel if e.g. the QM or QM/MM object has this information;
e.g. if an ORCA object has been defined with nprocs=8 then ORCA will run each Energy+Gradient evaluation with 8 cores using the OpenMPI parallelization of ORCA.
For numerical frequencies, it is usually much more efficient, however, to run the displacement jobs simutaneously in parallel fashion.
This is accomplished using runmode='parallel' and the parallelization will be linear scaling (almost always recommended).
As there are almost always many more displacements available than CPUs, the parallelization of the QM or QM/MM object is turned off and instead as many displacements
are run simultaneously as there are number of cores. For example, for a 30-atom system, there are 90 XYZ coordinates. For a 2-point Hessian, this means
that 180 displacements to be calculated. If 20 cores are available, then 20 displacements can be run simultaneously, fully utilizing all 20 cores.
This will require 9 runs in total (20*9=180).

*Full or partial Hessian*

A partial Hessian (NEEDS TO BE TESTED) can be easily performed instead of the full Hessian. This is an excellent approximation for vibrational modes with rather local character
and the quality of the approximation can be controlled. For a QM/MM model of a protein active site with an active region of a 1000 atoms, the full Hessian
of all 1000 atoms would typically not be doable; instead a partial Hessian job of the important atoms (e.g. the QM region) makes more sense.
A partial Hessian job is performed if a list of Hessian atoms (e.g. hessatoms=[0,1,2] ) is passed to the NumFreq function. In this case, the displacements
will only be calculated for the list of "hessatoms" and the result is a partial Hessian for the system.

*Final output*
Once the displacements are complete, the gradients for all displacements are combined to give the full (or partial) Hessian.
The Hessian is then mass-weighted and diagonalized. (Limitation: translational and rotational modes are currently not projected out).
This gives the frequencies as eigenvalues and the normal mode eigenvectors.
A normal mode composition factor analysis is automatically performed (NOT READY) as well as zero-point energy thermochemistry.


Example script below demonstrates a combined geometry optimization (using geomeTRIC).
The QM code used here is ORCA but any QM, MM or QM/MM object can be used.

.. code-block:: python

    from ash import *
    import sys
    settings_ash.init() #initialize

    #the total number of CPU cores available to Ash (should match the job-script)
    Ashnumcores=8

    orcadir='/opt/orca_4.2.1'
    orcasimpleinput="! HF-3c "
    orcablocks="%scf maxiter 200 end"

    reactstring="""
       C  -2.66064921   -0.44148342    0.02830018
       H  -2.26377685   -1.23173358    0.68710920
       H  -2.29485851   -0.62084858   -0.99570465
       H  -2.27350346    0.53131334    0.37379014
       F  -4.03235214   -0.44462811    0.05296388
    """
    Reactant=Fragment(coordsstring=reactstring)

    #Calculator object without frag. nprocs=8 is used here for parallelizing ORCA during optimization.
    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput=orcasimpleinput, orcablocks=orcablocks, nprocs=Ashnumcores)

    #Geometry optimization of Reactant object and ORCAcalc theory object.
    #Each Energy+Grad step is parallelized by ORCA.
    geomeTRICOptimizer(theory=ORCAcalc,fragment=Reactant)


    #Numfreq job. A 2-point Hessian is requested in runmode parallel (recommended).
    #Ash will use the number of cores given to run same number of displacments simultaneouslyu.
    #ORCA parallelization is turned off automatically.
    NumFreq(Reactant, ORCAcalc, npoint=2, runmode='parallel', numcores=Ashnumcores)



##################################
Nudged Elastic Band Calculations
##################################

Through an interface to an external code, nudged elastic band (NEB) calculations are possible.
Both regular NEB and CI-NEB calculations are possible.

Any QM or QM/MM Hamiltonian can be used.

.. code-block:: python

    from ash import *
    import sys
    settings_ash.init() #initialize
    import interface_knarr

    Reactant=Fragment(xyzfile="react.xyz")
    Product=Fragment(xyzfile="prod.xyz")

    #Calculator object without frag
    xtbcalc = xTBTheory(charge=0, mult=1, xtbmethod='GFN2', runmode='library')

    interface_knarr.NEB(reactant=Reactant, product=Product, theory=xtbcalc, images=10, CI=True)


###########################
Saddle-point optimization
###########################


###########################
Surface scans
###########################

**Unrelaxed scan**
TODO

**Relaxed scan**
TODO

###########################
Molecular Dynamics
###########################

Not yet ready