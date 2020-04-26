==========================
Job Types
==========================

There are a few different job-types in Yggdrasill: single-point , geometry optimization, nudged-elastic band optimization, numerical frequencies and molecular dynamics (not ready).
- Single-point energy/property jobs in Yggdrasill (instead of using the QM code directly) are useful for the purpose of doing electrostatically embedded QM/MM, running multiple energy/property calculations in parallel
- Geometry optimizations can be performed using a simple internal Optimizer or via more flexible external optimizers that can be easily installed.
- Numerical frequencies can be performed for any Hamiltonian (QM, MM or QM/MM). To be finished.
- Nudged elastic band calculations are available via an interface to an external code.
- Molecular dynamics (not ready)

The job-types can be used with any theory object available, e.g. one of the QMTheories in :doc:`QM-interfaces` or using
a QM/MM Theory object from :doc:`QM-MM`

###########################
Single-point calculation
###########################
A single-point calculation is the most basic job to perform.
After creating an Yggdrasill fragment, you create a Theory object, e.g. a QMTheory from: :doc:`QM-interfaces` an
MMTheory (TODO) or a QM/MMTheory (see XX).
The ORCATheory class is recommended for general jobs as this interface is more supported than others.
Below, the ORCASP object is created from the ORCATheory class, passing various ORCA-specific variables to it
(location of ORCA dir and specifying how the inputfile should look). An Yggdrasill fragment object is also passed onto the object (now alwasy required).

For a single-point calculation only then simply runs the Theory object via executing the internal run function of the
object. For the object run command to work, the fragment would have to be associated with the object (as in this example).

TODO: Create SinglePointEnergy function too?

.. code-block:: python

    from yggdrasill import *
    import sys
    settings_yggdrasill.init() #initialize

    HF_frag=Fragment(xyzfiles="hf.xyz")
    #ORCA
    orcadir='/opt/orca_4.2.1'
    orcasimpleinput="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    orcablocks="%scf maxiter 200 end"
    ORCASP = ORCATheory(orcadir=orcadir, fragment=HF_frag, charge=0, mult=1,
                        orcasimpleinput=orcasimpleinput, orcablocks=orcablocks)

    #Simple Energy SP calc
    ORCASP.run()

The flexible input-nature of the ORCA interface here allows one to use any method/basis/property inside ORCA for the
single-point job. Thus one can define any calculation one wants:
DFT job, coupled-cluster, TDDFT, CASSCF, multi-reference configuration interaction, NMR/EPR properties.
Only the total energy of the system would typically be picked up by Yggdrasill though.

It is also possible to request a gradient calculation :

.. code-block:: python

    #An Energy+Gradient calculation
    ORCASP.run(Grad=True)

While Yggdrasill will print out basic information about the run at runtime (e.g. the energy), the energy or gradient
(if requested) is also stored inside the object and can be accessed once the job is completed:

.. code-block:: python

    print(ORCASP.energy)
    print(ORCASP.grad)

By default, the files created by the Theory interface are not cleaned up. To have ORCA (in this example) clean up
temporary files (e.g. so they don't interfere with a future job), one can use the cleanup function.

.. code-block:: python

    #Clean up
    ORCASP.cleanup()



###########################
Geometry optimization
###########################
Geometry optimizations are easily performed in Yggdrasill due to availability of a few different optimization codes.

- An internal optimizer is available (called "Optimizer") that can optimize the system in Cartesian coordinates only using the LBFGS algorithm. While frozen atoms are supported, no other constraints are supported.

- An interface to the PyBerny optimization program (https://github.com/jhrmnn/pyberny) is available that allows efficient optimization in redundant internal coordinates. No frozen atoms or constraints are available currently. PyBerny requires installation via pip.

- The **recommended** optimizer is geomeTRIC (https://github.com/leeping/geomeTRIC) for which there is an Yggdrasill interface. geomeTRIC allows efficient optimization in multiple coordinate systems: TRIC, HDLC, DLC, Cartesian, redundant internals. Supports constraints as well as frozen atoms and the "ActiveRegion" feature inside Yggdrasill allows definition of an active region that aids QM/MM optimization (where most atoms are frozen).

.. code-block:: python

    from yggdrasill import *
    import sys
    settings_yggdrasill.init() #initialize

    HF_frag=Fragment(xyzfile="hf.xyz")
    #ORCA
    orcadir='/opt/orca_4.2.1'
    orcasimpleinput="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    orcablocks="%scf maxiter 200 end"
    ORCAcalc = ORCATheory(orcadir=orcadir, fragment=HF_frag, charge=0, mult=1,
                        orcasimpleinput=orcasimpleinput, orcablocks=orcablocks)

    #Geometry optimization of the ORCA using geomeTRIC optimizer
    #Note: if fragment is passed to optimizer it is not necessary to pass it to the QMtheory (here ORCAcalc) object
    geomeTRICOptimizer(fragment=HF_frag, theory=ORCAcalc, coordsystem='tric')
    #PyBerny example: BernyOpt(ORCAcalc,HF_frag)
    # Internal Cartesian-LBFGS Optimizer:
    #Opt_frag = Optimizer(fragment=HF_frag, theory=ORCAcalc, optimizer='KNARR-LBFGS', frozen_atoms=[])
    #Opt_frag.run()


###########################
Numerical frequencies
###########################


##################################
Nudged Elastic Band Calculations
##################################

Through an interface to an external code, nudged elastic band (NEB) calculations are possible.
Both regular NEB and CI-NEB calculations are possible.

Any QM or QM/MM Hamiltonian can be used.

.. code-block:: python

    from yggdrasill import *
    import sys
    settings_yggdrasill.init() #initialize
    import interface_knarr

    Reactant=Fragment(xyzfile="react.xyz")
    Product=Fragment(xyzfile="prod.xyz")

    #Calculator object without frag
    xtbcalc = xTBTheory(charge=0, mult=1, xtbmethod='GFN2', runmode='library')

    interface_knarr.NEB(reactant=Reactant, product=Product, theory=xtbcalc, images=10, CI=True)



###########################
Molecular Dynamics
###########################

