==========================
Job Types
==========================

The job-types available in ASH:

- Single-point energy/property jobs in Ash (instead of using the QM code directly) are useful for the purpose of doing electrostatically embedded QM/MM, running multiple energy/property calculations in parallel, creating advanced workflows etc.
- Geometry optimizations can be performed using a simple internal Optimizer or via more flexible external optimizers that can be easily installed.
- Numerical frequencies can be performed for any Hamiltonian (QM, MM or QM/MM). Analytical frequencies available for some theories.
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
Geometry optimizations are easily performed in Ash due to availability of the flexible optimizer: geomeTRIC (https://github.com/leeping/geomeTRIC)
geomeTRIC allows efficient optimization in multiple coordinate systems: TRIC, HDLC, DLC, Cartesian, redundant internals. Supports constraints as well as frozen atoms natively. Furthermore, the "ActiveRegion" feature inside Ash allows definition of an active region that allows efficient QM/MM optimizations of large systems (where most atoms are frozen). Only the active region coordinates are passed to geomeTRIC.
ASH features a full-featured interface to geomeTRIC that allows flexible constraint input, QM/MM optimizations, relaxed and unrelaxed 1D/2D surface scans and more.

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

Options to **geometricOptimizer**:

- coordsystem (Default: 'tric', other options: 'hdlc', 'dlc', 'cartesian', 'prim', 'tric-p')
- frozenatoms (default: None, provide list of atoms to be frozen in space)
- constraintsinputfile (default: None, provide name of constraints-inputfile according to geomeTRIC syntax.
- constraints (default: None, provide dictionary of constraint definitions, with or without the value of the constraint. Example: constraints = { 'bond' : [[0,1]]}
- constrainvalue (default: False, Boolean, whether constrain-value is provided in constraints dictionary or not)
- maxiter (default: 50, maximum number of iterations)
- ActiveRegion (default:False, whether to use an active region or not. Requires accompanying actatoms list.
- actatoms (default: None, list of atoms that are active during optimization, all others are frozen)
- convergence_setting (default: 'ORCA'. What type of convergence criteria to use. Valid options are: 'ORCA', 'Chemshell', 'ORCA_TIGHT', 'GAU', 'GAU_TIGHT', 'GAU_VERYTIGHT', 'SuperLoose'.

    - ORCA:    conv_criteria = {'convergence_energy' : 5e-6, 'convergence_grms' : 1e-4, 'convergence_gmax' : 3.0e-4, 'convergence_drms' : 2.0e-3, 'convergence_dmax' : 4.0e-3 }
    - Chemshell:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 3e-4, 'convergence_gmax' : 4.5e-4, 'convergence_drms' : 1.2e-3, 'convergence_dmax' : 1.8e-3 }
    - ORCA_TIGHT:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 3e-5, 'convergence_gmax' : 1.0e-4, 'convergence_drms' : 6.0e-4, 'convergence_dmax' : 1.0e-3 }
    - GAU:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 3e-4, 'convergence_gmax' : 4.5e-4, 'convergence_drms' : 1.2e-3, 'convergence_dmax' : 1.8e-3 }
    - GAU_TIGHT:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 1e-5, 'convergence_gmax' : 1.5e-5, 'convergence_drms' : 4.0e-5, 'convergence_dmax' : 6e-5 }
    - GAU_VERYTIGHT:    conv_criteria = {'convergence_energy' : 1e-6, 'convergence_grms' : 1e-6, 'convergence_gmax' : 2e-6, 'convergence_drms' : 4.0e-6, 'convergence_dmax' : 6e-6 }
    - SuperLoose:            conv_criteria = { 'convergence_energy' : 1e-1, 'convergence_grms' : 1e-1, 'convergence_gmax' : 1e-1, 'convergence_drms' : 1e-1, 'convergence_dmax' : 1e-1 }

- conv_criteria (Default: see ORCA setting above. Optionally provide your own dictionary of settings using syntax above.



Other optimizers:

- An internal optimizer is available (called "SimpleOpt") that can optimize the system in Cartesian coordinates only using the LBFGS algorithm. While frozen atoms are supported, no other constraints are supported.
- An interface to the PyBerny optimization program (https://github.com/jhrmnn/pyberny) is available that allows efficient optimizations in redundant internal coordinates. No frozen atoms or constraints are available currently. PyBerny requires installation via pip.

.. code-block:: python

    #PyBerny example:
    BernyOpt(ORCAcalc,HF_frag)

    # Internal Cartesian-LBFGS Optimizer:
    SimpleOpt(fragment=HF_frag, theory=ORCAcalc, optimizer='KNARR-LBFGS', frozen_atoms=[])


################################
Analytical frequencies (Hessian)
################################
Analytical frequencies can be requested in some cases if supported by the theory-level interface as well as the Hamiltonian inside program.
Currently analytical frequencies are supported in: ORCATheory



.. code-block:: python

    def AnFreq(fragment=None, theory=None, numcores=1, temp=298.15, pressure=1.0)


Example:

.. code-block:: python

    HF_frag=Fragment(xyzfile="hf.xyz")
    ORCAcalc = ORCATheory(orcadir='/opt/orca_4.2.1', charge=0, mult=1,
                        orcasimpleinput='BP def2-SVP def2/J', orcablocks=", nprocs=1)
    thermochem_dict = AnFreq(theory=ORCAcalc, fragment=HF_frag)

    print("Thermochem properties dict:", thermochem_dict)
    print("Vibrational frequencies (cm**-1) : ", thermochem_dict['freqs'])
    print("ZPVE (Eh) : ", thermochem_dict['ZPVE'])
    print("Gibbs energy corrections (Eh) : ", thermochem_dict['Gcorr'])

A dictionary containing various properties is returned (dictionary keys) from an AnFreq job:
(freqs, ZPVE, E_trans, E_rot, E_vib, E_tot, TS_trans, TS_rot, TS_vib, TS_el, vibenergycorr, Hcorr, Gcorr, TS_tot)

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


    #Numfreq job. A 1-point or 2-point Hessian can be requested.
    # Either serial or parallell runmode can be used.
    # For parallel: Ash will use the number of cores given to run same number of displacments simultaneouslyu.
    #ORCA parallelization is turned off automatically.

    #Serial mode:
    freqresult = NumFreq(fragment=Reactant, theory=ORCAcalc, npoint=2, runmode='serial')
    #Parallel mode:
    freqresult = NumFreq(fragment=Reactant, theory=ORCAcalc, npoint=2, runmode='parallel', numcores=Ashnumcores)

    print("freqresult:", freqresult)
    #The resulting object from a NumFreq calculation is a dictionary (here called freqresult)
    # It contains the calculated frequencies and results from the Thermochemical analysis.
    #Individual items from the dictionary can be accessed by specifying the dictionary key:
    # Available keys: frequencies, ZPVE, vibenergy, transenergy, rotenergy, vibenergy, vibenergycorr
    # TO BE FINISHED...
    print("Frequencies : ", freqresult['frequencies])
    print("ZPVE : ", freqresult['ZPVE])


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
Surface scans
###########################
Potential Energy Surfaces can be conveniently scanned in ASH using the **calc_surface function** that uses the **geometric** optimization library.
Both unrelaxed and relaxed scans be calculated, using either 1 and 2 reaction coordinates.

The calc_surface function takes a fragment object and theory object as input. The type of scan is specified ('Unrelaxed' or 'Relaxed') and
then either 1 or 2 reaction coordinates are specified via keyword arguments: RC1_type, RC1_range and RC1_indices (and RC2 versions if using two reaction coordinates).

- The RC1_type/RC2_type keyword can be: 'bond', 'angle' or 'dihedral'.
- The RC1_indices/RC2_indices keyword defines the atom indices for the bond/angle/dihedral. Note: ASH counts from zero.
- The RC1_range/RC2_range keyword species the start coordinate, end coordinate and the stepsize (Å units for bonds, ° for angles/dihedrals).

The resultfile keyword should be used to specify the name of the file that contains the results of the scan ( format: coord1 coord2 energy).
This file can be used to restart an incomplete/failed scan. If ASH finds this file in the same dir as the script, it will read the data and skip unneeded calculations.
Default name : 'surface_results.txt'


**calc_surface** returns a dictionary of total energies for each surface point. The key is a tuple of coordinate value and the value is the energy, i.e.
(RC1value,RC2value) : energy

**1D scan:**

.. code-block:: python

    surfacedictionary = calc_surface(fragment=frag, theory=ORCAcalc, type='Unrelaxed', resultfile='surface_results.txt', runmode='serial',
        RC1_range=[180,110,-10], RC1_type='angle', RC1_indices=[1,0,2])

**2D scan:**

If both RC1 and RC2 keywords are provided then a 2D scan will be calculated.

.. code-block:: python

    surfacedictionary = calc_surface(fragment=frag, theory=ORCAcalc, type='Unrelaxed', resultfile='surface_results.txt', runmode='serial',
        RC1_type='bond', RC1_range=[2.0,2.2,0.01], RC1_indices=[[0,1],[0,2]], RC2_range=[180,110,-10], RC2_type='angle', RC2_indices=[1,0,2])

NOTE: It is possible to have each chosen reaction coordinate apply to multiple sets of atom indices by specifying a list of lists.
In the 2D scan example above, the RC1_indices keyword (a 'bond' reaction coordinate) will apply to both atoms [0,1] as well as [0,2].
This makes sense when preserving symmetry of a system e.g. the O-H bonds in H2O.

NOTE: Currently the runmode is serial which means that one surface point is run after the other and only the theory level can be parallelized.
A future parallel runmode will become available where X surfacepoints can be run simultaneously using X available cores.

Other options to calc_surface:

- coordsystem  (for geomeTRICOptimizer, default: 'dlc'. Other options: 'hdlc' and 'tric')
- maxiter (for geomeTRICOptimizer,default : 50)
- extraconstraints (for geomeTRICOptimizer, default : None. dictionary of additional constraints. Same syntax as constraints in **geomeTRICOptimizer**)
- convergence_setting (for geomeTRICOptimizer, same syntax as in **geomeTRICOptimizer**)

**Working with a previous scan from collection of XYZ files**

If a surface scan has already been performed, it's possible to use the created XYZ-files and calculate energies for each surfacepoint with
e.g. a high-level of theory (CCSD(T) for instance).

We can use the **calc_surface_fromXYZ** function to read in previous XYZ-files (named like this: RC1_2.0-RC2_180.0.xyz for a 2D scan and like this: RC1_2.0.xyz for a 1D scan).
These files should have been created from **calc_surface** already (present in surface_xyzfiles results directory).
By providing a theory level object we can then easily perform single-point calculations for each surface point.
The results is a dictionary like before.

.. code-block:: python

    #Directory of XYZ files. Can be full path or relative path.
    surfacedir = '/users/home/ragnarbj/Fe2S2Cl4/PES/Relaxed-Scan-test1/SP-DLPNOCC/surface_xyzfiles'

    #Calculate surface from collection of XYZ files. Will read old surface-results.txt file if requested (resultfile="surface-results.txt")
    surfacedictionary = calc_surface_fromXYZ(xyzdir=surfacedir, theory=ORCAcalc, dimension=2, resultfile='surface_results.txt' )



**Plotting**

The final result of the scan is stored in a dictionary (named 'surfacedictionary' in the examples above) and can be easily
plotted by giving the dictionary as input to plotting functions (based on Matplotlib).
See :doc:`plotting`) page.

The dictionary has the format: (coord1,coord2) : energy  for a 2D scan  and (coord1) : energy for a 1D scan
where (coord1,coord2)/(coord1) is a tuple of floats and energy is the total energy as a float.

A dictionary using data from a previous job (stored e.g. in surface_results.txt) can be created via the **read_surfacedict_from_file** function:

.. code-block:: python

    surfacedictionary = read_surfacedict_from_file("surface_results.txt", dimension=1)



###########################
Saddle-point optimization
###########################

Not yet ready



###########################
Molecular Dynamics
###########################

Not yet ready