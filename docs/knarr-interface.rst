KNARR interface
======================================

Nudged elastic band (NEB) calculations with/without a climbing image (CI-NEB) are possible in ASH using any theory level via an interface to the KNARR code.

.. code-block:: python

    def NEB(reactant=None, product=None, theory=None, images=None, interpolation=None, CI=None, free_end=None, restart_file=None,
            conv_type=None, tol_scale=None, tol_max_fci=None, tol_rms_fci=None, tol_max_f=None, tol_rms_f=None,
            tol_turn_on_ci=None, ActiveRegion=False, actatoms=None, runmode='serial', printlevel=0,
            idpp_maxiter=None, charge=None, mult=None):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``reactant``
     - ASH Fragment
     - None
     - An ASH fragment for the reactant. Geometry should ideally be previously optimized.
   * - ``product``
     - ASH Fragment
     - None
     - An ASH fragment for the product. Geometry should ideally be previously optimized.
   * - ``theory``
     - ASH THeory
     - None
     - An ASH Theory. This will be used to supply energy+gradient information for each image during the NEB.
   * - ``images``
     - integer
     - None
     - The number of images to use during the NEB.
   * - ``interpolation``
     - string
     - None
     - The type of interpolation used. Default: 'IDPP'. Other options: 'linear'
   * - ``CI``
     - Boolean
     - None
     - Whether to use a climbing image or not during NEB.
   * - ``free_end``
     - Boolean
     - None
     - Whether initial images are free or frozen.
   * - ``restart_file``
     - 'string'
     - None
     - Name of restart-file (XYZ-file) to use as initial guess path instead of IDPP/linear.
   * - ``conv_type``
     - string
     - None
     - Convergence type. Default: 'all'
   * - ``tol_scale``
     - float
     - None
     - Scaling tolerance during NEB.
   * - ``tol_max_fci``
     - float
     - None
     - Max tolerance for when the climbing image is converged.
   * - ``tol_rms_fci``
     - float
     - None
     - RMS tolerance for when the climbing image is converged.
   * - ``tol_max_f``
     - float
     - None
     - Max tolerance for when images are converged (all except CI.)
   * - ``tol_rms_f``
     - float
     - None
     - RMS tolerance for when images are converged (all except CI.)
   * - ``tol_turn_on_ci``
     - float
     - None
     - The threshold for turning on the climbing image.
   * - ``ActiveRegion``
     - Boolean
     - False
     - Whether to use an Active Region during the NEB job. This requires setting the number of active atoms (actatoms list).
   * - ``actatoms``
     - list
     - None
     - List of atom indices that are active during the NEB job. All other atoms are frozen. 
   * - ``runmode``
     - string
     - 'serial'
     - Runmode: 'serial' or 'parallel'. 'parallel' is not yet available.
   * - ``printlevel``
     - integer
     - 0
     - The printlevel to use in NEB calculations. This printlevel is enforced in the Theory object.
   * - ``idpp_maxiter``
     - integer
     - None
     - Maximum number of iterations to use in the IDPP guess path calculations
   * - ``charge``
     - integer
     - None
     - Optional specification of the charge of the system (if QM).
   * - ``mult``
     - integer
     - None
     - Optional specification of the spin multiplicity of the system (if QM).

################################################################################
How to use
################################################################################

- Generally, the default settings are quite reliable and should not be touched.
- CI-NEB is almost always desired.
- The most important parameter is the number of images. Few images will result in fewer energy+gradient calculations but may hinder convergence. Many images will give a high resolution of the minimum energy path, will have fewer convergence problems but will result in an expensive calculation. A good number is usually 5-11 or so.
- When a partially converged NEB path reveals that there probably is an intermediate inbetween, it is best to cancel the calculation and split the job into 2 jobs, i.e. start a new job from reactant to intermediate and another from intermediate to product.
- It can be a good idea to do an initial NEB from a lower level of theory (e.g. xTB) before doing the higher level of theory (DFT). Use restart_file option to read in lower-level MEP as guess.
- If a CI-NEB calculation converges, the saddlepoint geometry can be confirmed as a saddlepoint via a NumFreq job.


################################################################################
Examples
################################################################################

10-image NEB calculation at the XTB level of theory

.. code-block:: python

    from ash import *

    numcores=2

    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    prod=Fragment(xyzfile="prod.xyz", charge=0, mult=1)


    #Theory to use for NEB
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library')

    #Run NEB to find saddlepoint. Returns saddlepoint as ASH fragment
    SP = NEB(reactant=react, product=prod, theory=xtbcalc, images=10)

    #Optional NumFreq job on saddlepoint to confirm that a saddlepoint was found.
    NumFreq(theory=xtbcalc, fragment=SP)


Restarting a calculation with user-defined path-file. 
Here, using the *restart_file* option to the NEB we read in a previous Knarr path-file ("knarr_MEP.xyz") instead of doing the regular IDPP interpolation
This file must contain the coordinates of the same number of images (here 10) as number of images specified.
The file can come from a previously unconverged NEB calculation or perhaps a converged MEP from a calculation at another level of theory.

.. code-block:: python

    from ash import *

    numcores=4

    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    prod=Fragment(xyzfile="prod.xyz", charge=0, mult=1)


    #Theory to use for NEB
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library', numcores=numcores)

    #Run NEB to find saddlepoint. Returns saddlepoint as ASH fragment
    SP = NEB(reactant=react, product=prod, theory=xtbcalc, images=10, restart_file="knarr_MEP.xyz")


################################################################################
Controlling printout
################################################################################

During a NEB calculation the theory code is called multiple times to calculate the energy and gradient.
As the printout can become excessive (especially if using a QMMMTheory object) it is usually desirable to reduce printout considerably for NEB calculations.
Currently this is performed by default by setting the printlevel in the NEB calculation as a keyword argument.
The NEB printlevel is then used to set the printlevel in the Theory objects.
The default printlevel is 0 (barely any output from other modules) but this can be increased to 1,2 or 3 to get more output.

Example:

.. code-block:: python

    #Run NEB to find saddlepoint. Returns saddlepoint as ASH fragment
    SP = NEB(reactant=react, product=prod, theory=xtbcalc, images=10, printlevel=1)


################################################################################
Parallelization
################################################################################

Currently, NEB calculations can only be parallelized via the Theory object, i.e. in the example above, xTB utilizes 4 CPU cores each time it is called.
In the future it will be possible to run multiple NEB images in parallel (more efficient parallelization).