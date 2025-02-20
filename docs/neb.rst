Nudged Elastic Band
======================================


Nudged elastic band (NEB) calculations with/without a climbing image (CI-NEB) are possible in ASH using any theory level via an interface to the state-of-the-art NEB code, KNARR code, Vilhjálmur Ásgeirsson.
This NEB code uses energy-weighted springs (instead of a single spring constant for all images), together with a global L-BFGS optimizer that results in a greatly improved NEB algorithm for molecular systems.

The Knarr implementation predates the ORCA implementation but is overall very similar to the one described there:

V. Ásgeirsson, B. Birgisson, R. Bjornsson, U. Becker, F. Neese, C: Riplinger,  H. Jónsson, J. Chem. Theory Comput. 2021,17, 4929–4945.
DOI: 10.1021/acs.jctc.1c00462

A NEB-TS implementation is also available. 


.. note:: While QM/MM NEB calculations are possible, running NEB in parallel with a QM/MM Hamiltonian is currently not possible. Only Theory-parallelization can be used.


.. code-block:: python

  def NEB(reactant=None, product=None, theory=None, images=8, CI=True, free_end=False, 
          conv_type="ALL", tol_scale=10, tol_max_fci=0.026, tol_rms_fci=0.013, tol_max_f=0.26, tol_rms_f=0.13,
          tol_turn_on_ci=1.0,  runmode='serial', numcores=1, 
          charge=None, mult=None,printlevel=0, ActiveRegion=False, actatoms=None,
          interpolation="IDPP", idpp_maxiter=300, 
          restart_file=None, TS_guess_file=None, mofilesdir=None):

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
     - 8
     - | The number of images to use during the NEB. For free_end=False (recommended) this means that this 
       | number of intermediate images will be active and calculated. For free_end=True, the reactant and
       | product images will also be included i.e. total_images = images+2 
   * - ``interpolation``
     - string
     - 'IDPP'
     - The type of interpolation used. Default: 'IDPP'. Other options: 'linear', 'geodesic'
   * - ``CI``
     - Boolean
     - True
     - Whether to use a climbing image or not during NEB. CI is always recommended.
   * - ``free_end``
     - Boolean
     - False
     - | Whether initial images are free or frozen. Warning: the number of active images will be images+2 
       | when this is True.
   * - ``restart_file``
     - 'string'
     - None
     - Name of restart-file (XYZ-file) to use as initial guess path instead of IDPP/linear.
   * - ``TS_guess_file``
     - 'string'
     - None
     - | Name of XYZ-file containing a guess for the saddlepoint("TS"). The TS guess will be used during the
       | interpolation step to provide a better guess for the overall guess-path.
   * - ``mofilesdir``
     - 'string'
     - None
     - | Path (absolute or relative path) to a directory with MO-files to be used as guess orbitals in the first
       | NEB iterations. For ORCATheory these should be GBW-files named: "current_image0.gbw" etc.
   * - ``conv_type``
     - string
     - 'ALL'
     - Convergence type. Options: 'ALL', 'CIONLY'. 
   * - ``tol_scale``
     - float
     - 10
     - Scaling tolerance during NEB. Units: eV/Å
   * - ``tol_max_fci``
     - float
     - 0.026
     - The Max Force convergence tolerance for the Climbing IMage. Units: eV/Å
   * - ``tol_rms_fci``
     - float
     - 0.013
     - The RMS Force convergence tolerance for the Climbing IMage. Units: eV/Å
   * - ``tol_max_f``
     - float
     - 0.26
     - The Max Force convergence tolerance for all other images except CI. Units: eV/Å
   * - ``tol_rms_f``
     - float
     - 0.13
     - The RMS Force convergence tolerance for all other images except CI. Units: eV/Å
   * - ``tol_turn_on_ci``
     - float
     - 1.0
     - The threshold on Max Force for turning on the climbing image. Units: eV/Å
   * - ``ActiveRegion``
     - Boolean
     - False
     - | Whether to use an Active Region during the NEB job. This requires setting the
       | number of active atoms (actatoms list) below.
   * - ``actatoms``
     - list
     - None
     - List of atom indices that are active during the NEB job. All other atoms are frozen. 
   * - ``runmode``
     - string
     - 'serial'
     - Runmode: 'serial' or 'parallel'. Runmode 'parallel' requires setting numcores keyword below.
   * - ``numcores``
     - integer
     - 1
     - | Numcores of cores to use in a parallel NEB job. 
       | Recommended: use same number of cores as there are active images.
   * - ``printlevel``
     - integer
     - 0
     - | The printlevel to use in NEB calculations. This same printlevel is enforced in the Theory object. Change to
       | 1,2 or 3 to get more printing.
   * - ``idpp_maxiter``
     - integer
     - 300
     - Maximum number of iterations to use in the IDPP guess path calculations
   * - ``idpp_springconst``
     - float
     - 5.0
     - Value of the spring constant used in IDPP interpolation. Units: eV/Å^2
   * - ``charge``
     - integer
     - None
     - Optional specification of the charge of the system (if QM).
   * - ``mult``
     - integer
     - None
     - Optional specification of the spin multiplicity of the system (if QM).

################################################################################
Recommendations and how to use
################################################################################

- Generally, the default settings are reasonably reliable and need not be touched. The number of images is the most important parameter.
- CI-NEB is almost always desired. A NEB job without CI, gives only a crude estimate of the barrier.
- Frozen endpoints (free_end=False, i.e. reactant and product frozen) is usually recommended (more efficient). Reactant and product should have been previously optimized.
- If you want to explore active endpoints (free_end=True), note that the number of active images will be images+2,
- The most important parameter is the number of images. Few images will result in fewer energy+gradient calculations but may hinder convergence. Many images will give a high resolution of the minimum energy path, will have fewer convergence problems but will result in an expensive calculation. A good number is usually 5-11 or so.
- In runmode = 'parallel' you should generally choose the number of active images to be equal to the number of CPU cores provided to NEB.
- If you activate parallelization of the theory level also, this will be the number of cores used per image. So if you do ORCATheory(...numcores=2) and NEB(...images=8,numcores=8) ASH will be attempting to use 2x8 = 16 cores. 
- When a partially converged NEB path reveals that there probably is an intermediate inbetween, it is best to cancel the calculation and split the job into 2 jobs, i.e. start a new job from reactant to intermediate and another from intermediate to product. A CI-NEB job would only converge to the higher energy saddlepoint in such a case.
- It can be a good idea to do an initial NEB from a lower level of theory (e.g. xTB) before doing the higher level of theory (DFT). Use restart_file option to read in lower-level MEP as guess.
- If you already know approximately what the saddlepoint geometry should look like you can provide such a geometry using the TS_guess_file option. The geometry will be used during the interpolation to provide a more accurate guess path. This could also be a previously obtained saddlepoint at another level of theory.
- In rare cases the IDPP interpolation goes wrong, you can either 1) try modify the idpp_springconst value or 2) switch to a simpler linear Cartesian interpolation (interpolation="linear" option) instead, perhaps in combination with a TS_guess_file (guides the linear interpolation).
- There is now also the option of using the 'GEODESIC' option which uses the geodesic_interpolate library to perform the interpolation.
- If a CI-NEB calculation converges, the saddlepoint geometry can be confirmed as a saddlepoint via a NumFreq job. NEB returns an ASH Fragment inside the ASH-Results object (saddlepoint_fragment attibute) of the saddlepoint geometry as well as an XYZ-file.
- Any ASH Theory level can in principle be used (although only ORCA and xTB have been tested). In practice you want to use a QM method and code with an analytical gradient available.


################################################################################
Examples
################################################################################

**8-image NEB calculation at the XTB level of theory (Theory parallelization):**

.. code-block:: python

    from ash import *

    numcores=8

    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    prod=Fragment(xyzfile="prod.xyz", charge=0, mult=1)

    #Theory to use for NEB. Setting number of cores for xTB.
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library', numcores=numcores)

    #Run NEB to find saddlepoint. Returns ASH Results object.
    #Note: the saddlepoint fragment
    NEB_result = NEB(reactant=react, product=prod, theory=xtbcalc, images=8)
    print(NEB_result) #Printout of the NEB_result object

    #Optional NumFreq job on saddlepoint to confirm that a saddlepoint was found.
    NumFreq(theory=xtbcalc, fragment=NEB_result.saddlepoint_fragment)


**Restarting a calculation with user-defined path-file.**

Here, using the *restart_file* option to the NEB we read in a previous Knarr path-file ("knarr_MEP.xyz") instead of doing the regular IDPP interpolation
This file must contain the coordinates of the same number of images (here 10) as number of images specified.
The file can come from a previously unconverged NEB calculation or perhaps a converged MEP from a calculation at another level of theory.

.. code-block:: python

    from ash import *

    numcores=1

    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    prod=Fragment(xyzfile="prod.xyz", charge=0, mult=1)


    #Theory to use for NEB
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library')

    #Run NEB to find saddlepoint. Returns an ASH Results object 
    NEB_result = NEB(reactant=react, product=prod, theory=xtbcalc, images=10, restart_file="knarr_MEP.xyz")


**A calculation with user-defined guess for the saddlepoint.**


Here, using the *TS_guess_file* option. This will influence the initial interpolation path generation by interpolating between reactant and guess_TS structure and guess_TS structure and product.

.. code-block:: python

    from ash import *

    numcores=1

    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    prod=Fragment(xyzfile="prod.xyz", charge=0, mult=1)

    #Theory to use for NEB
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library')

    #Run NEB to find saddlepoint. Returns an ASH Results object 
    NEB_result = NEB(reactant=react, product=prod, theory=xtbcalc, images=10, TS_guess_file="guess_TS_geometry.xyz")
    print(NEB_result)

################################################################################
Guess pathway and Interpolation
################################################################################

The initial guess pathway plays an important role in NEB calculations.
If you end up with NEB convergence problems that are never resolved or perhaps even SCF convergence problems in the very first NEB iteration,
it is likely that there is something wrong with the initial guess pathway.
Visualizing the guess pathway, present in the file initial_guess_path.xyz may reveal the problem.
A common issue is that the reactant and product geometries do not have atoms ordered in a consistent way which will lead to a problematic pathway.

However, it is also possible that the default IDPP interpolation fails to produce a good enough pathway for your system.
The problem can potentially be fixed by tweaking the idpp_maxiter (default value 700 )and idpp_springconst (default 5.0) parameters.
But there is also an alternative guess-option in ASH now, the 'GEODESIC' option which is based on geodesic interpolation by Todd Martínez and coworkers.
The algorithm is described in :
Xiaolei Zhu, Keiran C. Thompson, Todd J. Martínez, J. Chem. Phys. 2019, 150, 164103. `Article <https://pubs.aip.org/aip/jcp/article/150/16/164103/198363/Geodesic-interpolation-for-reaction-pathways>`_ 

In initial tests GEODESIC seems to improve upon IDPP for molecular reactions and maybe become the default in NEB and NEBTS jobs in ASH in the future

Use like this:

.. code-block:: python

    NEB_result = NEB(reactant=react, product=prod, theory=xtbcalc, images=10, interpolation="GEODESIC")

################################################################################
Controlling printout
################################################################################

During a NEB calculation the theory code is called multiple times to calculate the energy and gradient.
As the printout can become excessive (especially if using a QMMMTheory object) it is usually desirable to reduce printout considerably for NEB calculations.
This is performed by setting the printlevel in the NEB calculation as a keyword argument.
The NEB printlevel is then used to set the printlevel in the Theory objects.
The default printlevel is 0 (barely any output from other modules) but this can be increased to 1,2 or 3 to get more output, both from the NEB function and the Theory level etc.
Printlevel 1 is useful for getting more useful information from the NEB module (especially regarding guess orbital logic) as well as slightly more information from the Theory object. Printlevel 2 will probably give too much output in general but can be useful for diagnostic purposes.

Example:

.. code-block:: python

    #Run NEB to find saddlepoint. Returns an ASH Results object (NEB_result.saddlepoint_fragment is the saddlepoint ASH Fragment).
    NEB_result = NEB(reactant=react, product=prod, theory=xtbcalc, images=10, printlevel=1)
    print(NEB_result)

################################################################################
Controlling guess orbitals during SCF of Theory level
################################################################################

During the NEB job the Theory level object is called multiple times using each iteration. The Theory level object will handle what guess orbitals are used during this step and you can modify the Theory object as desired 
(e.g. for ORCATheory you can change autostart and moreadfile keywords as desired).

For a default NEB calculation in runmode='serial':
for e.g. ORCATheory, the first calculation in the NEB job (NEB iteration -1) will be on the reactant. ORCA will in this case use brand-new guess orbitals (from PModel guess typically). 
Once converged, the orca.gbw file will be copied and stored as current_image0.gbw by the NEB module.
Next calculation on the product will use the previous orca.gbw file (from reactant) since ORCA will by default try to read orbitals from that file (since the inputfile has the same basename) but once ORCA is finished we will store the file as e.g. current_image11.gbw
This is repeated for intermediate images: image1, image2, ..., image10.gbw in NEB iteration 0.
However, in the next NEB iterations, the code will find and use current_image1.gbw for image no. 1 etc. since these files have been stored. These files will be updated during the job, ensuring that each image has converged image-specific orbitals from the last iteration available.
In order to see detailed printout for what is going on w.r.t. ORCA GBW-file book-keeping during the NEB job, set the NEB printlevel to 1. 

For parallel NEB calculations with ORCATheory, things are just slightly different as there will be different directories for each Python multiprocessing worker, called e.g. 'Pooljob_image_9'. Orbitals inside file 'orca.gbw' from last NEB iteration for that image will be read each time.

Generally this behaviour works well as previously converged orbitals, specific to each image are being used.
If you require even more control over which orbitals should be used there are a few options.

**1. Reading in a single initial guess orbital-file (ORCATheory) :**

By doing ORCATheory(...moreadfile="test.gbw"), ORCA will read in orbitals from file "test.gbw" (make sure to copy file "test.gbw" to scratch or provide full path) in the first calculation by NEB(reactant calculation).
This option is primarily useful if the system is tough to converge (e.g. a BS-DFT job on a spin-coupled system).
Note: By default, the moreadfile option is turned off in the ORCATheory object after that so if you want to enforce moreadfile behaviour for every calculation during the NEB job, you could do: ORCATheory(...moreadfile="test.gbw", moreadfile_always=True).
This is probably unlikely to be useful though.


.. code-block:: python

  from ash import *

  numcores=1
  #SN2 reaction
  Reactant=Fragment(xyzfile="react.xyz", charge=-1, mult=1)
  Product=Fragment(xyzfile="prod.xyz",charge=-1, mult=1)

  #Calculator object without frag
  calc = ORCATheory(orcasimpleinput="!r2scan-3c tightscf CPCM", numcores=numcores, moreadfile="test.gbw")

  NEB_result = NEB(reactant=Reactant, product=Product, theory=calc, images=10, printlevel=0)
  print(NEB_result)

**2. Reading in guess orbitals for each image separately from a directory (ORCATheory):**

A better way to control the original guess is to provide to the NEB function, a mofilesdir keyword pointing to a directory-path that contains GBW files for each image.
The directory should contain GBW files for each image and should be called: 

*current_image0.gbw, current_image1.gbw, current_image2.gbw, ..., current_image11.gbw* etc.

This allows you more flexibility in choosing precisely what orbitals will be read in initially.

**Note:** Orbitals will only be read from the mofilesdir directory in NEB-iteration -1 (first reactant and product calcs) and NEB-iteration 0 (first intermediate image calculations). In the subsequent NEB iterations, the program will use image-specific GBW files from the previous iteration.

**Note:**  The mofilesdir path must either be a full path to a directory that is available to the computing node (e.g. /home/bjornsson/NEBjob1/mofilesdir or something) or a directory that is copied over to the the scratch
directory by your job-submission script.

.. code-block:: python

  from ash import *

  numcores=1
  #SN2 reaction
  Reactant=Fragment(xyzfile="react.xyz", charge=-1, mult=1)
  Product=Fragment(xyzfile="prod.xyz",charge=-1, mult=1)

  #Calculator object without frag
  calc = ORCATheory(orcasimpleinput="!r2scan-3c tightscf CPCM", numcores=numcores)

  NEB_result = NEB(reactant=Reactant, product=Product, theory=calc, images=10, printlevel=0, mofilesdir="/home/bjornsson/NEBjob1/mofiles_dir")
  print(NEB_result)

################################################################################
Controlling convergence
################################################################################

NEB convergence is controlled by a number of thresholds. Note that Knarr internally utilizes units of Å (distances and coordinates), eV (energies), eV/Å (forces), eV/Å^2 (force constants).
For now, the interface requires you to specify convergence tolerances in these units as well.

conv_type: 'ALL' or 'CIONLY' options specifies whether the NEB job should end when all the tolerances of the images have been met ('ALL') or only on the CI ('CIONLY')
The default is 'ALL' and is recommended. All 4 threshold belows have to be met in this case (only the first 2 in the case of 'CIONLY').

**Convergence tolerances:**

+------------------+---------------+-------------------------------------------------+
| **Tolerance**    | **Default**   | **Description**                                 |
+------------------+---------------+-------------------------------------------------+
| tol_max_fci      | 0.026 eV/Å    | when Max Force on the CI is met.                |
+------------------+---------------+-------------------------------------------------+
| tol_rms_fci      | 0.013 eV/Å    | when RMS Force on the CI is met                 |           
+------------------+---------------+-------------------------------------------------+
| tol_max_f        | 0.26 eV/Å     | when Max Force on all other images is met.      |
+------------------+---------------+-------------------------------------------------+
| tol_rms_f        | 0.13 eV/Å     | when RMS Force on all other images is met.      |           
+------------------+---------------+-------------------------------------------------+

**Other thresholds:**


+------------------+---------------+------------------------------------------------------+
| **Tolerance**    | **Default**   | **Description**                                      |
+------------------+---------------+------------------------------------------------------+
| tol_turn_on_ci   | 1.0 eV/Å      |  Specifies at which MaxF value, the CI is turned on  |
+------------------+---------------+------------------------------------------------------+
| tol_scale        | 10            |                                                      |           
+------------------+---------------+------------------------------------------------------+

################################################################################
Free-end NEB calculations
################################################################################

A recommended NEB job has endpoints (reactant and product) previously optimized at the same level of theory and are then kept frozen during the NEB job.
This usually results in a more efficient NEB job as it constrains the possibilities for the minimum energy path and saddlepoint search.

A free_end = True option where the endpoints are also minimized during the NEB is also possible but as there are more degrees of freedom, it can be trickier to converge.
This may be a good option when the endpoints have deliberately not been minimized in an effort to explore multiple potential reaction pathways.


################################################################################
NEB on systems with an active region (e.g. QM/MM)
################################################################################

For large systems, e.g. a QM/MM model of a protein active site, it is possible to perform a NEB calculation of only a selected group of atoms, with other atoms being frozen during the NEB iterations.
You should set ActiveRegion=True in this case and then specify the list of active atoms by their indices via e.g. actatoms=[17,18,19,20,21,22,23,24]
As a NEB calculation is a difficult minimization problem it is advised to keep the active region as small as possible, at least to begin with. For a QM/MM job it might be a good idea to first set actatoms = qmatoms. i.e. only allow the QM atoms to move during the NEB path minimization.
Future version of the code may further allow one to use weights 

Note: When an active region is used, the RMSD minimization for images is turned off automatically (used to superimpose images to avoid complicated NEB paths).

################################################################################
Parallelization
################################################################################

During each NEB iteration, X number of images are active and their energy+gradient needs to be calculated for each new geometry in each iteration.
As each E+G image calculation is independent from the others it is possible to utilize parallelization very effectively in a NEB job.
It is generally recommended to prioritize parallelization over images rather than the Theory level (QM parallelization never scales perfectly)
Theory parallelization is also available, however, and can be used to further speed up NEB job.

NEB-parallelization with a QM/MM Hamiltonian is currently not possible due to problems with the multiprocessing library and OpenMM.

**Example: 8-image NEB calculation at the XTB level of theory (NEB parallelization):**

If you are calculating 8 images then you should set runmode='parallel' and use numcores=8.

.. code-block:: python

    from ash import *

    numcores=8
    numimages=numcores
    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    prod=Fragment(xyzfile="prod.xyz", charge=0, mult=1)

    #Theory to use for NEB
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library', numcores=numcores)

    #Run NEB to find saddlepoint. Returns an ASH Results object
    NEB_result = NEB(reactant=react, product=prod, theory=xtbcalc, images=numimages, runmode='parallel', numcores=numcores)
    print(NEB_result)

    #Optional NumFreq job on saddlepoint to confirm that a saddlepoint was found.
    NumFreq(theory=xtbcalc, fragment=NEB_result.saddlepoint_fragment)

If you have additional CPU cores available on your computing node that you would like to use to speed up an NEB job you have 2 options:

- You could increase the number of images as well as CPU cores to e.g. 16. Such a 16-image/core-job would run each NEB iteration at the same speed as the 8 image/core job but since there are more images it may facilitate convergence and locate the saddlepoint more efficiently.
- Or you could active Theory parallelization by setting the numcores keyword for the Theory level. If you have 16 cores available on your node, you could set Theory parallelization to 2 which would result in each of the 8 images utilizing 2 CPU cores to speed up the E+G step, resulting in 16 cores being used. Note that if the Theory parallelization utilizes MPI it is possible that problems could occur.


**Example: 16-core job using 8-image NEB parallelization + Theory parallelization:**

This NEB job would run 8 active images simultaneously (via Python multiprocessing library) while parallelizing each xTB E+G calculation by 2 cores.
This job requires 16 available CPU cores.

.. code-block:: python

    from ash import *

    numcores=16 #Total number of CPU cores to be used (makes sure to submit a job with this number of slots)
    numimages=8 #Number of images in NEB job and the number of cores available to the NEB parallelization
    cores_theory=numcores/numimages #Number of cores used to parallelize the Theory level
    ################################################
    # Defining reactant and product ASH fragments
    #################################################
    react=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    prod=Fragment(xyzfile="prod.xyz", charge=0, mult=1)

    #Theory to use for NEB. Parallelizing
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library', numcores=cores_theory)

    #Run NEB to find saddlepoint. Returns an ASH Results object
    NEB_result = NEB(reactant=react, product=prod, theory=xtbcalc, images=numimages, runmode='parallel', numcores=numimages)
    print(NEB_result)

    #Optional NumFreq job on the NEB saddlepoint to confirm that a saddlepoint was found.
    NumFreq(theory=xtbcalc, fragment=NEB_result.saddlepoint_fragment)



################################################################################
NEB-TS : combining CI-NEB with TS-optimization
################################################################################

As discussed in the article:

V. Ásgeirsson, B. Birgisson, R. Bjornsson, U. Becker, F. Neese, C: Riplinger,  H. Jónsson, J. Chem. Theory Comput. 2021,17, 4929–4945.
DOI: 10.1021/acs.jctc.1c00462

a CI-NEB calculation is well suited to be combined with an eigenvector-following method for improved efficiency of a saddlepoint search.
The idea is to only partially converge a minimum energy path and saddlepoint via the CI-NEB method (that requires multiple images and a more complicated minimization)
but then use the approximate saddlepoint geometry to start an eigenvector-following optimization which can both make the overall saddlepoint search more efficient (as only a single image is calculated in the latter part) 
but can also ensure that a proper 1st-order saddlepoint is located via the use of exact/approximate Hessian information.

In a NEB-TS job in ASH, the Knarr library is used to perform a CI-NEB calculation while the geomeTRIC library is used to perform the eigenvector-following optimization. 

The NEBTS function is very similar to the NEB function:

.. code-block:: python

  def NEBTS(reactant=None, product=None, theory=None, images=8, CI=True, OptTS=True, free_end=False, maxiter=100,
          conv_type="ALL", tol_scale=10, tol_max_fci=0.10, tol_rms_fci=0.05, tol_max_f=1.03, tol_rms_f=0.51,
          tol_turn_on_ci=1.0,  runmode='serial', numcores=1, charge=None, mult=None, printlevel=0, ActiveRegion=False, actatoms=None,
          interpolation="IDPP", idpp_maxiter=300, restart_file=None, TS_guess_file=None, mofilesdir=None, 
          OptTS_maxiter=100, OptTS_print_atoms_list=None, OptTS_convergence_setting=None, OptTS_conv_criteria=None, OptTS_coordsystem='tric',
          hessian_for_TS=None, modelhessian='unit', tsmode_tangent_threshold=0.1):

with additional keywords: *OptTS_maxiter*, *OptTS_print_atoms_list*, *OptTS_convergence_setting*, *OptTS_conv_criteria* and *OptTS_coordsystem*  being keywords that belong to the Optimizer.
See :doc:`Geometry-optimization` for explanations.

An important option is the *hessian_for_TS* keyword which controls what type of Hessian should be used during the OptTS job.

Options to *hessian_for_TS* are:

.. list-table::
   :widths: 15 60
   :header-rows: 1

   * - hessian_for_TS value
     - Description
   * - ``first``
     - Optimizer calculates exact Hessian in the first step of the OptTS procedure.
   * - ``each``
     - Optimizer calculates exact Hessian in each step of the OptTS procedure (expensive).
   * - ``xtb``
     - Calculate an exact Hessian but at the cheap GFN1-xTB level of theory.
   * - ``model``
     - | Calculate a model Hessian (default: *modelhessian* ='unit') to be used as approximation to the exact Hessian. Requires ORCA.
       | *modelhessian* options: 'unit', 'Almloef', 'Lindh', 'Schegel'  
   * - ``partial``
     - | Calculate a partial exact Hessian using only the atoms that contribute the most to approximate TS-mode (from CI-NEB job).
       | Use *tsmode_tangent_threshold* to control the size of the partial Hessian.
       | Rest is approximated by a model Hessian or unit atrix. *modelhessian* options: 'unit','Almloef', 'Lindh', 'Schegel'  

*hessian_for_TS* ='xtb' is the currently recommended option. This will do an xTB NumFreq calculation at the saddlepoint geometry and this Hessian will then be used
as an initial Hessian in the eigenvector-following minimization. Unless the system is very large, this option is the most cost-effective. 
This requires an active xTB interface (xTB needs to installed on the computer).
If this option fails: 'first' will calculate an exact Hessian in the first step. A safe but very expensive option is to use 'each' (exact Hessian in every Opt step).


**Example:**

.. code-block:: python

    from ash import *

    Reactant=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    Product=Fragment(xyzfile="prod.xyz", charge=0, mult=1)
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf") #ORCATheory object creation

    #NEB-TS combines a CI-NEB job (note: looser thresholds than default CI-NEB) and a Optimizer(OptTS=True) job.
    NEBTS_result = NEBTS(reactant=Reactant, product=Product, theory=calc, images=12, printlevel=0, hessian_for_TS='xtb')
    print(NEBTS_result)


Parallelization of a **NEBTS** job can be controlled by the *numcores* keyword and for the CI-NEB part it will behave like in the **NEB** function.
However, once the CI-NEB part is complete, and the NEBTS job switches to performing the eigenvector-following minimization, ASH will automatically
change the number of cores available to the Theory object to use the maximum number of CPU cores provided to either NEBTS or the Theory object. 
This maximizes use of CPU cores during the job.

**Parallelization example:**

.. code-block:: python

    from ash import *

    numcores=16 #Total number of CPU cores to be used by ASH. OptTS will later use all of these.
    numimages=8 #Number of images in NEB job and the number of cores available to the NEB parallelization.
    cores_theory=numcores/numimages #Number of cores used to parallelize the Theory level during NEB.

    Reactant=Fragment(xyzfile="react.xyz", charge=0, mult=1)
    Product=Fragment(xyzfile="prod.xyz", charge=0, mult=1)
    ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf", numcores=cores_theory) #ORCATheory object creation

    #NEB-TS combines a CI-NEB job (note: looser thresholds than default CI-NEB) and a Optimizer(OptTS=True) job.
    NEBTS_result = NEBTS(reactant=Reactant, product=Product, theory=calc, numcores=numimages, images=numimages, printlevel=0, hessian_for_TS='xtb')
    print(NEBTS_result)