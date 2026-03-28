DL-FIND Optimizer
======================================

ASH also features an interface to the powerful `DL-FIND library <https://www.itheoc.uni-stuttgart.de/research/kaestner/research/dlfind/>`_ .
DL-FIND is developed by Prof. Johannes Kaestner and includes various powerful optimization algorithms.
If you use the interface, make sure to cite the `DL-FIND article <https://pubs.acs.org/doi/10.1021/jp9028968>`_
The ASH interface to DL-FIND utilizes the `libdlfind C/Python API <https://github.com/digital-chemistry-laboratory/libdlfind>`_ .

The interface supports most of the algorithms available in DL-FIND including :

- Geometry optimizations using Cartesian, HDLC internal coordinates using L-BFGS algorithm.
- Saddlepoint optimizations using P-RFO algorithm
- Dimer-method for saddle-point searches requiring only gradients
- Nudged elastic band calculations
- Instanton calculations
- Optimizations with active/frozen atoms as well as bond, angle and dihedral constraints.

The DL-FIND Optimizer can be used in ASH surface scans, see :doc:`surfacescan`
Since the interface supports active regions, **DLFIND_optimizer** can e.g. be used for QM/MM geometry optimizations in ASH.

As DL-FIND is written in Fortran, the execution speed of the optimizer is very fast compared to geomeTRIC.
DLFIND_optimizer in ASH can hence be particulary useful if the execution speed of the energy+gradient by the Theory object
is very fast (e.g. an OpenMMTheory, a semiempirical xTB-theory, an ML potential) as the speed of the Optimizer will 
then start to play a role for the overall wall-time and DL-FIND execution can be very fast.


################################
Installation
################################

To install DL-FIND one needs the DL-FIND Fortran library as well as the C-API/Python-API extension, `libdlfind <https://github.com/digital-chemistry-laboratory/libdlfind>`_
It is easiest to install both using:

.. code-block:: shell

    pip install geometric

For problems associated with installation,
see: `DL-FIND <https://chemshell.org/dl-find/>`_ and `libdlfind Github repository <https://github.com/digital-chemistry-laboratory/libdlfind>`_



################################
The DLFIND_optimizer function
################################

.. code-block:: python

    def DLFIND_optimizer(jobtype=None, theory=None, fragment=None, fragment2=None, charge=None, mult=None, 
                        maxcycle=250, tolerance=4.5E-4, tolerance_e=1E-6,
                        actatoms=None, frozenatoms=None, residues=None, constraints=None,
                        printlevel=2, NumGrad=False, delta=0.01,
                        icoord=None, iopt=None, nimage=None, 
                        hessian_choice="numfreq", inithessian=0, 
                        numfreq_npoint=1, numfreq_displacement=0.005, numfreq_hessatoms=None,
                        numfreq_force_projection=None, print_atoms_list=None,
                        force_noPBC=False, PBC_format_option='CIF'):

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
   * - ``fragment2``
     - ASH Fragment
     - None
     - Some DL-FIND jobtypes require 2 Fragment objects.
   * - ``jobtype``
     - string
     - None
     - Type of DL-FIND job. Options: 'opt', 'tsopt', 'neb', 'dimer', 'qts'. 
   * - ``iopt``
     - integer
     - None
     - DLFIND-code: Alternative to jobtype, choose code for optimization algorithm, e.g. iopt=3 (L-BFGS optimization)
   * - ``icoord``
     - integer
     - None
     - DLFIND-code: Alternative to jobtype, choose code for job-type, e.g. icoord=120 (NEB with frozen endpoints)
   * - ``nimage``
     - integer
     - None
     - Number of images for NEB job.
   * - ``maxcycle``
     - integer
     - 250
     - Max number of cycles(iterations).
   * - ``tolerance``
     - float
     - 4.5E-4
     - Tolerance on gradient in Eh/Bohr
   * - ``tolerance_e``
     - float
     - 1e-6
     - Tolerance on gradient in Eh.
   * - ``actatoms``
     - list
     - None
     - List of indices that should be active (rest is frozen)
   * - ``frozenatoms``
     - list
     - None
     - List of indices that should be active (rest is active)
   * - ``residues``
     - list of lists
     - None
     - For HDLCs, it is possible to provide HDLC residue definitions as a list of lists.
   * - ``constraints``
     - dictionary
     - None
     - Dictionary of constraints. See section.
   * - ``NumGrad``
     - Boolean
     - False
     - Whether to perform a numerical gradient optimization (uses NumGradClas)
   * - ``delta``
     - float
     - 0.01
     - Delta parameter used for dimer method.
   * - ``inithessian``
     - integer
     - 0
     - For P-RFO TS optimizations, the DL-FIND code for Hessian option. Options: 0 (external), 1(1-point FD), 2(2-point FD), 3(diagonal 1-point FD), 4(identity matrix)
   * - ``hessian_choice``
     - string or Numpy array
     - numfreq
     - The input Hessian, only if inithessian=0. If "numfreq" then an ASH Numfreq Hessian is calculated. Options : 'anfreq', 'xtb', 'file:hessianfile' or a Numpy array.
   * - ``numfreq_npoint``
     - integer
     - 1
     - For a NumFreq Hessian, whether to use 1-point or 2-point Hessian (see Numfreq documentation)
   * - ``numfreq_displacement``
     - float
     - 0.005
     - For a NumFreq Hessian, displacement (in Bohrs).
   * - ``numfreq_hessatoms``
     - list
     - None
     - For a NumFreq Hessian, list of atoms to use in a partial Hessian calculation.
   * - ``numfreq_force_projection``
     - Boolean
     - None
     - For a NumFreq Hessian, whether to do translation-rotation projection or not.
   * - ``numfreq_force_projection``
     - Boolean
     - None
     - For a NumFreq Hessian, whether to do translation-rotation projection or not.
   * - ``print_atoms_list``
     - list
     - None
     - What atoms to use when printing coordinates in output.        
   * - ``force_noPBC``
     - Boolean
     - False
     - Whether to force PBCs to not be activated.   
   * - ``PBC_format_option``
     - string
     - 'CIF'
     - For a PBC optimization, what type of PBC fileformat to print in the end.   

######################################################
Controlling DL-FIND
######################################################

*Controlling DL-FIND via jobtype*

It is easiest to use **DLFIND_optimizer** together with the jobtype keyword, which will select recommended optimization algorithm and coordinate system.
The available options are: 

- 'opt'. Selects HDLC coordinates and L-BFGS optimizer.
- 'tsopt'. Selects HDLC coordinates and P-RFO optimizer. 
- 'neb'. Selects a frozen endpoint NEB job with a L-BFGS optimizer.
- 'dimer'. Selects a dimer job with an L-BFGS optimizer. 
- 'qts' (or 'instanton'). Selects instanton optimization with an L-BFGS optimizer.

*Controlling DL-FIND via jobtype*

Alternatively is also possible to control DL-FIND behaviour by not using a jobtype (jobtype=None)
and specify both icoord and iopt keywords.
Consult the DL-FIND manual or libdlfind README in this case:

See `DL-FIND manual <https://github.com/digital-chemistry-laboratory/libdlfind/blob/4167998d16d8dac4a484ba9305f27d6325a7a28d/docs/documentation.pdf>`_
See also `libdlfind README <https://github.com/digital-chemistry-laboratory/libdlfind/blob/4167998d16d8dac4a484ba9305f27d6325a7a28d/docs/README.md>`_


######################################################
Defining constraints
######################################################

The DL-FIND library supports bond, angle, dihedral constraints as well as frozen atoms.

Similar to the geomeTRIC interface, a dictionary defining constraints (*constraints* keyword) should be provided.
Syntax to use for the constraints dictionary:

.. code-block:: python

    constraints={'bond':[[0,1]]} #This defines a bond/distance constraint between atoms 0 and 1
    constraints={'bond':[[0,1],[3,4]]} #This defines multiple bond constraints: between atoms 0 and 1 AND also between atoms 3 and 4
    constraints={'angle':[[98,99,100]]} #This defines a angle constraint between atoms 98,99 and 100
    constraints={'dihedral':[[98,99,100,101]]} #This defines a dihedral constraint between atoms 98,99,100 and 101.
    constraints={'bond':[[0,1],[3,4]], 'angle':[[98,99,100]]} #This defines 2 bond constraints and 1 angle constraint.
    constraints={'xyz':[5,6]} #This defines XYZ constraints for the indicated atoms (i.e. freeze atoms). Alternative to frozenatoms option.

Note that, unlike the geomeTRIC library, DL-FIND constraints act directly on the initial geometry provided. 
It is not possible to specify what the constraint should be (i.e. like the constrainvalue option in geomeTRIC).

If you want a specific value of a constraint to be set (e.g. a specific dihedral angle), and the initial geometry does 
not have that value of the constraint, you have to first change the geometry so that the geometry has that value, prior to starting the minimization.
This can be done in a few ways:

- Manually (e.g. in a GUI molecular editor).
- Use the geomeTRIC Optimizer with the *constrainvalue* option.
- Do a RestraintTheory minimization first. See more information on :doc:`Geometry-optimization` page.

To freeze atoms in space it is easiest to provide a list of frozen atoms by the *frozenatoms* keyword.

Alternatively, if a large number of atoms should be frozen it may be easier to provide a list of active atoms by the *actatoms* keyword.

######################################################
Examples
######################################################


**Example: Default geometry optimization in HDLC internal coordinates**

Here we use jobtype="opt" option.

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
  theory = xTBTheory()

  DLFIND_optimizer(theory=theory, fragment=frag,  jobtype="opt", maxcycle=300)

**Example: Geometry optimization in HDLC internal coordinates using icoord/iopt syntax and with constraints**

Here we manually select the DL-FIND options by specifying icoord=1 (HDLC coordinates) and iopt=3 (L-BFGS algorithm).
Also showing how constraints can be provided by providing a dictionary (same format as used in geomeTRIC interface).

.. code-block:: python

  from ash import *

  frag = Fragment(databasefile="h2o.xyz", charge=0, mult=1)
  theory = xTBTheory()

  DLFIND_optimizer(theory=theory, fragment=frag,  icoord=1, iopt=3, maxcycle=300, constraints={'bond':[[0,1]]})


**Example: Dimer saddlepoint optimization**

.. code-block:: python

  from ash import *

  frag = Fragment(xyzfile="knarr_saddle.xyz", charge=0, mult=1)
  frag2 = Fragment(xyzfile="system102-sp.xyz", charge=0, mult=1)

  theory = xTBTheory(xtbmethod="GFN2")

  DLFIND_optimizer(theory=theory, fragment=frag, fragment2=frag2, icoord=210, iopt=3, maxcycle=300)

**Example: P-RFO saddlepoint optimization**

A P-RFO saddlepoint job requires an input-Hessian.

*inithessian* controls what Hessian DL-FIND will use:

-  0: external program. if failure we go to 2-point FD
- 1: 1-point FD
- 2: 2-point FD
- 3: diagonal 1-point FD
- 4: identity matrix

For inithessian=0, ASH computes the Hessian in one of various ways.
The *hessian_choice* keyword can be set to "numfreq", "anfreq", "xtb", "file:Hessianfilename" (read from file) or defined as a 2d numpy-array.
For the "numfreq" option we can control the approximate Hessian calculated via *numfreq_npoint*, *numfreq_displacement*,
*numfreq_hessatoms* keywords.

.. code-block:: python

  from ash import *

  frag = Fragment(xyzfile="saddle_guess.xyz", charge=0, mult=1)

  theory = xTBTheory(xtbmethod="GFN2")

  # Read previously calculated Hessian from file
  hessian = np.loadtxt("Hessian")
  # Start P-RFO job with this input Hessian
  DLFIND_optimizer(theory=theory, fragment=frag, jobtype="tsopt", inithessian=0, 
            hessian_choice=hessian, maxcycle=300)

**Example: NEB**

jobtype="neb" selects a climbing-image NEB job with frozen endpoints (same as icoord=120).

.. code-block:: python

  from ash import *

  frag = Fragment(xyzfile="system10-react.xyz", charge=0, mult=1)
  frag2 = Fragment(xyzfile="system10-prod.xyz", charge=0, mult=1)

  theory = xTBTheory(xtbmethod="GFN2")

  DLFIND_optimizer(theory=theory, fragment=frag, fragment2=frag2, jobtype="neb", maxcycle=300, nimage=30)



**Example: QM/MM geometry optimization of an active-region around a metalloprotein active-site**

Here we provide a list of active atoms to the DL-FIND optimizer (everything else will be frozen).

.. code-block:: python


  from ash import *

  #Define number of cores variable
  numcores=1

  #Fe(SCH2)4 indices (inspect system_aftersolvent.pdb file to get indices)
  qmatoms=[93,94,95,96,133,134,135,136,564,565,566,567,604,605,606,607,755]

  #Defining fragment containing coordinates (can be read from XYZ-file, ASH fragment, PDB-file)
  lastpdbfile="final_MDfrag_laststep_imaged.pdb"
  fragment=Fragment(pdbfile=lastpdbfile)

  #Creating new OpenMM object from OpenMM XML files (built-in CHARMM36 and a user-defined one)
  omm = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "./specialresidue.xml"], pdbfile=lastpdbfile, periodic=True,
              platform='CPU', numcores=numcores, autoconstraints=None, rigidwater=False)

  #QM theory
  xtbobject = xTBTheory(xtbmethod="GFN1", numcores=numcores)
  #QM/MM theory
  qmmm = QMMMTheory(qm_theory=xtbobject, mm_theory=omm, fragment=fragment,
          embedding="Elstat", qmatoms=qmatoms, printlevel=1, qm_charge=-1, qm_mult=6)

  # QM/MM geometry optimization
  actatoms=read_intlist_from_file("active_atoms")

  DLFIND_optimizer(jobtype="opt", theory=qmmm, fragment=fragment, actatoms=actatoms, maxcycle=200)



######################################################
The DLFIND_optimizerClass
######################################################

The **DLFIND_optimizer** described above is actually a wrapper function around a class: **DLFIND_optimizerClass**.

It is generally recommended to use the function call above directly, most of the time.
However, if you do require more flexibility for your ASH script then it is also possible to 
create an object from the class directly and use the built-in *run* method.


.. code-block:: python

    class DLFIND_optimizerClass:
        def __init__(self,jobtype=None, fragment=None, fragment2=None, theory=None, charge=None, mult=None, 
                    maxcycle=250, tolerance=4.5E-4, tolerance_e=1E-6, 
                    printlevel=2, result_write_to_disk=True, actatoms=None, frozenatoms=None, residues=None, constraints=None,
                    icoord=None, iopt=None, nimage=None, delta=0.01, 
                    hessian_choice='numfreq', inithessian=None, 
                    numfreq_npoint=1,numfreq_displacement=0.005,numfreq_force_projection=None,
                    numfreq_hessatoms=None, print_atoms_list=None,
                    force_noPBC=False, PBC_format_option='CIF'):


Example on how to use:

.. code-block:: python

    #Create optimizer object
    optimizer = DLFIND_optimizerClass(theory=theory, fragment=fragment)
    #Run the optimizer object
    result = optimizer.run()
