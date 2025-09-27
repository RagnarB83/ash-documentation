Geometry optimization
======================================

Geometry optimizations in ASH are almost exclusively performed via an interface to the powerful `geomeTRIC optimizer library <https://github.com/leeping/geomeTRIC>`_.
In addition there is a very simple optimizer, *SimpleOpt*, that performs geometry optimizations exclusively using Cartesian coordinates,
by basic algorithms such as steepest descent and LBFGS.
There is also a preliminary interface to the powerful `DL-FIND library <https://www.itheoc.uni-stuttgart.de/research/kaestner/research/dlfind/>`_ , which will be improved in the near future.


######################################################
geomeTRICOptimizer
######################################################

The interface to the geomeTRIC optimization library allows efficient optimization in multiple coordinate systems: TRIC, HDLC, DLC, Cartesian, redundant internals. Constraints and frozen atoms are supported natively.
Any ASH theory object can be used using in principle any available Hamiltonian (implemented analytical gradient strongly recommended though).
Furthermore, the "ActiveRegion" feature allows definition of an active region that allows efficient QM/MM optimizations of large systems (where most atoms are typically frozen). 
Only the active region coordinates are in this case passed to geomeTRIC.
ASH features a full-featured interface to geomeTRIC that allows flexible constraint input, QM/MM optimizations (via ActiveRegion feature), 
relaxed and unrelaxed 1D/2D surface scans (see  :doc:`surfacescan`), saddlepoint optimizations and more.

If you use geometry optimizations in ASH using the geomeTRIC library, make sure to cite the article:

`Geometry optimization made simple with translation and rotation coordinates <https://doi.org/10.1063/1.4952956>`_
by Lee-Ping Wang, Chenchen Song, *J. Chem. Phys.* **2016**, *144*, 214108. 


The geomeTRICOptimizer function can also be called via the shorter aliases: 
**Optimizer** or **Opt**.

.. code-block:: python

    def geomeTRICOptimizer(theory=None, fragment=None, coordsystem='tric', frozenatoms=None, constraints=None, constraintsinputfile=None, 
                        constrainvalue=False, maxiter=50, ActiveRegion=False, actatoms=None, convergence_setting=None, conv_criteria=None,
                        TSOpt=False, hessian=None, print_atoms_list=None, charge=None, mult=None):


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
   * - ``coordsystem``
     - string
     - 'tric'
     - | Which coordinate system to use during optimization. Options: 'tric', 'hdlc', 'dlc', 'prim', 'cart'  
       | Default: 'tric' (TRIC: translation+rotation internal coordinates), for an active region 'hdlc' is used instead.
   * - ``maxiter``
     - integer
     - 100
     - Maximum number of optimization iterations before giving up.
   * - ``TSOpt``
     - Boolean
     - False
     - Whether to do saddlepoint/TS optimization or not. 
   * - ``hessian``
     - string
     - None
     - | Hessian option for geomeTRIC. Options: 'never', 'first', 'last', 'first+last', 'each', or 'file:<path>''
       | Keywords refer to when the exact Hessian is calculated or the path to an external Hessian-file.
       | Default: No Hessian ('never') for TSOpt=False; 
       | for TSOpt=True the Hessian is calculated in the first step ('first').
   * - ``partial_hessian_atoms``
     - list
     - None
     - | List of atom indices for which a partial numerical Hessian will be calculated.
   * - ``frozenatoms``
     - list
     - None
     - List of frozen atoms during the optimization. The specified atoms will be added as Cartesian constraints.
   * - ``constraints``
     - Dict
     - None
     - Dictionary that specifies bond, angle or dihedral constraints. See Constraints section for syntax.
   * - ``constrainvalue``
     - Boolean
     - False
     - | Whether the value of the constraint is specified in the constraint 
       | specification or not (constraints dict above).
   * - ``constraintsinputfile``
     - string
     - None
     - | Alternative to the constraints dictionary. The file name of a constraints input file. 
       | See Constraints section for syntax. 
   * - ``convergence_setting``
     - string
     - None.
     - | Specifies the type of convergence criteria. Options: 'ORCA', 'Chemshell', 'ORCA_TIGHT', 'GAU',
       | 'GAU_TIGHT', 'GAU_VERYTIGHT', 'SuperLoose'. See Convergence section for details.
   * - ``conv_criteria``
     - Dict
     - None
     - Alternative manual way of defining convergence criteria. See Convergence section for details.
   * - ``print_atoms_list``
     - list
     - None
     - | Optional list of atom indices for which the Cartesian coordinates will be printed out
       | in each optimization step. Defaults to all atoms in general, active-atoms if ActiveRegion
       | is True or QM-region for QM/MM theories.
   * - ``ActiveRegion``
     - Boolean
     - False
     - | Whether to use an Active Region during the optimization. This requires setting
       |  the number of active atoms (actatoms list) below.
   * - ``actatoms``
     - list
     - None
     - List of atom indices that are active during the optimization job. All other atoms are frozen. 
   * - ``printlevel``
     - integer
     - 0
     - The printlevel to use during optimization.
   * - ``charge``
     - integer
     - None
     - | Optional specification of the charge of the system (if QM)
       | if the information is not present in the fragment.
   * - ``mult``
     - integer
     - None
     - | Optional specification of the spin multiplicity of the system (if QM) 
       | if the information is not present in the fragment.



######################################################
Examples
######################################################

*Geometry optimization of an H2O fragment at the xTB level, with default options.*


.. code-block:: python

    from ash import *

    frag=Fragment(databasefile="h2o.xyz",charge=0, mult=1)
    xtbcalc=xTBTheory(xtbmethod='GFN1')

    Optimizer(theory=xtbcalc, fragment=frag)

*Geometry optimization of an H2O fragment at the BP86 DFT-level with ORCA, with default options.*

.. code-block:: python

    from ash import *

    frag=Fragment(databasefile="h2o.xyz",charge=0, mult=1)
    orcacalc=ORCATheory(orcasimpleinput='! BP86 def2-SVP def2/J tightscf')

    Optimizer(theory=orcacalc, fragment=frag)


*Geometry optimization of a QM/MM system with an active region:*

.. code-block:: python

    from ash import *

    #Fe(SCH2)4 indices (inspect system_aftersolvent.pdb file to get indices)
    qmatoms=[93,94,95,96,133,134,135,136,564,565,566,567,604,605,606,607,755]

    #Defining fragment containing coordinates (can be read from XYZ-file, ASH fragment, PDB-file)
    lastpdbfile="final_MDfrag_laststep_imaged.pdb"
    fragment=Fragment(pdbfile=lastpdbfile)
    #Creating new OpenMM object from OpenMM XML files (built-in CHARMM36 and a user-defined one)
    omm = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "./specialresidue.xml"], pdbfile=lastpdbfile, periodic=True,
                platform='CPU', numcores=numcores, autoconstraints=None, rigidwater=False)
    #QM theory: r2SCAN-3c DFT-composite method using ORCA
    orca = ORCATheory(orcasimpleinput="! r2SCAN-3c tightscf", numcores=numcores)
    #QM/MM theory
    qmmm = QMMMTheory(qm_theory=orca, mm_theory=omm, fragment=fragment,
            embedding="Elstat", qmatoms=qmatoms, printlevel=1)

    # QM/MM geometry optimization with an active region (here QM-region only)
    Optimizer(fragment=fragment, theory=qmmm, ActiveRegion=True, actatoms=qmatoms, maxiter=200, charge=-1, mult=6)

######################################################
Files created
######################################################

Once the Optimizer is done, the coordinates in the Fragment object are automatically updated (to be the optimized coordinates) so the Fragment could be immediately used for another job (e.g. a NumFreq job).

During the geometry optimization the following files are created and updated:

.. code-block:: text

  - initialxyzfiletric.xyz: The initial XYZ coordinates read into the geomeTRIC optimizer
  - geometric_OPTtraj.log : A logfile containing the optimizer settings and also the data for each Step (RMS/Max Gradient and Displacement values and Energy)
  - geometric_OPTtraj_optim.xyz: An XYZ trajectory containing the geometry of each optimization step. Can be visualized using e.g. VMD/Chemcraft.
  - Fragment-currentgeo.xyz: An XYZ-file containing the coordinates of the current optimization step.

If the geometry optimization converges without problems, the 'Fragment-optimized.xyz' file is available, which is an XYZ-file containing the optimized coordinates.

If the theory level is a QMMMTheory object then additional files are created for convenience:

.. code-block:: text

  - geometric_OPTtraj_Full.xyz : An XYZ trajectory file containing the full system (not just the active region).
  - geometric_OPTtraj_QMregion.xyz:  An XYZ trajectory file containing the QM-region only.
  - optimization_energies.log: A logfile containing the QM-energy, MM-energy and QM/MM-energy for each optimization step.

######################################################
Constraints
######################################################

Constraints can be provided to the Optimizer in two different ways: either via providing a dictionary definition of the constraints (*constraints* keyword) or alternatively by providing a valid constraint-parameter file (*constraintsinputfile* keyword) in geomeTRIC library syntax.
The former way is recommended.
Syntax to use for the constraints dictionary:

.. code-block:: python

    constraints_dict={'bond':[[0,1]]} #This defines a bond/distance constraint between atoms 0 and 1
    constraints={'bond':[[0,1],[3,4]]} #This defines multiple bond constraints: between atoms 0 and 1 AND also between atoms 3 and 4
    constraints={'angle':[[98,99,100]]} #This defines a angle constraint between atoms 98,99 and 100
    constraints={'dihedral':[[98,99,100,101]]} #This defines a dihedral constraint between atoms 98,99,100 and 101.
    constraints={'bond':[[0,1],[3,4]], 'angle':[[98,99,100]]} #This defines 2 bond constraints and 1 angle constraint.

*Example:*

.. code-block:: python

    from ash import *

    h2ostring="""
    O        1.586357512      0.000000000     -6.179217668
    H        1.586357512      0.759337000     -5.583174668
    H        1.586357512     -0.759337000     -5.583174668
    """
    frag=Fragment(coordsstring=h2ostring,charge=0, mult=1)
    
    xtbcalc=xTBTheory(xtbmethod='GFN1')

    #Defining constraints: 1 O-H bond is constrained to its current value 
    constraints_dict ={'bond':[[0,1]]}
    Optimizer(theory=xtbcalc, fragment=frag, constraints=constraints_dict)


When the above syntax is used, the constraint is applied according to the initial geometry provided (the O-H bond (between atoms 0 and 1) is constrained to 0.965 Å)). 
If one wants to constrain e.g. a bond distance to a specific value
then this can be done by providing an extra value to the list while also providing the *constrainvalue=True* option.

.. code-block:: python

    constraints_dict={'bond':[[0,1,0.97]]} #This defines a bond/distance constraint of 0.97 Å between atoms 0 and 1
    constraints={'bond':[[0,1,0.97],[3,4,0.97]]} #This defines multiple bond constraints of 0.97 Å.
    constraints={'angle':[[98,99,100,104.5]]} #This defines an angle constraint of 104.5° between atoms 98,99 and 100


It should be noted that geomeTRIC handles constraints a little different than many codes and the constraints are not fully enforced until the end of the constrained optimization.
There are also rare cases where the minimization stalls because constraints can not be fully satisfied.
Changing the convergence tolerance of the constraints are then necessary (see Convergence criteria section below).

*Example:*

.. code-block:: python

    from ash import *

    h2ostring="""
    O        1.586357512      0.000000000     -6.179217668
    H        1.586357512      0.759337000     -5.583174668
    H        1.586357512     -0.759337000     -5.583174668
    """
    frag=Fragment(coordsstring=h2ostring,charge=0, mult=1)
    
    xtbcalc=xTBTheory(xtbmethod='GFN1')

    #Defining constraints: 1 O-H bond is constrained to to 0.97 Å
    constraints_dict ={'bond':[[0,1,0.97]]}
    Optimizer(theory=xtbcalc, fragment=frag, constraints=constraints_dict, constrainvalue=True)


An alternative way of specifying constraint is to provide a file with the constraints defined according to the syntax of the geomeTRIC library.
See `geomeTRIC constraints file format <https://github.com/leeping/geomeTRIC/blob/master/examples/constraints.txt>`_ for more information.
The drawback of this approach is that atom indices will use 1-based indexing (unlike ASH in general), indices would have to be checked and modified in case of an Active Region,
and finally either a global path to this file needs to be provided (so that the computing node can access it) or the file copied over to the scratch on the node.

Format of the constraint file (*Warning: geomeTRIC counts from 1 (unlike ASH).*)

.. code-block:: text

    $freeze
    bond 5 6
    xyz 5 xyz
    xy 5-11,13,35
    $set
    angle 3 1 2 30.0
    z 36 10.0
    $scan
    dihedral 4 2 3 5 0.0 180.0 19


Finally it should be noted that the default constraint algorithm in geomeTRIC enforces constraints in a special way `as documented <https://geometric.readthedocs.io/en/latest/constraints.html#enforcing-constraint-satisfaction`_.
The exact constraints are not fully imposed until late in the optimization.
This behaviour can be controlled by enforcing a particular threshold when constraints are fully enforced.
In ASH this is controlled by the keyword *enforce_constraints*=X  where X is the desired threshold value for the constraint (Bohr for distance, radian for angles/dihedrals).


######################################################
Rigid optimization
######################################################

An alternative to user-supplied constraints is a rigid optimization as described in the `geomeTRIC library documentation <https://geometric.readthedocs.io/en/latest/constraints.html#rigid-optimizations>`_
In a rigid optimization, only the intermolecular positions and orientations are optimized while internal bonds etc. are kept fixed.
This option is available in ASH via the rigid Boolean keyword. A revised constraint algorithm is automatically used.

Example:

.. code-block:: python

  from ash import *

  coordsstring="""
  C       -0.911459798      0.000000000     -1.679453249
  O       -0.911459798      0.000000000     -2.895258249
  C       -0.911459798      1.292616000     -0.880357249
  C       -0.911459798     -1.292616000     -0.880357249
  H       -0.911459798      2.146342000     -1.559037249
  H       -0.911459798     -2.146342000     -1.559037249
  H       -0.031457798      1.340052000     -0.228371249
  H       -1.791461798      1.340052000     -0.228371249
  H       -1.791461798     -1.340052000     -0.228371249
  H       -0.031457798     -1.340052000     -0.228371249
  C        3.397788926     -2.560955125     -0.179308061
  H        3.397788926     -3.396837760     -0.883608867
  H        2.505939926     -1.948477943     -0.381012003
  H        4.289637926     -1.948477943     -0.381012003
  O        3.397788926     -3.120042498      1.123961354
  H        3.397788926     -2.390939481      1.756408337
  """
  frag = Fragment(coordsstring=coordsstring, charge=0, mult=1)
  theory = xTBTheory(xtbmethod="GFN1")

  # rigid opt
  Optimizer(theory=theory,fragment=frag, rigid=True)


######################################################
Convergence criteria
######################################################

The default convergence criteria of **geomeTRICOptimizer** are the same as used by the ORCA program by default. It is possible to change these default criteria by either specifying a string (*convergence_setting* keyword)
or manually setting all the criteria by providing a dictionary (*conv_criteria* keyword)

convergence_setting options (default: 'ORCA'). What type of convergence criteria to use. 

Valid options are: 'ORCA', 'ORCA_TIGHT', 'Chemshell', 'GAU', 'GAU_TIGHT', 'GAU_VERYTIGHT', 'SuperLoose'.



.. list-table::
   :widths: 15 15 15 15 15 15
   :header-rows: 1

   * - String keyword
     - convergence_energy
     - convergence_grms value
     - convergence_gmax
     - convergence_drms
     - convergence_dmax
   * - ``ORCA``
     - 5.0e-6
     - 1.0e-4
     - 3.0e-4
     - 2.0e-3
     - 4.0e-3
   * - ``ORCA_TIGHT``
     - 1.0e-6
     - 3.0e-5
     - 1.0e-4
     - 6.0e-4
     - 1.0e-3
   * - ``Chemshell``
     - 1.0e-6
     - 3.0e-4
     - 4.5e-4
     - 1.2e-3
     - 1.8e-3
   * - ``GAU``
     - 1.0e-6
     - 3.0e-4
     - 4.5e-4
     - 1.2e-3
     - 1.8e-3
   * - ``GAU_TIGHT``
     - 1.0e-6
     - 1.0e-5
     - 1.5e-5
     - 4.0e-5
     - 6.0e-5
   * - ``GAU_VERYTIGHT``
     - 1.0e-6
     - 1.0e-6
     - 2.0e-6
     - 4.0e-6
     - 6.0e-6
   * - ``SuperLoose``
     - 1.0e-1
     - 1.0e-1
     - 1.0e-1
     - 1.0e-1
     - 1.0e-1


.. note:: Additionally all settings above also use a tolerance for enforcing any present constraints (convergence_cmax) which is by default convergence_cmax=1.0e-2

*Example: Setting convergence criteria to GAU_TIGHT:*

.. code-block:: python

    from ash import *

    frag=Fragment(xyzfile="h2o.xyz",charge=0, mult=1)
    xtbcalc=xTBTheory(xtbmethod='GFN1')

    Optimizer(theory=xtbcalc, fragment=frag, convergence_setting='GAU_TIGHT')

*Example: Setting convergence criteria manually:*

.. code-block:: python

    from ash import *

    frag=Fragment(xyzfile="h2o.xyz",charge=0, mult=1)
    xtbcalc=xTBTheory(xtbmethod='GFN1')

    conv_criteria_dict = {'convergence_energy' : 1e-6, 'convergence_grms' : 1e-5, 'convergence_gmax' : 1.5e-5, 
        'convergence_drms' : 4.0e-5, 'convergence_dmax' : 6.0e-5, 'convergence_cmax' : 1.0e-1 }
    Optimizer(theory=xtbcalc, fragment=frag, conv_criteria=conv_criteria_dict)


geomeTRIC handles constraints a little different than many codes and the constraints are not fully enforced until the end of the constrained optimization.
There are also cases where the minimization stalls because constraints can not be fully satisfied.
Changing *convergence_cmax* to a smaller value than 1.0e-2 may be necessary in these cases.

######################################################
Transition-State/Saddlepoint Optimization
######################################################

A direct transition-state/Saddle-points optimization can be performed in the Optimizer via an eigenvector-following
algorithm as implemented in the geometric library. This option is actived by the *TSOpt=True* keyword as shown below:

.. code-block:: python

  from ash import *

  frag=Fragment(xyzfile="hf.xyz", charge=0, mult=1) #Fragment object creation
  ORCAcalc = ORCATheory(orcasimpleinput="! BP86 def2-SVP  tightscf") #ORCATheory object creation

  #TSOpt=True enables saddlepoint optimization in geomeTRIC. Note: Exact Hessian is calculated in the first step by default.
  Optimizer(fragment=frag, theory=ORCAcalc, coordsystem='tric', TSOpt=True)

It is important to realize that a direct TS-Optimization like this only makes sense when a good guess for the 
saddlepoint geometry is available, e.g. if the geometry has been estimated from a surface scan, aprevious NEB/NEB-CI job etc. 

Additionally, the algorithm requires a good initial approximation to the Hessian to be successful (unlike a regular minimization).
By default, if an Hessian-option (*hessian* keyword) is not specified for a *TSOpt=True* job, then an exact Hessian is estimated in the first step (*hessian='first'* option)
by a numerical-frequency calculation. The exact Hessian option can be expensive, especially if the system is large or the number of active atoms is large (often the case for QM/MM optimizations).

Options include:

- *hessian* = 'first'. Calculate the Hessian in the first step using geometric library.
- *hessian* = 'each'. Calculate the Hessian in each step (very expensive) using geometric library.
- *hessian* = 'partial'. Calculate an exact partial Hessian. Requires *partial_hessian_atoms* keyword to be defined.
- *hessian* = 'xtb'. Calculate an exact Hessian but using the cheap xTB level of theory. 
- *hessian* = '1point'. Calculate an exact Hessian using ASH using a cheap 1-point formula (requires 3*N energy+gradient displacement calculations).
- *hessian* = '2point'. Calculate an exact Hessian using ASH using a 2-point formula (requires 2*3*N energy+gradient displacement calculations).
- *hessian* = 'file\:/path/to/Hessianfile'. Read Hessian from file.
- *hessian* = <Numpy array>. Read Hessian from Numpy array.

The option to read in the Hessian from a file or a Numpy array offers a lot of flexibility.
Any Hessian (however calculated) can be read in from a file (or Numpy array) as long as it has the correct dimensions of the system (3*N, where N=numatoms).
See :doc:`module_freq` for options how an Hessian can be calculated numerically or analytically using various ASH Theories.
A Hessian-file is always written to disk (as text) following a successful NumFreq/AnFreq calculation.
The Hessian is also part of the Results object that is returned by NumFreq/AnFreq.


######################################################
Intrinsic Reaction Coordinate (IRC)
######################################################

The Intrinsic Reaction Coordinate (IRC) method is intended to find the minimum energy pathway, 
starting from a previously optimzied saddlepoint geometry and a Hessian. 
The primary intention is usually to validate that the located saddlepoint are connected to the assumed reaction and product minima.

The IRC method is implemented in the `geomeTRIC library <https://geometric.readthedocs.io/en/latest/irc.html>`_
and the ASH interface supports it. Note that version 1.1 of the geometric library is required for IRC.

IRC Example:

.. code-block:: python

  from ash import *

  coordsstring="""
  C   -0.1088783634   -0.6365101639    0.0043221742
  N   -0.6393457902    0.4205365638    0.0052498438
  H    0.7532976101    0.2173493463   -0.0090384631
  """
  frag = Fragment(coordsstring=coordsstring, charge=0, mult=1)
  theory = xTBTheory(xtbmethod="GFN1")
  # NumFreq to get the Hessian
  result = NumFreq(theory=theory, fragment=frag)
  # IRC with an input Hessian
  Optimizer(theory=theory,fragment=frag, maxiter=200, irc=True, hessian=result.hessian)



######################################################
The geomeTRICOptimizer class
######################################################

The **geomeTRICOptimizer** described above is actually a wrapper function around a class: **GeomeTRICOptimizerClass**.

It is strongly recommended to use the function described above, however, if you do require more flexibility for your 
ASH script then it is also possible to create an object from the class directly and use the built-in *run* method.


.. code-block:: python

    class GeomeTRICOptimizerClass:
            def __init__(self,theory=None, fragment=None, charge=None, mult=None, coordsystem='tric', frozenatoms=None, 
                        constraintsinputfile=None, constraints=None, constrainvalue=False, maxiter=50, print_atoms_list=None,
                        ActiveRegion=False, actatoms=None, convergence_setting=None, conv_criteria=None):


Example on how to use:

.. code-block:: python

    #Create optimizer object
    optimizer = GeomeTRICOptimizerClass(theory=theory, fragment=fragment, charge=0, mult=1))
    #Run the optimizer object
    result = optimizer.run()

######################################################
DL-FIND interface
######################################################

In addition to geomeTRIC, ASH features an interface to the powerful `DL-FIND library <https://www.itheoc.uni-stuttgart.de/research/kaestner/research/dlfind/>`_ .
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


.. code-block:: python


  def DLFIND_optimizer(jobtype=None, theory=None, fragment=None, fragment2=None, charge=None, mult=None, 
                      maxcycle=250, tolerance=4.5E-4, tolerance_e=1E-6,
                      actatoms=None, frozenatoms=None, residues=None, constraints=None,
                      printlevel=2, NumGrad=False, delta=0.01,
                      icoord=None, iopt=None, nimage=None, 
                      hessian_choice="numfreq", inithessian=0, 
                      numfreq_npoint=1, numfreq_displacement=0.005, numfreq_hessatoms=None,
                      numfreq_force_projection=None, print_atoms_list=None):


Since the interface supports active regions, **DLFIND_optimizer** can e.g. be used for QM/MM geometry optimizations in ASH.

Some of the specific DL-FIND options must be chosen via the icoord and iopt keywords.
See `DL-FIND manual <https://github.com/digital-chemistry-laboratory/libdlfind/blob/4167998d16d8dac4a484ba9305f27d6325a7a28d/docs/documentation.pdf>`_
See also `libdlfind README <https://github.com/digital-chemistry-laboratory/libdlfind/blob/4167998d16d8dac4a484ba9305f27d6325a7a28d/docs/README.md>`_



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
SimpleOpt
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
  #Optimization
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
