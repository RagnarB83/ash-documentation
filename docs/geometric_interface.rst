geomeTRIC Optimizer
======================================

The interface to the geomeTRIC optimization library allows efficient optimization in multiple coordinate systems: TRIC, HDLC, DLC, and redundant internals. 
Constraints and frozen atoms are supported natively.
Any ASH theory object can be used, in principle with any available Hamiltonian (implemented analytical gradient strongly recommended though).
Furthermore, the "ActiveRegion" feature allows definition of an active region that allows efficient QM/MM optimizations of large systems (where most atoms are typically frozen). 
Only the active region coordinates are in this case passed to geomeTRIC.
ASH features a full-featured interface to geomeTRIC that allows flexible constraint input, QM/MM optimizations (via ActiveRegion feature), 
relaxed and unrelaxed 1D/2D surface scans (see  :doc:`surfacescan`), saddlepoint optimizations and more.

If you use geometry optimizations in ASH using the geomeTRIC library, make sure to cite the article:

`Geometry optimization made simple with translation and rotation coordinates <https://doi.org/10.1063/1.4952956>`_
by Lee-Ping Wang, Chenchen Song, *J. Chem. Phys.* **2016**, *144*, 214108. 


################################
Installation
################################

To install the geometric library, it is easiest to use pip:

.. code-block:: shell

    pip install geometric

For problems associated with the geomeTRIC library see: `geomeTRIC documentation <https://geometric.readthedocs.io/en/latest/>`_ and `geomeTRIC Github repository <https://github.com/leeping/geomeTRIC>`_

################################
geomeTRICOptimizer function
################################

The *geomeTRICOptimizer* function can also be called via the shorter aliases: 
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
    constraints={'xyz':[5,6]} #This defines XYZ constraints for the indicated atoms (i.e. freeze atoms). Alternative to frozenatoms option.
    constraints={'x':[5,6]} #This defines a partial X Cartesian constraint for the indicated atoms.
    constraints={'xy':[5,6]} #This defines a partial XY Cartesian constraint for the indicated atoms.

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


Finally it should be noted that the default constraint algorithm in geomeTRIC enforces constraints in a special way `as documented <https://geometric.readthedocs.io/en/latest/constraints.html#enforcing-constraint-satisfaction>`_.
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
algorithm as implemented in the geometric library. 
Specific geomeTRIC uses a restricted-step partitioned rational function optimization (RS-P-RFO) method,  `geomeTRIC documentation <https://geometric.readthedocs.io/en/latest/transition.html>`_ 
for more information on the algorithm and its implementation in geomeTRIC.

This option is actived by the *TSOpt=True* keyword as shown below:

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
Periodic Boundary Conditions with geomeTRIC
######################################################

The geomeTRIC library is normally intended only for molecular systems (i.e. no translational symmetry), 
not systems with periodic boundary conditions.
In previous ASH versions, if the system was described by a Theory-interface supporting PBCs (e.g. CP2K)
but optimized with geomeTRICOptimizer, a frozen-lattice optimization was automatically performed (i.e. only atoms of the cell were optimized).

More recently, the ASH interface to geomeTRIC, supports a way of coaxing the geomeTRIC library to 
simultaneously minimize atom positions and cell vectors of a periodic system.
This is performed by converting the cell-vectors into dummy atoms that are simultaneously optimized in geomeTRIC's internal coordinates.
This option will currently only work for HDLC internal coordinates (automatically enforced).

To utilize this option, one only needs to provide a Theory object that has native support for PBCs and has enabled PBCs in object (usually via *periodic* = True keyword)
The geomeTRICOptimizer will then automatically perform an atom+cellvector optimization. 

If this behaviour is not desired one can turn it off (*force_noPBC* = False in geomeTRICOptimizer) which should then correspond to a frozen lattice calculation.

During the optimization ASH will write XYZ-files and XYZ trajectories as normal.
Once the optimization finishes will additionally write the coordinates in a file-format suited for PBCs.
This format can be chosen by keyword *PBC_format_option* which is by default set to 'CIF' (see `CIF file format <https://www.ccdc.cam.ac.uk/community/access-deposit-structures/deposit-a-structure/guide-to-cifs/>`_ but other other options are 'XSF' 
(see `XSF file format <http://www.xcrysden.org/doc/XSF.html>`_) or 'POSCAR' (see `POSCAR file format <https://www.vasp.at/wiki/POSCAR>`_)


See also: :doc:`Periodic-systems` .

Example below shows a way of using geomeTRIC to optimize an ammonia crystal using CP2K as a PBC-DFT theory.

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
    Optimizer(theory=theory, fragment=frag, coordsystem="hdlc", PBC_format_option="XSF")



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
