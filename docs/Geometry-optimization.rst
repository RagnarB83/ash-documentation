Geometry optimization
======================================

Geometry optimizations in ASH are almost exclusively performed via an interface to the powerful geomeTRIC optimizer library  (https://github.com/leeping/geomeTRIC).

The interface to the geomeTRIC optimization library allows efficient optimization in multiple coordinate systems: TRIC, HDLC, DLC, Cartesian, redundant internals. Constraints and frozen atoms are supported natively.
Any ASH theory object can be used using in principle any available Hamiltonian (implemented analytical gradient strongly recommended though).
Furthermore, the "ActiveRegion" feature allows definition of an active region that allows efficient QM/MM optimizations of large systems (where most atoms are frozen). 
Only the active region coordinates are in this case passed to geomeTRIC.
ASH features a full-featured interface to geomeTRIC that allows flexible constraint input, QM/MM optimizations (via ActiveRegion feature), 
relaxed and unrelaxed 1D/2D surface scans (see  :doc:`surfacescan`) and more.

If you use geometry optimizations in ASH using the geomeTRIC library make sure to cite the article:

*Geometry optimization made simple with translation and rotation coordinates*  by    Lee-Ping Wang, Chenchen Song, *J. Chem. Phys.* **2016**, *144*, 214108. 

######################################################
geomeTRICOptimizer
######################################################

The geomeTRICOptimizer function can also be called via the shorter aliases: 
**Optimizer** or **Opt**.

.. code-block:: python

    def geomeTRICOptimizer(theory=None, fragment=None, coordsystem='tric', frozenatoms=None, constraints=None, constraintsinputfile=None, 
                        constrainvalue=False, maxiter=50, ActiveRegion=False, actatoms=None, convergence_setting=None, conv_criteria=None,
                        print_atoms_list=None, charge=None, mult=None):


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


Finally an alternative way of specifying constraint is to provide a file with the constraints defined according to the syntax of the geomeTRIC library.
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

*Example: Setting convergence criteria to GAU_TIGHT:*

.. code-block:: python

    from ash import *

    frag=Fragment(xyzfile="h2o.xyz",charge=0, mult=1)
    xtbcalc=xTBTheory(xtbmethod='GFN1')

    geomeTRICOptimizer(theory=xtbcalc, fragment=frag, convergence_setting='GAU_TIGHT')

*Example: Setting convergence criteria manually:*

.. code-block:: python

    from ash import *

    frag=Fragment(xyzfile="h2o.xyz",charge=0, mult=1)
    xtbcalc=xTBTheory(xtbmethod='GFN1')

    conv_criteria_dict = {'convergence_energy' : 1e-6, 'convergence_grms' : 1e-5, 'convergence_gmax' : 1.5e-5, 
        'convergence_drms' : 4.0e-5, 'convergence_dmax' : 6.0e-5 }
    geomeTRICOptimizer(theory=xtbcalc, fragment=frag, conv_criteria=conv_criteria_dict)


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
    finalenergy = optimizer.run()