QM/MM Theory
==========================

QM/MM in ASH is highly flexible as one can combine any QM-theory in ASH (that supports pointcharge embedding) with an MMTheory object of class NonbondedTheory (see :doc:`MM-interfaces`) or OpenMMTheory (see :doc:`OpenMM-interface`).

To do QM/MM, one combines a defined QMtheory object (:doc:`QM-interfaces`) and an MMtheory object in QMMMTheory object
and then specifies which atoms are QM and which are MM and the type of QM-MM coupling (typically electrostating embedding).
Note that in contrast to a QMtheory object or an MMtheory object, we also pass a fragment object to a QMMMTheory object so that
the QMMMTheory object can define the division of QM-region and MM-regions.


######################################
QMMMTheory class
######################################

.. code-block:: python
 
  class QMMMTheory:
      def __init__(self, qm_theory=None, qmatoms=None, fragment=None, mm_theory=None, charges=None,
                  embedding="Elstat", printlevel=2, numcores=1, actatoms=None, frozenatoms=None, excludeboundaryatomlist=None,
                  unusualboundary=False, openmm_externalforce=False, TruncatedPC=False, TruncPCRadius=55, TruncatedPC_recalc_iter=50,
                  qm_charge=None, qm_mult=None):

**QMMMTheory** options:

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``qm_theory``
     - ASHTheory
     - None
     - | Required: The theory level for the QM-region. Must be a valid ASH Theory
       | that supports electrostatic embedding.
   * - ``qmatoms``
     - list
     - None
     - | Required: List of QM-atom indices that defined the QM-region. Atoms not
       | in list are treated as MM atoms.
   * - ``mm_theory``
     - ASHTheory
     - None
     - | Required: The theory level for the MM-region. Must be an object of class
       | OpenMMTheory or NonbondedTheory.
   * - ``fragment``
     - ASH Fragment
     - None
     - Required: ASH fragment, needed for setting up QM-region and MM-region.
   * - ``qm_charge``
     - integer
     - None
     - | Optional: Specify the charge of the QM-region. This takes precedence
       | over other charge specifications.
   * - ``qm_mult``
     - integer
     - None
     - | Optional: Specify the spin multiplicity of the QM-region. This takes
       | precedence over other mult specifications.
   * - ``charges``
     - list
     - None
     - | Optional: list of atom charges. If not defined then charges will be read
       | from mm_theory.
   * - ``printlevel``
     - integer
     - 2
     - | Optional: The printlevel setting. If printlevel >= 3 then more printing
       | and gradient files are written to disk.
   * - ``numcores``
     - integer
     - 1
     - | Optional: Number of CPU cores to use for qm_theory. If defined, takes
       | precedence over QMTheory setting.
   * - ``excludeboundaryatomlist``
     - list
     - None
     - Optional: List of atoms that are excluded from adding linkatoms to.
   * - ``unusualboundary``
     - Boolean
     - False
     - Optional: Boundary-option: overrides ASH from quitting if an unusual QM-MM boundary is found. 
   * - ``openmm_externalforce``
     - Boolean
     - False
     - Optional: Option for passing QM/MM force as an external force to OpenMMTheory.
   * - ``TruncatedPC``
     - Boolean
     - False
     - Optional: Truncated Pointcharge Option on or off.
   * - ``TruncPCRadius``
     - float
     - 55
     - Optional: Truncated PC option; Radius (Å) for the truncated PC region.
   * - ``TruncatedPC_recalc_iter``
     - integer
     - 50
     - Optional: Truncated PC option; frequency for recalculating with full PC field.
   * - ``actatoms``
     - list
     - None
     - | Optional: List of active atoms in QM/MM. NOTE: Only compatible if
       | mm_theory is of NonBondedTheory class.
   * - ``frozenatoms``
     - list
     - None
     - | Optional: List of frozen atoms in QM/MM, alternative to actatoms. 
       | NOTE: Only compatible if mm_theory is of NonBondedTheory class.


Dummy example showing how to combine a QMTheory and MMTheory object into a QMMMTheory object:

.. code-block:: python

    #QM theory: xTB
    qm = xTBTheory(xtbmethod='GFN1')

    #Creating new OpenMM object from OpenMM XML files (built-in CHARMM36 and a user-defined one)
    omm = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "./specialresidue.xml"], pdbfile="topology.pdb", 
              periodic=True, platform='CPU', numcores=numcores, autoconstraints=None, rigidwater=False)

    #QM/MM theory object. QM-region defined as atom indices 500,501,502 and 503
    qmmm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=fragment, embedding="Elstat", 
              qmatoms=[500,501,502,503], printlevel=2, qm_charge=-1, qm_mult=6)


##################################################################
Defining the charge and spin multiplicity of the QM-region
##################################################################


To define the charge and spin multiplicity of the QM-region in QM/MM calculations you can choose between 3 options:

\- Define qm_charge and qm_mult attributes when defining the QMMMTheory object (**recommended**):

.. code-block:: python

    qmmm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=frag, qm_charge=-1, qm_mult=6)

\- Define as input to the job-function (e.g. Singlepoint):

.. code-block:: python

    Singlepoint(theory=qmmm, fragment=frag, charge=-1, mult=6)

\- Provide the information in the fragment definition:

.. code-block:: python

    frag=Fragment(xyzfile="system.xyz", charge=-1, mult=6)

This information will be passed onto the QM-program when called. The qm_charge/qm_mult option takes precedence over the other options, followed by the job-type keyword.

Note that the specified charge and multiplicity of the QM-region needs to be consistent with what chemical groups are present in the QM-region. 

######################################
Defining QM-region and active region
######################################

The QM-region needs to be defined in the QMMMTheory object by specifying which atom indices (of the full system) should be QM-atoms (everything else is MM).

.. code-block:: python

    qmmm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=fragment, 
              qmatoms=[500,501,502,503], qm_charge=-1, qm_mult=6)

Similarly the active-region (when performing a geometry optimization) needs to be defined by specifying which atoms are allowed to move.
This information should be provided to the Optimizer.

.. code-block:: python

    Optimizer(fragment=fragment, theory=QMMMobject, ActiveRegion=True, actatoms=[400,450,500,501,502,503,550,600,700])

Both the QM-region and Active regions are thus defined as simple Python lists of integers (corresponding to atom indices).
This approach allows you considerable flexibility in defining the QM/MM job. The QM-regions and active regions can be the same or different (quite common).

Definition of the QM-region when part of a larger molecule (e.g. a protein) requires a bit of insight into the system and knowledge of how the QM/MM boundary works (see next section).
It is usually best to define the QM-region by manually creating the list of atoms. 
One can double-check whether the region is correct by using the **fragedit.py** script (see :doc:`coordinate-tools`) or check the QM-region coordinates printed in the ASH output.

The active region is typically much larger than the QM-region (for a protein, an active region of approx. 1000 atoms is common)
and it is usually inconvenient to define it manually. ASH provides a convenient function **actregiondefine** (see :doc:`coordinate-tools`) to define
such a large list of atom indices. This function can also be used to define the QM-region.

As these lists can be large it is convenient to read them from a file. 
ASH provides a function **read_intlist_from_file** (see :doc:`coordinate-tools`) to read a list of integers from a file and return a Python list.
The file should contain integers separated by spaces or newlines.

.. code-block:: python

    qmatoms = read_intlist_from_file("qmatoms")
    actatoms = read_intlist_from_file("active_atoms")

    qmmm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=fragment, 
                qmatoms=qmatoms, qm_charge=-1, qm_mult=6)
    Optimizer(fragment=fragment, theory=QMMMobject, ActiveRegion=True, actatoms=actatoms)

######################################
QM/MM boundary treatment
######################################

If the QMregion-MMregion boundary is between two bonded atoms, then a boundary correction needs to be applied.
In ASH this is treated by the popular linkatom method, combined with charge-shifting.
A hydrogen-linkatom is added to cap the QM-subsystem. The hydrogen linkatoms are only visible to the QM theory, not the MM theory.
Additionally to prevent overpolarization, the atom charge of the MMatom is shifted towards its neighbours and a dipole correction
applied by adding additional pointcharges (can be turned off by dipole_correction=False). These pointcharges are only visible to the QM theory.

The recommended way of using link atoms is to define the QM-MM boundary for two carbon atoms that are as non-polar as possible.
In the CHARMM forcefield one should additionally make sure that one does not make a QM-MM boundary through a charge-group (check topology file).
By default ASH will exit if you try to define a QM-MM covalent boundary between two atoms that are not carbon atoms (since this is almost never desired). 
To override this behaviour add "unusualboundary=True" as keyword argument when creating QMMMTheory object.

In rare cases you may want to prevent ASH from adding a linkatom for a specific QM-atom, e.g. if you are making unusual QM-MM boundaries. 
This can be accomplished like below. Note, however, that the QM-MM bonded terms will still be included.

.. code-block:: python

    #Excluding QM-atom 5785 from linkatom-creation.
   qmmmobject = QMMMTheory(qm_theory=orcaobject, mm_theory=openmmobject, fragment=frag, embedding="Elstat",
            qmatoms=qmatoms, excludeboundaryatomlist=[5785])


Special care should be taken when defining a QM-region for a biomolecular system

**General recommendations:**

- Always cut a C-C bond that is as nonpolar as possible.
- Focus on including nearby sidechains of residues that are charged (e.g. Arg, LYS, ASP, GLU) or are involved in important hydrogen bonding. 
- Amino acid sidechains are straighforward but make sure to not cut through CHARMM charge groups.
- Including protein backbone is more involved and needs careful inspection. The only good option is typically to cut the C-C bond between the C=O and the C-alpha.


See :doc:`QM-MM-boundary_tutorial` for more information on how to define a good QM/MM boundary for proteins.

############################################
How QM/MM works behind the scenes in ASH
############################################

During a QM/MM energy+gradient calculation in ASH the following steps take place:

1.	QM/MM program reads in the full model of the system, containing all atoms (no fake atoms,dummy atoms or linkatoms).
2.	ASH determines connectivity of the system, i.e. finds what atoms are bonded to each other. E.g. atom no. 17 and atom no. 18 may be close enough that ASH thinks they are bonded.
3.	The program parses the qmatoms list. The qmatoms list only contains real atoms (not linkatoms because they don‘t exist yet). The qmatoms list may e.g. contain atoms 1,2,3,4,14,15,16,17) and the program next checks if there is a covalent QM-MM boundary. Since atom no. 17 was (according to step 2)  bound to atom no. 18 (which is an MM-atom) then that means we have a covalent QM-MM boundary.
4.	ASH will next automatically calculate the need for a linkatom for all QM-MM boundaries. Any required linkatom (H) coordinates will be calculated and MM charges will be modified to account for the QM-MM boundary.
5.	The Cartesian coordinates of the QM-atoms are taken and passed to the QM-program. A hydrogen link atom is automatically added to this list of QM-coordinates (so that the QM-system will not have a dangling boundary). Additionally the MM pointcharges (also passed to the QM-program) are modified so that overpolarization will not occur (electron density at atom no. 17 and linkatom would be overpolarized by the closeness to MM-atom no. 18 ). Additional MM charges are also added so that the dipole is more realistic.
6.	The QM-program calculates the energy and gradient of the QM-atoms (with linkatoms) with the electrostatic effect of the MM-atoms included (enters the 1-electron Hamiltonian). The QM energy and gradient (of the QM-atoms and also the PC gradient of the MM atoms) is passed back to ASH.
7.	An MM calculation is performed for the whole system. The pointcharges of the atoms that have been labelled QM-atoms have been set to zero to avoid calculating the electrostatic energy (would double-count otherwise). Bonded MM terms for the same QM-atoms are removed (if present, to avoid doublecounting). Bonded terms at the QM-MM boundary are dealt with in a special way. The MM program never sees any linkatoms,charge-shifted MM charges or dipole charges. 
8.	The QM/MM energy is calculated by combining the contribution calculated by the QM-program and the MM-program. This will include the coupling energy of the QM and MM subsystems. Correction for artificial linkatom energy could be done here (not done in practice in ASH). The QM/MM gradient of the full system is assembled from the QM-gradient, MM-gradient and PC-gradient. The gradient calculated on the dummy linkatom during the QM-atoms (which does  not exist in the real system) is taken and it is projected onto the relevant MM and QM-atoms instead.
9.	The complete QM/MM gradient of the whole system is used to make a step in the relevant job.
10. if part of a geometry optimization then the step is taken so as to minimize the QM/MM gradient
11. or: if part of a dynamics simulation is taken according to Newton‘s equations  (the QM/MM gradient or force is used to calculate an acceleration which results in a change in velocity and positions of all the (real) atoms).

.. note:: Neither the geometry optimization or dynamics algorithms see or experience any linkatoms, only real atoms of the system.


######################################
QM/MM Truncated PC approximation
######################################

For large systems (e.g. > 50 000 atoms) the evaluation of the QM-pointcharge interaction (calculated by the QM-code) will start to dominate the cost of the calculation in each QM/MM calculation step.
The QM-pointcharge gradient calculation is the main culprit and it depends on the QM-code how efficiently this step is carried out for a large number of pointcharges.
ASH features a convenient workaround for this problem in QM/MM geometry optimizations. Instead of reducing the system size, ASH can temporarily reduce the size of the PC field (MM calculation size remains the same) during the geometry optimization which can speed up the calculation a lot.
The size of the truncated PC field is controlled by the TruncPCRadius variable (radius in Å) which results in a truncated spherical PC field.

The algorith works like this:

.. code-block:: text

    Opt cycle 1: 
        Calculate truncated and full pointcharge field. Calculate gradient and energy correction.
    Opt cycle n: 
        if Opt cycle n is a multiple of TruncatedPC_recalc_iter then: 
            Recalculate correction using both full pointcharge field and truncated.
        else: 
            Use truncated PC field (defined by TruncPCRadius) in each QM run. Combine with energy and gradient corrections.
    Final Opt cycle: 
        Recalculate final geometry using full pointcharge field.

In a typical truncated-PC QM/MM optimization, the full pointcharge field (e.g. 1 million PCs) is used in the 1st step (expensive) but in later steps an approximated spherical PC-region (cheap) is used during the QM-steps (e.g. a spherical 35 Å radius region) 
until step 50/100/150 etc. (if TruncatedPC_recalc_iter=50) where the full pointcharge field is recalculated. When the optimization converges, e.g step 80, a final energy evaluation is performed using the full PC field.
For such an 80-iteration job, the full PC gradient may be calculated only 3 times (instead of 80 times) that can result in considerable time savings.

Note that QM and QM/MM energies are approximate during the optimization steps where a truncated PC field is used. The final energy is always calculated using the full PC field.
The error from the approximation depends on the TruncPCRadius parameter (smaller values than 30 not recommended) and TruncatedPC_recalc_iter (how often the full PC field is used). If TruncatedPC_recalc_iter=1 then no truncation is performed.

.. code-block:: python

    #QM/MM theory object defined with the truncated PC approximation
    qmmm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=frag, embedding="Elstat", qmatoms=qmatoms, printlevel=2,
        TruncatedPC=True, TruncPCRadius=35, TruncatedPC_recalc_iter=50)

  
#############################################
Example: QM/MM with ORCA and NonbondedTheory
#############################################

Example for a H2O-MeOH system where the MeOH is described by QM and H2O by MM. 
Here we read in a forcefield-file containing a nonbonded forcefield (see :doc:`MM-interfaces`).
The files for this example are available in the `examples/QM-MM-examples/QM-MM-ORCA-nonbondedtheory <https://github.com/RagnarB83/ash/tree/master/ash/examples/QM-MM-examples/QM-MM-ORCA-nonbondedtheory>`_ directory of the ASH repository.  


.. code-block:: python

    from ash import *

    #H2O...MeOH fragment defined. Reading XYZ file
    H2O_MeOH = Fragment(xyzfile="h2o_MeOH.xyz")

    # Specifying the QM atoms (3-8) by atom indices (MeOH). The other atoms (0,1,2) is the H2O and MM.
    #IMPORTANT: atom indices begin at 0.
    qmatoms=[3,4,5,6,7,8]

    # Charge definitions for whole fragment. Charges for the QM atoms are not important (ASH will always set QM atoms to zero)
    atomcharges=[-0.8, 0.4, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    #Defining atomtypes for whole system
    atomtypes=['OT','HT','HT','CX','HX', 'HX', 'HX', 'OT', 'HT']

    #Read forcefield (here containing LJ-part only) from file
    MM_forcefield=MMforcefield_read('MeOH_H2O-sigma.ff')

    #QM and MM objects
    ORCAQMpart = ORCATheory(orcasimpleinput="!BP86 def2-SVP def2/J tightscf", orcablocks="")
    MMpart = NonBondedTheory(charges = atomcharges, atomtypes=atomtypes, forcefield=MM_forcefield, 
        LJcombrule='geometric', codeversion="py")
    QMMMobject = QMMMTheory(fragment=H2O_MeOH, qm_theory=ORCAQMpart, mm_theory=MMpart, qmatoms=qmatoms,
                            charges=atomcharges, embedding='Elstat')

    #Single-point energy calculation of QM/MM object
    result = Singlepoint(theory=QMMMobject, fragment=H2O_MeOH, charge=0, mult=1)

    print("Single-point QM/MM energy:", result.energy)

    #Geometry optimization of QM/MM object (this may not converge)
    result2 = Optimizer(fragment=H2O_MeOH, theory=QMMMobject, coordsystem='tric', ActiveRegion=True, actatoms=[3,4,5,6,7,8], charge=0, mult=1)
    print("Optimized QM/MM energy:", result2.energy)


##########################################
Example: QM/MM with ORCA and OpenMMTheory
##########################################

See also :doc:`QM-MM-protein`.

The files for this example (DHFR protein) are available in the `examples/QM-MM-examples/QM-MM-CHARMM-example <https://github.com/RagnarB83/ash/tree/master/ash/examples/QM-MM-examples/QM-MM-CHARMM-example>`_ directory of the ASH repository.  


.. code-block:: python

  from ash import *

  numcores=1

  #Defining path to dir containing forcefield files and coordinates
  forcefielddir="./"
  psffile=forcefielddir+"step3_pbcsetup.psf"
  topfile=forcefielddir+"top_all36_prot.rtf"
  prmfile=forcefielddir+"par_all36_prot.prm"
  xyzfile=forcefielddir+"coordinates.xyz"

  #Read coordinates from XYZ-file
  frag = Fragment(xyzfile=xyzfile)

  #Creating OpenMM object
  openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
      charmmprmfile=prmfile, periodic=True, periodic_cell_dimensions=[80.0, 80.0, 80.0, 90.0, 90.0, 90.0],
      do_energy_decomposition=True, autoconstraints=None, rigidwater=False)


  #Creating ORCATheory object
  ORCAinpline="! HF-3c tightscf"
  ORCAblocklines="""
  %maxcore 2000
  """
  #Create ORCA QM object. Attaching numcores so that ORCA runs in parallel
  orcaobject = ORCATheory(orcasimpleinput=ORCAinpline,
                          orcablocks=ORCAblocklines, numcores=numcores)

  #act and qmatoms lists. Defines QM-region (atoms described by QM) and Active-region (atoms allowed to move)
  #IMPORTANT: atom indices begin at 0.
  #Here selecting the side-chain of threonine
  qmatoms = [569,570,571,572,573,574,575,576]
  actatoms = qmatoms #Same active region as QM-region here


  # Create QM/MM OBJECT by combining QM and MM objects above
  qmmmobject = QMMMTheory(qm_theory=orcaobject, mm_theory=openmmobject, printlevel=2,
                          fragment=frag, embedding="Elstat", qmatoms=qmatoms)

  #Run geometry optimization using geomeTRIC optimizer and HDLC coordinates. Using active region.
  Optimizer(theory=qmmmobject, fragment=frag, ActiveRegion=True, actatoms=actatoms,
                      maxiter=500, coordsystem='hdlc', charge=0,mult=1)



