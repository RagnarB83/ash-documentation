==========================
QM/MM Theory
==========================

QM/MM in ASH is flexible as one can in principle combine various QM-theories with either NonbondedTheory or OpenMMTHeory.

To do QM/MM, one combines a defined QMtheory object (:doc:`QM-interfaces`) and an MMtheory object in QMMMTheory object
and then specifies which atoms are QM and which are MM and the type of QM-MM coupling.
Note that in contrast to a QMtheory object or an MMtheory object, we pass a fragment object to a QMMMTheory object so that
the QMMMTheory object can define the division of QM-region and MM-region.

######################################
QM/MM boundary treatment
######################################

If the QMregion-MMregion boundary is between two bonded atoms, then a boundary correction need to be applied.
In ASH this is treated by the popular linkatom method, combined with charge-shifting.
A hydrogen-linkatom is added to cap the QM-subsystem. The hydrogen linkatoms are only visible to the QM theory, not the MM theory.
Additionally to prevent overpolarization, the atom charge of the MMatom is shifted towards its neighbours and a dipole correction
applied by adding additional pointcharges. These pointcharges are only visible to the QM theory.


######################################
Example: QM/MM with ORCA and NonbondedTheory
######################################

Example for a H2O-MeOH system where the MeOH is described by QM and H2O by MM.
Here we read in a forcefield-file (see :doc:`MM-interfaces`)

.. code-block:: python

    from ash import *
    settings_ash.init()

    #H2O...MeOH fragment defined
    H2O_MeOH = Fragment(xyzfile="h2o_MeOH.xyz")

    # Specifying MeOH QM atoms. Rest: 0,1,2 is H2O and MM. Note: atom indices begin at 0.
    qmatoms=[3,4,5,6,7,8]

    # Charge definitions for whole fragment.
    atomcharges=[-0.8, 0.4, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    #Defining atomtypes for whole system
    atomtypes=['OT','HT','HT','CX','HX', 'HX', 'HX', 'OT', 'HT']

    #Read forcefield (LJ-part only) from file
    MM_forcefield=MMforcefield_read('MeOH_H2O.ff')

    #QM and MM objects
    ORCAQMpart = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput=orcasimpleinput, orcablocks=orcablocks)
    MMpart = NonBondedTheory(charges = atomcharges, atomtypes=atomtypes, forcefield=MM_forcefield, LJcombrule='geometric')
    QMMMobject = QMMMTheory(fragment=H2O_MeOH, qm_theory=ORCAQMpart, mm_theory=MMpart, qmatoms=qmatoms,
                            charges=atomcharges, embedding='Elstat')

    #Geometry optimzation of QM/MM object
    geomeTRICOptimizer(fragment=H2O_MeOH, theory=QMMMobject, coordsystem='tric', ActiveRegion=True, actatoms=[3,4,5,6,7,8])


######################################
Example: QM/MM with ORCA and OpenMMTheory
######################################

Simple example:

.. code-block:: python

    from ash import *

    #Cores to use
    numcores=16
    #Forcefield files
    forcefielddir="/home/bjornsson/ASH-vs-chemshell-protein/QM-MM/FeMoco-test1/forcefielddir/"
    topfile=forcefielddir+"top_all36_prot.rtf"
    parfile=forcefielddir+"par_all36_prot.prm"
    psffile=forcefielddir+"new-XPLOR-psffile.psf"

    #Read coordinates from XYZ-file
    frag = Fragment(xyzfile="system.xyz", conncalc=False)

    #act and qmatoms lists. Defines QM-region and Active-region
    qmatoms = [13,14,15,20,22]
    actatoms = [13,14,15,20,22,300,320,340]

    #Creating OpenMMobject using CHARMM forcefield files
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=parfile)

    #Creating ORCATheory object
    orcadir="/opt/orca_current"
    ORCAinpline="! TPSSh RIJCOSX  D3BJ SARC/J ZORA-def2-SVP ZORA tightscf slowconv"
    ORCAblocklines="""
    %maxcore 2000
    """
    #Create ORCA QM object. Attaching numcores so that ORCA runs in parallel
    orcaobject = ORCATheory(orcadir=orcadir, charge=0,mult=1, orcasimpleinput=ORCAinpline,
                            orcablocks=ORCAblocklines, nprocs=numcores)

    # Create QM/MM OBJECT by combining QM and MM objects above
    qmmmobject = QMMMTheory(qm_theory=orcaobject, mm_theory=openmmobject, printlevel=2
                            fragment=frag, embedding="Elstat", qmatoms=qmatoms)

    #Run geometry optimization using geomeTRIC optimizer and HDLC coordinates. Using active region.
    geomeTRICOptimizer(theory=qmmmobject, fragment=frag, ActiveRegion=True, actatoms=actatoms,
                        maxiter=500, coordsystem='hdlc')


Advanced example:

.. code-block:: python

    from ash import *
    import time

    #Cores to use
    numcores=16
    #Forcefield files
    forcefielddir="/home/bjornsson/ASH-vs-chemshell-protein/QM-MM/FeMoco-test1/forcefielddir/"
    topfile=forcefielddir+"top_all36_prot.rtf"
    parfile=forcefielddir+"par_all36_prot.prm"
    psffile=forcefielddir+"new-XPLOR-psffile.psf"

    #Fragment file
    #Read old-chemshell file
    #frag = Fragment(chemshellfile="system.c", conncalc=False)
    #Read XYZ-file
    frag = Fragment(xyzfile="system.xyz", conncalc=False)

    #act and qmatoms lists
    #Reading in from qmatoms and act files. Here offsetting indices by -1 (to switch from 1-based to 0-based indexing)
    qmatoms = read_intlist_from_file("qmatoms",offset=-1)
    actatoms = read_intlist_from_file("act",offset=-1)

    #Creating OpenMMobject via CHARMM files
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=parfile, printlevel=1, platform='CPU' )

    #ORCA
    orcadir="/opt/orca_current"
    ORCAinpline="! TPSSh RIJCOSX  D3BJ SARC/J ZORA-def2-SVP ZORA tightscf slowconv"
    ORCAblocklines="""
    %maxcore 2000

    %basis
    newgto Fe \"ZORA-def2-TZVP(-f)\" end
    newgto V \"ZORA-def2-TZVP(-f)\" end
    newgto S \"ZORA-def2-TZVP(-f)\" end
    end

    %basis
    NewGTO Mo  \"old-ZORA-TZVP\" end
    addGTO Mo
    F 1
      1   0.6554500000      1.0000000000
    end
    end

    %scf
    MaxIter 1500
    diismaxeq 20
    end

    """
    #Charge/mult
    charge=-5
    mult=4

    #Brokensym options
    brokensym=True
    HSmult=36
    #Atoms in system to flop
    atomstoflip=[17763,17764,17766]
    # Atoms to put special basis set on
    extrabasisatoms=[17778]
    #Create ORCA QM object
    orcaobject = ORCATheory(orcadir=orcadir, charge=charge,mult=mult, orcasimpleinput=ORCAinpline, orcablocks=ORCAblocklines,
                            brokensym=brokensym, HSmult=HSmult, atomstoflip=atomstoflip, nprocs=numcores, extrabasisatoms=extrabasisatoms,
                            extrabasis="ZORA-def2-TZVP")

    # Create QM/MM OBJECT
    qmmmobject = QMMMTheory(qm_theory=orcaobject, mm_theory=openmmobject,
        fragment=frag, embedding="Elstat", qmatoms=qmatoms, printlevel=2)

    #Run Single-point job
    Singlepoint(theory=qmmmobject, fragment=frag, Grad=True)

    #Run geometry optimization using geomeTRIC optimizer and HDLC coordinates
    #Only active-region passed to optimizer
    geomeTRICOptimizer(theory=qmmmobject, fragment=frag, ActiveRegion=True, actatoms=actatoms, maxiter=500, coordsystem='hdlc')