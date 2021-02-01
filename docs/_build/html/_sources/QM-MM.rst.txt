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

The recommended way of using link atoms is to define the QM-MM boundary for two carbon atoms that are as non-polar as possible.
In the CHARMM forcefield one should additionally make sure that one does not make a QM-MM boundary through a charge-group (check topology file).

In rare cases you may want to prevent ASH from adding a linkatom for a specific QM-atom, e.g. if you are making unusual
QM-MM boundaries. This can be accomplished like below. Note, however, that the QM-MM bonded terms will still be included.

.. code-block:: python

    #Excluding QM-atom 5785 from linkatom-creation.
   qmmmobject = QMMMTheory(qm_theory=orcaobject, mm_theory=openmmobject, fragment=frag, embedding="Elstat",
            qmatoms=qmatoms, excludeboundaryatomlist=[5785])


#############################################
Example: QM/MM with ORCA and NonbondedTheory
#############################################

Example for a H2O-MeOH system where the MeOH is described by QM and H2O by MM.
Here we read in a forcefield-file (see :doc:`MM-interfaces`)

.. code-block:: python

    from ash import *
    settings_ash.init()

    #H2O...MeOH fragment defined
    H2O_MeOH = Fragment(xyzfile="h2o_MeOH.xyz")

    # Specifying MeOH QM atoms. Rest: 0,1,2 is H2O and MM.
    #IMPORTANT: atom indices begin at 0.
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


##########################################
Example: QM/MM with ORCA and OpenMMTheory
##########################################

See also :doc:`QM-MM-protein`.

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
    #IMPORTANT: atom indices begin at 0.
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
