==========================
QM/MM Theory
==========================

QM/MM in ASH is flexible. Note that currently only nonbonded QM/MM is available.

To do QM/MM, one combines a defined QMtheory object (:doc:`QM-interfaces`) and an MMtheory object in QMMMTheory object
and then specifies which atoms are QM and which are MM and the type of QM-MM coupling.
Note that in contrast to a QMtheory object or an MMtheory object, we pass a fragment object to a QMMMTheory object so that
the QMMMTheory object can define the division of QM-region and MM-region.

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
                            atomcharges=atomcharges, embedding='Elstat')

    #Geometry optimzation of QM/MM object
    geomeTRICOptimizer(fragment=H2O_MeOH, theory=QMMMobject, coordsystem='tric', ActiveRegion=True, actatoms=[3,4,5,6,7,8])

