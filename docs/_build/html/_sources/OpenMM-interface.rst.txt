OpenMM interface
======================================

The OpenMMTheory class:

.. code-block:: python

    class OpenMMTheory:
        def __init__(self, pdbfile=None, platform='CPU', active_atoms=None, frozen_atoms=None,
                     CHARMMfiles=False, psffile=None, charmmtopfile=None, charmmprmfile=None,
                     GROMACSfiles=False, gromacstopfile=None, grofile=None, gromacstopdir=None,
                     Amberfiles=False, amberprmtopfile=None, printlevel=2, nprocs=1,
                     xmlfile=None, periodic=False, periodic_cell_dimensions=None, customnonbondedforce=False):


Example creation of an OpenMMtheory object with CHARMM-files:

.. code-block:: python

    forcefielddir="/path/to/dir"
    topfile=forcefielddir+"top_all36_prot.rtf"
    parfile=forcefielddir+"par_all36_prot.prm"
    psffile=forcefielddir+"new-XPLOR-psffile.psf"
    openmmobject = OpenMMTheory(CHARMMfiles=True, psffile=psffile, charmmtopfile=topfile,
                               charmmprmfile=parfile)

Example creation of an OpenMMtheory object with GROMACS-files:

.. code-block:: python

    openmmobject = OpenMMTheory(GROMACSfiles=True, gromacstopdir="/path/to/gromacstopdir",
                    gromacstopfile="gromacstopfile.top", grofile="grofile.gro")

Example creation of an OpenMMtheory object with AMBER files:

.. code-block:: python

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile="/path/to/amberprmtopfile")

Example creation of an OpenMMtheory object with OpenMM XML file:

.. code-block:: python

    openmmobject = OpenMMTheory(xmlfile="exampl.xml")


An openmmtheory object can then be used to create a QM/MM theory object. See :doc:`module_QM-MM` page.