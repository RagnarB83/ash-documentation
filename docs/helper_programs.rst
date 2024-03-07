Helper-programs interfaces
======================================

ASH contains interfaces to various little helper programs with simple interfaces that are documented below.


####################################################################
DRACO: scaling of solvation radii based on geometry and charges
####################################################################

The Grimme group has developed an interesting modification of standard continuum solvation models
that simply involves scaling the atomic radii that are used to create the cavity in continuum solvation calculations for improved accuracy.
The scaling procedure is based on both the coordination-number as well as the atomic charge on each atom.
The article describing the method was published in J. Phys. Chem. Lett: https://pubs.acs.org/doi/10.1021/acs.jpclett.3c03551

The DRACO program is available in an open-source repository on `Github <https://github.com/grimme-lab/DRACO>`_ :

ASH features a simple convenient interface that easily allows one to get the DRACO-scaled atomic radii within an ASH calculation 
that can subsequently be used for a QM-continuum solvation calculation in any QM-program (assuming that the program supports manual specification of atomic radii).

The ASH-interface to DRACO is a simple function, **get_draco_radii**,  that can be called from within an ASH-Python script.
The function takes either an ASH Fragment or an XYZ-file of a molecule as input, 
additionally the total charge of the molecule needs to be specified (or found within the ASH fragment).
The solvent is by default water, the radii-types are by default 'cpcm' (other options: 'cosmo' or 'smd') and the charge-model is by default 'ceh' (other option is 'eeq').

.. code-block:: python

    def get_draco_radii(fragment=None, xyzfile=None, charge=None, dracodir=None, 
                    radii_type='cpcm', solvent='water', chargemodel='ceh'):


To use, you must first download the DRACO binary and make sure that it is available in the PATH environment variable when ASH is run (or specify the dracodir).

.. code-block:: python

    from ash import *
    #Define fragment: Here finding glycine from the ASH database
    fragment = Fragment(databasefile="glycine.xyz")

    #Call Draco to get the scaled CPCM atomic radii assuming a water solvent and using a CEH charge model
    draco_radii = get_draco_radii(fragment=fragment, radii_type='cpcm', solvent='water', chargemodel='ceh')

    #These are the scaled atomic radii for each atom (in the same order as the atoms in the fragment)
    print("draco_radii:", draco_radii)

ASH will call Draco to calculate the scaled atomic radii, an outputfile (draco.out) is written out, which can be 
inspected and ASH then grabs the radii and return as a list of floats. 

To more conveniently use DRACO-radii automatically in a calculation, 
you can combine a **get_draco_radii** call with a QM-continuum calculation. 
The ORCA interface in ASH is flexible enough to allow this (using the *cpcm_radii* keyword).

.. code-block:: python
    
    from ash import *
    fragment=Fragment(databasefile="glycine.xyz")
    draco_radii = get_draco_radii(fragment=fragment, radii_type="cpcm", solvent="water")

    #Define ORCA-CPCM-DFT calculation using manual radii (from DRACO-step)
    qm = ORCATheory(orcasimpleinput="! r2scan-3c tightscf CPCM(water)", cpcm_radii=draco_radii)

    #Singlepoint calculation
    Singlepoint(theory=qm, fragment=fragment)

The ORCA input file created by ASH will contain the scaled atomic radii in the CPCM section and the ORCA output can also be inspected
to make sure the new radii are being used.


####################################################################
DFT-D4 dispersion correction
####################################################################

It is usually convenient to utilize dispersion corrections as they have been implemented in the respective QM-programs (e.g. specify the ORCA built-in dispersion correction when defining the ORCATheory) but
sometimes the respective QM program has not implemented any dispersion corrections. 
Or perhaps more flexibility in the choice of dispersion correction is desired. 

ASH features a simple interface to the `DFT-D4 program <https://github.com/dftd4/dftd4>`_ by the Grimme group for such cases.
To install, see the Github page. Best option is probably to install via conda/mamba like this:

.. code-block:: bash

    mamba install dftd4-python

Once installed in the ASH Python environment you can either use the **calc_DFTD4** function or the DFTD4Theory class.

.. code-block:: python

    def calc_DFTD4(fragment=None, functional=None, Grad=True):

The function **calc_DFTD4** takes a fragment as input and the functional name (string) that needs of course to match the functional used by the QM_program.
It returns the DFTD4 energy and gradient.

If one, however, wants to use the DFTD4 interface to correct a QM-calculation that will be used for geometry optimization, frequencies, molecular dynamics etc. (i.e. anything beyond a single-point calculation)
then, it is necessary to use the DFTD4Theory class and then to combine it with the QM-theory using the WrapTheory class, see :doc:`module_Hybrid_Theory`.

.. code-block:: python

    class DFTD4Theory:
        def __init__(self, functional=None, printlevel=2, numcores=1):


Example below shows how to perform a geometry optimization using an ORCATheory object (defining a PBE calculation without dispersion correction) and the DFTD4 dispersion correction via the DFTD4 program.

.. code-block:: python

    from ash import *

    #Glycine fragment from database
    frag = Fragment(databasefile="glycine.xyz")

    #PBE/def2-SVP via ORCA (no dispersion correction)
    orca = ORCATheory(orcasimpleinput="! PBE def2-SVP tightscf")
    #DFTD4 dispersion correction using DFTD4 library
    dftd4 = DFTD4Theory(functional="PBE")
    #Combining the two theories using WrapTheory
    dft_plus_dftd4_theory = WrapTheory(theory1=orca, theory2=dftd4)

    #Calling the Optimizer function using the WrapTheory object as theory 
    Optimizer(theory=dft_plus_dftd4_theory, fragment=frag)


