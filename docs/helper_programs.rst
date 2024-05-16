Helper-programs interfaces
======================================

ASH contains simple interfaces to various little helper programs that are documented below.

####################################################################
Packmol
####################################################################

`Packmol <https://m3g.github.io/packmol/userguide.shtml>`_ is a program used to generate initial configurations of a system of molecules.
It is particularly useful for preparing an initial box of a particular solvent (or mixture) that can
then be used for classical or QM/MM molecular dynamics simulations.

ASH has a simple interface to Packmol that allows one to easily use the program to create either an XYZ or PDB-file of a solvent system.

.. code-block:: python

    def packmol_solvate(packmoldir=None, inputfiles=None, num_mols=None, tolerance=2.0, result_file="final", shape="box",
                    min_coordinates=[0.0, 0.0, 0.0], max_coordinates=[40.0, 40.0, 40.0],density=None):


Packmol needs to be downloaded and compiled. See `Packmol releases <https://github.com/m3g/packmol/releases>`_ for the latest version.
Once the archive has been downloaded and extracted, you have to compile the program which should be as simply as entering the directory and typing 'make' (requires Fortran compiler on system).
For compilation problems, see the Packmol documentation.
Once compiled, you can either add the directory to the PATH environment variable or specify the directory when calling the ASH function.

To use the program you have to specify input-coordinate files containing a single molecule each 
(e.g. water.pdb and ethanol.pdb) and the number of molecules of each type that you want in the final box. Alternatively to the number of molecules it is possible to specify the density of the system.
The inputfiles can be either PDB or XYZ files. 
The number of molecules should be a list of integers corresponding to the number of molecules of each type.
If the density is specified it should be in units of g/ml
The tolerance (default 2.0 Angstrom) is the minimum distance between any two atoms in the final box.
The shape of the system can be chosen (e.g. box or sphere) and the min and max coordinates of the box should then be specified
which will control the size of the box (and where in Cartesian space).

Simple *packmol_solvate* examples:


.. code-block:: python

    #Create a box of 5281 ethanol molecules in a 80x80x80 Angstrom box
    packmol_solvate(inputfiles=["ethanol.pdb"], num_mols=[5281],
        min_coordinates=[0,0,0], max_coordinates=[80,80,80])

    #Create a box of ethanol molecules corresponding to a density of 0.789 g/ml (will result in 5281 molecules)
    packmol_solvate(inputfiles=["ethanol.pdb"], density=0.789,
        min_coordinates=[0,0,0], max_coordinates=[80,80,80])

    #Create an ethanol-water mixture box (1000 and 500) 40x40x40 Angstrom box
    packmol_solvate(inputfiles=["water.pdb", "ethanol.pdb"], num_mols=[1000, 500],
        min_coordinates=[0,0,0], max_coordinates=[40,40,40])

The final coordinate file will be written out as 'final.xyz' or 'final.pdb' depending on the extension of the input files.
The packmol.out file can be inspected to see additional Packmol output.

The coordinate file can next be used to create an OpenMMTheory object for the purpose of running a classical molecular dynamics simulation.
It is best to use a PDB-file in this case (as it can define the topology). An XML-file corresponding to the forcefield also needs to be created.


Note that Packmol is capable of many more features than the ASH interface allows, use Packmol directly if you require a specific feature.

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


