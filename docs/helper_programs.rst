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

   
    def packmol_solvate(packmoldir=None,ligand_files=None, num_mols_ligands=None, solvent_files=None, solvents_ratio=None, tolerance=2.0, result_file="final", shape="box",
                    min_coordinates=[0.0, 0.0, 0.0], max_coordinates=[40.0, 40.0, 40.0],total_density=None):


Packmol needs to be either downloaded and compiled (see `Packmol releases <https://github.com/m3g/packmol/releases>`_ ) or installed via conda/mamba (see `conda-forge packmol package <https://anaconda.org/conda-forge/packmol>`_). 

.. note::

    Packmol compilation: Once the archive has been downloaded and extracted, you have to compile the program which should be as simply as entering the directory and typing 'make' (requires Fortran compiler on system).
    For compilation problems, see the Packmol documentation.
    Once compiled, you can either add the directory to the PATH environment variable or specify the directory when calling the ASH function.

The packmol interface in ASH, offers a lot of flexibility in creating complicated systems easily.
It can create a box with a mixture of ligands and solvents in a very flexible way (e.g. Multiple Ligands and multiple solvents).
The interface allows you to control the number of molecules of each type (ligands and solvents), the density of the system, the shape of the box (###currently supports only box shape####) and the dimensions of the box.
For a case of mixture of solvents (e.g. water and ethanol) the ratio of the solvents (e.g. 2:1) can be specified instead of manually calculating the number of molecules.
The program will automatically calculate the number of molecules of each solvent based on the total density desired. The total density should be specified in g/ml units.
The input files can be either PDB or XYZ files, and the output file will be written as 'final.pdb' file.



Simple *packmol_solvate* examples:


.. code-block:: python

    #Create a box with 5:1 water-ethanol mixture in a 60x60x60 Angstrom box
    packmol_solvate(solvent_files=["water.xyz","ethanol.pdb"],solvents_ratio=[5,1],
        min_coordinates=[0.0, 0.0, 0.0], max_coordinates=[60.0, 60.0, 60.0],total_density=1.0)

    #Create a box with 1 ligand (in center of the box by default) with a solvent mixture of Acetonitrile and DMF (5:1)
    packmol_solvate(ligand_files=["DUM.pdb"], num_mols_ligands=[1],solvent_files=["ACN.xyz","DMF.pdb"],solvents_ratio=[5,1],
        min_coordinates=[0.0, 0.0, 0.0], max_coordinates=[60.0, 60.0, 60.0],total_density=1.0)

   #Multiple ligands with a solvent mixture of Acetonitrile and DMF (5:1)
   packmol_solvate(ligand_files=["DUM.pdb","MOL.pdb"], num_mols_ligands=[1,2],solvent_files=["ACN.xyz","DMF.pdb"],solvents_ratio=[5,1],
        min_coordinates=[0.0, 0.0, 0.0], max_coordinates=[60.0, 60.0, 60.0],total_density=1.0)


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
Basis Set Exchange
####################################################################

The `Basis Set Exchange website <http://basissetexchange.org>`_ website is well-known in the community as a database of Gaussian basis sets
that can be downloaded for various elements and is exportable in various formats.
Perhaps less known is that a Python API is also available that allows one to extract the basis set via a Python-library.
See `basis_set_exchange repository <https://github.com/MolSSI-BSE/basis_set_exchange>`_ for details but in short the library can be installed using pip: pip install basis_set_exchange

It can then be used like this in an ASH Python script.

.. code-block:: python
    
    from ash import *

    # Fragment to be calculcated
    frag = Fragment(databasefile="acetone.xyz")

    #import basis_set exchange
    import basis_set_exchange as bse
    #Getting the def2-mTZVPP basis set and def2-mTZVPP J auxiliary basis set
    basisname='def2-mTZVPP'
    auxbasisname='def2-mTZVPP-RIJ'
    # Getting the basis and aux bases in NWChem format (used by NWChem and pySCF) for the desired elements
    basis = bse.get_basis(basisname, elements=frag.elems, fmt='nwchem')
    auxbasis= bse.get_basis(auxbasisname, elements=frag.elems, fmt='nwchem')
    #Writing basis set strings to files
    with open(basisname,'w') as f: f.write(basis)
    with open(auxbasisname,'w') as f: f.write(auxbasis)

The basis-set files can then be used in the respective QM-program interface (assuming the ASH interface supports reading the basis set from file).
Below we show how the basis-set files created above can be read into the ASH PySCF interface.

.. code-block:: python

    pyscf_r2scan = PySCFTheory(scf_type="RKS", functional="r2scan", basis_file=basisname, 
                    densityfit=True, auxbasis_file=auxbasisname)


####################################################################
DFT-D4 dispersion correction
####################################################################

It is usually the most convenient to utilize dispersion corrections as they are implemented in the respective QM-programs (e.g. specify the ORCA built-in dispersion correction when defining the ORCATheory) but
sometimes the respective QM program has not implemented dispersion corrections, or perhaps more flexibility in the choice of dispersion correction is desired. 

ASH features a simple interface to the `DFT-D4 program <https://github.com/dftd4/dftd4>`_ by the Grimme group for such cases.
The interface is based on the Python API and so should have no execution drawbacks due to I/O.
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


Example below shows how to perform a geometry optimization using an ORCATheory object (defining a PBE calculation without dispersion correction) and the DFTD4 dispersion correction via the DFTD4 program
by combining it into a WrapTheory object.

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


####################################################################
Geometrical Counter-Poise correction (gCP)
####################################################################

The geometrical counterpoise correction by Grimme and coworkers has been found to be useful for reducing the basis set superposition error (BSSE)
in small-basis DFT calculations. 
Unlike the regular counterpoise correction (CP) that requires multiple DFT calculations and ghost atoms (see),
the gCP correction, depending only on geometry, has effectively no computational cost and is thus highly cost-effective for combining with small DFT-basis protocols.
The gCP correction is part of composite methods such as HF-3c, PBEh-3c, r2SCAN-3c.

ASH features a basic interface to the gCP (see https://github.com/grimme-lab/gcp). A Python API is not yet available and so the interface does have some I/O.
To install, see the Github page for latest instructions. Best option is probably to install via conda/mamba like this:

.. code-block:: bash

    mamba install gcp-correction


Once installed in the ASH Python environment you can either use the **calc_gcp** function or the gcpTheory class.

.. code-block:: python

    def calc_gcp(fragment=None, xyzfile=None, current_coords=None, elems=None, functional=None, Grad=True):

The function **calc_gcp** takes an ASH fragment as input (or xyzfile or coordinates-array) and the functional name (string) 
that needs of course to match the functional used by the QM_program.
It returns the gCP energy and gradient.

If one, however, wants to use the gCP interface to correct a QM-calculation that will be used for geometry optimization, frequencies, molecular dynamics etc. (i.e. anything beyond a single-point calculation)
then, it is necessary to use the gcpTheory class and then to combine it with the QM-theory using the WrapTheory class, see :doc:`module_Hybrid_Theory`.

.. code-block:: python

    class gcpTheory:
        def __init__(self, functional=None, printlevel=2, numcores=1):


Example below shows how to perform a geometry optimization using an ORCATheory object (defining a plan PBE) and the gcp correction via the gcp program
by combining it into a WrapTheory object.

Counter-poise corrected PBE/def2-SVP:

.. code-block:: python

    from ash import *

    #Glycine fragment from database
    frag = Fragment(databasefile="glycine.xyz")

    #PBE/def2-SVP via ORCA (no dispersion correction)
    orca = ORCATheory(orcasimpleinput="! PBE def2-SVP tightscf")
    #gcp correction
    gcp_corr = gcpTheory(functional="PBE")
    #Combining the two theories using WrapTheory
    dft_plus_gcp = WrapTheory(theory1=orca, theory2=gcp_corr)

    #Calling the Optimizer function using the WrapTheory object as theory 
    Optimizer(theory=dft_plus_gcp, fragment=frag)

####################################################################
Combining DFT with both D4 dispersion and gCP correction
####################################################################

Sometimes one would of course like to include both D4 dispersion and gCP correction.
This can also be accomplished in ASH using WrapTheory which is convenient if the QM-code does not have an implementation of neither D4 or gCP.

The example below shows how the r2SCAN-3c method (contains both D4 and gCP corrections) can be defined by WrapTheory 
where the pure DFT-part is calculated using either ORCATheory or PySCFTheory but the D4 and gCP corrections via DFTD4 and gcp interfaces.

Importantly, a WrapTheory object can be used as input to almost any ASH job-type, including Optimizer, NumFreq, MolecularDynamics etc.

*Manual r2SCAN-3c via ORCA, D4 and gCP interfaces*

Here we show how we can combine an ORCATheory DFT calculation-object with the DFTD4Theory and gCPTheory objects using WrapTheory.

.. code-block:: python

    from ash import *

    #Acetone fragment from database
    frag = Fragment(databasefile="acetone.xyz")

    #r2SCAN/def2-mTZVPP via ORCA
    orca_r2scan = ORCATheory(orcasimpleinput="! r2SCAN def2-mTZVPP def2-mTZVPP/J printbasis tightscf noautostart")
    # gcp correction
    gcp_corr = gcpTheory(functional="r2SCAN-3c", printlevel=3)
    # D4 correction
    d4_corr = DFTD4Theory(functional="r2SCAN-3c", printlevel=3)

    #Combining the 3 theories using WrapTheory
    r2scan3c = WrapTheory(theories=[orca_r2scan, gcp_corr,d4_corr])

.. note:: Normally it would of course be easier to use ORCA to do the whole r2SCAN-3c calculation using the built-in r2SCAN-3c keyword.


*Manual r2SCAN-3c definition via pySCF, D4 and gCP interfaces*

Since the basis and auxiliary basis set used in the r2SCAN-3c method (def2-mTZVPP and def2-mTZVPP/J) is not yet built into pySCF,
we first have to get the basis set. Here we show how this can be accomplished using the basis-set-exchange Python API.
We then combine the PySCFTheory object with DFTD4Theory and gcpTheory objects like before.

.. code-block:: python

    from ash import *

    #Acetone fragment from database
    frag = Fragment(databasefile="acetone.xyz")

    #Getting the basis set used by the r2SCAN-3c method
    import basis_set_exchange as bse
    basisname='def2-mTZVPP'
    auxbasisname='def2-mTZVPP-RIJ'
    basis = bse.get_basis(basisname, elements=frag.elems, fmt='nwchem')
    auxbasis= bse.get_basis(auxbasisname, elements=frag.elems, fmt='nwchem')
    with open(basisname,'w') as f: f.write(basis)
    with open(auxbasisname,'w') as f: f.write(auxbasis)

    #Defining a pySCF r2SCAN calculation with density fitting and the basis sets above
    pyscf_r2scan = PySCFTheory(scf_type="RKS", functional="r2scan", basis_file=basisname, densityfit=True, auxbasis_file=auxbasisname)

    # gcp correction
    gcp_corr = gcpTheory(functional="r2SCAN-3c")
    # D4 correction
    d4_corr = DFTD4Theory(functional="r2SCAN-3c")

    #Combining the 3 theories using WrapTheory
    r2scan3c = WrapTheory(theories=[pyscf_r2scan, gcp_corr,d4_corr])

    #Calling the Singlepoint function using the WrapTheory object as theory
    res = Singlepoint(theory=r2scan3c, fragment=frag, Grad=True)
    #Or you can do:  Optimizer(theory=r2scan3c, fragment=frag)