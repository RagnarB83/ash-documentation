Quantum Chemistry Program Interfaces based on integrals
==========================================================

Some quantum chemistry programs allow calculation of the SCF or a correlated wavefunction directly from integrals, avoiding the basis-set specific creation of integrals typically required.
Often this involves reading in the integrals from an FCIDUMP file.

Running quantum chemistry from integrals directly would allow some freedom in performing quantum chemistry with integrals from different sources.

ASH has some basic support for facilitating this.
Currently supported programs are:

- **MRCC** : Arbitrary order CC, explicit correlation CC (F12), local natural orbital CC etc.
- **pyscf** : open-source general QM program: SCF, CC, MCSCF etc. 
- **ccpy** : Specialized CC methods
- **Block2**: DMRG
- **Dice**: semi-stochastic heat-bath CI


######################################################
Read/write AO integrals
######################################################




########################################################################
Running SCF from integrals, get mf object (pySCF) and MO integrals
########################################################################

In order to run a correlated WF calculation we need integrals in the MO-basis which means that we need to first solve the MO problem, i.e. the SCF.
Here we utilize the open-source pySCF program for this purpose.
It is possible to create a pySCF mean-field object directly from AO integrals in pySCF, without ever defining the standard GTO basis set.
The ASH function *create_pyscf_mol_and_mf* is convenient for this purpose.

To do this, we need to specify how many electrons we have, the nuclear repulsion energy of the system and read in the integrals as Numpy arrays.
Here we read in the integrals from Numpy binary files.

.. code-block:: python

    from ash import *

    #Defining the basic parameters of the system
    nuc_repulsion_energy = 22.51817918808511 # Nuclear repulsion energy in Eh
    # Or use function: nuc_nuc_repulsion(coords, charges)  (coordinates in Angstrom)
    num_el=14

    scf_type="RHF" # RHF, UHF, ROHF
    mult = 1 # Spin multiplicity
    #num_corr_el = 2 # Number of active electrons (i.e. non-frozen) in post-HF.
    #num_corr_orbs = 22 # Number of active orbitals in post-HF
    #numocc = num_corr_el/2 #Number of occupied orbitals

    #AO integrals read as Numpy arrays from disk (Numpy binary format)
    one_el_integrals = np.load("../AO_CO_1el.npy")
    two_el_integrals = np.load("../AO_CO_2el.npy")
    overlap = np.load("../AO_CO_overlap.npy")

    num_basis=one_el_integrals.shape[0] #Number of basis functions

    #Create mol and mf objects from integrals
    mol, mf = create_pyscf_mol_and_mf(numel=num_el, mult=mult, nuc_repulsion_energy=nuc_repulsion_energy,
        one_el_integrals=one_el_integrals, two_el_integrals=two_el_integrals, overlap=overlap )

    #Run the mf object (i.e. solve the SCF) to get the MOs
    mf.kernel()
    # MO coefficients : mf.mo_coeff
    # MO occupations: mf.mo_occ
    # MO energies: mf.mo_energy

Once the meanfield object has been run and we have the MOs, we can perform the AO->MO integral transformation. This can also be accomplished by tools inside pySCF.

.. code-block:: python

    #2-el and 1-el integral transformation from MO to AO
    twoel_MObas = ao2mo.kernel(two_el_integrals, mf.mo_coeff)
    oneel_MObas = np.einsum("ap, ab, bq -> pq", mf.mo_coeff, one_el_integrals, mf.mo_coeff)

######################################################
Write MO integrals to disk
######################################################

Now that we have the integrals in the MO-basis it is convenient to write them to disk.
A common standard is the FCIDUMP format.

ASH features the function *ASH_write_integralfile* for this purpose. 

.. code-block:: python

    def ASH_write_integralfile(two_el_integrals=None, one_el_integrals=None, nuc_repulsion_energy=None, header_format="MRCC",
                                num_corr_el=None, filename=None, int_threshold=1e-16, scf_type="RHF", mult=None,
                                symmetry_option=0, orbsym=None):



For MRCC calculations (see next) we want the FCIDUMP file to have an MRCC-specific header and have the filename be fort.55 and so we would run the function like this: 

.. code-block:: python

    nuc_repulsion_energy = 22.51817918808511 # Nuclear repulsion energy in Eh
    num_el=14
    scf_type="RHF"
    mult = 1
    ASH_write_integralfile(two_el_integrals=twoel_MObas, one_el_integrals=oneel_MObas,
        nuc_repulsion_energy=nuc_repulsion_energy, header_format="MRCC",
        num_corr_el=num_corr_el, filename="fort.55", int_threshold=1e-16, scf_type=scf_type, mult=mult)


######################################################
Running CC from integrals (MRCC)
######################################################

To run a MRCC calculation directly from MO integrals we need 2 files: the integral-file called fort.55 and a special basic inputfile named fort.56.
We can write the inputfile like below where we have specified the excitation level to be 4 (corresponding to CCSDTQ), scf_type to be RHF etc.
To request the calculation of reduced density matrices we specify dens = 1.
We also need to specify the occupations to use for the WF calculation and define how many electrons should be correlated etc.

.. code-block:: python

    num_corr_el = 2 # Number of active electrons (i.e. non-frozen) in post-HF.
    numocc = num_corr_el/2 #Number of occupied orbitals
    num_basis=one_el_integrals.shape[0] #Number of basis functions
    occupations = [2.0 if c < numocc else 0.0 for i,c in enumerate(range(num_basis))]

    MRCC_write_basic_inputfile(occupations=occupations, filename="fort.56", scf_type="RHF",
                               ex_level=4, nsing=1, ntrip=0, rest=0, CC_CI=1, dens=1, CS=1,
                               spatial=1, HF=1, ndoub=0, nacto=0, nactv=0, tol=9, maxex=0,
                               sacc=0, freq=0.0000, symm=0, conver=0, diag=0, dboc=0, mem=1024)

The inputfile is written to disk as fort.56 and looks like this:

.. code-block:: text


It can be modified as desired.
Once both fort.56 and fort.55 have been created we can run MRCC directly.

.. code-block:: python

    from ash import *

    numcores = 1
    run_mrcc("/path/to/mrccdir","mrcc.out", "OMP", numcores)

The density matrices are then available.




######################################################
Running CC calculations from integrals (pySCF)
######################################################

pySCF has support for various correlated wavefunctions, including the ability to get RDM1 and RDM2 from a CCSD(T) wavefunction.
If we already have a pySCF mean-field object created, we can use the ASH wrapper around pySCF to create a PySCFTheory object and request a CC calculation.

.. code-block:: python

    from ash import *
    
    #Fragment
    coordsstring="""C 0 0 0
    O 0 0 1.128
    """
    frag = Fragment(coordsstring=coordsstring, charge=0, mult=1)


    # Create ASH PySCFTheory object via previously created mf object (see earlier)
    pyscfobject = PySCFTheory(scf_type="RHF", mf_object=mf, CC=True, CCmethod="CCSD(T)", do_pop_analysis=False, symmetry=None)

    Singlepoint(theory=pyscfobject, fragment=frag)





######################################################
Running DMRG from integrals (Block2)
######################################################

A DMRG calculation with the Block2 interface can be run from integrals in a few different ways:

1. Via PySCFTheory object.
If the pySCF mean-field object has been created (see earlier) we can create an ASH pySCFTheory object which is conveniently passed over to BlockTheory.

.. code-block:: python

    from ash import *
    
    #Fragment
    coordsstring="""C 0 0 0
    O 0 0 1.128
    """
    frag = Fragment(coordsstring=coordsstring, charge=0, mult=1)


    # Create ASH PySCFTheory object via previously created mf object (see earlier)
    pyscfobject = PySCFTheory(scf_type="RHF", mf_object=mf, do_pop_analysis=False, symmetry=None)

    #Block2 DMRG calculation via input pyscftheory object. Here we define a DMRG-CASCI(2,10) calculation from the input SCF orbitals
    blocktheory = BlockTheory(pyscftheoryobject=pyscfobject, active_space=[2,10],  macroiter=0,
        numcores=1, memory=50000, tol=1e-8, initial_orbitals='HF', block_parallelization='OpenMP',
        maxM=1000, singlet_embedding=True, DMRG_DoRDM=True, DMRG_DoRDM2=True)

    # Run
    Singlepoint(theory=blocktheory, fragment=frag)

2. Via MO-basis integrals in FCIDUMP format.

If we have the 1- and 2-electron integrals available in the MO-basis we can write an FCIDUMP file to disk and start a DMRG calculation directly from the FCIDUMP file.
We first need to get the integrals in the MO-basis. Here we first run an SCF via pySCF as before and then convert the integrals from AO to MO basis.

.. code-block:: python

    from ash import *

    
    #Fragment
    coordsstring="""C 0 0 0
    O 0 0 1.128
    """
    frag = Fragment(coordsstring=coordsstring, charge=0, mult=1)

    # Get 1- and 2-integrals in AO basis
    one_el_integrals = np.load("../AO_CO_1el.npy")
    two_el_integrals = np.load("../AO_CO_2el.npy")
    overlap = np.load("../AO_CO_overlap.npy")

    #Create mol and mf objects from integrals
    mol, mf = create_pyscf_mol_and_mf(numel=14, mult=1,
        nuc_repulsion_energy=22.51817918808511,
        one_el_integrals=one_el_integrals, two_el_integrals=two_el_integrals,
        overlap=overlap )


    # AO->MO basis transformation via pyscf
    #2-el and 1-el integral transformation from MO to AO
    from pyscf import ao2mo
    twoel_MObas = ao2mo.kernel(two_el_integrals, mf.mo_coeff)
    oneel_MObas = np.einsum("ap, ab, bq -> pq", mf.mo_coeff, one_el_integrals, mf.mo_coeff)

    ###########################################################
    # Write MO-basis integrals to disk as FCIDUMP-style file
    ###########################################################
    ASH_write_integralfile(two_el_integrals=twoel_MObas, one_el_integrals=oneel_MObas,
        nuc_repulsion_energy=22.51817918808511, header_format="FCIDUMP",
        num_corr_el=14, filename="system.fcidump", int_threshold=1e-16,
        scf_type='RHF', mult=1)


    #Block2 DMRG calculation via input FCIDUMP file. Here we define a DMRG-CASCI(2,10) calculation from the input SCF orbitals
    blocktheory = BlockTheory(fcidumpfile="system.fcidump", active_space=[2,10],  macroiter=0,
        numcores=1, memory=50000, tol=1e-8, initial_orbitals='HF', block_parallelization='OpenMP',
        maxM=1000, singlet_embedding=True, DMRG_DoRDM=True, DMRG_DoRDM2=True)

    # Run
    Singlepoint(theory=blocktheory, fragment=frag)


Once a DMRG calculation has been run with RDMs requested, the RDMs are accessible from the BlockTheory object. 

.. code-block:: python

    #Grab RDMs in MO basis and do RDM1-MO to AO conversion
    rdm1_MO = blocktheory.properties["rdm1_MO"]
    print("RDM1-MO:", rdm1_MO)
    from ash.functions.functions_elstructure import DM_MO_to_AO
    # MO coefficients come from mf object
    rdm1_AO = DM_MO_to_AO(rdm1_MO, mf.mo_coeff)
    print("rdm1_AO:", rdm1_AO)


######################################################
Running SHCI from integrals (Dice)
######################################################

A semistochastic heatbath CI calculation with the Dice interface can be run from integrals almost identically to the DMRG-Block route above.
Simple follow the examples above but define a DiceTheory like this:

.. code-block:: python

    #From pyscf object
    dicetheory = DiceTheory(pyscftheoryobject=pyscftheoryobject, SHCI_active_space=[2,10], 
        numcores=1, memory=50000, initial_orbitals='HF',
        SHCI_DoRDM=True, SHCI_DoRDM2=True)

    # Or FCIDUMP-file:
    dicetheory = DiceTheory(fcidumpfile="system.fcidump", SHCI_active_space=[2,10], 
        numcores=1, memory=50000, initial_orbitals='HF',
        SHCI_DoRDM=True, SHCI_DoRDM2=True)

Once a SHCI calculation has been run with RDMs requested, the RDMs are accessible from the BlockTheory object. 

.. code-block:: python

    #Grab RDMs in MO basis and do RDM1-MO to AO conversion
    rdm1_MO = dicetheory.properties["rdm1_MO"]
    print("RDM1-MO:", rdm1_MO)
    from ash.functions.functions_elstructure import DM_MO_to_AO
    # MO coefficients come from mf object
    rdm1_AO = DM_MO_to_AO(rdm1_MO, mf.mo_coeff)
    print("rdm1_AO:", rdm1_AO)


######################################################
Running CC from integrals (ccpy)
######################################################

ccpy is a coupled cluster program in Python (with important routines in Fortran) that includes a number of interesting specialized CC methods.
A coupled cluster calculation with the interface to the ccpy program can be run from integrals like this:

*From a PySCFTheory object:*

First create a mf object as shown before, then create a PySCFTheory object.

.. code-block:: python

    #Fragment
    coordsstring="""C 0 0 0
    O 0 0 1.128
    """
    frag = Fragment(coordsstring=coordsstring, charge=0, mult=1)

    # Create ASH PySCFTheory object via previously created mf object (see earlier)
    pyscfobject = PySCFTheory(scf_type="RHF", mf_object=mf, do_pop_analysis=False, symmetry=None)

    # A CCSD ccpyTheory object using pyscfobj as input
    ccpy_theory = ccpyTheory(method="CCSD", pyscftheoryobject=pyscfobj, frozencore=True,
                cc_tol=1e-10, numcores=1, cc_maxiter=300)

    result = Singlepoint(theory=ccpy_theory, fragment=frag)

*From an FCIDUMP file:*

If an FCIDUMP file in MO-basis is already available (see earlier) we can start a ccpy calculation directly from it.

.. code-block:: python

    #Fragment
    coordsstring="""C 0 0 0
    O 0 0 1.128
    """
    frag = Fragment(coordsstring=coordsstring, charge=0, mult=1)

    # A CCSD ccpyTheory object using FCIDUMP-file (MO-baiss) as input
    theory = ccpyTheory(method="CCSD", fcidumpfile="FCIDUMP-file", frozencore=True,
                cc_tol=1e-10, numcores=1, cc_maxiter=300)
    
    result = Singlepoint(theory=ccpy_theory, fragment=frag)