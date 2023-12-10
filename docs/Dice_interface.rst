Dice interface
======================================

Dice is a specialized WF program developed by Prof. Sandeep Sharma.
See: `Dice Github page <https://github.com/sanshar/Dice>`_  and `Dice documentation <https://sanshar.github.io/Dice/>`_

Dice contains various methods such as SHCI, VMC, GFMC, DMC, FCIQMC, stochastic MRCI and SC-NEVPT2, and AFQMC calculations.
The program is best used together with PySCF (see :doc:`PySCF-interface`).
Primarily it is useful for conveniently performing semi-stochastic heat-bath CI calculations (SHCI) on large systems, allowing 
CAS-type calculations up to an active space of approximately 100 orbitals.

The program needs to be compiled and set-up together with PySCF.

**DiceTheory class:**

.. code-block:: python

    class DiceTheory:
        def __init__(self, dicedir=None, pyscftheoryobject=None, filename='input.dat', printlevel=2, numcores=1, 
                    moreadfile=None, initial_orbitals='MP2', CAS_AO_labels=None,memory=20000, frozencore=True,
                    SHCI=False, NEVPT2=False, AFQMC=False, label="Dice",
                    SHCI_stochastic=True, SHCI_PTiter=200, SHCI_sweep_iter= [0,3],
                    SHCI_DoRDM=False, SHCI_sweep_epsilon = [ 5e-3, 1e-3 ], SHCI_macroiter=0,
                    SHCI_davidsonTol=5e-05, SHCI_dE=1e-08, SHCI_maxiter=9, SHCI_epsilon2=1e-07, 
                    SHCI_epsilon2Large=1000, SHCI_targetError=0.0001, SHCI_sampleN=200, 
                    SHCI_nroots=1, SHCI_cas_nmin=1.999, SHCI_cas_nmax=0.0, SHCI_active_space=None, 
                    SHCI_active_space_range=None, SHCI_active_all_but_core=None,
                    Dice_SHCI_direct=None, fcidumpfile=None, Dice_refdeterminant=None,
                    QMC_trialWF=None, QMC_SHCI_numdets=1000, QMC_dt=0.005, QMC_nsteps=50, 
                    QMC_nblocks=1000, QMC_nwalkers_per_proc=5):


.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``dicedir``
     - string
     - None
     - Path to MRCC directory.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``label``
     - string
     - 'Dice'
     - String-label used for some output
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores that Dice will use (MPI parallelization).
   * - ``filename``
     - string
     - 'input.dat'
     - Name used for Dice inputfile
   * - ``pyscftheoryobject``
     - PySCFTheory object
     - None
     - A PySCFTheory object defining a mean-field calculation.
   * - ``moreadfile``
     - string
     - None
     - Name of a PySCF checkpoint file (.chk) to be used as input orbitals for Dice.
   * - ``initial_orbitals``
     - string
     - 'MP2'
     - Type of input-orbitals for Dice. Options: ['RKS', 'UKS', 'RHF', 'UHF', 'MP2', 'CCSD','CCSD(T)', 
            'SHCI', 'AVAS-CASSCF', 'DMET-CASSCF','CASSCF'] .
   * - ``CAS_AO_labels``
     - list
     - None
     - For input-orbitals options (AVAS-CASSCF, DMET-CASSCF) this list selects the active space (see PySCFTheory docs).
   * - ``frozencore``
     - Boolean
     - True
     - For AFQMC and NEVPT2 calculations (not SHCI!), this defines a frozen-core.
   * - ``SHCI``
     - Boolean
     - False
     - Whether to do an SHCI calculation or not
   * - ``NEVPT2``
     - Boolean
     - False
     - Whether to do a NEVPT2 calculation or not on top of the SHCI CAS WF.
   * - ``AFQMC``
     - Boolean
     - False
     - Whether to do an AFQMC calculation or not.
   * - ``SHCI_stochastic``
     - Boolean
     - True
     - For SHCI: Whether to do the stochastic PT contribution or not
   * - ``SHCI_PTiter``
     - integer
     - 200
     - For SHCI: How many PT iterations to do. Set to 0 to turn PT off.
   * - ``SHCI_sweep_iter``
     - list
     - [0,3]
     - Control the SHCI sweep iterations. See Dice documentation.
   * - ``SHCI_sweep_epsilon``
     - list
     - [ 5e-3, 1e-3 ]
     - The epsilon values to use in the beginning and end of SHCI sweep iterations.
   * - ``SHCI_DoRDM``
     - Boolean
     - False
     - Whether to calculate the density matrices or not.
   * - ``SHCI_macroiter``
     - integer
     - 0
     - SHCI macro iterations. If > 0 then SHCI-CASSCF is performed (orbital optimization).
   * - ``SHCI_davidsonTol``
     - float
     - 5e-05
     - SHCI tolerance for the final Davidson step
   * - ``SHCI_dE``
     - float
     - 1e-08
     - SHCI energy convergence tolerance for the variational part.
   * - ``SHCI_maxiter``
     - float
     - 1e-08
     - SHCI max iterations for the variational part. 
   * - ``SHCI_epsilon2``
     - float
     - 1e-07
     - SHCI: Lower limit for accepting determinants in the selection. Also controls stochastic PT part.
   * - ``SHCI_epsilon2Large``
     - float
     - 1000
     - SHCI: Activates stochastic PT. Specifies the limit for what PT part will be done deterministically (rest stochastically).
   * - ``SHCI_targetError``
     - float
     - 0.0001
     - SHCI: Target standard deviation of the (semi-)stochastic-corrected energy.
   * - ``SHCI_sampleN``
     - integer
     - 200
     - SHCI: Number of times the determinants outside variational space are sampled in each stochastic iteration.
   * - ``SHCI_nroots``
     - integer
     - 1
     - SHCI: Number of states to solve for.
   * - ``SHCI_active_space``
     - list
     - None
     - SHCI: Defining active space as n electrons in m orbitals: e.g. [2,4] for 2 electrons in 4 orbitals.
   * - ``SHCI_active_space_range``
     - list
     - None
     - SHCI: Defining active space as a range of orbital indices: e.g. [2,40].
   * - ``SHCI_active_all_but_core``
     - list
     - None
     - SHCI: Experimental: 
   * - ``SHCI_cas_nmin``
     - float
     - 1.999
     - SHCI: Upper limit for selecting natural orbitals for the active space. Requires initial_orbitals to have natural occupations.
   * - ``SHCI_cas_nmax``
     - float
     - 0.0
     - SHCI: Lower limit for selecting natural orbitals for the active space. Requires initial_orbitals to have natural occupations.
   * - ``Dice_SHCI_direct``
     - Boolean
     - False
     - SHCI: Run Dice directly without pySCF or SHCI plugin. Requires FCIDUMP file and ref determinant (see below)
   * - ``fcidumpfile``
     - string
     - None
     - SHCI: Name of FCIDUMP file containing orbitals and integrals.
   * - ``Dice_refdeterminant``
     - string
     - None
     - SHCI: String specifying reference determinant
   * - ``QMC_trialWF``
     - string
     - None
     - QMC: Trial WF for QMC. Set to 'SHCI' to use a SHCI trial WF.
   * - ``QMC_SHCI_numdets``
     - integer
     - 1000
     - QMC: Number of determinants to use in the SHCI trial WF.
   * - ``QMC_dt``
     - integer
     - 0.005
     - QMC: Time step for QMC.
   * - ``QMC_nsteps``
     - integer
     - 50
     - QMC: Number of steps per block in QMC.
   * - ``QMC_nblocks``
     - integer
     - 1000
     - QMC: Number of blocks used in QMC.
   * - ``QMC_nwalkers_per_proc``
     - integer
     - 5
     - QMC: Number of walkers per MPI process for QMC.


################################
Installing Dice
################################

You need to download the Dice source code from `Dice Github page <https://github.com/sanshar/Dice>`_  
and compile it according to the Github instructions.
You also need to have installed pyscf (see :doc:`PySCF-interface`) and install the SHCI plugin:
https://github.com/pyscf/shciscf
After some additional settings modification (ASH will prompt you) you should be ready to go.

################################
Using the interface
################################

Typically you first create a pySCFTheory object and then a DiceTheory object pointing to the pySCFTheory object.
The default settings for SHCI are mostly sensible with epsilon being the most important parameter.

See the Dice documentation for more details on various settings:
https://sanshar.github.io/Dice/gettingstarted.html
https://sanshar.github.io/Dice/keywords.html

################################
Parallelization
################################

Parallelization of Dice is possible via MPI if you compiled it with MPI. 
Just provide the numcores keyword and the MPI environment will need to have been set (PATH, LD_LIBRARY_PATH to the MPI program).

################################
Examples
################################


**Example 1: Dice semi-stochastic heat-bath CI calculation**

.. code-block:: python

    from ash import *

    numcores=10
    #Fragment
    fragment = Fragment(xyzfile="al2h2_mp2.xyz", charge=0, mult=1)
    #PySCF object: RHF/cc-pVTZ mean-field calculation
    pyscfobject = PySCFTheory(basis="cc-pVTZ", numcores=numcores, scf_type='RHF', conv_tol=1e-9,memory=50000)
    
    #DiceTheory object for an SHCI calculation.
    #The pySCFTheory mean-field object has to be provided. The active space is generated using CCSD natural orbitals as input and 
    #selecting natural orbitals with occupations between 1.999 and 0.0.
    eps=1e-4
    dicecalc = DiceTheory(pyscftheoryobject=pyscfobject, numcores=numcores, SHCI=True, memory=50000,
                initial_orbitals='CCSD', SHCI_cas_nmin=1.999, SHCI_cas_nmax=0.0, SHCI_stochastic=True, 
                SHCI_PTiter=400, SHCI_sweep_iter= [0,3,6],SHCI_sweep_epsilon = [ 4*eps,2*eps,eps ], 
                SHCI_davidsonTol=1e-8, SHCI_epsilon2=1e-8, SHCI_epsilon2Large=1e-5, SHCI_macroiter=0)

    #Now running Singlepoint job
    result = Singlepoint(fragment=fragment, theory=dicecalc)
    print(f"Dice eps={eps}: Energy: {result.energy}”)

**Example 2: Dice SHCI calculation with decreasing epsilon parameter**

.. code-block:: python

    from ash import *

    numcores=10
    #Fragment
    fragment = Fragment(xyzfile="al2h2_mp2.xyz", charge=0, mult=1)
    #PySCF object: RHF/cc-pVTZ mean-field calculation
    pyscfobject = PySCFTheory(basis="cc-pVTZ", numcores=numcores, scf_type='RHF', conv_tol=1e-9)
    
    #Looping over epsilon values
    for eps in [1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6]:
        dicecalc = DiceTheory(pyscftheoryobject=pyscfobject, numcores=numcores, SHCI=True, memory=50000,
                    initial_orbitals='CCSD', SHCI_cas_nmin=1.999, SHCI_cas_nmax=0.0, SHCI_stochastic=True, 
                    SHCI_PTiter=400, SHCI_sweep_iter= [0,3,6], SHCI_sweep_epsilon = [ 4*eps,2*eps,eps ], 
                    SHCI_davidsonTol=1e-8, SHCI_epsilon2=1e-8, SHCI_epsilon2Large=1e-5, SHCI_macroiter=0)

        #Now running Singlepoint job for each epsilon
        result = Singlepoint(fragment=fragment, theory=dicecalc)
        print(f"Dice eps={eps}: Energy: {result.energy}”)
