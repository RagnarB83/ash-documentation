Block2 interface
======================================

Block2 is a specialized WF program dedicated to density matrix renormalization group (DMRG) calculations, developed by the Garnet Chan group at Caltech
See: `Block2 Github page <https://github.com/block-hczhai/block2-preview>`_  and `Block2 documentation <https://block2.readthedocs.io>`_

Block2 contains an efficient and highly scalable implementation of DMRG based on the Matrix Product Operator formalism.
Its predecessors were StackBlock and the older Block code.

The interface in ASH intended to be used together with pySCF.

**BlockTheory class:**

.. code-block:: python

  class BlockTheory:
      def __init__(self, blockdir=None, pyscftheoryobject=None, blockversion='Block2', filename='input.dat', printlevel=2,
                  moreadfile=None, initial_orbitals='MP2', memory=20000, frozencore=True, fcidumpfile=None, 
                  active_space=None, active_space_range=None, cas_nmin=None, cas_nmax=None, macroiter=0,
                  Block_direct=False, maxM=1000, tol=1e-10, scratchdir=None, singlet_embedding=False,
                  block_parallelization='OpenMP', numcores=1, hybrid_num_mpi_procs=None, hybrid_num_threads=None,
                  FIC_MRCI=False, SC_NEVPT2_Wick=False, IC_NEVPT2=False, DMRG_DoRDM=False,
                  SC_NEVPT2=False, SC_NEVPT2_Mcompression=None, label="Block"):


.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``blockdir``
     - string
     - None
     - Path to Block2 directory.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``label``
     - string
     - 'Block'
     - String-label used for some output
   * - ``numcores``
     - integer
     - 1
     - Number of CPU cores that Block2 will use (OpenMP or MPI parallelization).
   * - ``memory``
     - integer
     - 20000
     - Memory in MB for Block2.
   * - ``scratchdir``
     - string
     - '.'
     - The path to the scratch directory used by Block2. If None then the current directory is used.
   * - ``filename``
     - string
     - 'input.dat'
     - Name used for Block2 inputfile
   * - ``pyscftheoryobject``
     - PySCFTheory object
     - None
     - A PySCFTheory object defining a mean-field calculation.
   * - ``moreadfile``
     - string
     - None
     - Name of a PySCF checkpoint file (.chk) to be used as input orbitals for Block2.
   * - ``initial_orbitals``
     - string
     - 'MP2'
     - Type of input-orbitals for Block2. Options: ['RKS', 'UKS', 'RHF', 'UHF', 'MP2', 'CCSD','CCSD(T)', 'SHCI', 'AVAS-CASSCF', 'DMET-CASSCF','CASSCF'] .
   * - ``CAS_AO_labels``
     - list
     - None
     - For input-orbitals options (AVAS-CASSCF, DMET-CASSCF) this list selects the active space (see PySCFTheory docs).
   * - ``active_space``
     - list
     - None
     - Defining active space as n electrons in m orbitals: e.g. [2,4] for 2 electrons in 4 orbitals.
   * - ``active_space_range``
     - list
     - None
     - Defining active space as a range of orbital indices: e.g. [2,40] for orbitals 2-40.
   * - ``cas_nmin``
     - float
     - 1.999
     - Upper limit for selecting natural orbitals for the active space. Requires initial_orbitals to have natural occupations.
   * - ``cas_nmax``
     - float
     - 0.0
     - Lower limit for selecting natural orbitals for the active space. Requires initial_orbitals to have natural occupations.
   * - ``macroiter``
     - integer
     - 0
     - Max number of macro iterations. If > 0 then DMRG-CASSCF is performed (orbital optimization).
   * - ``maxM``
     - integer
     - 1000
     - Max number of renormalized states in the DMRG calculation.
   * - ``tol``
     - float
     - 1e-10
     - The tolerance used in the calculation
   * - ``singlet_embedding``
     - Boolean
     - False
     - Whether to use singlet embedding or not (see Block2 manual).
   * - ``block_parallelization``
     - string
     - 'OpenMP
     - The type of Block2 parallelization to use. Options: ['OpenMP', 'MPI', 'Hybrid'].
   * - ``hybrid_num_mpi_procs``
     - integer
     - None
     - Number of MPI processes to use in hybrid parallelization.
   * - ``hybrid_num_threads``
     - integer
     - None
     - Number of OpenMP threads per MPI process to use in 'Hybrid' parallelization.
   * - ``FIC_MRCI``
     - Boolean
     - False
     - Do FIC-MRCI calculation on top of DMRG-SCF reference WF.
   * - ``SC_NEVPT2``
     - Boolean
     - False
     - Do SC_NEVPT2 calculation on top of DMRG-SCF reference WF.
   * - ``SC_NEVPT2_Mcompression``
     - integer
     - None
     - Do M-compression of reference WF for SC_NEVPT2 calculation.
   * - ``SC_NEVPT2_Wick``
     - Boolean
     - False
     - Do SC-NEVPT2-Wick calculation on top of DMRG-SCF reference WF.
   * - ``IC_NEVPT2``
     - Boolean
     - False
     - Do IC_NEVPT2 calculation on top of DMRG-SCF reference WF.
   * - ``DMRG_DoRDM``
     - Boolean
     - False
     - Calculate 1- and 2-RDMs for DMRG-SCF reference WF.


################################
Installing Block2
################################

See up-to-date information on:
https://block2.readthedocs.io/en/latest/user/installation.html

A simple installation of Block2 using OpenMP threading parallelization can be performed using pip:

.. code-block:: text

    #Simple OpenMP parallelization (simplest to setup )
    pip install block2
    #Possible MPI or hybrid OpenMP/MPI parallelization as well (more complicated to setup)
    pip install block2-mpi

A pySCF installation is required (see :doc:`PySCF-interface`) and additionally the `DMRGSCF interaface <https://github.com/pyscf/dmrgscf>`_   plugin must be installed:

This can be accomplished like this:

.. code-block:: text

    pip install git+https://github.com/pyscf/dmrgscf


After some additional settings modification (ASH will prompt you) you should be ready to go.

################################
Using the interface
################################

Typically you first create a pySCFTheory object and then a BlockTheory object pointing to the pySCFTheory object.
The default settings for DMRG are mostly sensible with maxM being the most important parameter.

See the Block2 documentation for details on the theory and various options (not all may be implemented in the interface).
https://block2.readthedocs.io/en/latest/index.html

################################
Parallelization
################################

Parallelization of Block2 is possible via either OpenMP (easiest to use), MPI or hybrid OpenMP/MPI.
You need to provide the *numcores* keyword when creating the BlockTheory object for OpenMP parallelization or MPI parallelization.
For hybrid OpenMP/MPI parallelization you need to provide the *numcores*, *hybrid_num_mpi_procs* and *hybrid_num_threads* keywords.
numcores must be equal to hybrid_num_mpi_procs*hybrid_num_threads.


################################
Examples
################################

**Example 1: Block2 DMRG CI calculation**

.. code-block:: python

    from ash import *

    numcores=10
    #Fragment
    fragment = Fragment(xyzfile="al2h2_mp2.xyz", charge=0, mult=1)
    #PySCF object: RHF/cc-pVTZ mean-field calculation
    pyscfobject = PySCFTheory(basis="cc-pVTZ", numcores=numcores, scf_type='RHF', conv_tol=1e-9,memory=50000)
    
    #Block2 DMRG calculation
    blockcalc = BlockTheory(pyscftheoryobject=pyscfobject, cas_nmin=1.999, cas_nmax=0.0, macroiter=0,
        numcores=numcores, memory=50000, tol=1e-8, initial_orbitals='CCSD', block_parallelization='OpenMP', 
        maxM=1000, singlet_embedding=True, DMRG_DoRDM=False)
    #Now running Singlepoint job
    result = Singlepoint(fragment=fragment, theory=blockcalc)
    print(f"Block DMRG-M=1000: Energy: {result.energy}")

**Example 2: Block2 DMRG calculation with increasing M states**

.. code-block:: python

    from ash import *

    numcores=10
    #Fragment
    fragment = Fragment(xyzfile="al2h2_mp2.xyz", charge=0, mult=1)
    #PySCF object: RHF/cc-pVTZ mean-field calculation
    pyscfobject = PySCFTheory(basis="cc-pVTZ", numcores=numcores, scf_type='RHF', conv_tol=1e-9)
    
    #Looping over M values
    for M in [100,200,300,400,500,600,700,800,900,1000]:
        blockcalc = BlockTheory(pyscftheoryobject=pyscfobject, cas_nmin=1.999, cas_nmax=0.0, macroiter=0,
          numcores=numcores, memory=50000, tol=1e-8, initial_orbitals='CCSD', block_parallelization='OpenMP', 
          maxM=M, singlet_embedding=True, DMRG_DoRDM=False)

        #Now running Singlepoint job for each epsilon
        result = Singlepoint(fragment=fragment, theory=blockcalc)
        print(f"Block-DMRG M={M}: Energy: {result.energy}")
