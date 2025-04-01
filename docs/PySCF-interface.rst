PySCF interface
======================================

`PySCF <https://pyscf.org>`_ is a very powerful open-source quantum chemistry program (or rather library) with the entire outer interface written in Python and everything else in C, 
including the very powerful libcint integral library. More recently the interface supports GPU-based computations.

ASH features a pretty good interface to PYSCF that allows one to conveniently use the various powerful DFT and WFT based features in the program 
that can be combined with the geometry optimization, surface scans, NEB, numerical frequencies, QM/MM,  MD and metadynamics features of ASH.
Due to the nature of PySCF as a Python library (with an essentially unlimited number of options) it is difficult to extensively support 
every PySCF feature and ASH instead takes the approach of writing wrapper-code around the most useful features that makes it suitable for ASH workflows, QM/MM etc.
This makes it very easy to use PySCF within ASH for the most basic features but the drawback being that every single PySCF method cannot be supported.
The pySCF interface in ASH is also used as part of the interface to Block2 and Dice for DMRG, SHCI and QMC calculations.
If the use of special features inside PySCF are desired, you may have to use PySCF directly (i.e. outside ASH) or contact us about adding the feature in the ASH interface.

**List of features:**

- One can define DFT, TDDFT, coupled cluster, CASSCF calculations fairly easily with support for various options.
- Control over special options such as BS-DFT, stability analysis, fractional occupationnoncollinear DFT (GHF/GKS), x2C relativistic Hamiltonians.
- Maximum overlap delta-SCF (OO-DFT) calculations for finding excited SCF states.
- Coupled cluster with multiple orbital references. Frozen natural orbital option (FNO).
- CASSCF starting orbital options: natural orbitals via MP2, CCSD and CCSD(T). Automatic CAS via DMET_CAS and AVAS.
- Read and write checkpointfiles, Molden files and Cubefiles
- PySCFTheory interface compatible with BlockTheory and DiceTheory interfaces for DMRG and stochastic heat-bath CI calculations.
- Support for `LOSC PySCF plugin <https://github.com/Yang-Laboratory/losc>`_ and MCPDFT
- Dispersion corrections via either pyscf-dispersion library (D3,D4) or via `vdw-wrapper <https://github.com/ajz34/vdw>`_ (D3, D4, TS and MBD).
- Electrostatically embedded QM/MM including pointcharge gradient is enabled
- Support for GPU4PYSCF: running HF/DFT on the GPU instead of CPU.

**Limitations:**

- Not all PySCF features are supported by the interface.
- post-SCF gradient currently not yet available in the interface



**PySCFTheory class:**


.. code-block:: python

  class PySCFTheory:
      def __init__(self, printsetting=False, printlevel=2, numcores=1, label="pyscf", platform='CPU', GPU_pcgrad=False,
                    scf_type=None, basis=None, basis_file=None, cartesian_basis=None, ecp=None, functional=None, gridlevel=5, symmetry='C1',
                    guess='minao', dm=None, moreadfile=None, write_chkfile_name='pyscf.chk', do_pop_analysis=True,
                    noautostart=False, autostart=True, fcidumpfile=None,
                    soscf=False, damping=None, diis_method='DIIS', diis_start_cycle=0, level_shift=None,
                    fractional_occupation=False, scf_maxiter=50, scf_noiter=False, direct_scf=True, GHF_complex=False, collinear_option='mcol',
                    NMF=False, NMF_sigma=None, NMF_distribution='FD', stability_analysis=False,
                    BS=False, HSmult=None, atomstoflip=None,
                    TDDFT=False, tddft_numstates=10, NTO=False, NTO_states=None,
                    mom=False, mom_occindex=0, mom_virtindex=1, mom_spinmanifold=0,
                    dispersion=None, densityfit=False, auxbasis=None, auxbasis_file=None, sgx=False, magmom=None,
                    pe=False, potfile=None, filename='pyscf', memory=3100, conv_tol=1e-8, verbose_setting=4,
                    CC=False, CCmethod=None, CC_direct=False, frozen_core_setting='Auto', cc_maxcycle=200, cc_diis_space=6,
                    CC_density=False, cc_conv_tol_normt=1e-06, cc_conv_tol=1e-07,
                    MP2=False,MP2_DF=False,MP2_density=False, DFMP2_density_relaxed=False,
                    CAS=False, CASSCF=False, CASSCF_numstates=1, CASSCF_weights=None, CASSCF_mults=None,
                    CASSCF_wfnsyms=None, active_space=None, casscf_maxcycle=200,
                    frozen_virtuals=None, FNO=False, FNO_orbitals='MP2', FNO_thresh=None, x2c=False,
                    AVAS=False, DMET_CAS=False, CAS_AO_labels=None, APC=False, apc_max_size=(2,2),
                    cas_nmin=None, cas_nmax=None, losc=False, loscfunctional=None, LOSC_method='postSCF',
                    loscpath=None, LOSC_window=None,
                    mcpdft=False, mcpdft_functional=None,
                    PBC_lattice_vectors=None,rcut_ewald=8, rcut_hcore=6, radii=None):

.. list-table::
  :widths: 15 15 15 60
  :header-rows: 1

  * - Keyword
    - Type
    - Default value
    - Details
  * - ``printsetting``
    - string
    - None
    - Printsetting. if True: printing to standard output otherwise disk.
  * - ``label``
    - string
    - None
    - Optional label.
  * - ``printlevel``
    - integer
    - 2
    - Printlevel
  * - ``numcores``
    - integer
    - 1
    - Number of CPU cores used by PySCF
  * - ``platform``
    - string
    - 'CPU'
    - Whether to use CPU-platform or GPU-platform (requires GPU4PySCF)
  * - ``scf_type``
    - string
    - None
    - Type of SCF-determinant. Options: 'RHF','UHF','RKS','UKS'.
  * - ``basis``
    - string
    - None
    - Name of basis set, can be a string (e.g. 'def2-SVP', must be valid PySCF basis-name) or a dict with element-specific keys and value-strings (basis-set name).
  * - ``ecp``
    - string
    - None
    - Name of ECP, can be a string ('e.g. 'def2-SVP', must be valid PySCF ECP-name) or a dict with element-specific keys and value-strings (ECP name).
  * - ``functional``
    - string
    - None
    - Name of DFT functional (must be valid PySCF functional-name).
  * - ``symmetry``
    - Boolean
    - False
    - Use of point-group symmetry or not.
  * - ``guess``
    - string
    - 'minao'
    - SCF guess options: 'minao', 'atom', 'huckel', 'vsap','1e'
  * - ``gridlevel``
    - string
    - ''
    - Name of DFT functional (must be valid PySCF functional-name).
  * - ``soscf``
    - Boolean
    - False
    - Second-order SCF algorithm active or not.
  * - ``damping``
    - float
    - None
    - Value of damping during SCF.
  * - ``diis_method``
    - string
    - 'DIIS'
    - DIIS method option: 'DIIS', 'ADIIS', 'EDIIS'
  * - ``diis_start_cycle``
    - integer
    - 0
    - In which SCF iteration to start the DIIS.
  * - ``level_shift``
    - float
    - None
    - Value of level-shift (Eh) during SCF.
  * - ``fractional_occupation``
    - Boolean
    - False
    - Whether fractional occupation is active or not.
  * - ``fractional_occupation``
    - Boolean
    - False
    - Whether fractional occupation is active or not.
  * - ``scf_maxiter``
    - integer
    - 50
    - Max number of SCF iterations.
  * - ``direct_scf``
    - Boolean
    - True
    - Whether direct SCF algorithm (recalculation of integrals in each iteration) is active or not. False is faster for small systems.
  * - ``densityfit``
    - Boolean
    - False
    - Whether to use density-fitting (RI) for Coulomb integrals. Use with auxbasis keyword.
  * - ``auxbasis``
    - string
    - None
    - Name of auxiliary basis set to use in density-fitting approximation. Example: 'def2-universal-jfit'.
  * - ``auxbasis_file``
    - string
    - None
    - Name of file containing auxiliary basis set (NWChem format).
  * - ``sgx``
    - Boolean
    - False
    - Whether to use semi-numerical exchange approximation for HF-exchange integrals. Note: gradient is not available
  * - ``stability_analysis``
    - Boolean
    - False
    - Whether SCF stability_analysis (calculation of orbital Hessian) is active or not.
  * - ``dispersion``
    - string
    - None
    - Dispersion correction to use. Options: 'D3', 'D4', 'TS', 'MBD'. Requires pyvdw package.
  * - ``do_pop_analysis``
    - Boolean
    - False
    - Whether to do automatic Mulliken population analysis after each SCF or not.
  * - ``moreadfile``
    - string
    - None
    - Name of PySCF checkpoint-file to read in as orbital guess.
  * - ``write_chkfile_name``
    - string
    - None
    - Name of the checkpointfile to write after SCF converges.
  * - ``noautostart``
    - Boolean
    - False
    - If True, then orbitals are not read in from a checkpoint-file.
  * - ``magmom``
    - list
    - None
    - If scf_type is 'GHF' or 'GKS', choose magnetic moment: list of the initial collinear spins of each atom.
  * - ``GHF_complex``
    - Boolean
    - False
    - If scf_type is 'GHF' or 'GKS', whether complex orbitals are used or not.
  * - ``collinear_option``
    - string
    - 'mcol'
    - If scf_type is 'GHF' or 'GKS', collinear option: col, ncol, mcol           
  * - ``GHF_complex``
    - Boolean
    - False
    - If scf_type is 'GHF' or 'GKS', whether complex orbitals are used or not.
  * - ``BS``
    - Boolean
    - False
    - Whether to find broken-symmetry solution by spin-flipping. Requires HSmult, and atomstoflip.
  * - ``HSmult``
    - integer
    - None
    - BS option: High-spin multiplicity to flip spin from.
  * - ``atomstoflip``
    - list of integers
    - None
    - What atom indices to flip spin on. Example: atomstoflip=[0,5] for atoms 0 and 5
  * - ``TDDFT``
    - Boolean
    - False
    - Whether to TDDFT on top of SCF solution or not.
  * - ``tddft_numstates``
    - integer
    - 10
    - Number of TDDFT states calculated.
  * - ``x2c``
    - Boolean
    - False
    - Whether to use the X2C scalar relativistic Hamiltonian or not.
  * - ``CAS``
    - Boolean
    - False
    - Whether to use a complete active space (CAS) or not. See also CASSCF and active_space keywords below.
  * - ``CASSCF``
    - Boolean
    - False
    - For CAS: Whether CASSCF orbital optimization is active. If False, then CAS-CI.
  * - ``active_space``
    - list of integers
    - None
    - Active space definition (electrons in orbitals), e.g. active_space=[3,2] (3 electrons in 2 orbitals).
  * - ``casscf_maxcycle``
    - integer
    - 200
    - Maximum number of CASSCF iterations.
  * - ``mcpdft``
    - Boolean
    - False
    - Whether multiconfigurational pair density functional theory (MCPDFT) method is active or not. Requires CAS keywords.
  * - ``mcpdft_functional``
    - string
    - None
    - Name of MCPDFT functional.
  * - ``AVAS``
    - Boolean
    - False
    - Whether to use the AVAS method to find CAS active space. Requires CAS_AO_labels keyword.
  * - ``DMET_CAS``
    - Boolean
    - False
    - Whether to use the DMET_CAS method to find CAS active space. Requires CAS_AO_labels keyword.
  * - ``CAS_AO_labels``
    - list of strings
    - None
    - List of atom-orbital label strings to use in AVAS/DMET_CAS selection.  Example: ['Fe 3d', 'Fe 4d', 'C 2pz']
  * - ``cas_nmin/cas_nmax``
    - float
    - None
    - If selecting active space from MP2 natural orbitals cas_nmin/cas_nmax tresholds determine active space.
  * - ``pe``
    - Boolean
    - False
    - Whether to use polarizable embedding in PySCF via CPPE library.
  * - ``potfile``
    - string
    - ''
    - Name of potential file for in PySCF CPPE polarizable embedding
  * - ``filename``
    - string
    - 'pyscf'
    - Filename used for PySCF output
  * - ``memory``
    - integer
    - 3100
    - Memory (in MB) used by PySCF .
  * - ``conv_tol``
    - float
    - 1e-8
    - Convergence tolerance in Eh .
  * - ``verbose_setting``
    - int
    - 4
    - How verbose PySCF output is.
  * - ``CC``
    - Boolean
    - False
    - Whether to do coupled-cluster on top of SCF or not.
  * - ``CCmethod``
    - string
    - None
    - Type of CCSD-method. Options:'CCSD', 'CCSD(T)'. More options will be available.
  * - ``CC_direct``
    - Boolean
    - False
    - Whether to use integral-direct CC or not.
  * - ``cc_maxcycle``
    - integer
    - 20
    - Maximum number of CC iterations.
  * - ``frozen_core_setting``
    - string
    - 'Auto'
    - How frozen core is handled. The ASH-default option is 'Auto' which means that frozen core settings are chosen by ASH (mimics ORCA-settings).
  * - ``frozen_virtuals``
    - list
    - None
    - Optionally freeze selected virtual orbitals in CC calculation.
  * - ``FNO``
    - Boolean
    - False
    - Do frozen natural orbital coupled cluster using MP2 natural orbitals.
  * - ``FNO_thresh``
    - float
    - None
    - Optional threshold to choose virtual natural orbitals to be skipped, based on natural occupation (from MP2 occupations).
  * - ``losc``
    - Boolean
    - False
    - Whether to do localized orbital scaling correction or not.
  * - ``loscfunctional``
    - string
    - None
    - The functional used (affect parameters chosen)
  * - ``LOSC_method``
    - string
    - None
    - LOSC correction post-SCF or full SCF. Options: 'postSCF' or 'SCF'
  * - ``LOSC_window``
    - list of floats.
    - None
    - LOSC energy window, e.g. [-30,-10].
  * - ``loscpath``
    - string
    - None
    - Path to losc package.
  * - ``mom``
    - Boolean
    - False
    - Whether to enable the maximum overlap method for delta-SCF calculations.
  * - ``mom_virtindex``
    - integer
    - 1
    - Which relative virtual orbital index to move electron from HOMO into. Default is 1 (LUMO); choose 2 for LUMO+1 etc.
  * - ``mom_spinmanifold``
    - integer
    - 0
    - What spin manifold to do MOM-deltaSCF calculations in. Default is 0 (i.e. alpha)

################################################################################
Advanced: PySCFTheory methods
################################################################################

The PySCFTheory class includes several methods that can also be called on their own (if you know what you are doing!)

.. code-block:: python

  def create_mol(self, qm_elems, current_coords, charge, mult):

  def define_basis(self,basis_string_from_file=None):

  def create_mf(self):

  def determine_frozen_core(self,elems):

  def set_numcores(self,numcores):

  def cleanup(self):

  def print_orbital_en_and_occ(self,mo_energies=None, mo_occ=None):

  def write_orbitals_to_Moldenfile(self,mol, mo_coeffs, occupations, mo_energies=None, label="orbs"):

  #Write Cube files for orbital, density or MEP
  def cubegen_orbital(self, mol, name, coeffs, nx=60,ny=60,nz=60):
  def cubegen_density(self, mol, name, dm, nx=60,ny=60,nz=60):
  def cubegen_mep(self, mol, name, dm, nx=60,ny=60,nz=60):

  def calculate_natural_orbitals(self,mol, mf, method='MP2', CAS_AO_labels=None, elems=None, relaxed=False, numcores=1):

  def calculate_CCSD_natorbs(self,ccsd=None, mf=None):

  def calculate_CCSD_T_natorbs(self,ccsd=None, mf=None):

  def run_population_analysis(self, mf, unrestricted=True, dm=None, type='Mulliken', label=None, verbose=3):

  def run_stability_analysis(self):

  def stability_analysis_loop(self,mf,mos,maxcyc=10):

  def read_chkfile(self,chkfile):

  def setup_guess(self):

  def calc_losc(self):

  def run_SCF(self,mf=None, dm=None, max_cycle=None):

  def run_MP2(self,frozen_orbital_indices=None, MP2_DF=None):

  def run_MP2_density(self, mp2object, MP2_DF=None, DFMP2_density_relaxed=None):

  def run_CC(self,mf, frozen_orbital_indices=None, CCmethod='CCSD(T)', CC_direct=False, mo_coefficients=None):

  def run_CC_density(self,ccobject,mf):

  def get_dipole_moment(self, dm=None, label=None):

  def get_polarizability_tensor(self):

  def set_mf_scfconv_options(self):

  def set_mf_smearing(self):

  def set_dispersion_options(self):

  def set_DF_mf_options(self):

  def set_DFT_options(self):

  def set_printing_option_mf(self):

  def set_collinear_option(self):

  def set_frozen_core_settings(self, elems):

  def set_embedding_options(self, PC=False):

  def density_potential_inversion(self, dm, lambda_par=8, method='ZMP', DF=True):

  def run(self, current_coords=None, current_MM_coords=None, MMcharges=None, qm_elems=None,
          elems=None, Grad=False, PC=False, numcores=None, pe=False, potfile=None, restart=False, label=None,
          charge=None, mult=None):
  def prepare_run(self, current_coords=None, current_MM_coords=None, MMcharges=None, qm_elems=None,
            elems=None, Grad=False, PC=False, numcores=None, pe=False, potfile=None, restart=False, label=None,
            charge=None, mult=None):
  def actualrun(self, current_coords=None, current_MM_coords=None, MMcharges=None, qm_elems=None,
          elems=None, Grad=False, PC=False, numcores=None, pe=False, potfile=None, restart=False, label=None,
          charge=None, mult=None,pyscf=None ):

################################################################################
PySCF installation
################################################################################

The PySCF interface is library-based and requires a PySCF installation inside the Python environment, typically via Pip (pip install pyscf).

################################################################################
Parallelization
################################################################################

The PySCF parallelization is OpenMP thread-based. The numcores keyword is used to specify the number of threads available to PySCF.


################################################################################
Running pySCF on the GPU
################################################################################

Recently GPU acceleration has become available in pySCF through the addition of a plugin package GPU4PySCF.
See `pySCF documentation <https://pyscf.org/user/gpu.html>`_. 
GPU acceleration requires the installation of the gpu4pyscf Python package, see `Github repository <https://github.com/pyscf/gpu4pyscf>`_
and requires CUDA-enabled hardware (i.e. Nvidia GPUs).

Installation may be as simple as:

.. code-block:: shell

  #For CUDA platform 12 (change 12 to 11 if only CUDA 11.x is available)
  pip3 install gpu4pyscf-cuda12x
  pip3 install cutensor-cu12

but see Github repo instructions if you have problems.

Once installed, a PySCFTheory object with platform="GPU" will automatically run on the GPU.
HF and DFT calculations are primarily useful for running on the GPU and the speedup on a single GPU vs. a multicore CPU can be quite impressive.
See articles on GPU4PySCF: `Introducing GPU Acceleration into the Python-Based Simulations of Chemistry Framework <https://pubs.acs.org/doi/10.1021/acs.jpca.4c05876>`_ , 
`Enhancing GPU-acceleration in the Python-based Simulations of Chemistry Framework <https://arxiv.org/abs/2404.09452>`_


################################################################################
Using the interface
################################################################################

Typicall the pySCFTheory theory object is simply used as an input-theory object

**Regular pySCF on the CPU:**

.. code-block:: python

  from ash import *
  n2_singlet= Fragment(diatomic="N2", bondlength=1.09, charge=0, mult=1)
  #Initialization of the PySCFTheory object (restricted HF here)
  pyscf_object = PySCFTheory(basis="cc-pVDZ", scf_type='RHF')
  #Calling Singlepoint function
  Singlepoint(theory=pyscf_object, fragment=n2_singlet)

**Using pySCF on the GPU (requires GPU4PySCF):**

.. code-block:: python

  from ash import *
  n2_singlet= Fragment(diatomic="N2", bondlength=1.09, charge=0, mult=1)
  #Initialization of the PySCFTheory object (restricted HF here)
  pyscf_object = PySCFTheory(basis="cc-pVDZ", scf_type='RHF', platform="GPU")
  #Calling Singlepoint function
  Singlepoint(theory=pyscf_object, fragment=n2_singlet)



In more advanced usage of the interface you can also call individual methods of the PySCFTheory object.
This is considered expert-territory and is typically not recommended.

.. code-block:: python

  from ash import *

  frag  = Fragment(diatomic="N2", bondlength=1.09, charge=0, mult=1)

  #Initialization of the PySCFTheory object
  pyscf_object = PySCFTheory(basis="cc-pVDZ", scf_type='RHF')

  #Prepare pySCFTheory object for run: This defines the pyscf mol and mf objects internally
  #Also sets various options inside mf and mol object previously defined
  pyscf_object.prepare_run(elems=frag.elems, current_coords=frag.coords, charge=frag.charge, mult=frag.mult)
  #Setup guess for SCF
  pyscf_object.setup_guess()
  #Run SCF with optional density-matrix input (dm) and max-cycle input (here 0, i.e. no SCF)
  pyscf_object.run_SCF(dm=None, max_cycle=0) #HF-SCF
  #Print orbitals, population analysis and dipole
  pyscf_object.print_orbital_en_and_occ() #HF-SCF
  pyscf_object.run_population_analysis(pyscf_object.mf)
  pyscf_object.get_dipole_moment()
  #Run CC using frozen core
  fc_indices=pyscf_object.set_frozen_core_settings(frag.elems)
  pyscf_object.run_CC(frozen_orbital_indices=fc_indices, CCmethod='CCSD(T)')

################################################################################
Controlling restart and guess 
################################################################################

How an SCF calculations begins can be controlled in a different ways.
Internally the SCF guess is handled by the setup_guess method which can be called on its own (see above for example).
First it is checked whether the PySCFTheory object already contains a density matrix (dm) and if so, then this is used as the guess.
Next it is checked whether the moreadfile keyword has been specified (should contain the name of a pySCF checkpointfile, something.chk) 
and if so, then the orbitals from the checkpoint-file are used as the guess.
Next it is checked whether Auto-Start has been disabled (either noautostart=True, or autostart=False). Autostart is on by default which means that it will try to read a checkpoint file in the directory with the default filename ("pyscf.chk").
If so then a new orbital-guess is used (based on the guess keyword, defaults to 'minao'). Guess options are: ['minao', 'atom', 'huckel', 'vsap','1e'].

.. code-block:: python

  #Reading in a density matrix. some_dm should here be a Numpy array
  pyscf_obj = PySCFTheory(scf_type="RHF", basis="def2-SVP", dm=some_dm)
  #Reading in a checkpoint file using moreadfile
  pyscf_obj = PySCFTheory(scf_type="RHF", basis="def2-SVP", moreadfile="previous.chk")
  #Disabling autostart by autostart=False
  pyscf_obj = PySCFTheory(scf_type="RHF", basis="def2-SVP", autostart=False)
  #Changing guess to huckel
  pyscf_obj = PySCFTheory(scf_type="RHF", basis="def2-SVP", autostart=False, guess="huckel")
 

The SCF-control functionality above can be utilized to do special things such as performing non-selfconsistent calculations using
some energy functional (HF or KS-DFT) on some other set of orbitals or density matrix. 
This requires one to i) read in the orbitals (or the density matrix) and ii) turn off SCF iterations.
Performing a non-selfconsistent DFT calculation using HF orbitals/density is called HF-DFT (or sometimes density-corrected DFT, DC-DFT) in the literature.
An example for this is shown below.

*Non-selfconsistent calculation using another set of orbitals (here HF-DFT)*

.. code-block:: python
  
  #Here we do a HF-DFT calculation by running first a HF calculation 
  #and then using the HF density matrix as a guess for the DFT calculation
  from ash import *
  frag = Fragment(databasefile="h2o.xyz")
  #Run HF calculation from scratch 
  pySCF_HF = PySCFTheory(scf_type="RHF", basis="def2-SVP", autostart=False)
  Singlepoint(fragment=frag, theory=pySCF_HF)
  #Create DFT object and reading in HF density matrix, also setting scf_maxiter=0 to avoid SCF
  pyscf_DFT_HF = PySCFTheory(scf_type="RHF", basis="def2-SVP", autostart=False, functional="PBE", dm=pySCF_HF.dm, scf_maxiter=0)
  Singlepoint(fragment=frag, theory=pyscf_DFT_HF)


Sometimes in unrestricted SCF calculations, one wants to guide the SCF procedure to find a symmetry-broken solution.
This is typically performed in the context of broken-symmetry DFT to describe spin-coupled antiferromagnetic states.
This can be performed in the PySCF interface by specifying BS=True, setting the spin multiplicity of the high-spin state (HSmult) 
and specifying the atom(s) to flip (atomstoflip is a list of atom indices ).

*Broken-symmetry solution via spin-flipping a spin-center from the high-spin solution*

.. code-block:: python

  #Here we do a HF-DFT calculation by running first a HF calculation 
  #and then using the HF density matrix as a guess for the DFT calculation
  from ash import *

  #Specify a BS-DFT calculation by setting BS=True and HSmult=3 (high-spin multiplicity)
  pySCF_HF = PySCFTheory(scf_type="RHF", basis="def2-SVP", functional='PBE', 
      autostart=False, BS=True, HSmult=3, atomstoflip=[0])
  Singlepoint(fragment=frag, theory=pySCF_HF, charge=0, mult=1)

################################################################################
Controlling basis set and ECP
################################################################################

Typically it is easiest to specify the basis set as in the examples above using a string (e.g. 'def2-SVP').
This option will work as long as the basis set is available inside pySCF and you specify the correct name of the basis set (see pySCF documentation).

If one needs more flexibility, such as using different basis sets for different elements, 
then one can use a dictionary to specify the basis set for each element.

.. code-block:: python

  from ash import *

  frag= Fragment(diatomic="HF", bondlength=1.09, charge=0, mult=1)
  #Initialization of the PySCFTheory object (restricted HF here)
  basis_dict={'H':'cc-pVDZ', 'F':'def2-SVP'}
  pyscf_object = PySCFTheory(basis=basis_dict, scf_type='RHF')
  #Calling Singlepoint function
  Singlepoint(theory=pyscf_object, fragment=frag)

If one wants to use a custom basis set, e.g. something not available inside pySCF, then one can use the basis_file keyword to specify a file containing the basis set.
The basis-set file needs to contain the basis set for all elements and needs to be in the NWChem basis set format.
It is best to go to the `Basis Set Exchange <https://www.basissetexchange.org/>`_ and download the basis set in NWChem format.

.. code-block:: python

  from ash import *

  frag= Fragment(diatomic="HF", bondlength=1.09, charge=0, mult=1)
  #pySCF object using a basis-set file
  pyscf_object = PySCFTheory(scf_type="UHF", basis_file="bla.basis")
  #Calling Singlepoint function
  Singlepoint(theory=pyscf_object, fragment=frag)


################################################################################
SCF convergence 
################################################################################

In case of SCF convergence problems there are a few options available.
One involves modifying the initial guess (see above) or reading in orbitals from a previous calculation (see also above).

If that does not work there are a few other options available such as turning on second-order SCF (SOSCF), 
using damping, modifying DIIS start-cycle, using level-shifting, enabling fractional occupation as well as increasing max iterations.

Shown below are the relevant keywords with their default values:

.. code-block:: python

  PySCFTheory(...,soscf=False, damping=None, diis_method='DIIS', diis_start_cycle=0, level_shift=None,
                  fractional_occupation=False, scf_maxiter=50)


################################################################################
Controlling integral approximation for Coulomb and HF Exchange
################################################################################

Density fitting for Coulomb and HF Exchange integrals is implemented in pySCF, it is not on by default in the interface.
For HF and hybrid-DFT it is also possible to use semi-numerical exchange approximation for HF exchange integrals (similar to RIJCOSX in ORCA).

See https://pyscf.org/user/df.html for more details on what is available in pySCF.

.. code-block:: python

  #Density fitting for Coulomb integrals only (recommended for non-hybrid DFT)
  #Note: Selecting the efficient Coulomb-only auxiliary basis set here
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RKS', functional='BLYP',
        densityfit=True, auxbasis='weigend')
  #RIJK: i.e. Density fitting for both Coulomb and HF Exchange (applies if HF or hybrid functional).
  #Note: Here we let pySCF automatically choose the RIJK auxiliary basis set (which hopefully exists for the basis set)
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RKS', functional='BLYP',
        densityfit=True)
  #Density fitting for Coulomb and + semi-numerical Exchange for HF Exchange integrals
  #Note: Here choosing again the more efficient Coulomb-only auxiliary basis set by Weigend
  #Warning: no analytical gradient available for this option
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RKS', functional='BLYP',
        densityfit=False, auxbasis='weigend', sgx=True)


################################################################################
Typical Examples
################################################################################

**HF-SCF example:**

.. code-block:: python

  from ash import *

  n2_singlet= Fragment(diatomic="N2", bondlength=1.09, charge=0, mult=1)

  #Minimal PySCFTheory definitino: RHF calculation
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RHF')
  Singlepoint(theory=PySCFcalc, fragment=n2_singlet)

**DFT-SCF example:**

.. code-block:: python

  from ash import *

  n2_singlet= Fragment(diatomic="N2", bondlength=1.09, charge=0, mult=1)

  #Define PySCF theory: RKS-PBE0 hybrid-DFT calculation
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RKS', functional="PBE0", gridlevel=6,
    numcores=2, memory=3000, filename='pyscf', printsetting=False)

  Singlepoint(theory=PySCFcalc, fragment=n2_singlet)


**Unrestricted CCSD(T) example:**

.. code-block:: python

  from ash import *

  o2_triplet= Fragment(diatomic="O2", bondlength=1.2075, charge=0, mult=3)

  #PySCF with UHF SCF and CCSD(T) on top
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", numcores=2, scf_type="UHF", CC=True,
    CCmethod='CCSD(T)', memory=3000, filename='pyscf', printsetting=False)

  Singlepoint(theory=PySCFcalc, fragment=o2_triplet)

################################################################################
Dispersion corrected DFT in pySCF
################################################################################

The ASH interface to pySCF features 3 different ways of including a dispersion correction,
utilizing different wrapper libraries to dispersion corrections: either `pyscf-dispersion <https://github.com/pyscf/dispersion>`_ or 
`pyvdw <https://github.com/ajz34/vdw>`_ or ASH's own interface to the DFTD4 library.
It is important to be aware that all dispersion corrections need to be tuned for each density functional and this present possibilities for things to go wrong.

**Option 1.** 
The easiest and recommended option is to directly use the correct functional keyword that includes the dispersion-part.
This option requires first installation of the pyscf-dispersion library: pip install pyscf-dispersion which is the recommended D3/D4 dispersion option for the pySCF object.
Example valid *functional* keywords are: 

- 'wb97x-d4' (i.e. wB97X-V functional with D4 dispersion correction)
- 'wb97m-d4' (i.e. wB97M-V functional with D4 dispersion correction)
- 'b3lyp-d3bj' (i.e. B3LYP with the D3BJ dispersion correction, the recommended D3 version)
- 'b3lyp-d3zero' (i.e. B3LYP with the older D3Zero dispersion correction)

One can also use the much more expensive nonlocal dispersion correlation functionals based on VV10 (built-into pyscf):

- 'wb97x-v' (i.e. wB97X-V functional)
- 'wb97m-v' (i.e. wB97M-V functional)

**Option 2.**
This option makes use of the `pyvdw <https://github.com/ajz34/vdw>`_ library which supports additional dispersion options.
To use this option, the dispersion keyword must be set to either: 'D3', 'D4', 'TS' or 'MBD'.
A few libraries must first be installed: 

.. code-block:: shell

  pip install pyvdw
  pip install toml
  mamba install simple-dftd3
  mamba install dftd3-python
  mamba install libmd

This option is less recommended since it relies on more third-party libraries, and is primarily useful for supporting the specific TS or MBD dispersion corrections.
This option will only work for some supported functionals and one should specify the functional name without dispersion-part, i.e. functional='b3lyp' instead of 'b3lyp-d3'.

**Option 3.**

Since ASH features it's own interface to DFTD4, see :doc:`helper_programs` it is also possible to combine a PySCFTheory object defining a DFT calculation without dispersion with a DFTD4Theory object
in a WrapTheory object, see :doc:`module_Hybrid_Theory`. 
This is useful for flexibility in defining e.g. composite methods (such as r2SCAN-3c) but is a more manual approach and offers more possibilities of selecting the wrong method.
It is thus only recommended if you carefully check the correctness of the results.


################################################################################
Natural orbital calculations from various WF methods
################################################################################

Natural orbitals are defined as the orbitals that diagonalize a one-particle reduced density matrix (1-RDM).
If the 1-RDM is available from a WF calculation (e.g. MP2, CCSD, CASSCF) then it is possible to calculate the natural orbitals, which is a convenient orbital representation of a many-particle WF.
The naturals orbitals can be visualized (from its MO coefficients) and the natural occupations can be analyzed.

The natural orbitals can be calculated from MP2, CCSD, CCSD(T), CASSCF methods (any method in principle as long as the 1-RDM is available). Natural orbitals can also be calculated for DMRG and SHCI WFs (see Block and Dice interfaces for more information).
For MP2 and CC methods the 1-RDM is not automatically available, you have to request it (e.g. by *MP2_density=True* or *CC_density=True*), since it takes additional computational effort.
If a CC_density or MP2_density calculation is requested then natural orbitals are automatically calculated and written to a Molden file. Otherwise, the natural orbitals can also be written to a Molden file manually as shown below.


################################################################################
Write Molden files of orbitals
################################################################################

To get access to the orbitals or wavefunction of a pySCF calculation it is easiest to write the orbitals to a Molden file.
For WFT calculations it is best to calculate the natural orbitals (see above) and then write the MO-coefficients associated with the natural orbitals to a Molden file.
This can be accomplished in a few different ways.

**Option 1:** 

If we are running a simple DFT-calculation (or HF) like below, we can simply call the **pySCF_write_Moldenfile** function that takes the pySCFTheory object as input and writes the orbitals to a Molden file.
Note that the pySCFTheory object needs to have been run (i.e. the Singlepoint calculation is necessary) before.

.. code-block:: python
    
  from ash import *

  #Fragment
  n2_singlet= Fragment(diatomic="N2", bondlength=1.09, charge=0, mult=1)
  #pyscf object
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RKS', functional="PBE0", gridlevel=6,
      numcores=2, memory=3000, filename='pyscf', printsetting=False)
  #Singlepoint  calculation
  Singlepoint(theory=PySCFcalc, fragment=n2_singlet)

  #Write moldenfile
  pySCF_write_Moldenfile(pyscfobject=PySCFcalc, label="orbs")


**Option 2:** 

Another option is to use the **write_orbitals_to_Moldenfile** method of the PySCFTheory object.
Here we specify the internal mol and mf objects of the PySCFTheory object as arguments to the write_orbitals_to_Moldenfile method.
This allows some more flexibility as we could in principle change the input MO-coefficients, MO-occupations and MO-energies (simple numpy arrays).
This could e.g. be used to write out natural orbitals if we have diagonalized a correlated density matrix and we have access to the natural orbital MO-coefficients and natural orbital occupations.

.. code-block:: python

  from ash import *
  n2_singlet= Fragment(diatomic="N2", bondlength=1.09, charge=0, mult=1)
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RKS', functional="PBE0", gridlevel=6,
    numcores=2, memory=3000, filename='pyscf', printsetting=False)
  Singlepoint(theory=PySCFcalc, fragment=n2_singlet)

  #Using the write_orbitals_to_Moldenfile method of the PySCTheory object
  PySCFcalc.write_orbitals_to_Moldenfile(PySCFcalc.mol, PySCFcalc.mf.mo_coeff, PySCFcalc.mf.mo_occ, PySCFcalc.mf.mo_energy, label="orbs")

**Option 3:**

If you have already run a calculation and you have a PySCF checkpoint file available (only for HF or DFT calculations) then you can create a Molden-file directly from the checkpoint file
using the **make_molden_file_PySCF_from_chkfile** function (note: you need to import it first). Note that the ASH fragment still needs to be loaded and you have to specify the same basis set as used previously.

.. code-block:: python

  from ash import *
  from ash.interfaces.interface_pyscf import make_molden_file_PySCF_from_chkfile

  #Define fragment again (geometry needed)
  n2_singlet= Fragment(diatomic="N2", bondlength=1.09, charge=0, mult=1)
  make_molden_file_PySCF_from_chkfile(fragment=n2_singlet, basis="cc-pVDZ",
      chkfile="pyscf.chk",label="orbs2")



################################################################################
Multireference calculations (CASSCF, MCPDFT etc.)
################################################################################

CASSCF calculations are possible in the interface.
Calculations are controlled by the CAS keyword (CAS=True or False) and the CASSCF keyword (CASSCF=True or False).
If CAS=True but CASSCF=False then a CAS-CI calculation is performed (only CI, no orbital optimization).
If CAS=True and CASSCF=True then a CASSCF calculation is performed (both CI and orbital optimization).
The active space is selected by providing a list of n electrons in m orbitals: active_space=[n,m].
Additionally one can solve for multiple states (controlled by CASSCF_numstates keyword) 
and it is also possible to specify the multiplicities for each state (CASSCF_mults), weights of the states (CASSCF_weights keyword).

.. code-block:: python

  #CASSCF calculation for a single-state
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='RHF', 
          CAS=True, CASSCF=True, CASSCF_numstates=1, active_space=[6,5], casscf_maxcycle=200)

  #State-averaged CASSCF calculations for 3-roots with equal weights
  PySCFcalc = PySCFTheory(basis="cc-pVDZ", scf_type='UHF', CAS=True, CASSCF=True, 
      CASSCF_numstates=3, active_space=[6,5], CASSCF_mults=[1,3,5], CASSCF_weights=[0.33,0.33,0.33])


A regular HF-SCF-calculation is currently automatically performed and can currently not be avoided. 
However, the HF-orbital guess for the CASSCF calculation can be controlled in a few different ways.
The options are: i) reading in a checkpoint-file (moreadfile keyword), ii) use the AVAS automatic active space method (AVAS=True),
iii) use the DMET_CAS automatic active space method (DMET_CAS=True), iv) use the APC automatic active space method v) use automatic MP2 natural orbitals.

AVAS and DMET_CAS requires one to set CAS_AO_labels keyword which is a list of atom-orbital labels (e.g. ['Fe 3d', 'Fe 4d', 'C 2pz']).

MC-PDFT calculations are also possible (mcpdft and mcpdft_functional keywords) but has not been tested.


################################################################################
Excited state calculation examples
################################################################################

**TDDFT calculations with NTO analysis**

.. code-block:: python

  from ash import *

  cstring="""
  O 0.0 0.0  0.0
  H 0.0 -0.757 0.587
  H 0.0 0.757 0.587
  """
  frag = Fragment(coordsstring=cstring, charge=0, mult=1)
  pyscf = PySCFTheory(scf_type='RKS', basis='6-31G', functional='b3lyp', 
    TDDFT=True, tddft_numstates=10, NTO=True, NTO_states=[1,2])
  Singlepoint(theory=pyscf, fragment=frag)


The relevant TDDFT output is shown in the main ASH output like below.
Also note that additional output will be present in the pySCF outputfile (by default: pyscf.out)

.. code-block:: text

  postSCF is True
  Now running TDDFT (Num states: 10)
  ----------------------------------------
  TDDFT RESULTS
  ----------------------------------------
  TDDFT transition energies (eV): [ 7.81984875  9.9212168   9.95812916 12.38331843 14.75956804 18.1889349
  27.77290941 28.15925452 29.1502703  30.1015163 ]
  Transition dipoles: [[-2.45304512e-01  2.68057788e-15  6.69547081e-16]
  [-2.01237402e-16 -1.21055864e-14  6.29424552e-01]
  [ 2.25211670e-15 -5.66428336e-15  2.04238232e-14]
  [ 5.34022012e-16 -5.35950517e-01 -7.01298803e-15]
  [ 1.12422599e-16  1.06732201e+00  2.04454308e-15]
  [-8.19417866e-16  2.28946438e-14  7.35926479e-01]
  [ 3.14351405e-14  2.32109432e-15 -7.66443771e-16]
  [-6.61659079e-16  1.36249903e-15  1.55571253e-01]
  [-3.49120535e-01  2.89921400e-15  1.48016888e-15]
  [-5.21074382e-15 -4.48622759e-01 -1.51123214e-14]]
  Oscillator strengths (length): [1.15283538e-02 9.62964261e-02 1.10832607e-28 8.71453664e-02
  4.11929077e-01 2.41342606e-01 6.76438154e-28 1.66969721e-02
  8.70464866e-02 1.48425612e-01]
  Oscillator strengths (velocity): [4.05425305e-02 1.70258256e-01 1.76701202e-28 1.13797326e-01
  3.86383743e-01 2.15322494e-01 1.72972290e-28 1.70834214e-02
  2.70699046e-02 1.03027463e-01]

  NTO analysis for state 1
  Now doing NTO analysis for states: [1, 2]
  See pySCF outputfile (pyscf.out) for the NTO analysis
  Doing NTO for state: 1
  Writing
  Doing NTO for state: 2
  Writing


pyscf.out contains the following NTO output:

.. code-block:: text

  State 1: 7.8198 eV  NTO largest component 0.9998830310985499
      occ-NTO: 1.000000 (MO #5)
      vir-NTO: 0.999752 (MO #6)
  State 2: 9.92117 eV  NTO largest component 0.986473412324918
      occ-NTO: 0.999699 (MO #4)
      vir-NTO: 0.999874 (MO #6)

The NTO-orbitals can be visualized using the Molden-files created: here nto-td-1.molden, nto-td-2.molden

**delta-SCF calculation using Maximum Overlap Method:**

PySCF includes the maximum overlap method that can be used to perform orbital-optimized SCF calculations of excited states (sometimes called delta-SCF approach).
You simply specify the SCF-type, functional and basis set as usual and then specify mom=True and optionally mom_virtindex and mom_spinmanifold keywords.

PySCF will first calculated the ground-state SCF with a regular Aufbau electron configuration and will then modify the guess to move an electron
from the HOMO to the specified virtual orbital index (default is mom_virtindex=1 which corresponds to the LUMO) of spin-manifold 0 (alpha).
If the SCF-type is restricted (RKS/RHF/ROHF/ROKS) then a ROHF/ROKS calculation will be carried out for the excited SCF calculations.
If the SCF type is unrestricted (UKS/UHF) then a UKS/UHF calculation will be carried out.

.. code-block:: python

  from ash import *

  cstring="""
  O 0.0 0.0  0.0
  H 0.0 -0.757 0.587
  H 0.0 0.757 0.587
  """
  frag = Fragment(coordsstring=cstring, charge=0, mult=1)
  pyscf = PySCFTheory(scf_type='RKS', basis='6-31G', functional='b3lyp', mom=True, mom_virtindex=1, mom_spinmanifold=0)
  Singlepoint(theory=pyscf, fragment=frag)

The output will look like this:

.. code-block:: text

  ----------------------------------------
  DELTA-SCF RESULTS
  ----------------------------------------

  Ground-state SCF energy -76.34781084088975 Eh
  Excited-state SCF energy -76.06068587471486 Eh

  delta-SCF transition energy 7.812957454584829 eV

  Alpha electron occupation pattern of ground state : [1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0.]
  Beta electron occupation pattern of ground state : [1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0.]

  Alpha electron occupation pattern of excited state : [1. 1. 1. 1. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
  Beta electron occupation pattern of excited state : [1. 1. 1. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0.]


**delta-SCF calculation using Maximum Overlap Method:**

