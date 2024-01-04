CP2K interface
======================================

`CP2K <https://www.cp2k.org>`_  is a very popular and powerful periodic electronic structure program, particulary known for its unique mixed
Gaussian and plane waves (GPW) approach and efficient parallelization which enables large-scale DFT MD simulations to be carried out.

ASH features a reasonably flexible interface to CP2K allowing the use of CP2K both as a QM-only code within ASH or as part of QM/MM Theory in ASH.
Energies and gradients are available in the interface so the CP2KTheory class in ASH can be used for single-point energies, geometry optimizations, 
numerical frequencies, surface scans, NEB and molecular dynamics. Neither the CP2K optimizer,MD or NEB routines are used, instead ASH handles all of this.
Both periodic and non-periodic systems are supported to some extent.
Additionally pointcharge-embedding is supported by the interface via the GEEP (Gaussian Expansion of Electrostatic Potential) approach in CP2K. 
This allows QM/MM calculations with CP2K as QM-code and OpenMM as MM-code within ASH.

If the purpose is to primarily carry out periodic DFT MD simulations, running CP2K via the ASH interface will not offer many benefits over CP2K directly, 
however, running geometry optimizations, surface scans and NEB via ASH may be more convenient than using CP2K directly.
Furthermore the QM/MM capabilities within ASH and the flexible forcefield support by OpenMM may be preferable to the CP2K options.

**CP2KTheory class:**

.. code-block:: python
    
    class CP2KTheory:
        def __init__(self, cp2kdir=None, filename='cp2k', cp2k_bin_name=None, printlevel=2, basis_dict=None, potential_dict=None, label="CP2K",
                    periodic=False, periodic_type='XYZ', qm_periodic_type=None,cell_dimensions=None, cell_vectors=None,
                    qm_cell_dims=None, qm_cell_shift_par=6.0, wavelet_scf_type=40,
                    functional=None, psolver='wavelet', potential_file='POTENTIAL', basis_file='BASIS',
                    basis_method='GAPW', ngrids=4, cutoff=250, rel_cutoff=60,
                    method='QUICKSTEP', numcores=1, center_coords=True, scf_convergence=1e-6,
                    coupling='COULOMB', GEEP_num_gauss=6, MM_radius_scaling=1, mm_radii=None):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``cp2kdir``
     - string
     - None
     - Directory where CP2K binaries are.
   * - ``cp2k_bin_name``
     - string
     - None
     - Name of CP2K binary to use. If not provided, ASH will search for cp2k.ssmp, cp2k.sopt, cp2k.psmp in the cp2kdir directory or PATH.
   * - ``printlevel``
     - integer
     - 2
     - Printlevel
   * - ``filename``
     - string
     - 'cp2k'
     - Name of CP2K inputfile created by ASH. 
   * - ``numcores``
     - integer
     - 1
     - The number of CPU cores used.
   * - ``basis_dict``
     - dict
     - None
     - Dictionary defining basis set information for each element.
   * - ``potential_dict``
     - dict
     - None
     - Dictionary defining pseudopotential information for each element.
   * - ``label``
     - string
     - 'CP2K'
     - Label for the CP2KTheory object.
   * - ``periodic``
     - Boolean
     - False
     - Whether to use periodic boundary conditions or not.
   * - ``periodic_type``
     - string
     - 'XYZ'
     - Type of PBC used. Options: 'XYZ', 'XY', 'XZ', 'YZ' etc.
   * - ``cell_dimensions``
     - list
     - None
     - The cell dimensions defining the box (applies to both PBC and no-PBC). List of cell-lengths in Angstroms and cell-angles in degrees.
   * - ``cell_vectors``
     - list of lists
     - None
     - Alternative way of specifying the cell dimensions. List of lists of cell-vectors in Angstroms.
   * - ``wavelet_scf_type``
     - integer
     - 40
     - Wavelet-only: Scaling function. Possible values: 8,14,16,20,24,30,40,50,60,100
   * - ``functional``
     - string
     - None
     - Name of DFT functional to use.
   * - ``psolver``
     - string
     - 'wavelet'
     - Name of Poisson solver to use. Options: 'PERIODIC', 'MT', 'WAVELET', 'MULTIPOLE'.
   * - ``potential_file``
     - string
     - 'POTENTIAL'
     - Name of potential file.
   * - ``basis_file``
     - string
     - 'BASIS'
     - Name of basis set file.
   * - ``basis_method``
     - string
     - 'GAPW'
     - Type of CP2K basis-set method to use. Options: 'GPW' (Gaussian-planewave with pseudopoentials), 'GAPW' (Gaussian augmented planewave). 
   * - ``ngrids``
     - integer
     - 4
     - Number of real-space grids to use.
   * - ``cutoff``
     - integer
     - 250
     - Planewave cutoff (in Ry) used.
   * - ``rel_cutoff``
     - integer
     - 60
     - Relative cutoff (in Ry) used. Controls which product Gaussians are mapped onto which level of the multi-grid.
   * - ``center_coords``
     - Boolean
     - True
     - Whether CP2K centers coordinates or not. Usually necessary.
   * - ``scf_convergence``
     - float
     - 1e-6
     - SCF convergence in Hartree.
   * - ``method``
     - string
     - 'QUICKSTEP'
     - The CP2K module to use.
   * - ``qm_periodic_type``
     - string
     - None
     - QM/MM only: Type of PBC used for QM-part. Options: 'XYZ', 'XY', 'XZ', 'YZ' etc.
   * - ``qm_cell_dims``
     - list
     - None
     - QM/MM only: List of cell-lengths in Angstroms for the QM region. Cell dimensions estimated if not provided.
   * - ``qm_cell_shift_par``
     - float
     - None
     - QM/MM only: Parameter used to shift the QM cell dimensions based on molecule size. Only used if qm_cell_dims is not provided.
   * - ``coupling``
     - string
     - 'GAUSSIAN'
     - QM/MM only: The type of QM/MM coupling used. Only 'GAUSSIAN' is supported for DFT. 'COULOMB' available for semi-empiricical systems.
   * - ``GEEP_num_gauss``
     - integer
     - 6
     - QM/MM only: Number of Gaussians used to expand each MM point charge.
   * - ``MM_radius_scaling``
     - float
     - 1.0
     - QM/MM only: Optional scaling factor of the MM radii.
   * - ``mm_radii``
     - dict
     - None
     - QM/MM only: Optional dictionary of MM radii for each element. If not provided, radii are estimated internal table.



################################################################################
CP2K installation
################################################################################

CP2K can be installed in several different ways, see https://www.cp2k.org/download
It is easiest to either download binaries (see link) or install via conda (see https://anaconda.org/conda-forge/cp2k).
Alternatively you can compile CP2K from source: https://github.com/cp2k/cp2k/blob/master/INSTALL.md

Note that downloaded or compiled CP2K binaries may come in a few different forms: e.g. cp2k.ssmp, cp2k.sopt, ccp2k.popt, cp2k.psmp 
where sopt means serial-optimized, ssmp means single-process with OpenMP, 
popt means parallel-optimized with MPI and psmp means parallel-optimized with MPI and OpenMP.

ASH will find a CP2K binary to use according to this logic:

1. if cp2kdir variable provided (containing path to where the binaries are) and cp2k_bin_name provided: use that binary in that directory
2. if cp2kdir variable provided but cp2k_bin_name NOT provided: search for cp2k.X executables in the cp2kdir directory
3. if cp2kdir variable NOT provided but cp2k_bin_name provided: search for cp2k_bin_name in PATH
4. if cp2kdir variable NOT provided and cp2k_bin_name not provided: search for cp2k.X executables in PATH

Note that the search for executables will only work if the binaries are named: cp2k.X where X is one of ssmp, sopt, popt, psmp.
ASH will search for executables in this order: ["cp2k.psmp", "cp2k.popt", "cp2k.ssmp","cp2k.sopt"]


################################################################################
Parallelization
################################################################################

CP2K binaries differ in their parallelization:
- sopt: no parallelization
- ssmp: uses OpenMP parallelization
- popt: uses MPI parallelization
- psmp: mixed MPI and OpenMP parallelization. Primarily useful for massive parallelization (>10K cores). 

The CP2K manual advises to use the cp2k.psmp executable as it is the most flexible.
The number of cores that CP2K will use for either MPI-runs or OpenMP runs is controlled by the numcores keyword in the CP2KTheory object.

Warning: Massively parallel CP2K within ASH has not been tested much.


################################################################################
Controlling the basis set
################################################################################

The primary purpose of using CP2K is probably to take advantage of the efficient mixed Gaussian and plane wave (GPW) approach where Gaussians are used to calculate
the 1-electron integrals and plane waves are used to calculate the 2-electron integrals.
Furthermore the user should specify whether the standard GPW (Gaussian and planewaves) or GAPW (Gaussian augmented GPW) method should be used.
Pseudopotential-based calculations can be performed with both methods, however, all-electron calculations can only be performed with GAPW.
GAPW may have more stable forces and require reduced cutoff but may be more expensive.

Depending on whether GPW or GAPW is used, suitable basis set and pseudopotential information should be provided.
This is controlled by defining the basis_dict and potential_dict keywords in the CP2KTheory object.
The chosen basis sets and pseudopotentials must be available in the specified basis and potential files.
For all-electron GAPW calculations one should set value for each element in the potential_dict to 'ALL'.

.. code-block:: python

    
    #Defining MOLOPT basis sets and GTH pseudopotentials for each element
    basis_dict={'C':'DZVP-MOLOPT-SR-GTH','O':'DZVP-MOLOPT-SR-GTH','H':'DZVP-MOLOPT-SR-GTH'}
    potential_dict={'C':'GTH-PBE-q4','O':'GTH-PBE-q6','H':'GTH-PBE-q1'}
    cp2k_object = CP2KTheory(basis_method='GPW', basis_dict=basis_dict,potential_dict=potential_dict, 
            potential_file='POTENTIAL', basis_file='BASIS',)

Note that if the specified basis-file or potential-file is not in the current dir (or parent dir) then ASH will automatically
copy a file containing GTH pseudopotentials (renamed from GTH_POTENTIALS to POTENTIAL) and MOLOPT basis sets (renamed from BASIS_MOLOPT to BASIS).
This will only work if MOLOPT basis sets are being used. For all other basis sets, then the user must provide the basis and potential files.

For the planewave part of the basis set, the cutoff and rel_cutoff keywords can be used to control the cutoffs used.
The number of grids also play a role in the accuracy of the calculation and can be controlled by the ngrids keyword (default=4).
Suitable cutoff values and grids require some experience or testing.
See https://www.cp2k.org/howto:converging_cutoff for some information on how to choose cutoffs and grids.
A reasonable value for the Cutoff is 250 Ry and a good value for the rel_cutoff is usually 60 Ry. These cutoff should be varied simultaneously.
These are the ASH defaults but we don't have a lot of experience with CP2K. 
Some system setups (depends on elements, basis set and pseudopotential) may require larger values and other systems will run more efficiently with smaller values.


################################################################################
Periodic vs. non-periodic calculations
################################################################################

CP2K is a code first and foremost developed for the purpose of periodic calculations. 
It is nonetheless possible to perform non-periodic calculations and this is probably preferable for calculations on molecules in vacuum (to avoid PBC artifacts) 
and may also be beneficial for some QM/MM applications within ASH.

Regardless of whether the system is periodic or not, the system cell needs to be specified.
The *cell_dimensions* (e.g. *cell_dimensions=[10.0,10.0,10.0,90.0,90.0,90.0]* or *cell_vectors* (*cell_vectors=[[10.0,0.0,0.0],[0.0,10.0,0.0],[0.0,0.0,10.0]]*) keywords 
should be used to define the box size. For non-periodic calculations this is necessary as the basis set and solver are based on the box dimensions.
However, if cell information is not provided, then by default a cell size will be automatically estimated (by ASH) based on the molecule size and the *qm_cell_shift_par* parameter 
will extend the box by an additional amount (6.0 Angstrom by default).

**Non-periodic calculations**

For non-periodic calculations, the CP2KTheory object should be defined with *periodic=False*, this is the default.
The cell should be specified as described above. Poisson solver should also be specified with the psolver keyword. 
The default is 'wavelet' which is probably the most efficient for non-periodic calculations. The accuracy of the solver can be controlled by the wavelet_scf_type keyword (see `CP2K-manual-wavelet <https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/POISSON/WAVELET.html>`_ ).
Another Poisson solver option is 'MT' (`CP2K-manual-MT <https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/POISSON/MT.html>`_ ).
In the case of the MT solver the cell should be at least 2 as large as the charge density (i.e. the molecule). The cell can be smaller for the wavelet solver.

**Periodic calculations**

For periodic calculations, the CP2KTheory object should be defined with *periodic=True*. The *periodic_type* is by default 'XYZ' (i.e. PBC in all directions).
The cell size should be specified as described above.
Poisson solver options are : 'PERIODIC', 'WAVELET', 'MULTIPOLE' or 'IMPLICIT'. The PERIODIC solver is recommended (only available for full 3D periodicity).



################################################################################
QM/MM
################################################################################

QM/MM calculations are possible in the ASH interface to CP2K. 
Unlike most other QM-codes, however, regular electrostatic embedding is not available for DFT-methods in CP2K so instead we use the 
GEEP (Gaussian Expansion of Electrostatic Potential) approach available in CP2K. This approach expands the MM pointcharges as Gaussians.
The GEEP approach is overall an improvement over traditional electrostatic embedding as it should prevent charge-leakage onto MM atoms (electron spill-out effect).
The GEEP approach, however, requires definition of radii on the MM-atoms which control the width of the Gaussians used to expand the MM pointcharges.

To use CP2K as QM-code in an ASH QM/MM calculation one needs be aware of a few things:

- The cell size must be specified (either *cell_dimension* or *cell_vectors*) but counterintuitively it needs to be specified for the whole system (QM+MM) and not just the QM-part as CP2K needs this information.
- Additionally the QM-cell size should be specified (where the electrons and basis sets are) and this should be a box encompassing the whole QM-region (slightly larger).
  The *qm_cell_dims* keyword can be used to specify this or alternatively ASH can also estimate the QM-cell size based on the QM-region size and the *qm_cell_shift_par* extension parameter (default 6).
  If the Poisson solver is wavelet, the QM-cell needs to be cubic (automatically done if the QM-cell size is estimated from the QM-region).
- A QM/MM job with CP2K in ASH can either be periodic or non-periodic. For non-periodic calculations it is recommended to use the wavelet Poisson solver.
- For periodic QM/MM calculations, one should typically set: *periodic=True* and *psolver='PERIODIC'*


One then should specify the QM/MM electrostatic coupling. For DFT only the Gaussian-based GEEP approach is available (*coupling='GAUSSIAN'*) while *coupling='COULOMB'* is available for semi-empirical systems.
GEEP can only be used with the wavelet or periodic Poisson solver (not 'MT')
The number of Gaussians used to expand each MM-center is controlled by the *GEEP_num_gauss keyword* (default=6). 
The width of the Gaussians depends on the defined MM-radius for each MM site which should vary according to the element. 
Element information of the MM-region is automatically passed onto CP2K and default MM-radii will be used:

.. code-block:: python

    #Element radii in Angstrom (will be converted to Bohrs by CP2K)
    element_radii_for_cp2k = {'H':0.44,'He':0.44,'Li':0.6,'Be':0.6,'B':0.78,'C':0.78,'N':0.78,'O':0.78,'F':0.78,'Ne':0.78,
                        'Na':1.58,'Mg':1.58,'Al':1.67,'Si':1.67,'P':1.67,'S':1.67,'Cl':1.67,'Ar':1.67,
                        'K':1.52,'Ca':1.6,'Sc':1.6,'Ti':1.6,'V':1.6,'Cr':1.6,'Mn':1.6,'Fe':1.6,'Co':1.6,
                        'Ni':1.6,'Cu':1.6,'Zn':1.6,'Br':1.6, 'Mo':1.7}

If the user wants to use different MM radii for each element, then this can be specified with the *mm_radii* keyword which should point to a dictionary containing radii
for each element present in the system. It is also possible to use the *MM_radius_scaling* keyword to scale the radii by a factor (default=1.0).

ASH handles the creation of linkatoms and charge-shifting at a QM-MM boundary and this information is provided to CP2K as a modified XYZ-file.
It is unclear whether the automatic dipole-correction (addition of charges to maintain dipole), commonly employed in charge-shifted electrostatic embedding QM/MM is useful when combined
with the GEEP procedure of CP2K. It is thus possible to turn it off with the *dipole_correction* keyword in the QMMMTheory object.


################################################################################
Examples
################################################################################

In the examples below it is assumed that the CP2K binaries are already in PATH (no need to use cp2kdir)

**Minimal non-periodic geometry optimization of MeOH in vacuum:**

.. code-block:: python

    from ash import *

    numcores=2
    frag = Fragment(xyzfile="MeOH.xyz",charge=0, mult=1)

    #Basis set and pseudopotential information per element
    basis_dict={'C':'DZVP-MOLOPT-SR-GTH','O':'DZVP-MOLOPT-SR-GTH','H':'DZVP-MOLOPT-SR-GTH'}
    potential_dict={'C':'GTH-PBE-q4','O':'GTH-PBE-q6','H':'GTH-PBE-q1'}

    #Minimal CP2KTheory definition: no periodicity, psolver=wavelet by default, basis_method='GAPW', cutoff=250,rel_cutoff=60
    #Cell dimensions are estimated from molecule size
    qm = CP2KTheory(cp2k_bin_name="cp2k.ssmp",basis_dict=basis_dict,potential_dict=potential_dict,functional='PBE',numcores=numcores)
    #Geometry optimization
    Optimizer(theory=qm, fragment=frag)

**Periodic geometry optimization of MeOH in vacuum:**

.. code-block:: python

    from ash import *

    numcores=2
    frag = Fragment(xyzfile="MeOH.xyz",charge=0, mult=1)

    #Basis set and pseudopotential information per element
    basis_dict={'C':'DZVP-MOLOPT-SR-GTH','O':'DZVP-MOLOPT-SR-GTH','H':'DZVP-MOLOPT-SR-GTH'}
    potential_dict={'C':'GTH-PBE-q4','O':'GTH-PBE-q6','H':'GTH-PBE-q1'}

    #Periodic CP2KTheory definition with specified cell dimensions
    qm = CP2KTheory(cp2k_bin_name="cp2k.ssmp",basis_dict=basis_dict,potential_dict=potential_dict,functional='PBE',numcores=numcores,
                    periodic=True,cell_dimensions=[10,10,10,90,90,90], psolver='periodic',basis_method='GPW', ngrids=4, cutoff=450, rel_cutoff=50)
    #Geometry optimization
    Optimizer(theory=qm, fragment=frag)


**Non-periodic QM/MM geometry optimization and frequencies of a protein**

Here using the simple solvated lysozyme protein as a test system with a threonine sidechain in the QM-region.

.. code-block:: python

    from ash import *

    numcores=4

    #Defining path to dir containing forcefield files and coordinates
    forcefielddir="./"
    psffile=forcefielddir+"step3_pbcsetup.psf"
    topfile=forcefielddir+"top_all36_prot.rtf"
    prmfile=forcefielddir+"par_all36_prot.prm"
    xyzfile=forcefielddir+"coordinates.xyz"

    #Read coordinates from XYZ-file
    frag = Fragment(xyzfile=xyzfile)

    #Creating OpenMM object from CHARMM-files
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=prmfile, periodic=True, charmm_periodic_cell_dimensions=[80.0, 80.0, 80.0, 90.0, 90.0, 90.0],
        autoconstraints=None, rigidwater=False)
    #CP2KTheory object
    basis_dict={'C':'DZVP-MOLOPT-SR-GTH','N':'DZVP-MOLOPT-SR-GTH','O':'DZVP-MOLOPT-SR-GTH','H':'DZVP-MOLOPT-SR-GTH'}
    potential_dict={'C':'GTH-PBE-q4','N':'GTH-PBE-q5', 'O':'GTH-PBE-q6','H':'GTH-PBE-q1'}
    functional='PBE'
    #cell_dimensions are for full system (slight expansion was necessary
    #QM-cell dimensions here defined manually
    qm = CP2KTheory(basis_dict=basis_dict,potential_dict=potential_dict,functional=functional, psolver='wavelet', coupling='GAUSS',
        periodic=False, cell_dimensions=[82.0, 82.0, 82.0, 90.0, 90.0, 90.0], qm_cell_dims=[12.0,12.0,12.0], numcores=numcores)

    #act and qmatoms lists. Defines QM-region (atoms described by QM) and Active-region (atoms allowed to move)
    #IMPORTANT: atom indices begin at 0.
    #Here selecting the side-chain of a threonine residue
    qmatoms = [569,570,571,572,573,574,575,576]
    actatoms = qmatoms

    # Create QM/MM OBJECT by combining QM and MM objects above. Dipole-correction turned off.
    qmmmobject = QMMMTheory(qm_theory=qm, mm_theory=openmmobject, printlevel=2, dipole_correction=False,
                            fragment=frag, embedding="Elstat", qmatoms=qmatoms, qm_charge=0, qm_mult=1)

    #Run geometry optimization using geomeTRIC optimizer and HDLC coordinates. Using QM-region as active-region.
    Optimizer(theory=qmmmobject, fragment=frag, ActiveRegion=True, actatoms=actatoms,
                        maxiter=200, coordsystem='hdlc', charge=0,mult=1)

    #Partial numerical Hessian calculation
    NumFreq(theory=qmmmobject, fragment=frag, hessatoms=actatoms)

**Periodic QM/MM MD simulation of a protein**