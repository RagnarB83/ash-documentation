ORCA interface
======================================

ORCATheory class:

.. code-block:: python
    
  class ORCATheory:
      def __init__(self, orcadir=None, orcasimpleinput='', printlevel=2, extrabasisatoms=None, extrabasis=None, 
                  TDDFT=False, TDDFTroots=5, FollowRoot=1,
                  orcablocks='', extraline='', first_iteration_input=None, 
                  brokensym=None, HSmult=None, atomstoflip=None, numcores=1, nprocs=None, label=None, 
                  moreadfile=None,moreadfile_always=False, autostart=True,
                  propertyblock=None, keep_each_run_output=False, print_population_analysis=False, 
                  filename="orca", check_for_errors=True, check_for_warnings=True,
                  fragment_indices=None):

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``orcadir``
     - string
     - None
     - Path to ORCA directory.
   * - ``orcasimpleinput``
     - string
     - ''
     - Definition of the ORCA simple-input line
   * - ``orcablocks``
     - string (multiline)
     - ''
     - Used for block-input in the ORCA inputfile.
   * - ``extraline``
     - string
     - ''
     - Additional inputfile-option for ORCA.
   * - ``printlevel``
     - integer
     - 2
     - How much output printed by the ORCA module.
   * - ``extrabasisatoms``
     - list
     - None
     - What atomindices should have a different basis set (gets added to coordinate block)
   * - ``extrabasis``
     - string
     - None
     - What the basis set on extrabasisatoms should be
   * - ``TDDFT``
     - Boolean
     - False
     - | Whether to do TDDFT or not. If part of a Gradient job or Optimization job then the excited state
       | gradient is calculated and used.
   * - ``TDDFTroots``
     - integer
     - 5
     - How many TDDFT roots to calculate if TDDFT=True
   * - ``FollowRoot``
     - integer
     - 1
     - What excited state root to calculate gradient for if TDDFT=True.
   * - ``brokensym``
     - Boolean
     - False
     - | Whether to do a Flipspin ORCA calculation to find a BS solution. 
       | Requires HSmult and atomstoflip options.
   * - ``HSmult``
     - integer
     - None
     - What high-spin multiplicity to use in a brokensym=True job.
   * - ``atomstoflip``
     - list
     - None
     - What atom indices to spin-flip.
   * - ``moreadfile``
     - string
     - None
     - Name of file or path to file of a GBWfile to read in to the ORCA calculation
   * - ``moreadfile_always``
     - Boolean
     - False
     - | Whether moreadfile option is constantly applied for all runs using this ORCATheory
       | object or only for first run. Default: False meaning moreadfile is only used for 
       | first run using ORCATHeory object.
   * - ``autostart``
     - Boolean
     - True
     - Whether ORCA will automatically try to read orbitals from a GBW file with same basename.
   * - ``numcores``
     - integer
     - 1
     - Number of cores to use for ORCA
   * - ``filename``
     - string
     - 'orca'
     - Name of inputfile and outputfile
   * - ``label``
     - string
     - None
     - Label for ORCA object. Useful if working with many ORCATheory objects to distinguish them.
   * - ``propertyblock``
     - string
     - None
     - String containing ORCA-block input (e.g. %eprnmr) that must come after the coordinates.
   * - ``keep_each_run_output``
     - Boolean
     - False
     - Whether to keep copy of each ORCA outputfile from each run-call (e.g. each Opt-step).
   * - ``print_population_analysis``
     - Boolean
     - False
     - Whether to print Mulliken population analysis for each step
   * - ``print_population_analysis``
     - Boolean
     - False
     - Whether to print Mulliken population analysis for each step
   * - ``check_for_errors``
     - Boolean
     - True
     - Whether to check for errors in ORCA output once ORCA calculation is done.
   * - ``check_for_warnings``
     - Boolean
     - True
     - Whether to check for warnings in ORCA output once ORCA calculation is done.
   * - ``fragment_indices``
     - list of lists
     - None
     - | Optional: list of lists of atom indices that specify whether atoms belong 
       | to a specific ORCA fragment (e.g. for ORCA multi-level PNO calculations). 
       | Example: [[1,2,3],[10,11,12],[13,14,15]]. Will affect the coordinate-block
       | in the ORCA inputfile. For QM/MM: atom indices must be in QM-region. 
   * - ``cpcm_radii``
     - list of floats
     - None
     - | By providing a list of radii (in Å) for each atom in the molecule, 
       | the CPCM radii will manually be changed in the ORCA inputfile. 
       | Typically used with DRACO-radii :doc:`helper_programs`

################################
Finding the ORCA program
################################

ASH can find the ORCA program in a few different ways.

- ASH will first check if the orcadir argument has been set which should be a string that points to the directory where the orca program is located, e.g. "orcadir=/path/to/orca_5_0_2". This option takes precedence.
- If the orcadir argument has not been provided ASH will next see if orcadir has been provided in the ASH settings (~/ash_user_settings.ini file): See :doc:`basics`
- If orcadir has also not been defined at all, ASH will next search the operating systems's PATH environment variable for an executable "orca" and if found, will set the orcadir accordingly and use that ORCA version.  This can be a convenient option if you make sure to define your shell environments carefully in your jobscript or shell-startup file. Be careful, however, if you have multiple versions of the program available.


.. warning:: The ORCA program binaries are nowadays often provided as a small-size shared version (has dynamically linked binaries). This means that for ORCA to run using the shared-library version, both the PATH and LD_LIBRARY_PATH needs to be set in the shell environment (should point to the ORCA directory).
  ASH can not set the LD_LIBRARY_PATH (must be done in the shell environment beforehand) and thus if LD_LIBRARY_PATH has not been set properly in the shell, ORCA will crash when called by ASH.
  This means that it is usually best to set the PATH and LD_LIBRARY_PATH to ORCA in your jobscript or login shell-file (.bashrc, .bash_profile etc.) and ASH will then be able to find ORCA like that.


################################################################################
Parallelization
################################################################################

ORCA parallelization is handled by OpenMPI. By specifying the numcores=X as a keyword when creating the ORCATheory object,
a *%pal numcores X end block* will be added to the ORCA inputfile that ASH creates. ORCA then handles its own parallelization, 
will call the OpenMPI mpirun binary when needed which does requires the correct OpenMPI version to be installed and available in PATH.
Make sure the recommended OpenMPI version for the ORCA version you are using is available. This typically requires
setting (in the shell or jobscript):

.. code-block:: text

  export PATH=/path/to/openmpi/bin:$PATH
  export LD_LIBRARY_PATH=/path/to/openmpi/lib:$LD_LIBRARY_PATH

or alternatively loading the appropriate module (if the computer is using modules). 
Set these variables in the job-script (see :doc:`basics`) that you are using.

################################################################################
Examples
################################################################################

The ORCA interface is quite flexible. orcasimpleinput and orcablocks keyword arguments (accepts single or multi-line strings) have to be provided and these keywords define what the ORCA-inputfile looks like. 
This means that you can completely control what type of electronic structure method should be used by ORCA including choosing aspects such as basis set, convergence and grid settings etc.
The geometry block will be added to the inputfile by ASH.
Note that ASH handles aspects such as telling ORCA what orbitals to read as well as parallelization.

.. warning:: Do not put parallelization information (! Pal4 or %pal nprocs 4 end)or job-type keywords such as "! Opt" "!Freq" to the orcasimpleinput and orcablocks variables. 
  Such functionality is handled by ASH separately.

.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz', charge=0, mult=1)
    #ORCA
    input="! BP86 def2-SVP tightscf"
    blocks="""
    %scf
    maxiter 200
    end
    %basis
    newgto F "ma-def2-SVP" end
    end
    """

    ORCAcalc = ORCATheory(orcasimpleinput=input, orcablocks=blocks, numcores=8)

    #Run a single-point energy job
    Singlepoint(theory=ORCAcalc, fragment=HF_frag)
    #An Energy+Gradient calculation
    Singlepoint(theory=ORCAcalc, fragment=HF_frag, Grad=True)


Here a fragment (here called HF_frag with a defined charge and multiplicity) is defined (from an XYZ file) and passed to the Singlepoint function along with an ORCAtheory object (called ORCAcalc). The input, and blocks string variables are defined and passed onto the ORCA object via keyword arguments. 
By default, the ORCA autostart feature is active, meaning that if an inputfile with name "orca-input.inp" is run, ORCA will
try to read orbitals from "orca-input.gbw" file if present. This is utilized automatically during geometry optimizations, numerical frequencies as well
as multiple single-point calculations sequentially. It is possible to turn this off by adding "!Noautostart" in the simple-inputline of the orcasimpleinput variable or by setting autostart=False when defining ORCATheory object.
It is also possible to have each ORCA-calculation read in orbitals from another source by using the: moreadfile keyword argument option:

.. code-block:: python

    ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput=input,
                        orcablocks=blocks, numcores=8, moreadfile="orbitals.gbw")


Note: For parallel-ASH calculations (ASH in parallel, ORCA in serial). The full path to the moreadfile may be required.


The ORCA object is then used by passing it to a function: e.g. Singlepoint, an optimizer, a QM/MM object, NumFreq function etc.
When the ORCA object is run (e.g. by the Singlepoint function, an optimizer etc.) it will create an ORCA inputfile
that will always be called orca-input.inp. This inputfile will look familiar to any ORCA user as it will contain a "Simpleinput line", Block-input
a coordinate block etc. (cordinates in Å). ASH will then tell ORCA to run the inputfile and an outputfile called orca-input.out will be created.
Once the ORCA calculation is done the outputfile (or other files) is read for information (usually the energy and gradient) by ASH
and ASH will continue. The ORCA inputfile , "orca-input.inp" may be replaced later (e.g. if an optimization job" and ORCA
will be run again.

################################################################################
Broken-symmetry DFT example
################################################################################

ORCA is quite a convenient program for finding broken-symmetry SCF solutions and within an ORCA inputfile one can 
easily tell ORCA to find a broken-symmetry solution within the %scf block (Flipspin or Brokensym options). 
While this could in principle simply be specified by the user in the orcablocks variable, this would have the drawback of ORCA
attempting a broken-symmetry search everytime the program is called, e.g. in every ASH optimization step or an ASH MD run. 
This is almost never what we want since we simply want to find the broken-symmetry SCF solution once and then reuse those orbitals in a subsequent step (and so on).
This is why ASH features a way to control the broken-symmetry search by the *brokensym* keyword in the ORCATheory object as shown below.

.. code-block:: python

    #Create fragment object from XYZ-file. Here a hypothetical Fe dimer complex
    frag=Fragment(xyzfile='fedimer.xyz', charge=0, mult=1)
    
    #ORCA settings
    inputline="! BP86 def2-SVP tightscf"
    #Here we specify a broken-symmetry solution with a high-spin multiplicity of 11 and flipping atoms no. 0
    ORCAcalc = ORCATheory(orcasimpleinput=inputline, brokensym=True, HSmult=11, atomstoflip=[0])

    #Run a broken-symmetry DFT geometry optimization
    Optimizer(theory=ORCAcalc, fragment=frag)

Running this job would have the effect of ASH initially writing an ORCA inputfile containing broken-symmetry settings (Flipspin and FinalMS keywords, high-spin multiplicity etc.)
but this would only apply to the first step of the geometry optimization. 
Once the ORCATheory object has been run once with broken-symmetry settings, the broken-symmetry feature is automatically turned off.
The next time the ORCATheory object is run (the next geometry optimization step of the above example), ASH creates an ORCA inputfile
with regular SCF inputsettings with the spin multiplicity being the low-spin BS multiplicity. 
Since the broken-symmetry SCF orbitals are available in the GBW file they are automatically loaded.



################################################################################
ORCA_External_Optimizer
################################################################################

It is possible to use ORCA as an external optimizer for ASH. 
This means that the ORCA geometry optimizer will be used with an ASH Theory level as input.
This functionality has not been tested much.

.. code-block:: python

  def ORCA_External_Optimizer(fragment=None, theory=None, orcadir=None, charge=None, mult=None):


################################################################################
Wrapper around ORCA helper programs.
################################################################################

ASH features wrappers around useful ORCA programs such as orca_plot, orca_mapspc and orca_2mkl.

**run_orca_plot**

.. code-block:: python

  # Simple Wrapper around orca_plot for creating Cube-files of MOs or densitities.
  def run_orca_plot(filename, option, orcadir=None, gridvalue=40,densityfilename=None, mo_operator=0, mo_number=None):

**run_orca_mapspc**
  
  .. code-block:: python
    
    # Simple Wrapper around orca_mapspc to create a broadened spectrum from a ORCA outputfile (creates .dat and .stk files)
    def run_orca_mapspc(filename, option, start=0.0, end=100, unit='eV', broadening=1.0, points=5000, orcadir=None):

**make_molden_file_ORCA**

  .. code-block:: python

    #Make a Molden file from ORCA GBW file (uses orca_2mkl)
    def make_molden_file_ORCA(GBWfile, orcadir=None):

################################################################################
ORCA fragment guess
################################################################################

It is possible to use the function **orca_frag_guess** to divide an ASH fragment into two fragments, run an ORCA calculation on each fragment
using an ORCATheory level and then combine the orbitals from the two fragments into a single GBW file (uses orca_mergefrag). 
This could be utilized to make a more accurate guess of the whole system.

  .. code-block:: python

    #Make an ORCA fragment guess. Returns name of GBW-file created ("orca_frag_guess.gbw")
    def orca_frag_guess(fragment=None, theory=None, A_indices=None, B_indices=None, A_charge=None, B_charge=None, A_mult=None, B_mult=None):


################################################################################
ORCA_JSON 
################################################################################

Since ORCA 6.0, ORCA-JSON feature has become more powerful, allowing for extracting MOs and all integrals from an ORCA GBW-file.
ASH features a few functions for conveniently creating or reading ORCA-JSON files.

.. code-block:: python

  #Wrapper around orca_2json to create JSON file from ORCA GBW file
  def create_ORCA_json_file(file, orcadir=None, format="json", basis_set=True, mo_coeffs=True, one_el_integrals=True,
                            two_el_integrals=False, two_el_integrals_type="ALL", dipole_integrals=False, full_int_transform=False):

  #Read ORCA json file: MO-coefficients, MO-energies, basis set, H,S,T matrices, 2-electron ints, densities etc.
  #Returns a dictionary with all information
  def read_ORCA_json_file(file):

  #Read ORCA MSPack (JSON-like binary format) file
  #Returns a dictionary with all information
  def read_ORCA_msgpack_file(file):

  #Read ORCA BSON (JSON-like binary format) file
  #Returns a dictionary with all information
  def read_ORCA_bson_file(file):

  #Get densities from data dictionary (from read_ORCA_json_file)
  def get_densities_from_ORCA_json(data):

  #Grab ORCA wfn from jsonfile or data-dictionary
  def grab_ORCA_wfn(data=None, jsonfile=None, density=None):

  #Reverse JSON to GBW
  def create_GBW_from_json_file(jsonfile, orcadir=None):

.. warning::  Do note that if the GBW-file contains a ROHF wavefunction then this will most likely not work due to the lack of ORCA-JSON handling for ROHF.

################################################################################
Creating FCIDUMP file from ORCA
################################################################################

The ORCA-JSON functionality can be utilized to create FCIDUMP files using the function **create_ORCA_FCIDUMP**.

.. code-block:: python

  def create_ORCA_FCIDUMP(gbwfile, header_format="FCIDUMP", filename="FCIDUMP_ORCA", orca_json_format="msgpack",
                          int_threshold=1e-16,  mult=1, full_int_transform=False,
                          convert_UHF_to_ROHF=True):

Examples:

.. code-block:: python

  # Create standard FCIDUMP file from ORCA GBW-file
  create_ORCA_FCIDUMP("orca.gbw", header_format="FCIDUMP", filename="FCIDUMP_ORCA",
                        int_threshold=1e-16, scf_type="RHF", mult=1)
  # Create MRCC-style FCIDUMP-file (fort.55) from ORCA GBW-file
  create_ORCA_FCIDUMP("orca.gbw", header_format="MRCC", int_threshold=1e-16, scf_type="RHF", mult=1)

.. warning::  Do note that if the GBW-file contains a ROHF wavefunction then this will most likely not work due to the lack of ORCA-JSON handling for ROHF.

.. warning:: If a UHF/UKS WF is found, then this is currently not handled. However, the convert_UHF_to_ROHF keyword can be set to True to make a naive conversion of UHF/UKS to ROHF.
  
################################################################################
Workflow to automate ORCA-orbital creation
################################################################################

ORCA is capable of producing various type of orbitals such as SCF-orbitals (RHF,UHF,ROHF etc.), MP2 natural orbitals, CC natural orbitals,
MRCI natural orbitals. The latter type of orbitals require a bit of know-how.
To automate the creation of these orbitals, ASH features a function called **ORCA_orbital_setup**.

.. code-block:: python

  #Function to prepare ORCA orbitals for another ORCA calculation
  def ORCA_orbital_setup(orbitals_option=None, fragment=None, basis=None, basisblock="", extrablock="", extrainput="", label="frag",
          MP2_density=None, MDCI_density=None, memory=10000, numcores=1, charge=None, mult=None, moreadfile=None,
          gtol=2.50e-04, nmin=1.98, nmax=0.02, CAS_nel=None, CAS_norb=None,CASCI=False, natorb_iterations=None,
          FOBO_excitation_options=None, MRCI_natorbiterations=0, MRCI_tsel=1e-6,
          ROHF=False, ROHF_case=None, MP2_nat_step=False, MREOMtype="MR-EOM",
          NMF=False, NMF_sigma=None):

Example on how to get CCSD natural orbitals from an unrelaxed CCSD density:

.. code-block:: python

  newmofile, nat_occupations = ORCA_orbital_setup(orbitals_option="CCSD", fragment=frag, label="CCSD" 
                basis="def2-SVP", MDCI_density="unrelaxed", charge=0, mult=1)
  # Returns name of the MO-file (here called CCSD_orca.mdci.nat)

################################################################################
Useful ORCA functions
################################################################################

In addition to the ORCATheory class, there are a number of built-in functions in ASH that are useful for ORCA functionlaity.
For example functions to grab specific information from an ORCA outputfile etc.
To use most these functions, the module has to be loaded first: 

.. code-block:: python

  from ash.interfaces.interface_ORCA.py import *


Functions for grabbing information from ORCA outputfiles:

.. code-block:: python

  #Simple function that grabs elements and coordinates from ORCA outputfile
  def grab_coordinates_from_ORCA_output(filename):

  #Grab Final single point energy. Ignoring possible encoding errors in file
  def ORCAfinalenergygrab(file, errors='ignore'):

  #Grab multiple Final single point energies in output. e.g. new_job calculation
  def finalenergiesgrab(file):

  #Grab SCF energy (non-dispersion corrected)
  def scfenergygrab(file):

  #Grab HF and correlation energies from ORCA output
  def grab_HF_and_corr_energies(file, DLPNO=False, F12=False):

  #Grab energies from unrelaxed scan in ORCA (paras block type)
  def grabtrajenergies(filename):

  #Grab ORCA timings. Return dictionary
  def ORCAtimingsgrab(file):

  #Grab gradient from ORCA engrad file
  def ORCAgradientgrab(engradfile):

  #Grab pointcharge gradient from ORCA pcgrad file
  def ORCApcgradientgrab(pcgradfile):

  #Grab XES state energies and intensities from ORCA output
  def xesgrab(file):

  #Grab TDDFT state energies from ORCA output
  def tddftgrab(file):

  #Grab TDDFT state intensities from ORCA output
  def tddftintens_grab(file):

  #Grab TDDFT orbital pairs from ORCA output
  def tddft_orbitalpairs_grab(file):

  #Grab molecular orbital energies from ORCA outputfile
  def MolecularOrbitalGrab(file):

  #Grab QRO energies from ORCA outputfile
  def QRO_occ_energies_grab(filename):

  #Grab <S**2> expectation values from outputfile
  def grab_spin_expect_values_ORCA(file):

  #Grab MP2 natural occupations from ORCA outputfile
  def MP2_natocc_grab(filename):

  #Grab SCF FOD occupations from ORCA outputfile
  def SCF_FODocc_grab(filename):

  #Grab CASSCF natural occupations from ORCA outputfile
  def CASSCF_natocc_grab(filename):

  #Find localized orbitals in ORCA outputfile for a given element. Returns orbital indices (to be fed into run_orca_plot)
  def orblocfind(outputfile, atomindex_strings=None, popthreshold=0.1):

  #Grab spin populations from ORCA outputfile
  def grabspinpop_ORCA(chargemodel,outputfile):

  #Grab atomic charges from ORCA outputfile
  def grabatomcharges_ORCA(chargemodel,outputfile):

  #Grab IPs from an EOM-IP calculation and also largest singles amplitudes.
  def grabEOMIPs(file):

  #Grab electric field gradients from ORCA outputfile
  def grab_EFG_from_ORCA_output(filename):

  #Grab ICE-WF info from CASSCF job
  def ICE_WF_size(filename):

  #Grab ICE-WF CFG info from CI job
  def ICE_WF_CFG_CI_size(filename):

  #Reading stability analysis from output. Returns true if stab-analysis good, otherwise falsee
  def check_stability_in_output(file):

Functions related to ORCA Hessian files:

.. code-block:: python

  #write ORCA-style Hessian file
  def write_ORCA_Hessfile(hessian, coords, elems, masses, hessatoms,outputname):

  #Function to grab Hessian from ORCA-Hessian file. Returns 2d Numpy array
  def Hessgrab(hessfile):

  #Grab coordinates from ORCA-Hessian file. Returns elements and coordinates.
  def grabcoordsfromhessfile(hessfile):

  #Function to grab masses and elements from an ORCA Hessian file
  def masselemgrab(hessfile):

  #Read ORCA Hessian-file and return Hessian, elems, coords and masses
  def read_ORCA_Hessian(hessfile):

  #Grab frequencies from ORCA-Hessian file
  def ORCAfrequenciesgrab(hessfile):


Functions for creating ORCA inputfiles:

.. code-block:: python

  #Create PC-embedded ORCA inputfile from elems,coords, input, charge, mult,pointcharges
  def create_orca_input_pc(name,elems,coords,orcasimpleinput,orcablockinput,charge,mult, Grad=False, extraline='',
                          HSmult=None, atomstoflip=None, Hessian=False, extrabasisatoms=None, extrabasis=None,
                          moreadfile=None, propertyblock=None, fragment_indices=None):

  #Create simple ORCA inputfile from elems,coords, input, charge, mult,pointcharges
  def create_orca_input_plain(name,elems,coords,orcasimpleinput,orcablockinput,charge,mult, Grad=False, Hessian=False, extraline='',
                              HSmult=None, atomstoflip=None, extrabasis=None, extrabasisatoms=None, moreadfile=None, propertyblock=None,
                              ghostatoms=None, dummyatoms=None,fragment_indices=None):

  # Create ORCA pointcharge file based on provided list of elems and coords (MM region elems and coords) and list of point charges of MM atoms
  def create_orca_pcfile(name,coords,listofcharges):

  # Chargemodel select. Creates ORCA-inputline with appropriate keywords
  def chargemodel_select(chargemodel):


Functions for other ORCA functionality:

.. code-block:: python

  #Print gradient in ORCA format to disk
  def print_gradient_in_ORCAformat(energy,gradient,basename):


################################################################################
Useful ORCA workflows
################################################################################

Examples of useful ways to automate various ORCA calculations.


**Plot ORCA-calculated spectra (using orca_mapspc) and normalize**

Uses ASH functions: **grab_coordinates_from_ORCA_output**, **run_orca_mapspc**, **read_datafile**, **write_datafile**

.. code-block:: python

  from ash import *
  import glob

  #Simple ASH script to plot XES spectra from multiple ORCA XES-job outputfiles and normalize w.r.t. to number of absorber elements
  absorber_element="Fe"

  #orca_mapspc settings
  orca_mapspc_option='XESQ'
  broadening=1.0
  numpoints=5000
  start_value=0
  end_value=8000
  unit='eV'

  #Loop over ORCA outputfiles and run orca_mapspc
  for outfile in glob.glob("*.out"):
      print("Outfile:", outfile)
      #Get number of absorber elements in molecule from outputfile
      elems,coords = grab_coordinates_from_ORCA_output(outfile)
      elementcount = elems.count(absorber_element)
      print(f"Number of {absorber_element} atoms in file:", elementcount)
      #Get XES .at and .stk files via orca_mapspc
      run_orca_mapspc(outfile, orca_mapspc_option, start=start_value, end=end_value, unit=unit, broadening=broadening, points=numpoints)
      #Read .dat file. Get x and y values as numpy arrays
      x, y = read_datafile(outfile+".xesq.dat")
      #Scale y-values
      scalingfactor=elementcount
      write_datafile(x,y/scalingfactor, filename=outfile+f"_SCALED_by_{scalingfactor}.xesq.dat")
      #Read .stk file
      x, y = read_datafile(outfile+".xesq.stk")
      #Scale y-values
      write_datafile(x,y/scalingfactor, filename=outfile+f"_SCALED_by_{scalingfactor}.xesq.stk")
  #
