ORCA interface
======================================

ORCATheory class:

.. code-block:: python
    
  class ORCATheory:
      def __init__(self, orcadir=None, fragment=None, orcasimpleinput='', printlevel=2, extrabasisatoms=None, extrabasis=None, TDDFT=False, TDDFTroots=5, FollowRoot=1,
                  orcablocks='', extraline='', first_iteration_input=None, brokensym=None, HSmult=None, atomstoflip=None, numcores=1, nprocs=None, label=None, moreadfile=None, 
                  autostart=True, propertyblock=None, keep_each_run_output=False, print_population_analysis=False):



**ORCATheory** options:

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
   * - ``fragment``
     - ASH Fragment
     - None
     - ASH Fragment object.
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
     - Whether to do TDDFT or not. If part of a Gradient job or Optimization job then the excited state gradient is calculated and used.
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
     - Whether to do a Flipspin ORCA calculation to find a BS solution. Requires HSmult and atomstoflip options.
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
   * - ``autostart``
     - string
     - 'Plottyplot'
     - X
   * - ``numcores``
     - integer
     - 1
     - Number of cores to use for ORCA
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
Parallelization
################################################################################

ORCA parallelization is handled by OpenMPI. By specifying the numcores=X, a *%pal numcores X end block* will be added to the
ORCA inputfile created by Ash. ORCA will then call the OpenMPI mpirun binary when needed and this requires the
correct OpenMPI version to be available.
Make sure the recommended OpenMPI version for the ORCA version you are using is available. This typically requires
setting (in the shell or jobscript):

export PATH=/path/to/openmpi/bin:$PATH and export LD_LIBRARY_PATH=/path/to/openmpi/lib:$LD_LIBRARY_PATH

or alternatively loading the appropriate module. Set these variables in the job-script (see :doc:`basics`) that you are using.



################################################################################
Useful ORCA functions
################################################################################

In addition to the ORCATheory class, there are a number of built-in functions inside that can help grab specific information from an ORCA outputfile etc.
To use these functions, the module has to be loaded first: import interfaces.interface_ORCA.py

**run__orca_plot**

.. code-block:: python

  def run_orca_plot(filename, option, orcadir=None, gridvalue=40,densityfilename=None, mo_operator=0, mo_number=None):





To be documented:

grab_EFG_from_ORCA_output(filename)

ICE_WF_size(filename)

QRO_occ_energies_grab(filename)

CASSCF_natocc_grab(filename)

SCF_FODocc_grab(filename)

MP2_natocc_grab(filename)

check_stability_in_output(file)

grabEOMIPs(file)

grabatomcharges_ORCA(chargemodel,outputfile)

chargemodel_select(chargemodel)

read_ORCA_Hessian(hessfile)

ORCAfrequenciesgrab(hessfile)

write_ORCA_Hessfile(hessian, coords, elems, masses, hessatoms,outputname)

grab_spin_expect_values_ORCA(file)

MolecularOrbitalGrab(file)

grabtrajenergies(filename)

tddftgrab(file)

xesgrab(file)

grab_HF_and_corr_energies(file, DLPNO=False, F12=False)

scfenergygrab(file)

finalenergiesgrab(file)

checkORCAfinished(file)

grab_coordinates_from_ORCA_output(outfile)

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
