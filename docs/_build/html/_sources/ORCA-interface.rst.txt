ORCA interface
======================================

ORCATheory class:

.. code-block:: python
    
    class ORCATheory:
        def __init__(self, orcadir=None, fragment=None, charge=None, mult=None, orcasimpleinput='', orcablocks='', 
        printlevel=2, extrabasisatoms=None, extrabasis=None, TDDFT=False, TDDFTroots=5, FollowRoot=1, extraline='', 
        brokensym=None, HSmult=None, atomstoflip=None, nprocs=1, label=None, moreadfile=None, autostart=True, propertyblock=None):

Options:
- orcadir: string. Path to ORCA
- fragment: ASH Fragment object.
- charge: integer. Charge of molecule.
- mult: integer. Spin multiplicity of molecule
- orcasimpleinput: string. Definition of the ORCA simple-input line
- orcablocks: string (can be multiline). Used for block-input in the ORCA inputfile.
- extraline: string. Additional inputfile-option for ORCA.
- printlevel: integer. how much output by the ORCA module.
- extrabasisatoms: list. What atomindices should have a different basis set (gets added to coordinate block)
- extrabasis: string. What the basis set on extrabasisatoms should be
- TDDFT: Boolean. Whether to do TDDFT or not. If part of a Gradient job or Optimization job then the excited state gradient is calculated and used.
- TDDFTroots: integer. How many TDDFT roots to calculate.
- FollowRoot: integer. What excited state root to calculate gradient for.
- brokensym: Boolean. Whether to do a Flipspin ORCA calculation to find a BS solution. Requires HSmult and atomstoflip options.
- HSmult: integer. What high-spin multiplicity to use in a Flipspin job.
- atomstoflip: list. What atom indices to spin-flip.
- moreadfile: string. Name of file or path to file of a GBWfile to read in to the ORCA calculation
- autostart: Boolean. Whether to turn Autostart on or off (default: True)
- nprocs: integer. Number of cores to use for ORCA
- label: string. Label for ORCA object. Useful if working with many.
- propertyblock: string. String containing ORCA-block input (e.g. %eprnmr) that comes after the coordinates.

The ORCA interface is quite flexible. It currently requires the path to the ORCA installation to be passed on as a keyword
argument when creating object (alternatively it can be set in the ~/ash_user_settings.ini file). orcasimpleinput and orcablocks keyword arguments (accepts single or multi-line strings) have to be provided and these keywords define what the ORCA-inputfile looks like. The geometry block would be added to the inputfile by ASH.
Functionality such as telling ORCA what orbitals to read and parallelization would be handled by ASH as well.


.. code-block:: python

    #Create fragment object from XYZ-file
    HF_frag=Fragment(xyzfile='hf.xyz')
    #ORCA
    orcadir='/opt/orca_4.2.1'
    input="! BP86 def2-SVP Grid5 Finalgrid6 tightscf"
    blocks="""
    %scf
    maxiter 200
    end
    %basis
    newgto F "ma-def2-SVP" end
    end
    """

    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput=input, orcablocks=blocks, nprocs=8)

    #Run a single-point energy job
    Singlepoint(theory=ORCAcalc, fragment=HF_frag)
    #An Energy+Gradient calculation
    Singlepoint(theory=ORCAcalc, fragment=HF_frag, Grad=True)



Here a fragment (here called HF_frag) is defined (from an XYZ file) and passed to the Singlepoint function along with an
ORCAtheory object (called ORCAcalc). The orcadir, input, and blocks string variables are defined and passed onto the ORCA object via keywords, as
are charge and spin multiplicity. By default, the ORCA autostart feature is active, meaning that if an inputfile with name "orca-input.inp" is run, ORCA will
try to read orbitals from "orca-input.gbw" file if present. This is utilized automatically during geometry optimizations, numerical frequencies as well
as multiple single-point calculations sequentially. It is possible to turn this off by adding "!Noautostart" in the simple-inputline of the orcasimpleinput variable or by setting autostart=False when defining ORCATheory object.
It is also possible to have each ORCA-calculation read in orbitals from another source by using the: moreadfile keyword argument option:

.. code-block:: python

    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput=input,
                        orcablocks=blocks, nprocs=8, moreadfile="orbitals.gbw")


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

ORCA parallelization is handled by OpenMPI. By specifying the nprocs=X, a *%pal nprocs X end block* will be added to the
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
To be documented:

grab_EFG_from_ORCA_output(filename)

ICE_WF_size(filename)

QRO_occ_energies_grab(filename)

CASSCF_natocc_grab(filename)

SCF_FODocc_grab(filename)

MP2_natocc_grab(filename)

check_stability_in_output(file)

grabEOMIPs(file)

run_orca_plot(orcadir, filename, option, gridvalue=40,densityfilename=None, mo_operator=0, mo_number=None)

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