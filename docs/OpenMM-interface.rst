======================================
OpenMM interface
======================================

`OpenMM <https://openmm.org>`_ is an open-source molecular mechanics library written in C++. ASH features a flexible interface to the Python API of the OpenMM library. 
OpenMM has been designed to run on both CPU and GPU codes with the GPU code being particulary fast.



######################################
The OpenMMTheory class 
######################################

OpenMM can run either on the CPU or on the GPU platform by specifying platform as either: 'CPU', 'OpenCL' or 'CUDA' depending on whether there is a GPU available on the computer. 
For platform='CPU' the numcores can additionally be specified in order to control the number of CPU cores. Alternatively the shell environment variable OPENMM_CPU_THREADS can be set 
to control the number of CPU cores that OpenMM will use. If neither numcores keyword is provided or OPENMM_CPU_THREADS variable set, OpenMM will use the number of physical cores present.

The OpenMMTheory class:

.. code-block:: python

    class OpenMMTheory:
        def __init__(self, printlevel=2, platform='CPU', numcores=None, topoforce=False, forcefield=None, topology=None,
                     CHARMMfiles=False, psffile=None, charmmtopfile=None, charmmprmfile=None,
                     GROMACSfiles=False, gromacstopfile=None, grofile=None, gromacstopdir=None,
                     Amberfiles=False, amberprmtopfile=None,
                     cluster_fragment=None, ASH_FF_file=None, PBCvectors=None,
                     xmlfiles=None, pdbfile=None, use_parmed=False,
                     xmlsystemfile=None,
                     do_energy_decomposition=False,
                     periodic=False, charmm_periodic_cell_dimensions=None, customnonbondedforce=False,
                     periodic_nonbonded_cutoff=12, dispersion_correction=True,
                     switching_function_distance=10,
                     ewalderrortolerance=5e-4, PMEparameters=None,
                     delete_QM1_MM1_bonded=False, applyconstraints_in_run=False,
                     constraints=None, restraints=None, frozen_atoms=None, dummy_system=False, fragment=None, 
                     autoconstraints='HBonds', hydrogenmass=1.5, rigidwater=True):


**OpenMMTheory** options:

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``printlevel``
     - integer
     - 2
     - The printlevel.
   * - ``platform``
     - string
     - 'CPU'
     - Run on CPU or GPU. Options: 'CPU', 'OpenCL', 'CUDA'
   * - ``numcores``
     - integer
     - None
     - The number of CPU cores to use for 'CPU' platform.
   * - ``topoforce``
     - Boolean
     - False
     - | Whether to read in topology and forcefield objects created by 
       | OpenMM manually. Requires setting topology and forcefield also.
   * - ``forcefield``
     - OpenMM forcefield object
     - None
     - Read in OpenMM forcefield object defined manually by OpenMM.
   * - ``topology``
     - OpenMM forcefield object
     - None
     - Read in OpenMM topology object defined by manually by OpenMM.
   * - ``CHARMMfiles``
     - Boolean
     - False
     - Read in CHARMM forcefiles or not.
   * - ``psffile``
     - string
     - None
     - Name of CHARMM protein structure file.
   * - ``charmmtopfile``
     - string
     - None
     - YYYYYYYYYYYY
   * - ``charmmprmfile``
     - string
     - None
     - Name of CHARMM parameter file.
   * - ``GROMACSfiles``
     - Boolean
     - False
     - Read in GROMACS forcefiles or not.
   * - ``gromacstopfile``
     - string
     - None
     - Name of the GROMACS topology file.
   * - ``grofile``
     - string
     - None
     - Name of Gromacs coordinate file (.gro extension)
   * - ``gromacstopdir``
     - string
     - None
     - Path to the GROMACS topology directory.
   * - ``Amberfiles``
     - Boolean
     - False
     - Read in Amber forcefiles or not.
   * - ``amberprmtopfile``
     - string
     - None
     - Name of the Amber PRMTOP file.
   * - ``cluster_fragment``
     - Fragment
     - None
     - ASH Fragment objected created by Molcrys.
   * - ``ASH_FF_file``
     - string
     - None
     - Name of ASH forcefield file.
   * - ``PBCvectors``
     - list
     - None
     - List of lists of floats (nm) or list of OpenMM Vec3 objects. Units are nm.
   * - ``xmlfiles``
     - list
     - None
     - List of XML-files to read forcefield from.
   * - ``pdbfile``
     - string
     - None
     - | Name of PDB-file. Used for reading topology for xmlfiles, xmlsystemfile 
       | and ASH_FF_file options.
   * - ``use_parmed``
     - Boolean
     - False
     - | Whether to use Parmed library to help read in CHARMM/Amber/GROMACS
       | files. Requires install of Parmed Python library.
   * - ``xmlsystemfile``
     - string
     - None
     - Name of XML system file to read in.
   * - ``do_energy_decomposition``
     - Boolean
     - False
     - | Do energy decomposition of each energy evaluation 
       | (when called by Singlepoint or optimizer).
   * - ``periodic``
     - Boolean
     - False
     - Periodic boundary conditions or not.
   * - ``charmm_periodic_cell_dimensions``
     - None
     - None
     - | Periodic cell dimension for CHARMM-setup of system. Example: 
       | charmm_periodic_cell_dimensions= [200, 200, 200, 90, 90, 90]
   * - ``customnonbondedforce``
     - Boolean
     - False
     - | Expert option: whether CustomNonbondedForce is used instead 
       | of NonbondedForce.
   * - ``periodic_nonbonded_cutoff``
     - float
     - 12.0
     - Cutoff for the periodic nonbonding interaction.
   * - ``dispersion_correction``
     - Boolean
     - True
     - Dispersion correction in periodic nonbonding evaluation.
   * - ``switching_function_distance``
     - float
     - 10.0
     - Switching function distance in Å units.
   * - ``ewalderrortolerance``
     - float
     - 5e-4
     - Error tolerance for the periodic electrostatics Ewald algorithm.
   * - ``PMEparameters``
     - list
     - None
     - | Optional manual parameters for the Particle Mess Ewald algorithm. 
       | Alternative to ewalderrortolerance keyword.
   * - ``delete_QM1_MM1_bonded``
     - Boolean
     - False
     - For QM/MM job, whether QM1-MM1 are deleted or not.
   * - ``applyconstraints_in_run``
     - Boolean
     - False
     - Exper option: Whether constraints are applied in run method. Should be False.
   * - ``constraints``
     - list of lists
     - None
     - | List of lists of constraint definitions based on atom indices. Either 
       | [[atom_i,atom_j]] or [[atom_i,atom_j, d]], e.g. [[700,701],[703,704]]
       | or [[700,701, 1.05],[702,703, 1.14]], where d: distance (Å))
   * - ``restraints``
     - list of lists
     - None
     - | List of lists of restraint definitions ([[atom_i,atom_j, d, k ]], e.g.
       | [[700,701, 1.05, 5.0 ]], d: distance (Å) k: force constant (kcal/mol*Å^-2))
   * - ``frozen_atoms``
     - list
     - None
     - List of atom indices to keep frozen during MD (particle mass set to 0).
   * - ``dummy_system``
     - Boolean
     - False
     - | If True, OpenMM will set up a dummy MM system based on provided fragment
       | (see below). Used for QM dynamics option in OpenMM_MD.
   * - ``fragment``
     - ASH Fragment
     - None
     - ASH fragment to provide when dummy_system is True.
   * - ``autoconstraints``
     - string
     - 'HBonds'
     - | Type of automatic constraints to apply to system. Options: 'HBonds' 
       | (constrain all X-H bonds), 'AllBonds' (constrain all bonds), 'HAngles'
       | (constrain all bonds and  H-X-H and H-O-X angles).
   * - ``hydrogenmass``
     - float
     - 1.5
     - | Hydrogen mass repartioning value. 1.5 is OpenMM and ASH default. 
       | Improves numerical stability.
   * - ``rigidwater``
     - Boolean
     - True
     - | Whether to automatically apply rigid water constraints for recognized 
       | water models (e.g. TIP3P) found in system. Note: needs to be turned off for 
       | Singlepoint/Optimizations.







It is possible to read in multiple types of forcefield files: AmberFiles, CHARMMFiles, GROMACSFiles or an OpenMM XML forcefieldfile.
Note: In rare cases OpenMM fails to read in Amber/CHARMM/GROMACS files correctly. In those cases the Parmed library may be more successful (use_parmed=True). Requires ParMed (pip install parmed).

*Example creation of an OpenMMtheory object with CHARMM-files:*

.. code-block:: python

    forcefielddir="/path/to/dir"
    topfile=forcefielddir+"top_all36_prot.rtf"
    parfile=forcefielddir+"par_all36_prot.prm"
    psffile=forcefielddir+"new-XPLOR-psffile.psf"
    openmmobject = OpenMMTheory(CHARMMfiles=True, psffile=psffile, charmmtopfile=topfile,
                               charmmprmfile=parfile)

*Example creation of an OpenMMtheory object with GROMACS-files:*

.. code-block:: python

    openmmobject = OpenMMTheory(GROMACSfiles=True, gromacstopdir="/path/to/gromacstopdir",
                    gromacstopfile="gromacstopfile.top", grofile="grofile.gro")

*Example creation of an OpenMMtheory object with AMBER files:*

.. code-block:: python

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile="/path/to/amberprmtopfile")

*Example creation of an OpenMMtheory object with OpenMM XML file:*

.. code-block:: python

    openmmobject = OpenMMTheory(xmlfiles=["example.xml"]) #File example.xml should be in dir
    #or
    openmmobject = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "specialresidue.xml"]) 
    #Here the "charmm36.xml" and "charmm36/water.xml" files are found in the OpenMM library.



Any Openmmtheory object can used to create a QM/MM theory object. See :doc:`module_QM-MM` page.

**Periodic boundary conditions:**

- If periodic boundary conditions are chosen (periodic=True) then the PBC box parameters are automatically found in the Amber PRMTOP file or the GROMACS Grofile or in the case of CHARMM-files they need to be provided: charmm_periodic_cell_dimensions
- The Ewald error tolerance (ewalderrortolerance) can be modified (default: 5e-4)
- PME parameters can be modified: PMEparameters=[alpha_separation,numgridpoints_X,numgridpoints_Y,numgridpoints_Z] 
- The periodic nonbonded cutoff can be modified. Default: 12 Å
- Long-range dispersion correction can be turned on or off. Default: True
- The switching function distance can be changed. Default: 10 Å. Used for CHARMM and XML files.
- The box dimensions can also be modified by PBCvectors= keyword argument:
    Example: PBCvectors=[[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]] with values in Å.

######################################
Molecular Dynamics via OpenMM
######################################

It is possible to run molecular dynamics of an ASH system using the MD algorithms present in the OpenMM library.
The OpenMM_MD function takes as argument an ASH fragment, a theory object and the user then selects an integrator of choice, simulation temperature, simulation length, timestep, optional additional thermostat, barostat etc.
The theory level can be OpenMMTheory, QMMMTheory or even a simple QMTheory.
Some options are only available for OpenMMTheory.

See `OpenMM documentation page <http://docs.openmm.org/latest/userguide/application.html#integrators>`_  for details about the integrators, thermostats, barostats etc.

- Available Integrators: Langevin, LangevinMiddleIntegrator, NoseHooverIntegrator, VerletIntegrator, VariableLangevinIntegrator, VariableVerletIntegrator
- Available Barostat: MonteCarloBarostat
- Optional additional thermostat: Anderson

.. code-block:: python

    def OpenMM_MD(fragment=None, theory=None, timestep=0.004, simulation_steps=None, simulation_time=None,
                traj_frequency=1000, temperature=300, integrator='LangevinMiddleIntegrator',
                barostat=None, pressure=1, trajectory_file_option='DCD', trajfilename='trajectory',
                coupling_frequency=1, charge=None, mult=None,
                anderson_thermostat=False,
                enforcePeriodicBox=True, dummyatomrestraint=False, center_on_atoms=None, solute_indices=None,
                datafilename=None, dummy_MM=False, plumed_object=None, add_center_force=False,
                center_force_atoms=None, centerforce_constant=1.0, barostat_frequency=25, specialbox=False):


**OpenMM_MD** options:

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``fragment``
     - ASH Fragment
     - None
     - The ASH fragment.
   * - ``theory``
     - ASH Theory
     - None
     - The ASH Theory object.
   * - ``timestep``
     - float
     - 0.004
     - | The timestep . Default: 0.004 ps (suitable for LangevinMiddleIntegrator 
       | dynamics with frozen X-H bonds)
   * - ``simulation_steps``
     - integer
     - None
     - Number of simulation steps to take. Alternative to simulation_time below
   * - ``simulation_time``
     - float
     - None
     - Length of simulation in picoseconds. Alternative to simulation_time above.
   * - ``temperature``
     - integer
     - 300
     - The temperature in Kelvin.
   * - ``integrator``
     - string
     - LangevinMiddleIntegrator
     - | The integrator to use. Options: 'Langevin', 'LangevinMiddleIntegrator', 
       | 'NoseHooverIntegrator', 'VerletIntegrator', 'VariableLangevinIntegrator',
       | 'VariableVerletIntegrator'
   * - ``coupling_frequency``
     - integer
     - 1
     - The coupling frequency of thermostat (in ps^-1 for Nosé-Hoover and Langevin-type)
   * - ``barostat``
     - string
     - None
     - Barostat to use for NPT simulations. Options: 'MonteCarloBarostat'
   * - ``barostat_frequency``
     - int
     - 25
     - Frequency of barostat update.
   * - ``pressure``
     - int
     - 1
     - Pressure to enforce by barostat
   * - ``anderson_thermostat``
     - Boolean
     - False
     - Whether to use Andersen Thermostat
   * - ``trajectory_file_option``
     - string
     - 'DCD'
     - | Type of trajectory file. Options: 'DCD' (compressed), 'PDB', 'NetCDFReporter' 
       | (compressed), 'HDF5Reporter' (compressed). Applies only to pure MM simulations.
   * - ``traj_frequency``
     - integer
     - 1000
     - Frequency of writing trajectory (every Xth timestep).
   * - ``trajfilename``
     - string
     - None
     - Name of trajectory file (without suffix).
   * - ``enforcePeriodicBox``
     - Boolean
     - True
     - Enforce PBC image during simulation. Fixes PBC-image artifacts in trajectory.
   * - ``center_on_atoms``
     - list
     - None
     - Expert options: Center system on these atoms.
   * - ``dummyatomrestraint``
     - Boolean
     - False
     - Expert options: Dummy atom restraints.
   * - ``solute_indices``
     - list
     - None
     - Expert options: solute_indices
   * - ``add_center_force``
     - Boolean
     - False
     - Whether to add a spherical force that pushes atoms to the center.
   * - ``center_force_atoms``
     - list
     - None
     - List of atom indices that the center force acts on.
   * - ``centerforce_constant``
     - float
     - None
     - Value of the spherical center force in kcal/mol/Ang^2.
   * - ``specialbox``
     - Boolean
     - False
     - Expert option: Special box for QM/MM.
   * - ``plumed_object``
     - ASH-Plumed object
     - None
     - Expert option: Plumed object for biased dynamics.



Note that constraints, autoconstraints, restraints and frozen_atoms must be defined in the OpenMMTHeory object before.



Example:

.. code-block:: python

    from ash import *

    #Forcefield parameters
    forcefielddir="/home/bjornsson/ASH-DEV_GIT/testsuite/OpenMM-files-for-tests/dhfr/charmm/"
    psffile=forcefielddir+"step3_pbcsetup.psf"
    topfile=forcefielddir+"top_all36_prot.rtf"
    prmfile=forcefielddir+"par_all36_prot.prm"

    #Defining fragment
    xyzfile=forcefielddir+"file.xyz"
    frag = Fragment(xyzfile=xyzfile, conncalc=False)

    #Defining OpenMM theory object: CHARMM forcefield with periodic boundary conditions
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=prmfile, periodic=True, charmm_periodic_cell_dimensions=[80, 80, 80, 90, 90, 90],
        dispersion_correction=False, periodic_nonbonded_cutoff=12, switching_function_distance=10,
        PMEparameters=[1.0/0.34, 90, 90, 90])

    #Launching a molecular dynamics simulation
    OpenMM_MD(fragment=frag, theory=openmmobject, timestep=0.001, simulation_steps=20, traj_frequency=1, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')


**General constraints or H-mass modification:**

- In order to allow shorter timesteps in MD simulations it is common to utilize some general constraints in biomolecular simulations, e.g. all X-H bonds, all bonds or even all-bond and some angles. This can be accomplished  via the autoconstraints option (NOTE: an option to OpenMMTheory). autoconstraints can be set to: 'HBonds' (X-H bonds constrained), 'AllBonds' (all bonds constrained), 'HAngles' (all bonds and H-X-H and H-O-X angles constrained) or None (default)
- An alternative (or addition) is to change the masses of the hydrogen atoms (fastest-moving atoms). This is also an option to OpenMMTheory. hydrogenmass keyword takes an integer and can e.g. be 2 (mass of deuterium) or heavier. hydrogenmass=1.5 is default (default in OpenMM) .



General X-H constraints and deuterium-mass example:

.. code-block:: python

    from ash import *

    #Forcefield parameters
    forcefielddir="/home/bjornsson/ASH-DEV_GIT/testsuite/OpenMM-files-for-tests/dhfr/charmm/"
    psffile=forcefielddir+"step3_pbcsetup.psf"
    topfile=forcefielddir+"top_all36_prot.rtf"
    prmfile=forcefielddir+"par_all36_prot.prm"

    #Defining fragment
    xyzfile=forcefielddir+"file.xyz"
    frag = Fragment(xyzfile=xyzfile, conncalc=False)

    #Defining OpenMM theory object: CHARMM forcefield with periodic boundary conditions
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=prmfile, periodic=True, charmm_periodic_cell_dimensions=[80, 80, 80, 90, 90, 90], autoconstraints='HBonds', hydrogenmass=2)

    #Launching a molecular dynamics simulation
    OpenMM_MD(fragment=frag, theory=openmmobject, timestep=0.001, simulation_steps=20, traj_frequency=1, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')



Dealing with PBC image problems in trajectory. See `OpenMM FAQ <https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#how-do-periodic-boundary-conditions-work>`_
To obtain a more pleasing visualization of the trajectory you can "reimage" the trajectory afterwards using the program mdtraj (requires installation of mdtraj: pip install mdtraj)

Example:

.. code-block:: python

    from ash import *
    #Provide trajectory file, PDB topology file and final format of trajectory
    MDtraj_imagetraj("output_traj.dcd", "final_MDfrag_laststep.pdb", format='DCD')
    
    #If periodic box info is missing from trajectory file (can happen with CHARMM files):
    MDtraj_imagetraj("out", pdbtopology, format='DCD', unitcell_lengths=[100.0,100.0,100.0], unitcell_angles=[90.0,90.0,90.0])


######################################
PBC box relaxation via NPT 
######################################

This function allows one to conveniently run multiple NPT simulations (constant pressure and temperature) in order to converge the periodic box dimensions
of the system.


.. code-block:: python

    def OpenMM_box_relaxation(fragment=None, theory=None, datafilename="nptsim.csv", numsteps_per_NPT=10000,
                              volume_threshold=1.0, density_threshold=0.001, temperature=300, timestep=0.004,
                              traj_frequency=100, trajfilename='relaxbox_NPT', trajectory_file_option='DCD', coupling_frequency=1):
        """NPT simulations until volume and density stops changing

        Args:
            fragment ([type], optional): [description]. Defaults to None.
            theory ([type], optional): [description]. Defaults to None.
            datafilename (str, optional): [description]. Defaults to "nptsim.csv".
            numsteps_per_NPT (int, optional): [description]. Defaults to 10000.
            volume_threshold (float, optional): [description]. Defaults to 1.0.
            density_threshold (float, optional): [description]. Defaults to 0.001.
            temperature (int, optional): [description]. Defaults to 300.
            timestep (float, optional): [description]. Defaults to 0.004.
            traj_frequency (int, optional): [description]. Defaults to 100.
            trajectory_file_option (str, optional): [description]. Defaults to 'DCD'.
            coupling_frequency (int, optional): [description]. Defaults to 1.
        """


######################################
Simple minimization via OpenMM
######################################

A classical system setup typically requires a minimization to get rid of large initial forces related to non-ideal atom positions.
The simple minimizer in the OpenMM library works well for this purpose although achieving convergence can be difficult.
Typically a few 100-1000 steps of minimization is sufficient to get rid of the major forces.

Example:

.. code-block:: python

    from ash import *

    #Forcefield parameters
    forcefielddir="/home/bjornsson/ASH-DEV_GIT/testsuite/OpenMM-files-for-tests/dhfr/charmm/"
    psffile=forcefielddir+"step3_pbcsetup.psf"
    topfile=forcefielddir+"top_all36_prot.rtf"
    prmfile=forcefielddir+"par_all36_prot.prm"

    #Defining fragment
    xyzfile=forcefielddir+"file.xyz"
    frag = Fragment(xyzfile=xyzfile, conncalc=False)

    #Defining OpenMM theory object: CHARMM forcefield with periodic boundary conditions
    openmmobject = OpenMMTheory(psffile=psffile, CHARMMfiles=True, charmmtopfile=topfile,
        charmmprmfile=prmfile, periodic=True, charmm_periodic_cell_dimensions=[80, 80, 80, 90, 90, 90])

    #Launching a minimization
    OpenMM_Opt(fragment=frag, theory=openmmobject, maxiter=1000, tolerance=1)
    #After minimization, the ASH fragment is updated, a PDB-file is written out: frag-minimized.pdb
    #Alternative XYZ write-out:
    frag.write_xyzfile(xyzfilename="frag_afteropt.xyz")


If you want to do a simple minimization of only the H-atoms of your system (e.g. your protein with newly added H-atoms),
you can do this by freezing all non-H atoms. An ASH fragment can conveniently give you lists of atom indices by the built-in functions:

- fragment.get_atomindices_for_element('C')   #List of atom-indices for carbon atoms in the system
- fragment.get_atomindices_except_element('H')   #List of atom-indices for all atoms except the chosen element (here H).

Note: all constraints in the OpenMM object needs to be turned off for (autoconstraints=None, rigidwater=False) for this many frozen atoms (frozen atoms can not have constraints).

.. code-block:: python

    from ash import *


    pdbfile = "ash_inp.pdb"
    prmtopfile = "prmtop"

    frag = Fragment(pdbfile=pdbfile)

    #List of all non-H atoms
    allnonHatoms=frag.get_atomindices_except_element('H')

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile=prmtopfile, periodic=True,
            platform='CPU', autoconstraints=None, rigidwater=False, frozen_atoms=allnonHatoms)

    OpenMM_Opt(fragment=frag, theory=openmmobject, maxiter=1000, tolerance=1)



######################################
System setup via OpenMM: Modeller
######################################

OpenMM features a convenient a `PDBfixer program <https://github.com/openmm/pdbfixer>`_ and  `Modeller tool <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.modeller.Modeller.html>`_
that together are capable of setting up a new biomolecular system from scratch. See also `OpenMM-Model-building and editing <http://docs.openmm.org/7.2.0/userguide/application.html#model-building-and-editing>`_
As ASH features a highly convenient interface to these programs and OpenMM itself this allows near-automatic system-setup for biomolecular systems.

.. code-block:: python

    def OpenMM_Modeller(pdbfile=None, forcefield=None, xmlfile=None, waterxmlfile=None, watermodel=None, pH=7.0,
                        solvent_padding=10.0, solvent_boxdims=None, extraxmlfile=None, residue_variants=None,
                        ionicstrength=0.1, pos_iontype='Na+', neg_iontype='Cl-', platform='CPU'):


The OpenMM_Modeller function returns an ASH OpenMMTheory object and ASH fragment object that can be used directly as theory level for future calculations.
OpenMM_Modeller will also print various PDB-files associated with each step of the setup (H-addition, solvation, ionization etc.) that can be visualized for correctness.
An XML file associated with the system is created that can be used to create future OpenMMtheory objects from.

Lysozyme example (simple, no modifications required):

.. code-block:: python

    from ash import *

    #Original raw PDB-file (no hydrogens, nosolvent)
    #Download from https://www.rcsb.org/structure/1AKI
    pdbfile="1aki.pdb"


    #Defining residues with special user-wanted protonation states for residues in each indicated chain
    #Dictionary of dictionaries with the chainname (e.g. 'A','B') acting as keys for the outer dictionary 
    #and the resids being keys for the inner dictionary
    #Example: residue_variants={'A':{0:'LYN', 17:'CYX', 18:'ASH', 19:'HIE', 20:'HID', 21:'GLH' }}
    #resid 1: neutral LYS, resid 17, deprotonated CYS, resid 18 protonated ASP, 
    #resid 19 epsilon-protonated HIS, resid 20 delta-protonated HIS, 21 protonated GLU.
    residue_variants={}

    #Setting up new system, adding hydrogens, solvent, ions and defining forcefield, topology
    openmmobject, ashfragment = OpenMM_Modeller(pdbfile=pdbfile, forcefield='CHARMM36', watermodel="tip3p", pH=7.0, 
        solvent_padding=10.0, ionicstrength=0.1, residue_variants=residue_variants)

    #MM minimization for 100 steps
    OpenMM_Opt(fragment=ashfragment, theory=openmmobject, maxiter=100, tolerance=1)

    #Classical MD simulation for 10 ps
    OpenMM_MD(fragment=ashfragment, theory=openmmobject, timestep=0.001, simulation_time=10, traj_frequency=100, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')


If the protein contains nonstandard residues (e.g. metallocofactors) that are not present in a typical protein forcefield (OpenMM_Modeller will exit with errors),
then these need to be provided using the extraxmlfile option.

.. code-block:: python

    openmmobject, ashfragment = OpenMM_Modeller(pdbfile=pdbfile, forcefield='CHARMM36', watermodel="tip3p", pH=7.0, 
        solvent_padding=10.0, ionicstrength=0.1, residue_variants=residue_variants, extraxmlfile="cofactor.xml")


The cofactor.xml file needs to define a forcefield (a nonbonded one at least) for the residue. 
Here defining a simple Fe(III) ion:

.. code-block:: text

    <ForceField>
    <AtomTypes>
    <Type name="FEX" class="Fe" element="Fe" mass="55.84700"/>
    </AtomTypes>
    <Residues>
    <Residue name="FE">
    <Atom name="FE" type="FEX"/>
    </Residue>
    </Residues>
    <NonbondedForce coulomb14scale="1.0" lj14scale="1.0">
    <Atom type="FEX" charge="3.0" sigma="1.3" epsilon="0.0"/>
    </NonbondedForce>
    <LennardJonesForce lj14scale="1.0">
    <Atom type="FEX" sigma="0.3" epsilon="0.00000"/>
    </LennardJonesForce>
    </ForceField>


See e.g. `Molecular Mechanics Tools <https://education.molssi.org/mm-tools/01-introduction/index.html>`_ for information on the format of the XML file.
See also information on the **write_nonbonded_FF_for_ligand** function on this page.

See :doc:`OpenMM-interface` for details and the :doc:`Metalloprotein-I` and :doc:`Metalloprotein-II` for step-by-step tutorials on the rubredoxin and ferredoxin metalloproteins.

Common error messages encountered when using OpenMM_Modeller on PDB-files:

-**ValueError: No template found for residue X (YYY).  This might mean your input topology is missing some atoms or bonds, or possibly that you are using the wrong force field.**

*This means that the parser encountered a completely unknown residue. You might have forgotten to read in the XML file to OpenMM_Modeller or the resname is not the same in the PDBfile as in the XML file. The atomnames and residue name in PDB-file must match the atomnames and residue name in the XML file. Also, element information (column 77-78) must be present in the PDB-file. 
It is also possible that PDB-file does not contain a valid N- or C-terminus for each peptide chain. Note that for a C-terminus, a terminal oxygen atom with atomname OXT is required.*

- **ValueError: Found multiple definitions for atom type: X**  :  

*This means that the atomtypes you defined inside <AtomTypes> </AtomTypes> are not unique. Either you have accidentally two lines with the same name for atomtype or the general forcefield (e.g. CHARMM) already contains an atomtype definition with this name. In that case, choose a unique name.*

- **KeyError: 'SXM'**  :  

*Possibly means that your atomname definition points to an atomtype-name that does not exist*


- **ValueError: No template found for residue X (YYY).  The set of atoms matches YYY, but the bonds are different.  Perhaps the chain is missing a terminal group?'**  :  

*This means there is some mismatch between the information present in the PDB-file and the information in the XML-file you provided.
It's possible that the PDB-file contains connectivity statements at the bottom of the PDB-file (CONE lines) but no bond information is present in the XML file.
Solution: Either add the missing bond to the residue definition so that it matches the CONE lines or simply delete the CONE information that you don't need.*



Valid alternative residue names for alternative protonation states of titratable residues:

- LYN instead of LYS: deprotonated lysine residue (NH2 instead of NH3)
- CYX instead of CYS: deprotonated cysteine residue (S- instead of SH)
- ASH instead of ASP: protonated aspartate residue (COOH instead of COO-)
- GLH instead of GLU: protonated glutamate residue (COOH instead of COO-)
- HID instead of HIS: histidine protonated at delta nitrogen
- HIE instead of HIS: histidine protonated at epsilon nitrogen

.. note:: Note: these names should not be used in the PDB-file. Only in the residue_variants dictionary that you provide to OpenMM_Modeller.



######################################
Small molecule solvation
######################################

**WORK IN PROGRESS**

ASH also features a function to solvate a small molecule automatically. This also makes use of the Modeller functionality of OpenMM but is intended to be used for molecules 
for where forcefield parameters are typically not available: e.g. metal complexes. Instead of regular forcefield parameters, nonbonded parameters (charges and Lennard-Jones parameters) 
are defined for the solute (used for classical and QM/MM simulations) which can be used to perfrom classical MM dynamics or QM/MM dynamics.

See also :doc:`Explicit-solvation` workflow page.


.. code-block:: python

    def solvate_small_molecule(fragment=None, charge=None, mult=None, watermodel=None, solvent_boxdims=[70.0,70.0,70.0], 
                               nonbonded_pars="CM5_UFF", orcatheory=None, numcores=1):

The solvate_small_molecule function reads in an ASH fragment, as well as charge and multiplicity, name of watermodel (e.g. "TIP3P"), size of solvent box, option for 
how the nonbonded parameters should be prepared, an optional ORCATheory object and optional numcores.

Options:

- watermodel (string). Can be: 'TIP3P' only for now
- solvent_boxdims (list of floats). Cubic box dimensions in Angstrom.
- nonbonded_pars (string). Options: 'CM5_UFF', 'DDEC3', 'DDEC6' or 'xtb_UFF'
- orcatheory (ORCATheory object). Optional ORCAtheory object defining the theory for deriving charges/LJ parameters
- numcores (integer). Number of cores used in ORCA/xTB calculations

nonbonded_pars options:

- 'CM5_UFF' derives CM5 charges (scaled Hirshfeld charges) from an ORCA calculation of the molecule and uses UFF Lennard-Jones parameters
- 'DDEC3' and 'DDEC6' derive both charges and LJ parameters from an ORCA calculation. Uses the Chargemol program.
- 'xtb_UFF' performs an xTB calculation to derive charges and uses UFF for LJ.



Example:

.. code-block:: python

    from ash import *

    numcores=4
    #Molecule definition
    mol=Fragment(xyzfile="3fgaba.xyz")
    mol.charge=0;mol.mult=1

    #Solvate molecule (70x70x70 Å TIP3P box)
    forcefield, topology, ashfragment = solvate_small_molecule(fragment=mol, charge=mol.charge, 
        mult=mol.mult, watermodel='tip3p', solvent_boxdims=[70,70,70], nonbonded_pars="CM5_UFF", 
        numcores=numcores)


The output of the solvate_small_molecule function are files: "system_aftersolvent.pdb", "newfragment.ygg", "newfragment.xyz" that can be used to inspect the coordinates of the system.

Additionally the function returns an OpenMM forcefield object, an OpenMM topology and an ASH fragment. These can be used in a next step to create an OpenMMTheory object:

.. code-block:: python

    from ash import *

    #Creating new OpenMM object from forcefield, topology and and fragment
    openmmobject =OpenMMTheory(numcores=numcores, Modeller=True, forcefield=forcefield, topology=topology, 
                    periodic=True, autoconstraints='HBonds', rigidwater=True)


The OpenMMTheory object can then be used on its own or can be combined with a QM theory to define a QM/MM theory object etc.
See :doc:`Explicit-solvation` workflow for more information on how to use solvate_small_molecule in a multi-step workflow.


##############################################
Create nonbonded forcefield file for ligand
##############################################


ASH features a function (**write_nonbonded_FF_for_ligand**) that allows one to quickly
create an XML forcefield file for any residue based on xTB-derived or DFT-derived atomic charges (CM5 charges) together with element-specific
Lennard-Jones parameters.

DFT-example:

.. code-block:: text

    from ash import *

    frag=Fragment(xyzfile="ligand.xyz")

    #Script to get nonbonded model parameters for a ligand
    orcatheory=ORCATheory(orcasimpleinput="!r2scan ZORA ZORA-def2-TZVP tightscf CPCM", numcores=8)

    write_nonbonded_FF_for_ligand(fragment=frag, resname="MCMtest", charge=0, mult=1,
        coulomb14scale=1.0, lj14scale=1.0, charge_model="CM5_ORCA", theory=orcatheory, LJ_model="UFF", charmm=True)

xTB-example:

.. code-block:: text

  from ash import *

  frag=Fragment(xyzfile="ligand.xyz")
  #Script to get nonbonded model parameters for a ligand
  write_nonbonded_FF_for_ligand(fragment=frag, charge=-3, mult=1, resname="MCMtest",
      coulomb14scale=1.0, lj14scale=1.0, charge_model="xTB", LJ_model="UFF", charmm=True)

**Options:**

- coulomb14scale and lj14scale parameters can be changed, depending on what other forcefield this ligand-forcefield will be combined with.
- charmm=True keyword writes the forcefield file so that it is compatible with the CHARMM36 forcefield (i.e. containing both a NonbondedForce and LennardJonesForce block)


**NOTES**

- Parameters will be derived for each atom in the XYZ-file. Symmetry is currently not incorporated and this means that very 
  similar atoms in the structure will have their own charge/LJ parameters. Since this is not always desired, the user
  should take care to combine and symmetrize the parameters in the XML-file manually.
- For a ligand bound to the protein, special care must be taken. Charges are best derived from a ligand structure with all metal ions
  coordinated (e.g. including an amino acid side chain) but then the calculation will contain those extra atoms.
  This requires manual tweaking of the final charges (make sure that the sum of atom charges add up to the correct total charge).