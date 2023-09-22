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
Note: OpenMM_box_relaxation is an alias for penMM_box_equilibration

.. code-block:: python

    def OpenMM_box_equilibration(fragment=None, theory=None, datafilename="nptsim.csv", numsteps_per_NPT=10000,
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
These large initial forces are usually responsible for the system blowing up in the beginning (error messages of e.g. 'Particle number is NaN' etc.).
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
Gentle WarmupMD
######################################
A function to gently warm up a newly setup system to a target temperature. 
Can also be used to help diagnose why an MD simulation crashes (reports initial high atomic forces as well as root-mean-square fluctuations).

.. code-block:: python

  def Gentle_warm_up_MD(theory=None, fragment=None, time_steps=[0.0005,0.001,0.004], steps=[10,50,10000], 
      temperatures=[1,10,300], check_gradient_first=True, gradient_threshold=100, use_mdtraj=True)

The minimization algorithm in **OpenMM_Opt** described above can occasionally fail to reduce the main problematic forces
present in a newly setup system. It can even occasionally crash during the minimization without revealing the cause.
Starting MD simulations directly can also lead to crashes without helpful error messages.
The reason for these crashes is usually due to these large forces resulting in high atom velocities (or similar problems in the minimization) 
which causes the system to blow up (error messages such as 'Particle number is NaN' etc.).
Furthermore, the OpenMM minimization algorithm currently does not report any progress on the minimization (`see Github issue <https://github.com/openmm/openmm/issues/1155>`_)

An alternative (or addition) to a minimization is to instead start MD simulations using a very low temperature and small timesteps and then gradually increase the temperature and timestep.
Such a protocol can work where a minimization fails or at the very least it can provide information about what part of the system has these large forces.

ASH provides a convenient function, Gentle_warm_up_MD, that can be called to do such a gentle warmup MD in a few steps.
In addition, the function reports the largest atom forces present in the initial geometry and will report atoms with the largest root-mean-square fluctuations
after each MD simulation it performs (requires mdtraj to be installed).

To use it, you simple call the function with the OpenMMTheory object and Fragment object as input.

.. code-block:: python

  Gentle_warm_up_MD(theory=openmmobject, fragment=frag, time_steps=[0.0005,0.001,0.004], 
                    steps=[10,50,10000], temperatures=[1,10,300])

By default, the function will perform a warmup protocol consisting of:

- 10-step MD simulation with a 0.5 fs timestep (0.0005 ps) at temperature 1 K
- 50-step MD simulation with a 1.0 fs timestep (0.001 ps) at temperature 10 K
- 10000-step MD simulation with a 4.0 fs timestep (0.004 ps) at temperature 300 K

This protocol may be sufficient to warm up your system without it blowing up but the protocol can also be modified in any way you like.
By adding values to the lists above you add extra simulations, change the steps, change the temperatures, timesteps etc.
A DCD trajectory is written for each MD simulation and each snapshot is written to disk (traj_frequency=1) which can be visualized in VMD.

Gentle_warm_up_MD will by default use `mdtraj <https://www.mdtraj.org>`_ to image trajectories
for better visualization as well as calculate root-mean-square fluctuations.  mdtraj can be installed like this: pip install mdtraj

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
    <Atom type="FEX" charge="3.0" sigma="1.0" epsilon="0.0"/>
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


#######################################################
Setting up a protein system with implicit solvation
#######################################################

It is also possible to use OpenMM_Modeller to setup a protein system with an implicit solvent instead of explicit.
Note that the protein-forcefield must be compatible with the chosen implicit solvent.
See `Open MM documentation <http://docs.openmm.org/latest/userguide/application/02_running_sims.html#implicit-solvent>`_ for more information.

Below is an example for setting up a protein using Amber14 and the OBC2 implicit solvation model.
.. code-block:: python

  from ash import *

  OpenMM_Modeller(pdbfile="combmol.pdb", forcefield="Amber14", implicit=True, implicit_solvent_xmlfile="implicit/obc2.xml")

The system will be setup as usual using the steps in Open_Modeller but no explicit solvent or counterions will be added.
Additionally periodicity will not be assumed during the creation of the files as implicit solvation calculations should be run without PBC.

.. code-block:: python

  from ash import *

  #Read in previous system from OpenMM_Modeller
  frag = Fragment(pdbfile="finalsystem.pdb")

  #Create an OpenMMTheory object without PBC. Here using a cutoff of 20 Angstroms for nonbonded interactions
  omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "implicit/obc2.xml", "gaff_ligand.xml"], pdbfile="finalsystem.pdb",
      periodic=False, platform='OpenCL', nonbondedMethod_noPBC='CutoffNonPeriodic', nonbonded_cutoff_noPBC=20)

  #NVT MD simulation for 1000 ps = 1 ns
  OpenMM_MD(fragment=frag, theory=omm, timestep=0.004, simulation_time=10, traj_frequency=100, temperature=300,
      integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajfilename='NVTtrajectory',trajectory_file_option='DCD')

MD simulations with an implicit solvation model can have their advantages as they should run considerably quicker.
While the implicit solvation model is not as accurate as explicit solvation, it can be a good starting point for a system that is later simulated with explicit solvent.


#######################################################
Create forcefield for ligand / small molecule
#######################################################

Often one wants to perform a classical or QM/MM simulation of a small molecule in solution (either as part of a biomolecular system or on its own)
but one lacks forcefield parameters to do so. One has typically 2 options for how to proceed in this case:

- Create only a nonbonded forcefield (charges and Lennard-Jones parameters) for the small molecule.
- Create a full forcefield for the small molecule (bonded and nonbonded parameters).

The first option (nonbonded only) is sufficient if one primarily intends to perform QM/MM simulations where the molecule will always be in the QM-region.
This may also be the only easy option if the molecule is inorganic (e.g. a metal complex) where forcefield parameterization is less straightforward. 
The nonbonded forcefield can also be used in classical simulation if one makes sure the ligand is rigid (all bonds constrained, possibly angles and dihedrals as well).
See next section below: **write_nonbonded_FF_for_ligand**

The second option (full forcefield) is generally better and is required if one wants to perform classical simulations where the molecule is flexible.
ASH features a function (**small_molecule_parameterizor**) that allows one to expedite this process with the help of the `openmm-forcefields <https://github.com/openmm/openmmforcefields>`_, 
that provides a convenient way of getting forcefield parameters from the `GAFF <https://ambermd.org/antechamber/gaff.html>`_ and `OpenFF <https://openforcefield.org>`_ projects. 
The limitation is that this option is primarily available for organic or drug-like molecules.
Additionally these small-molecule forcefields are intended to be only used together with Amber biomolecular forcefields (if your system also includes protein/DNA).


##############################################
write_nonbonded_FF_for_ligand
##############################################

.. code-block:: python

  def write_nonbonded_FF_for_ligand(fragment=None, charge=None, mult=None, coulomb14scale=1.0, lj14scale=1.0, 
    ff_type="CHARMM", charge_model="CM5", theory=None, LJ_model="UFF", resname="LIG", numcores=1):


ASH features a function (**write_nonbonded_FF_for_ligand**) that allows one to quickly create an OpenMM-style XML forcefield file for any ligand/molecule
with only nonbonded parameters specified which can be sufficient for QM/MM simulations or classical simulations where the ligand/molecule is rigid (all bonds constrained).

One can choose to derive the atom charges from either an xTB-calculation (using the xTB interface) or a DFT-calculation (ORCA interface).
The charge_model options are: CM5 charges or DDEC3/DDEC6 charges (requires DDEC3/DDEC6).
The Lennard-Jones parameters can either come from UFF (very crude: element-specific LJ parameters) or via DDEC3/DDEC6 population analysis.


.. warning:: It is up to you the user to make sure that the nonbonded parameters from this procedure are sensible and compatible with other molecules present in your system (described by another forcefield).
  You may have to change the parameters manually 

*Example:*

.. code-block:: text

    from ash import *

    frag=Fragment(xyzfile="ligand.xyz")

    #Script to get nonbonded model parameters for a ligand
    orcatheory=ORCATheory(orcasimpleinput="!r2scan ZORA ZORA-def2-TZVP tightscf CPCM", numcores=8)

    write_nonbonded_FF_for_ligand(fragment=frag, resname="MCMtest", charge=0, mult=1,
        coulomb14scale=1.0, lj14scale=1.0, charge_model="CM5_ORCA", theory=orcatheory, LJ_model="UFF", ff_type="CHARMM")


**Options:**

- charge_model: Options are 'CM5', 'xTB', 'DDEC3', 'DDEC6'
- LJ_model: Options are 'UFF', 'DDEC3', 'DDEC6'
- The ff_type keyword (options: 'CHARMM', 'AMBER', 'None'), writes the forcefield file so that it is compatible with the CHARMM, Amber biomolecular forcefields. Choose 'None' if not needed.
- coulomb14scale and lj14scale parameters can be changed, depending on what other forcefield this ligand-forcefield will be combined with  (OpenMM requires compatibility)

**NOTES**

- Parameters will be derived for each atom in the XYZ-file. Symmetry is currently not incorporated and this means that very 
  similar atoms in the structure will have their own charge/LJ parameters. Since this is not always desired, the user
  should take care to combine and symmetrize the parameters in the XML-file manually.
- For a ligand bound to the protein, special care must be taken. Charges are best derived from a ligand structure with all metal ions
  coordinated (e.g. including an amino acid side chain) but then the calculation will contain those extra atoms.
  This requires manual tweaking of the final charges (make sure that the sum of atom charges add up to the correct total charge).
- DDEC3/DDEC6: Both atom charges and LJ parameters can be determined from a DFT-calculation and a DDEC3/DDEC6 population analysis using the Chargemodel. This options has not been well tested and requires external programs (Chargemol and mol2aim)


##############################################
small_molecule_parameterizor
##############################################

.. code-block:: python

  def small_molecule_parameterizor(xyzfile=None, pdbfile=None, molfile=None, sdffile=None, smiles_string=None,
                                  forcefield_option='GAFF', gaffversion='gaff-2.11',
                                  output_xmlfile="ligand.xml", openff_file="openff-2.0.0.offxml",
                                  expected_coul14=0.8333333333333334, expected_lj14=0.5):

**small_molecule_parameterizor** allows you to quickly create an OpenMM XML forcefield file with bonded and nonbonded parameters for your molecule.
You can choose between two general forcefields: `GAFF <https://ambermd.org/antechamber/gaff.html>`_  or `OpenFF <https://openforcefield.org>`_. 
Different GAFF and OpenFF versions are also available. The limitation is that creating the small-molecule forcefield from these general forcefields can only be done for "organic" chemical elements (H,C,N,O,S,P,F,Cl,Br,I; also ions such as 
Li+, Na+, K+, Rb+, F-, Cl-, Br-, and I-).
These small-molecule forcefields are intended to be only used together with Amber biomolecular forcefields (if your system also includes protein/DNA).

**small_molecule_parameterizor** is very easy to use. You simply need to provide a molecular structure. This can be an XYZ-file, PDB-file, MDL Mol-file, SDF-file but it can also
be a `SMILES string <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system>`_ .

The program depends on a few Python libraries that have to be installed when prompted: `openmmforcefields <https://github.com/openmm/openmmforcefields>`_, `parmed <https://github.com/ParmEd/ParmEd>`_, `openff-toolkit <https://github.com/openforcefield/openff-toolkit>`_ and `OpenBabel <http://openbabel.org/wiki/Main_Page>`_
ASH will tell you which libraries are missing and how to install them when you try to use the function.

*Example using GAFF*

.. code-block:: python

  from ash import *
  #Creating forcefield for nitrate using GAFF. Here providing a SMILES string as input
  small_molecule_parameterizor(forcefield_option="GAFF", smiles_string="[N+](=O)([O-])[O-]")

*Example using OpenFF*

.. code-block:: python

  from ash import *
  #Creating forcefield for nitrate using OpenFF. Here providing xyz-file as input
  small_molecule_parameterizor(forcefield_option="OpenFF", xyzfile="no3.xyz"


The output is an XML-file that can then be used as input to **OpenMMTheory**, **OpenMM_Modeller** or **solvate_small_molecule** functions (see below).

.. warning:: The XML-file created by this function will contain bonded parameters and it is thus important that the topology of the molecule is available when using the XML-file
  together with OpenMM. Otherwise, the pairing of molecule and small-molecule forcefield in the XML-file will not work. As OpenMM will typically get the topology from a PDB-file you must ensure 
  to have a PDB-file that contains CONECT lines at the bottom of the PDB-file that describes the connectivity of the small molecule. A PDB-file with connectivity is automatically created if you read in an XYZ-file
  to small_molecule_parameterizor above. You can also use the  **xyz_to_pdb_with_connectivity** function.


######################################
Small molecule solvation
######################################

.. code-block:: python

  def solvate_small_molecule(fragment=None, charge=None, mult=None, watermodel=None, solvent_boxdims=[70.0, 70.0, 70.0],
                            xmlfile=None):

ASH also features a function to solvate a small molecule automatically. This also makes use of the Modeller functionality of OpenMM but is a bit simpler.
It requires reading an ASH fragment, selection of a water model and an XML-file containing the small-molecule forcefield.
The XML-file can come from either **write_nonbonded_FF_for_ligand** or **small_molecule_parameterizor**
The size of the solvent box can be modified as required (default 70x70x70 Angstrom).

Options:

- watermodel (string). Can be: 'TIP3P' only for now
- xmlfile (string). Name of the XML-file containing either a nonbonded or full forcefield of the molecule.
- solvent_boxdims (list of floats). Cubic box dimensions in Angstrom.


*Example:*

.. code-block:: python

    from ash import *

    numcores=4
    #Molecule definition
    mol=Fragment(xyzfile="3fgaba.xyz", charge=0, mult=1)

    #Solvate molecule (70x70x70 Å TIP3P box)
    forcefield, topology, ashfragment = solvate_small_molecule(fragment=mol, watermodel='tip3p', solvent_boxdims=[70,70,70])


The output of **solvate_small_molecule**  are coordinate files: "system_aftersolvent.pdb" and "system_aftersolvent.xyz" .

ASH will print information about how to create an OpenMMTheory for the system but typically it would look like this:

.. code-block:: python

    from ash import *
    #Read in coordinates: either XYZ-file or PDB-file
    fragment = Fragment(xyzfile="system_aftersolvent.xyz")
    #Create an OpenMMTheory object based on PDB-file and XML-files for water and small-molecule
    #Note: that the XML-file for the solvent may be different (CHARMM-style, Amber-style or OpenMM-style)
    openmmobject =OpenMMTheory(xmlfiles=["molecule.xml", "amber/tip3p_standard.xml"], pdbfile="system_aftersolvent.pdb", 
            periodic=True, rigidwater=True, autoconstraints='HBonds')


Additionally the function returns an OpenMM forcefield object, an OpenMM topology and an ASH fragment. These could also be used to create an OpenMMTheory object, 
but would have to be performed in the same script as **solvate_small_molecule**

.. code-block:: python

    #Creating new OpenMM object from forcefield, topology and and fragment
    openmmobject =OpenMMTheory(numcores=numcores, topoforce=True, forcefield=forcefield, topology=topology, 
                    periodic=True, autoconstraints='HBonds', rigidwater=True)

The OpenMMTheory object can then be used on its own or can be combined with a QM theory to define a QM/MM theory object etc.
See :doc:`Explicit-solvation` workflow for more information on how to use **solvate_small_molecule** in a multi-step workflow.