======================================
OpenMM interface
======================================

OpenMM is an open-source molecular mechanics library written in C++. It comes with a handy Python interface that was easily adapted for use with ASH. It has been designed for both CPU and GPU codes.
The GPU code is particulary fast.



######################################
The OpenMMTheory theory level
######################################

OpenMM can run either on the CPU or on the GPU platform by specifying platform as either: 'CPU', 'OpenCL' or 'CUDA' depending on whether there is a GPU available on the computer. For platform='CPU' the numcores can additionally be specified in order to control the number of CPU cores. Alternatively the shell environment variable OPENMM_CPU_THREADS can be set to control the number of CPU cores that OpenMM will use. If neither numcores keyword is provided or OPENMM_CPU_THREADS variable set, OpenMM will use the number of physical cores present.

The OpenMMTheory class:

.. code-block:: python

    class OpenMMTheory:
        def __init__(self, printlevel=2, platform='CPU', numcores=None, Modeller=False, forcefield=None, topology=None,
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
                     ewalderrortolerance=1e-5, PMEparameters=None,
                     delete_QM1_MM1_bonded=False, applyconstraints_in_run=False,
                     constraints=None, restraints=None, frozen_atoms=None, fragment=None,
                     autoconstraints='HBonds', hydrogenmass=1.5, rigidwater=True):




It is possible to read in multiple types of forcefield files: AmberFiles, CHARMMFiles, GROMACSFiles or an OpenMM XML forcefieldfile.
Note: In rare cases OpenMM fails to read in Amber/CHARMM/GROMACS files correctly. In those cases the Parmed library may be more successful (use_parmed=True). Requires ParMed (pip install parmed).

Example creation of an OpenMMtheory object with CHARMM-files:

.. code-block:: python

    forcefielddir="/path/to/dir"
    topfile=forcefielddir+"top_all36_prot.rtf"
    parfile=forcefielddir+"par_all36_prot.prm"
    psffile=forcefielddir+"new-XPLOR-psffile.psf"
    openmmobject = OpenMMTheory(CHARMMfiles=True, psffile=psffile, charmmtopfile=topfile,
                               charmmprmfile=parfile)

Example creation of an OpenMMtheory object with GROMACS-files:

.. code-block:: python

    openmmobject = OpenMMTheory(GROMACSfiles=True, gromacstopdir="/path/to/gromacstopdir",
                    gromacstopfile="gromacstopfile.top", grofile="grofile.gro")

Example creation of an OpenMMtheory object with AMBER files:

.. code-block:: python

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile="/path/to/amberprmtopfile")

Example creation of an OpenMMtheory object with OpenMM XML file:

.. code-block:: python

    openmmobject = OpenMMTheory(xmlfiles=["example.xml"]) #File example.xml should be in dir
    #or
    openmmobject = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "specialresidue.xml"]) 
    #Here the "charmm36.xml" and "charmm36/water.xml" files are found in the OpenMM library.



An openmmtheory object can then be used to create a QM/MM theory object. See :doc:`module_QM-MM` page.

**Periodic boundary conditions:**

- If periodic boundary conditions are chosen (periodic=True) then the PBC box parameters are automatically found in the Amber PRMTOP file or the GROMACS Grofile or in the case of CHARMM-files they need to be provided: charmm_periodic_cell_dimensions
- PME parameters can be modified: PMEparameters=[alpha_separation,numgridpoints_X,numgridpoints_Y,numgridpoints_Z] 
- The ewalderrortolerance can be modified (default: 1e-5)
- The periodic nonbonded cutoff can be modified. Default: 12 Å
- Long-range dispersion correction can be turned on or off.
- The switching function distance can be changed. Default: 10 Å. Used for CHARMM and XML files.
- The box dimensions can also be modified by PBCvectors= keyword argument:
    Example: PBCvectors=[[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]]

######################################
Molecular Dynamics via OpenMM
######################################

It is possible to run MM molecular dynamics of system using the OpenMMTheory object created.
This is accomplished directly via the MD algorithms present in the OpenMM library.
The OpenMM_MD function takes as argument an ASH fragment, a teory object and then the user can select an integrator of choice, simulation temperature, simulation length, timestep, optional additional thermostat, barostat etc.

Most general options available in OpenMM are available in this interface. 
See OpenMM documentation page: http://docs.openmm.org/latest/userguide/application.html#integrators  for details about the integrators, thermostats, barostats etc.

- Available Integrators: Langevin, LangevinMiddleIntegrator, NoseHooverIntegrator, VerletIntegrator, VariableLangevinIntegrator, VariableVerletIntegrator
- Available Barostat: MonteCarloBarostat
- Optional additional thermostat: Anderson

.. code-block:: python

    def OpenMM_MD(fragment=None, theory=None, timestep=0.004, simulation_steps=None, simulation_time=None,
                  traj_frequency=1000, temperature=300, integrator='LangevinMiddleIntegrator',
                  barostat=None, pressure=1, trajectory_file_option='PDB', trajfilename='trajectory',
                  coupling_frequency=1,
                  anderson_thermostat=False,
                  enforcePeriodicBox=True, 
                  datafilename=None, dummy_MM=False, plumed_object=None, add_center_force=False,
                  center_force_atoms=None, centerforce_constant=1.0):


Options:

- fragment: ASH Fragment object.
- theory: should either be an ASH OpenMMTheory object, ASH QMMMTheory object (with mm_theory=OpenMMTheoryobject) or an ASH QMtheory.
- timestep: float (default: 0.001 ps). Size of timestep in picoseconds.
- simulation_steps: integer. Number of steps to take. (Use either simulation_steps or simulation_time)
- simulation_time: integer. Length of simulation time in ps. (Use either simulation_steps or simulation_time)
- temperature: integer (default:300). Temperature in Kelvin
- integrator: string (regular integrator or integrator+thermostat, e.g. 'LangevinMiddleIntegrator')
- barostat: string (e.g. 'MonteCarloBarostat'). Whether to add barostat to simulation for NPT simulations.
- coupling_frequency: frequency (ps^-1) to update thermostat/integrator. Applies to Nose-Hoover/Langevin.
- anderson_thermostat: Boolean (default: False)
- trajectory_file_option: 'PDB' or 'DCD'. Creates an ASCII PDB-trajectory or a compressed DCD trajectory.
- trajfilename : 'string'. Name of trajectory file (without extension)
- traj_frequency: integer (default: 1000). How often to write coordinates to trajectory file (every nth step)
- enforcePeriodicBox: Boolean (default: False). Option to fix PBC-image artifacts in trajectory.

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



Dealing with PBC image problems in trajectory. See https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#how-do-periodic-boundary-conditions-work
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

This function allows one to run multiple NPT simulations (constant pressure and temperature) in order to relax the periodic box dimensions
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

    numcores = 4

    pdbfile = "ash_inp.pdb"
    prmtopfile = "prmtop"

    frag = Fragment(pdbfile=pdbfile)

    #List of all non-H atoms
    allnonHatoms=frag.get_atomindices_except_element('H')

    openmmobject = OpenMMTheory(Amberfiles=True, amberprmtopfile=prmtopfile, periodic=True,
            platform='CPU', autoconstraints=None, rigidwater=False, frozen_atoms=allnonHatoms)



    OpenMM_MD(fragment=frag, theory=openmmobject, timestep=0.001, simulation_steps=100,
            traj_frequency=1, temperature=300, integrator="LangevinIntegrator",
            coupling_frequency=1, trajectory_file_option="PDB")


######################################
System setup via OpenMM: Modeller
######################################

OpenMM features a convenient PDBfixer program (https://github.com/openmm/pdbfixer) and a Modeller tool (http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.modeller.Modeller.html)
that is capable of setting up a new biomolecular system from scratch. See also: http://docs.openmm.org/7.2.0/userguide/application.html#model-building-and-editing . ASH features a highly convenient interface to these programs and allows near-automatic system-setup for favorable systems.

.. code-block:: python

    def OpenMM_Modeller(pdbfile=None, forcefield=None, xmlfile=None, waterxmlfile=None, watermodel=None, pH=7.0,
                        solvent_padding=10.0, solvent_boxdims=None, extraxmlfile=None, residue_variants=None,
                        ionicstrength=0.1, pos_iontype='Na+', neg_iontype='Cl-'):


The OpenMM_Modeller function returns an ASH OpenMMTheory object that can be used directly as theory level for future calculations.
OpenMM_Modeller will also print various PDB-files associated with each step of the setup (H-addition, solvation, ionization etc.).
And an XML file associated with the system that can be used to create future OpenMMtheory objects from.

Lysozyme example (simple, no modifications required):

.. code-block:: python

    from ash import *

    #Original raw PDB-file (no hydrogens, nosolvent)
    #Download from https://www.rcsb.org/structure/1AKI
    pdbfile="1aki.pdb"


    #Defining residues with special user-wanted protonation states for residues in each indicated chain
    #Dictionary of dictionaries with the chainname (e.g. 'A','B') acting as keys for the outer dictionary and the resids being keys for the inner dictionary
    #Example: residue_variants={'A':{0:'LYN', 17:'CYX', 18:'ASH', 19:'HIE', 20:'HID', 21:'GLH' }}
    #resid 1: neutral LYS, resid 17, deprotonated CYS, resid 18 protonated ASP, 
    #resid 19 epsilon-protonated HIS, resid 20 delta-protonated HIS, 21 protonated GLU.
    residue_variants={}

    #Setting up new system, adding hydrogens, solvent, ions and defining forcefield, topology
    openmmobject = OpenMM_Modeller(pdbfile=pdbfile, forcefield='CHARMM36', watermodel="tip3p", pH=7.0, 
        solvent_padding=10.0, ionicstrength=0.1, residue_variants=residue_variants)

    #MM minimization for 100 steps
    OpenMM_Opt(fragment=ashfragment, theory=openmmobject, maxiter=100, tolerance=1)

    #Classical MD simulation for 10 ps
    OpenMM_MD(fragment=ashfragment, theory=openmmobject, timestep=0.001, simulation_time=10, traj_frequency=100, temperature=300,
        integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')


If the protein contains nonstandard residues (e.g. metallocofactors) that are not present in a typical protein forcefield (OpenMM_Modeller will exit with errors),
then these need to be provided using the extraxmlfile option.

.. code-block:: python

    openmmobject = OpenMM_Modeller(pdbfile=pdbfile, forcefield='CHARMM36', watermodel="tip3p", pH=7.0, 
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


See e.g. https://education.molssi.org/mm-tools/01-introduction/index.html for information on the format of the XML file.

Common error messages encountered when reading in user-defined XML-files:

-**ValueError: No template found for residue X (YYY).  This might mean your input topology is missing some atoms or bonds, or possibly that you are using the wrong force field.**

*This means that the parser encountered a completely unknown residue. You might have forgotten to read in the XML file to OpenMM_Modeller or the resname is not the same in the
PDBfile as in the XML file.*

- **ValueError: Found multiple definitions for atom type: X**  :  

*This means that the atomtypes you defined inside <AtomTypes> </AtomTypes> are not unique. Either you have accidentally two lines with the same name for atomtype or the general forcefield (e.g. CHARMM) already contains an atomtype definition with this name. In that case, choose a unique name.*

- **KeyError: 'SXM'**  :  

*Possibly means that your atomname definition points to an atomtype-name that does not exist*


- **ValueError: No template found for residue X (YYY).  The set of atoms matches YYY, but the bonds are different.  Perhaps the chain is missing a terminal group?'**  :  

*This means there is some mismatch between the information present in the PDB-file and the information in the XML-file you provided.
It's possible that the PDB-file contains connectivity statements at the bottom of the PDB-file (CONE lines) but no bond information is present in the XML file.
Solution: Either add the missing bond to the residue definition so that it matches the CONE lines or simply delete the CONE information that you don't need.*


Advanced example (additional forcefield parameters required):

.. code-block:: python

    from ash import *


    #TODO


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

ASH also features a function to solvate a small molecule automatically. This also makes use of the Modeller functionality of OpenMM but is intended to be used for molecules for where forcefield parameters are typically not available: e.g. metal complexes. Instead of regular forcefield parameters, nonbonded parameters (charges and Lennard-Jones parameters) are defined for the solute (used for classical and QM/MM simulations) which can be used to perfrom classical MM dynamics or QM/MM dynamics.

See also :doc:`Explicit-solvation` workflow page.


.. code-block:: python

    def solvate_small_molecule(fragment=None, charge=None, mult=None, watermodel=None, solvent_boxdims=[70.0,70.0,70.0], 
                               nonbonded_pars="CM5_UFF", orcatheory=None, numcores=1):

The solvate_small_molecule function reads in an ASH fragment, as well as charge and multiplicity, name of watermodel (e.g. "TIP3P"), size of solvent box, option for how the nonbonded parameters should be prepared, an optional ORCATheory object and optional numcores.

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


