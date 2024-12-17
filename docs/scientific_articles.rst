Scientific articles using ASH
================================

ASH is a relatively new code but has been used for a few different research projects that led to publications. 
Below you can find various examples, including example ASH scripts that were used in the publications.

If you have used ASH in your research and you would like to mention it here, contact me!


###################################
MM and QM/MM protein studies
###################################

**Micro-second classical dynamics simulation of an artificial Mo-Cu hydrogenase**

`2023 J. Am. Chem. Soc. article <https://pubs.acs.org/doi/10.1021/jacs.3c01350>`_ 

This article used the OpenMM interface in ASH (:doc:`ORCA-interface`) to set up new solvated models of the Orange protein (Orp) and to explore possible binding
sites of an anionic S2MoS2CuS2MoS2 cluster to the protein. Long-time scale classical MD simulations were carried out, up to 1 microsecond (1000 nanoseconds) on a workstation with an Nvidia GPU card (platform='CUDA'). 
The X-ray structure of the apoprotein (`PDB-ID : 2WFB <https://www.rcsb.org/structure/2WFB>`_) together with a DFT-optimized (r2SCAN/def2-TZVP) structure of the metal cluster was used to set up the system which was fully solvated.
The metal cluster was described as a semi-rigid entity (bond constraints) with a nonbonded model (DFT-derived  CM5 charges and UFF Lennard-Jones parameters).
Periodic boundary conditions were used, all X-H bonds frozen and all water molecules were rigid. Before production NVT simulations, the system was minimized for 100 steps, gently heated up using 3 short NVT simulations
before running an NPT simulation until volume and density became stable. Production NVT simulations at 300 K were run for up to 1 microsecond using an integration timestep of 4 fs, Langevin friction coefficient of 1 ps-1 and hydrogen mass repartioning of 1.5 amu.
The simulations were found to be consistent with the anionic cluster binding to a positively charged Arg,Lys-containing pocket.

*Example setup script:*

.. toggle::
        
    .. code-block:: python


        from ash import *

        #Original raw PDB-file (no hydrogens, nosolvent)
        pdbfile="2wfb-modRB.pdb"
        #XML-file to deal with cofactor
        extraxmlfile="MCM.xml"
        #Setting possible manual protonation states.
        residue_variants={}
        # Setting up system via Modeller
        OpenMM_Modeller(pdbfile=pdbfile, forcefield='CHARMM36',
            extraxmlfile=extraxmlfile, watermodel="tip3p", pH=7.0, solvent_padding=10.0,
            ionicstrength=0.1, pos_iontype='Na+', neg_iontype='Cl-', residue_variants=residue_variants, use_higher_occupancy=True)

    where MCM.xml is the forcefield file for the Mo2CuS8 cluster:

    .. code-block:: text

        <ForceField>
        <AtomTypes>
        <Type name="MOX" class="Mo" element="Mo" mass="95.96"/>
        <Type name="CUX" class="Cu" element="Cu" mass="63.546"/>
        <Type name="SXB" class="S" element="S" mass="32.065"/>
        <Type name="SXT" class="S" element="S" mass="32.065"/>
        </AtomTypes>
        <Residues>
        <Residue name="MCM">
        <Atom name="MO1" type="MOX"/>
        <Atom name="MO2" type="MOX"/>
        <Atom name="CU1" type="CUX"/>
        <Atom name="S1" type="SXT"/>
        <Atom name="S2" type="SXT"/>
        <Atom name="S3" type="SXB"/>
        <Atom name="S4" type="SXB"/>
        <Atom name="S5" type="SXB"/>
        <Atom name="S6" type="SXB"/>
        <Atom name="S7" type="SXT"/>
        <Atom name="S8" type="SXT"/>
        </Residue>
        </Residues>
        <NonbondedForce coulomb14scale="1.0" lj14scale="1.0">
        <Atom type="MOX" charge="0.909" sigma="0.0" epsilon="0.0"/>
        <Atom type="CUX" charge="0.074" sigma="0.0" epsilon="0.0"/>
        <Atom type="SXT" charge="-0.691" sigma="0.0" epsilon="0.0"/>
        <Atom type="SXB" charge="-0.532" sigma="0.0" epsilon="0.0"/>
        </NonbondedForce>
        <LennardJonesForce lj14scale="1.0">
        <Atom type="MOX" sigma="0.2719022887764316" epsilon="0.234304"/>
        <Atom type="CUX" sigma="0.3113691019900486" epsilon="0.02092"/>
        <Atom type="SXT" sigma="0.3594776327696269" epsilon="1.146416"/>
        <Atom type="SXB" sigma="0.3594776327696269" epsilon="1.146416"/>
        </LennardJonesForce>
        </ForceField>


*Example MD script (opt,gentle heating, NPT and production NVT for 2.5 ns):*

.. toggle::

    .. code-block:: python

        from ash import *

        numcores=1

        #Defining list of lists of bond-constraints for cluster
        bondconstraints=[[1708,1711],[1708,1712],[1708,1714],[1708,1713],[1708,1710],[1710,1713],[1710,1714],[1710,1716],[1710,1715],[1710,1709],[1709,1715],[1709,1716],[1709,1718],[1709,1717]]

        #PDB-file to read topology from (and also initial coordinates)
        pdbfile="finalsystem.pdb"

        #Read coordinates from PDB-file only this time
        fragment=Fragment(pdbfile=pdbfile)

        #OpenMM object with constraints
        omm = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml", "MCM.xml"], pdbfile=pdbfile, periodic=True,
            platform='CUDA', numcores=numcores, autoconstraints='HBonds', constraints=bondconstraints, rigidwater=True)


        #MM minimization for 100 steps
        OpenMM_Opt(fragment=fragment, theory=omm, maxiter=100, tolerance=1)

        #Gentle heating up protocol
        OpenMM_MD(fragment=fragment, theory=omm, timestep=0.0005, simulation_steps=10, traj_frequency=1, temperature=1,
            integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajfilename='NVTtrajectorystepA', trajectory_file_option='DCD')
        OpenMM_MD(fragment=fragment, theory=omm, timestep=0.001, simulation_steps=50, traj_frequency=1, temperature=10,
            integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajfilename='NVTtrajectorystepB', trajectory_file_option='DCD')
        OpenMM_MD(fragment=fragment, theory=omm, timestep=0.004, simulation_steps=10000, traj_frequency=1, temperature=300,
            integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajfilename='NVTtrajectorystepC', trajectory_file_option='DCD')

        #NPT simulation until density and volume converges
        OpenMM_box_relaxation(fragment=fragment, theory=omm, datafilename="nptsim.csv", numsteps_per_NPT=10000,
                            volume_threshold=1.0, density_threshold=0.001, temperature=300, timestep=0.004,
                            traj_frequency=100, trajfilename='relaxbox_NPT', trajectory_file_option='DCD', coupling_frequency=1)

        #Classical NVT MD simulation for 2500 ps at 300 K
        OpenMM_MD(fragment=fragment, theory=omm, timestep=0.004, simulation_time=2500, traj_frequency=1000, temperature=300,
            integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD', trajfilename='NVTtrajectory')

        #Re-image trajectory so that protein is in middle
        MDtraj_imagetraj("NVTtrajectory.dcd", "final_MDfrag_laststep.pdb", format='DCD')




**QM/MM modelling of a CN-inhibited state of FeFe hydrogenase**

`2023 Chem. Sci. article <https://pubs.rsc.org/en/content/articlelanding/2023/sc/d2sc06098a>`_ 

This article used the QM/MM module of ASH together with the ORCA interface (:doc:`ORCA-interface`)
for the QM part and the OpenMM interface (:doc:`OpenMM-interface`) to setup
The OpenMMTheory interface used CHARMM-style forcefield files.

**QM/MM modelling of dinitrogen binding to redox states of nitrogenase**

`2023 Inorg. Chem. article <https://doi.org/10.1021/acs.inorgchem.2c03967>`_

This work, exploring dinitrogen binding to multiple redox states of the complex iron-molybdenum cofactor of nitrogenase 
used the QM/MM module of ASH together with the ORCA interface (:doc:`ORCA-interface`)
for the QM part and the OpenMM interface (:doc:`OpenMM-interface`) for the MM part (CHARMM36 forcefield with CHARMM files).
Broken-symmetry solutions in the QM-part were controlled by the ORCA interface (brokensym, HSmult, atomstoflip keywords, see :doc:`ORCA-interface`). 

**QM/MM modelling of Cu proteins**

Uncovering primary and secondary coordination sphere effects in S-nitrosylating azurin:
`2023 JACS article <https://pubs.acs.org/doi/10.1021/jacs.3c07399>`_

Second sphere variants of Type 1 Cu site in azurin:
`2024 J Phys. Chem. B article <https://pubs.acs.org/doi/full/10.1021/acs.jpcb.3c08194>`_



###################################
Highlevel WFT workflows
###################################

**High-level correlated WF density and ELF analysis using many ASH interfaces**

`Paper: The diradicaloid electronic structure of dialumenes: a benchmark study at the Full-CI limit <https://pubs.rsc.org/en/content/articlelanding/2024/cp/d4cp03005b>`_ 

In this article we performed high-level WFT calculations using both single-reference (MP2, CC) and multi-reference methods (CASSCF, MRCI)
as well as near-Full-CI calculations on a small dialumene, Al2H2.
ASH interfaces to ORCA, MRCC, CFour, pySCF, Dice (SHCI) and Block2 (DMRG) were used to conveniently carry out the calculations.
Energy calculations were carried out but also CCSD(T) and CCSDT geometry optimizations (using CFour and MRCC interfaces),
and difference density and ELF analysis was carried out using various correlated WF methods.
ASH was used to conveniently creating Molden files of the natural orbitals (derived from the correlated 1-particle density matrix) 
and an interfaces to Multiwfn was used to carry out the density and ELF analysis.

Note that a tutorial on high-level density analysis can be found here which covers most of the functionality used in the article:
:doc:`Highlevel-density-analysis`


**Multistep DLPNO-CCSD(T)/CBS workflow for a transition metal complex**

`2023 PCCP article <https://pubs.rsc.org/en/content/articlelanding/2023/cp/d2cp04715b>`_ 

This article used ORCA_CC_CBS_Theory (:doc:`module_highlevel_workflows`) functionality in ASH.
Below is a script that describes a recommended DLPNO-CCSD(T)/CBS workflow that worked well for this class of metallocenes
and should be reasonably reliable in general (assuming coupled cluster is reliable).
It uses a CBS(3/4) basis set extrapolation using the cc-pVnZ-DK basis set family, together with BP86 reference orbitals, 
DKH scalar relativistic Hamiltonian, PNO extrapolation using the cheaper approach and the cheaper T1 correction described
in th article.

*Example script:*

.. toggle::
        
    .. code-block:: python

        from ash import *
        numcores=24 #Number of cores reserved
        actualcores=16 #Number of cores used
        #Defining molecular fragments
        cpco0=Fragment(xyzfile="CpCo_0_gas.xyz", charge=0, mult=2)
        cpcoI=Fragment(xyzfile="CpCo_I_gas.xyz", charge=1, mult=1)
        # Defining species, stoichiometry and reaction specieslist=[cpco0,cpcoI]
        stoichiometry=[-1, 1]
        specieslist=[cpco0,cpcoI]
        reaction = Reaction(fragments=specieslist, stoichiometry=stoichiometry, unit='eV')
        #Defining a ORCA_CC_CBS_Theory object
        cc = ORCA_CC_CBS_Theory(elements=cpco0.elems, cardinals=[3,4], basisfamily="cc-dk", DFTreference="BP86", 
            DLPNO=True, CVSR=False, T1correction=True, T1corrbasis_size='Small', T1corrpnosetting='NormalPNOreduced', 
            numcores=actualcores, pnosetting="extrapolation", pnoextrapolation=[1e-6,3.33e-7,2.38,'NormalPNO'], 
            memory=20000, scfsetting="Verytightscf", relativity='DKH', SCFextrapolation=False)
        #Running reaction
        Singlepoint_reaction(theory=cc, reaction=reaction)



