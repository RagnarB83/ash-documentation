Tutorial: Running fast protein QM/MM MD simulations in ASH
==================================================================

In many research projects the question sometimes arises whether it might be feasible to
perform a long enough MD simulation (and perhaps enhanced sampling MD) of the system at some kind of quantum level of theory. 
Often such simulations are simply unfeasible at the direct QM-level, however, hybrid QM/MM methodology
can make such simulations possible, as one can focus the expensive QM calculation on a small important region of the system.

However, what is the best way to perform such a QM/MM MD simulation? 

- What program to use ? 
- What QM theory level? 
- Is fast communication between QM and MM program important ?
- Is there an interface already available between QM and MM program?
- Does the MM calculation need to be parallelized ?
- How to think about CPU parallelization and scaling?
- What about using GPUs instead of CPUs ?

This tutorial intends to give insight into how to choose a good QM/MM MD protocol
and to demonstrate how ASH is ideally suited for performing such simulations due to
the flexibility offered by the general QM/MM approch available and the many QM-code interfaces available.
The flexibility offered by ASH means that it is easy to switch from running classical MM (OpenMM), 
semi-empirical QM/MM (MNDO, xTB), Gaussian-basis DFT/MM (ORCA, pySCF, NWChem etc.), GPU-based codes (TeraChem, QUICK, pyscf), 
mixed Gaussian-planewave/MM (CP2K), and even post-HF/MM (ORCA,CFour, MRCC) etc.


################################################################################
0. Lysozyme system
################################################################################

We will use a previously set-up system of the solvated lysozyme enzyme as an example here.
The files can be found in the `Github repository <https://github.com/RagnarB83/ash/tree/master/examples/OpenMM_Modeller-setups/lysozyme>`_ 
The system contains 25825 atoms in total (1960 protein atoms, 23865 solvent/ions atoms) in a cubic simulation box.
In the QM/MM MD simulations we will use a QM-region containing the sidechains of GLU35, ASP52 and GLN57 which roughly corresponds
to the active site residues of this enzyme. This is a QM-region of 26 atoms +3 H-linkatoms = 29 QM atoms (14 of which are non-H).

In the examples we will typically run a 10 ps QM/MM MD simulations using a 1 fs timestep, i.e. 10000 steps.
Note that 1 fs timestep requires the use of X-H constraints (for MM region).
 
Neither the simulation lengths nor the QM-regions are intended to be realistic in terms of solving some interesting chemical problem. 
We are only using this setup to demonstrate the flexibility of ASH and what can be achieved in terms of speed of simulations for 
a typical closed-shell system with a fairly small QM-region. 
Furthermore the timings from these benchmarks should be taken with a grain of salt as they will depend on specific hardware
and software-setup on each computer.

The simulations will mostly use the ASH-script as shown below, where only the QM-theory object definition changes.

.. toggle::

    .. code-block:: python

        from ash import *

        #Cores to use for QM-program (and MM-program if platform='CPU')
        numcores=1
        #MM-platform options: 'CPU', 'OpenCL', 'CUDA'
        MM_platform='CUDA'
        #Traj frequency
        traj_frequency=1000

        #PDB-file and fragment
        pdbfile="Setup_equilibrated.pdb"
        frag = Fragment(pdbfile=pdbfile)

        #OpenMM object
        omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "amber14/tip3pfb.xml"], pdbfile=pdbfile, periodic=True,
            autoconstraints='HBonds', rigidwater=True, platform=MM_platform, numcores=numcores)

        #QM-theory object
        qm = MNDOTheory(method="AM1", numcores=numcores)

        #QM-region definition
        GLU35=list(range(539, 547+1))
        ASP52=list(range(778,783+1))
        GLN57=list(range(856,866+1))
        qmatoms = GLU35+ASP52+GLN57

        #QM/MM object
        qm_mm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=frag,
                qmatoms=qmatoms, printlevel=2, qm_charge=-2, qm_mult=1)

        #Classical MD simulation for 1000 ps
        OpenMM_MD(fragment=frag, theory=qm_mm, timestep=0.001, simulation_time=10, traj_frequency=traj_frequency, temperature=300,
            integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')


When reporting timings we use the time for the **OpenMM_MD run** step, printed at the end of the ASH output, 
(this time does not include the setup and creation of the Theory objects, done only once per job).

################################################################################
1. Lysozyme: Classical simulations
################################################################################

It is useful to first see the speed that can be obtained at the classical MM level of theory using the OpenMM library.
Since the simulations run so fast at the MM level and to get more reliable statistics, 
we ran these MM simulations for both 10 ps and 100 ps.

Simulations were run on a Linux computing node:
24-core Intel(R) Xeon(R) Gold 5317 CPU @ 3.00GHz with a Nvidia A100 80GB VRAM GPU.

Script used:

.. toggle::

    .. code-block:: python

        from ash import *

        #Cores to use for QM-program (and MM-program if platform='CPU')
        numcores=1
        #MM-platform options: 'CPU', 'OpenCL', 'CUDA'
        MM_platform='CUDA'
        #traj_frequency
        traj_frequency=10000

        #PDB-file and fragment
        pdbfile="Setup_equilibrated.pdb"
        frag = Fragment(pdbfile=pdbfile)

        #OpenMM object
        #CHARMM:
        #omm = OpenMMTheory(xmlfiles=["charmm36.xml", "charmm36/water.xml"], pdbfile=pdbfile, periodic=True,
        #    autoconstraints='HBonds', rigidwater=True, platform=MM_platform, numcores=numcores)
        #Amber:
        omm = OpenMMTheory(xmlfiles=["amber14-all.xml", "amber14/tip3pfb.xml"], pdbfile=pdbfile, periodic=True,
            autoconstraints='HBonds', rigidwater=True, platform=MM_platform, numcores=numcores)

        #Classical MD simulation for 1000 ps
        OpenMM_MD(fragment=frag, theory=omm, timestep=0.001, simulation_time=100, traj_frequency=traj_frequency, temperature=300,
            integrator='LangevinMiddleIntegrator', coupling_frequency=1, trajectory_file_option='DCD')

--------------------------------------------------------------
CPU vs. GPU
--------------------------------------------------------------

We used the Amber14 forcefield and ran simulations by changing the *platform* keyword between 'CPU', 'OpenCL' and 'CUDA'.
For 'CPU' we also changed the number of cores in the OpenMMTheory object (*numcores* keyword).

==============  ========================= =============================== ================================
Hardware          Time (sec) for 10 ps     Time (sec) for 100 ps           Ave. time (sec) per timestep
==============  ========================= =============================== ================================
 CPU(1 cores)        1434                        13866                      0.1387
 CPU(4 cores)        562                         5371                       0.0537
 CPU(8 cores)        344                         4100                       0.0410
 CPU(16 cores)       345                         3072                       0.0307 
 CPU(24 cores)       388                         3758                       0.0376
 OpenCL              11                          42                         0.0004
 CUDA                8                           35                         0.0003
==============  ========================= =============================== ================================

As the results show, there is a massive speed difference between running classical simulations on the CPU vs. GPU. 
While using multiple CPU cores help speed up the simulation (up to approx. 8-16 cores here), we can't approach 
the speed of running GPU-tailored code on a single GPU (~100 times faster).
There is additionally a small advantage in utilizing the CUDA GPU-code in OpenMM (only for Nvidia GPUs) 
rather than the more general OpenCL GPU-code.

.. note:: It should be noted that OpenMM is primarily designed for running fast on the GPU. There are MM codes that have faster CPU execution than OpenMM but generally GPU MM implementations are faster than CPU implementations and OpenMM is one of the fastest.

--------------------------------------------------------------
Forcefield: Amber vs. CHARMM
--------------------------------------------------------------

It is easy to switch between CHARMM and Amber forcefields (see script above) for this simple protein setup (no ligand or cofactor present)
and so we can compare whether there is a difference in speed when using a CHARMM forcefield vs. an Amber forcefield. 
Calculations below ran either on the GPU(CUDA) or CPU(8 cores).

=================  =============================== =============================== ================================
Forcefield          Time (sec) for 10 ps             Time (sec) for 100 ps            Ave. time (sec) per timestep
=================  =============================== =============================== ================================
 Amber  (GPU)              9                             35                               0.00035
 Amber  (8 CPUs)         336                           3143                                0.03143
 CHARMM (GPU)              9                             44                               0.00044
 CHARMM (8 CPUs)           843                           8757                                0.08757
=================  =============================== =============================== ================================

It turns out that there is a speed difference, with Amber being here a bit faster than CHARMM at the GPU-level (~25 %) and quite a bit faster at CPU-level (2.5x).
The reason for this is likely due to differences in the Lennard-Jones terms in CHARMM vs. Amber and how they are implemented and 
handled in OpenMM (see `Github discussion <https://github.com/openmm/openmm/issues/4311>`_ ). 
This difference is more severe on the CPU platform and thus may be worth taking into account when choosing a forcefield for a classical simulation, if CPU is the only hardware option.


--------------------------------------------------------------
Trajectory writing
--------------------------------------------------------------

It is also good to be aware of in classical simulations, that because each timestep is so fast to compute, any
additional procedure during each timestep may strongly affect the performance by increasing I/O (reading and writing to disk).
Here we show how the act of writing the geometry to a trajectory file after each timestep, affects the overall speed.
Trajectory-writing is always active but the frequency of writing is controlled by the *traj_frequency* keyword.
If *traj_frequency=1* then we write a frame to trajectory every single step (this is slow would produces very large trajectory files, possibly filling up the scratch),
while if *traj_frequency=10000* we write to the trajectory every 10000 steps (little cost and smaller files).

The table below shows that as long as *traj_frequency* is 1000 or larger then no severe speed-penalty is encountered.
Calculations used Amber forcefield and ran on the GPU(CUDA).

===================  ================================ ================================ ================================
traj_frequency          Time (sec) (for 10 ps)          Time (sec) (for 100 ps)          Ave. time (sec) per timestep
===================  ================================ ================================ ================================
 1                     951                             8978                               0.0898
 2                     458                             4516                               0.0452
 5                     184                             1890                               0.0189
 10                    97                              971                                0.0097
 100                   17                              131                                0.0013
 1000                  8                               72                                 0.0007
 10000                 7                               63                                 0.0006
===================  ================================ ================================ ================================


--------------------------------------------------------------
Timestep
--------------------------------------------------------------

Finally we note that in pure MM MD simulations it is easy to use larger timesteps than possible in QM or QM/MM MD simulations
as long as appropriate constraints are employed at the MM level.
We can typically gain a speed-up of approx. 3-4 by going from a 1 fs timestep to a 3-4 fs timestep.
This is essentially without loss in accuracy as long as the water model is rigid (rigidwater=True) and 
all X-H bonds are constrained (autoconstraints='HBonds') and HMR (increased hydrogen mass) is being utilized.
These tricks are typically not possible at the QM/MM level.

===================  ================================
Timestep (fs)          Time (sec) (for 100 ps)        
===================  ================================
 1                     34.94                            
 2                     20.27                              
 3                     15.07                              
 4                     12.97                               
===================  ================================

Calculations in table used Amber forcefield and ran on the GPU(CUDA).

################################################################################
2. Lysozyme: QM/MM MD using semi-empirical methods
################################################################################

When you switch from MM to QM/MM you should expect a massive drop in speed. This is because of 2 factors:

A. The slower speed of the QM energy+gradient calculation that has to be performed in each simulation step.
B. A regular MM simulation in ASH keeps positions and velocities in memory while running efficient C++/OpenCl/CUDA code; meanwhile a QM/MM simulation will by necessity do the simulation step-by-step, with data exchange in each step, calling the QM and MM program and even having some I/O (reading and writing to disk).

We can see some of this speed-drop from factor B that occurs if we switch from running a regular MM MD (with all positions and velocities in memory
and simulation proceeding by the compiled code) to a simulation where each simulation step is iterated at the Python-level and each MM-call is performed explicitly by ASH.
The latter can be accomplished by providing the *dummy_MM* keyword to the **OpenMM_MD** function.

===================  ================================ ================================ ================================
QM-method             Time (sec) (for 10 ps)           Time (sec) (for 100 ps)           Ave. time (sec) per timestep
===================  ================================ ================================ ================================
 Regular                          7                              35                         0.00035
 Dummy-MM                         23                            192                         0.00192
===================  ================================ ================================ ================================


This speed-penalty factor of ~5 (going from 0.35 ms 1.92 ms, per step) is caused by the need for data-exchange between the Python, 
C++ and CUDA/OpenCL layers of OpenMM and ASH. While this looks at first glance like an issue, it is actually a very small penalty compared to the cost of the QM-step (factor A) that will be added on top of this cost.

The fastest way to run a QM step in QM/MM is using a semi-empirical QM Hamiltonian that avoids calculating all the expensive 
2-electron integrals that are otherwise present in regular HF or DFT theory. Here we compare the old-school semi-empirical AM1 method from 1985 (using the fast MNDO code), PM3
and the tightbinding GFN1-xTB method (using the xTB code).

Switching from an MM simulation to a QM/MM simulation is very simple in ASH. 
We simply have to combine a QM-theory object with an MM-theory object
in a QMMMTheory object.

The MNDOTheory object (see :doc:`MNDO-interface` documentation ) is created like this:

.. code-block:: python
    
    #Note: MNDO does not run in parallel
    qm = MNDOTheory(method="AM1", numcores=1)

while the xTBTheory object (see :doc:`xTB-interface`) is created like this:

.. code-block:: python
    
    #QM object
    qm = xTBTheory(xtbmethod="GFN1", numcores=numcores)

We then simply provide either object (MNDOTheory or xTBTheory) as argument to the the *qm_theory* keyword in the QMMMTheory object.

.. code-block:: python

    #QM/MM object
    qm_mm = QMMMTheory(qm_theory=qm, mm_theory=omm, fragment=frag,
            qmatoms=qmatoms, printlevel=2, qm_charge=-2, qm_mult=1)

where we defined the QM-region to look like:

.. code-block:: python

    #QM-region
    GLU35=list(range(539, 547+1))
    ASP52=list(range(778,783+1))
    GLN57=list(range(856,866+1))
    qmatoms = GLU35+ASP52+GLN57


We note that while the MM Hamiltonian is still calculated with periodic boundary conditions (using either CPU or GPU), 
the QM-Hamiltonian is calculated here without periodic boundary conditions. This is an approximation which is reliable as long as the QM-region is approximately in the center of the box.

The results are shown in table below.

===================  ================================ ================================
QM-method             Time (sec) (for 10 ps)           Ave. time (sec) per timestep
===================  ================================ ================================
 MNDO-PM3                        4650                               0.465
 ORCA-PM3                        8949                               0.895
 OM2 or OM3                         X                                  X
 ODM2 or ODM3                         X                                  X
 GFN0-xTB (1 CPU)             X                                  0.311
 GFN1-xTB (1 CPU)             6551                               0.655
 GFN1-xTB (8 CPUs)            7194                                  X
 GFN2-XTB  (1 CPU)                    5454                                  0.545
 GFN2-XTB  (8 CPU)                    X                                  X
===================  ================================ ================================

As we can see, going from pure MM to the fastest semi-empirical QM/MM results in a major drop in speed of simulation (~800 times slower). 
The speed of the MNDO program and the PM3 semiempirical method is still quite impressive with each QM/MM step 
here taking only 0.465 s (465 ms) (really fast for a QM method!).
Each MNDO step takes approx. 376 ms, the rest is data exchange + MM-step (180 ms). 
As we will see, this is obviously much faster than what most regular QM-methods (e.g. DFT) will take.

The speed of the GFN-xTB methods is also very good, with each step taking less than a second.

TODO: xTB on GPU ??

xTB via CP2K??

Sparrow?

.. note:: xTB in ASH has 2 runmodes : runmode='library' and runmode='inputfile'. The runmode='inputfile' creates input on disk and calls the xtb binary while runmode='library' calls a Python-C-API and runs calculations as a library. The latter involves no I/O which can make calculations quicker. Unfortunately, pointcharge-embeded xTB calculations are not yet available with this option.

################################################################################
3. Lysozyme: QM/MM MD using non-hybrid DFT and composite methods
################################################################################

It should be clear by now that the speed of a QM/MM MD simulation is dominated by the QM-step and as we go beyond semi-empirical methods to DFT and others, the possible obtainable speed will drop further.
The data-exchange and communication between QM-program, MM-program and ASH is usually negligible. Only for the superfast MNDO-AM1 case does it reach ~XX % of the total step-time while it is ~X % for xTB. 
This essentially means that we can freely choose whatever QM-program we want that has an interesting electronic structure method implemented (in terms of speed or accuracy), as long as ASH features an interface to that program.

The disadvantage of semiempirical methods in the MNDO and xTB programs is the limited accuracy and they typically struggle with e.g. transition metal chemistry, open-shell electronic structure etc. even when parameters are available for all elements.
We now explore how the performance changes as we go to composite methods and non-hybrid DFT methods implemented in other QM programs.

We first consider the "3c composite" methods from the Grimme research group (HF-3c, B97-c, r2SCAN-3c and PBEh-3c).
These are HF or DFT methods with an accompanying small-ish tailored basis set and featuring some corrections (dispersion, basis set correction). 
These methods are designed to be reasonably fast (enabled by a small tailored basis set) yet reasonably accurate (achieved by the corrections).
These methods have perhaps a minor semi-empirical component, however, they do not make any approximations to the 2-electron integrals (like MNDO and xTB), instead utilize smaller basis sets (at least smaller than a well-polarized triple-zeta basis set) 
so there are fewer integrals to calculate.
These methods can thus be expected to be more accurate and robust but also considerably more expensive than traditional semiempirical methods.
All the "3c composite" methods are available in the ORCA program and we will use the ORCA interface in ASH to run these simulations.

.. code-block:: python
    
    #QM object example
    qm = ORCATheory(orcasimpleinput="! HF-3c tightscf", numcores=numcores)

=========================  ================================ ================================
QM-method                       Time (sec) (for 10 ps)        Ave. time (sec) per timestep
=========================  ================================ ================================
 ORCA5-HF-3c   (10 CPU)                 49848                               4.98
 ORCA5-B97-3c (10 CPU)                  164402                              16.4
 ORCA5-r2SCAN-3c (10 CPU)                185043                               18.5
 ORCA5-PBEh-3c (10 CPU)                  170138                                17.1
=========================  ================================ ================================

The results show that there is a considerably time-penalty per timestep now when running these "low-cost" QM calculations.
Using 10 CPU cores and the cheapest HF-3c method has each timestep taking 5 seconds and the cheapest DFT-method taking around 16 seconds.
HF-3c is thus an interesting choice if it is accurate enough for the system (closed-shell organic molecules primarily) while the DFT methods are likely necessary for sufficient accuracy of transition metal systems.
B97-3c is here the cheapest option. Interestingly the PBEh-3c method is faster than the r2SCAN-3c method due to smaller DZ basis set that is defined as part of this method (compared to the TZ basis sets in B97-3c and r2SCAN-3c).

If the "3c composite" methods are not a good option for some reason we can of course also utilize regular DFT methods, 
perhaps employing a DFT protocol that is known to be a good compromise between speed and accuracy (often system-dependent).
DFT programs come in many flavors and can roughly be grouped into programs that utilize local Gaussian-basis sets (e.g. ORCA, Gaussian, NWChem, Q-Chem, Molpro, Psi4 etc.) 
and programs that utilize plane-wave basis sets (e.g. CP2K, VASP, Quantum Espresso etc.).
The advantage of ASH is that we can in principle choose any of these programs and run them in a QM/MM simulation.
ASH at the moment has more interfaces to Gaussian-based DFT-codes (ORCA, pySCF, NWChem, Psi4, Gaussian) than to plane-wave DFT-codes (e.g. CP2K, VASP, Quantum Espresso).
Here we compare 2 diffent QM-codes: a state-of-the-art Gaussian-basis DFT code (ORCA) and a mixed planewave-Gaussian code (CP2K).
For ORCA we compare a GGA-functional (PBE) and the newer meta-GGA functionals TPSS and r2SCAN.

===========================================  ================================ ================================
QM-method                                        Time (sec) (for 10 ps)        Ave. time (sec) per timestep
===========================================  ================================ ================================
 ORCA-PBE/def2-SVP   (10 CPU)                       108770.294                      10.9
 ORCA-PBE/def2-SV(P) (10 CPU)                       98097.007                      9.8
 ORCA-TPSS/def2-SVP (10 CPU)                        123286.881                      12.3
 ORCA-r2SCAN/def2-SVP (10 CPU)                      124380.528                      12.4
===========================================  ================================ ================================

As can be seen, PBE, being a regular GGA functional runs marginally faster than TPSS and r2SCAN which are meta-GGA functionals.
These calculations take less time than before because we utilize the small def2-SVP basis set in the case of ORCA which is smaller than the tailored triple-zeta basis set used in the B97-3c and r2SCAN-3c methods.
A double-zeta basis set (like def2-SVP) is typically the minimum basis set that can be used for reliable DFT calculations but the accuracy of such small basis sets should be checked for each system.
Reducing the basis set size down to def2-SV(P) instead (fewer polarization functions) results in only a tiny speed-up here which may or may not be worth it.

Parallelization of the DFT calculation is an important aspect of running these types of calculations. More and more CPU cores does not necessarily mean much faster performance and can even in some cases be detrimental to performance.
The parallel performance depends strongly on the QM-code in question.
The benchmarking below reveals that parallelizing the ORCA-PBE calculation results in decent speed-up up to about 8-10 cores and then we get diminishing returns.

=========================  ================================ ================================
QM-method                       Time (sec) (for 10 ps)        Ave. time (sec) per timestep
=========================  ================================ ================================
  ORCA-PBE   (1 CPU)                 (not completed)                       63.2
  ORCA-PBE   (4 CPU)                 201493.058                            20.1
  ORCA-PBE   (8 CPU)                  124299.836                          12.4
  ORCA-PBE   (10 CPU)                  108770.294                        10.9
  ORCA-PBE   (16 CPU)                  88595.223                         8.9
  ORCA-PBE   (20 CPU)                  82957.761                         8.2
  ORCA-PBE   (24 CPU)                  79830.670                         7.9
=========================  ================================ ================================


CP2K is a program designed for fast DFT calculations, especially DFT MD simulations. 
It uses a hybrid Gaussian-planewaves approaches unlike the purely local atom-centered Gaussian basis-function approach of ORCA.
It is designed to scale up to a large number of CPU cores.


===========================================  ================================ ================================
QM-method                                        Time (sec) (for 10 ps)        Ave. time (sec) per timestep
===========================================  ================================ ================================
 CP2K-PBE/GPW (1 CPU,MPI)                                   (not completed)                   32.996 
 CP2K-PBE/GPW (4 CPU,MPI)                                   (not completed)                   12.451 
 CP2K-PBE/GPW (8 CPU,MPI)                                   (not completed)                   XX 
 CP2K-PBE/GPW (10 CPU,MPI)                                    81790.886                         8.2 
 CP2K-PBE/GPW (10 CPU,MPI)                                    81790.886                         8.2 
 CP2K-PBE/GPW (12 CPU,MPI)                                    76127.111                         7.6
 CP2K-PBE/GPW (24 CPU,MPI)                                    70305.023                         7.0
 CP2K-PBE/GPW (12 CPU,OMP)                                   93554.211                        9.3
 CP2K-PBE/GAPW (12 CPU,MPI)                             85668.633                               8.5
===========================================  ================================ ================================



################################################################################
4. Lysozyme: QM/MM MD using a special tailored minimal basis set
################################################################################

https://github.com/grimme-lab/qvSZP?tab=readme-ov-file
https://pubs.aip.org/aip/jcp/article/159/16/164108/2918302/An-atom-in-molecule-adaptive-polarized-valence



################################################################################
5. Lysozyme: QM/MM MD using hybrid-DFT
################################################################################

The HF exchange integrals in hybrid-DFT calcualtions typically dominate the cost of a hybrid-DFT calculation and this makes
hybrid-DFT ill-suited for dynamics studies as each timestep simply will be too expensive too compute.

However, hybrid-DFT is nevertheless typically the more accurate flavor of DFT (with exceptions) and for some systems 
hybrid-DFT may actually be necessary for even a qualitatively correct description of the system.
We here discuss options for running efficient hybrid-DFT QM/MM simulations.

ORCA RIJONX 
ORCA RIJK 
ORCA-RIJCOSX, 
CP2K
Turbomole ?

################################################################################
6. Lysozyme: QM/MM MD using GPU-based DFT-programs
################################################################################

NOTE: CP2K on GPU

While QM-programs are primarily written for the CPU, a few QM programs feature GPU implementation of HF/DFT. 
We will explore the performance of these here.
TeraChem, QUICK, pyscf-GPU, CP2K-GPU.


################################################################################
7. Lysozyme: QM/MM MD using WFT methods
################################################################################

Typically QM/MM MD simulations are limited to semi-empirical or DFT-based Hamiltonians.
MD simulations based on correlated wavefunction methods are typically too expensive and often lack gradients.

We will run the cheapest correlated WF method, MP2, as implemented in ORCA for comparison.


################################################################################
8. Lysozyme: Machine-learning potentials
################################################################################

################################################################################
9. Lysozyme: Truncated PC gradient approximation
################################################################################


################################################################################
10. Lysozyme: QM/MM MD SUMMARY
################################################################################

We have now gone through a variety of ways to perform QM/MM MD in ASH. As should be clear, ASH offers obvious advantages to performing these kinds
of studies as we can utilize a common system setup, common QM/MM coupling and a common dynamics propagation to run the calculations.
We only have to change the QM-program component to whatever method fits the particular study or what program is available.

The benchmarking numbers on this page should of course be taken with a pinch of salt as the relative performance of methods and programs
heavily depends on the nature of the system e.g. closed-shell vs. open-shell, the size of QM-region, computer hardware, how the QM-program was compiled/installed etc.
The speed of the QM-method/QM-program is also not the only variable to consider in a QM/MM MD study.

Nonetheless, the numbers here reveal some generally useful trends that may be worth considering when planning your next QM/MM MD study.
For example, a semiempirical protocol (using MNDO or xTB) is probably worth considering if your QM-region only contains "organic" elements
as the computational cost is much lower than for DFT methods.

TODO

Comparison of the most interesting methods from the previous sections:

=========================  ================================ ================================
QM-method                       Time (sec) (for 10 ps)        Ave. time (sec) per timestep
=========================  ================================ ================================
 XX   (1 CPU)                       X                               X

=========================  ================================ ================================

In this tutorial we have until now not considered accuracy.
TODO