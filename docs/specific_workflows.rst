Specific workflows
======================================

This page contains information on very specific workflow functionality for carrying out specific things such as:

- counterpoise corrections
- finding non-Aufbau SCF solutions
- finding an active space automatically
- calculating redox density differences

###################################
Counter-poise correction (ORCA)
###################################

.. code-block:: python
	
	def counterpoise_calculation_ORCA(fragments=None, theory=None, monomer1_indices=None, monomer2_indices=None):

ASH can perform Boys-Bernardi counterpoise corrections (single-point energy level only) together with ORCA in a convenient way.
All that is required are geometries (previously optimized) for the AB dimer as well as monomers A and B respectively, a theory level definition and lists of atom indices that specify which atoms in the AB dimer belong to monomer A and B, respectively. 

.. code-block:: python

	from ash import *

	#Define ASH fragments for the A-B adduct (dimer) and monomers from XYZ-files
	#Dimer: H2O...MeOH H-bonded complex
	dimer=Fragment(xyzfile="h2o_meoh.xyz", charge=0, mult=1)
	#H2O monomer
	h2o=Fragment(xyzfile="h2o.xyz", charge=0, mult=1)
	#MeOH monomer
	meoh=Fragment(xyzfile="meoh.xyz", charge=0, mult=1)
	#Combine fragments in a list
	all_fragments=[dimer, h2o, meoh]

	#Define ORCA theory
	simple=" ! RI-MP2 def2-SVP def2-SVP/C RIJCOSX def2/J tightscf "
	blocks="""
	%scf
	maxiter 300
	end
	"""
	orcacalc = ORCATheory(orcasimpleinput=simple, orcablocks=blocks)


	#Run counterpoise_calculation giving fragment-list, orcacalculation and atom-indices as input
	# monomer1_indices and monomer2_indices specify which atoms in the dimer correspond to monomer1 and monomer2
	counterpoise_calculation_ORCA(fragments=all_fragments, theory=orcacalc, monomer1_indices=[0,1,2], monomer2_indices=[3,4,5,6,7,8])


The final output looks like :


.. code-block:: text

	                #######################################
	                #                                     #
	              #     COUNTERPOISE CORRECTION JOB     #
	                #                                     #
	                #######################################



	 Boys-Bernardi counterpoise correction

	monomer1_indices: [0, 1, 2]
	monomer2_indices: [3, 4, 5, 6, 7, 8]

	Monomer 1:
	--------------------
	Defined coordinates (Å):
	   O  -0.52532979   -0.05097108   -0.31451686
	   H  -0.94200663    0.74790163    0.01125282
	   H   0.40369652    0.05978598   -0.07356837
	Monomer 1 indices in dimer: [0, 1, 2]

	Monomer 2:
	--------------------
	Defined coordinates (Å):
	   O   2.31663329    0.04550085    0.07185839
	   H   2.68461611   -0.52657655    0.74938672
	   C   2.78163836   -0.42612907   -1.19030072
	   H   2.35082127    0.22496462   -1.94341475
	   H   3.86760205   -0.37533621   -1.26461265
	   H   2.45329574   -1.44599856   -1.38938136
	Monomer 2 indices in dimer: [3, 4, 5, 6, 7, 8]

	Dimer:
	--------------------
	0   O -0.525329794 -0.050971084 -0.314516861   Monomer1
	1   H -0.942006633 0.747901631 0.011252816   Monomer1
	2   H 0.403696525 0.059785981 -0.073568368   Monomer1
	3   O 2.316633291 0.045500849 0.071858389   Monomer2
	4   H 2.684616115 -0.526576554 0.749386716   Monomer2
	5   C 2.781638362 -0.426129067 -1.190300721   Monomer2
	6   H 2.350821267 0.224964624 -1.943414753   Monomer2
	7   H 3.867602049 -0.375336206 -1.264612649   Monomer2
	8   H 2.453295744 -1.445998564 -1.389381355   Monomer2


	----LOTS OF CALCULATION OUTPUT---

	COUNTERPOISE CORRECTION RESULTS
	==================================================

	Monomer 1 energy: -76.162192724532 Eh
	Monomer 2 energy: -115.290878785879 Eh
	Sum of monomers energy: -191.453071510411 Eh
	Dimer energy: -191.465349252819 Eh

	Monomer 1 at dimer geometry: -115.290878793717 Eh
	Monomer 2 at dimer geometry: -76.162192727048 Eh
	Sum of monomers at dimer geometry energy: -191.45307152076498 Eh

	Monomer 1 at dimer geometry with dimer basis: -115.29491810198 Eh
	Monomer 2 at dimer geometry with dimer basis: -76.163483336908 Eh
	Sum of monomers at dimer geometry with dimer basis: -191.45840143888802 Eh
	counterpoise_corr: 3.344574118169517 kcal/mol

	Uncorrected interaction energy: -7.704399681128008 kcal/mol
	Corrected interaction energy: -4.359825562958491 kcal/mol



###################################
Automatic non-Aufbau calculator
###################################

Excited SCF configurations can be tricky to converge to without falling back to the ground-state. 
While various different algorithms have recently been suggested in the literature to help locating such excited SCF configurations, often the methods are  only available in specific QM codes.
The `STEP <https://doi.org/10.1021/acs.jctc.0c00502>`_ algorithm by Carter-Fenk and Herbert is a much simpler algorithm and can be used with any QM program with level-shifting implemented (a common SCF convergence aid).
The idea is simply to choose to change the MO occupation as desired (e.g. swap the HOMO and LUMO orbitals) and then choose a specific levelshift and hopefully converge to the desired SCF configuration.
The levelshift is chosen based on the occupied-virtual orbital-energy gap together with an extra parameter, epsilon (0.1 by default).

ASH allows one to utilize the STEP algorithm in a convenient way together with ORCA (only QM-program supported so far) using the **AutoNonAufbau** function.

.. code-block:: python

	def AutoNonAufbau(fragment=None, theory=None, num_occ_orbs=1, num_virt_orbs=3, spinset=[0], stability_analysis_GS=False, 
					TDDFT=False, epsilon=0.1, maxiter=500, manual_levelshift=None):


**AutoNonAufbau** options:

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
     - ASH Fragment object.
   * - ``theory``
     - ASH Theory
     - None
     - An ASH theory level. Currently only ORCATheory is supported.
   * - ``num_occ_orbs``
     - integer
     - 1
     - Number of occupied orbitals to include the in orbital rotation procedure. A value of 3 would include the HOMO, HOMO-1 and HOMO-2
   * - ``num_virt_orbs``
     - integer
     - 3
     - Number of virtual orbitals to include the in orbital rotation procedure. A value of 3 would include the LUMO, LUMO+1 and LUMO+2
   * - ``spinset``
     - list
     - [0]
     - What spin manifold to use. Alpha: [0] or  Beta: [1] or Both: [0,1]
   * - ``TDDFT``
     - Boolean
     - False
     - Whether to do TDDFT on the ground-state in order to select orbitals to rotate. Experimental feature
   * - ``epsilon``
     - float
     - 0.1
     - The value of the epsilon parameter in the STEP algorithm. 
   * - ``maxiter``
     - integer
     - 500
     - Maximum number of ORCA SCF iterations.
   * - ``manual_levelshift``
     - float
     - None
     - Manual levelshift instead of the automatic levelshift (based on the gap and epsilon parameter)


**How to use:**

One reads in an ASH fragment, an ORCATheory object and optionally specifies how many occupied and virtual orbitals, what spin manifold to use etc.
ASH will then tell ORCA to calculate the ground-state SCF in a first step.
ASH will then go through all possible SCF configurations that involve the highest-energy occupied MOs and the lowest-energy virtual MOs based on the user selection,
rotate the respective occupied and virtual orbital pair (e.g. the HOMO and LUMO+1), apply a levelshift based on the orbital-energy gap and the epsilon parameter (0.1 Eh by default)
and then attempt to converge the SCF for each possible guess configuration.


Example: Finding excited SCF states of the water molecule.

.. code-block:: python

	from ash import *
	#
	frag=Fragment(databasefile="h2o.xyz", charge=0, mult=1)

	orcacalc=ORCATheory(orcasimpleinput="! UHF def2-QZVPPD usesym")

	AutoNonAufbau(theory=orcacalc, fragment=frag, num_occ_orbs=2, num_virt_orbs=16, spinset=[0])

**Output:**

.. code-block:: text

						#########################
						#                       #
						#     AutoNonAufbau     #
						#                       #
						#########################


	Spin orbital sets to choose: [0]
	Number of occupied orbitals allowed in MO swap: 2
	Number of virtual orbitals allowed in MO swap: 16
	Total number of states: 32
	TDDFT: False
	stability_analysis_GS: False
	Epsilon: 0.1
	manual_levelshift: None
	Cleaning up old ORCA files
	Now doing initial state SCF calculation
	...
	Energy:  -76.06701160263
	...
	==========================================================================================
	Now running excited state SCF no. 0 with multiplicity: 1 and spinset 0
	==========================================================================================
	Simple MO selection scheme
	homo_lumo_gap: -0.5747076804099635
	Will rotate orbital 4 (HOMO) and orbital 5
	lshift: 0.6747076804099634
	Now doing SCF calculation with rotated MOs and levelshift
	...
	Energy:  -75.835914096332
	GS/ES state energy gap: 6.29 eV
	Found something different than ground state
	Converged SCF energy higher than ground-state SCF. Found new excited state SCF solution !

	==========================================================================================
	Now running excited state SCF no. 1 with multiplicity: 1 and spinset 0
	==========================================================================================
	Simple MO selection scheme
	homo_lumo_gap: -0.58996231689521
	Will rotate orbital 4 (HOMO) and orbital 6
	lshift: 0.6899623168952099
	Now doing SCF calculation with rotated MOs and levelshift
	...
	Energy:  -75.771691604422
	GS/ES state energy gap: 8.04 eV
	Found something different than ground state
	Converged SCF energy higher than ground-state SCF. Found new excited state SCF solution !
	==========================================================================================
	Now running excited state SCF no. 2 with multiplicity: 1 and spinset 0
	==========================================================================================
	Simple MO selection scheme
	homo_lumo_gap: -0.6157346044574922
	Will rotate orbital 4 (HOMO) and orbital 7
	lshift: 0.7157346044574922
	Now doing SCF calculation with rotated MOs and levelshift
	...
	Energy:  -76.067008807825
	GS/ES state energy gap: 0.00 eV
	GS/ES state energy gap smaller than 0.04 eV. Presumably found the original SCF again.


The final output will print a list of all the states found.

.. code-block:: text

	Ground-state SCF energy: -76.06701160263 Eh
	-----
	Excited state index: 0
	Spin multiplicity: 1
	Excited-state SCF energy 0: -75.835914096332 Eh
	Excited-state SCF transition energy: 6.29 eV
	Excited-state SCF orbital rotation: Occ:[4] Virt:[5]
	Rotation in spin manifold:  [0]
	Excited-state SCF HOMO-LUMO gap: -0.5747076804099635 Eh (-15.638599999999999) eV:
	Excited-state SCF Levelshift chosen: 0.6747076804099634
	-----
	Excited state index: 1
	Spin multiplicity: 1
	Excited-state SCF energy 1: -75.771691604422 Eh
	Excited-state SCF transition energy: 8.04 eV
	Excited-state SCF orbital rotation: Occ:[4] Virt:[6]
	Rotation in spin manifold:  [0]
	Excited-state SCF HOMO-LUMO gap: -0.58996231689521 Eh (-16.0537) eV:
	Excited-state SCF Levelshift chosen: 0.6899623168952099
	-----
	Excited state index: 2
	Spin multiplicity: 1
	Excited-state SCF energy 2: -76.067008807825 Eh
	Excited-state SCF transition energy: 0.00 eV
	Excited-state SCF orbital rotation: Occ:[4] Virt:[7]
	Rotation in spin manifold:  [0]
	Excited-state SCF HOMO-LUMO gap: -0.6157346044574922 Eh (-16.755) eV:
	Excited-state SCF Levelshift chosen: 0.7157346044574922

ORCA outputfiles and GBW files for each state is kept and can be inspected or used for future calculations.

- orca_GS.out
- orca_GS.gbw
- orcaES_SCF0_mult1_spinset0.out
- orcaES_SCF0_mult1_spinset0.gbw
- orcaES_SCF1_mult1_spinset0.out
- orcaES_SCF1_mult1_spinset0.gbw
- orcaES_SCF2_mult1_spinset0.out
- orcaES_SCF2_mult1_spinset0.gbw



###################################
Automatic active-space selection
###################################

.. code-block:: python

	def auto_active_space(fragment=None, orcadir=None, basis="def2-SVP", scalar_rel=None, charge=None, mult=None, 
    initial_orbitals='MP2', functional='TPSS', smeartemp=5000, tgen=1e-1, selection_thresholds=[1.999,0.001],
    numcores=1):

Workflow to guess a good active space for CASSCF calculation based on a 2-step procedure:
1. Calculate MP2-natural orbitals (alternative Fractional occupation DFT orbitals)
2. ICE-CI on top of MP2-natural orbitals using a large active-space but with small tgen threshold


Example on ozone:

.. code-block:: python

	from ash import *

	fragstring="""
	O       -2.219508975      0.000000000     -0.605320629
	O       -1.305999766     -0.913250049     -0.557466332
	O       -2.829559171      0.140210894     -1.736132689
	"""

	fragment=Fragment(coordsstring=fragstrin, charge=0, mult=1)

	activespace_dictionary = auto_active_space(fragment=fragment, basis="def2-TZVP", charge=0, mult=1,
	    initial_orbitals='MP2', tgen=1.0)
	#Returns dictionary with various active_spaces based on thresholds

Output:

.. code-block:: text

	ICE-CI step done
	Note: New natural orbitals from ICE-CI density matrix formed!

	Wavefunction size:
	Tgen: 1.0
	Tvar: 1e-07
	Orbital space of CAS(18,37) used for ICE-CI step
	Num generator CFGs: 4370
	Num CFGS after S+D: 4370

	Table of natural occupation numbers

	Orbital   MP2natorbs ICE-nat-occ
	----------------------------------------
	0            2.0000    2.0000
	1            2.0000    2.0000
	2            2.0000    2.0000
	3            1.9859    1.9898
	4            1.9809    1.9869
	5            1.9747    1.9836
	6            1.9637    1.9791
	7            1.9607    1.9787
	8            1.9360    1.9665
	9            1.9223    1.9631
	10           1.9197    1.9603
	11           1.8522    1.9371
	12           0.1868    0.0779
	13           0.0680    0.0349
	14           0.0612    0.0318
	15           0.0241    0.0122
	16           0.0171    0.0093
	17           0.0146    0.0081
	18           0.0117    0.0076
	19           0.0106    0.0067
	20           0.0105    0.0064
	...

	Recommended active spaces based on ICE-CI natural occupations:
	Minimal (1.95,0.05): CAS(2,2)
	Medium1 (1.98,0.02): CAS(12,9)
	Medium2 (1.985,0.015): CAS(14,10)
	Medium3 (1.99,0.01): CAS(18,13)
	Medium4 (1.992,0.008): CAS(18,15)
	Large (1.995,0.005): CAS(18,19)
	Orbital file to use for future calculations: orca.gbw
	Note: orbitals are new natural orbitals formed from the ICE-CI density matrix

