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



###################################
Redox difference density
###################################

.. image:: figures/cocene-redox-diffdens-300.png
   :align: center
   :width: 700

To understand the nature of a redox process it can be useful to calculate and visualize the density difference between the two redox species.
This can only cleanly be done for a vertical redox process.
First one requires Cube files of the electron density for both species that should have the same geometry and same number of grid points (and ideally quite fine).
Then one can call the built-in ASH function: **read_cube** and **write_cube_diff**  to output the density difference as a Cube file.
The calculations, Cube-file creation and the density-difference can be done in a single script as shown below if using ORCA (orca_plot used to generate the Cubefiles.)

**Example: Vertical ionization of Cobaltocene**

.. code-block:: python

	from ash import *
	import shutil

	string="""
	Co       6.344947000     -1.560817000      5.954256000
	C        6.026452000     -0.546182000      7.802276000
	C        5.793563000     -1.965872000      7.908267000
	C        7.027412000     -2.657637000      7.660083000
	C        7.976024000     -1.681400000      7.254253000
	C        7.360141000     -0.371346000      7.359546000
	H        5.287260000      0.234934000      7.955977000
	H        4.853677000     -2.430623000      8.196105000
	H        7.174908000     -3.733662000      7.680248000
	H        7.845323000      0.573837000      7.130600000
	C        7.003675000     -1.909758000      4.002507000
	C        6.025582000     -2.892708000      4.310293000
	C        4.831240000     -2.191536000      4.690416000
	C        5.029751000     -0.780020000      4.468693000
	C        6.380790000     -0.601338000      4.083942000
	H        8.038293000     -2.097559000      3.727053000
	H        6.179063000     -3.967316000      4.350882000
	H        3.905308000     -2.652860000      5.025216000
	H        6.875603000      0.346060000      3.887076000
	H        4.297348000      0.002948000      4.644348000
	H        8.999375000     -1.871686000      6.940915000
	"""

	#Defining fragment for redox reaction
	Co_neut=Fragment(coordsstring=string, charge=0, mult=2)
	Co_ox=Fragment(coordsstring=string, charge=1, mult=1)
	label="Cocene_"+'_'
	#Defining QM theory as ORCA here
	qm=ORCATheory(orcasimpleinput="! BP86 def2-SVP tightscf notrah")

	#Run neutral species with ORCA
	e_neut=Singlepoint(theory=qm, fragment=Co_neut)
	shutil.copyfile(qm.filename+'.gbw', label+"neut.gbw") # Copy GBW file
	#Run orca_plot to request electron density creation from ORCA gbw file
	run_orca_plot(label+"neut.gbw", "density", gridvalue=80)
	
	#Run oxidized species with ORCA
	e_ox=Singlepoint(theory=qm, fragment=Co_ox)
	shutil.copyfile(qm.filename+'.gbw', label+"ox.gbw")  # Copy GBW file
	#Run orca_plot to request electron density creation from ORCA gbw file
	run_orca_plot(label+"ox.gbw", "density", gridvalue=80)

	#Read Cubefiles from disk. 
	neut_cube_data = functions.functions_elstructure.read_cube(label+"neut.eldens.cube")
	ox_cube_data = functions.functions_elstructure.read_cube(label+"ox.eldens.cube")
	#Write out difference density as a Cubefile
	functions.functions_elstructure.write_cube_diff(neut_cube_data, ox_cube_data, label+"diffence_density.cube")

The script will output the files Cocene_neut.eldens.cube and Cocene_ox.eldens.cube that are here generated by orca_plot. 
The file Cocene_diffence_density.cube is generated by **write_cube_diff**.
