Workflows
======================================

ASH includes a number of convenient built-in functionality to carry out multi-step workflows. See below.
See also :doc:`workflows-examples` for a tutorial on how you can build your own workflows. 


#####################
ReactionEnergy
#####################

.. code-block:: python

	def ReactionEnergy(list_of_energies=None, stoichiometry=None, list_of_fragments=None, unit='kcal/mol', label=None, reference=None, silent=False):

Options:

- list_of_energies (list of floats). List of total energies in Eh (hartrees).
- stoichiometry (list of integers). Integers that determine the stoichiometry of the reaction. Order must match list_of_energies or list_of_fragments.
- list_of_fragments (list of ASH Fragments). ASH fragments must have a set energy attribute.
- unit (string). String that indicates the final unit to convert reaction energy to. Options: 'kcal/mol', 'eV', 'kJ/mol', 'cm-1'. Default: 'kcal/mol'
- reference (float). If set, this signifies the reference energy and ReactionEnergy will print both energy and error w.r.t. reference. Default=None
- silent (Boolean). If True, ReactionEnergy will not print to standard-output, only return the result.

The simple ReactionEnergy function is a convenient way to calculate the reaction energy for a reaction from a list of energies and the stoichiometry associated with the reaction.
The function prints to standard output the reaction energy (unless silent=True) and returns the relative energy converted into a unit of choice (default: kcal/mol).

Simple example for Haber-Bosch reaction:  N\ :sub:`2` \  + 3H\ :sub:`2`\  → 2NH\ :sub:`3`\

.. code-block:: python

	from ash import *

	#Haber-Bosch reaction: N2 + 3H2 => 2NH3
	N2=Fragment(diatomic="N2", diatomic_bondlength=1.0975, charge=0, mult=1)
	H2=Fragment(diatomic="H2", diatomic_bondlength=0.741, charge=0, mult=1)
	NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)
	specieslist=[N2, H2, NH3] #An ordered list of ASH fragments.
	stoichiometry=[-1, -3, 2] #Using same order as specieslist.
	xtbcalc=xTBTheory(xtbmethod='GFN1') # GFN1-xTB theory-level
	energies = Singlepoint_fragments(theory=xtbcalc, fragments=specieslist) #Calculating list of energies

	#Calculating reaction-energy using list and stoichiometry
	reaction_energy, unused = ReactionEnergy(stoichiometry=stoichiometry, list_of_energies=energies, unit='kcal/mol', label='ΔE')

.. code-block:: text

	Reaction_energy(Δ):  -136.6723479900558 kcal/mol


If there is an energy attribute associated with each fragment it is also possible to just provide ReactionEnergy with a list of the fragments involved.
This will only work if the energy attribute of the fragment has been definded. Some ASH functions will do this: **Singlepoint**, **Singlepoint_fragments**, **geomeTRICOptimizer**

.. code-block:: python

	#Calculating reaction-energy using list_of_fragments and stoichiometry
	specieslist=[N2, H2, NH3]
	reaction_energy, unused = ReactionEnergy(stoichiometry=stoichiometry, list_of_fragments=specieslist, unit='kcal/mol', label='ΔE')

#####################
Thermochemprotocols
#####################


The **thermochemprotocol_reaction** and **thermochemprotocol_single** functions can be used to
perform a multi-step Opt+Freq+HL-single-point protocol on either a reaction or a single species.


The thermochemprotocol_reaction is used for chemical reactions by giving a list of ASH fragments, stoichiometry and theory levels.

.. code-block:: python

    def thermochemprotocol_reaction(fraglist=None, stoichiometry=None, Opt_theory=None, SP_theory=None, orcadir=None, numcores=None, memory=5000,
                       analyticHessian=False, temp=298.15, pressure=1.0)

while thermochemprotocol_single is used for a single fragment (thermochemprotocol_reaction calls thermochemprotocol_single).

.. code-block:: python

    def thermochemprotocol_single(fragment=None, Opt_theory=None, SP_theory=None, orcadir=None, numcores=None, memory=5000,
                       analyticHessian=True, temp=298.15, pressure=1.0):


The reaction is defined via a list of defined fragments and stoichiometry, a theory object for Opt+Freq steps is defined (Opt_theory)
and then a theory for the high-level single-point level is chosen (SP_theory). Can be any ASH Theory including ORCATheory, CC_CBS_Theory etc.

.. code-block:: python

    from ash import *

    #
    orcadir='/opt/orca_5.0.2'
    numcores=4

    N2=Fragment(xyzfile="n2.xyz", charge=0, mult=1)
    H2=Fragment(xyzfile="h2.xyz", charge=0, mult=1)
    NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)

    # List of species from reactant to product
    specieslist=[N2, H2, NH3] #Use same order as stoichiometry
    #Equation stoichiometry : negative integer for reactant, positive integer for product
    # Example: N2 + 3H2 -> 2NH3  reaction should be:  [-1,-3,2]
    stoichiometry=[-1, -3, 2] #Use same order as specieslist

    #Opt+Freq theory
	B3LYP_opt=ORCATheory(orcasimpleinput="! B3LYP D3BJ def2-TZVP def2/J tightscf", numcores=numcores)
	#HL theory
	DLPNO_CC_calc = CC_CBS_Theory(elements=["N", "H"], cardinals = [2,3], basisfamily="def2", DLPNO=True, 
                  pnosetting='extrapolation', pnoextrapolation=[6,7], numcores=numcores)

	#Example: Thermochemistry protocol on the single N2 species
    thermochemprotocol_single(fragment=N2, stoichiometry=stoichiometry, orcadir=orcadir,
                        numcores=numcores, Opt_theory=None, SP_theory=DLPNO_CC_calc)

    #Alternative: Thermochemistry protocol on the whole N2 + 3 H2 => 2 NH3 reaction
    thermochemprotocol_reaction(fraglist=specieslist, stoichiometry=stoichiometry, orcadir=orcadir,
                        numcores=numcores, Opt_theory=B3LYP_opt, SP_theory=DLPNO_CC_calc)

###############################################################
calc_xyzfiles: Run calculations on a collection of XYZ-files
###############################################################

calc_xyzfiles is similar to Singlepoint_fragments (:doc:`singlepoint`) but saves you the step of defining fragments manually if you already have XYZ-files collected in a directory.


.. code-block:: python

	def calc_xyzfiles(xyzdir=None, theory=None, Opt=False, Freq=False, charge=None, mult=None, xtb_preopt=False):


If you have a collection of XYZ-files that you wish to run calculations on (either single-point energy evalutation or geometry optimizations) 
then this can be easily accomplished using the calc_xyzfiles function. 
Charge and multiplicities for each XYZ-file need to be given in the description-line (2nd line) of each XYZ-file like this:

HCl.xyz example:

.. code-block:: text

	2
	0 1
	H 0.0 0.0 0.0
	Cl 0.0 0.0 0.0

Alternatively, if all molecules are e.g. neutral singlets then one can give charge=0, mult=1 keyword arguments to **calc_xyzfiles()**

Example script:

.. code-block:: python

	from ash import *

	numcores=24
	#Directory of XYZ files. Can be full path or relative path (dir needs to be copied to scratch location in this case).
	dir = '/home/bjornsson/FeCO4_N2/r2scan-opt/xyzfiles_temp'

	#Defining theory. Charge/mult is skipped here
	ORCAcalc = ORCATheory(orcasimpleinput="! r2SCAN-3c", orcablocks="%scf maxiter 500 end", numcores=numcores)

	#Call calc_xyzfiles giving xyzdir and theory. 
	#Geometry optimizations for each XYZ-file can be requested via Opt=True (default False, i.e. singlepoint) 
	calc_xyzfiles(xyzdir=dir, theory=ORCAcalc, Opt=True)

	# Same but with an xTB pre-optimization (requires xtb to be installed)
	#calc_xyzfiles(xyzdir=dir, theory=ORCAcalc, Opt=True, xtb_preopt=True)



The ASH script then runs through and gives a table at the end with the energies. 
In the case of Opt=True, a geometry optimization is performed for each molecule and a final directory of XYZ-files with optimized coordinates is created.


.. code-block:: text

	XYZ-file             Charge     Mult           Energy(Eh)
	----------------------------------------------------------------------
	no.xyz                     0       2      -129.8755914784
	no_plus.xyz                1       1      -129.5232460574
	h2.xyz                     0       1        -1.1693816161
	n2.xyz                     0       1      -109.5070757384
	hbr.xyz                    0       1     -2574.7361724856


	XYZ-files with optimized coordinates can be found in: optimized_xyzfiles


###################################
High-level single-point workflows
###################################

See :doc:`module_highlevel_workflows`

#######################################################################
confsampler_protocol : Automatic Crest+DFTopt+DLPNO-CCSD(T) workflow
#######################################################################

See :doc:`crest-interface`



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

	activespace_dictionary = auto_active_space(fragment=fragment, orcadir=None, basis="def2-TZVP", charge=0, mult=1,
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




