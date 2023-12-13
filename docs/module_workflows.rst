Workflow functionality
======================================

ASH includes a number of convenient built-in functionality to help the user carry out various multi-step workflows.

- The page :doc:`singlepoint` discusses various ways of performing singlepoint calculations with multiple fragments or multiple theories.
- The page :doc:`module_highlevel_workflows` introduces high-level multistep workflows for doing CCSD(T)/CBS calculations
- The page :doc:`workflows-examples` is a tutorial on how you can build your own workflows.
- The page :doc:`specific_workflows` includes various highly specific workflows implemented in ASH.
- The page :doc:`crest-interface` introduces an automatic conformation sampling + optimization + highlevel SP workflow.

This page documents the **Reaction** class, the **ReactionEnergy** function, **thermochemprotocol_single** and **thermochemprotocol_reaction** functions, 
and the **calc_xyzfiles** function 


#####################
Reaction class
#####################

The Reaction class is a simple convenient ASH object to define a molecular reaction object that can be used
as input in various multistep workflows.


.. code-block:: python

	class Reaction:
		def __init__(self, fragments, stoichiometry, label=None):

To create a Reaction object one first defines all molecular species of reaction as ASH fragments (with charge and multiplicity defined).
The ASH fragments are collected in a list and a stoichiometry is defined (same order as list of fragments):

**Example:**

.. code-block:: python

	N2=Fragment(xyzfile="n2.xyz", charge=0, mult=1)
	H2=Fragment(xyzfile="h2.xyz", charge=0, mult=1)
	NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)
	# List of species from reactant to product
	specieslist=[N2, H2, NH3] #Use same order as stoichiometry
	#Equation stoichiometry : negative integer for reactant, positive integer for product
	# Example: N2 + 3H2 -> 2NH3  reaction should be:  [-1,-3,2]
	stoichiometry=[-1, -3, 2] #Use same order as specieslist
	HB_reaction = Reaction(fragments=specieslist, stoichiometry=stoichiometry)
	#Theory object
	dft = ORCATheory(orcasimpleinput="! BP86 def2-SVP def2/J tightscf")
	#Running a Singlepoint_reaction job that supports a Reaction object as input
	Singlepoint_reaction(reaction=HB_reaction, theory=dft)

Job-types that can take a Reaction as inputobject:

- **Singlepoint_reaction** (See :doc:`singlepoint`)
- **thermochemprotocol_reaction** (see below)


#####################
ReactionEnergy
#####################

.. code-block:: python

	def ReactionEnergy(list_of_energies=None, stoichiometry=None, list_of_fragments=None, unit='kcal/mol', 
		label=None, reference=None, silent=False):


**ReactionEnergy** options:

.. list-table::
   :widths: 15 15 15 60
   :header-rows: 1

   * - Keyword
     - Type
     - Default value
     - Details
   * - ``list_of_energies``
     - List of floats
     - None
     - List of energies as floats for reaction. Order must match stoichiometry list.
   * - ``stoichiometry``
     - ASH Theory
     - list of integers
     - Integers for stoichiometry of the reaction. Order must match list_of_energies or list_of_fragments.
   * - ``list_of_fragments``
     - list of ASH Fragments
     - None
     - List of ASH fragments that must have a set energy attribute. Alternative to list_of_energies.
   * - ``unit``
     - string
     - 'kcal/mol'
     - String for final unit to convert reaction energy to. Options: 'kcal/mol', 'eV', 'kJ/mol', 'cm-1', 'Eh', 'mEh', 'meV'.
   * - ``reference``
     - float
     - None
     - If set, will print both energy and error w.r.t. reference.
   * - ``silent``
     - Boolean
     - False
     - Whether function prints to stdout (True) or not (False)



The simple ReactionEnergy function is a convenient way to calculate the reaction energy for a reaction from a list of energies and the stoichiometry associated with the reaction.
The function prints to standard output the reaction energy (unless silent=True) and returns the relative energy converted into a unit of choice (default: kcal/mol).

Simple example for Haber-Bosch reaction:  N\ :sub:`2` \  + 3H\ :sub:`2`\  → 2NH\ :sub:`3`\

.. code-block:: python

	from ash import *

	#Haber-Bosch reaction: N2 + 3H2 => 2NH3
	N2=Fragment(diatomic="N2", bondlength=1.0975, charge=0, mult=1)
	H2=Fragment(diatomic="H2", bondlength=0.741, charge=0, mult=1)
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
This will only work if the energy attribute of the fragment has been defined. Some ASH functions will do this: **Singlepoint**, **Singlepoint_fragments**, **geomeTRICOptimizer**

.. code-block:: python

	#Calculating reaction-energy using list_of_fragments and stoichiometry
	specieslist=[N2, H2, NH3]
	reaction_energy, unused = ReactionEnergy(stoichiometry=stoichiometry, list_of_fragments=specieslist, unit='kcal/mol', label='ΔE')

#####################
Thermochemprotocols
#####################


The **thermochemprotocol_reaction** and **thermochemprotocol_single** functions can be used to
perform a multi-step Opt+Freq+HL-single-point protocol on either a reaction or a single species.


The **thermochemprotocol_reaction** is used for chemical reactions by providing multiple theory level (for Opt+Freq and High-level singlepoint)
and an ASH Reaction object.

.. code-block:: python

	def thermochemprotocol_reaction(Opt_theory=None, SP_theory=None, reaction=None, fraglist=None, stoichiometry=None, numcores=1, memory=5000,
						analyticHessian=True, temp=298.15, pressure=1.0, unit='kcal/mol'):

while **thermochemprotocol_single** is used for a single fragment (**thermochemprotocol_reaction** calls **thermochemprotocol_single**).

.. code-block:: python

    def thermochemprotocol_single(fragment=None, Opt_theory=None, SP_theory=None, orcadir=None, numcores=None, memory=5000,
                       analyticHessian=True, temp=298.15, pressure=1.0):


The reaction must first be defined for a list of defined fragments and stoichiometry, a theory object for Opt+Freq steps is defined (Opt_theory)
and then a theory for the high-level single-point level is chosen (SP_theory). Can be any ASH Theory including ORCATheory, CC_CBS_Theory etc.

**thermochemprotocol_reaction example:**

.. code-block:: python

	from ash import *

	#
	numcores=4

	N2=Fragment(xyzfile="n2.xyz", charge=0, mult=1)
	H2=Fragment(xyzfile="h2.xyz", charge=0, mult=1)
	NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)

	# List of species from reactant to product
	specieslist=[N2, H2, NH3] #Use same order as stoichiometry
	#Equation stoichiometry : negative integer for reactant, positive integer for product
	# Example: N2 + 3H2 -> 2NH3  reaction should be:  [-1,-3,2]
	stoichiometry=[-1, -3, 2] #Use same order as specieslist
	#ASH reaction object
	HB_reaction = Reaction(fragments=specieslist, stoichiometry=stoichiometry)

	#Opt+Freq theory
	B3LYP_opt=ORCATheory(orcasimpleinput="! B3LYP D3BJ def2-TZVP def2/J tightscf", numcores=numcores)
	#HL theory
	DLPNO_CC_calc = ORCA_CC_CBS_Theory(elements=["N", "H"], cardinals = [2,3], basisfamily="def2", DLPNO=True, 
					pnosetting='extrapolation', pnoextrapolation=[6,7], numcores=numcores)
	#Alternative: Thermochemistry protocol on the whole N2 + 3 H2 => 2 NH3 reaction
	thermochemprotocol_reaction(fraglist=specieslist, stoichiometry=stoichiometry,
						numcores=numcores, Opt_theory=B3LYP_opt, SP_theory=DLPNO_CC_calc, unit='kcal/mol')


**thermochemprotocol_single example:**

.. code-block:: python

	from ash import *

	#Fragment
	N2=Fragment(xyzfile="n2.xyz", charge=0, mult=1)
	#Theories
	B3LYP_opt=ORCATheory(orcasimpleinput="! B3LYP D3BJ def2-TZVP def2/J tightscf", numcores=numcores)
	DLPNO_CC_calc = ORCA_CC_CBS_Theory(elements=["N", "H"], cardinals = [2,3], basisfamily="def2", DLPNO=True, 
					pnosetting='extrapolation', pnoextrapolation=[6,7], numcores=1)
	#Job
	thermochemprotocol_single(fragment=N2, Opt_theory=B3LYP_opt, SP_theory=DLPNO_CC_calc)

###############################################################
calc_xyzfiles: Run calculations on a collection of XYZ-files
###############################################################

**calc_xyzfiles** is similar to **Singlepoint_fragments** (:doc:`singlepoint`) but saves you the step of defining fragments manually if you already have XYZ-files collected in a directory.


.. code-block:: python

	def calc_xyzfiles(xyzdir=None, theory=None, Opt=False, Freq=False, charge=None, mult=None, xtb_preopt=False):


If you have a collection of XYZ-files that you wish to run calculations on (either single-point energy evalutation or geometry optimizations) 
then this can be easily accomplished using the **calc_xyzfiles** function. 
Charge and multiplicities for each XYZ-file need to be given in the description-line (2nd line) of each XYZ-file like this:

HCl.xyz example:

.. code-block:: text

	2
	0 1
	H 0.0 0.0 0.0
	Cl 0.0 0.0 1.3

Alternatively, if all molecules are e.g. neutral singlets then one can give charge=0, mult=1 keyword arguments to **calc_xyzfiles()**

Example script:

.. code-block:: python

	from ash import *

	numcores=24
	#Directory of XYZ files. Can be full path or relative path (dir needs to be copied to scratch location in this case).
	dir = '/home/bjornsson/FeCO4_N2/r2scan-opt/xyzfiles_temp'

	#Defining theory.
	ORCAcalc = ORCATheory(orcasimpleinput="! r2SCAN-3c", orcablocks="%scf maxiter 500 end", numcores=numcores)

	#Call calc_xyzfiles giving xyzdir and theory. 
	#Geometry optimizations for each XYZ-file can be requested via Opt=True (default False, i.e. singlepoint) 
	calc_xyzfiles(xyzdir=dir, theory=ORCAcalc, Opt=True)

	# Same but with an xTB pre-optimization (requires xtb to be installed)
	#calc_xyzfiles(xyzdir=dir, theory=ORCAcalc, Opt=True, xtb_preopt=True)



The ASH script then runs through and gives a table at the end with the energies. 
In the case of Opt=True, a geometry optimization is performed for each molecule at the chosen theory-level instead of a singlepoint calculations 
and a final directory of XYZ-files with optimized coordinates is created.


.. code-block:: text

	XYZ-file             Charge     Mult           Energy(Eh)
	----------------------------------------------------------------------
	no.xyz                     0       2      -129.8755914784
	no_plus.xyz                1       1      -129.5232460574
	h2.xyz                     0       1        -1.1693816161
	n2.xyz                     0       1      -109.5070757384
	hbr.xyz                    0       1     -2574.7361724856


	XYZ-files with optimized coordinates can be found in: optimized_xyzfiles



