Workflows module
======================================


#####################
Thermochemprotocols
#####################

The **thermochemprotocol_reaction** and **thermochemprotocol_single** functions can be used to
perform a multi-step Opt+Freq+HL-single-point protocol on either a reaction or a single species.


The thermochemprotocol_reaction is used for chemical reactions by giving a list of ASH fragments, stoichiometry and theory levels.

.. code-block:: python

    def thermochemprotocol_reaction(fraglist=None, stoichiometry=None, Opt_theory=None, SP_theory=None, orcadir=None, numcores=None, memory=5000,
                       workflow_args=None, analyticHessian=False)

while thermochemprotocol_single is used for a single fragment (thermochemprotocol_reaction calls thermochemprotocol_single).

.. code-block:: python

    def thermochemprotocol_single(fragment=None, Opt_theory=None, SP_theory=None, orcadir=None, numcores=None, memory=5000,
                       workflow_args=None, analyticHessian=True, temp=298.15, pressure=1.0):


The reaction is defined via a list of defined fragments and stoichiometry, a theory object for Opt+Freq steps is defined (Opt_theory)
and then a protocol for the high-level single-point level is chosen (SP_theory).
The available high-level single-point calculations are defined later.

###################################
High-level single-point workflows
###################################

See :doc:`module_highlevel_workflows`


###################################
Automatic active-space selection
###################################

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

	fragment=Fragment(coordsstring=fragstring)

	activespace_dictionary = auto_active_space(fragment=fragment, orcadir=None, basis="def2-TZVP", charge=0, mult=1,
	    initial_orbitals='MP2', tgen=1.0)
	#Returns dictionary with various active_spaces based on thresholds

Output:

.. code-block:: shell

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




