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
	Tgen: 0.01
	Tvar: 1e-09
	Orbital space of CAS(25,37) used for ICE-CI step
	Num generator CFGs: 28982
	Num CFGS after S+D: 29197

	Table of natural occupation numbers

	Orbital   MP2natorbs ICE-nat-occ
	----------------------------------------
	0            2.0000    2.0000
	1            2.0000    2.0000
	2            2.0000    2.0000
	3            2.0000    2.0000
	4            1.9854    1.9901
	5            1.9825    1.9882
	6            1.9791    1.9866
	7            1.9756    1.9845
	8            1.9720    1.9820
	9            1.9696    1.9810
	10           1.9666    1.9786
	11           1.9662    1.9777
	12           1.9639    1.9764
	13           1.9626    1.9754
	14           1.9584    1.9732
	15           1.9560    1.9723
	16           0.9950    0.9988
	17           0.0270    0.0197
	18           0.0249    0.0184
	19           0.0231    0.0176
	20           0.0210    0.0161
	21           0.0200    0.0158
	22           0.0174    0.0145
	23           0.0172    0.0142
	24           0.0166    0.0137
	25           0.0162    0.0136
	26           0.0143    0.0126
	27           0.0133    0.0120
	28           0.0119    0.0089
	29           0.0099    0.0060
	...

	Recommended active spaces based on ICE-CI natural occupations:
	Minimal (1.95,0.05): CAS(1,1)
	Medium1 (1.98,0.02): CAS(13,7)
	Medium2 (1.985,0.015): CAS(19,15)
	Medium3 (1.99,0.01): CAS(23,23)
	Medium4 (1.992,0.008): CAS(25,25)
	Large (1.995,0.005): CAS(25,30)
	Orbital file to use for future calculations: orca.gbw
	Note: orbitals are new natural orbitals formed from the ICE-CI density matrix





