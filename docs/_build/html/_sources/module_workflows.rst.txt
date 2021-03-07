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