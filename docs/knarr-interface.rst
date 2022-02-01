KNARR interface
======================================

State of the art nudged elastic band (NEB) calculations are possible in ASH using any theory level via an interface to the KNARR code.

.. code-block:: python

    def NEB(reactant=None, product=None, theory=None, images=None, interpolation=None, CI=None, free_end=None, restart_file=None,
            conv_type=None, tol_scale=None, tol_max_fci=None, tol_rms_fci=None, tol_max_f=None, tol_rms_f=None,
            tol_turn_on_ci=None, ActiveRegion=False, actatoms=None, runmode='serial', printlevel=1,
            idpp_maxiter=None):

TODO: Document options

################################################################################
Examples
################################################################################

10-image NEB calculation at the XTB level of theory

.. code-block:: python

    from ash import *

    numcores=2

    ################################################
    # Defining reactant and product ASH fragment
    #################################################
    charge=0
    mult=1
    react=Fragment(xyzfile="react.xyz", charge=charge, mult=mult)
    prod=Fragment(xyzfile="prod.xyz", charge=charge, mult=mult)


    #Theory to use for NEB
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library')

    #Run NEB to find saddlepoint. Returns saddlepoint as ASH fragment
    SP = NEB(reactant=react, product=prod, theory=xtbcalc, images=10, charge=charge


Restarting a calculation with user-defined path-file. 
Here, using the *restart_file* option to the NEB we read in a previous Knarr path-file ("knarr_MEP.xyz") instead of doing the regular IDPP interpolation
This file must have contain the coordinates of the same number of images (here 10) as number of images specified.
The file from a previously unconverged NEB calculation can be used or a converged MEP from a calculation at another level of theory.

.. code-block:: python

    from ash import *

    numcores=2

    ################################################
    # Defining reactant and product ASH fragment
    #################################################
    charge=0
    mult=1
    react=Fragment(xyzfile="react.xyz", charge=charge, mult=mult)
    prod=Fragment(xyzfile="prod.xyz", charge=charge, mult=mult)


    #Theory to use for NEB
    xtbcalc = xTBTheory(xtbmethod='GFN2', runmode='library', restart_file="knarr_MEP.xyz")

    #Run NEB to find saddlepoint. Returns saddlepoint as ASH fragment
    SP = NEB(reactant=react, product=prod, theory=xtbcalc, images=10)
