

High-level Thermochemistry in ASH
======================================

ASH includes various high-level thermochemistry protocols (a.k.a. composite methods) and include convenient functionality to use them automatically on a set of compounds.

See aslo Benchmarking page and Workflows page.


##############################################################################
thermochemprotocol: Automatic Opt+Freq+HL-theory
##############################################################################

The thermochemprotocol function performs a geometry optimization, frequency calculation and high-level single-point protocol.
While intended for calculation of a molecular reaction it can be used for single molecules as well.

The reaction is defined via a list of defined fragments and stoichiometry, a theory object for Opt+Freq steps is defined (Opt_theory)
and then a protocol for the high-level single-point level is chosen (SPprotocol).
The available high-level single-point calculations are defined below.

Reaction example:

.. code-block:: python

    from ash import *
    settings_ash.init() #initialize

    orcadir='/opt/orca_4.2.1'
    numcores=4
    N2=Fragment(xyzfile="n2.xyz", charge=0, mult=1)
    H2=Fragment(xyzfile="h2.xyz", charge=0, mult=1)
    NH3=Fragment(xyzfile="nh3.xyz", charge=0, mult=1)

    ##Defining reaction##
    # List of species from reactant to product
    specieslist=[N2, H2, NH3] #Use same order as stoichiometry

    #Equation stoichiometry : negative integer for reactant, positive integer for product
    # Example: N2 + 3H2 -> 2NH3  reaction should be:  [1,3,-2]
    stoichiometry=[-1, -3, 2] #Use same order as specieslist
    ##
    #Defining theory for Opt+Freq step in thermochemprotocol
    orcadir='/Applications/orca_4.2.1'
    simpleinput="! B3LYP D3BJ def2-TZVP TightSCF Grid5 Finalgrid6"
    blockinput="""
    %scf maxiter 200 end
    """
    ORCAopt = ORCATheory(orcadir=orcadir, orcasimpleinput=simpleinput, orcablocks=blockinput, nprocs=numcores)

    thermochemprotocol(Opt_theory=ORCAopt, SPprotocol='W1', fraglist=specieslist, stoichiometry=stoichiometry, orcadir=orcadir, numcores=numcores)

Single fragment example:

.. code-block:: python

    H2=Fragment(xyzfile='h2.xyz')
    #Defining theory for Opt+Freq step in thermochemprotocol
    orcadir='/Applications/orca_4.2.1'
    simpleinput="! B3LYP D3BJ def2-TZVP TightSCF Grid5 Finalgrid6"
    blockinput="""
    %scf maxiter 200 end
    """
    ORCAobject = ORCATheory(orcadir=orcadir, orcasimpleinput=simpleinput, orcablocks=blockinput, nprocs=numcores)

    thermochemprotocol(Opt_theory=ORCAobject, SPprotocol='DLPNO-W1', fraglist=[H2], stoichiometry=[1], orcadir=orcadir, numcores=numcores)

##############################################################################
Available High-level SinglePoint Protocols
##############################################################################

Available protocols are (all use ORCA):

TO BE DOCUMENTED

- W1theory_SP
- W1F12theory_SP
- DLPNO_W1theory_SP
- DLPNO_W1F12theory_SP
- DLPNO_F12_SP
- DLPNO_W2theory_SP
- DLPNO_CC_CBS_SP