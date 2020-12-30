

High-level Thermochemistry in ASH
======================================

ASH includes various high-level thermochemistry protocols (a.k.a. composite methods) and include convenient functionality to use them automatically on a set of compounds.

See also Benchmarking page and Workflows page.


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

    thermochemprotocol(Opt_theory=ORCAobject, SPprotocol='DLPNO-W1', workflow_args=fraglist=[H2], stoichiometry=[1], orcadir=orcadir, numcores=numcores)


Example with additional SPprotocol workflow arguments:

.. code-block:: python

    H2=Fragment(xyzfile='h2.xyz')
    #Defining theory for Opt+Freq step in thermochemprotocol
    orcadir='/Applications/orca_4.2.1'
    simpleinput="! B3LYP D3BJ def2-TZVP TightSCF Grid5 Finalgrid6"
    blockinput="""
    %scf maxiter 200 end
    """
    ORCAobject = ORCATheory(orcadir=orcadir, orcasimpleinput=simpleinput, orcablocks=blockinput, nprocs=numcores)
    DLPNO_CC_CBS_SP_args = {'cardinals' : [2,3], "basisfamily" : "def2", 'stabilityanalysis' : True, 'pnosetting' : 'extrapolation', 'pnoextrapolation' : [5,6], 'CVSR' : True,
                    'memory' : 5112, 'extrablocks' : "%scf\ndirectresetfreq 1\nend\n", 'extrainputkeyword' : 'Slowconv'}
    thermochemprotocol(Opt_theory=ORCAobject, SPprotocol='DLPNO-W1', workflow_args=DLPNO_CC_CBS_SP_args, fraglist=[H2], stoichiometry=[1], orcadir=orcadir, numcores=numcores)






##############################################################################
Available High-level SinglePoint Protocols
##############################################################################
These high-level singlepoint energy protocols can either be called on their own (see below) or used in the SPprotocol keyword argument in thermochemprotocol (see above).


Available protocols are (all use ORCA):

**DLPNO_CC_CBS_SP**

.. code-block:: python

    def DLPNO_CC_CBS_SP(fragment=None, cardinals = [2,3], basisfamily="def2", charge=None, orcadir=None, mult=None, stabilityanalysis=False,
        numcores=1, CVSR=False, memory=5000, pnosetting='NormalPNO', pnoextrapolation=[5,6], T1=False, scfsetting='TightSCF',
        extrainputkeyword='', extrablocks='', **kwargs):

Example:

.. code-block:: python

    H2=Fragment(xyzfile='h2.xyz')
    DLPNO_CC_CBS_SP(fragment=H2, cardinals = [2,3], basisfamily="def2", charge=0 orcadir=/opt/orca_4.2.1, mult=1, stabilityanalysis=False,
    numcores=1, CVSR=False, memory=5000, pnosetting='extrapolation', pnoextrapolation=[5,6], T1=False, scfsetting='TightSCF')

The code above defines an H2 fragment (from h2.xyz) and runs a single-point DLPNO-CCSD(T)/CBS calculation, utilizing basis-set and PNO extrapolation.
Cardinals=[2,3] and basisfamily="def2" means that the def2-SVP and def2-TZVPP basis sets will be used.
pnosetting="extrapolation" and pnoextrapolation=[5,6] means that the DLPNO-calculations will be run using 2 different TCutPNO cutoffs and then extrapolated to the PNO limit.

- Cardinals can be set to [2,3] or [3,4].
- basisfamily can be set to "def2" (Ahlrichs basis sets) or "cc" (correlation consistent basis sets).
- If a heavy element is chosen (heavier than Kr) then the cc-pVnZ-PP and corresponding ECP will be used for the heavy element.
- T1 option utilizes iterative triples, i.e. DLPNO-CCSD(T1) instead (more accurate, more expensive).
- CVSR adds a Core-Valence-Scalar-Relativistic correction (more accurate, more expensive).


TO BE DOCUMENTED:

- **W1theory_SP**
- **W1F12theory_SP**
- **DLPNO_W1theory_SP**
- **DLPNO_W1F12theory_SP**
- **DLPNO_F12_SP**
- **DLPNO_W2theory_SP**


