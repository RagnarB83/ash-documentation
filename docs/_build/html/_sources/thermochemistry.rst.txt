High-level Thermochemistry in ASH
======================================

ASH contains various high-level thermochemistry protocols (a.k.a. composite methods) and convenient functionality to use them automatically on a set of compounds.
See also :doc:`module_benchmarking` and :doc:`module_workflows` page.

ASH can also conveniently calculate thermodynamic corrections via the RRHO approximaation

##############################################################################
thermochemistry corrections
##############################################################################

Thermochemistry corrections are automatically calculated (return object) when either a Numfreq or Anfreq job is requested.

.. code-block:: python

    thermochem_an = AnFreq(theory=ORCAcalc, fragment=HF_frag)
    thermochem_num = NumFreq(theory=ORCAcalc, fragment=HF_frag)

    print("Thermochem property dict:", thermochem_an)
    print("ZPVE (Eh) : ", thermochem_an['ZPVE'])

A dictionary containing various properties is returned (dictionary keys):
(frequencies, ZPVE, E_trans, E_rot, E_vib, E_tot, TS_trans, TS_rot, TS_vib, TS_el, vibenergycorr, Hcorr, Gcorr, TS_tot)

Alternatively, the thermochemcalc function can be called directly.

.. code-block:: python

    def thermochemcalc(vfreq,atoms,fragment, multiplicity, temp=298.18,pressure=1.0):

This function calculates the thermodynamic corrections from a list of available frequencies, number of atoms, ASH fragment object and spin multiplicity.
The temperature (default: 298.15 K) and pressure (default: 1.0 atm) can be specified.

.. code-block:: python

    h2o_frag = Fragment(xyzfile="h2o.xyz")
    #Manually defined frequencies for system
    frequencies=[1600.1, 2300.2, 2400.3]
    thermochemcalc(frequencies,3,h2o_frag, 1, temp=298.18, pressure=1.0)


##########################################################################################
thermochemprotocol functions: Automatic Opt+Freq+HL-theory
##########################################################################################

Reaction example:

.. code-block:: python

    from ash import *

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
    ORCAopt = ORCATheory(orcadir=orcadir, orcasimpleinput=simpleinput, orcablocks=blockinput, numcores=numcores)

    thermochemprotocol_reaction(Opt_theory=ORCAopt, SP_theory=DLPNO_CC_CBS, 
    fraglist=specieslist, stoichiometry=stoichiometry, orcadir=orcadir, numcores=numcores)

Single fragment example:

.. code-block:: python

    H2=Fragment(xyzfile='h2.xyz')
    #Defining theory for Opt+Freq step in thermochemprotocol
    orcadir='/Applications/orca_4.2.1'
    simpleinput="! B3LYP D3BJ def2-TZVP TightSCF Grid5 Finalgrid6"
    blockinput="""
    %scf maxiter 200 end
    """
    ORCAobject = ORCATheory(orcadir=orcadir, orcasimpleinput=simpleinput, orcablocks=blockinput, numcores=numcores)s
    thermochemprotocol_single(fragment=H2, Opt_theory=ORCAobject, SP_theory=DLPNO_CC_CBS, orcadir=orcadir, numcores=numcores)


Example with additional SP_theory workflow arguments:

.. code-block:: python

    H2=Fragment(xyzfile='h2.xyz')
    #Defining theory for Opt+Freq step in thermochemprotocol
    orcadir='/Applications/orca_4.2.1'
    simpleinput="! B3LYP D3BJ def2-TZVP TightSCF Grid5 Finalgrid6"
    blockinput="""
    %scf maxiter 200 end
    """
    ORCAobject = ORCATheory(orcadir=orcadir, orcasimpleinput=simpleinput, orcablocks=blockinput, numcores=numcores)
    DLPNO_CC_CBS_SP_args = {'cardinals' : [2,3], "basisfamily" : "def2", 'stabilityanalysis' : True, 
    'pnosetting' : 'extrapolation', 'pnoextrapolation' : [5,6], 'CVSR' : True,
                    'memory' : 5112, 'extrablocks' : "%scf\ndirectresetfreq 1\nend\n", 
                    'extrainputkeyword' : 'Slowconv'}
    thermochemprotocol_reaction(fraglist=[H2], stoichiometry=[1], Opt_theory=ORCAobject, 
    SP_theory=DLPNO_CC_CBS, workflow_args=DLPNO_CC_CBS_SP_args, orcadir=orcadir, numcores=numcores)



