

High-level Thermochemistry in ASH
======================================

ASH contains various high-level thermochemistry protocols (a.k.a. composite methods) and convenient functionality to use them automatically on a set of compounds.
See also :doc:`benchmarking` and :doc:`workflows` page.

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
(freqs, ZPVE, E_trans, E_rot, E_vib, E_tot, TS_trans, TS_rot, TS_vib, TS_el, vibenergycorr, Hcorr, Gcorr, TS_tot)

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

    thermochemprotocol_reaction(Opt_theory=ORCAopt, SP_theory=DLPNO_CC_CBS, fraglist=specieslist, stoichiometry=stoichiometry, orcadir=orcadir, numcores=numcores)

Single fragment example:

.. code-block:: python

    H2=Fragment(xyzfile='h2.xyz')
    #Defining theory for Opt+Freq step in thermochemprotocol
    orcadir='/Applications/orca_4.2.1'
    simpleinput="! B3LYP D3BJ def2-TZVP TightSCF Grid5 Finalgrid6"
    blockinput="""
    %scf maxiter 200 end
    """
    ORCAobject = ORCATheory(orcadir=orcadir, orcasimpleinput=simpleinput, orcablocks=blockinput, nprocs=numcores)s
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
    ORCAobject = ORCATheory(orcadir=orcadir, orcasimpleinput=simpleinput, orcablocks=blockinput, nprocs=numcores)
    DLPNO_CC_CBS_SP_args = {'cardinals' : [2,3], "basisfamily" : "def2", 'stabilityanalysis' : True, 'pnosetting' : 'extrapolation', 'pnoextrapolation' : [5,6], 'CVSR' : True,
                    'memory' : 5112, 'extrablocks' : "%scf\ndirectresetfreq 1\nend\n", 'extrainputkeyword' : 'Slowconv'}
    thermochemprotocol_reaction(fraglist=[H2], stoichiometry=[1], Opt_theory=ORCAobject, SP_theory=DLPNO_CC_CBS, workflow_args=DLPNO_CC_CBS_SP_args, orcadir=orcadir, numcores=numcores)

##############################################################################
Available High-level SinglePoint Protocols
##############################################################################
These high-level singlepoint energy protocols can either be called on their own (see below) or used as SP_theory keyword argument in thermochemprotocol (see above).
All of these protocols use the ORCA quantum chemistry code and give the 0 K electronic energy.


**DLPNO_CC_CBS_SP**

This workflow carries out multiple ORCA calculations for the given geometry and extrapolates to the DLPNO-CCSD(T)/CBS limit.
This workflow is flexible and features multiple ways of approaching the CBS limit and the PNO limit.
Various options affecting the accuracy, efficiency and robustness of the protocol can be chosen.
The basis set families: cc-pVnZ ('cc') and Ahlrichs def2 ('def2') can be chosen that are available for most of the periodic table (cc-pVnZ-PP for heavy elements).
A corevalence+scalar-relativistic correction (CVSR option) can be included upon request (important for high-accuracy atomization energies).
Atomic spin-orbit coupling is automatically included if system is an atom.

.. code-block:: python

    def DLPNO_CC_CBS(fragment=None, cardinals = [2,3], basisfamily="def2", charge=None, orcadir=None, mult=None, stabilityanalysis=False,
        numcores=1, CVSR=False, CVbasis="W1-mtsmall", memory=5000, pnosetting='NormalPNO', pnoextrapolation=[5,6], T1=False, scfsetting='TightSCF',
        extrainputkeyword='', extrablocks='', **kwargs):

Example:

.. code-block:: python

    N2=Fragment(xyzfile='n2.xyz')
    DLPNO_CC_CBS(fragment=N2, cardinals = [2,3], basisfamily="def2", charge=0, orcadir='/opt/orca_4.2.1', mult=1, stabilityanalysis=False,
    numcores=1, CVSR=False, memory=5000, pnosetting='extrapolation', pnoextrapolation=[5,6], T1=False, scfsetting='TightSCF')

The example above defines an N2 fragment (from n2.xyz) and runs multiple DLPNO-CCSD(T) calculations, utilizing basis-set and PNO extrapolation to give a final CCSD(T)/CBS estimate.
Cardinals=[2,3] and basisfamily="def2" means that the def2-SVP and def2-TZVPP basis sets will be used and extrapolated to the basis set limit.
pnosetting="extrapolation" and pnoextrapolation=[5,6] means that the DLPNO-calculations will be run using 2 different TCutPNO cutoffs and then extrapolated to the PNO limit.

- Cardinals can be set to [2,3] or [3,4].
- basisfamily can be set to "def2" (Ahlrichs basis sets) or "cc" (correlation consistent basis sets).
- If a heavy element is chosen (heavier than Kr) and basisfamily="cc", then the cc-pVnZ-PP and corresponding ECP will be used for the heavy element.
- T1 option utilizes iterative triples, i.e. DLPNO-CCSD(T1) instead (more accurate, more expensive).
- CVSR adds a Core-Valence-Scalar-Relativistic correction (more accurate, more expensive). The correction is performed at the DLPNO-CCSD(T) level (hardcoded to NormalPNO) using the W1-mtsmall basis set.
  Note: If "W1-mtsmall" is not available for the element involved you may have to provide an appropriate core-valence basis set keyword via CVbasis="basisname" .

TO BE DOCUMENTED:

- **W1theory**
- **W1F12theory**
- **DLPNO_W1theory**
- **DLPNO_W1F12theory**
- **DLPNO_F12**
- **DLPNO_W2theory**


