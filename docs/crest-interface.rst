CREST interface
======================================

ASH features a simple interface to the powerful conformational sampling program `crest <https://xtb-docs.readthedocs.io/en/latest/crest.html>`_ by the Grimme group.


By providing an ASH fragment object to the **call_crest** function, the Crest xTB-metadynamics-based conformational sampling procedure is invoked.
The output from crest is written to standard output. If successful, Crest will create a file crest_conformers.xyz
that can be directly read into ASH for further processing or further calculations.
This allows one to write a multi-step workflow of which the crest-procedure is one of many steps.

.. code-block:: python

    #Function to call crest
    def call_crest(fragment=None, xtbmethod=None, crestdir=None,charge=None, mult=None, solvent=None, energywindow=6, numcores=1,
                   constrained_atoms=None, forceconstant_constraint=0.5)
    #Function to grab conformers. Returns list of conformers and list of xtb energies
    def get_crest_conformers()


**call_crest** requires one to specify: an ASH fragment, xtbmethod (GFN1-xTB or GFN2-xTB ), location of crest directory, charge, multiplicity.

Optional keywords are: solvent, energywindow (default 6), numcores (default 1), constrained_atoms (list of integers) and the value of the force-consstant.

If you specify a list of constrained atoms then ASH will create an .xcontrol file that defines the constraints according to `crest Example Applications <https://xtb-docs.readthedocs.io/en/latest/crestxmpl.html>`_.




################################################################################
Example workflow 1. Call crest to get low-energy conformers as ASH fragments.
################################################################################
.. code-block:: python

    from ash import *

    crestdir='/opt/crest'
    numcores=24

    #0. Starting structure and charge and mult
    molecule = Fragment(xyzfile="ethanol.xyz", charge=0, mult=1)

    #1. Calling crest and getting list of conformer fragments and energies
    list_conformer_frags, xtb_energies = call_crest(fragment=molecule, xtbmethod='GFN2-xTB', crestdir=crestdir, 
        numcores=numcores)

    print("list_conformer_frags:", list_conformer_frags)
    print("")
    print("Crest Conformer Searches done. Found {} conformers".format(len(xtb_energies)))
    print("xTB energies: ", xtb_energies)


################################################################################
confsampler_protocol : Automatic Crest+DFTopt+DLPNO-CCSD(T) workflow
################################################################################

It is also possible to call the **confsampler_protocol** function that carries out an automatic multi-step workflow
at various levels of theory.

1. conformational sampling using crest and GFN-xTB (**low-level** theory).
2. Geometry optimizations for each low-energy conformer at a **medium-level** of theory (typically DFT using e.g. ORCATheory)
3. **High-level** single-point calculation (e.g. DLPNO-CCSD(T)/CBS using e.g. ORCA_CC_CBS_Theory)

.. code-block:: python

    def confsampler_protocol(fragment=None, crestdir=None, xtbmethod='GFN2-xTB', MLtheory=None, 
                            HLtheory=None, numcores=1, charge=None, mult=None, crestoptions=None,
                            optimizer_maxiter=200):


.. code-block:: python

    from ash import *

    #
    crestdir='/opt/crest'
    numcores=4
    #Fragment to define
    frag=Fragment(xyzfile="ethanol.xyz", charge=0, mult=1)

    #Defining MLTheory: DFT optimization
    MLsimpleinput="! B3LYP D3BJ def2-TZVP TightSCF "
    MLblockinput="""
    %scf maxiter 200 end
    """
    ML_B3LYP = ORCATheory(orcasimpleinput=MLsimpleinput, orcablocks=MLblockinput, numcores=numcores)
    #Defining HLTheory: DLPNO-CCSD(T)/CBS
	HL_CC = ORCA_CC_CBS_Theory(elements=frag.elems, cardinals = [2,3], basisfamily="def2", DLPNO=True, 
                  pnosetting='extrapolation', pnoextrapolation=[6,7], numcores=numcores)

    #Call confsampler_protocol
    confsampler_protocol(fragment=frag, crestdir=crestdir, xtbmethod='GFN2-xTB', MLtheory=ML_B3LYP,
                             HLtheory=HL_CC, orcadir=orcadir, numcores=numcores)

Final result table of calculated conformers at 3 different theory levels:

.. code-block:: text

    =================
    FINAL RESULTS
    =================

     Conformer   xTB-energy    DFT-energy    HL-energy (Eh)
    ----------------------------------------------------------------
             0 -25.8392205500 -346.2939482921 -345.2965932205
             1 -25.8377914500 -346.2884905132 -345.2911748671
             2 -25.8358803400 -346.2818766960 -345.2848279253
             3 -25.8313250600 -346.2788608396 -345.2815202116
             4 -25.8307377800 -346.2788662649 -345.2815419285
             5 -25.8303374700 -346.2775476223 -345.2792917601
             6 -25.8300128900 -346.2776089771 -345.2794648759

     Conformer   xTB-energy    DFT-energy    HL-energy (kcal/mol)
    ----------------------------------------------------------------
             0  0.0000000000  0.0000000000  0.0000000000
             1  0.8967737821  3.4248079602  3.4000680178
             2  2.0960134034  7.5750408530  7.3828340833
             3  4.9544947374  9.4675192805  9.4584557521
             4  5.3230184983  9.4641148891  9.4448282319
             5  5.5742168139 10.2915756050 10.8568301896
             6  5.7778938373 10.2530749008 10.7481984235