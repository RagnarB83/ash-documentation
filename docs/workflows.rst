

Workflows in ASH
======================================

As an ASH-script is pure Python, this allows one to easily create advanced workflows in a single script.

For example, a geometry optimization of a structure in on QM-program can easily be combined with a subsequent frequency job and this
can be followed by a subsequent higher-level single-point energy job using another QM-program even.

Simple for-loops can also be created to run multiple jobs with slightly different parameters (different theory level, different geometry etc.).

##############################################################################
Example: Running multiple single-point energies with different functionals
##############################################################################


.. code-block:: python

    from ash import *
    import sys
    import PES
    settings_ash.init() #initialize

    h2string="""
    H 0 0 0
    H 0 0 0.7
    """

    h2=Fragment(coordsstring=h2string)

    #List of functional keywords (strings) to loop over
    functionals=['BP86', 'B3LYP', 'TPSS', 'TPSSh', 'PBE0', 'BHLYP', 'CAM-B3LYP']

    #Dictionary to keep track of results
    energies_dict={}

    for functional in functionals:
        print("FUNCTIONAL: ", functional)
        orcadir='/opt/orca_4.2.1'
        #Appending functional keyword to the string-variable that contains the ORCA inputline
        input="! def2-SVP Grid5 Finalgrid6 tightscf slowconv " + functional
        blocks="""
        %scf
        maxiter 200
        end
        """
        #Defining/redefining ORCA theory. Does not need charge/mult keywords.
        ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput=input, orcablocks=blocks, nprocs=8, charge=0, mult=1)

        # Run single-point job
        energy = Singlepoint(theory=ORCAcalc, fragment=h2)
        #Adding energy to dictionary
        energies_dict[functional] = energy

        #Cleaning up after each job (not always necessary)
        ORCAcalc.cleanup()
        print("=================================")

    print("Dictionary with results:", energies_dict)
    print("")
    #Pretty formatted printing:
    print("")
    print(" Functional   Energy (Eh)")
    print("----------------------------")
    for func, e in energies_dict.items():
        print("{:10} {:13.10f}".format(func,e))


Producing a nice table of results:

.. code-block:: shell

     Functional   Energy (Eh)
    ----------------------------
    BP86       -1.1689426849
    B3LYP      -1.1642632249
    TPSS       -1.1734355861
    TPSSh      -1.1729787552
    PBE0       -1.1610065506
    BHLYP      -1.1624650247
    CAM-B3LYP  -1.1625896338


###########################################################################################
Example: Running conformer-sampling, geometry optimizations and High-level single-points
###########################################################################################
This example utilizes the interface to Crest to perform metadynamics-based conformational sampling from a starting geometry at a semi-empirical level of theory.
This is then followed by DFT geometry optimizations for each conformer found by the Crest procedure.
Finally high-level coupled cluster single-point calculations (here DLPNO-CCSD(T)) are performed for each conformer.


.. code-block:: python

    from ash import *
    import sys
    import PES
    from interface_crest import *
    settings_ash.init() #initialize

    orcadir='/opt/orca_4.2.1/'
    crestdir='/opt/crest'
    numcores=24

    #0. Starting structure and charge and mult
    molecule = Fragment(xyzfile="dmp.xyz")
    charge=0
    mult=1

    #1. Calling crest
    #call_crest(fragment=molecule, xtbmethod='GFN2-xTB', crestdir=crestdir, charge=charge, mult=mult, solvent='H2O', energywindow=6 )
    call_crest(fragment=molecule, xtbmethod='GFN2-xTB', crestdir=crestdir, charge=charge, mult=mult, numcores=numcores)

    #2. Grab low-lying conformers from crest_conformers.xyz as list of ASH fragments.
    list_conformer_frags, xtb_energies = get_crest_conformers()

    print("list_conformer_frags:", list_conformer_frags)
    print("")
    print("Crest Conformer Searches done. Found {} conformers".format(len(xtb_energies)))
    print("xTB energies: ", xtb_energies)

    #3. Run DFT geometry optimizations for each crest-conformer
    #ML Theory level. TODO: Run in ASH parallel instead of ORCA parallel?
    MLorcasimpleinput="! BP86 def2-SVP def2/J Grid5 Finalgrid6 tightscf"
    MLorcablocks="%scf maxiter 200 end"
    MLORCATheory = ORCATheory(orcadir=orcadir, charge=charge, mult=mult,
                        orcasimpleinput=MLorcasimpleinput, orcablocks=MLorcablocks, nprocs=numcores)

    DFT_energies=[]
    print("")
    for index,conformer in enumerate(list_conformer_frags):
        print("")
        print("Performing DFT Geometry Optimization for Conformer ", index)
        geomeTRICOptimizer(fragment=conformer, theory=MLORCATheory, coordsystem='tric')
        DFT_energies.append(conformer.energy)
        #Saving ASH fragment and XYZ file for each DFT-optimized conformer
        os.rename('Fragment-optimized.ygg', 'Conformer{}_DFT.ygg'.format(index))
        os.rename('Fragment-optimized.xyz', 'Conformer{}_DFT.xyz'.format(index))

    print("")
    print("DFT Geometry Optimization done")
    print("DFT_energies: ", DFT_energies)

    #4.Run high-level DLPNO-CCSD(T). Ash should now have optimized conformers
    HLorcasimpleinput="! DLPNO-CCSD(T) def2-SVP def2-SVP/C tightscf"
    HLorcablocks="""
    %scf
    maxiter 200
    end
    %mdci
    maxiter 100
    end
    """

    HLORCATheory = ORCATheory(orcadir=orcadir, charge=charge, mult=mult,
                        orcasimpleinput=HLorcasimpleinput, orcablocks=HLorcablocks, nprocs=numcores)
    HL_energies=[]
    for index,conformer in enumerate(list_conformer_frags):
        print("")
        print("Performing High-level calculation for DFT-optimized Conformer ", index)
        HLenergy = Singlepoint(theory=HLORCATheory, fragment=conformer)
        HL_energies.append(HLenergy)


    print("")
    print("=================")
    print("FINAL RESULTS")
    print("=================")

    #Printing total energies
    print("")
    print(" Conformer   xTB-energy    DFT-energy    HL-energy (Eh)")
    print("----------------------------------------------------------------")

    min_xtbenergy=min(xtb_energies)
    min_dftenergy=min(DFT_energies)
    min_HLenergy=min(HL_energies)

    for index,(xtb_en,dft_en,HL_en) in enumerate(zip(xtb_energies,DFT_energies, HL_energies)):
        print("{:10} {:13.10f} {:13.10f} {:13.10f}".format(index,xtb_en, dft_en, HL_en))

    print("")
    #Printing relative energies
    min_xtbenergy=min(xtb_energies)
    min_dftenergy=min(DFT_energies)
    min_HLenergy=min(HL_energies)
    harkcal = 627.50946900
    print(" Conformer   xTB-energy    DFT-energy    HL-energy (kcal/mol)")
    print("----------------------------------------------------------------")
    for index,(xtb_en,dft_en,HL_en) in enumerate(zip(xtb_energies,DFT_energies, HL_energies)):
        rel_xtb=(xtb_en-min_xtbenergy)*harkcal
        rel_dfT=(dft_en-min_dftenergy)*harkcal
        rel_HL=(HL_en-min_HLenergy)*harkcal
        print("{:10} {:13.10f} {:13.10f} {:13.10f}".format(index,rel_xtb, rel_dfT, rel_HL))

    print("")
    print("Workflow done!")

Final result table of calculated conformers at 3 different theory levels:

.. code-block:: shell

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
