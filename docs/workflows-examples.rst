Workflow examples in ASH
======================================

As an ASH-script is pure Python and the user has access to various functions for manipulating coordinates, create fragments,
call QM-code and calculate energy, etc. this allows one to easily create advanced workflows in a single script.


As a basic example, a geometry optimization of a structure with a QM method can easily be combined with a subsequent frequency job and this
can be followed by a subsequent higher-level single-point energy job (even with different QM programs).

- See Example 1 for such a workflow.

Another workflow might involve calculating all species of a chemical reaction with the reaction energy being the final result.

- Example 2a shows an example of this using a simple for-loop.
- Example 2b shows how the reaction energy could be calculated with a high-level thermochemistry protocol (thermochemprotocol function).

It can also be advantageous to run multiple jobs with slightly different parameters (different theory level, different geometry etc.)

- Example 3a shows how multiple single-point energies with different functionals on 1 geometry can be easily run.
- Example 3b shows how the same can be accomplished in a completely parallel fashion (Singlepoint_parallel)
- Example 4a shows how multiple single-point energies on multiple XYZ-files can be easily run.
- Example 4b shows how the same can be accomplished in a completely fashion (Singlepoint_parallel)
- Example 5 shows how one can calculate localized orbitals and create Cube files for a collection of XYZ-files or a multi-XYZ file.


An even more advanced workflow combines metadynamics-based conformational sampling (Crest procedure by Grimme) from a starting structure,
automatically performs DFT geometry optimizations for each conformer and finally evaluates a high-level single-point energy.

- See Example 6.


##############################################################################
Example 1 : Optimization + Frequency + HL-singlepoint
##############################################################################

.. code-block:: python

    from ash import *

    #Defining molecular fragment
    molstring="""
    N 0 0 0
    N 0 0 0.9
    """
    molecule=Fragment(coordsstring=molstring)
    #Defining ORCA object for optimization
    numcores=2
    orcadir='/opt/orca_4.2.1'
    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", nprocs=numcores)

    #Geometry optimization of molecule and ORCAcalc theory object.
    geomeTRICOptimizer(theory=ORCAcalc,fragment=molecule)

    #Numfreq job of molecule (contains optimized coordinates). A 2-point Hessian is requested in runmode serial.
    thermochem = NumFreq(fragment=molecule, theory=ORCAcalc, npoint=2, runmode='serial')

    #Single-point HL job on optimized geometry
    HLORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput="! DLPNO-CCSD(T) Extrapolate(2/3,def2) def2-QZVPP/C",
                orcablocks="", nprocs=numcores)
    HLenergy = Singlepoint(theory=HLORCAcalc, fragment=molecule)

    print("DLPNO-CCSD(T)/CBS//BP86/def2-SVP energy: ", HLenergy, "Eh")
    print("ZPVE: ", thermochem['ZPVE'], "Eh")

.. code-block:: shell

    DLPNO-CCSD(T)/CBS//BP86/def2-SVP energy:  -109.421012242536 Eh

#######################################################################################################
Example 2a : Direct calculation of Reaction Energy:  N\ :sub:`2` \  + 3H\ :sub:`2`\  → 2NH\ :sub:`3`\
#######################################################################################################

.. code-block:: python

    from ash import *

    #Defining all reaction species as ASH objects from XYZ-files
    N2=Fragment(xyzfile="n2.xyz")
    H2=Fragment(xyzfile="h2.xyz")
    NH3=Fragment(xyzfile="nh3.xyz")

    ##Defining reaction##
    # List of species from reactant to product
    specieslist=[N2, H2, NH3] #Use same order as stoichiometry

    #Equation stoichiometry : negative integer for reactant, positive integer for product
    # Example: N2 + 3H2 -> 2NH3  reaction should be:  [-1,-3,2]
    stoichiometry=[-1, -3, 2] #Use same order as specieslist
    ##
    numcores=1
    orcadir='/opt/orca_4.2.1'

    FinalEnergies=[]
    for molecule in specieslist:
        #Defining ORCA object.
        ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", nprocs=numcores)
        energy = Singlepoint(theory=ORCAcalc, fragment=molecule)
        #Storing energy as list. Energy is also stored as part of fragment.
        FinalEnergies.append(energy)
        ORCAcalc.cleanup()

    #Reaction Energy via list of total energies:
    ReactionEnergy(stoichiometry=stoichiometry, list_of_fragments=specieslist, list_of_energies=FinalEnergies)

    ##Reaction Energy via internal energies of fragment objects:
    #ReactionEnergy(stoichiometry=stoichiometry, list_of_fragments=specieslist)


.. code-block:: shell

      Reaction_energy: -65.12668956189346 kcalpermol


#######################################################################################################
Example 2b : Direct calculation of Reaction Energy with an Automatic Thermochemistry Protocol
#######################################################################################################

A more advanced feature is to run each fragment with a high-level thermochemistry protocol (using ORCA) and get the final
reaction energy with chemical accuracy. Here the coupled-cluster based W1 method is used as part of the
thermochemprotocol function. The protocol will run a DFT opt+Freq job (as defined via the ORCA-inputline string shown)
and then do the high-level W1 theory protocol on top (multiple CCSD, CCSD(T) jobs with extrapolation, core-valence, scalar relativistic and atomic spin-orbit corrections etc.)
This feature is in progress and will be made more userfriendly soon. Note that W1 is only doable for really small molecules (1-4 heavy atom systems are doable).

See :doc:`thermochemistry` for more information.

.. code-block:: python

    from ash import *

    #
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

    #ORCA theory inputline for Opt+Freq
    Opt_protocol_inputline="! B3LYP D3BJ def2-TZVP TightSCF Grid5 Finalgrid6"

    #Thermochemistry protocol
    thermochemprotocol(SP_theory='W1', fraglist=specieslist, stoichiometry=stoichiometry, orcadir=orcadir,
                        numcores=numcores, Opt_protocol_inputline=Opt_protocol_inputline)


Final output:

.. code-block:: shell

     Reaction_energy(ΔSCF):  -33.980155385058865
     Reaction_energy(ΔCCSD):  -6.937247193220541
     Reaction_energy(Δ(T)):  1.4333499904116154
     Reaction_energy(ΔCV+SR):  -0.07653672690344188
     Reaction_energy(ΔSO):  0.0
     Reaction_energy(ΔZPVE):  20.455727327700334
    ----------------------------------------------
     Reaction_energy(Total ΔE):  -19.104861987083417

The output shows the total reaction energy (0 K enthalpy) and the contribution from Hartree-Fock (SCF), singles-doubles excitations (ΔCCSD),
perturbative triples (Δ(T)), core-valence + scalar-relativistics (CV+SR), atomic spin-orbit coupling (ΔSO, here none), and zero-point
vibrational energy (ΔZPVE).
The agreement with experiment (-18.4 kcal/mol) is excellent.


############################################################################################
Example 3a : Running multiple single-point energies with different functionals (sequential)
############################################################################################


.. code-block:: python

    from ash import *

    h2string="""
    H 0 0 0
    H 0 0 0.7
    """

    h2=Fragment(coordsstring=h2string)

    #List of functional keywords (strings) to loop over. Need to be valid ORCA keywords.
    functionals=['BP86', 'B3LYP', 'TPSS', 'TPSSh', 'PBE0', 'BHLYP', 'CAM-B3LYP']

    #Dictionary to keep track of energies
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
        ORCAcalc = ORCATheory(orcadir=orcadir, orcasimpleinput=input, orcablocks=blocks, nprocs=4, charge=0, mult=1)

        # Run single-point job
        energy = Singlepoint(theory=ORCAcalc, fragment=h2)

        #Keep ORCA outputfile for each functional
        os.rename('orca-input.out', functional+'_orcajob.out')

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


############################################################################################
Example 3b : Running multiple single-point energies with different functionals (in parallel)
############################################################################################
The example in 3a ran each job sequentially, one after the other, according to the list of functional strings.
While ORCA parallelization was utilized, it may be more economical to run the jobs simultaneously instead, especially if there are lot of jobs to go through.
This can be accomplished using the Singlepoint_parallel function inside ASH.
Here Python multiprocessing (pool.map) is utilized.
In this case ORCA parallelization must be turned off as the parallelization strategies are not compatible.

.. code-block:: python

    from ash import *
    #Fragment
    h2string="""
    H 0 0 0
    H 0 0 0.7
    """
    h2=Fragment(coordsstring=h2string)

    #Single-point job parallelization
    #Case: Multiple theories
    orcadir='/opt/orca_4.2.1'
    #Creating multiple ORCA objects and storing in list: orcaobjects
    #Important: use a label (here functional-name)for the created ORCA object to distinguish jobs
    orcaobjects=[]
    for functional in ['B3LYP', 'BP86', 'PBE0']:
        ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput="! def2-SVP def2/J "+functional, orcablocks="", label=functional)
        orcaobjects.append(ORCAcalc)

    #Calling the Singlepoint_parallel function and providing list of fragments and list of theories:
    results = Singlepoint_parallel(fragments=[h2], theories=orcaobjects, numcores=4)

    #results is a dictionary of energies
    print("results :", results)

###########################################################################################
Example 4a : Running single-point energies on a collection of XYZ files (sequential)
###########################################################################################

.. code-block:: python

    from ash import *
    import glob
    #
    orcadir='/opt/orca_4.2.1'
    numcores=1
    #Directory of XYZ files. Can be full path or relative path.
    dir = './xyz_files'
    #Changing to dir
    os.chdir(dir)

    energies=[]
    for file in glob.glob('*.xyz'):
        print("XYZ-file:", file)
        mol=Fragment(xyzfile=file)
        ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", nprocs=1)
        energy = Singlepoint(theory=ORCAcalc, fragment=mol)
        print("Energy of file {} : {} Eh".format(file, energy))
        ORCAcalc.cleanup()
        energies.append(energy)
        print("")
    #Pretty print
    print(" XYZ-file             Energy (Eh)")
    print("-----------------------------------------------")
    for xyzfile, e in zip(glob.glob('*.xyz'),energies):
        print("{:20} {:>13.10f}".format(xyzfile,e))


Output:

.. code-block:: python

     XYZ-file             Energy (Eh)
    -----------------------------------------------
    h2.xyz               -1.1715257797
    h2o_MeOH.xyz         -192.0023991603
    O-O-dimer.xyz        -149.8555328055
    butane.xyz           -158.3248873844
    nh3.xyz              -56.5093301286
    n2.xyz               -109.4002969311
    hi.xyz               -298.3735362292
    h2o_strained.xyz     -76.2253312246


############################################################################################
Example 4b : Running single-point energies on a collection of XYZ files (parallel)
############################################################################################
The example in 4a ran each job sequentially, one after the other, according to the list of XYZ-files available.
While ORCA parallelization was utilized, it may be more economical to run the jobs simultaneously instead, especially if there are lot of XYZ-files.
This can be accomplished using the Singlepoint_parallel function inside ASH.
Here Python multiprocessing (pool.map) is utilized.
In this case ORCA parallelization must be turned off as the parallelization strategies are not compatible.

.. code-block:: python

    from ash import *
    import glob
    #
    orcadir='/opt/orca_4.2.1'
    ORCAcalc = ORCATheory(orcadir=orcadir, charge=0, mult=1, orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks="", nprocs=1)
    #Directory of XYZ files. Can be full path or relative path.
    dir = './xyz_files'

    molecules=[]
    #Creating list of ASH fragments from XYZ files. Using filename as label
    for file in glob.glob(dir+'/*.xyz'):
        print("XYZ-file:", file)
        basename=os.path.basename(file)
        label=os.path.splitext(basename)[0]
        molecule=Fragment(xyzfile=file,label=label)
        molecules.append(molecule)

    #Calling the Singlepoint_parallel function and providing list of fragments and list of theories:
    results = Singlepoint_parallel(fragments=molecules, theories=[ORCAcalc], numcores=4)

    #results is a dictionary of energies
    print("results :", results)



###########################################################################################################
Example 5 : Calculate localized orbitals and create Cube files for multiple XYZ files or an XYZ-trajectory
###########################################################################################################

Analyzing electronic structure along a reaction path (e.g. a NEB or IRC path) or a trajectory (optimization or MD)
can be useful to understand the nature of the reaction. The code below shows how this can be accomplished in ASH
via a workflow involving single-point DFT, orbital localization and Cube-file creation (via orca_plot).

TODO: Add centroid analysis

Using a collection of XYZ-files:

.. code-block:: python

    from ash import *
    import glob
    #
    orcadir='/opt/orca_4.2.1'
    numcores=1
    #Directory of XYZ files. Can be full path or relative path.
    dir = '/home/bjornsson/ASH-DEV_GIT/testsuite/localized-orbital-IRC-workflow/calcs/images'
    #Changing to dir
    #os.chdir(dir)
    #Localization block in ORCA inputfile
    blockinput="""
    %loc
    LocMet IAOIBO
    end
    """

    #Looping over XYZ-files in directory, creating ASH fragments, running ORCA and calling orca_plot
    for file in sorted(glob.glob(dir+'/*.xyz')):
        basefile=os.path.basename(file)
        print("XYZ-file:", basefile)
        mol=Fragment(xyzfile=file)
        ORCAcalc = ORCATheory(orcadir=orcadir, charge=-1, mult=1, orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks=blockinput, nprocs=1)
        energy = Singlepoint(theory=ORCAcalc, fragment=mol)
        print("Energy of file {} : {} Eh".format(basefile, energy))
        locfile=basefile.split('.')[0]+'_calc.loc'
        os.rename('orca-input.loc', locfile)
        #Call ORCA_plot and create Cube file for specific MO in locfile: here alpha-MOs 13 and 17
        run_orca_plot(orcadir, locfile, 'mo', gridvalue=30, mo_operator=0, mo_number=13)
        run_orca_plot(orcadir, locfile, 'mo', gridvalue=30, mo_operator=0, mo_number=17)

        ORCAcalc.cleanup()
        print("")


Using a multi-XYZ file containing multiple sets of geometries (could be a NEB path, MD/Opt trajectory, XYZ animation etc.)

.. code-block:: python

    from ash import *
    import glob
    #
    orcadir='/opt/orca_4.2.1'
    numcores=1

    #Name of trajectory file containing multiple geometries (could be optimization traj, MD traj, NEB-path traj, Hessian XYZ animation etc.)
    #File should be in dir
    trajectoryfile="neb-ts_MEP_trj.xyz"

    blockinput="""
    %loc
    LocMet IAOIBO
    end
    """

    fraglist = get_molecules_from_trajectory(trajectoryfile)

    for index,frag in enumerate(fraglist):
        print("Frag :", index)
        ORCAcalc = ORCATheory(orcadir=orcadir, charge=-1, mult=1, orcasimpleinput="! BP86 def2-SVP def2/J", orcablocks=blockinput, nprocs=1)
        energy = Singlepoint(theory=ORCAcalc, fragment=frag)
        print("Energy of frag {} : {} Eh".format(index, energy))
        locfile='frag{}_calc.loc'.format(index)
        os.rename('orca-input.loc', locfile)
        #Call ORCA_plot and create Cube file for specific MO in locfile: here alpha-MOs 13 and 17
        run_orca_plot(orcadir, locfile, 'mo', gridvalue=30, mo_operator=0, mo_number=13)
        run_orca_plot(orcadir, locfile, 'mo', gridvalue=30, mo_operator=0, mo_number=17)

        ORCAcalc.cleanup()


###########################################################################################
Example 6 : Running conformer-sampling, geometry optimizations and High-level single-points
###########################################################################################
This example utilizes the interface to Crest to perform metadynamics-based conformational sampling from a starting geometry at a semi-empirical level of theory.
This is then followed by DFT geometry optimizations for each conformer found by the Crest procedure.
Finally high-level coupled cluster single-point calculations (here DLPNO-CCSD(T)/CBS extrapolations) are performed for each conformer.


.. code-block:: python

    from ash import *
    from interface_crest import *

    orcadir='/opt/orca_4.2.1/'
    crestdir='/opt/crest'
    numcores=24

    #0. Starting structure and charge and mult
    molecule = Fragment(xyzfile="ethanol.xyz")
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
    MLorcasimpleinput="! BP86 D3 def2-TZVP def2/J Grid5 Finalgrid6 tightscf"
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
    HLorcasimpleinput="! DLPNO-CCSD(T) Extrapolate(2/3,def2) def2-QZVPP/C tightscf TightPNO"
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


The manually defined workflow above can also be more conveniently run like this:

.. code-block:: python

    from ash import *

    #
    crestdir='/opt/crest'
    orcadir='/opt/orca_4.2.1'
    numcores=4
    #Fragment to define
    frag=Fragment(xyzfile="ethanol.xyz", charge=0, mult=1)

    #Defining MLTheory: DFT optimization
    orcadir='/opt/orca_4.2.1'
    MLsimpleinput="! B3LYP D3BJ def2-TZVP TightSCF Grid5 Finalgrid6"
    MLblockinput="""
    %scf maxiter 200 end
    """
    ML_B3LYP = ORCATheory(orcadir=orcadir, orcasimpleinput=MLsimpleinput, orcablocks=MLblockinput, nprocs=numcores, charge=frag.charge, mult=frag.mult)
    #Defining HLTheory: DLPNO-CCSD(T)/CBS
    HLsimpleinput="! DLPNO-CCSD(T) Extrapolate(2/3,def2) def2-QZVPP/C TightSCF"
    HLblockinput="""
    %scf maxiter 200 end
    """
    HL_CC = ORCATheory(orcadir=orcadir, orcasimpleinput=HLsimpleinput, orcablocks=HLblockinput, nprocs=numcores, charge=frag.charge, mult=frag.mult)

    #Call confsampler_protocol
    confsampler_protocol(fragment=frag, crestdir=crestdir, xtbmethod='GFN2-xTB', MLtheory=ML_B3LYP,
                             HLtheory=HL_CC, orcadir=orcadir, numcores=numcores, charge=frag.charge, mult=frag.mult)

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
